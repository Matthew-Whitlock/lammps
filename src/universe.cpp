// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "universe.h"
#include "input.h"

#include "error.h"
#include "memory.h"

#include <cstring>

using namespace LAMMPS_NS;

static constexpr int MAXLINE = 256;

/* ----------------------------------------------------------------------
   create & initialize the universe of processors in communicator
------------------------------------------------------------------------- */

Universe::Universe(LAMMPS *lmp, MPI_Comm communicator) : Pointers(lmp)
{
  uworld = uorig = communicator;
  external_comm = MPI_COMM_NULL;

  MPI_Comm_rank(uworld,&me);
  MPI_Comm_size(uworld,&nprocs);

  uscreen = stdout;
  ulogfile = nullptr;

  existflag = 0;
  nworlds = 0;
  procs_per_world = nullptr;
  root_proc = nullptr;

  memory->create(uni2orig,nprocs,"universe:uni2orig");
  for (int i = 0; i < nprocs; i++) uni2orig[i] = i;

  screenarg = logarg = partscreenarg = partlogarg = nullptr;
}

/* ---------------------------------------------------------------------- */

Universe::~Universe()
{
  if (screen && screen != stdout) {
    fclose(screen);
    screen = nullptr;
  }
  if (logfile) {
    fclose(logfile);
    logfile = nullptr;
  }
  if (nworlds > 1 && ulogfile){
    fclose(ulogfile);
    ulogfile = nullptr;
  }
  if (world != uworld) MPI_Comm_free(&world);
  if (uworld != uorig) MPI_Comm_free(&uworld);
  if (external_comm != MPI_COMM_NULL) MPI_Comm_free(&uorig);
  memory->destroy(procs_per_world);
  memory->destroy(root_proc);
  memory->destroy(uni2orig);
}

/* ----------------------------------------------------------------------
   Parse argument and its parameters if a valid universe arg, and increment
   start value based on number of parsed strings from arg.
   Returns true if this argument is handled by universe
------------------------------------------------------------------------- */

bool Universe::parse_arg(int* start, int narg, char** arg){
  int& iarg = *start;

  if (strcmp(arg[iarg],"-mpicolor") == 0 ||
      strcmp(arg[iarg],"-m") == 0) {
    if (iarg+2 > narg)
      error->universe_all(FLERR,"Invalid command-line argument");
    if (iarg != 1) error->universe_all(FLERR,"Invalid command-line argument");
    int color = std::stoi(arg[iarg+1]);
    this->multiapp_split(color);

    iarg += 2;

  } else if (strcmp(arg[iarg],"-partition") == 0 ||
             strcmp(arg[iarg],"-p") == 0) {
    if (iarg+2 > narg)
      error->universe_all(FLERR,"Invalid command-line argument");
    iarg++;
    while (iarg < narg && arg[iarg][0] != '-') {
      this->add_world(arg[iarg]);
      iarg++;
    }

  } else if (strcmp(arg[iarg],"-reorder") == 0 ||
             strcmp(arg[iarg],"-ro") == 0) {
    if (iarg+3 > narg)
      error->universe_all(FLERR,"Invalid command-line argument");
    if (nworlds>0)
      error->universe_all(FLERR,"Cannot use -reorder after -partition");
    this->reorder(arg[iarg+1],arg[iarg+2]);
    iarg += 3;

  } else if (strcmp(arg[iarg],"-screen") == 0 ||
             strcmp(arg[iarg],"-sc") == 0) {
    if (iarg+2 > narg)
      error->universe_all(FLERR,"Invalid command-line argument");
    screenarg = arg[iarg + 1];
    iarg += 2;

  } else if (strcmp(arg[iarg],"-log") == 0 ||
             strcmp(arg[iarg],"-l") == 0) {
    if (iarg+2 > narg)
      error->universe_all(FLERR,"Invalid command-line argument");
    logarg = arg[iarg + 1];
    iarg += 2;

  } else if (strcmp(arg[iarg],"-plog") == 0 ||
             strcmp(arg[iarg],"-pl") == 0) {
    if (iarg+2 > narg)
     error->universe_all(FLERR,"Invalid command-line argument");
    partlogarg = arg[iarg + 1];
    iarg += 2;

  } else if (strcmp(arg[iarg],"-pscreen") == 0 ||
             strcmp(arg[iarg],"-ps") == 0) {
    if (iarg+2 > narg)
     error->universe_all(FLERR,"Invalid command-line argument");
    partscreenarg = arg[iarg + 1];
    iarg += 2;

  } else {
    return false;
  }

  return true;
}

/* ----------------------------------------------------------------------
  Run if multiple apps were launched with one mpirun command.
    If so, the passed communicator (e.g. MPI_COMM_WORLD) is bigger than LAMMPS
    universe should be, so shrink it based on passed color.
  syntax: -mpicolor color
    color = integer for this app, different than any other app(s)
    Must be the first argument passed
  do the following:
    perform an MPI_Comm_split() to create a new LAMMPS-only subcomm
    NOTE: this assumes other app(s) make same call, else will hang!
    store comm that all apps belong to in external_comm
------------------------------------------------------------------------- */
void Universe::multiapp_split(int color)
{
  external_comm = uorig;
  MPI_Comm_split(external_comm,color,me,&uorig);
  uworld = uorig;
  MPI_Comm_rank(uworld, &me);
  MPI_Comm_size(uworld, &me);
}

/* ----------------------------------------------------------------------
   reorder universe processors
   create uni2orig as inverse mapping
   re-create uworld communicator with new ordering via Comm_split()
   style = "nth", arg = N
   move every Nth proc to end of rankings
   style = "custom", arg = filename
   file has nprocs lines with I J
   I = universe proc ID in original communicator uorig
   J = universe proc ID in reordered communicator uworld
------------------------------------------------------------------------- */

void Universe::reorder(char *style, char *arg)
{
  char line[MAXLINE] = {'\0'};

  if (uworld != uorig) MPI_Comm_free(&uworld);

  if (strcmp(style,"nth") == 0) {
    int n = utils::inumeric(FLERR,arg,false,lmp);
    if (n <= 0)
      error->universe_all(FLERR,"Invalid -reorder N value");
    if (nprocs % n)
      error->universe_all(FLERR,"Nprocs not a multiple of N for -reorder");
    for (int i = 0; i < nprocs; i++) {
      if (i < (n-1)*nprocs/n) uni2orig[i] = i/(n-1) * n + (i % (n-1));
      else uni2orig[i] = (i - (n-1)*nprocs/n) * n + n-1;
    }

  } else if (strcmp(style,"custom") == 0) {

    if (me == 0) {
      FILE *fp = fopen(arg,"r");
      if (fp == nullptr)
        error->universe_one(FLERR,fmt::format("Cannot open -reorder "
                                              "file {}: {}",arg,
                                              utils::getsyserror()));

      // skip header = blank and comment lines

      char *ptr;
      if (!fgets(line,MAXLINE,fp))
        error->one(FLERR,"Unexpected end of -reorder file");
      while (true) {
        if ((ptr = strchr(line,'#'))) *ptr = '\0';
        if (strspn(line," \t\n\r") != strlen(line)) break;
        if (!fgets(line,MAXLINE,fp))
          error->one(FLERR,"Unexpected end of -reorder file");
      }

      // read nprocs lines
      // uni2orig = inverse mapping

      int me_orig,me_new,rv;
      rv = sscanf(line,"%d %d",&me_orig,&me_new);
      if (me_orig < 0 || me_orig >= nprocs ||
          me_new < 0 || me_new >= nprocs || rv != 2)
        error->one(FLERR,"Invalid entry '{} {}' in -reorder "
                                     "file", me_orig, me_new);
      uni2orig[me_new] = me_orig;

      for (int i = 1; i < nprocs; i++) {
        if (!fgets(line,MAXLINE,fp))
          error->one(FLERR,"Unexpected end of -reorder file");
        rv = sscanf(line,"%d %d",&me_orig,&me_new);
        if (me_orig < 0 || me_orig >= nprocs ||
            me_new < 0 || me_new >= nprocs || rv != 2)
          error->one(FLERR,"Invalid entry '{} {}' in -reorder "
                                       "file", me_orig, me_new);
        uni2orig[me_new] = me_orig;
      }
      fclose(fp);
    }

    // bcast uni2org from proc 0 to all other universe procs

    MPI_Bcast(uni2orig,nprocs,MPI_INT,0,uorig);

  } else error->universe_all(FLERR,"Invalid command-line argument");

  // create new uworld communicator

  int ome,key;
  MPI_Comm_rank(uorig,&ome);
  for (int i = 0; i < nprocs; i++)
    if (uni2orig[i] == ome) key = i;

  MPI_Comm_split(uorig,0,key,&uworld);
  MPI_Comm_rank(uworld,&me);
  MPI_Comm_size(uworld,&nprocs);
}

/* ----------------------------------------------------------------------
   add 1 or more worlds to universe
   str == nullptr -> add 1 world with all procs in universe
   str = NxM -> add N worlds, each with M procs
   str = P -> add 1 world with P procs
------------------------------------------------------------------------- */

void Universe::add_world(const char *str)
{
  int n,nper;

  n = 1;
  nper = 0;

  if (str != nullptr) {
    // flag that this universe is partitioned (ie 'exists')

    existflag = 1;

    // check for valid partition argument

    bool valid = true;

    // str may not be empty and may only consist of digits or 'x'

    std::string part(str);
    if (part.size() == 0) valid = false;
    if (part.find_first_not_of("0123456789x") != std::string::npos) valid = false;

    if (valid) {
      std::size_t found = part.find_first_of('x');

      // 'x' may not be the first or last character

      if ((found == 0) || (found == (part.size() - 1))) {
        valid = false;
      } else if (found == std::string::npos) {
        nper = std::stoi(part);
      } else {
        n = std::stoi(part.substr(0,found));
        nper = std::stoi(part.substr(found+1));
      }
    }

    // require minimum of 1 partition with 1 processor

    if (n < 1 || nper < 1) valid = false;

    if (!valid)
      error->universe_all(FLERR, fmt::format("Invalid partition string '{}'", str));
  } else nper = nprocs;

  memory->grow(procs_per_world,nworlds+n,"universe:procs_per_world");
  memory->grow(root_proc,(nworlds+n),"universe:root_proc");

  for (int i = 0; i < n; i++) {
    procs_per_world[nworlds] = nper;
    if (nworlds == 0) root_proc[nworlds] = 0;
    else
      root_proc[nworlds] = root_proc[nworlds-1] + procs_per_world[nworlds-1];
    if (me >= root_proc[nworlds]) iworld = nworlds;
    nworlds++;
  }
}

/* ----------------------------------------------------------------------
   check if total procs in all worlds = procs in universe
------------------------------------------------------------------------- */

int Universe::consistent()
{
  int n = 0;
  for (int i = 0; i < nworlds; i++) n += procs_per_world[i];
  if (n == nprocs) return 1;
  else return 0;
}

/* ---------------------------------------------------------------------- */

void Universe::create(int helpflag)
{
  // if no partition command-line switch, universe is one world with all procs

  if(nworlds == 0) this->add_world(nullptr);

  // sum of procs in all worlds must equal total # of procs

  if(!this->consistent())
    error->universe_all(FLERR,"Processor partitions do not match number of allocated processors");

  // if no partition command-line switch, cannot use -pscreen option

  if(nworlds == 1 && partscreenarg)
    error->universe_all(FLERR,"Can only use -pscreen with multiple partitions");

  // if no partition command-line switch, cannot use -plog option

  if(nworlds == 1 && partlogarg)
    error->universe_all(FLERR,"Can only use -plog with multiple partitions");

  // create communicator(s) for lammps->world

  this->create_worlds();

  // set universe-level screen and logfile

  uscreen = screenarg ? nullptr : stdout;
  if (me == 0 && screenarg && strcmp(screenarg,"none") != 0) {
    uscreen = fopen(screenarg,"w");
    if (uscreen == nullptr)
      error->universe_one(FLERR,fmt::format("Cannot open universe screen file {}: {}",
                                            screenarg,utils::getsyserror()));
  }

  ulogfile = nullptr;
  const char* logname = logarg == nullptr ? "log.lammps" : logarg;
  if (me == 0 && strcmp(logname,"none") != 0) {
    ulogfile = fopen(logname, "w");
    bool skip_check = logarg == nullptr && helpflag;
    if(!skip_check && ulogfile == nullptr)
      error->universe_one(FLERR,fmt::format("Cannot open universe log file {}: {}",
                                            logname,utils::getsyserror()));
  }

  // set world-level screen and logfile

  int world_rank;
  MPI_Comm_rank(world, &world_rank);

  if (nworlds == 1) {
    //Simply inherit from universe
    screen = uscreen;
    logfile = ulogfile;
  } else {
    screen = logfile = infile = nullptr;

    const char* screen_prefix = partscreenarg ? partscreenarg : screenarg;
    if(!screen_prefix) screen_prefix = "screen";
    if(world_rank == 0 && strcmp(screen_prefix, "none") != 0) {
      std::string str = fmt::format("{}.{}", screen_prefix, iworld);
      screen = fopen(str.c_str(), "w");
      if (screen == nullptr)
        error->one(FLERR,"Cannot open screen file {}: {}",str,utils::getsyserror());
    }

    const char* log_prefix = partlogarg ? partlogarg : logname;
    if(world_rank == 0 && strcmp(log_prefix, "none") != 0) {
      std::string str = fmt::format("{}.{}", log_prefix, iworld);
      logfile = fopen(str.c_str(), "w");
      if (logfile == nullptr)
        error->one(FLERR,"Cannot open logfile {}: {}", str, utils::getsyserror());
    }
  }
}

/* ----------------------------------------------------------------------
   Construct individual worlds if multiple partitions requested
------------------------------------------------------------------------- */

void Universe::create_worlds(){
  if(nworlds == 1){
    world = uworld;
  } else {
    MPI_Comm_split(uworld, iworld, me, &world);
  }
}

/* ----------------------------------------------------------------------
   Run the input file across the universe.
------------------------------------------------------------------------- */

void Universe::run(){
  input->file();
}
