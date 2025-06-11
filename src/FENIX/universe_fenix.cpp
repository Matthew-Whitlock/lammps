/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "universe_fenix.h"
#include "memory.h"
#include "error.h"
#include "utils.h"
#include "input.h"
#include "comm.h"

#include "fenix.hpp"

#include <fmt/ranges.h>

namespace LAMMPS_NS {

/* ---------------------------------------------------------------------- */

UniverseFenix::UniverseFenix(LAMMPS *lmp, MPI_Comm comm) : Universe(lmp, comm){
  role = FENIX_ROLE_INITIAL_RANK;
  status = FENIX_SUCCESS;
}

/* ---------------------------------------------------------------------- */

UniverseFenix::~UniverseFenix() {
  if(full_world != MPI_COMM_NULL && full_world != uworld)
    MPI_Comm_free(&full_world);

  if(fenix_config != MPI_INFO_NULL) MPI_Info_free(&fenix_config);
  if(full_world != MPI_COMM_NULL && full_world != uorig && full_world != world)
    MPI_Comm_free(&full_world);
  memory->destroy(spares_per_world);
}

/* ---------------------------------------------------------------------- */

bool UniverseFenix::parse_arg(int* start, int narg, char** arg){
  int& iarg = *start;

  if(strcmp(arg[iarg], "-spares") == 0){
    if(++iarg >= narg)
      error->universe_all(FLERR,"Invalid command-line argument");

    std::string spares = arg[iarg++];
    if(!utils::is_double(spares))
      error->universe_all(FLERR,"Invalid command-line argument");

    default_spares = stof(spares);
    if(default_spares<0 ||
       (default_spares>1 && default_spares != std::trunc(default_spares)))
      error->universe_all(FLERR,"Invalid command-line argument, spares value not valid");
  
  } else if(strcmp(arg[iarg], "-restart-file") == 0){
    if(++iarg >= narg)
      error->universe_all(FLERR,"Invalid command-line argument");
    restart_file = arg[iarg++];

  } else {
    return Universe::parse_arg(start, narg, arg);
  }

  return true;
}

/* ----------------------------------------------------------------------
   Handle partitions that include a spare ranks specifier following an 's'
------------------------------------------------------------------------- */

void UniverseFenix::add_world(const char* str) {
  using namespace std;

  string spec = str ? str : "";
  string base_spec = spec;

  int spares = 0;
  double frac_spares = 0;

  size_t spares_idx = spec.find('s');
  string spares_spec;
  if(spares_idx != string::npos) {
    spares_spec = spec.substr(spares_idx+1);
    base_spec = spec.substr(0, spares_idx);
  }
  
  if(spares_spec.empty()){
    if(default_spares != -1) frac_spares = default_spares;
    else spares = -1;
  } else if(utils::is_double(spares_spec)){
    frac_spares = stof(spares_spec);
  } else {
    // Failed to parse, try just passing the whole spec to Universe.
    base_spec = spec;
    spares = -1;
  }

  if(frac_spares != 0){
    if(frac_spares < 0)
      error->universe_all(FLERR, fmt::format(
        "Invalid partition string '{}': spares must be >= 0",
        spec
      ));
    
    bool is_int = frac_spares == std::trunc(frac_spares);
    if(is_int){
      spares = frac_spares;
      frac_spares = 0;
    }
    if(frac_spares >= 1){
      error->universe_all(FLERR, fmt::format(
        "Invalid partition string '{}': fractional spares must be < 1",
        spec
      ));
    }
  }


  int old_nworlds = nworlds;

  // allow e.g. "-partition s5" to function as if no partition was passed
  // while still setting the number of spares.
  const char* ptr = base_spec.empty() ? nullptr : base_spec.c_str();
  Universe::add_world(ptr);
  
  memory->grow(spares_per_world, nworlds, "universe/fenix:spares_per_world");
  for(int i = old_nworlds; i < nworlds; i++){
    if(spares != 0) spares_per_world[i] = spares;
    else spares_per_world[i] = frac_spares * procs_per_world[i];

    if(spares_per_world[i] >= procs_per_world[i])
      error->universe_all(FLERR, fmt::format(
        "Invalid partition string '{}': {} spares is too many for world size {}",
        spec, spares_per_world[i], procs_per_world[i]
      ));

    if(frac_spares > 0 && spares_per_world[i] == 0 && me == 0)
      error->universe_warn(FLERR, fmt::format(
        "Requested {}% spares is 0 on world size {}",
        frac_spares*100, procs_per_world[i]
      ));
  }
}

/* ----------------------------------------------------------------------
   For worlds with spare ranks, initialize Fenix and replace the world.
------------------------------------------------------------------------- */

void UniverseFenix::create_worlds(){
  Universe::create_worlds();

  full_world = world;
  // TODO: update uworld to only have active ranks.
  //   Will also need to update uni2orig

  // Update procs per world to no longer include spares
  for(int i = 0; i < nworlds; i++){
    if(spares_per_world[i] > 0){
      procs_per_world[i] -= spares_per_world[i];
    }
  }

  // If my world had no spares specification, don't initialize Fenix on it
  if(spares_per_world[iworld] < 1) return;

  MPI_Info fenix_config;
  MPI_Info_create(&fenix_config);
  MPI_Info_set(fenix_config, "FENIX_RESUME_MODE", "THROW");

  Fenix_Init(&role, full_world, &res_world, nullptr, nullptr,
             spares_per_world[iworld], 0, fenix_config, &status);

  MPI_Info_free(&fenix_config);
  fenix_config = MPI_INFO_NULL;

  world = res_world;
}

/* ----------------------------------------------------------------------
   Catch any process failures during execution then rebuild and run the
   restart file (or repeat the input file if restart file not specified)
------------------------------------------------------------------------- */

void UniverseFenix::run(){
  bool should_throw = role == FENIX_ROLE_RECOVERED_RANK;
  while(true){ try {
    if(should_throw){
      should_throw = false;
      throw Fenix::CommException(world, status);
    }

    if(role == FENIX_ROLE_INITIAL_RANK ||
        restart_file.empty()){
      input->file();
    } else {
      input->file(restart_file.c_str());
    }

    Fenix_Finalize();
    break;
  } catch (const Fenix::CommException& e){
    try{
      world = res_world;

      if(comm->me == 0){
        int *faults, n_faults = Fenix_Process_fail_list(&faults);
        utils::logmesg(lmp, "Fenix recovering from rank failure(s): [{}]\n",
            fmt::join(faults, faults+n_faults, ", ")
        );
      }

      if(status != FENIX_SUCCESS){
        error->one(FLERR, "Fenix unable to recover, status: {}", status);
      }

      if(role == FENIX_ROLE_SURVIVOR_RANK){
        // Survivors recreate everything
        lmp->destroy();
        if(infile && infile != stdin){
          rewind(infile);
        }

        int narg = input->narg;
        char** arg = input->arg;
        delete input;
        input = new Input(lmp, narg, arg);

        lmp->create();
        lmp->post_create();
      }
    } catch (const Fenix::CommException& f){
      should_throw = true;
    }
  }}
}

}
