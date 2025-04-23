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

#ifndef LMP_UNIVERSE_H
#define LMP_UNIVERSE_H

#include "pointers.h"

namespace LAMMPS_NS {

class Universe : protected Pointers {
 public:
  MPI_Comm uworld;    // communicator for entire universe
  int me, nprocs;     // my place in universe

  FILE *uscreen;     // universe screen output
  FILE *ulogfile;    // universe logfile

  int existflag;           // 1 if universe exists due to -partition flag
  int nworlds;             // # of worlds in universe
  int iworld;              // which world I am in
  int *procs_per_world;    // # of procs in each world
  int *root_proc;          // root proc in each world

  MPI_Comm uorig;    // unordered communicator passed to LAMMPS instance
  int *uni2orig;     // proc I in universe uworld is
                     // proc uni2orig[I] in original communicator

  MPI_Comm external_comm;    // Original communicator passed to lammps,
                             //  if lammps was launched alongside multiple
                             //  MPI apps. Else MPI_COMM_NULL.

  char* screenarg;       // arguments for various configuration flags
  char* logarg;          //  or nullptr if not specified.
  char* partscreenarg;
  char* partlogarg;

  Universe(class LAMMPS *, MPI_Comm);
  ~Universe() override;

  void reorder(char *, char *);

  virtual void add_world(const char *);
  virtual int consistent();

  virtual bool parse_arg(int*, int, char**);   // handle one argument and its parameters
  virtual void create(int);                    // create universe and world communicators, screens, and logs
  virtual void run();                          // run all input

 protected:
  virtual void multiapp_split(int);            // split off just LAMMPS' portion of uorig if -mpicolor passed
  virtual void create_worlds();                // create world communicator(s) from universe
};

}    // namespace LAMMPS_NS

#endif
