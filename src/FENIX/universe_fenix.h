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

#ifndef LMP_UNIVERSE_FENIX_H
#define LMP_UNIVERSE_FENIX_H

#include "universe.h"

namespace LAMMPS_NS {

class UniverseFenix : public Universe {
 public:
  MPI_Info fenix_config = MPI_INFO_NULL;
  int *spares_per_world = nullptr;
  int role;
  int status;
  
  // world including spare ranks
  MPI_Comm full_world = MPI_COMM_NULL;
  // resilient world not including spare ranks
  MPI_Comm res_world = MPI_COMM_NULL;

  double default_spares = -1;

  //If set, this is run instead of infile when restarting after failures
  std::string restart_file;

  UniverseFenix(class LAMMPS *, MPI_Comm);
  ~UniverseFenix() override;

  void add_world(const char *) override;

  bool parse_arg(int*, int, char**) override;

  void run() override;

 protected:
  void create_worlds() override;
};

}

#endif
