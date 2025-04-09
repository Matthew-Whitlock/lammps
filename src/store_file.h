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

#ifndef LMP_STORE_FILE_H
#define LMP_STORE_FILE_H

#include "store_impl.h"
#include <cstdio>

namespace LAMMPS_NS {

class StoreFile : public StoreImpl {
 public:
  StoreFile(FILE* fileptr);

  virtual
  size_t write(const void*, size_t, size_t) override final;

  virtual
  size_t read(void*, size_t, size_t) override final;

  virtual
  void sread(const char*, int, void*, size_t, size_t, const char*, Error*) override final;

  virtual
  bigint tell() override final;

  virtual
  int seek(bigint) override final;

  virtual
  int truncate(bigint) override final;

  virtual
  bool is_file() override final;

  virtual
  int close() override final;

  virtual
  int error() override final;

  FILE* file;
};

}

#endif
