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

#ifndef LMP_STORE_IMPL_H
#define LMP_STORE_IMPL_H

#include "error.h"

#include <cstdio>

namespace LAMMPS_NS {

class StoreImpl {
 public:
  /* As fwrite */
  virtual
  size_t write(const void*, size_t, size_t) = 0;

  /* As fread */
  virtual
  size_t read(void* buf, size_t, size_t) = 0;

  /* As utils::sfread */
  virtual
  void sread(const char*, int, void*, size_t, size_t, const char*, Error*) = 0;

  /* As platform::ftell */
  virtual
  bigint tell() = 0;

  /* As platform::fseek */
  virtual
  int seek(bigint) = 0;

  /* As platform::ftruncate */
  virtual
  int truncate(bigint) = 0;

  /* If this StoreImpl actually has a file pointer for IO */
  virtual
  bool is_file() { return false; }

  /* As fclose */
  virtual
  int close() { return 0; }

  /* As ferror */
  virtual
  int error() { return 0; }
};

}

#endif
