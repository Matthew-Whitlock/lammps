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

#ifndef LMP_STORE_H
#define LMP_STORE_H

#include <memory>
#include <cstdio>
#include <vector>

#include "store_impl.h"
#include "error.h"
#include "safe_pointers.h"

namespace LAMMPS_NS {

/* A Store works with the FILE* API, but allows different backends
 * to store data in other structures (e.g. directly in memory).
 */
class Store {
 public:
  Store() = default;
  Store(FILE*);
  Store(const SafeFilePtr&);
  Store(std::shared_ptr<StoreImpl> impl);

  static
  Store Sizer(LAMMPS*, bigint&);
  static
  Store Buffer(LAMMPS*, std::vector<char>&);
  static
  Store Buffer(LAMMPS*, char*, bigint);

  Store& operator=(FILE*);
  bool operator==(FILE* const&);
  bool operator!=(FILE* const&);
  bool operator==(std::nullptr_t const&);
  bool operator!=(std::nullptr_t const&);
  bool operator==(Store const&);
  bool operator!=(Store const&);

  explicit operator bool() const;

  std::shared_ptr<StoreImpl> impl;
};


//Overloads for the FILE* functions

namespace utils {
/* See utils.h */
void sfread(const char*, int, void*, size_t, size_t, Store, const char*, Error*);
}

namespace platform {
/* See platform.h */
bigint ftell(Store);
int fseek(Store, bigint);
int ftruncate(Store, bigint);
}

size_t fwrite(const void*, size_t, size_t, Store);
size_t fread(void*, size_t, size_t, Store);
int fclose(Store);
int ferror(Store);

}

#endif
