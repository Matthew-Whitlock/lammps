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

#include "store_sizer.h"
#include "utils.h"
#include "platform.h"

namespace LAMMPS_NS {

StoreSizer::StoreSizer(LAMMPS* lmp, bigint& output)
  : Pointers(lmp), max_pos(output)
{
  cur_pos = max_pos = 0;
};

size_t StoreSizer::write(const void*, size_t size, size_t count){
  cur_pos += size*count;
  max_pos = std::max(max_pos, cur_pos);
  return size == 0 ? 0 : count;
}

size_t StoreSizer::read(void*, size_t, size_t){
  lmp->error->one(FLERR, "Attempting to read from StoreSizer");
  return 0;
}

void StoreSizer::sread(
  const char* srcname, int srcline, void*, size_t, size_t,
  const char* filename, Error* err
) {
  err->one(
    srcname, srcline, "Attempting to read from StoreSizer named {}", filename
  );
}

bigint StoreSizer::tell(){
  return cur_pos;
}

int StoreSizer::seek(bigint pos){
  if(pos == platform::END_OF_FILE) pos = max_pos;
  cur_pos = pos;
  max_pos = std::max(max_pos, cur_pos);
  return 0;
}

int StoreSizer::truncate(bigint length){
  cur_pos = length;
  max_pos = std::max(max_pos, cur_pos);
  return 0;
}

}
