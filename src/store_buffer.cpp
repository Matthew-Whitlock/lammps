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

#include "store_buffer.h"
#include "utils.h"
#include "platform.h"

#include <vector>
#include <cstring>
#include <algorithm>

namespace LAMMPS_NS {

StoreBuffer::StoreBuffer(LAMMPS* lmp, std::vector<char>& buf_vec)
  : Pointers(lmp), buf(buf_vec.data()), buf_len(buf_vec.size())
{
  cur_pos = cur_len = 0;
  eflag = false;
};

StoreBuffer::StoreBuffer(LAMMPS* lmp, char* buf_ptr, bigint len)
  : Pointers(lmp), buf(buf_ptr), buf_len(len)
{
  cur_pos = cur_len = 0;
  eflag = false;
}

size_t StoreBuffer::write(const void* data, size_t size, size_t count){
  if(size == 0 || count == 0) return 0;

  bigint avail = buf_len - cur_pos;
  bigint store_count = std::min(count, avail / size);
  if(store_count != count) eflag = true;

  bigint store_bytes = store_count * size;
  memcpy(buf+cur_len, data, store_bytes);

  cur_pos += store_bytes;
  cur_len = std::max(cur_len, cur_pos);

  if(eflag) lmp->error->one(FLERR, "Could not write {} items of size {}, could only write {}\n", count, size, store_count);
  return store_count;
}

size_t StoreBuffer::read(void* out, size_t size, size_t count){
  if(size == 0 || count == 0) return 0;

  bigint avail = cur_len - cur_pos;
  bigint read_count = std::min(count, avail / size);
  if(read_count != count) eflag = true;

  bigint read_bytes = read_count * size;
  memcpy(out, buf+cur_pos, read_bytes);

  cur_pos += read_bytes;
  
  if(eflag) lmp->error->one(FLERR, "Could not read {} items of size {}, could only rad {}\n", count, size, read_count);
  return read_count;
}

void StoreBuffer::sread(
  const char* srcname, int srcline, void* out, size_t size, size_t count,
  const char* filename, Error* err
) {
  size_t read_count = read(out, size, count);
  if(read_count != count){
    err->one(
      srcname, srcline,
      "Attempting to read {} items of size {} from StoreBuffer named {}, but "
      "only {} could be read", count, size, filename, read_count
    );
  }
}

bigint StoreBuffer::tell(){
  return cur_pos;
}

int StoreBuffer::seek(bigint pos){
  if(pos == platform::END_OF_FILE) pos = cur_len;
  if(pos > buf_len){
    if(eflag) lmp->error->one(FLERR, "Could not seek to {} when buf_len is {}\n", pos, buf_len);
    return -1;
  }
  if(pos > cur_len) cur_len = pos;
  cur_pos = pos;
  return 0;
}

int StoreBuffer::truncate(bigint length){
  if(length > buf_len){
    eflag = true;
    return -1;
  }
  cur_pos = cur_len = length;
  return 0;
}

int StoreBuffer::error() {
  return eflag;
}

}
