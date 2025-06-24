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

#include "store.h"
#include "store_file.h"
#include "store_sizer.h"
#include "store_buffer.h"

namespace LAMMPS_NS {

Store::Store(FILE* file) : impl(std::make_shared<StoreFile>(file)) {};
Store::Store(const SafeFilePtr& file) : Store(static_cast<FILE*>(file)) {};
Store::Store(std::shared_ptr<StoreImpl> impl_in) : impl(impl_in) {};

Store Store::Sizer(LAMMPS* lmp, bigint& output){
  return Store(std::make_shared<StoreSizer>(lmp, output));
}

Store Store::Buffer(LAMMPS* lmp, std::vector<char>& buf){
  return Store(std::make_shared<StoreBuffer>(lmp, buf));
}
Store Store::Buffer(LAMMPS* lmp, char* buf, bigint buf_len){
  return Store(std::make_shared<StoreBuffer>(lmp, buf, buf_len));
}

Store& Store::operator=(FILE* file){
  if(file == nullptr) impl = nullptr;
  else impl = std::make_shared<StoreFile>(file);
  return *this;
}

bool Store::operator==(FILE* const& file){
  if(impl == nullptr) return file == nullptr;
  return impl->is_file() &&
    file == std::static_pointer_cast<StoreFile>(impl)->file;
}
bool Store::operator!=(FILE* const& file){
  return !(*this == file);
}

bool Store::operator==(std::nullptr_t const& ptr){
  return impl == nullptr || (impl->is_file() &&
    std::static_pointer_cast<StoreFile>(impl)->file == nullptr);
}
bool Store::operator!=(std::nullptr_t const& ptr){
  return !(*this == ptr);
}

bool Store::operator==(Store const& other){
  return impl == other.impl;
}
bool Store::operator!=(Store const& other){
  return impl != other.impl;
}

Store::operator bool() const {
  return impl != nullptr;
}


namespace utils {
void sfread(const char* srcname, int srcline, void* s, size_t size, size_t num,
            Store store, const char* filename, Error* error){
  store.impl->sread(srcname, srcline, s, size, num, filename, error);
}
}

namespace platform {
bigint ftell(Store store){
  return store.impl->tell();
}

int fseek(Store store, bigint pos){
  return store.impl->seek(pos);
}

int ftruncate(Store store, bigint length){
  return store.impl->truncate(length);
}
}

size_t fwrite(const void* buf, size_t size, size_t count, Store store){
  return store.impl->write(buf, size, count);
}

size_t fread(void* buf, size_t size, size_t count, Store store){
  return store.impl->read(buf, size, count);
}

int fclose(Store store){
  return store.impl->close();
}

int ferror(Store store){
  return store.impl->error();
}

}
