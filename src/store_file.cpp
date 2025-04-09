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

#include "store_file.h"
#include "utils.h"
#include "platform.h"

namespace LAMMPS_NS {

StoreFile::StoreFile(FILE* fileptr) : file(fileptr) {};

size_t StoreFile::write(const void* buf, size_t size, size_t count){
  return fwrite(buf, size, count, file);
}

size_t StoreFile::read(void* buf, size_t size, size_t count){
  return fread(buf, size, count, file);
}

void StoreFile::sread(const char* srcname, int srcline, void* s, size_t size, size_t num,
                      const char* filename, Error* error){
  utils::sfread(srcname, srcline, s, size, num, file, filename, error);
}

bigint StoreFile::tell(){
  return platform::ftell(file);
}

int StoreFile::seek(bigint pos){
  return platform::fseek(file, pos);
}

int StoreFile::truncate(bigint length){
  return platform::ftruncate(file, length);
}

bool StoreFile::is_file(){
  return true;
}

int StoreFile::close(){
  return fclose(file); 
}

int StoreFile::error(){
  return ferror(file);
}

}
