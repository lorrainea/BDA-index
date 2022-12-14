/**
 * @file    em_sparse_phi_src/utils.cpp
 * @section LICENCE
 *
 * This file is part of EM-SparsePhi v0.1.0
 * See: http://www.cs.helsinki.fi/group/pads/
 *
 * Copyright (C) 2016
 *   Juha Karkkainen <juha.karkkainen (at) cs.helsinki.fi>
 *   Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cerrno>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string>
#include <fstream>
#include <algorithm>

#include "utils.hpp"


namespace em_sparse_phi_private {
namespace utils {

long double wclock() {
  timeval tim;
  gettimeofday(&tim, NULL);
  return tim.tv_sec + (tim.tv_usec / 1000000.0L);
}

std::FILE *file_open(std::string filename, std::string mode) {
  std::FILE *f = std::fopen(filename.c_str(), mode.c_str());
  if (f == NULL) {
    std::perror(filename.c_str());
    std::exit(EXIT_FAILURE);
  }
  return f;
}

std::FILE *file_open_nobuf(std::string filename, std::string mode) {
  std::FILE *f = std::fopen(filename.c_str(), mode.c_str());
  if (f == NULL) {
    std::perror(filename.c_str());
    std::exit(EXIT_FAILURE);
  }
  if(std::setvbuf(f, NULL, _IONBF, 0) != 0) {
    perror("setvbuf failed");
    std::exit(EXIT_FAILURE);
  }
  return f;
}

std::uint64_t file_size(std::string filename) {
  std::FILE *f = file_open_nobuf(filename, "r");
  std::fseek(f, 0, SEEK_END);
  long size = std::ftell(f);
  if (size < 0) {
    std::perror(filename.c_str());
    std::exit(EXIT_FAILURE);
  }
  std::fclose(f);
  return (std::uint64_t)size;
}

bool file_exists(std::string filename) {
  std::FILE *f = std::fopen(filename.c_str(), "r");
  bool result = (f != NULL);
  if (f != NULL)
    std::fclose(f);
  return result;
}

void file_delete(std::string filename) {
  int res = std::remove(filename.c_str());
  if (res != 0) {
    std::perror(filename.c_str());
    std::exit(EXIT_FAILURE);
  }
}

std::string absolute_path(std::string filename) {
  char path[1 << 12];
  bool created = false;
  if (!file_exists(filename)) {
    std::fclose(file_open(filename, "w"));
    created = true;
  }
  if (!realpath(filename.c_str(), path)) {
    std::perror(filename.c_str());
    std::exit(EXIT_FAILURE);
  }
  if (created)
    file_delete(filename);
  return std::string(path);
}

void drop_disk_pages(std::string filename) {
  int fd = open(filename.c_str(), O_RDWR);
  if (fd == -1) {
    std::perror(filename.c_str());
    std::exit(EXIT_FAILURE);
  }
  off_t length = lseek(fd, 0, SEEK_END);
  lseek(fd, 0L, SEEK_SET);
  posix_fadvise(fd, 0, length, POSIX_FADV_DONTNEED);
  close(fd);
}

std::int32_t random_int32(std::int32_t p, std::int32_t r) {
  return p + rand() % (r - p + 1);
}

std::int64_t random_int64(std::int64_t p, std::int64_t r) {
  std::int64_t x = random_int32(0, 1000000000);
  std::int64_t y = random_int32(0, 1000000000);
  std::int64_t z = x * 1000000000L + y;
  return p + z % (r - p + 1);
}

void fill_random_string(std::uint8_t* &s, std::uint64_t length, std::uint64_t sigma) {
  for (std::uint64_t i = 0; i < length; ++i)
    s[i] = random_int32(0, sigma - 1);
}

void fill_random_letters(std::uint8_t* &s, std::uint64_t length, std::uint64_t sigma) {
  fill_random_string(s, length, sigma);
  for (std::uint64_t i = 0; i < length; ++i)
    s[i] += 'a';
}

std::string random_string_hash() {
  uint64_t hash = (uint64_t)rand() * RAND_MAX + rand();
  std::stringstream ss;
  ss << hash;
  return ss.str();
}

std::uint64_t log2ceil(std::uint64_t x) {
  std::uint64_t pow2 = 1, w = 0;
  while (pow2 < x) { pow2 <<= 1; ++w; }
  return w;
}

std::uint64_t log2floor(std::uint64_t x) {
  std::uint64_t pow2 = 1, w = 0;
  while ((pow2 << 1) <= x) { pow2 <<= 1; ++w; }
  return w;
}

}  // namespace utils
}  // namespace em_sparse_phi_private
