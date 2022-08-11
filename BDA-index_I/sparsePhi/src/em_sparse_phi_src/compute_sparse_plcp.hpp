/**
 * @file    em_sparse_phi_src/compute_sparse_plcp.hpp
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

#ifndef __EM_SPARSE_PHI_SRC_COMPUTE_SPARSE_PHI_HPP_INCLUDED
#define __EM_SPARSE_PHI_SRC_COMPUTE_SPARSE_PHI_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <cinttypes>
#include <algorithm>
#include <parallel/algorithm>
#include <omp.h>

#include "io/async_stream_reader.hpp"
#include "io/async_stream_writer.hpp"
#include "io/async_multi_stream_writer.hpp"
#include "utils.hpp"


namespace em_sparse_phi_private {

class text_accessor {
  public:
    text_accessor(std::string filename,
        std::uint64_t buf_size = (1UL << 20)) {
      m_buf_size = buf_size;
      m_file_size = utils::file_size(filename);
      m_file = utils::file_open_nobuf(filename, "r");
      m_buf = new std::uint8_t[m_buf_size];
      m_buf_pos = 0;
      m_buf_filled = 0;
      m_bytes_read = 0;
    }

    inline std::uint8_t access(std::uint64_t i) {
      if (!(m_buf_pos <= i && i < m_buf_pos + m_buf_filled)) {
        if (m_buf_pos + m_buf_filled != i)
          std::fseek(m_file, i, SEEK_SET);
        m_buf_pos = i;
        m_buf_filled = std::min(m_buf_size, m_file_size - m_buf_pos);
        utils::read_from_file(m_buf, m_buf_filled, m_file);
        m_bytes_read += m_buf_filled;
      }
      return m_buf[i - m_buf_pos];
    }

    inline std::uint64_t bytes_read() const {
      return m_bytes_read;
    }

    ~text_accessor() {
      std::fclose(m_file);
      delete[] m_buf;
    }

  private:
    std::uint64_t m_bytes_read;
    std::uint64_t m_file_size;
    std::uint64_t m_buf_size;
    std::uint64_t m_buf_pos;
    std::uint64_t m_buf_filled;
    std::uint8_t *m_buf;
    std::FILE *m_file;
};

template<typename T>
struct key_val_pair {
  key_val_pair() {}
  key_val_pair(T key, T val)
    : m_key(key), m_val(val) {}

  inline bool operator < (const key_val_pair<T> &x) const {
    return m_key < x.m_key;
  }

  T m_key;
  T m_val;
};

// Note: the only assumption about saidx_t we have is that it is
// capable of storing values in the range [0..text_length). text_length
// cannot be assigned to variable of type saidx_t!

// Returns the position, for which Phi is undefined
template<typename saidx_t>
std::uint64_t compute_sparse_phi(std::string sa_filename, std::uint64_t text_length,
    std::uint64_t ram_use, std::uint64_t max_halfsegment_size,
    std::uint64_t **halfsegment_pair_count, saidx_t *sparse_phi,
    std::uint64_t plcp_sampling_rate, std::uint64_t &total_io_volume) {
  fprintf(stderr, "  Compute sparse Phi: ");
  long double start = utils::wclock();

  std::uint64_t phi_undefined_position = 0;
  std::uint64_t n_halfsegments = (text_length + max_halfsegment_size - 1) / max_halfsegment_size;
  std::uint64_t sparse_phi_size = (text_length + plcp_sampling_rate - 1) / plcp_sampling_rate;
  std::uint64_t sparse_phi_ram = sparse_phi_size * sizeof(saidx_t);

#ifdef _OPENMP
  // Allocate space for auxiliary data structures (buffers, etc.).
  static const std::uint64_t opt_in_buf_ram = (32UL << 20);
  static const std::uint64_t opt_par_buf_ram = (4UL << 20);

  std::uint64_t in_buf_ram = opt_in_buf_ram;
  std::uint64_t par_buf_ram = opt_par_buf_ram;

  {
    std::uint64_t total_buf_ram = in_buf_ram + par_buf_ram;

    if (sparse_phi_ram + total_buf_ram > ram_use) {
      std::uint64_t ram_budget = ram_use - sparse_phi_ram;
      long double shrink_factor = (long double)ram_budget / (long double)total_buf_ram;
      in_buf_ram = (std::uint64_t)((long double)in_buf_ram * shrink_factor);
      par_buf_ram = (std::uint64_t)((long double)par_buf_ram * shrink_factor);
    }
  }

  typedef async_stream_reader<saidx_t> sa_reader_type;
  sa_reader_type *sa_reader = new sa_reader_type(sa_filename, in_buf_ram, std::max(4UL, in_buf_ram / (2UL << 20)));

  std::uint64_t par_buf_size = par_buf_ram / sizeof(saidx_t);
  saidx_t *sa_buffer = new saidx_t[par_buf_size];

  std::uint64_t sa_items_read = 0;
  std::uint64_t prev_sa = 0;
  std::uint64_t step_count = 0;

  while (sa_items_read < text_length) {
    ++step_count;
    if (step_count == (1UL << 28) / par_buf_size) {
      step_count = 0;
      std::uint64_t io_vol = sa_items_read * sizeof(saidx_t);
      long double elapsed = utils::wclock() - start;
      fprintf(stderr, "\r  Compute sparse Phi: %.1Lf%%, time = %.2Lfs, I/O = %.2LfMiB/s",
          (100.L * sa_items_read) / text_length, elapsed, (1.L * io_vol / (1L << 20)) / elapsed);
    }

    std::uint64_t sa_buffer_filled = std::min(text_length - sa_items_read, par_buf_size);
    sa_reader->read(sa_buffer, sa_buffer_filled);

    // Process first element separatelly.
    if (sa_items_read == 0) {
      phi_undefined_position = sa_buffer[0];
    } else {
      std::uint64_t sa_val = sa_buffer[0];
      if ((sa_val % plcp_sampling_rate) == 0)
        sparse_phi[sa_val / plcp_sampling_rate] = prev_sa;
      if (halfsegment_pair_count != NULL)
        ++halfsegment_pair_count[sa_val / max_halfsegment_size][prev_sa / max_halfsegment_size];
    }

    // Process remaining elements in the buffer.
    std::uint64_t max_threads = omp_get_max_threads();
    std::uint64_t max_range_size = (sa_buffer_filled + max_threads - 1) / max_threads;
    std::uint64_t n_ranges = (sa_buffer_filled + max_range_size - 1) / max_range_size;
    #pragma omp parallel num_threads(n_ranges)
    {
      std::uint64_t thread_id = omp_get_thread_num();
      std::uint64_t range_beg = thread_id * max_range_size;
      std::uint64_t range_end = std::min(range_beg + max_range_size, sa_buffer_filled);
      std::uint64_t **local_halfsegment_pair_count = NULL;
      if (halfsegment_pair_count != NULL) {
        local_halfsegment_pair_count = new std::uint64_t*[n_halfsegments];
        for (std::uint64_t j = 0; j < n_halfsegments; ++j) {
          local_halfsegment_pair_count[j] = new std::uint64_t[n_halfsegments];
          std::fill(local_halfsegment_pair_count[j], local_halfsegment_pair_count[j] + n_halfsegments, 0UL);
        }
      }

      std::uint64_t prev_segment_id = ((std::uint64_t)sa_buffer[std::max(1UL, range_beg) - 1]) / max_halfsegment_size;
      for (std::uint64_t i = std::max(1UL, range_beg); i < range_end; ++i) {
        std::uint64_t sa_val = sa_buffer[i];
        if ((sa_val % plcp_sampling_rate) == 0)
          sparse_phi[sa_val / plcp_sampling_rate] = sa_buffer[i - 1];

        if (halfsegment_pair_count != NULL) {
          std::uint64_t segment_id = sa_val / max_halfsegment_size;
          ++local_halfsegment_pair_count[segment_id][prev_segment_id];
          prev_segment_id = segment_id;
        }
      }

      if (halfsegment_pair_count != NULL) {
        #pragma omp critical
        {
          for (std::uint64_t i = 0; i < n_halfsegments; ++i)
            for (std::uint64_t j = 0; j < n_halfsegments; ++j)
              halfsegment_pair_count[i][j] += local_halfsegment_pair_count[i][j];
        }
      }

      if (halfsegment_pair_count != NULL) {
        for (std::uint64_t j = 0; j < n_halfsegments; ++j)
          delete[] local_halfsegment_pair_count[j];
        delete[] local_halfsegment_pair_count;
      }
    }

    // Update the number of read SA items.
    prev_sa = sa_buffer[sa_buffer_filled - 1];
    sa_items_read += sa_buffer_filled;
  }

  delete[] sa_buffer;
#else
  // Allocate space for auxiliary data structures (buffers, etc.).
  static const std::uint64_t opt_in_buf_ram = (32UL << 20);
  std::uint64_t in_buf_ram = opt_in_buf_ram;
  if (sparse_phi_ram + in_buf_ram > ram_use)
    in_buf_ram = ram_use - sparse_phi_ram;

  typedef async_stream_reader<saidx_t> sa_reader_type;
  sa_reader_type *sa_reader = new sa_reader_type(sa_filename, in_buf_ram, std::max(4UL, in_buf_ram / (2UL << 20)));

  saidx_t prev_sa = 0;
  std::uint64_t prev_sa_halfsegment_id = 0;
  std::uint64_t step_count = 0;

  for (std::uint64_t i = 0; i < text_length; ++i) {
    ++step_count;
    if (step_count == (8UL << 20)) {
      step_count = 0;
      std::uint64_t io_vol = i * sizeof(saidx_t);
      long double elapsed = utils::wclock() - start;
      fprintf(stderr, "\r  Compute sparse Phi: %.1Lf%%, time = %.2Lfs, I/O = %.2LfMiB/s",
          (100.L * i) / text_length, elapsed, (1.L * io_vol / (1L << 20)) / elapsed);
    }

    saidx_t sa_i = sa_reader->read();
    std::uint64_t sa_i_uint64 = sa_i;

    if (i == 0) {
      phi_undefined_position = sa_i_uint64;
      prev_sa = sa_i;
      prev_sa_halfsegment_id = sa_i_uint64 / max_halfsegment_size;
      continue;
    }

    if ((sa_i_uint64 % plcp_sampling_rate) == 0)
      sparse_phi[sa_i_uint64 / plcp_sampling_rate] = prev_sa;

    if (halfsegment_pair_count != NULL) {
      std::uint64_t halfsegment_id = sa_i_uint64 / max_halfsegment_size;
      ++halfsegment_pair_count[halfsegment_id][prev_sa_halfsegment_id];
      prev_sa_halfsegment_id = halfsegment_id;
    }

    prev_sa = sa_i;
  }
#endif

  // halfsegment_pair_count == NULL means than we
  // do now employ the text-partitioning.
  if (halfsegment_pair_count != NULL) {
    for (std::uint64_t i = 0; i < n_halfsegments; ++i) {
      for (std::uint64_t j = 0; j < n_halfsegments; ++j) {
        if (i > j) {
          halfsegment_pair_count[j][i] += halfsegment_pair_count[i][j];
          halfsegment_pair_count[i][j] = 0;
        }
      }
    }
  }

  long double total_time = utils::wclock() - start;
  std::uint64_t io_vol = sa_reader->bytes_read();
  total_io_volume += io_vol;
  fprintf(stderr, "\r  Compute sparse Phi: 100.0%%, time = %.2Lfs, I/O = %.2LfMiB/s, total I/O vol = %.2Lfn\n",
      total_time, (1.L * io_vol / (1L << 20)) / total_time, (1.L * total_io_volume) / text_length);

  delete sa_reader;
  return phi_undefined_position;
}

// Returns the sparse PLCP array.
template<typename saidx_t>
saidx_t* compute_sparse_plcp(std::string text_filename, std::string sa_filename, std::string output_filename,
    std::uint64_t text_length, std::uint64_t max_halfsegment_size, std::uint64_t **halfsegment_pair_count,
    std::uint64_t plcp_sampling_rate, std::uint64_t ram_use, std::uint64_t &total_io_volume) {
  std::uint64_t sparse_plcp_size = (text_length + plcp_sampling_rate - 1) / plcp_sampling_rate;
  std::uint64_t sparse_phi_size = sparse_plcp_size;
  std::uint64_t sparse_plcp_ram = sparse_plcp_size * sizeof(saidx_t);
  std::uint64_t sparse_phi_ram = sparse_plcp_ram;

  fprintf(stderr, "Compute sparse PLCP:\n");
  long double start = utils::wclock();

  // STEP 1: Compute sparse Phi from SA.
  saidx_t *sparse_phi = (saidx_t *)malloc(sparse_phi_ram);
  std::uint64_t phi_undefined_position = compute_sparse_phi(sa_filename, text_length, ram_use,
      max_halfsegment_size, halfsegment_pair_count, sparse_phi, plcp_sampling_rate, total_io_volume);

  std::uint64_t max_segment_size = 0;
  std::uint64_t n_segments = 0;

  // STEP 2: Compute sparse PLCP and write into disk.
  {
#ifdef _OPENMP
    std::uint64_t aux_sparse_plcp_sampling_rate = plcp_sampling_rate * std::max(10000UL, ((64UL << 20) / plcp_sampling_rate));
    std::uint64_t aux_sparse_plcp_size = (text_length + aux_sparse_plcp_sampling_rate - 1) / aux_sparse_plcp_sampling_rate;
    std::uint64_t aux_sparse_plcp_ram = aux_sparse_plcp_size * sizeof(std::uint64_t);
    std::uint64_t *aux_sparse_plcp = (std::uint64_t *)malloc(aux_sparse_plcp_ram);

    // Compute aux sparse PLCP.
    {
      fprintf(stderr, "  Compute mini sparse PLCP: ");
      long double mini_plcp_start = utils::wclock();

      std::uint64_t *aux_sparse_phi = new std::uint64_t[aux_sparse_plcp_size];
      for (std::uint64_t j = 0; j < aux_sparse_plcp_size; ++j)
        aux_sparse_phi[j] = sparse_phi[j * (aux_sparse_plcp_sampling_rate / plcp_sampling_rate)];

      static const std::uint64_t opt_txt_acc_buf_ram = (1UL << 20);
      std::uint64_t txt_acc_buf_ram = opt_txt_acc_buf_ram;

      {
        std::uint64_t total_txt_acc_buf_ram = 2 * txt_acc_buf_ram;
        if (sparse_phi_ram + aux_sparse_plcp_ram + total_txt_acc_buf_ram > ram_use) {
          total_txt_acc_buf_ram = ram_use - sparse_phi_ram - aux_sparse_plcp_ram;
          txt_acc_buf_ram = total_txt_acc_buf_ram / 2;
        }
      }

      text_accessor *a1 = new text_accessor(text_filename, txt_acc_buf_ram);
      text_accessor *a2 = new text_accessor(text_filename, txt_acc_buf_ram);
      std::uint64_t lcp = 0;
      std::uint64_t prev_i = 0;
      std::uint64_t prev_lcp = 0;
      for (std::uint64_t j = 0; j < aux_sparse_plcp_size; ++j) {
        std::uint64_t i = j * aux_sparse_plcp_sampling_rate;
        std::uint64_t phi_i = aux_sparse_phi[j];
        if (i == phi_undefined_position) lcp = 0;
        else {
          lcp = (std::uint64_t)std::max(0L, (std::int64_t)(prev_i + prev_lcp) - (std::int64_t)i);
          while (i + lcp < text_length && phi_i + lcp < text_length && a1->access(i + lcp) == a2->access(phi_i + lcp))
            ++lcp;
        }
        aux_sparse_plcp[j] = lcp;
        prev_i = i;
        prev_lcp = lcp;
      }

      long double mini_plcp_time = utils::wclock() - mini_plcp_start;
      std::uint64_t io_vol = a1->bytes_read() + a2->bytes_read();
      total_io_volume += io_vol;
      fprintf(stderr, "%.2Lfs, I/O = %.2LfMiB/s, total I/O vol = %.2Lfn\n", mini_plcp_time,
          (1.L * io_vol / (1L << 20)) / mini_plcp_time, (1.L * total_io_volume) / text_length);
      delete a1;
      delete a2;
      delete[] aux_sparse_phi;
    }
#endif

    max_segment_size = (std::uint64_t)((long double)ram_use * 0.9L);  // minimal acceptable values of max_segment_size

    // Shrink buffers or enlarge segment to fill in the available RAM.
#ifndef _OPENMP
    static const std::uint64_t opt_in_buf_ram = (32UL << 20);
    static const std::uint64_t opt_out_buf_ram = (4UL << 20);
    static const std::uint64_t opt_txt_acc_buf_ram = (1UL << 20);
    static const std::uint64_t opt_overflow_size = (1UL << 20);

    std::uint64_t in_buf_ram = opt_in_buf_ram;
    std::uint64_t out_buf_ram = opt_out_buf_ram;
    std::uint64_t txt_acc_buf_ram = opt_txt_acc_buf_ram;
    std::uint64_t overflow_size = opt_overflow_size;

    {
      std::uint64_t total_buf_size = in_buf_ram + out_buf_ram + txt_acc_buf_ram + overflow_size;

      if (max_segment_size + total_buf_size > ram_use) {
        std::uint64_t ram_budget = ram_use - max_segment_size;
        long double shrink_factor = (long double)ram_budget / (long double)total_buf_size;
        in_buf_ram = (std::uint64_t)((long double)in_buf_ram * shrink_factor);
        out_buf_ram = (std::uint64_t)((long double)out_buf_ram * shrink_factor);
        txt_acc_buf_ram = (std::uint64_t)((long double)txt_acc_buf_ram * shrink_factor);
        overflow_size = (std::uint64_t)((long double)overflow_size * shrink_factor);
      } else max_segment_size = ram_use - total_buf_size;
    }
#else
    static const std::uint64_t opt_in_buf_ram = (2UL << 20);
    static const std::uint64_t opt_out_buf_ram = (4UL << 20);
    static const std::uint64_t opt_txt_acc_buf_ram = (1UL << 20);
    static const std::uint64_t opt_local_buf_ram = (1UL << 20);
    static const std::uint64_t opt_overflow_size = (1UL << 20);

    std::uint64_t in_buf_ram = opt_in_buf_ram;
    std::uint64_t out_buf_ram = opt_out_buf_ram;
    std::uint64_t txt_acc_buf_ram = opt_txt_acc_buf_ram;
    std::uint64_t local_buf_ram = opt_local_buf_ram;
    std::uint64_t overflow_size = opt_overflow_size;

    {
      std::uint64_t max_threads = omp_get_max_threads();
      std::uint64_t total_buf_size = out_buf_ram + overflow_size + max_threads * (in_buf_ram + local_buf_ram + txt_acc_buf_ram);

      if (max_segment_size + aux_sparse_plcp_ram + total_buf_size > ram_use) {
        std::uint64_t ram_budget = ram_use - aux_sparse_plcp_ram - max_segment_size;
        long double shrink_factor = (long double)ram_budget / (long double)total_buf_size;
        in_buf_ram = (std::uint64_t)((long double)in_buf_ram * shrink_factor);
        out_buf_ram = (std::uint64_t)((long double)out_buf_ram * shrink_factor);
        txt_acc_buf_ram = (std::uint64_t)((long double)txt_acc_buf_ram * shrink_factor);
        local_buf_ram = (std::uint64_t)((long double)local_buf_ram * shrink_factor);
        overflow_size = (std::uint64_t)((long double)overflow_size * shrink_factor);
      } else max_segment_size = ram_use - aux_sparse_plcp_ram - total_buf_size;
    }
#endif

    n_segments = (text_length + max_segment_size - 1) / max_segment_size;

    fprintf(stderr, "  Max segment size = %lu (%.2LfMiB)\n", max_segment_size, (1.L * max_segment_size / (1UL << 20)));
    fprintf(stderr, "  Number of segments = %lu\n", n_segments);

    std::uint64_t max_possible_sparse_phi_for_segment_size = (max_segment_size + plcp_sampling_rate - 1) / plcp_sampling_rate;
    std::uint64_t max_possible_sparse_phi_for_segment_ram = 2UL * sizeof(saidx_t) * max_possible_sparse_phi_for_segment_size;
    std::uint64_t extra_ram_needed = 0;
#ifdef _OPENMP
    extra_ram_needed += aux_sparse_plcp_ram;
#endif

    // Two cases, depending on whether we can sort sparse_phi_for_segment in RAM.
    if (max_possible_sparse_phi_for_segment_ram + extra_ram_needed > ram_use) {
      // Write every pair (i, Phi[i]) into file
      // corresponding to segment containing Phi[i].
      {
        fprintf(stderr, "  Distribute Phi into buckets: ");
        long double phi_distr_start = utils::wclock();

        static const std::uint64_t free_distr_bufs = 4;

        // Compute ram_for_distr_buf.
#ifdef _OPENMP
        static const std::uint64_t opt_distr_buf_ram = (4UL << 20);
        std::uint64_t distr_buf_ram = opt_distr_buf_ram;
        {
          std::uint64_t total_distr_bufs_ram = (n_segments + free_distr_bufs) * distr_buf_ram;
          if (total_distr_bufs_ram + sparse_phi_ram + aux_sparse_plcp_ram > ram_use) {
            total_distr_bufs_ram = ram_use - sparse_phi_ram - aux_sparse_plcp_ram;
            distr_buf_ram = total_distr_bufs_ram / (n_segments + free_distr_bufs);
          }
        }
#else
        static const std::uint64_t opt_distr_buf_ram = (4UL << 20);
        std::uint64_t distr_buf_ram = opt_distr_buf_ram;

        {
          std::uint64_t total_distr_bufs_ram = (n_segments + free_distr_bufs) * distr_buf_ram;
          if (total_distr_bufs_ram + sparse_phi_ram > ram_use) {
            total_distr_bufs_ram = ram_use - sparse_phi_ram;
            distr_buf_ram = total_distr_bufs_ram / (n_segments + free_distr_bufs);
          }
        }
#endif

        // Actual writing follows.
        typedef async_multi_stream_writer<saidx_t> phi_pair_multi_stream_writer_type;
        phi_pair_multi_stream_writer_type *phi_pair_multi_stream_writer = new phi_pair_multi_stream_writer_type(distr_buf_ram, free_distr_bufs);
        for (std::uint64_t segment_id = 0; segment_id < n_segments; ++segment_id)
          phi_pair_multi_stream_writer->add_file(output_filename + ".phi_pairs." + utils::intToStr(segment_id));

        std::uint64_t step_counter = 0;
        for (std::uint64_t i = 0, j = 0; j < sparse_phi_size; ++j, i += plcp_sampling_rate) {
          ++step_counter;
          if (step_counter == (16UL << 20)) {
            step_counter = 0;
            long double elapsed = utils::wclock() - phi_distr_start;
            std::uint64_t io_volume = 2UL * j * sizeof(saidx_t);
            fprintf(stderr, "\r  Distribute Phi into buckets: %.1Lf%%, time = %.1Lfs, I/O = %.2LfMiB/s",
                (100.L * j) / sparse_phi_size, elapsed, (1.L * io_volume / (1L << 20)) / elapsed);
          }

          if (i == phi_undefined_position)
            continue;

          saidx_t phi_i = sparse_phi[j];
          std::uint64_t segment_id = phi_i / max_segment_size;
          phi_pair_multi_stream_writer->write_to_ith_file(segment_id, i);
          phi_pair_multi_stream_writer->write_to_ith_file(segment_id, phi_i);
        }

        // Print summary.
        long double phi_distr_time = utils::wclock() - phi_distr_start;
        std::uint64_t io_vol = phi_pair_multi_stream_writer->bytes_written();
        total_io_volume += io_vol;
        fprintf(stderr, "\r  Distribute Phi into buckets: 100.0%%, time = %.1Lfs, I/O = %.2LfMiB/s, total I/O vol = %.2Lfn\n",
            phi_distr_time, ((1.L * io_vol) / (1L << 20)) / phi_distr_time, (1.L * total_io_volume) / text_length);

        // Clean up.
        delete phi_pair_multi_stream_writer;
      }

      free(sparse_phi);

      {
        // Load every segment and compute PLCP for pairs
        // (i, Phi[i]), where Phi is inside the segment.
        for (std::uint64_t segment_id = 0; segment_id < n_segments; ++segment_id) {
          std::uint64_t segment_beg = segment_id * max_segment_size;
          std::uint64_t segment_end = std::min(segment_beg + max_segment_size, text_length);
          std::uint64_t ext_segment_end = std::min(segment_end + overflow_size, text_length);
          std::uint64_t ext_segment_size = ext_segment_end - segment_beg;

          fprintf(stderr, "  Process segment %lu/%lu [%lu..%lu):\n", segment_id + 1, n_segments, segment_beg, segment_end);

          std::string phi_pairs_filename = output_filename + ".phi_pairs." + utils::intToStr(segment_id);
          std::string lcp_pairs_filename = output_filename + ".lcp_pairs." + utils::intToStr(segment_id);

          typedef async_stream_reader<saidx_t> phi_pair_reader_type;
          typedef async_stream_writer<saidx_t> lcp_pair_writer_type;
          typedef async_stream_reader<std::uint8_t> text_reader_type;

          // Allocate the segment.
          std::uint8_t *segment = (std::uint8_t *)malloc(ext_segment_size);

          // Read the segment from disk.
          fprintf(stderr, "    Read segment: ");
          long double read_start = utils::wclock();
          utils::read_at_offset(segment, segment_beg, ext_segment_size, text_filename);
          total_io_volume += ext_segment_size;
          long double read_time = utils::wclock() - read_start;
          fprintf(stderr, "%.2Lfs, I/O = %.2LfMiB/s, total I/O vol = %.2Lfn\n", read_time,
              ((1.L * ext_segment_size) / (1L << 20)) / read_time, (1.L * total_io_volume) / text_length);

          // Process the segment
          fprintf(stderr, "    Process (i, Phi[i]) pairs: ");
          long double process_start = utils::wclock();

          // Initialize readers and writers.
#ifndef _OPENMP
          lcp_pair_writer_type *lcp_pair_writer = new lcp_pair_writer_type(lcp_pairs_filename, out_buf_ram, std::max(4UL, out_buf_ram / (2UL << 20)));
          text_accessor *accessor = new text_accessor(text_filename, txt_acc_buf_ram);
          phi_pair_reader_type *phi_pair_reader = new phi_pair_reader_type(phi_pairs_filename, in_buf_ram / 2, std::max(4UL, in_buf_ram / (2UL << 20)));
          text_reader_type *text_reader = new text_reader_type(text_filename, in_buf_ram / 2, std::max(4UL, in_buf_ram / (2UL << 20)));

          if (phi_pair_reader->empty() == false) {
            std::uint64_t i = (std::uint64_t)phi_pair_reader->read();
            std::uint64_t phi_i = (std::uint64_t)phi_pair_reader->read();
            std::uint64_t lcp = 0;

            std::uint64_t buf_beg = 0;
            while (buf_beg < text_length) {
              text_reader->receive_new_buffer();
              std::uint64_t buf_filled = text_reader->get_buf_filled();
              const std::uint8_t *buffer = text_reader->get_buf_ptr();

              while (true) {
                while (phi_i + lcp < text_length && i + lcp < buf_beg + buf_filled) {
                  std::uint8_t next_char = (phi_i + lcp < ext_segment_end) ?
                    segment[(phi_i + lcp) - segment_beg] : accessor->access(phi_i + lcp);
                  if (next_char == buffer[(i + lcp) - buf_beg]) ++lcp;
                  else break;
                }

                if (i + lcp < buf_beg + buf_filled || phi_i + lcp == text_length || buf_beg + buf_filled == text_length) {
                  lcp_pair_writer->write((saidx_t)i);
                  lcp_pair_writer->write((saidx_t)lcp);
                  if (phi_pair_reader->empty() == false) {
                    std::uint64_t next_i = (std::uint64_t)phi_pair_reader->read();
                    std::uint64_t next_phi_i = (std::uint64_t)phi_pair_reader->read();
                    lcp = std::max(0L, (std::int64_t)(i + lcp) - (std::int64_t)next_i);
                    i = next_i;
                    phi_i = next_phi_i;
                  } else {
                    i = text_length;
                    break;
                  }
                } else break;
              }

              if (i == text_length)
                break;

              buf_beg += buf_filled;
            }
          }

          // Print summary.
          long double process_time = utils::wclock() - process_start;
          std::uint64_t io_vol = text_reader->bytes_read() + accessor->bytes_read() +
          lcp_pair_writer->bytes_written() + phi_pair_reader->bytes_read();
          total_io_volume += io_vol;
          fprintf(stderr, "\r    Process (i, Phi[i]) pairs: 100.0%%, time = %.2Lfs, I/O = %.2LfMiB/s, total I/O vol = %.2Lfn\n",
              process_time, (1.L * io_vol / (1L << 20)) / process_time, (1.L * total_io_volume) / text_length);

          delete lcp_pair_writer;
          delete accessor;
          delete phi_pair_reader;
          delete text_reader;
#else
          // Parallel computation of sparse PLCP values assuming each
          // thread handles the same amount of (i, Phi[i]) pairs (but
          // different threads can stream different amount of text).
          lcp_pair_writer_type *lcp_pair_writer = new lcp_pair_writer_type(lcp_pairs_filename, out_buf_ram, 2);
          std::uint64_t n_pairs = utils::file_size(phi_pairs_filename) / (2 * sizeof(saidx_t));
          std::uint64_t max_threads = omp_get_max_threads();
          std::uint64_t max_range_size = (n_pairs + max_threads - 1) / max_threads;
          std::uint64_t n_ranges = (n_pairs + max_range_size - 1) / max_range_size;
          std::uint64_t io_vol = 0;

          #pragma omp parallel num_threads(n_ranges)
          {
            std::uint64_t thread_id = omp_get_thread_num();
            std::uint64_t range_beg = thread_id * max_range_size;
            std::uint64_t range_end = std::min(n_pairs, range_beg + max_range_size);
            std::uint64_t range_size = range_end - range_beg;

            text_accessor *accessor = new text_accessor(text_filename, txt_acc_buf_ram);
            std::uint64_t local_buf_size = local_buf_ram / (2 * sizeof(saidx_t));
            saidx_t *local_buf = new saidx_t[2 * local_buf_size];
            phi_pair_reader_type *phi_pair_reader = new phi_pair_reader_type(phi_pairs_filename,
                in_buf_ram / 2, 2, 2 * range_beg * sizeof(saidx_t));

            std::uint64_t first_i = (std::uint64_t)phi_pair_reader->peek();
            std::uint64_t sample_plcp_addr = first_i / aux_sparse_plcp_sampling_rate;
            std::uint64_t sample_plcp_pos = sample_plcp_addr * aux_sparse_plcp_sampling_rate;
            std::uint64_t sample_plcp_val = aux_sparse_plcp[sample_plcp_addr];
            std::uint64_t first_i_plcp_lower_bound = (std::uint64_t)std::max(0L, (std::int64_t)sample_plcp_val -
                ((std::int64_t)first_i - (std::int64_t)sample_plcp_pos));
            std::uint64_t text_pos = first_i + first_i_plcp_lower_bound;
            text_reader_type *text_reader = new text_reader_type(text_filename, in_buf_ram / 2, 2, text_pos);

            std::uint64_t lcp = 0;
            std::uint64_t pairs_processed = 0;
            while (pairs_processed < range_size) {
              std::uint64_t local_buf_filled = std::min(local_buf_size, range_size - pairs_processed);
              for (std::uint64_t j = 0; j < local_buf_filled; ++j) {
                std::uint64_t i = (std::uint64_t)phi_pair_reader->read();
                std::uint64_t phi_i = (std::uint64_t)phi_pair_reader->read();
                lcp = (std::uint64_t)std::max(0L, (std::int64_t)text_pos - (std::int64_t)i);
                if (i + lcp != text_pos) {
                  text_reader->skip((i + lcp) - text_pos);
                  text_pos = i + lcp;
                }

                while (text_pos < text_length && phi_i + lcp < text_length) {
                  std::uint8_t next_char = (phi_i + lcp < ext_segment_end) ?
                    segment[(phi_i + lcp) - segment_beg] : accessor->access(phi_i + lcp);
                  if (next_char == text_reader->peek()) {
                    ++lcp;
                    ++text_pos;
                    text_reader->read();
                  } else break;
                }

                local_buf[2 * j] = (saidx_t)i;
                local_buf[2 * j + 1] = (saidx_t)lcp;
              }

              #pragma omp critical
              {
                lcp_pair_writer->write(local_buf, 2 * local_buf_filled);
              }
              pairs_processed += local_buf_filled;
            }

            #pragma omp critical
            {
              io_vol += accessor->bytes_read();
              io_vol += text_reader->bytes_read();
              io_vol += phi_pair_reader->bytes_read();
            }

            delete[] local_buf;
            delete accessor;
            delete text_reader;
            delete phi_pair_reader;
          }

          // Print summary.
          long double process_time = utils::wclock() - process_start;
          io_vol += lcp_pair_writer->bytes_written();
          total_io_volume += io_vol;
          fprintf(stderr, "\r    Process (i, Phi[i]) pairs: 100.0%%, time = %.2Lfs, I/O = %.2LfMiB/s, total I/O vol = %.2Lfn\n",
              process_time, (1.L * io_vol / (1L << 20)) / process_time, (1.L * total_io_volume) / text_length);

          delete lcp_pair_writer;
#endif
          utils::file_delete(phi_pairs_filename);
          free(segment);
        }
      }
    } else {
      // Processing that streams the text buffer-by-buffer.

      // Write every pair (i, Phi[i]) into file
      // corresponding to segment containing max(i, Phi[i]).
#ifdef _OPENMP
      // Bucket size used to partition (i, Phi[i]) pairs, so that
      // each thread processes roughly equal number of pairs.
      // Note: this array is assumed to be small, so we don't take it into
      // consideration when computing exact RAM usages. E.g., for 256GiB text
      // and 3.5GiB of RAM it takes only about 9MiB.
      static const std::uint64_t bucket_size_log = 24UL;
      static const std::uint64_t bucket_size = (1UL << bucket_size_log);
      std::vector<std::uint64_t> **type_1_bucket_sizes = new std::vector<std::uint64_t>*[n_segments];
      for (std::uint64_t segment_id = 0; segment_id < n_segments; ++segment_id) {
        std::uint64_t segment_beg = segment_id * max_segment_size;
        std::uint64_t segment_end = std::min(text_length, segment_beg + max_segment_size);
        std::uint64_t n_buckets = (segment_end + bucket_size - 1) / bucket_size;
        type_1_bucket_sizes[segment_id] = new std::vector<std::uint64_t>(n_buckets, 0UL);
      }
#endif

      typedef key_val_pair<saidx_t> pair_type;

      {
        fprintf(stderr, "  Distribute Phi into buckets: ");
        long double phi_distr_start = utils::wclock();

        // Initialize distr_buf_ram.
        static const std::uint64_t free_distr_bufs = 4;
        static const std::uint64_t opt_distr_buf_ram = (4UL << 20);
        std::uint64_t distr_buf_ram = opt_distr_buf_ram;

        // Shrink distr_buf_ram if necessary to fit in RAM budget.
        {
          std::uint64_t total_distr_bufs_ram = (2 * n_segments + free_distr_bufs) * distr_buf_ram;
#ifdef _OPEMP
          if (total_distr_bufs_ram + sparse_phi_ram + aux_sparse_plcp_ram > ram_use) {
            total_distr_bufs_ram = ram_use - sparse_phi_ram - aux_sparse_plcp_ram;
            distr_buf_ram = total_distr_bufs_ram / (2 * n_segments + free_distr_bufs);
          }
#else
          if (total_distr_bufs_ram + sparse_phi_ram > ram_use) {
            total_distr_bufs_ram = ram_use - sparse_phi_ram;
            distr_buf_ram = total_distr_bufs_ram / (2 * n_segments + free_distr_bufs);
          }
#endif
        }

        typedef async_multi_stream_writer<pair_type> phi_pair_multi_stream_writer_type;
        phi_pair_multi_stream_writer_type *phi_pair_multi_stream_writer = new phi_pair_multi_stream_writer_type(distr_buf_ram, free_distr_bufs);
        for (std::uint64_t segment_id = 0; segment_id < n_segments; ++segment_id) {
          phi_pair_multi_stream_writer->add_file(output_filename + ".phi_pairs." + utils::intToStr(segment_id) + ".1");
          phi_pair_multi_stream_writer->add_file(output_filename + ".phi_pairs." + utils::intToStr(segment_id) + ".2");
        }

        std::uint64_t step_counter = 0;
        for (std::uint64_t i = 0, j = 0; j < sparse_phi_size; ++j, i += plcp_sampling_rate) {
          ++step_counter;
          if (step_counter == (16UL << 20)) {
            step_counter = 0;
            long double elapsed = utils::wclock() - phi_distr_start;
            std::uint64_t io_volume = phi_pair_multi_stream_writer->bytes_written();
            fprintf(stderr, "\r  Distribute Phi into buckets: %.1Lf%%, time = %.1Lfs, I/O = %.2LfMiB/s",
                (100.L * j) / sparse_phi_size, elapsed, (1.L * io_volume / (1L << 20)) / elapsed);
          }

          if (i == phi_undefined_position)
            continue;

          saidx_t phi_i = sparse_phi[j];
          if (i < (std::uint64_t)phi_i) {
            std::uint64_t segment_id = (std::uint64_t)phi_i / max_segment_size;
            phi_pair_multi_stream_writer->write_to_ith_file(segment_id * 2, pair_type(i, phi_i));
#ifdef _OPENMP
            std::uint64_t bucket_id = (i >> bucket_size_log);
            (*type_1_bucket_sizes[segment_id])[bucket_id] += 1;
#endif
          } else {
            std::uint64_t segment_id = i / max_segment_size;
            phi_pair_multi_stream_writer->write_to_ith_file(segment_id * 2 + 1, pair_type(phi_i, i));
          }
        }

        // Print summary.
        long double phi_distr_time = utils::wclock() - phi_distr_start;
        std::uint64_t io_vol = phi_pair_multi_stream_writer->bytes_written();
        total_io_volume += io_vol;
        fprintf(stderr, "\r  Distribute Phi into buckets: 100.0%%, time = %.1Lfs, I/O = %.2LfMiB/s, total I/O vol = %.2Lfn\n",
            phi_distr_time, ((1.L * io_vol) / (1L << 20)) / phi_distr_time, (1.L * total_io_volume) / text_length);

        // Clean up.
        delete phi_pair_multi_stream_writer;
      }

      free(sparse_phi);

      {
        // Step 3, load every segment and compute PLCP for pairs
        // (i, Phi[i]), where Phi is inside the segment.
        for (std::uint64_t segment_id = 0; segment_id < n_segments; ++segment_id) {
          std::uint64_t segment_beg = segment_id * max_segment_size;
          std::uint64_t segment_end = std::min(segment_beg + max_segment_size, text_length);
          std::uint64_t ext_segment_end = std::min(segment_end + overflow_size, text_length);
          std::uint64_t ext_segment_size = ext_segment_end - segment_beg;

          fprintf(stderr, "  Process segment %lu/%lu [%lu..%lu):\n", segment_id + 1, n_segments, segment_beg, segment_end);

          std::string phi_pairs_filename_1 = output_filename + ".phi_pairs." + utils::intToStr(segment_id) + ".1";
          std::string phi_pairs_filename_2 = output_filename + ".phi_pairs." + utils::intToStr(segment_id) + ".2";
          std::string lcp_pairs_filename = output_filename + ".lcp_pairs." + utils::intToStr(segment_id);

          std::uint64_t n_pairs_type_1 = utils::file_size(phi_pairs_filename_1) / sizeof(pair_type);
          std::uint64_t n_pairs_type_2 = utils::file_size(phi_pairs_filename_2) / sizeof(pair_type);

          typedef async_stream_reader<pair_type> phi_pair_reader_type;
          typedef async_stream_reader<std::uint8_t> text_reader_type;
          typedef async_stream_writer<saidx_t> lcp_pair_writer_type;

#ifdef _OPENMP
          std::uint64_t n_buckets = (segment_end + bucket_size - 1) / bucket_size;
          std::vector<std::uint64_t> type_2_bucket_sizes(n_buckets, 0UL);
#endif

          // Read the file contaiing pairs (i, Phi[i]) such that Phi[i] < i
          // sort by Phi[i] and write to disk (overwriting the old file).
          {
            fprintf(stderr, "    Preprocess (i, Phi[i]) pairs: ");
            long double preprocess_start = utils::wclock();
            std::uint64_t sparse_phi_for_segment_size = n_pairs_type_2;
            std::uint64_t sparse_phi_for_segment_ram = sparse_phi_for_segment_size * sizeof(pair_type);
            pair_type *sparse_phi_for_segment = (pair_type *)malloc(sparse_phi_for_segment_ram);
            utils::read_from_file(sparse_phi_for_segment, sparse_phi_for_segment_size, phi_pairs_filename_2);
            total_io_volume += sparse_phi_for_segment_size * sizeof(sparse_phi_for_segment[0]);
            std::sort(sparse_phi_for_segment, sparse_phi_for_segment + sparse_phi_for_segment_size);

#ifdef _OPENMP
            for (std::uint64_t i = 0; i < sparse_phi_for_segment_size; ++i) {
              std::uint64_t bucket_id = ((std::uint64_t)(sparse_phi_for_segment[i].m_key)) >> bucket_size_log;
              type_2_bucket_sizes[bucket_id] += 1;
            }
#endif

            // Write sorted sparse_phi_for_segment back to disk (overwrite old file).
            utils::write_to_file(sparse_phi_for_segment, sparse_phi_for_segment_size, phi_pairs_filename_2);
            total_io_volume += sparse_phi_for_segment_size * sizeof(sparse_phi_for_segment[0]);

            // Clean up.
            free(sparse_phi_for_segment);
            long double preprocess_time = utils::wclock() - preprocess_start;
            fprintf(stderr, "%.2Lfs, total I/O vol = %.2Lfn\n", preprocess_time,
                (1.L * total_io_volume) / text_length);
          }

#ifdef _OPENMP
          std::uint64_t max_threads = omp_get_max_threads();
          std::uint64_t max_range_size = (segment_end + max_threads - 1) / max_threads;
          std::uint64_t n_ranges = (segment_end + max_range_size - 1) / max_range_size;
          std::vector<std::uint64_t> range_begs(n_ranges);
          std::vector<std::uint64_t> type_1_ptrs(n_ranges);
          std::vector<std::uint64_t> type_2_ptrs(n_ranges);

          std::uint64_t n_pairs = n_pairs_type_1 + n_pairs_type_2;
          std::uint64_t ideal_range_size = n_pairs / n_ranges;

          // Make the size of each bucket as close as possible to ideal_range_size.
          // Alternative to this would be to find the smallest upper bound on the
          // size of all buckets. It could be binary searched (to check a single
          // value during binary search, we greedily partition the buckets and then
          // check if the number of groups is <= n_ranges).
          std::uint64_t bucket_range_beg = 0;
          std::uint64_t type_1_buckets_total_size = 0;
          std::uint64_t type_2_buckets_total_size = 0;
          for (std::uint64_t range_id = 0; range_id < n_ranges; ++range_id) {
            range_begs[range_id] = bucket_range_beg * bucket_size;
            type_1_ptrs[range_id] = type_1_buckets_total_size;
            type_2_ptrs[range_id] = type_2_buckets_total_size;

            std::uint64_t bucket_range_end = bucket_range_beg;
            std::uint64_t cur_range_size = 0;

            // Add buckets to the current range as long as adding a new bucket
            // brings the range size closer to the 'ideal_range_size'.
            while (bucket_range_end < n_buckets && (range_id + 1 == n_ranges || std::abs((std::int64_t)ideal_range_size -
                    (std::int64_t)cur_range_size) >= std::abs((std::int64_t)ideal_range_size - (std::int64_t)(cur_range_size +
                        (*type_1_bucket_sizes[segment_id])[bucket_range_end] + type_2_bucket_sizes[bucket_range_end])))) {
              cur_range_size += (*type_1_bucket_sizes[segment_id])[bucket_range_end] + type_2_bucket_sizes[bucket_range_end];
              type_1_buckets_total_size += (*type_1_bucket_sizes[segment_id])[bucket_range_end];
              type_2_buckets_total_size += type_2_bucket_sizes[bucket_range_end];
              ++bucket_range_end;
            }

            bucket_range_beg = bucket_range_end;
          }
#endif

          // Allocate the segment.
          std::uint8_t *segment = (std::uint8_t *)malloc(ext_segment_size);

          // Read the segment from disk.
          fprintf(stderr, "    Read segment: ");
          long double read_start = utils::wclock();
          utils::read_at_offset(segment, segment_beg, ext_segment_size, text_filename);
          total_io_volume += ext_segment_size;
          long double read_time = utils::wclock() - read_start;
          fprintf(stderr, "%.2Lfs, I/O = %.2LfMiB/s, total I/O vol = %.2Lfn\n",
              read_time, ((1.L * ext_segment_size) / (1L << 20)) / read_time,
              (1.L * total_io_volume) / text_length);

          // Process the segment
          fprintf(stderr, "    Process (i, Phi[i]) pairs: ");
          long double process_start = utils::wclock();

          // Initialize readers and writers.
#ifndef _OPENMP
          text_reader_type *text_reader = new text_reader_type(text_filename, in_buf_ram / 2, std::max(4UL, (in_buf_ram / (4UL << 20))));
          lcp_pair_writer_type *lcp_pair_writer = new lcp_pair_writer_type(lcp_pairs_filename, out_buf_ram, std::max(4UL, out_buf_ram / (2UL << 20)));

          text_accessor *accessor_1 = new text_accessor(text_filename, txt_acc_buf_ram / 2);
          text_accessor *accessor_2 = new text_accessor(text_filename, txt_acc_buf_ram / 2);
          phi_pair_reader_type *phi_pair_reader_1 = new phi_pair_reader_type(phi_pairs_filename_1, in_buf_ram / 4, std::max(4UL, in_buf_ram / (4UL << 20)));
          phi_pair_reader_type *phi_pair_reader_2 = new phi_pair_reader_type(phi_pairs_filename_2, in_buf_ram / 4, std::max(4UL, in_buf_ram / (4UL << 20)));

          std::uint64_t i_1 = text_length;
          std::uint64_t phi_i_1 = 0;
          std::uint64_t lcp_1 = 0;
          if (n_pairs_type_1 > 0) {
            pair_type pp = phi_pair_reader_1->read();
            i_1 = (std::uint64_t)pp.m_key;
            phi_i_1 = (std::uint64_t)pp.m_val;
          }

          std::uint64_t i_2 = text_length;
          std::uint64_t phi_i_2 = 0;
          std::uint64_t lcp_2 = 0;
          if (n_pairs_type_2 > 0) {
            pair_type pp = phi_pair_reader_2->read();
            phi_i_2 = (std::uint64_t)pp.m_key;
            i_2 = (std::uint64_t)pp.m_val;
          }

          std::uint64_t buf_beg = 0;
          while (buf_beg < text_length && (i_1 != text_length || i_2 != text_length)) {
            text_reader->receive_new_buffer();
            std::uint64_t buf_filled = text_reader->get_buf_filled();
            const std::uint8_t *buffer = text_reader->get_buf_ptr();

            while (i_1 != text_length) {
              while (phi_i_1 + lcp_1 < text_length && i_1 + lcp_1 < buf_beg + buf_filled) {
                std::uint8_t next_char = (phi_i_1 + lcp_1 < ext_segment_end) ?
                  segment[(phi_i_1 + lcp_1) - segment_beg] : accessor_1->access(phi_i_1 + lcp_1);
                if (next_char == buffer[(i_1 + lcp_1) - buf_beg]) ++lcp_1;
                else break;
              }

              if (i_1 + lcp_1 < buf_beg + buf_filled || phi_i_1 + lcp_1 == text_length || buf_beg + buf_filled == text_length) {
                lcp_pair_writer->write((saidx_t)i_1);
                lcp_pair_writer->write((saidx_t)lcp_1);
                if (phi_pair_reader_1->empty() == false) {
                  pair_type pp = phi_pair_reader_1->read();
                  std::uint64_t next_i = (std::uint64_t)pp.m_key;
                  std::uint64_t next_phi_i = (std::uint64_t)pp.m_val;
                  lcp_1 = std::max(0L, (std::int64_t)(i_1 + lcp_1) - (std::int64_t)next_i);
                  i_1 = next_i;
                  phi_i_1 = next_phi_i;
                } else i_1 = text_length;
              } else break;
            }

            while (i_2 != text_length) {
              while (i_2 + lcp_2 < text_length && phi_i_2 + lcp_2 < buf_beg + buf_filled) {
                std::uint8_t next_char = (i_2 + lcp_2 < ext_segment_end) ?
                  segment[(i_2 + lcp_2) - segment_beg] : accessor_2->access(i_2 + lcp_2);
                if (next_char == buffer[(phi_i_2 + lcp_2) - buf_beg]) ++lcp_2;
                else break;
              }

              if (phi_i_2 + lcp_2 < buf_beg + buf_filled || i_2 + lcp_2 == text_length || buf_beg + buf_filled == text_length) {
                lcp_pair_writer->write((saidx_t)i_2);
                lcp_pair_writer->write((saidx_t)lcp_2);
                if (phi_pair_reader_2->empty() == false) {
                  pair_type pp = phi_pair_reader_2->read();
                  std::uint64_t next_phi_i = (std::uint64_t)pp.m_key;
                  std::uint64_t next_i = (std::uint64_t)pp.m_val;
                  lcp_2 = std::max(0L, (std::int64_t)(phi_i_2 + lcp_2) - (std::int64_t)next_phi_i);
                  i_2 = next_i;
                  phi_i_2 = next_phi_i;
                } else i_2 = text_length;
              } else break;
            }

            buf_beg += buf_filled;
          }

          // Print summary.
          long double process_time = utils::wclock() - process_start;
          std::uint64_t io_vol = text_reader->bytes_read() + accessor_1->bytes_read() + accessor_2->bytes_read() +
            lcp_pair_writer->bytes_written() + phi_pair_reader_1->bytes_read() + phi_pair_reader_2->bytes_read();
          total_io_volume += io_vol;
          fprintf(stderr, "\r    Process (i, Phi[i]) pairs: 100.0%%, time = %.2Lfs, I/O = %.2LfMiB/s, total I/O vol = %.2Lfn\n",
              process_time, (1.L * io_vol / (1L << 20)) / process_time, (1.L * total_io_volume) / text_length);

          delete lcp_pair_writer;
          delete accessor_1;
          delete accessor_2;
          delete phi_pair_reader_1;
          delete phi_pair_reader_2;
          delete text_reader;
#else
          lcp_pair_writer_type *lcp_pair_writer = new lcp_pair_writer_type(lcp_pairs_filename, out_buf_ram, std::max(4UL, out_buf_ram / (2UL << 20)));
          std::uint64_t io_vol = 0;

          // Parallel computation of sparse PLCP values assuming each
          // thread handles the same amount of text (but different
          // threads can handle different amount of (i, Phi[i]) pairs.
          // Each thread processes a range of text, in a sense that
          // processed are all pairs with the beginning position
          // inside that range of text.

          #pragma omp parallel num_threads(n_ranges)
          {
            std::uint64_t range_id = omp_get_thread_num();
            std::uint64_t range_beg = range_begs[range_id];
            std::uint64_t range_end = 0;
            if (range_id + 1 == n_ranges) range_end = segment_end;
            else range_end = range_begs[range_id + 1];

            std::uint64_t type_1_ptr = type_1_ptrs[range_id];
            std::uint64_t type_2_ptr = type_2_ptrs[range_id];
            if (type_1_ptr != n_pairs_type_1 || type_2_ptr != n_pairs_type_2) {
              // At least one pair (i, Ph[i]) with i >= range_beg was found.
              std::uint64_t local_buf_size = local_buf_ram / (2 * sizeof(saidx_t));
              saidx_t *local_buf = new saidx_t[2 * local_buf_size];
              std::uint64_t local_buf_filled = 0;

              text_reader_type *text_reader = new text_reader_type(text_filename, in_buf_ram / 2, 2, range_beg);

              text_accessor *accessor_1 = new text_accessor(text_filename, txt_acc_buf_ram / 2);
              text_accessor *accessor_2 = new text_accessor(text_filename, txt_acc_buf_ram / 2);

              phi_pair_reader_type *phi_pair_reader_1 = new phi_pair_reader_type(phi_pairs_filename_1, in_buf_ram / 4, 2, type_1_ptr * 2 * sizeof(saidx_t));
              phi_pair_reader_type *phi_pair_reader_2 = new phi_pair_reader_type(phi_pairs_filename_2, in_buf_ram / 4, 2, type_2_ptr * 2 * sizeof(saidx_t));

              std::uint64_t i_1 = text_length;
              std::uint64_t phi_i_1 = 0;
              std::uint64_t lcp_1 = 0;
              if (type_1_ptr != n_pairs_type_1 && (std::uint64_t)(phi_pair_reader_1->peek().m_key) < range_end) {
                pair_type pp = phi_pair_reader_1->read();
                i_1 = (std::uint64_t)pp.m_key;
                phi_i_1 = (std::uint64_t)pp.m_val;

                // Compute lower bound for PLCP[i_1] using mini sparse PLCP array.
                std::uint64_t sample_plcp_addr = i_1 / aux_sparse_plcp_sampling_rate;
                std::uint64_t sample_plcp_pos = sample_plcp_addr * aux_sparse_plcp_sampling_rate;
                std::uint64_t sample_plcp_val = aux_sparse_plcp[sample_plcp_addr];
                std::uint64_t first_i_plcp_lower_bound = (std::uint64_t)std::max(0L, (std::int64_t)sample_plcp_val -
                  ((std::int64_t)i_1 - (std::int64_t)sample_plcp_pos));
                lcp_1 = first_i_plcp_lower_bound;
              }

              std::uint64_t i_2 = text_length;
              std::uint64_t phi_i_2 = 0;
              std::uint64_t lcp_2 = 0;
              if (type_2_ptr != n_pairs_type_2 && (std::uint64_t)(phi_pair_reader_2->peek().m_key) < range_end) {
                pair_type pp = phi_pair_reader_2->read();
                phi_i_2 = (std::uint64_t)pp.m_key;
                i_2 = (std::uint64_t)pp.m_val;

                // Compute lower bound for PLCP[i_1] using mini sparse PLCP array.
                std::uint64_t sample_plcp_addr = i_2 / aux_sparse_plcp_sampling_rate;
                std::uint64_t sample_plcp_pos = sample_plcp_addr * aux_sparse_plcp_sampling_rate;
                std::uint64_t sample_plcp_val = aux_sparse_plcp[sample_plcp_addr];
                std::uint64_t first_i_plcp_lower_bound = (std::uint64_t)std::max(0L, (std::int64_t)sample_plcp_val -
                  ((std::int64_t)i_2 - (std::int64_t)sample_plcp_pos));
                lcp_2 = first_i_plcp_lower_bound;
              }

              std::uint64_t buf_beg = range_beg;
              while (buf_beg < text_length && (i_1 != text_length || i_2 != text_length)) {
                text_reader->receive_new_buffer();
                std::uint64_t buf_filled = text_reader->get_buf_filled();
                const std::uint8_t *buffer = text_reader->get_buf_ptr();

                while (i_1 != text_length) {
                  while (phi_i_1 + lcp_1 < text_length && i_1 + lcp_1 < buf_beg + buf_filled) {
                    std::uint8_t next_char = (phi_i_1 + lcp_1 < ext_segment_end) ?
                      segment[(phi_i_1 + lcp_1) - segment_beg] : accessor_1->access(phi_i_1 + lcp_1);
                    if (next_char == buffer[(i_1 + lcp_1) - buf_beg]) ++lcp_1;
                    else break;
                  }

                  if (i_1 + lcp_1 < buf_beg + buf_filled || phi_i_1 + lcp_1 == text_length || buf_beg + buf_filled == text_length) {
                    local_buf[2 * local_buf_filled] = (saidx_t)i_1;
                    local_buf[2 * local_buf_filled + 1] = (saidx_t)lcp_1;
                    ++local_buf_filled;

                    if (local_buf_filled == local_buf_size) {
                      #pragma omp critical
                      {
                        lcp_pair_writer->write(local_buf, 2 * local_buf_filled);
                      }
                      local_buf_filled = 0;
                    }

                    if (phi_pair_reader_1->empty() == false && (std::uint64_t)(phi_pair_reader_1->peek().m_key) < range_end) {
                      pair_type pp = phi_pair_reader_1->read();
                      std::uint64_t next_i = (std::uint64_t)pp.m_key;
                      std::uint64_t next_phi_i = (std::uint64_t)pp.m_val;
                      lcp_1 = std::max(0L, (std::int64_t)(i_1 + lcp_1) - (std::int64_t)next_i);
                      i_1 = next_i;
                      phi_i_1 = next_phi_i;
                    } else i_1 = text_length;
                  } else break;
                }

                while (i_2 != text_length) {
                  while (i_2 + lcp_2 < text_length && phi_i_2 + lcp_2 < buf_beg + buf_filled) {
                    std::uint8_t next_char = (i_2 + lcp_2 < ext_segment_end) ?
                      segment[(i_2 + lcp_2) - segment_beg] : accessor_2->access(i_2 + lcp_2);
                    if (next_char == buffer[(phi_i_2 + lcp_2) - buf_beg]) ++lcp_2;
                    else break;
                  }

                  if (phi_i_2 + lcp_2 < buf_beg + buf_filled || i_2 + lcp_2 == text_length || buf_beg + buf_filled == text_length) {
                    local_buf[2 * local_buf_filled] = (saidx_t)i_2;
                    local_buf[2 * local_buf_filled + 1] = (saidx_t)lcp_2;
                    ++local_buf_filled;

                    if (local_buf_filled == local_buf_size) {
                      #pragma omp critical
                      {
                        lcp_pair_writer->write(local_buf, 2 * local_buf_filled);
                      }
                      local_buf_filled = 0;
                    }

                    if (phi_pair_reader_2->empty() == false && (std::uint64_t)(phi_pair_reader_2->peek().m_key) < range_end) {
                      pair_type pp = phi_pair_reader_2->read();
                      std::uint64_t next_phi_i = (std::uint64_t)pp.m_key;
                      std::uint64_t next_i = (std::uint64_t)pp.m_val;
                      lcp_2 = std::max(0L, (std::int64_t)(phi_i_2 + lcp_2) - (std::int64_t)next_phi_i);
                      i_2 = next_i;
                      phi_i_2 = next_phi_i;
                    } else i_2 = text_length;
                  } else break;
                }

                buf_beg += buf_filled;
              }

              if (local_buf_filled > 0) {
                #pragma omp critical
                {
                  lcp_pair_writer->write(local_buf, 2 * local_buf_filled);
                }
              }

              #pragma omp critical
              {
                io_vol += accessor_1->bytes_read();
                io_vol += accessor_2->bytes_read();
                io_vol += text_reader->bytes_read();
                io_vol += phi_pair_reader_1->bytes_read();
                io_vol += phi_pair_reader_2->bytes_read();
              }

              delete[] local_buf;
              delete text_reader;
              delete accessor_1;
              delete accessor_2;
              delete phi_pair_reader_1;
              delete phi_pair_reader_2;
            }
          }

          // Print summary.
          long double process_time = utils::wclock() - process_start;
          io_vol += lcp_pair_writer->bytes_written();
          total_io_volume += io_vol;
          fprintf(stderr, "\r    Process (i, Phi[i]) pairs: 100.0%%, time = %.2Lfs, I/O = %.2LfMiB/s, total I/O vol = %.2Lfn\n",
              process_time, (1.L * io_vol / (1L << 20)) / process_time, (1.L * total_io_volume) / text_length);

          delete lcp_pair_writer;
#endif


          utils::file_delete(phi_pairs_filename_1);
          utils::file_delete(phi_pairs_filename_2);
          free(segment);
        }
      }

#ifdef _OPENMP
      // Clean up.
      for (std::uint64_t segment_id = 0; segment_id < n_segments; ++segment_id)
        delete type_1_bucket_sizes[segment_id];
      delete[] type_1_bucket_sizes;
#endif
    }

#ifdef _OPENMP
    free(aux_sparse_plcp);
#endif
  }

  // Allocate sparse PLCP.
  saidx_t *sparse_plcp = (saidx_t *)malloc(sparse_plcp_ram);

  {
    // Load pairs (i, PLCP[i]) from disk and permute in RAM.
    fprintf(stderr, "  Load pairs (i, PLCP[i]) and permute in RAM: ");
    long double lcp_pair_permute_start = utils::wclock();
    std::uint64_t io_vol = 0;

    for (std::uint64_t segment_id = 0; segment_id < n_segments; ++segment_id) {
      std::string lcp_pairs_filename = output_filename + ".lcp_pairs." + utils::intToStr(segment_id);
      std::uint64_t n_pairs = utils::file_size(lcp_pairs_filename) / (2 * sizeof(saidx_t));

#ifdef _OPENMP
      static const std::uint64_t opt_in_buf_ram = (32UL << 20);
      static const std::uint64_t opt_par_buf_ram = (4UL << 20);

      std::uint64_t in_buf_ram = opt_in_buf_ram;
      std::uint64_t par_buf_ram = opt_par_buf_ram;

      {
        std::uint64_t total_buf_ram = in_buf_ram + par_buf_ram;

        if (total_buf_ram + sparse_plcp_ram > ram_use) {
          std::uint64_t ram_budget = ram_use - sparse_plcp_ram;
          long double shrink_factor = (long double)ram_budget / (long double)total_buf_ram;
          in_buf_ram = (std::uint64_t)((long double)in_buf_ram * shrink_factor);
          par_buf_ram = (std::uint64_t)((long double)par_buf_ram * shrink_factor);
        }
      }

      typedef async_stream_reader<saidx_t> lcp_pair_reader_type;
      lcp_pair_reader_type *lcp_pair_reader = new lcp_pair_reader_type(lcp_pairs_filename, in_buf_ram, std::max(4UL, in_buf_ram / (2UL << 20)));

      std::uint64_t par_buf_size = par_buf_ram / (2 * sizeof(saidx_t));
      saidx_t *pair_buffer = new saidx_t[2 * par_buf_size];

      std::uint64_t pairs_read = 0;
      while (pairs_read < n_pairs) {
        std::uint64_t buf_filled = std::min(n_pairs - pairs_read, par_buf_size);
        lcp_pair_reader->read(pair_buffer, buf_filled * 2);

        #pragma omp parallel for
        for (std::uint64_t j = 0; j < buf_filled; ++j) {
          std::uint64_t i = pair_buffer[2 * j];
          std::uint64_t lcp = pair_buffer[2 * j + 1];

          // Invariant: PLCP[i] = lcp.
          sparse_plcp[i / plcp_sampling_rate] = lcp;
        }

        pairs_read += buf_filled;
      }

      delete[] pair_buffer;
#else
      std::uint64_t opt_in_buf_ram = (32UL << 20);
      std::uint64_t in_buf_ram = opt_in_buf_ram;
      if (sparse_plcp_ram + in_buf_ram > ram_use) {
        in_buf_ram = ram_use - sparse_plcp_ram;
      }

      typedef async_stream_reader<saidx_t> lcp_pair_reader_type;
      lcp_pair_reader_type *lcp_pair_reader = new lcp_pair_reader_type(lcp_pairs_filename, in_buf_ram, std::max(4UL, in_buf_ram / (2UL << 20)));

      for (std::uint64_t j = 0; j < n_pairs; ++j) {
        std::uint64_t i = (std::uint64_t)lcp_pair_reader->read();
        std::uint64_t lcp = (std::uint64_t)lcp_pair_reader->read();

        // Invariant: PLCP[i] = lcp.
        sparse_plcp[i / plcp_sampling_rate] = lcp;
      }
#endif

      io_vol += lcp_pair_reader->bytes_read();
      delete lcp_pair_reader;
      utils::file_delete(lcp_pairs_filename);
    }

    if ((phi_undefined_position % plcp_sampling_rate) == 0)
      sparse_plcp[phi_undefined_position / plcp_sampling_rate] = 0;

    long double lcp_pair_permute_time = utils::wclock() - lcp_pair_permute_start;
    total_io_volume += io_vol;
    fprintf(stderr, "\r  Load pairs (i, PLCP[i]) and permute in RAM: 100.0%%, time = %.1Lfs, I/O = %.2LfMiB/s, total I/O vol = %.2Lfn\n",
        lcp_pair_permute_time, ((1.L * io_vol) / (1L << 20)) / lcp_pair_permute_time, (1.L * total_io_volume) / text_length);
  }

  long double total_time = utils::wclock() - start;
  fprintf(stderr, "  Total time: %.2Lfs, total I/O vol = %.2Lfn\n\n", total_time,
      (1.L * total_io_volume) / text_length);

  return sparse_plcp;
}

}  // namespace em_sparse_phi_private

#endif  // __EM_SPARSE_PHI_SRC_COMPUTE_SPARSE_PHI_HPP_INCLUDED
