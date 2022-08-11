/**
 * @file    em_sparse_phi_src/em_sparse_phi.hpp
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

#ifndef __EM_SPARSE_PHI_SRC_EM_SPARSE_PHI_HPP_INCLUDED
#define __EM_SPARSE_PHI_SRC_EM_SPARSE_PHI_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <string>
#include <algorithm>
#include <parallel/algorithm>
#include <omp.h>
#include <unistd.h>

#define PRODUCTION_VER
#define USE_VBYTE_ENCODING

#include "compute_sparse_plcp.hpp"
#include "compute_final_lcp.hpp"
#include "compute_lcp_lower_bounds.hpp"
#include "compute_lcp_delta.hpp"
#include "utils.hpp"


namespace em_sparse_phi_private {

template<typename saidx_t>
void write_sparse_plcp_to_file(saidx_t *sparse_plcp, std::uint64_t text_length,
    std::uint64_t sampling_rate, std::string filename, std::uint64_t &total_io_volume) {
  std::uint64_t size = (text_length + sampling_rate - 1) / sampling_rate;
  fprintf(stderr, "  Write sparse PLCP to disk: ");
  long double start = utils::wclock();

  // Space-efficient encoding using at most n(1+1/q) bits,
  // i.e., assuming q >= 5 at most 0.15n bytes.
  typedef async_stream_writer<std::uint64_t> writer_type;
  writer_type *writer = new writer_type(filename, (1UL << 20), 2);
  std::uint64_t buffer = 0, filled = 0, prev = 0;
  for (std::uint64_t j = 0; j < size; ++j) {
    std::uint64_t cur = (std::uint64_t)sparse_plcp[j];
    std::uint64_t diff = cur - std::max(0L, (std::int64_t)prev - (std::int64_t)sampling_rate);
    prev = cur;
    while (diff > 0) {
      ++filled;
      if (filled == 64) {
        writer->write(buffer);
        buffer = 0;
        filled = 0;
      }
      --diff;
    }
    buffer |= (1UL << filled);
    ++filled;
    if (filled == 64) {
      writer->write(buffer);
      buffer = 0;
      filled = 0;
    }
  }
  if (filled > 0)
    writer->write(buffer);
  long double write_time = utils::wclock() - start;
  std::uint64_t io_vol = writer->bytes_written();
  total_io_volume += io_vol;
  fprintf(stderr, "time = %.1Lfs, I/O = %.1LfMiB/s, total I/O vol = %.2Lfn\n", write_time,
      ((1.L * io_vol) / (1L << 20)) / write_time, (1.L * total_io_volume) / text_length);
  delete writer;
}

template<typename saidx_t>
void read_sparse_plcp_from_file(saidx_t *sparse_plcp, std::uint64_t text_length,
    std::uint64_t sampling_rate, std::string filename, std::uint64_t &total_io_volume) {
  std::uint64_t size = (text_length + sampling_rate - 1) / sampling_rate;
  fprintf(stderr, "  Read sparse PLCP array from disk: ");
  long double start = utils::wclock();

  // Space-efficient encoding using at most n(1+1/q) bits,
  // i.e., assuming q >= 5 at most 0.15n bytes.
  typedef async_stream_reader<std::uint64_t> reader_type;
  reader_type *reader = new reader_type(filename, (1UL << 20), 2);
  std::uint64_t buffer = 0, filled = 0, prev = 0;
  for (std::uint64_t j = 0; j < size; ++j) {
    std::uint64_t diff = 0;
    while (true) {
      if (filled == 0) {
        buffer = reader->read();
        filled = 64;
      }
      if (buffer & 1) {
        --filled;
        buffer >>= 1;
        break;
      } else {
        --filled;
        buffer >>= 1;
        ++diff;
      }
    }
    std::uint64_t cur = std::max(0L, (std::int64_t)prev - (std::int64_t)sampling_rate) + diff;
    prev = cur;
    sparse_plcp[j] = (saidx_t)cur;
  }
  long double read_time = utils::wclock() - start;
  std::uint64_t io_vol = reader->bytes_read();
  total_io_volume += io_vol;
  fprintf(stderr, "time = %.1Lfs, I/O = %.1LfMiB/s, total I/O vol = %.2Lfn\n", read_time,
      ((1.L * io_vol) / (1L << 20)) / read_time, (1.L * total_io_volume) / text_length);
  delete reader;
}

template<typename saidx_t>
void em_sparse_phi(std::string text_filename, std::string sa_filename,
    std::string output_filename, std::uint64_t ram_use, bool near_inplace_mode = false) {
  srand(time(0) + getpid());
  utils::drop_disk_pages(text_filename);
  utils::drop_disk_pages(sa_filename);
  std::string sparse_plcp_filename = output_filename + ".sparse_plcp";
  long double start = utils::wclock();
  std::uint64_t total_io_volume = 0;

  // Compute some basic parameters of computation.
  std::uint64_t text_length = utils::file_size(text_filename);
  long double text_to_ram_ratio = (long double)text_length / (long double)ram_use;
  enum partitioning_t { text_partitioning, lex_partitioning };
  partitioning_t partitioning_type = lex_partitioning;
  if (text_to_ram_ratio > 5.0L && near_inplace_mode == true)
    partitioning_type = text_partitioning;

  // Compute max halfsegment size and buffer
  // sizes used when computing lcp deltas.
  static const std::uint64_t opt_overflow_size_compute_lcp_delta = (1UL << 20);
  static const std::uint64_t opt_in_buf_ram_compute_lcp_delta = (32UL << 20);
  static const std::uint64_t opt_out_buf_ram_compute_lcp_delta = (4UL << 20);
  static const std::uint64_t opt_local_buf_ram_compute_lcp_delta = (64UL << 20);

  std::uint64_t max_halfsegment_size = (std::uint64_t)((long double)ram_use * 0.45L);  // min acceptable value
  std::uint64_t overflow_size_compute_lcp_delta = opt_overflow_size_compute_lcp_delta;
  std::uint64_t in_buf_ram_compute_lcp_delta = opt_in_buf_ram_compute_lcp_delta;
  std::uint64_t out_buf_ram_compute_lcp_delta = opt_out_buf_ram_compute_lcp_delta;
  std::uint64_t local_buf_ram_compute_lcp_delta = opt_local_buf_ram_compute_lcp_delta;

  {
    std::uint64_t total_buf_size_compute_lcp_delta = 2 * overflow_size_compute_lcp_delta +
      in_buf_ram_compute_lcp_delta + out_buf_ram_compute_lcp_delta + local_buf_ram_compute_lcp_delta;
    if (2 * max_halfsegment_size + total_buf_size_compute_lcp_delta > ram_use) {
      std::uint64_t ram_budget = ram_use - 2 * max_halfsegment_size;
      long double shrink_factor = (long double)ram_budget / (long double)total_buf_size_compute_lcp_delta;
      overflow_size_compute_lcp_delta = (std::uint64_t)((long double)overflow_size_compute_lcp_delta * shrink_factor);
      in_buf_ram_compute_lcp_delta = (std::uint64_t)((long double)in_buf_ram_compute_lcp_delta * shrink_factor);
      out_buf_ram_compute_lcp_delta = (std::uint64_t)((long double)out_buf_ram_compute_lcp_delta * shrink_factor);
      local_buf_ram_compute_lcp_delta = (std::uint64_t)((long double)local_buf_ram_compute_lcp_delta * shrink_factor);
    } else max_halfsegment_size = (ram_use - total_buf_size_compute_lcp_delta) / 2;
  }

  std::uint64_t n_halfsegments = (text_length + max_halfsegment_size - 1) / max_halfsegment_size;
  std::uint64_t n_different_halfseg_pairs = (n_halfsegments * (n_halfsegments + 1)) / 2;

  // Compute PLCP sampling rate and size of buffer used for each halfsegment pair.
  static const std::uint64_t opt_in_buf_ram = (32UL << 20);
  static const std::uint64_t opt_out_buf_ram = (4UL << 20);
  std::uint64_t opt_local_buf_ram = 0;
#ifdef _OPENMP
  if (partitioning_type == lex_partitioning)
    opt_local_buf_ram = (100UL << 20);
  else opt_local_buf_ram = (140UL << 20);
#else
  opt_local_buf_ram = (64UL << 20);
#endif

  std::uint64_t sparse_plcp_and_halfseg_buffers_ram = (std::uint64_t)((long double)ram_use * 0.9L);  // min acceptable value
  std::uint64_t sparse_plcp_ram = (std::uint64_t)((long double)ram_use * 0.2L);                      // min acceptable value
  std::uint64_t halfseg_buffers_ram = n_different_halfseg_pairs * (4UL << 20);
  std::uint64_t in_buf_ram = opt_in_buf_ram;
  std::uint64_t out_buf_ram = opt_out_buf_ram;
  std::uint64_t local_buf_ram = opt_local_buf_ram;

  {
    std::uint64_t total_buf_ram = in_buf_ram + out_buf_ram + local_buf_ram;
    if (sparse_plcp_and_halfseg_buffers_ram + total_buf_ram > ram_use) {
      // Shrink reader/writer/local buffers.
      {
        std::uint64_t ram_budget = ram_use - sparse_plcp_and_halfseg_buffers_ram;
        long double shrink_factor = (long double)ram_budget / (long double)total_buf_ram;
        in_buf_ram = (std::uint64_t)((long double)in_buf_ram * shrink_factor);
        out_buf_ram = (std::uint64_t)((long double)out_buf_ram * shrink_factor);
        local_buf_ram = (std::uint64_t)((long double)local_buf_ram * shrink_factor);
      }

      if (sparse_plcp_ram + halfseg_buffers_ram > sparse_plcp_and_halfseg_buffers_ram) {
        // Decrese RAM for halfseg buffers.
        halfseg_buffers_ram = sparse_plcp_and_halfseg_buffers_ram - sparse_plcp_ram;
      } else {
        // Increase the RAM for sparse PLCP.
        sparse_plcp_ram = sparse_plcp_and_halfseg_buffers_ram - halfseg_buffers_ram;
      }
    } else {
      // Increase the space for sparse PLCP and halfseg buffers.
      sparse_plcp_and_halfseg_buffers_ram = ram_use - total_buf_ram;

      if (sparse_plcp_ram + halfseg_buffers_ram > sparse_plcp_and_halfseg_buffers_ram) {
        // Decrese RAM for halfseg buffers.
        halfseg_buffers_ram = sparse_plcp_and_halfseg_buffers_ram - sparse_plcp_ram;
      } else {
        // Increase the RAM for sparse PLCP.
        sparse_plcp_ram = sparse_plcp_and_halfseg_buffers_ram - halfseg_buffers_ram;
      }
    }
  }

  std::uint64_t sparse_plcp_size = std::max(1UL, sparse_plcp_ram / sizeof(saidx_t));
  std::uint64_t plcp_sampling_rate = (text_length + sparse_plcp_size - 1) / sparse_plcp_size;

  // Turn paths absolute.
  text_filename = utils::absolute_path(text_filename);
  sa_filename = utils::absolute_path(sa_filename);
  output_filename = utils::absolute_path(output_filename);

  // Print summary of basic parameters.
  fprintf(stderr, "Text filename = %s\n", text_filename.c_str());
  fprintf(stderr, "SA filename = %s\n", sa_filename.c_str());
  fprintf(stderr, "Output filename = %s\n", output_filename.c_str());
  fprintf(stderr, "Integer size = %lu bytes\n", sizeof(saidx_t));
  fprintf(stderr, "Text length = %lu (%.2LfMiB)\n", text_length, 1.L * text_length / (1 << 20));
  fprintf(stderr, "Text size / ram_use = %.2Lf\n", text_to_ram_ratio);
  fprintf(stderr, "RAM use = %lu (%.2LfMiB)\n", ram_use, ram_use / (1024.L * 1024));
  fprintf(stderr, "Max halfsegment size = %lu (%.2LfMiB)\n", max_halfsegment_size, 1.L * max_halfsegment_size / (1L << 20));
  fprintf(stderr, "Number of halfsegments = %lu\n", n_halfsegments);
  fprintf(stderr, "Number of halfsegment pairs = %lu\n", n_different_halfseg_pairs);
  fprintf(stderr, "RAM for halfsegment buffers = %lu (%.2LfMiB)\n", halfseg_buffers_ram, 1.L * halfseg_buffers_ram / (1L << 20));
  fprintf(stderr, "RAM for sparse PLCP = %lu (%.2LfMiB)\n", sparse_plcp_ram, 1.L * sparse_plcp_ram / (1L << 20));
  fprintf(stderr, "PLCP sampling rate = %lu\n", plcp_sampling_rate);
#ifdef _OPENMP
  fprintf(stderr, "Max number of threads = %d\n", omp_get_max_threads());
#endif
  fprintf(stderr, "Near inplace mode = %s\n", near_inplace_mode ? "TRUE" : "FALSE");
  fprintf(stderr, "\n");

  // Allocate array with info about the number of items in each paris of halfsegments.
  std::uint64_t **pairs_count = NULL;
  if (partitioning_type == text_partitioning) {
    pairs_count = new std::uint64_t*[n_halfsegments];
    for (std::uint64_t j = 0; j < n_halfsegments; ++j) {
      pairs_count[j] = new std::uint64_t[n_halfsegments];
      std::fill(pairs_count[j], pairs_count[j] + n_halfsegments, 0UL);
    }
  }

  // Step 1: compute sparse PLCP array.
  saidx_t *sparse_plcp = NULL;
  sparse_plcp = compute_sparse_plcp<saidx_t>(text_filename, sa_filename,
      output_filename, text_length, max_halfsegment_size, pairs_count,
      plcp_sampling_rate, ram_use, total_io_volume);

  std::uint64_t lcp_sum = 0;
  std::uint64_t max_lcp = 0;

  if (partitioning_type == lex_partitioning) {
    // Step 2: compute lcp-delta values.
    std::uint64_t max_sa_parts = (near_inplace_mode ? 2 : 1);
    std::uint64_t max_sa_part_size = (text_length + max_sa_parts - 1) / max_sa_parts;
    std::uint64_t n_sa_parts = (text_length + max_sa_part_size - 1) / max_sa_part_size;
    for (std::uint64_t sa_part_id = 0; sa_part_id < n_sa_parts; ++sa_part_id) {
      std::uint64_t sa_part_beg = sa_part_id * max_sa_part_size;
      std::uint64_t sa_part_end = std::min(sa_part_beg + max_sa_part_size, text_length);
      long double process_part_start = utils::wclock();
      fprintf(stderr, "Process lex-part %lu/%lu:\n", sa_part_id + 1, n_sa_parts);

      // Step 2.1: compute LCP lower bounds
      compute_lcp_lower_bounds_lex_partitioning(sa_filename, output_filename, text_length,
          ram_use, sparse_plcp, plcp_sampling_rate, max_halfsegment_size, sa_part_beg,
          sa_part_end, halfseg_buffers_ram, in_buf_ram, local_buf_ram, overflow_size_compute_lcp_delta,
          total_io_volume);

      // Step 2.2: write sparse PLCP to disk.
      write_sparse_plcp_to_file(sparse_plcp, text_length,
          plcp_sampling_rate, sparse_plcp_filename, total_io_volume);
      free(sparse_plcp);

      // Step 2.3: compute LCP delta-values.
      compute_lcp_delta_lex_partitioning<saidx_t>(text_filename, output_filename, text_length,
          overflow_size_compute_lcp_delta, max_halfsegment_size, (sa_part_id == 0),
          in_buf_ram_compute_lcp_delta, out_buf_ram_compute_lcp_delta, local_buf_ram_compute_lcp_delta,
          total_io_volume);

      // Step 2.4: load sparse PLCP from disk.
      sparse_plcp = (saidx_t *)malloc(sparse_plcp_size * sizeof(saidx_t));
      read_sparse_plcp_from_file(sparse_plcp, text_length,
          plcp_sampling_rate, sparse_plcp_filename, total_io_volume);
      fprintf(stderr, "  Total time: %.2Lfs, total I/O vol = %.2Lfn\n", utils::wclock() - process_part_start,
          (1.L * total_io_volume) / text_length);
    }
    utils::file_delete(sparse_plcp_filename);

    // Step 3: compute the final LCP array.
    compute_final_lcp_lex_partitioning(sa_filename, output_filename, text_length,
        sparse_plcp, plcp_sampling_rate, max_halfsegment_size, halfseg_buffers_ram,
        in_buf_ram, out_buf_ram, local_buf_ram, overflow_size_compute_lcp_delta, lcp_sum,
        max_lcp, total_io_volume);
    free(sparse_plcp);
  } else {  // Use text partitioning.
    static const std::uint64_t n_parts = 2;
    std::uint64_t items_per_part = (text_length - 1) / n_parts;

    // Compute the number of pairs to be processed in each of
    // the halfsegments pairs for all text parts.
    std::uint64_t ***items_per_halfseg_pair = new std::uint64_t**[n_parts];
    for (std::uint64_t part_id = 0; part_id < n_parts; ++part_id) {
      items_per_halfseg_pair[part_id] = new std::uint64_t*[n_halfsegments];
      for (std::uint64_t i = 0; i < n_halfsegments; ++i) {
        items_per_halfseg_pair[part_id][i] = new std::uint64_t[n_halfsegments];
        std::fill(items_per_halfseg_pair[part_id][i], items_per_halfseg_pair[part_id][i] + n_halfsegments, 0UL);
      }

      std::uint64_t items_in_this_part = items_per_part;
      if (part_id + 1 == n_parts)
        items_in_this_part = text_length - 1 - (n_parts - 1) * items_per_part;

      for (std::uint64_t diff = 0; diff < n_halfsegments; ++diff) {
        for (std::uint64_t j = n_halfsegments; j > diff; --j) {
          std::uint64_t i = (j - 1) - diff;
          std::uint64_t count = std::min(items_in_this_part, pairs_count[i][j - 1]);
          items_per_halfseg_pair[part_id][i][j - 1] = count;
          pairs_count[i][j - 1] -= count;
          items_in_this_part -= count;
          if (items_in_this_part == 0) break;
        }
        if (items_in_this_part == 0) break;
      }
    }

    for (std::uint64_t part_id = 0; part_id < n_parts; ++part_id) {
      long double process_part_start = utils::wclock();
      fprintf(stderr, "Process text-part %lu/%lu:\n", part_id + 1, n_parts);

      // Step 2.1: compute LCP lower bounds.
      compute_lcp_lower_bounds_text_partitioning(sa_filename, output_filename, text_length,
          ram_use, sparse_plcp, plcp_sampling_rate, max_halfsegment_size, halfseg_buffers_ram,
          in_buf_ram, local_buf_ram, overflow_size_compute_lcp_delta, items_per_halfseg_pair,
          part_id, total_io_volume);

      // Step 2.2: write sparse PLCP to disk.
      write_sparse_plcp_to_file(sparse_plcp, text_length,
          plcp_sampling_rate, sparse_plcp_filename, total_io_volume);
      free(sparse_plcp);

      // Step 2.3: compute LCP delta-values.
      compute_lcp_delta_text_partitioning<saidx_t>(text_filename, output_filename, text_length,
          overflow_size_compute_lcp_delta, max_halfsegment_size, in_buf_ram_compute_lcp_delta,
          out_buf_ram_compute_lcp_delta, local_buf_ram_compute_lcp_delta, part_id, total_io_volume);

      // Step 2.4: load sparse PLCP from disk.
      sparse_plcp = (saidx_t *)malloc(sparse_plcp_size * sizeof(saidx_t));
      read_sparse_plcp_from_file(sparse_plcp, text_length,
          plcp_sampling_rate, sparse_plcp_filename, total_io_volume);
      fprintf(stderr, "  Total time: %.2Lfs, total I/O vol = %.2Lfn\n", utils::wclock() - process_part_start,
          (1.L * total_io_volume) / text_length);
    }
    utils::file_delete(sparse_plcp_filename);

    // Step 3: compute the final LCP array.
    compute_final_lcp_text_partitioning(sa_filename, output_filename, text_length,
        sparse_plcp, plcp_sampling_rate, max_halfsegment_size, halfseg_buffers_ram,
        in_buf_ram, out_buf_ram, local_buf_ram, overflow_size_compute_lcp_delta, n_parts,
        items_per_halfseg_pair, lcp_sum, max_lcp, total_io_volume);

    // Clean up.
    for (std::uint64_t part_id = 0; part_id < n_parts; ++part_id) {
      for (std::uint64_t i = 0; i < n_halfsegments; ++i)
        delete[] items_per_halfseg_pair[part_id][i];
      delete[] items_per_halfseg_pair[part_id];
    }
    delete[] items_per_halfseg_pair;
    free(sparse_plcp);
  }

  // Clean up.
  if (pairs_count != NULL) {
    for (std::uint64_t j = 0; j < n_halfsegments; ++j)
      delete[] pairs_count[j];
    delete[] pairs_count;
  }

  // Print summary.
  long double total_time = utils::wclock() - start;
  long double avg_lcp = (long double)lcp_sum / text_length;
  fprintf(stderr, "\n\nComputation finished. Summary:\n");
  fprintf(stderr, "  elapsed time = %.2Lfs (%.3Lfs/MiB of text)\n", total_time, total_time / (1.L * text_length / (1L << 20)));
  fprintf(stderr, "  speed = %.2LfMiB of text/s\n", (1.L * text_length / (1L << 20)) / total_time);
  fprintf(stderr, "  I/O volume = %lu (%.2Lfbytes/input symbol)\n", total_io_volume, (1.L * total_io_volume) / text_length);
  fprintf(stderr, "  average LCP = %.2Lf\n", avg_lcp);
  fprintf(stderr, "  maximal LCP = %lu\n", max_lcp);
}

}  // namespace em_sparse_phi_private

template<typename saidx_t>
void em_sparse_phi(std::string text_filename, std::string sa_filename,
    std::string output_filename, std::uint64_t ram_use) {
  em_sparse_phi_private::em_sparse_phi<saidx_t>(text_filename,
      sa_filename, output_filename, ram_use);
}

#endif  // __EM_SPARSE_PHI_SRC_EM_SPARSE_PHI_HPP_INCLUDED
