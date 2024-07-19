BDA-index: Text Indexing for Long Patterns
===

Description
-----------

This repository maintains a time- and space-efficient construction algorithm of the <b>BDA-index</b>, a text index for long patterns introduced by [Loukides, Pissis, and Sweering](https://doi.org/10.1109/TKDE.2022.3231780).
This new construction relies on a linear-time algorithm for computing the bd-anchors and on a semi-external memory implementation to
construct the index in small space. These algorithms are now outperformed by the [rrBDA-index](https://github.com/lorrainea/rrBDA-index).

Requirements
-----------
* A GNU/Linux system
* A modern C++11 ready compiler such as g++ version 4.9 or higher

How to use
----------

The <b>BDA-index</b> comes in two flavours: BDA-index I and BDA-index II. BDA-index I has provably near-optimal queries but relies on a 2D range reporting data structure (based on wavelet trees), which may be slow in practice.
BDA-index II is a heuristic that is considerably faster, especially when the number of occurrences is high; and it is also smaller.

<b>INPUT</b>: A file containing the text and a file containing a set of patterns seperated by a new line.

<b>OUTPUT</b>: A file containing the set of patterns and the starting position of their occurrences within the text.

### Installation

```
cd BDA-index_I
./pre-install.sh
make -f Makefile.32-bit.gc
```

```
cd BDA-index_II
./pre-install.sh
make -f Makefile.32-bit.gc
```

### Usage

```
./bda-index_I <text_file> <ell> <pattern_file> <output_filename> <ram_use> <block_size> <index_filename>
./bda-index_II <text_file> <ell> <pattern_file> <output_filename> <ram_use> <block_size> <index_filename>

<text_file> - name of input text file.
<ell> - lower bound on the length of input patterns to consider. 
<pattern_file> - name of input file containing the patterns.
<output_filename> - name of output file, where pattern occurrences will be output.
<ram_use> - RAM usage for external SA and LCP array construction (MiB).
<block_size> - size of block to use for constructing the bd-anchors (bytes).
<index_filename> - name of the index file to be used (if it exists) otherwise to be created.
```

### Examples

```
 $ ./bda-index_I ./data/text 3 ./data/patterns out 1024 10 index
 $ ./bda-index_II ./data/text 3 ./data/patterns out 1024 10 index
```

Datasets
--------

The Pizza&Chili datasets and a sample of the patterns used in the experimental analysis can be found at https://bit.ly/3pdViRs.

Citation
--------

Lorraine A. K. Ayad, Grigorios Loukides, and Solon P. Pissis. 2023. Text Indexing for Long Patterns: Anchors are All you Need. Proc. VLDB Endow. 16, 9 (May 2023), 2117–2131. https://doi.org/10.14778/3598581.3598586

License
--------

GNU GPLv3 License; Copyright (C) 2023 Lorraine A. K. Ayad, Grigorios Loukides and Solon P. Pissis.
