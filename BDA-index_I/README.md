BDA-index I
===

How to use
----------

### Installation

```
cd BDA-index_I
./pre-install.sh
make -f Makefile.32-bit.gc
```

### Usage

```
./bda-index_I <text_file> <ell> <pattern_file> <output_filename> <ram_use> <block_size> <index_filename>

<text_file> - name of input text file.
<ell> - lower bound on the length of input patterns to consider. 
<pattern_file> - name of input file containing the patterns.
<output_filename> - name of output file, where pattern occurrences will be output.
<ram_use> - RAM usage for external SA and LCP array construction (Mbits).
<block_size> - size of block to use for constructing the bd-anchors (bytes).
<index_filename> - name of the index file to be used (if it exists) otherwise to be created.
```

### Examples

```
 $ ./bda-index_I ./data/text 3 ./data/patterns out 1024 10 index
```

Work-in-progress
----------------

This index is currently considerably slower in practice than BDA-index II, especially when the number of occurrences is high.
We are in the process of re-implementing the 2D range reporting algorithm to make it more competitive to BDA-index II.

License
-------

GNU GPLv3 License; Copyright (C) 2023 Lorraine A. K. Ayad, Grigorios Loukides and Solon P. Pissis.
