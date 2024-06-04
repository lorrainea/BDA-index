BDA-index_I
===

<b>Installation</b>: To install and compile BDA-index_I, read the INSTALL file.

<b>INPUT</b>: A file containing a single text and a file containing a set of patterns seperated by a new line.

<b>OUTPUT</b>: A file containing the set of patterns and the starting position of their occurrences within the text.


```
Usage:
./bda-index_I <text_file> <ell> <pattern_file> <output_filename> <ram_use> <block_size> <index_filename>

<text_file> - name of input text file
<ell> - minimum size of pattern to consider searching for within text.
<pattern_file> - name of input file containing patterns
<output_filename> - name of output file where pattern occurrences will be placed.
<ram_use> - ram usage for external SA and LCP
<block_size> - size of block size b to use.
<index_filename> - name of the index file to be used if they exist otherwise to be created.
```

<b>Example</b>
```
 $ ./bda-index_I ./data/text 3 ./data/patterns out 150 10 index
```

<b>WORK-IN-PROGRESS</b>

This index is currently much slower than BDA-index_II, especially when the number of occurrences is high.
We are in the process of re-implementing the 2D range reporting algorithm to make it more competitive to BDA-index_II.

<b>License</b>

GNU GPLv3 License; Copyright (C) 2023 Lorraine A. K. Ayad, Grigorios Loukides and Solon P. Pissis.
