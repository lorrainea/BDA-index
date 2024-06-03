BDA-index_I and BDA-index_II
===

<b>Installation</b>: To install and compile BDA-index_I or BDA-index_II, read the INSTALL file within the BDA-index_I or BDA-index_II folders.

<b>INPUT</b>: A file containing a single text and a file containing a set of patterns seperated by a new line.

<b>OUTPUT</b>: A file containing the set of patterns and the starting position of their occurrences within the text.


```
Usage: 
./bda-index_I <text_file> <ell> <pattern_file> <output_filename> <ram_use> <block_size> <index_filename>
./bda-index_II <text_file> <ell> <pattern_file> <output_filename> <ram_use> <block_size> <index_filename>

<text_file> - name of input text file.
<ell> - minimum size of pattern to consider searching for within text. 
<pattern_file> - name of input file containing patterns.
<output_filename> - name of output file where pattern occurrences will be placed.
<ram_use> - ram usage for external SA and LCP (Mbit).
<block_size> - size of block size b to use (bytes).
<index_filename> - name of the index file to be used if they exist otherwise to be created.
```

<b>Examples</b>
```
 $ ./bda-index_I ./data/text 3 ./data/patterns out 150 10 index
 $ ./bda-index_II ./data/text 3 ./data/patterns out 150 10 index
```

<b>Data</b>

The Pizza&Chili datasets and patterns used in the experimental analysis can be found at https://bit.ly/3pdViRs.

<b>Citation</b>

Lorraine A. K. Ayad, Grigorios Loukides, and Solon P. Pissis. 2023. Text Indexing for Long Patterns: Anchors are All you Need. Proc. VLDB Endow. 16, 9 (May 2023), 2117–2131. https://doi.org/10.14778/3598581.3598586

<b>License</b>

GNU GPLv3 License; Copyright (C) 2023 Lorraine A. K. Ayad, Grigorios Loukides and Solon P. Pissis.
