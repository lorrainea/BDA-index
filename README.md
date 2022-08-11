BDA-index_I and BDA-index_II
===

<b>Installation</b>: To install and compile BDA-index_I or BDA-index_II, read the INSTALL file within the BDA-index_I or BDA-index_II folders.

<b>INPUT</b>: A file containing a single text and a file containing a set of patterns seperated by a new line.

<b>OUTPUT</b>: A file containing the set of patterns and the starting position of their occurrences within the text.


```
Usage: 
./index <text_file> <ell> <pattern_file> <output_filename> <ram_use> <block_size>

<text_file> - name of input text file
<ell> - minimum size of pattern to consider searching for within text. 
<pattern_file> - name of input file containing patterns
<output_filename> - name of output file.
<ram_use> - ram usage for external SA and LCP
<block_size> - size of block size b to use.
```

<b>Examples</b>
```
 $ ./bda-index_I ./data/text 3 ./data/patterns out 150 10
 $ ./bda-index_II ./data/text 3 ./data/patterns out 150 10
```

<b>Data</b>
All data used in the experimental analysis can be found at https://www.dropbox.com/home/bd-index_ALENEX23.
