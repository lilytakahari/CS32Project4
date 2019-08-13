## UCLA CS 32 Winter 2019 Project 4 - Gee-nomics

```Correctness score: 82/95```

This program processes genetic data. It maintains a library of genomes (in memory) from multiple organisms; the user can add new genomes to this library. The user can then search the library for a specified DNA sequence and identify all genomes in the library that contain that DNA sequence or any SNiP of it. The user can also present the genome/DNA of a new organism and quickly identify all genomes in the library that have a high percentage of matching DNA with the presented genome.

My implementation of what was required in the spec is contained in the `Trie.h, Genome.cpp`, and `GenomeMatcher.cpp` files.

### How to Run the Program
Download the `skeleton.zip` file. Replace `Trie.h, Genome.cpp`, and `GenomeMatcher.cpp` with the files above. Edit the string literal in `main.cpp` to be the file path to the folder that contains the genome text files. Build the program and run the executable (from the command line if necessary). 

These are the commands that the test harness accepts:

```
c - create new genome library               s - find matching SNiPs
a - add one genome manually                 r	- find related genomes (manual)
l - load one data file                      f	- find related genomes (file)
d - load all provided data files            ?	- show this menu
e - find matches exactly                    q	- quit
```
