# MOVE-r
This is an optimized and parallelized implementation of the modified r-index described in [1] ([arxiv.org](https://arxiv.org/abs/2006.05104)).

## External Dependencies
- [OpenMP](https://www.openmp.org/)
- [intel TBB](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onetbb.html)

## Included Dependencies
- [ips4o](https://github.com/ips4o)
- [concurrentqueue](https://github.com/cameron314/concurrentqueue)
- [malloc_count](https://github.com/bingmann/malloc_count)
- [libsais](https://github.com/IlyaGrebnov/libsais)
- [sdsl-lite](https://github.com/simongog/sdsl-lite)

## Build Instructions
This implementation has been tested on Ubuntu 20.04 with GCC 9.4.0, libtbb-dev and libomp-dev installed.
```
clone https://github.com/LukasNalbach/MOVE-r.git
mkdir build
cd build
cmake ..
make
```
This creates five executeables in the build folder:
- MOVE-r-build
- MOVE-r-count
- MOVE-r-locate
- MOVE-r-revert
- MOVE-r-genpatterns

## Usage
### MOVE-r-build: builds the r-index (extension .MOVE-r is automatically).
```
usage: MOVE-r-build [options] <text_file>
   -o <basename>     names the index file basename.MOVE-r (default: text_file)
   -s <mode>         support mode: revert, revert-count or revert-count-locate (default:
                     revert-count-locate)
   -mu <mode>        memory usage mode: low, high or auto (low stores some data structures on
                     the disk to reduce the peak memory usage (-25% max.); auto only does this
                     if the ram capacity would be exceeded; high does not store anything on the 
                     disk; default: auto)
   -p <integer>      number of threads to use during index construction (default: all threads)
   -pr <integer>     restricts, how many threads can be used when reverting the index afterwards
                     (default: 256)
   -pd <integer>     restricts, how many threads can be used when encoding and reconstructing
                     the index (default: 256)
   -a <integer>      balancing parameter; a must be an integer and at least 2 (default: 8)
   -e <integer>      epsilon, 0 < epsilon <= 1 (default: 0.125)
   -m_idx <m_file>   file to write measurement data of the index construction to
   -m_mds <m_file>   file to write measurement data of the move-datastructures construction to
   <text_file>       input text file
```

### MOVE-r-count: count all occurrences of the input patterns.
```
usage: MOVE-r-count <index_file> <patterns_file>
   -m <m_file> <text_name>    m_file is the file to write measurement data to, text_name should be
                              the name of the original file
   -pd <integer>              the number of threads to use when reconstructing the index (default:
                              greatest possible)
   <index_file>               index file (with extension .MOVE-r)
   <patterns_file>            file in pizza&chili format containing the patterns
```

### MOVE-r-locate: locate all occurrences of the input patterns.
```
usage: MOVE-r-locate [options] <index_file> <patterns>
   -c <text_file>             check correctness of each pattern occurrence on this text file (must be
                              the indexed text file)
   -m <m_file> <text_name>    m_file is the file to write measurement data to, text_name should be the
                              name of the original file
   -o <output_file>           write pattern occurrences to this file (ASCII)
   -pd <integer>              the number of threads to use when reconstructing the index (default:
                              greatest possible)
   <index_file>               index file (with extension .MOVE-r)
   <patterns_file>            file in pizza&chili format containing the patterns
```

### MOVE-r-revert: reconstruct the original file.
```
usage: MOVE-r-revert [options] <index_file>
   -pr <integer>              number of threads to use while reverting, must not exceed the maximum
                              number of threads to use that has been specified when the index was
                              built (default: greatest possible)
   -m <m_file> <text_name>    m_file is the file to write measurement data to, text_name should be
                              the name of the original file
   -o <text_file>             output text file
   -pd <integer>              the number of threads to use when reconstructing the index (default:
                              greatest possible)
   <index_file>               index file (with extension .MOVE-r)
```

### MOVE-r-genpatterns: generate patterns from a file.
```
usage: MOVE-r-genpatterns <file> <length> <number> <patterns file> <forbidden>
       randomly extracts <number> substrings of length <length> from <file>,
       avoiding substrings containing characters in <forbidden>.
       The output file, <patterns file> has a first line of the form:
       # number=<number> length=<length> file=<file> forbidden=<forbidden>
       and then the <number> patterns come successively without any separator
```

## References
[1] Takaaki Nishimoto and Yasuo Tabei. Optimal-time queries on bwt-runs compressed indexes.
In 48th International Colloquium on Automata, Languages, and Programming (ICALP 2021),
volume 198, page 101. Schloss Dagstuhl–Leibniz-Zentrum für Informatik, 2021.
