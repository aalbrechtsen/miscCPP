# refFinder

A small, fast C program to extract bases from a fasta file. It can work as a stand-alone program, or you can use its API to look up reference bases from your own program.

It uses an indexed fasta file. If there is no `.fai` index next to the fasta file, refFinder builds one the first time it runs.

## Install

```bash
git clone https://github.com/aalbrechtsen/miscCPP.git
cd miscCPP/refFinder
make
```

It needs a C/C++ compiler and zlib.

## Stand-alone

Pipe a list of positions (chromosome and position, tab-separated) into refFinder, and give it the fasta file:

```bash
./refFinder ref.fa [full] [inputIsZero] < positions.txt
```

There are two optional flags:

- **full**: print the chromosome and position as well as the base. Without it, only the base is printed.
- **inputIsZero**: the input positions are 0-based. By default they are 1-based, as in samtools.

### Example

Make a list of chromosome, position and reference base with samtools:

```bash
samtools mpileup -b smallBam.filelist -f hg19NoChr.fa | cut -f1-3 > small.sam
```

Use refFinder to find the bases for each position in `small.sam`, and check they match:

```bash
cut -f1-2 small.sam | ./refFinder hg19NoChr.fa full > tst
cmp tst small.sam
```

Only the bases:

```
cut -f1-2 small.sam | ./refFinder hg19NoChr.fa | head
a
g
c
t
a
c
t
c
g
g
```

With the chromosome and position:

```
cut -f1-2 small.sam | ./refFinder hg19NoChr.fa full | head
1	13999902	a
1	13999903	g
1	13999904	c
1	13999905	t
1	13999906	a
1	13999907	c
1	13999908	t
1	13999909	c
1	13999910	g
1	13999911	g
```

If the positions are 0-based instead of 1-based:

```
cut -f1-2 small.sam | ./refFinder hg19NoChr.fa full inputIsZero | head
1	13999902	g
1	13999903	c
1	13999904	t
1	13999905	a
1	13999906	c
1	13999907	t
1	13999908	c
1	13999909	g
1	13999910	g
1	13999911	g
```

## API

```cpp
#include "refFinder.h"

perFasta *pf = init("hg19.fa");
char refbase = getchar("chr20", 130224101, pf);
// refbase now contains the reference base for chr20 at position 130,224,101
destroy(pf);
```

Link with `refFinder.o` and `-lz`:

```bash
g++ sampleProg.cpp refFinder.o -lz
```

## Known bugs

- It does not check whether the reference file exists.

Moved here in 2026 from the popgen.dk software wiki.
