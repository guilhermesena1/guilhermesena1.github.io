---
layout: post
math: true
title:  "The most naive phylogenetic reconstruction algorithm"
date:   '2022-07-18'
categories: phylogenetics C++ UCSC-genome-browser
---

The full source code for this post can be found on
[GitHub](https://github.com/guilhermesena1/phylogenetic-by-genome).

This post was motivated by the following question: *what would be the
simplest problem we could formulate to introduce a computer science
student to the field of bioinformatics?* This is the answer I came up
with, and is a (hopefully) fun attempt to analyze easily accessible
biological data with very little code or effort. The question we will
address is the following:

**given the genome sequences of a set of species, can we reconstruct
their speciation history?**

Below we will provide a brief introduction to the problem, then
formulate it more formally and solve it with some code, after which we
will discuss the code's correctness, efficiency and potential
limitations.

# Introduction

Two species that are phenotypically very similar are expected to have
had a recent common ancestor. In other words, not too long ago (which
in "evolution time" may mean hundreds of thousands of years ago), a
speciation event occurred, and the two resulting species evolved
independently, resulting in the observable differences. The longer ago
the speciation event, the more different species are today.

In theory we can (and have) reconstruct speciation events by observing
physical traits (phenotypes) of a pair of species. For example, we do
not need to look at DNA to have a strong conviction that the Chinese
hamster and the mouse are close siblings evolutionarily. Just look at
them! (and take a second to appreciate how adorable they are)

![a Chinese hamster (left) and a mouse
(right)](https://i.ibb.co/L5SGKd5/post-1.png)

Conversely, relying purely on traits can be misleading. Bats and birds
have wings and other similar physical traits, but bats are more
closely related to cats than to birds. Similarly, mammalian aquatic
animals like whales and dolphins are closer to wolves and alpacas than
to fish (by "closer" I mean that the speciation event of the most
recent common ancestor occurred more recently). These are examples of
[convergent
evolution](https://en.wikipedia.org/wiki/Convergent_evolution), where
two species have similar traits not necessarily because they are
related by recent common ancestry, but because they are under similar
selective pressure to acquire those traits.

With recent genome sequencing technologies, we can objectively measure
similarity between two species: simply measure the "similarity"
between the genome sequences, where the sequences are strings of
characters (A, C, G, T) representing chromosomes. This is easier said
than done. First, we need to define what "similarity" of sequences
means. When speciation occurs, independent genome divergence occurs
through a sequence of transformations: mutations, nucleotide
insertions, deletions, duplications of certain sequences, and possible
translocations and inversions of large portions of the genome. Simply
put, pieces of DNA can move around, mutate, copy-paste themselves, and
other small operations that can significantly change the resulting
sequence. Can you think of a similarity metric between two strings
that would identify similar genomic sequences under these potential
change operations? Think if edit distance would be enough. Yes? No? To
me at least it is not immediately obvious. Simple translocations of
genome sequences can significantly change the edit distance.

In most bioinformatics methods used today, comparison of sequences is
done through (possibly some modification of) the [local alignment of
two
sequences](https://www.sciencedirect.com/science/article/abs/pii/0022283681900875).
Given sequences $A$ and $B$, of sizes $|A|$ and $|B|$, the local
alignment of $A$ and $B$ is the most similar substring of $A$ to a
substring of $B$. Objectively, we devise a scoring scheme that rewards
similar letters in $A$ and penalize letter differences, insertions and
deletions. The local alignment can be found through the Smith-Waterman
algorithm in $O(|A| |B|)$ by filling a $|A| \times |B|$ matrix. The
local alignment score $S[i, j]$ of the first $i$ characters of $A$ and
the first $j$ characters of $B$ can be found by the following dynamic
programming iteration:  $S[i, 0] = S[0, j] = 0$ and for $i, j > 0$:

$$
S[i, j] =
\begin{cases}
S[i - 1, j - 1] + s(A[i], B[j]) \\
S[i - 1], j] - \delta \\
S[i, j - 1] - \delta \\
0
\end{cases}
$$

where $s(x, y) = 1$ if $x = y$ or $-\mu$ if $x \neq y$. Here $\mu$ is
the penalty of mismatch of two letters, and $\delta$ is the
penalty to insert or delete a character in either $A$ or $B$.

This is not a post about local alignments, so this extremely simple
definition of local alignment will suffice for our purposes. We can
even assume $\mu = \delta = 1$, which is possibly the simplest
alignment scheme. I should, however point out that the resulting score
matrix $S$ contains a lot of information about the evolutionary
relationship between $A$ and $B$. The highest-scoring location in $S$
is the most similar subsequence of the two, but other high-scoring
cells in the matrix can potentially identify subsequences that
translocated across the genome and diverged independently.  Possible
inversions can be found by aligning $A$ to the reverse-complement of
$B$ and vice-versa. By studying the alignment matrix, we can learn, in
great detail, what likely happened to certain genome subsequences, and
quantify how exactly they diverged between the two species.

**The limitation of alignments**: Quantifying similarities between
species through local alignments is not immediately obvious. The
highest alignment score (which, recall, is the most similar
subsequence) may not be representative of the global similarity, and
this becomes more true when comparing more distant species. We can use
other high scores to reconstruct the divergence of subsequences, but
it is not entirely obvious how to combine high scores in the matrix
into one global metric of similarity between sequences. Furthermore,
many genomes, including most mammalians, are billions of sequences
long. A $O(|A| |B|)$ algorithm is prohibitive both in time and memory.

My proposition to solve this problem is to take a step back and think
about the absolutely most naive way to compare two sequences: Consider
small sequences of length $k$ (for example, $k = 12$), which we will
call $k$-mers. What if we just counted $k$-mers in two sequences and
calculated the sum of squares of their differences? This is analogous
to measuring the similarity of two books based on their word
frequencies. Two dictionaries written by different people would be
very similar, and so would two math textbooks that would constantly
use words like "function", "variable", "number"< etc.

Note: we are shifting our thought process from "similarity"  -a number
that is higher the more similar two sequences are - to "distance".  A
distance function between two strings is a function $d$ such that, for
strings $A$ and $B$, $d(A, B) = 0$ if and only if $A = B$, $d(A,B) =
d(B,A)$ and for any three strings $A, B, C$, $d(A,B) \leq d(A,C) +
d(B,C)$. The edit distance is an example of a distance metric between
strings, and we can also prove that our $k$-mer distance also is,
since our $k$-mer counts induces a $4^k$-dimensional vector for which
Euclidean geometry properties apply.

This approach is not very rigorous, and it may not necessarily work,
but it is a start when thinking of a complex problem. Take a second to
think about how this could break: Could two highly similar sequences
have highly different $k$-mer counts?  Conversely, could two very
different sequences have similar $k$-mer counts? If you cannot
convince yourself that it can break easily, then the approach may have
some merit.

In fact, hopefully I will be able to show you that this extremely naive
approach leads to somewhat satisfiable results. The take-home message
should be that these simple approaches should often not be overlooked.

By the way, this is the most basic example of an [alignment-free
genome comparison
method](https://academic.oup.com/bib/article-abstract/15/3/343/182355).
Alignment-free sequence comparison research a beautiful field in and
of itself with various applications in metagenomics, virology and
population genetics.

# Problem formulation

We are given a set of $n$ species and their reference genome sequences
$S_1, \cdots, S_n$, which is the concatenation of all of their
chromosomes. We wish to obtain the following.

 * A $n \times n$ distance matrix between all pairs of species
 * A phylogenetic tree reconstructing the ancestral speciation events
   that led to the observed species, based on the similarity matrix.

For this post, we will use the following species: human, gorilla, mouse,
Chinese hamster, cat, alpaca and whale. Quick pause for some more
adorable pictures.

![our adorable species set](https://i.ibb.co/Rcwck5C/datasets.png)

Genomes of each species are given as [FASTA
files](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/). Those are
human-readable text files, where lines indicating chromosome names
start with the `>` character (e.g. `>chr1`), and below each such line
we see the genome sequence, with some line breaks after at most 80
characters (so the sequence can be read in a terminal).  For example,
below are the first lines of the cat genome (`felCat9.fa`):

```bash
$ head genomes/felCat9.fa
>chrA1
atcaggagatctagatgcctggagaggagtggagaaaacgggaaaccctc
ttATGggaagaggtaatatgtatttctccttcgaatataaaaaaagtaaa
aagaaggaaaacttaccaaattcacttatgagccattcattaccctgata
ccaaaaccagataaagccctccactaaaaccaaaactgcagcggcgcctt
gtgggctcggtcggttttactgtccaactcttaatttcagattaggaaat
aatcttgcggtgcatgggttcaagtcccacgttggaccctgccatgacag
tgtggggaatggctaggattctctctctccctgtctctctgcccctccct
cacttttttgtactctaaggaaagaaataaacatttaaaaaaatgttgaa
aattttttaaataaaactgcataccaatagccttgatgagtatgtatgcc
```

For our problem, we will download some FASTA files, inflate them,
count the k-mers and save the k-mer frequencies into a table, which
will then be turned into a distance matrix and a phylogenetic tree in
R.

# Obtaining the data

Put this in a text file and call it `genome_urls.txt`. These are the
URLs for the FASTA files of our species.

```
$ cat genome_urls.txt
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
https://hgdownload.soe.ucsc.edu/goldenPath/gorGor6/bigZips/gorGor6.fa.gz
https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
https://hgdownload.soe.ucsc.edu/goldenPath/criGriChoV2/bigZips/criGriChoV2.fa.gz
https://hgdownload.soe.ucsc.edu/goldenPath/felCat9/bigZips/felCat9.fa.gz
https://hgdownload.soe.ucsc.edu/goldenPath/vicPac2/bigZips/vicPac2.fa.gz
https://hgdownload.soe.ucsc.edu/goldenPath/balAcu1/bigZips/balAcu1.fa.gz
```

To download everything, on the same directory of your
`genome_urls.txt`, use `wget` on each line of the file. Then create a
`genomes` directory, unzip all the `.fa` files and move them to the
`genomes` folder as shown below.

```bash
for i in $(cat genome_urls.txt); do echo $i; wget ${i}; done
gunzip *.fa.gz;
mkdir genomes
mv *.fa genomes
```

Now that we have the FASTA files, we will read them in C++ and create
a $k$-mer frequency matrix.

# The C++ code

Let us first write a struct that will store $k$-mer frequencies for a
genome. We will use $k = 12$, and pre-allocate counts for all $4^{12}$
possible $k$-mers. We will assume counts fit into a 32-bit integer,
since $2^{32}$ is larger than any of the genomes we are using.

```cpp
struct KmerStats {
  KmerStats() { kmer_count = vector<uint32_t>(num_kmers, 0); }
  void count_kmers(const string &chrom);
  vector<uint32_t> kmer_count;
  static const uint32_t kmer_size = 12;
  static const uint32_t num_kmers = (2 << (2*kmer_size));
};
```

The function `count_kmers` takes a chromosome and increments the
vector `kmer_count`. We will encode $k$-mers as 24-bit numbers, with 2
bits per letter. Here we have `00 = A`, `01 = C`, `10 = G`, and `11 =
T`, so the number `0b110011001010010100011011` is the sequence
`TATAGGCCACGT`.

Here is the implementation of `count_kmers`:

```cpp
void
KmerStats::count_kmers(const string &chrom) {
  uint32_t mer = 0;

  auto itr(begin(chrom));
  const auto lim(end(chrom));
  for (size_t i = 0; itr != lim && i < kmer_size - 1; ++i, ++itr)
    shift_hash_key(*itr, mer);

  for (; itr != lim; ++itr) {
    shift_hash_key(*itr, mer);
    ++kmer_count[mer];
  }
}
```

The function `shift_hash_key` simply shifts a k-mer two bits, appends
a new letter to the end of it, then discards the first number by
applying an AND (`&`) operation to the number
`0b111111111111111111111111` (twenty-four 1s), which is `1 <<
(2*KmerStats::kmer_size) - 1`.

```cpp
inline void
shift_hash_key(const uint8_t c, uint32_t &mer) {
  static const uint32_t hash_mask = KmerStats::num_kmers - 1;
  mer = ((mer << 2) | encode_char[c]) & hash_mask;
}
```

Finally, we write a function to take a file name, open it, read the
chromosomes and process them onto a `KmerStats` object. Here we use a
pre-allocation size of 250 million for the `chrom` string, which is
larger than the largest chromosome we expect. This allows us to re-use
memory each time we read a chromosome.

```cpp
void
process_species(const string &file, KmerStats &v) {
  // get file size
  ifstream in(file);

  if (!in)
    throw runtime_error("cannot open file " + file);

  static const size_t RESERVE_SIZE = 250000000;
  string chrom;
  chrom.reserve(RESERVE_SIZE);

  string line;
  while (getline(in, line)) {
    if (line[0] != '>')
      copy(begin(line), end(line), back_inserter(chrom));
    else
      process_chrom(chrom, v);
  }

  process_chrom(chrom, v); // last chromosome after EOF
  in.close();
}
```

The `process_chrom` function simply encodes the string in two bits per
letter, replacing any non-ACGT character to a random character among
the two. This ensures that, when counting $k$-mers, we only see ACGTs
in our chromosome.

```cpp
inline void
process_chrom(string &chrom, KmerStats &v) {
  // makes sure letters are only ACGT:
  for (auto it(begin(chrom)); it != end(chrom); ++it)
    *it = ((encode_char[static_cast<uint8_t>(*it)] == 4) ?  "ACGT"[rand()%4] : *it);

  v.count_kmers(chrom);
  chrom.clear();
}
```

and the `encode_char` is a static look-up conversion between ASCII
characters and two-bit representations. We set 0 for the characters
`a` and `A`, 1 for `c` and `C`, 2 for `g` and `G` and 3 for `t` and
`T`. The rest we set to 4, and if we read a character that is a 4 in
our look-up, we replace with a nucleotide:

```cpp
static const uint8_t encode_char[256] = {
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //4
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //17
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //33
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //49
  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, //@,A-O
  4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //P-Z
  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, //`,a-o
  4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //p-z
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};
```

Finally, our `main` function takes two arguments, the first is a "file
of files", which is a two-column file showing the species name in the
first column and the path to its FASTA genome in the second. For
example, this is a file called `genomes.txt`:

```bash
$ cat genomes.txt
chinese_hamster genomes/criGriChoV2.fa
cat     genomes/felCat9.fa
gorilla genomes/gorGor6.fa
human   genomes/hg38.fa
mouse   genomes/mm39.fa
alpaca  genomes/vicPac2.fa
whale genomes/balAcu1.fa
```

And our `main` function simply reads this file and calls the functions
we created. When done, our output will be a table with $4^{12}$ rows and
$n$ columns, where $n$ is the number of species. The element in row
`i` and column `j` is the frequency of $k$-mer `i` (as its binary
representation) in species `j`. This is a data frame that we can read
into R:

```cpp
int
main(int argc, const char **argv) {
  if (argc != 2) {
    cout << "usage: ./phylo <input-species.txt>" << endl;
    return 0;
  }

  // ensures the k-mer size used fits in a 32-bit number
  static_assert(KmerStats::kmer_size <= 16);

  vector<string> species;
  vector<string> files;

  string tmp1, tmp2;
  ifstream in(argv[1]);
  while (in >> tmp1 >> tmp2) {
    species.push_back(tmp1);
    files.push_back(tmp2);
  }
  in.close();

  // print the headers: the species names, separated by tabs
  copy(begin(species), end(species), std::ostream_iterator<string>(cout, "\t"));
  cout << "\n";

  vector<KmerStats> v(species.size());

  omp_set_num_threads(8); // comment if OpenMP not used
#pragma omp parallel for
  for (size_t i = 0; i < species.size(); ++i) {
#ifdef VERBOSE

#pragma omp critical
    {
      cerr << "processing " << species[i] << "...\n";
    }
#endif
    process_species(files[i], v[i]);
  }

#ifdef VERBOSE
  cerr << "writing output\n";
#endif
  const size_t num_species = species.size();
  for (size_t i = 0; i < KmerStats::num_kmers; ++i) {
    for (size_t j = 0; j < num_species; ++j)
      printf("%d\t", v[j].kmer_count[i]);
    printf("\n");
  }
  return EXIT_SUCCESS;
}
```

## Compiling the code

We will create a `Makefile` to compile the program. If you do not have
OpenMP you can remove the `-fopenmp` flag below and the code will run
single-thread. It will use less memory and more time, but it should
still finish in a few minutes.

```make
all : phylo

phylo: src/phylo.cpp
	g++ -O3 -Wall -std=c++11 -o phylo src/phylo.cpp -fopenmp

clean:
	rm phylo
```

We can compile by simply running

```bash
make all
```

and run the program by writing

```bash
./phylo genomes.txt >kmer-counts.tsv
```

## Profiling the code

We can further profile using `/usr/bin/time` (GNU time) to see how
much time and memory it takes. Using 8 cores we get the following:

```bash
$ /usr/bin/time -v ./phylo genomes.txt >kmer-counts.tsv
processing chinese_hamster...
processing mouse...
processing whale...
processing gorilla...
processing cat...
processing alpaca...
processing human...
writing output
        Command being timed: "./phylo genomes.txt"
        User time (seconds): 525.60
        System time (seconds): 17.25
        Percent of CPU this job got: 524%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 1:43.54
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 1653768
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 0
        Minor (reclaiming a frame) page faults: 412911
        Voluntary context switches: 659062
        Involuntary context switches: 91400
        Swaps: 0
        File system inputs: 0
        File system outputs: 805360
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4096
        Exit status: 0
```

So we used a little under 2 minutes and 1.6 GB to create our table. Note
that the high memory use is because we used OpenMP to count the k-mers
of all 7 species in parallel, so all `KmerStats` objects were loaded.
Each `KmerStats` object allocates a vector of size $4^{12}$ with 32
bits, so each takes 64 MB. Then each thread uses an additional 250 MB
as pre-allocation for the chromosomes.

We can further profile the code using
[perf](https://www.man7.org/linux/man-pages/man1/perf.1.html).

```bash
perf record -v ./phylo genomes.txt >kmer-counts.tsv
perf report
```

# The R code

For the final part, we will construct a hierarchical clustering based
on the $n \times n$ distance matrix, defined as the sum of squares of
the differences between $k$-mer frequencies for any two species. We
will load the output of the C++ program into a matrix `x`, calculate
the pairwise distances using the `dist` function (we need to transpose
because `dist` is between all pairs of rows), then use `heatmap` and
`hclust` to make plots.

```r
> x <- read.table('kmer-counts.tsv', header = T, row.names=NULL)
> the.dist <- dist(t(x)) # distance between all pairs of species
> heatmap(the.dist) # heatmap of distance matrix
> plot(hclust(the.dist), hang = -1) # phylogenetic tree
```

This is the heatmap result:
![heatmap of distances](https://i.ibb.co/JpqbJ8R/heatmap.png)

And this is the resulting tree

![our resulting phylogenetic tree](https://i.ibb.co/Syhv1hb/phylo.png)

If you compare where the species lie in the [UCSC genome browser
tree](https://genome.ucsc.edu/cgi-bin/hgGateway), you will see that
our tree is consistent with theirs. Human and gorilla cluster
together, as do mouse and hamster. Whale, cat and alpaca cluster
together, with whale closer to alpaca. Hooray!

# Limitations of the method

This algorithm is admittedly naive. First, it treats the $k$-mer
frequencies in isolation. We have not accounted for co-occurences of
$k$-mers. Consider the $k$-mer `AAAAAAAAAAAA` and three genomes. In
genomes A and B, we have a large run of 1000 consecutive A, resulting
in 989 occurrences of this $k$-mer. Then, in genome C, the sequence
`AAAAAAAAAAAA` appears uniformly at random in the middle of larger
more complex sequences. The contribution of `AAAAAAAAAAAA` is
identical in all three genomes, but the correlation of $k$-mers in A
and B implies that they are more similar. We can incorporate
co-occurrences of sequences by modeling them as Markov chains, and
comparing the Markov chain parameters instead of the $k$-mers
directly. There is certainly a lot to improvement avenues to explore,
but our idea is a start.

# Follow-up questions

We can use a similar approach to guess, given a sequencing dataset
(e.g. Illumina, PacBio or Oxford Nanopore), which species it most
likely comes from. Can you think about how to do this? What
adjustments are required in this method?
