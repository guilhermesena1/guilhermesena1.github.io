---
layout: post
math: true
title:  "Wavefront alignments and the Myers edit distance algorithm"
date:   '2022-07-22'
categories: C++ alignment-algorithms
---

This is post explores the basis for the [wavefront alignment
algorithm](https://academic.oup.com/bioinformatics/article/37/4/456/5904262),
a 2020 paper that finds the alignment between two sequences $A$ and
$B$ in $O((|A| + |B|)s)$, where $s$ is the alignment "distance". We
are using "distance" instead of "score" because these algorithms only
make sense when the match score is 0 and the mismatch, indel and
possibly gap-open scores are non-negative.

Instead of exploring the generalized version, in this post we will
discuss the
[algorithm](https://link.springer.com/article/10.1007/BF01840446)
introduced by Prof. Gene Myers in 1986. Besides being the foundation
of the wavefront alignment, this algorithm is also the basis of the
`diff` tool in Linux, used to find the differences between two files.
I think this says volumes about the algorithm's relevance. Personally
I can only dream of writing an algorithm that will be adopted by a GNU
tool.

# Motivation

Many fields in computational biology revolve around finding
similarities and differences between extremely large sequences. In
most cases, these differences are formalized through approximate
string matching problem formulations.
[Local](https://www.sciencedirect.com/science/article/abs/pii/0022283681900875)
and
[global](https://www.sciencedirect.com/science/article/abs/pii/0022283670900574)
alignments are ubiquitous, so a lot of interest exists in creating
efficient implementations of alignment algorithms.

The Needleman-Wunsch and Smith-Waterman alignment algorithms are the
only ones guaranteed to find optimal global and local alignment,
respectively, of sequences $A$ and $B$ for arbitrary scoring schemes.
They are also $\Theta(|A| |B|)$, meaning that, whether $A$ and $B$ are
very similar or very different, the number of operations is the same:
We have to fill out the entire traceback matrix to find our answer,
either as the score at the bottom-right element in the matrix (global
alignment) or the maximum element in the matrix (local alignment).

In most applications in molecular biology, however, we only reach the
point of aligning two sequences when there is already some evidence
that the two sequences are likely to be similar. This is either
because we have prior knowledge of their similarity or because, prior
to comparison, we "filtered" a set of candidate sequences for which
there is some evidence of similarity. For example, in sequence
database search algorithms such as
[BLAST](https://www.sciencedirect.com/science/article/abs/pii/S0022283605803602)
or [Kraken](https://genomebiology.biomedcentral.com/articles/10.1186),
we use smaller subsequence in our query (called "seeds") and only
fully compare our sequence to "hits" in the database that match the
selected subsequence exactly. This means that the hits already have,
by definition, a significant similarity to the sequence. The same
procedure of "seeding" is done to map short sequences to a reference
genome in mapping algorithms like [minimap2](https://github.com/lh3/minimap2),
[abismal](https://github.com/smithlabcode/abismal) and literally every other
mapping tool ever written in the last 5 years or so.

For both the purposes of database search and read mapping (which in
some sense is also a database search), we want two types of alignment
functions. The first type  is a "fast" function that only computes the
alignment score and nothing else. We want these to perform as fast as
possible, and it is in our interest to minimize the number of
operations to obtain the score, after all we are only interested in
whichever sequence attains the highest score until we actually have to
report how the query and the best matching sequence differ. The second
"slow" function computes both the score and the exact operations (substitutions,
insertions and deletions) that are necessary to transform our query to
our best match. This information is used downstream to study *how*
they differ, but we only call this function to the highest scoring
candidate among all our hits.

It is very easy to verify (see
[appendix](#appendix-the-time-mappers-spend-aligning-reads)) that most
of the run time for read mapping algorithms is spent on alignment
algorithms to fully compare reads to hits (about 60-65% of the time
from the mappers we tested). In fact, we can even see that the most
time-consuming step are alignments that *do not* compute the
traceback, in other words, the false-positive alignments. We will
therefore focus on alignment algorithms that perform well for highly
similar sequences. We will begin with what is possibly the simplest
form of alignment of all: the *edit distance*.

# Edit distance problem formulation

The *edit distance* $D(A, B)$ between strings $A$ and $B$ is the
number of substitutions, deletions and insertions of characters in $A$
to make it identical to $B$. Obviously edit distance is a symmetric,
that is, $D(A,B) = D(B,A)$. In fact, edit distance is a metric (can
you think of a proof for the triangle inequality?). Similar to
alignments, edit distances can be computed in $O(|A| |B|)$ time. Let
$A = a_1, \dots a_m$ and $B = b_1, \dots, b_n$ be two strings and for
string $S$, $S[i]$ denote the prefix of $S$ with the first $i$
characters, with $S[0]$ being the empty string, then the edit distance
$D[i, j]$ of prefix $A[i]$ with prefix $B[j]$ is given by $D[i, 0] =
i$, $D[0, j] = j$ and $D[i, j] = D[i - 1, j - 1]$ if $A[i] = B[j]$,
or $D[i, j] = 1 + \mathrm{min}(D[i - 1, j - 1], D[i - 1, j],
D[i, j-1])$ otherwise. Those who have seen the [Smith-Waterman
recursion](https://guilhermesena1.github.io/posts/possibly-the-most-naive-phylogenetic-reconstruction-algorithm)
should find this recursion familiar, with the exception that we always
use $D[i-1, j-1]$ if there is a match, but sometimes we also take
diagonals for substitutions of characters.

**A quick note on our definition of edit distance:**
In the next section we will describe an algorithm to compute edit
distance based on the paper by Prof. Gene Myers. In the original 1986
formulation, he described that the longest common subsequence and the
shortest edit distance are dual problems. This is only true if the
only possible you allow are insertions and deletions. Here we are also
considering substitutions to be a single edit, which makes our
upcoming description of the algorithm a bit different to his. For
example, Myers would consider the edit distance of strings `ABA` and
`AAA` to be 2 (delete `B` and insert `A`), whereas we consider it to
be 1 (substitute `B` with `A`).

# The Myers algorithm for edit distance calculation

If we are simply interested in the edit distance between $A$ and $B$,
we only care about the value $D[m, n]$. In the recursion above, note
that we always take the diagonal when $A[i] = B[j]$. In other words,
the traceback of $D[m, n]$ is composed of vertical lines, horizontal
lines, and a series of diagonals where there are huge runs of matches
between substrings of $A$ and $B$. We can take advantage of these
large diagonals if we rethink the edit distance problem appropriately.

From the Myers paper:
>"In practical situations, it is usually the parameter $d$ that is small.
Programmers wish to know how they have altered a text file. Biologists
wish to know how one DNA stand has mutated into another."

The idea of the Myers algorithm is to take advantage of these runs of
diagonals by parametrizing the problem not by the prefix of the two
strings, but by the diagonals of the edit distance matrix $D$. First,
we denote $m + n - 1$ diagonals of the matrix by how far off they are
from the main diagonal, so $k = 0$ is the main diagonal, $k = 1$ is
one diagonal below the main, and $k = -1$ is one diagonal to the right
of the main. In other words, the cells in diagonal $k$ are the cells
$(x, y)$ where $x - y = k$.

We will show that, with this parametrization, we can find the edit
distance of $A$ and $B$ in $\Theta((m + n) d)$, where $m$ is the
length of $A$, $n$ is the length of $B$ and $d$ is the edit distance
between $A$ and $B$. This is much faster than the traditional $\Theta (mn)$
when $A$ and $B$ are expected to be similar, and not asymptotically
different when that is the case since $d \leq \mathrm{max}(m, n)$ when
substitutions only count as one edit and $d \leq m + n$ when we only
allow insertions and deletions (the Myers definition).

### The furthest-reaching point of diagonals

For diagonal $k$, let $V[d, k]$ be the *furthest reaching point* from
the origin of diagonal $k$ when we are allowed $d$ non-diagonal
changes (e.g. only insertions and deletions). In other words, $V$ is a
pair of coordinates $(x, y)$ where $x$ (or $y$) is maximum and $(x, y)$ can
be reached from $(k, 0)$ (if $k \geq 0)$ or $(0, -k)$ (if $k < 0$)
with $d$ edits. Then the edit distance of $A$ and $B$ is the minimum
value of $d$ for which $V[d, k] = (m, n)$ for some $k$.

The key idea for the algorithm is that, for any diagonal the furthest
reaching point with $d$ edits can be constructed easily from the
furthest reaching points with $d - 1$ edits.  First, if $d = 0$, $V[0,
k]$ is found by the largest series of matches between either the
suffix of $A$ starting at $k$ and $B$ or the suffix of $B$ starting at
$k$ and $A$ (depending on the sign of $k$).  For $d > 0$, we can
proceed as follows.  Start with the furthest reaching points of
$V[d-1, k-1]$, $V[d-1, k]$ and $V[d-1, k+1]$, both of which are one
cell away from the diagonal $k$. Then we have a starting point at $k$
either by going one cell to the left, one cell down, or one diagonal
cell from these two points, since these will increment one extra edit
to get to diagonal $k$. Pick whichever of the two is furthest in
diagonal $k$ to the origin, then go down diagonal $k$ on the series of
matches between $A$ and $B$ until the first mismatch is found. When we
find the point $(m, n)$ as the furthest reaching point for some pair
$(d, k)$, the edit distance is $d$.

The algorithm takes $O((m+n)d)$  time and, if we are only interested
in $d$, it takes $O(m+n)$ space. For $0 \leq d' \leq d$ (i.e. all
possible edit distances until our answer), we compute the furthest
reaching point when we are allowed $d'$ edits. When $d' = d$, we reach
the endpoint $(m, n)$ and stop.

# A short implementation of the algorithm

Here is a C++ program that computes the edit distance between
two strings passed on two lines of STDIN. The `edit_distance` function
implemented below also takes a third parameter `MAX_D`, and stops
trying to compute the edit distance if the answer goes above `MAX_D`.
This has useful applications in the problems we stated above. For
example, in read mapping, we will compare read to many sequences in
the reference genome to which we are mapping the reads. In the seeding
step, many of the retrieved sequences are false positives, and we want
to stop comparing as soon as we realize that we will not get a
satisfactory answer. This is normally done by setting a maximum
acceptable number of edits, which we can use to accelerate our
comparison.

A quick note about this code, which is easier to visualize with pen
and paper: Suppose you are at diagonal `k` and your horizontal offset
to the start of the diagonal is `x`. If you go one cell below to
diagonal `k+1` your horizontal offset to this new diagonal is `x+1`.
However, if you go one cell to the right to diagonal `k-1`, your
offset at this new diagonal will still be `x`. Finally, if you go
diagonally at diagonal `k`, you are still at diagonal `k`, but your
horizontal offset is now `x+1`. This is the rational for lines 30 to
38 in the code below. The last thing to mention is that we do not need
to keep a pair `(x, y)` as the furthest reaching point for diagonal
`k`. If we know `x`, then `y = x - k`, so everything is parametrized
by the horizontal offset `x`.

```cpp
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using std::string;
using std::vector;
using std::cin;
using std::cout;
using std::endl;

// we can use std::max from algorithm too, but this is faster
inline uint32_t
max32(const uint32_t a, const uint32_t b) { return (a>b) ? a: b; }

// edit distance function that stop if the distance is > MAX_D
uint32_t
edit_distance(const string &a, const string &b, const uint32_t MAX_D) {
  const uint32_t n = a.size();
  const uint32_t m = b.size();

  vector<uint32_t> V(2*MAX_D + 1, 0);
  vector<uint32_t> Vp(2*MAX_D + 1, 0);

  int32_t k = 0;
  uint32_t x = 0, y = 0;
  for (int32_t D = 0; D <= static_cast<int32_t>(MAX_D); ++D) {
    for (k = -D; k <= D; ++k) {

      // from the diagonal: increment one position
      x = (D == 0) ? 0 : Vp[MAX_D + k] + 1;

      // from above: we are also one point farther from
      // the diagonal above us when we move one position down
      x = (k == -D) ? x : max32(x, Vp[MAX_D + k - 1] + 1);

      // from the left: we are at the same diagonal position
      x = (k == D) ? x : max32(x, Vp[MAX_D + k + 1]);

      // extend the reaching point with run of matches
      for (y = x - k; x < n && y < m && a[x] == b[y]; ++x, ++y);
      if (x == n && y == m) return D;
      V[MAX_D + k] = x;
    }
    // copies the 2*D answers we found to use it in the next iteration
    if (D != static_cast<int32_t>(MAX_D))
      copy(begin(V) + MAX_D - D, begin(V) + MAX_D + D + 1, begin(Vp) + MAX_D - D);
  }

  // edit distance > MAX, reject comparison.
  return MAX_D + 1;
}

// general-purpose edit distance. If MAX_D not defined, then
// we use the maximum, which is |A| + |B|
uint32_t edit_distance(const string &a, const string &b) {
  return edit_distance(a.size() > b.size() ? a : b,
                       a.size() > b.size() ? b: a,
                       max32(a.size(), b.size()));
}

int main(int argc, const char **argv) {
  string a, b;
  cin >> a >> b;
  cout << "edit distance: " << edit_distance(a, b) << endl;
  return EXIT_SUCCESS;
}
```

# Benchmarking the algorithm

I ran the edit distance calculation to find the difference between
the human chromosome X in the human genome versions 37 (hg19) and 38
(hg38). For simplicity, I compared the first 1M characters of the two
assemblies, then the first 10M characters. I also removed all Ns from
both assemblies and converted letters to uppercase prior to
comparison.

This is the result for the first 1M characters:

```bash
$ /usr/bin/time -v ./wavefront <in-1m.txt
edit distance: 101430
  Command being timed: "./wavefront"
  Elapsed (wall clock) time (h:mm:ss or m:ss): 0:55.59
  Maximum resident set size (kbytes): 36512
  Exit status: 0
```

And this is the result for the first 10M characters:

```bash
$ /usr/bin/time -v ./wavefront <in-10m.txt
edit distance: 285326
  Command being timed: "./wavefront"
  Elapsed (wall clock) time (h:mm:ss or m:ss): 7:55.84
  Maximum resident set size (kbytes): 335364
  Exit status: 0
```

The difference is a bit steep. We have 10% on the first 1M characters
then less than 3% for 1M characters. Let us do 40M characters now!

```bash
edit distance: 330050
  Command being timed: "./wavefront"
  Elapsed (wall clock) time (h:mm:ss or m:ss): 10:53.26
  Maximum resident set size (kbytes): 1331500
  Exit status: 0
```

Looks like it's less than 1% now! The more we compare, the more
similar they become (though the number of edits is still increasing,
which is encouraging to confirm that the algorithm is working).

Do note that 40M reads would require $\approx 10^{15}$ operations in
the traditional edit distance algorithm proposed in the start of this
post, so we made good progress for similar sequences apparently :)

In the next post we will extend this idea of "furthest reaching point"
to affine local alignments to fully describe and implement the
beautiful algorithm shown by Marco-Sola and colleagues!

# Appendix: The time mappers spend aligning reads

To measure the alignment time, I took one million reads in a FASTQ
file and ran both [minimap2](https://github.com/lh3/minimap2) and
[abismal](https://github.com/smithlabcode/abismal) to map them to
`hg38` (it doesn't matter much how much error reads have, as long as
it's uniformly sampled from the genome).  Using `perf record <command>`,
we can store the time spent on each function in the
program on a summary file, then the `perf report` reads the file and
sorts the program
functions by time spent, finally summarizes it to us on screen.
It is an amazing profiling tool!

For `minimap2`, here is the result:

![results of perf for minimap2](https://i.ibb.co/JzrWGkp/minimap2.jpg)

For `abismal`, here is the result:

![results of perf for abismal](https://i.ibb.co/kqdXvZz/abismal.jpg)

 * In `minimap2`, we can see that 36.86% of the time is spent at function
`ksw_extd2_sse41`. This is the beautiful SSE implementation of
Smith-Waterman in the mapper, which we will cover in more detail in
another post.
 * In `abismal` we have 33.67% of the time on
`process_seeds` and 16.32% of the time on `AbismalAlign`. The
`process_seeds` step compares reads to hits using only Hamming
distance, and the `AbismalAlign::align<false>` function is the banded
Smith-Waterman algorithm (the "false" in the template means to only
calculate score, not traceback/CIGAR string).

In both cases, we can see that the bulk of time is not in collecting
seeds, finding them in the genome, ecoding or any other secondary
operation. Instead, the alignment itself comprises the majority of the
mapping effort, which greatly motivates fast implementations of the
algorithms we have available today.
