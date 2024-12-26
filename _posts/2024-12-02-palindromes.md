---
layout: post
math: true
title:  "Palindromes palindromes palindromes!"
date:   '2024-12-02'
categories: C++ strings palindromes
---

This post explores the solution of a problem I proposed to my
undergrad students when I was a TA in a computational molecular
biology course. It involves palindromes so readers who suffer from
[aibohphobia](https://en.wiktionary.org/wiki/aibohphobia) should
probably sit this one out.

I made this problem up during a discussion session and, when I did, I
actually had no idea how to solve it. We had just covered various
exact string matching algorithms, so they knew I was expecting a
linear solution.  In computational biology, given the length of the
sequences we work with, we often cannot afford to do worse than
$$O(n)$$.

**Problem**: Given a DNA sequence of length $$n$$, what is the minimum number of
letters that must be appended to the end of the sequence to make it a
DNA palindrome? What is the resulting palindrome?

Note that a "biological" palindrome is not the same as a lexicographic
palindrome. In english, a palindrome is a word or sentence that reads
the same when reversed. Examples include "ana", "racecar", "a man, a
plan, a canal: Panama", "Did Hannah see bees? Hannah did" and, of
course, the 224-word palindromic poem [Dammit, I'm
mad!](https://www.nku.edu/~longa/classes/mat385_resources/docs/Palindromic-Poem.pdf).
In biology it is almost the same but with a small difference. A
palindromic DNA sequence is one that reads the same when
*reverse-complemented*. Examples include GAATTC, ACGT, AGGCCT,
AGGCCCT, and so on. Funnily enough none of these are lexicographic
palindromes and, conversely, simple sequences like AAAA are
lexicographic palindromes but not biological palindromes.

# Examples

1. `ATTGCTT` $$\to$$ `ATTGCTTAAGCAAT`. This is a "worst" case scenario
and an example of why there is always a solution. If you
reverse-complement the DNA sequence and append it to itself, it
becomes a palindrome.

2. `ATTGCA` $$\to$$ `ATTGCAAT`. This example is an almost palindrome,
and we only needed to add two bases to make it palindromic. Note that
we could also have appended its reverse-complement, but the resulting
sequence `ATTGCATGCAAT` is much longer than the optimal answer.

3. `ACGT` $$\to$$ `ACGT`. This one is already a palindrome, so no need
to add bases.

4. `AAAA` $$\to$$ `AAAATTTT`

# Why care?

Palindromic DNA sequences are recurring themes in molecular biology.  Consensus
sequences that enzymes bind to are often palindromic.  Examples include the [CpG
dinucleotide](https://en.wikipedia.org/wiki/CpG_site) that is the primary target
for methylation, and the GAATTC motif that is bound by the [EcoRI restriction
enzyme](https://en.wikipedia.org/wiki/EcoRI) in *E. coli*.  Palindromic RNA
transcripts can also form secondary structures by binding ends to form stable
geometries. In the enzymatic case, this is not by chance, as a palindromic DNA
suggest that it matters not in which strand of the DNA the enzyme binds to,
making it easier for it to find the motif of interest.

Please note, however, that this specific problem does not have
practical applications to my knowledge. We will not be building
synthetic cut sites from *in vivo* DNA or anything like that! :)

Read more about palindromes in genetics
[here](https://en.wikipedia.org/wiki/Palindromic_sequence).

# A linear solution

We should look at some examples that deviate from the worst and best
case scenarios to build an intuition for an algorithm. Let's look at
example 2 above: `ATTGCA`. If we reverse-complement the sequence, we
get `TGCAAT`. Here we notice that the sequence and its
reverse-complement partially overlap, and it just so happens that
merging the sequence and its rc to the partial overlap results in our
answer:

```
ATTGCA
  TGCAAT
========
ATTGCAAT

```

So all we need to do is to find the largest prefix of the sequence's
reverse-complement that matches a suffix of the original sequence.
This resembles string matching algorithms very closely. We can think
of prefix-suffix overlaps as follows. Call the original sequence $S$
and the reverse-complement $T$. Think of $T$ as a "needle" we want to
find in the "haystack" $S$, just like in classic string matching
problems. The intermediate steps of looking for $T$ in $S$ often keep
track of the largest match for a prefix in $T$ to any location in $S$.
If we don't stop the process and let $S$ be searched for to its end,
we will end up with a partial solution of an exact match between a
prefix of $T$ and a suffix of $S$. We can write this as a simplified
version of the [KMP
algorithm](https://en.wikipedia.org/wiki/Knuth%E2%80%93Morris%E2%80%93Pratt_algorithm),
so the major steps of our solution will be:

1. Reverse-complement the input string $S$ to make string $T$
2. Build the [partial match
table](https://en.wikipedia.org/wiki/Knuth%E2%80%93Morris%E2%80%93Pratt_algorithm#%22Partial_match%22_table_(also_known_as_%22failure_function%22))
for $T$
3. Look for $T$ in $S$ using the KMP procedure, but don't stop until
we reach the end of $S$
4. Append the suffix of $T$ after the partial prefix-suffix match to
$S$ and return it.

# C++ code

I tried to write a solution with the fewest lines. Obviously there is
room for optimization, such as reverse-complementing using a look-up
table instead of conditions. I believe we can do the matching
procedure without having to make a copy of the string to
reverse-complement but I think this would make the solution hard to
read.


```cpp
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

using namespace std; // don't do this ever

// this would be better as a look-up ASCII array
char comp(const char c) {
  if (c == 'A') return 'T';
  if (c == 'C') return 'G';
  if (c == 'G') return 'C';
  if (c == 'T') return 'A';

  // we will assume the input is always well-formatted and we will not
  // reach this part
  return 'N';
}

// low-effort revcomp function
void
revcomp(string &s) {
  reverse(begin(s), end(s));
  transform(begin(s), end(s), begin(s), comp);
}

string
palindromize(const string &s) {
  const size_t len = s.size();

  // s is already a palindrome if 0 or 1 characters
  if (len <= 1) return s;

  string t = s;
  revcomp(t);

  // build longest-prefix-suffix table for t
  vector<int> lps(len, 0);
  lps[0] = 0;
  for (size_t i = 1, j = 0; i < len;) {
    if (t[i] == t[j]) lps[i++] = ++j;
    else if (j != 0) j = lps[j - 1];
    else { j = 0; ++i; }

  }

  // j stores the longest prefix of T in S at any time
  size_t j = 0;
  for (size_t i = 0; i < len; ) {
    if (s[i] == t[j]) { ++i; ++j; }
    else if (j != 0) j = lps[j - 1];
    else ++i;
  }

  return s + t.substr(j);
}


int main() {
  string s;
  cin >> s;
  s = palindromize(s);
  cout << s << endl;

  // sanity check: Do we produce the same sequence if we revcomp it?
  // revcomp(s);
  // cout << s << endl;
  return 0;
}

```

# Results
```
| input              | output                |
|--------------------|-----------------------|
|AAAA                |AAAATTTT               |
|ATTATATTGC          |ATTATATTGCAATATAAT     |
|CGCGCCGGCGGC        |CGCGCCGGCGGCCGCCGGCGCG |
|AATTCCGG            |AATTCCGGAATT           |
|AACCTTGG            |AACCTTGGCCAAGGTT       |
|AACCGGTT            |AACCGGTT               |
|AACCGGGTT           |AACCGGGTTAACCCGGTT     |

```
I think it works! :)


# Follow-up question

1. (Easy) Is there a sequence that is simultaneously a lexicographic palindrome
and a biological one? If so, show it, if not, prove it.

2. (Hard) What if we could add letters to any part of the sequence, not just
append to the end? What would be the shortest palindrome we could make?
