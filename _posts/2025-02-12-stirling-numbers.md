---
layout: post
math: true
title:  "Stirling numbers and a useful combinatorics problem for bioinformatics"
date:   '2025-02-12'
categories: combinatorics bioinformatics R
---

# Balls, boxes bioinformatics

It never ceases to amaze me that sometimes in combinatorics the
difference between an elementary school problem and a lifelong career
in math is one additional sentence in a problem statement.

Consider the following combinatorics problem, which we will call the
"box problem".

<span style="color:#DDDDDD">
There are $$k$$ empty boxes. At each turn, a marble is dropped into one
of the $$k$$ boxes with equal probability. What is the minimum number
$$C$$ of marbles that need to be dropped to ensure that, with
probability at least $$\alpha$$, there are at least $$r$$ marbles in all
boxes?
</span>

For example if $$\alpha=0.99$$, $$k=2$$ and $$r=1$$. There are two
boxes and we want at least one marble in each box. Obviously it must
be that $$C \geq k \times r = 2$$ given the problem constraints. For
$C=2$, with 50% probability, both balls are dropped in one of the
boxes. With $$C=3$$, we only fail if all 3 marbles are dropped in
either box, which happens with probability $$2/(2^3)$$. Continuing
this logic we need the smallest $$C$$ for which $$2/(2^C) \leq 1 -
\alpha = 0.05$$, which happens for $$C \geq 6$$. Try doing this
exercise for $$k=3$$ and $$r=1$$ and your brain may start to hurt a
bit. Try for $$k=3$$ and $$r=2$$ and you may, like I almost did, toss
your notebook out the window trying to find a generalizable formula.

## The "box problem" in bioinformatics

This problem arises when we want to reconstruct the sequence of a
locus that may contain multiple copies in a targeted assay, such as a
retrotransposon or a gene with variable copy number. A common approach
to do this shown in the figure below. The genome is represented in
lilac and the four copies of a locus in green. You can cut the DNA using
restriction enzymes or guide RNAs in the ends of your sequence of
interest, represented as the light green corners of the locus. Each
copy can be slightly different from each other, as represented by the
orange and purple SNPs whose combination is unique to each copy.

When we sequence the targeted data that is cut and isolated, each read
can come from any copy of the locus with equal probability. Reads can
contain sequencing error, so in order to reconstruct the sequence of
each copy, it is often important to have at least $$r$$ reads for each
of the $$k$$ copies of the locus. We can assume that the error is
small enough that we will always cover sufficient variation to assign
reads to a given copy unambiguously. Knowing the number $$k$$ of
copies we can expect (or at least the maximum possible value of
$$k$$), what is the number $$C$$ of reads we need to guarantee, with
high probability $$\alpha$$, that we will have at least $$r$$ reads
per locus we can use to reconstruct their sequence? 

![multiple copies of a locus present on a reference
genome](https://i.ibb.co/pBH7Wb7V/coverage.png)

## How is this not a trivial problem?!

[Fermat's last
theorem](https://en.wikipedia.org/wiki/Fermat%27s_Last_Theorem) is the
living proof that a simple problem statement does not equal a simple
solution. I have no intuition as to why this problem is so difficult
besides the fact that I was trying to parametrize as a function of
$$k$$, $$r$$ and $$C$$ but there was no efficient strategy that did
not double count some set configuration.

In what follows, let us forget the $$\alpha$$ part for a second and
focus on trying to count how many ways we can distribute the marbles
in the $$k$$ boxes with at least $$r$$ marbles per box. If we succeed,
we can calculate both the numerator and the denominator for the
probability.

# Set partitions and numbers named after smart people

We can think of each of the $$C$$
marbles as elements of a set, and distributing the marbles across
boxes is a way of partitioning the set into at most $$k$$ subsets,
some of which can be potentially the empty set.  Intuitively this
analogy makes sense, also because the number of partitions that
distributes the elements "somewhat evenly" is much larger than the
number of partitions that, for example, puts all marbles in the first
box. This matches our intuition that we should always have at least
some marbles in most boxes given a sufficient number of marbles

## Bell number

The number of ways to partition a set with $n$-set $$\{1, 2, \dots,
n\}$$ is the [Bell number](https://en.wikipedia.org/wiki/Bell_number)
$$B_n$$.The first Bell numbers are $$(B_1, B_2, B_3, B_4, B_5) = (1,
2, 5, 15, 52)$$.  For example, $$B_2 = 2$$ because the only ways to
partition the set $$\{1, 2\}$$ are $$\{\{1\}, \{2\}\}$$ and $$\{\{1,
2\}\}$$. 

There is
no closed formula for it, but it obeys the following recurrence:

$$
B_{n+1} = \sum_{m=0}^{n} {n \choose m} B_m.
$$

This stems from the fact that, for any of the $B_{n+1}$ partitions,
removing the partition containing the first element leaves a partition
of a smaller number $m$ of elements. There are $${n \choose m}$$ ways
to decide which elemens are part of the smaller set (all the other
elements are part of the removed set) and $$B_{m}$$ ways of
partitioning them.


## Stirling numbers of the second kind

An analogous number to the bell number is [Stirling's number of the
second
kind](https://en.wikipedia.org/wiki/Stirling_numbers_of_the_second_kind),
$$S(n, k)$$ which counts the number of ways to partition the $$n$$-set
into $$k$$ non-empty subsets. It immediately follows form this
definition that $$B_n = \sum_{k=0}^{n} S(n, k)$$. Furthermore, a very
simple recursion can be used to find $$S(n, k)$$. Obviously $$S(n, 1)
= 1$$ for every $$n$$

$$S(n+1,k) = kS(n, k) + S(n, k-1).$$

The rationale for this recurrence is: When adding a new number to an
$$n-set$$, you can either add a new partition to any of the $$S(n,
k-1)$$ existing $$k-1$$ partitions or add the new number to an
existing $$k$$-partition, for which there are $$k$$ possible sets to
add it to. This procedure also allows all partitions to be generated.

Therefore, if $$r = 1$$, our problem is solved. The number of ways to
add $$C$$ marbles to $$k$$ boxes such that none are empty is
$$N=S(C,k)$$. The total number of ways to partition the marbles into the
boxes is $$D=\sum_{i=1}^{k} S(C,k)$$ and our answer is the minimum
value of $$C$$ for which $$N/D \geq \alpha$$.

Unfortunately the general cass for $$r > 1$$ significantly cranks up
the difficulty of this problem.

## Restricted & associated Stirling numbers of the second kind



## Calculating restricted Stirling numbers in R

