---
id: THM-4113
title: "Maximal noncrossing half-Kempe atlas and the 2--3 partition law"
status: >
  PROVED + VERIFIED-EXACT. The nonredundant planar partial matchings produced
  by Algorithm B.1.2 of arXiv:2608.22870v1 are exactly the maximal
  noncrossing partial matchings on a cyclic boundary. Rooting the boundary
  gives a bijection with either one noncrossing partition into blocks of
  sizes 2 and 3, or one root singleton together with an ordered pair of such
  partitions. If u marks singleton half-chains, their bivariate generating
  function is M=C+uxC^2 where C=1+x^2C^2+ux^3C^3. This gives an explicit
  coefficient formula, a direct output-sensitive generator, the total
  sequence 1,1,1,3,4,10,20,42,98,..., and the exact shifted signed
  A321197 identity. At boundary size five, coarse pair compatibility is the
  Petersen/A4-root graph of THM-261: planarity removes a crossing C5 from
  its 15 edges and leaves precisely the 10 half-Kempe atlas states.
source: codex-snark-apex-260822870-20260825
depends_on: []
related:
  - THM-261-petersen-root-orthogonality
  - THM-2290-context-selected-colored-pair-kernel-is-hafnian-complete
  - THM-2672-slope-seven-carry-nerve-exact-eleven-simplex-and-root-zero-cap
  - THM-3990-componentwise-harmonic-obstruction-and-repair-quotient
  - THM-4104-selected-order-eleven-strong-ear-solid-interval
external:
  - "Yuta Inoue, Ken-ichi Kawarabayashi, Ritarou Matsuo, Atsuyuki Miyashita, Bojan Mohar, and Tomohiro Sonobe, Three-edge-coloring apex cubic graphs, arXiv:2608.22870v1, Appendix B.1."
  - "OEIS A001005, noncrossing 2--3-block partitions."
  - "OEIS A321197, A-sequence of a Riordan matrix."
script: 04-computation/snark_half_kempe_atlas_thm4113.py
output: 05-knowledge/results/snark_half_kempe_atlas_thm4113.out
script_sha256: 87419c3c74782ffa9a220098306e665124ffc11a242be2bf20b1efa8ad8a214a
output_sha256: ce0157c903e68ec842ed1fb72e84017c2530a2866ba2721a91ceaf67186fde48
hash_basis: raw working-tree bytes
audit: >
  PASS. A literal set-valued transcription of Algorithm B.1.2 agrees through
  boundary size 12 with an independent Motzkin recursion followed by direct
  maximality testing, and through size 14 with the closed coefficient law.
  Explicit maps in both directions audit 25,203 maximal matching encodings,
  11,223 type-C partitions, and 13,980 rooted type-X objects. The five-point
  Petersen split is checked independently. The script contains no assert
  statements; normal and optimized runs reproduce the frozen transcript.
---

# THM-4113 -- maximal noncrossing half-Kempe atlas

**PROVED + VERIFIED-EXACT.** The apex-cubic proof introduces singleton
half-Kempe chains because a two-colour alternating path can terminate at a
degree-two vertex. Its Appendix B.1.2 generates a nonredundant family of
planar partial matchings on each boundary. That family has an exact closed
description and an algebraic counting law.

The closest proved repo mechanism is THM-2290: a matching calculation is
exact only while its endpoint selectors and contraction family are retained.
The canonical hostile here occurs already at five boundary positions:
disjointness alone admits 15 pair-of-pair states, whereas cyclic planarity
admits only 10. The least-used relevant sidecar is THM-261's identification
of that coarse disjointness graph with the Petersen graph and the positive
roots of `A_4`.

## 1. The Appendix algorithm generates precisely the maximal atlas

Put `n` labelled points in cyclic order. A **partial matching** is a set of
pairwise vertex-disjoint unordered pairs. It is **noncrossing** when no two
pairs have cyclic order

```text
a<c<b<d.                                                     (1)
```

It is **maximal noncrossing** when no pair of currently unmatched points can
be added without producing a crossing. Unmatched points are the singleton
half-chains.

Let `H_n` be the output of Algorithm B.1.2 in arXiv:2608.22870v1. The paper
defines it for `n>=1`; put `H_0={emptyset}` as a harmless empty-boundary
convention.

> **Theorem 1 (algorithmic identification).** For every `n>=0`, `H_n` is
> exactly the set of maximal noncrossing partial matchings on the cyclically
> ordered set `[n]`.

### Proof

The algorithm processes the number `q` of pairs in decreasing order. For a
chosen `2q`-subset it generates every Catalan noncrossing perfect matching on
that subset. It rejects a candidate if two unmatched points are cyclically
adjacent. It then calls the candidate redundant exactly when adding some pair
of unmatched points produces an earlier generated candidate.

Every maximal noncrossing partial matching has no adjacent unmatched points:
their boundary edge could otherwise be added. Hence every maximal matching
passes the algorithm's first filter. If a candidate with `q` pairs can be
extended, the extension has `q+1` pairs and was already put in the algorithm's
set `A`, because pair counts are processed downward. Conversely, membership
of such an extension in `A` says precisely that the added pair does not cross
the old pairs. Therefore the algorithm rejects exactly the extendable
candidates and retains exactly the maximal ones. The `n=0,1` cases are
immediate. **QED.**

This also explains why the algorithm's cyclic-adjacency test is necessary
for pruning candidates but carries no extra condition in the final theorem:
maximality already forces it.

## 2. Root-region bijection with 2--3-block partitions

Distinguish the boundary gap between `n-1` and `0`, thereby choosing a root
region of the chord diagram. The noncrossing chords cut the disk into
regions. Each nonroot region has a unique **parent chord**, namely the chord
separating it from the root region.

Two unmatched boundary points can be joined without crossing the matching if
and only if they lie in the same region. Consequently a noncrossing partial
matching is maximal exactly when

```text
each complementary region contains at most one unmatched point.             (2)
```

For every parent chord `{a,b}` do the following:

```text
no unmatched point in its child region:       make the block {a,b};
one unmatched point c in its child region:    make the block {a,b,c}.        (3)
```

The blocks in `(3)` form a noncrossing partition into blocks of size two or
three. There are two cases.

1. If the root region contains no unmatched point, `(3)` partitions all `n`
   boundary points.
2. If the root region contains its unique unmatched point `r`, no chord can
   straddle `r`. Deleting `r` therefore leaves an ordered pair of independent
   2--3-block partitions, one before and one after `r` in the rooted order.

The construction is invertible. In each 2-block, match its two points. In a
3-block `{a<b<c}`, match `{a,c}` and leave `b` unmatched. In the second case
also leave the distinguished root point unmatched. Noncrossing of the block
partition makes the resulting pairs noncrossing, and every unmatched point
occupies a different region, so `(2)` makes the matching maximal. These two
maps undo each other.

> **Theorem 2 (root-region bijection).** Maximal noncrossing partial
> matchings are in bijection with the disjoint union
>
> ```text
> P_23  disjoint_union  (one root point) * P_23 * P_23,       (4)
> ```
>
> where `P_23` is the class of noncrossing partitions into blocks of sizes
> two and three. A 3-block records one singleton half-chain; the distinguished
> root records one additional singleton in the second summand.

The boundary root is a gauge, not discarded information. Moving it changes
which singleton, if any, is designated as the outer one, while leaving the
underlying cyclic semi-matching unchanged.

## 3. Bivariate generating function and explicit coefficients

Let `C(x,u)` count `P_23`, where `x` marks boundary points and `u` marks
3-blocks. Looking at the block containing the first point gives the unique
formal-power-series solution with constant term one of

```text
C = 1 + x^2 C^2 + u x^3 C^3.                               (5)
```

Let `h_(n,s)` be the number of members of `H_n` having exactly `s`
singleton half-chains, and put

```text
M(x,u)=sum_(n,s>=0) h_(n,s) x^n u^s.                        (6)
```

The bijection `(4)` gives

```text
boxed: M(x,u)=C(x,u)+u x C(x,u)^2.                          (7)
```

In particular `h_(n,s)=0` unless `n` and `s` have the same parity. Lagrange
inversion makes every coefficient explicit. Put

```text
a=(n-3s)/2,              b=(n+2-3s)/2.                      (8)
```

Interpret a term as zero when its displayed integer is negative or
nonintegral, or when a factorial has a negative argument. Then

```text
boxed:
h_(n,s)
 = n!/(a! s! (n+1-a-s)!)
   + 2 n!/(b! (s-1)! (n+2-b-s)!).                           (9)
```

The first term is the `C` case and the second is the rooted `uxC^2` case.
For example,

```text
h_(10,0)=42,       h_(10,2)=420,       h_(10,4)=30.          (10)
```

On setting `u=1`, write `C(x)=C(x,1)` and `M(x)=M(x,1)`. If `Y=xC`, then

```text
Y=x(1+Y^2+Y^3),        M=(Y+Y^2)/x.                          (11)
```

Therefore

```text
[x^n]M
 = 1/(n+1) ([t^n]+2[t^(n-1)])(1+t^2+t^3)^(n+1).            (12)
```

Eliminating `Y` gives the cubic equation

```text
x^4 M^3-2x^2 M^2+(1+x^2)M-(1+x)=0.                         (13)
```

The resulting totals, beginning with the empty convention, are

```text
n:       0  1  2  3  4   5   6   7   8    9   10    11
|H_n|:   1, 1, 1, 3, 4, 10, 20, 42, 98, 210, 492, 1122.    (14)
```

The coefficients of `C` are OEIS A001005, the noncrossing 2--3-block
partition numbers. There is also an exact, previously unrecorded-in-this-repo
Riordan identity. If `a_k` denotes OEIS A321197, then

```text
boxed: |H_n|=(-1)^(n+3) a_(n+2),          n>=0.             (15)
```

Indeed, let `f(z)=z/(1+z^2-z^3)` and let `A(z)=z/f^[-1](z)` be the
A-sequence series of the associated Riordan matrix. Equation `(11)` says
`f(-Y)=-x`, so `A(-x)=1/C(x)`. Dividing `(5)` by `C` gives

```text
A(-x)=1/C=1-x^2(C+xC^2)=1-x^2M,                            (16)
```

and coefficient comparison is `(15)`. Thus the sequence match is a theorem,
not an identification from finitely many terms.

## 4. The five-boundary Petersen firewall

At `n=5`, a maximal state has two pairs and one singleton, and there are ten
such states. Regard each possible pair as a vertex. Two pair-vertices are
coarsely compatible exactly when their endpoint sets are disjoint. By
THM-261 this ten-vertex compatibility graph is

```text
K(5,2) = the Petersen graph = A_4 positive-root orthogonality. (17)
```

The Petersen graph has 15 edges. The five edges whose two chords cross form
a `C_5` on the five diagonal pair-vertices. Removing that crossing pentagon
leaves exactly the ten edges corresponding to `H_5`; the surviving graph is
connected and has degree profile

```text
1^5, 3^5.                                                    (18)
```

Hence disjointness/root orthogonality is a real but incomplete carrier:

```text
source:       cyclic boundary semi-matching;
coarse map:   endpoint pairs -> vertices of K(5,2);
preserved:    pair identity and endpoint disjointness;
destroyed:    cyclic crossing, boundary face, colour transposition;
sidecar:      cyclic order + face label + chosen colour pair;
hostile:      15 Petersen edges versus 10 planar states.      (19)
```

This is also the exact reason not to replace the boundary system by a
tournament. The native binary observable is compatibility of two chain
endpoints, and even that graph needs a rotation sidecar. An arbitrary
orientation adds a gauge while failing to restore the missing `C_5` test.

## 5. Direct generation and multi-boundary products

Equation `(4)` gives a direct generator with no generate-and-prune table:

1. recursively generate a 2-block or 3-block containing the first point;
2. recurse independently in the two or three intervening intervals;
3. generate either one such partition, or select one root singleton and an
   ordered pair of partitions;
4. apply the inverse map following Theorem 2.

With standard persistent output representation this is output-sensitive up
to the linear cost of writing each state. It can replace the redundancy loop
of Algorithm B.1.2 without changing the atlas.

The paper forms a multi-boundary semi-matching by concatenating the choices
on its boundary faces. Therefore boundaries of sizes `n_1,...,n_r` have
exactly

```text
product_(i=1)^r |H_(n_i)|                                  (20)
```

nonredundant topology choices. The singleton profile is the coefficient
profile of the product of the fixed-size polynomials
`sum_s h_(n_i,s)u^s`. This makes boundary count a literal state-space
coordinate, rather than an informal complexity warning.

## 6. Relation to earlier repo mechanisms

1. **THM-2290 (hafnian carrier).** Both results retain endpoint selectors
   and the matching contraction family. THM-2290 sums over unrestricted
   perfect pairings; this theorem uses maximal **noncrossing partial**
   pairings with singleton half-chains. Replacing one contraction family by
   the other is invalid. The five-point Petersen hostile is the cheapest
   witness.
2. **THM-3990 (componentwise repair).** A boundary region is a genuine
   component coordinate: maximality permits at most one singleton per
   region. This is a combinatorial analogy to retaining componentwise repair
   classes, not an application of a Laplacian quotient.
3. **THM-2672 and multi-boundary nerves.** Splitting the root singleton into
   two interval partitions is an exact instance of boundary components
   becoming part of the carrier. It does not identify the arithmetic carry
   nerve with a Kempe chain.
4. **THM-4104 (selected insertion response).** A full boundary response is
   compositional while an aggregate mean is not. Here `(7)` retains the
   complete singleton-refined response; the scalar total `|H_n|` alone does
   not decide which ring coloring a topology can repair.

## 7. Exact audit and scope

Reproduce with

```bash
python3 -B 04-computation/snark_half_kempe_atlas_thm4113.py
python3 -B -O 04-computation/snark_half_kempe_atlas_thm4113.py
```

The first path literally generates the paper's descending Catalan candidate
sets. The second path independently generates every Motzkin partial matching
through `n=12` and directly tests whether any unmatched pair is addable. The
third path evaluates `(9)` through `n=14`. Explicit encoders and decoders test
both sides of `(4)`, and a separate graph calculation checks `(17)--(18)`.
The normal and optimized transcripts agree with the frozen output.

This theorem audits and improves one finite generator in the v1 proof. It
does not independently reproduce the paper's 915-configuration discharging
check, its 109,501-island reducibility check, or the theorem that every
2-connected apex cubic graph is three-edge-colourable. It does not count
ring **colorings**, which add ternary labels and semi-consistency closure, and
it does not make the mutable external verifier repositories a pinned proof
artifact. The paper's global theorem remains a cited preprint claim under its
own source and computation audits.
