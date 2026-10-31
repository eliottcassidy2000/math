# Independent audit: Session124 finite tournament path-count spectrum

**PASS — bounded construction theorem only.** Every positive odd integer at most
999999, except 7 and 21, is attained by a finite tournament, conditional only on
the elementary constructor identities proved below (not a conjectural global
spectrum assumption). This audit does not prove that all larger allowed values
are attained. It does not independently re-prove the global exclusions 7,21.
No novelty claim is made.

## Independent computations

No audit executable imports code from `count_constructions`.

* `independent.cpp` enumerates every labelled tournament of orders 1 through 5,
  every external cut, every vertex permutation, and every split of that
  permutation into two paths. It gives exactly **273** states `(h,F,q,t)`.
* Its recursive word-tree traversal uses **no global state deduplication**.
  It visits **37,079,906** tree nodes, including over-ceiling leaves, and gives
  exactly **499,980** distinct values below 1000001. Its support equals all
  positive odd values in that range except 7,21 and the 18 stated misses.
  The emitted support is a 1000001-byte membership bitmap, not an external
  input or a presumed coverage interval.
* `check.py` independently computes coordinates for each of the 273 stored
  actual seed witnesses, checks equality with the independently enumerated
  universe, verifies all 14 factorizations and membership of every factor in
  the independent bitmap, and constructs four literal adjacency matrices.
* `adjacency_dp.cpp` uses a **backward endpoint** recurrence, unlike the producer's
  forward extension recurrence, on those full adjacency matrices. It obtains:

  | Vertices | Actual directed Hamiltonian paths |
  |---:|---:|
  |16|103061|
  |20|594793|
  |19|774073|
  |20|949777|

  No transfer-coordinate computation is used in these four counts. All graphs
  are checked to have exactly one directed edge per distinct vertex pair and
  no loops. Unsigned 64-bit arithmetic is safe: each partial path count is at
  most 20!, which is less than 2^64. All four full adjacency files are retained.
* `mod3.py` verifies the relevant matrix identities over F_3 directly, rather
  than repeating the producer's 81-input state loop.
* Independent traversal compiled with `-O0` and `-O2` has byte-identical output.
  Seed/factor/adjacency preparation under Python normal and `-O` is identical.
  No validation depends on `assert`.

The exact allowed misses before repairs are:

34869,85885,103061,116587,153451,245923,535405,594793,595543,
609873,618723,767135,774073,782103,824119,826373,868367,949777.

## Why the finite traversal covers all word lengths

For an actual parent tournament A and external vertex x, let h count its full
paths, F count paths of A+x, t count full paths ending at a vertex that beats x,
and q count ordered path covers (P,Q) with Q nonempty, P possibly empty, and
x beating Q's first vertex. A permutation plus a split represents each such
cover uniquely; the literal seed enumeration thus counts the stated objects.
Appending a parent sink v gives

* x→v: `(h,F,q,t) -> (h,F,2q+h,0)`;
* v→x: `(h,F,q,t) -> (h,F+q+h-t,2q,h)`.

For q, delete v from its path if that path is not the singleton Q=(v): placing
v at the end of P or at the end of Q gives the two old-cover terms. The
singleton Q case contributes h exactly when x→v. If x→v, v is an augmented
sink and F is unchanged. If v→x, paths ending at v contribute F-t; remaining
paths have v immediately before x and contribute q+h. These are exhaustive
and disjoint. Thus each generated state retains an actual constructor.

Always h≥1, q≥0, and 0≤t≤h. F and q are nondecreasing. If F>N no descendant
returns below N. If q>N, any subsequent zero increases F by at least q and
overshoots N; an all-one suffix leaves F unchanged. Therefore our DFS emits
the present F and stops at either bound without losing any value below N.
For q>0, each step at least doubles q, hence bounded branches are finite. For
q=0 a zero either is an exact self-loop or reaches a state where the next zero
is a self-loop; a one makes q=h>0. Deleting only exact self-loop edges preserves
all reachable states and makes the traversal finite. This is an all-finite-word
argument, not a depth cutoff assumption.

The producer's additional suppression of a one-child with q>N is also safe:
that child's F equals its parent's already emitted F, and its remaining
all-one descendants add no new value.

## Product repairs and the final finite conclusion

For an ordinal join A▷B, a directed Hamiltonian path cannot move from B back
to A. Consequently it is uniquely a full path of A followed by a full path of
B, so H(A▷B)=H(A)H(B). Iteration supplies all 14 verified factorizations using
actually attained factors. In particular the repair for 618723 uses
49·9·23·61, not an unattained factor7. These 14 repaired values and the four
full-adjacency witnesses partition the 18 misses exactly. Together with the
499980 direct values they supply all 499998 allowed positive odds ≤999999.
The stored support is not a separate adjacency certificate for every value:
existence follows inductively from the exact finite state/word constructors.

## Modulo-three restriction: correct scope

On `(F,q,g,h)` with g=h-t, a nonconstant triangle action C_j, j=1,2,
reduces modulo3 to `((3-j)g+2h,0,jh,0)`. Constant-cut triangle actions annihilate
all four coordinates. After any nonconstant triangle followed by a singleton,
q,g,h vanish modulo3; they remain zero through further singletons. Another
triangle therefore annihilates F as well. Three consecutive nonconstant
triangles also annihilate everything. Zero residues stay zero thereafter.

Hence any non-3-divisible output in this specified singleton/triangle language
uses at most two triangle blocks, both nonconstant, adjacent if there are two.
This holds for arbitrary initial integer vectors, not merely sampled seeds.
It is a necessary language restriction, **not** an impossibility statement
about unrestricted tournaments, and not a proof that the language misses
infinitely many primes. The singleton q>N pruning is **not** transferred to
arbitrary blocks: a p=0 block can multiply F while ignoring q. The four repair
certificates need no expanded-language completeness claim.

## Reproduction

From this package's root, run:

```sh
python3 run_checks.py --work-dir /tmp/spectrum-check
```

The root README gives the individual commands. Generated executables and the
1,000,001-byte support bitmap remain in that explicit work directory; neither
is a publication input. `SHA256SUMS` freezes package sources, actual witness
data, and fresh receipts. The normal/optimized Python modes and O0/O2 C++
modes produce identical relevant output and data. The two independent
algorithms have not been replaced by shared constructor imports.

Additional hostile control: reversing seed edge (0,1) separately in each
of the four full adjacency witnesses changes all four counts and is rejected
by the DP checker. Counts become 101247,769467,754839,759405 respectively.
See `mutation.py` and `../receipts/mutation.normal.out`. Malformed adjacency
input is also explicitly rejected; no validation depends on Python asserts.
