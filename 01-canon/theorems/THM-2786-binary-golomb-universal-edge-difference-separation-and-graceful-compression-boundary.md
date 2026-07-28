---
id: THM-2786
title: "Binary Golomb universal edge-difference separation and graceful compression boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The labels
  2^i-1 give every finite simple graph pairwise distinct
  absolute edge differences, with an exact 2-adic decoder for both
  endpoints.  More generally every radix q works.  A self-contained
  quadratic prime construction compresses universal all-pair separation
  below 8n^2.  Conversely any m-edge distinct-difference labeling has span
  at least m, and any labeling separating all vertex pairs has span at least
  binom(n,2).  Therefore a graceful n-vertex tree must deliberately allow
  collisions involving at least one nonedge difference: universal Golomb
  separation over-solves the problem and cannot be compressed to the
  graceful span n-1.  The
  already-proved rooted Nullstellensatz theorem THM-2765 does exploit this
  distinction to reach span 3n-5; the remaining debt is the exact
  factor-three-to-one compression.  This does not prove Graceful Tree.
source: root/binary-golomb-graceful-boundary-2026-07-28
external_input: >
  Bertrand's postulate, used only for the corollary that an odd prime
  n<=p<2n exists for n>=2 and hence the displayed prime ruler has span
  <8n^2.  The n=1 case is the trivial one-point ruler of span zero.
depends_on:
  - THM-2761-graph-edge-sum-discriminant-codegree-factorization-and-graceful-sign-gauge
  - THM-2765-rooted-nullstellensatz-linear-range-distinct-edge-labeling
related:
  - THM-2783-weighted-long-wall-binary-null-avoidance-and-ternary-state-reconstruction
script: 04-computation/binary_golomb_graceful_boundary_thm2786.py
output: 05-knowledge/results/binary_golomb_graceful_boundary_thm2786.out
script_sha256: 61ca79bae3dc9317833f2310745a7270d017a4880ffbba1d5ce64960628b0f9a
output_sha256: 3210f27f8a9d4083e295d212cfcedeebd64fe6e48b9d344959a9d4b9b9bc9f1e
hash_basis: LF-normalized bytes
---

# THM-2786 -- a universal ruler solves the wrong graceful problem

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

There is an elementary collision-free labeling for every finite simple graph:

```text
vertex i |-> 2^i-1.                                      (1)
```

It is far stronger than graceful labeling, because it separates the
differences of every pair of vertices, not only the edges.  That strength is
exactly why it cannot reach the graceful range.  The theorem below turns
this familiar powers-of-two observation into a source/target map, compresses
it to quadratic range, and proves that any linear compression must use the
tree's missing edges.  THM-2765 already does so and reaches `3n-5`; the
theorem here explains why that change of category is load-bearing.

## 1. Exact binary endpoint decoder

Label ordered vertices `0,...,n-1` by

```text
a_i=2^i-1.                                                (2)
```

For `i<j`,

```text
a_j-a_i=2^i(2^(j-i)-1).                                  (3)
```

The second factor is odd.  Therefore

```text
i=v_2(a_j-a_i),
j=i+log_2((a_j-a_i)/2^i+1).                              (4)
```

Thus the positive difference determines the ordered endpoint pair exactly.
In particular, every simple graph on these vertices has pairwise distinct
absolute edge differences.

The same proof works for every integer radix `q>=2`:

```text
a_i=q^i-1,
a_j-a_i=q^i(q^(j-i)-1),                                 (5)
```

where the second factor is not divisible by `q`.  The largest power of `q`
dividing the difference recovers `i`, and the quotient recovers `j-i`.
Among these positional rulers, binary has the smallest span

```text
2^(n-1)-1.                                               (6)
```

This is the pair-difference counterpart of THM-2783's binary radix wall.
It is not yet a bounded graceful labeling.

## 2. Quadratic universal compression by a prime ruler

The exponential span in `(6)` is not intrinsic to all-pair separation.  Let
`p` be any odd prime with `p>=n`, and for `0<=i<n` let

```text
r_i = the least residue of i^2 modulo p,
b_i = 2pi+r_i.                                           (7)
```

The `b_i` are strictly increasing, because a consecutive linear increment
`2p` dominates any residue drop.  They form a Sidon set.  Suppose

```text
b_i+b_j=b_k+b_l,             i<=j, k<=l.                 (8)
```

Then

```text
2p(i+j-k-l)=r_k+r_l-r_i-r_j.                             (9)
```

The right side has absolute value below `2p`, so both sides vanish:

```text
i+j=k+l,                  r_i+r_j=r_k+r_l.              (10)
```

Modulo `p`, `(10)` and equality of square sums give

```text
ij=kl mod p.                                             (11)
```

Since `p` is odd, the two pairs are the roots of the same quadratic over
`F_p`.  All indices lie in `{0,...,p-1}`, hence

```text
{i,j}={k,l}.                                             (12)
```

This proves Sidon sum uniqueness.  If two positive differences were equal,

```text
b_j-b_i=b_l-b_k,
```

then `b_j+b_k=b_l+b_i`; applying `(12)` gives the same endpoint pair.
Therefore `(7)` is a universal Golomb ruler as well.

Its span satisfies

```text
b_(n-1)<2p^2.                                            (13)
```

For `n>=2`, Bertrand's postulate permits an odd prime `n<=p<2n`, with the
trivial choice `p=3` at `n=2`.  Hence every `n`-vertex finite simple graph
has one explicit integral all-pair-separating labeling of span

```text
<8n^2.                                                   (14)
```

The formula and proof `(7)--(13)` are self-contained; Bertrand is used only
for the displayed uniform constant in `(14)`.  For `n=1`, the one-point
ruler has span zero.

## 3. The exact graceful compression obstruction

Let a graph have `m` edges, translate any integral labeling so its minimum is
zero, and let its span be `L`.  If its `m` positive edge differences are
distinct, then they lie in `{1,...,L}`, so

```text
L>=m.                                                     (15)
```

Equality in `(15)` is exactly the bounded-difference part of gracefulness:
the differences must be `{1,...,m}`.  For an injective labeling into
`{0,...,m}`, THM-2761 proves that this is equivalent to nonvanishing of the
squared-difference discriminant.

Now impose the stronger, graph-independent requirement that **all**
`binom(n,2)` vertex-pair differences be distinct.  The same count gives

```text
L>=binom(n,2).                                            (16)
```

For a tree, `m=n-1`.  When `n>=3`,

```text
binom(n,2)>n-1.                                          (17)
```

Therefore no universal all-pair Golomb ruler—binary, prime-quadratic, or
otherwise—can be compressed to the graceful tree span while retaining
all-pair separation.  A graceful-tree argument must use the tree incidence
and permit collisions involving at least one nonedge difference.

This is the sharp conceptual boundary:

```text
universal ruler:   separates K_n, needs at least binom(n,2) span;
graceful tree:     separates only E(T), seeks exactly n-1 span.          (18)
```

The missing operation is not generic collision avoidance.  It is
tree-specific selective collision: preserve the `n-1` edge differences
while folding differences involving the much larger nonedge family.

THM-2765 proves that this distinction already has force.  Its rooted
coefficient-one construction labels every tree in `{0,...,3n-5}` while
keeping only the tree-edge absolute differences distinct.  For `n>=6`,

```text
binom(n,2)>3n-5,                                           (18a)
```

so the counting bound `(16)` proves that every such linear-range route must
allow a collision involving at least one nonedge difference.  This need not
be a nonedge--nonedge collision: the graceful `P3` labels `(0,2,1)` have
edge differences two and one, while their sole nonedge difference one
collides with the latter edge difference.  The remaining graceful problem
is therefore not the quadratic-to-linear transition supplied by `(7)`.
It is the sharper compression

```text
rooted 1+2 channel range 3n-5  ->  graceful range n-1,    (18b)
```

while retaining vertex injection and the two mirror signs of THM-2765.

## 4. Small complete-graph controls

For a complete graph, `(16)` is the elementary counting floor for a Golomb
ruler.  Exact exhaustive search gives the first optimal spans

```text
n=2,3,4,5,6:          1,3,6,11,17.                      (19)
```

The counting floor `binom(n,2)` is attained through `n=4` and already fails
at `n=5,6`.  First witnesses are

```text
(0,1),
(0,1,3),
(0,1,4,6),
(0,1,4,9,11),
(0,1,4,10,12,17).                                      (20)
```

Only the finite cases in `(19)--(20)` are asserted here.  They are controls
showing that even complete-graph compression has structure beyond counting;
they are not used in the all-`n` theorem.

## 5. Consequences and stopping boundary

For every tree, `(1)` and `(7)` make THM-2761's absolute-difference
discriminant nonzero.  This supplies a universal rational/integral point
away from the graceful collision hypersurfaces, first at exponential and
then at quadratic height.  It proves that the discriminant polynomial is
not identically zero without constructing a point in the small graceful
box.  It is not a better range theorem than THM-2765; its value is the
all-pairs decoder and the impossibility boundary `(16)--(18b)`.

The transfer ledger is:

```text
source:       a radix or prime-quadratic Sidon ruler;
target:       a nonzero tree edge-difference discriminant;
preserved:    every edge difference, indeed every vertex-pair difference;
destroyed
 by desired
 compression: identities involving at least one nonedge may and must collide;
needed
 sidecar:      improve THM-2765's rooted 1+2 channel load without losing
               vertex injection or either mirror sign;
 cheapest test: lower one rooted exponent layer while preserving the
                coefficient-one monomial and the range bound.
                                                                    (21)
```

No tournament orientation, modular action, LRC owner/current, Keller map, or
Jacobian information is present.

## 6. Exact verification

Run

```bash
python 04-computation/binary_golomb_graceful_boundary_thm2786.py
python -O 04-computation/binary_golomb_graceful_boundary_thm2786.py
```

The exact companion uses explicit exceptions, integer arithmetic, and no
truth-bearing Python assertions.  It checks all binary endpoint decoders
through `n=13`; all `20,825` prime-ruler pairs through `n=50`; the sum-Sidon,
difference, ordering, prime-window, and span claims; and exhaustively proves
the optimal values `(19)` by testing every smaller normalized ruler.  Normal
and optimized runs byte-match the stored transcript.

```text
PROVED HERE:              binary and general-radix endpoint decoder;
                          universal finite-simple-graph edge-difference
                          separation;
                          explicit prime-quadratic Sidon/Golomb ruler;
                          span <2p^2 and, by Bertrand, <8n^2;
                          edge-count and all-pair span lower bounds;
                          necessity of a nonedge-involving collision for
                          graceful span;
                          THM-2765 linear-range comparison and remaining
                          factor-three compression boundary;
                          complete-graph optima through six vertices.

NOT PROVED:               an improvement to THM-2765's 3n-5 tree range;
                          a graceful labeling for any previously open tree;
                          the Graceful Tree Conjecture;
                          an optimal general Golomb-ruler bound;
                          a tournament, modular-group, Keller, or LRC result;
                          JC(2), DC(2), or LRC(14).                       (22)
```

QED.
