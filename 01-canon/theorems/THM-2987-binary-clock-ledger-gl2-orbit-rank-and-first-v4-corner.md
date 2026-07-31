---
id: THM-2987
title: "Binary-clock ledger GL2 orbit rank and first V4 corner"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For the three
  GL(2,F_2)-conjugate padding gauges of the THM-2976
  homogeneous forced-parity ledger, dyadic checkpoints vanish. At every
  non-dyadic checkpoint their span has rank one for depth zero, rank two and
  an exact affine V4 carrier for positive power-of-two depth, and rank three
  otherwise. Consequently a THM-2976 corner-timed odd-unit ladder has ranks
  1,2,3,3,... by 2-adic clock depth: V4 occurs only at the first nontrivial
  clock. The physical padding line has only a C2 stabilizer, and no
  entrywise-nonnegative integral projective order-three lift exists; no
  physical cell, capacity, owner, or phase transport is claimed.
source: codex-binary-clock-ledger-gl2-orbit-2026-07-30
depends_on:
  - THM-2976-binary-clock-parity-for-critical-run-deficit-ledgers
related:
  - THM-2595-modular-v4-affine-lift-dichotomy-and-six-vertex-tournament-no-go
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
script: 04-computation/amm12592_binary_clock_gl2_orbit_rank_thm2987.py
output: 05-knowledge/results/amm12592_binary_clock_gl2_orbit_rank_thm2987.out
script_sha256: 45f542fbb5b6c0e2fab924cc7d0e64458ec4dfc6aa0ed9dae6b92e3cc5ca9dd2
output_sha256: 8f5fca5d1fa10aa93bf832408d389f4209780e07edb56fa87e9e55c8080a434f
hash_basis: LF-normalized bytes
audit: >
  thm2987-hostile-audit-2026-07-30 (independent Frobenius/Lucas and
  set-polynomial rank derivation; THM-2976 corner, pointed-line, and K4 audit;
  PGL(2,Z) no-nonnegative-projective-order-three proof and bounded hostile
  census; normal/-O/stored LF-hash replay and docs: ACCEPT)
---

# THM-2987 -- binary-clock ledger GL2 orbit rank and first V4 corner

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2976 proves the binary clock of one physically distinguished homogeneous
padding gauge. The abstract `V4 semidirect S3` anatomy is already canonical
in THM-2595/2606. The new statement here is narrower and ledger-specific:
exactly when do the three conjugate padding gauges of the THM-2976 polynomial
form a `V4`, and does that structure persist along its corner-timed ladder?

The answer is a sharp trichotomy. The `V4` occurs at one clock depth only.
Moreover, the order-three permutation of the three gauges is not a symmetry
of the pointed physical ledger and admits no entrywise-nonnegative integral
lift remaining projectively order three.

## 1. The homogeneous factorization

Work in the domain

```text
R=F_2[p,q],                 l=p+q.
```

Put `N=M+1` and `d=d_M`. Homogenizing THM-2976 equation (B) gives

```text
B_l=l^(N+d)+(p^N+q^N)l^d
   =S_N l^d,                                            (1)

S_N=p^N+q^N+l^N.                                       (2)
```

The subscript on `B_l` records that the physical THM-2976 construction pads
with the distinguished line `l=p+q`. Define its two conjugate gauges

```text
B_p=S_N p^d,             B_q=S_N q^d.                  (3)
```

The group `G=GL(2,F_2)` permutes the three nonzero linear forms

```text
L={p,q,l}.                                               (4)
```

Consequently `S_N`, the sum of their `N`th powers, is `G`-invariant, and
`G` permutes the three polynomials in `(3),(1)`.

The checkpoint factor has the exact Frobenius converse

```text
S_N=0             iff             N is a power of two. (5)
```

The forward direction for dyadic `N` is Frobenius. Conversely, if `N` is
not a power of two, its binary expansion has at least two nonzero bits.
Choose a non-leading set bit `2^a`. Lucas gives

```text
binom(N,2^a)=1 mod 2,
```

so `(p+q)^N` has a nonzero interior term and cannot equal `p^N+q^N`.
This proves `(5)` and recovers THM-2976 T1 in the homogeneous gauge.

## 2. Exact orbit-span trichotomy

Let

```text
W_(N,d)=span_F2{B_p,B_q,B_l}.                            (6)
```

Then:

```text
N dyadic:                 dim W_(N,d)=0 for every d;

N non-dyadic, d=0:       dim W_(N,d)=1;

N non-dyadic, d=2^a>0:   dim W_(N,d)=2 and
                          {0,B_p,B_q,B_l}=W_(N,d)=V4;

N non-dyadic, d>0 not
a power of two:           dim W_(N,d)=3.                 (7)
```

*Proof.* The dyadic case is `(5)`. Suppose `N` is non-dyadic, so `S_N` is
nonzero. Multiplication by `S_N` is injective because `R` is a domain. It is
therefore enough to compute the rank of

```text
p^d, q^d, l^d.                                          (8)
```

For `d=0` all three equal one. Let `d>0`. No singleton relation is possible.
No pair relation is possible either: `p^d` and `q^d` are distinct monomials,
while `l^d` contains both endpoint monomials `p^d` and `q^d`. Thus the only
remaining possible nonzero relation is the three-term one, and

```text
p^d+q^d+l^d=0        iff        d is a power of two     (9)
```

by the same Frobenius/Lucas argument as `(5)`. This proves every rank in
`(7)`. In the rank-two case the three nonzero, pairwise distinct vectors sum
to zero, so adjoining zero gives the full two-dimensional vector group
`V4`. QED.

For `d=2^a`, Frobenius makes the identification with the natural module
literal:

```text
(u p+v q)^d=u p^d+v q^d,             u,v in F_2.        (10)
```

Multiplication by the invariant nonzero `S_N` then identifies the natural
`F_2^2` with the carrier in `(7)`. Hence `G=S3` acts faithfully on its three
nonzero elements. For other positive `d`, the rank-three object is the
three-letter permutation module in characteristic two; it is not an affine
`V4`.

## 3. Corner timing meets V4 only once

At a non-dyadic THM-2976 corner-timed level, T3 says

```text
v=v_2(N),                  d=2^v-1.                     (11)
```

Substitution into `(7)` gives

```text
v=0:        d=0,                         rank one;
v=1:        d=1,                         rank two, V4;
v>=2:       d=2^v-1 odd and greater 1,   rank three.     (12)
```

Thus `V4` is not a persistent carrier of the corner clock. It occurs only at
the first nontrivial `2`-adic clock `v=1`.

More explicitly, THM-2976 C5 gives for every odd `J>=3` the corner ladder

```text
N=J 2^v,                   d=2^v-1.                     (13)
```

Its orbit ranks are exactly

```text
1,2,3,3,3,...                                                (14)
```

as `v=0,1,2,...`. The `J=1` classical ladder is different: every `N` is
dyadic and the entire forced-parity carrier vanishes by `(5)`. In particular,
the repeated corner timing at rate `1/3` does not carry one same `V4` through
all scales.

## 4. Pointed-line guardrail

The physical polynomial is `B_l`, not the unpointed set in `(6)`. Its line
stabilizer in `G` has order two:

```text
Stab_G(l)=C2.                                             (15)
```

Either order-three element cycles

```text
p -> q -> l -> p                                         (16)
```

and therefore cycles three gauge-conjugate ledgers. It does not act on one
same pointed physical ledger. The algebraic orbit preserves homogeneous
degree, forced parity up to variable relabelling, and the ranks in `(7)`. It
forgets:

```text
the distinguished padding line l;
which exponent is the THM-2976 q/corner coordinate;
the integer signs and magnitudes before reduction mod two;
row cells, binomial capacities, and admissible scheme choices;
owner, phase, and endpoint data.                          (17)
```

This is the same kind of quotient discipline emphasized independently by
THM-2984: retaining a primitive direction can make a finite cell test strict,
but that theorem supplies no transport between the three padding gauges here.

## 5. No positive integral order-three lift

The marked order-three matrix over `F_2` may be written

```text
Tbar=[0 1]
     [1 1].                                               (18)
```

The entrywise-nonnegative integral matrix with the same entries is not of
projective order three:

```text
Tplus^2=[1 1],                 det(Tplus)=-1,             (19)
        [1 2]
```

and its powers follow the unbounded Fibonacci recurrence. A genuine signed
integral lift is

```text
Tsigned=[0 -1],                Tsigned^3=I,               (20)
        [1 -1]
```

so the missing sign is load-bearing.

In fact there is no entrywise-nonnegative `A in GL(2,Z)` whose class in
`PGL(2,Z)` has order three. If `[A]^3=1`, then `A^3=lambda I`. Taking
determinants gives `det A=1` and `lambda=+-1`. In the nontrivial case:

```text
A^3= I  forces characteristic polynomial X^2+X+1,
         hence trace(A)=-1, impossible for A>=0;

A^3=-I  forces characteristic polynomial X^2-X+1,
         hence trace(A)=1.                               (21)
```

For a nonnegative integral two-by-two matrix of trace one, the two diagonal
entries have product zero, so

```text
det A=ad-bc<=0,
```

contradicting `det A=1`. This proves the no-lift statement.

Therefore the mod-two order-three gauge cycle cannot be promoted to a
positive integral basis/cell/capacity transport. Any physical use would need
new labelled data carrying all three gauges into one common owner/capacity
chart. No such sidecar is constructed here.

## 6. K4 and ordered-basis context only

In the `V4` case, each nonzero translation of

```text
X={0,B_p,B_q,B_l}                                        (22)
```

pairs the four points into two edges. These are the three perfect matchings
of the abstract `K4`, exactly the generic frame in THM-2606. Likewise the six
ordered bases of `X` carry the regular `S3` action already implicit in
THM-2595/2606 and visible in other six-sheet canon.

These are related interpretations, not new physical identifications. There
is no canonical map here to four quartic roots, a Keller source, an LRC
runner packet, or a THM-2966 scheme. In particular, the abstract zero in
`(22)` does not select a physical quartic section or repair the pointed-line
loss in `(17)`.

## 7. Exact evidence and audited scope

The dependency-free companion represents homogeneous `F_2[p,q]`
polynomials as q-exponent bitsets. It independently checks:

```text
the THM-2976 factorization (1);
S_N=0 iff N is dyadic;
all 8,320 cases 1<=N<=128, 0<=d<=64 of (7);
the corner classification through N=256;
every odd J=3,5,...,31 through clock depth eight;
all six GL(2,F_2) matrices and the C2 line stabilizer;
the three K4 perfect matchings;
the positive and signed integral lifts (19)--(20).         (23)
```

Run

```text
python 04-computation/amm12592_binary_clock_gl2_orbit_rank_thm2987.py
python -O 04-computation/amm12592_binary_clock_gl2_orbit_rank_thm2987.py
```

Normal and optimized transcripts LF-normalize byte-for-byte to the declared
stored output, and both declared hashes match.

```text
PROVED CONCLUSIONS:
  exact ledger-orbit rank trichotomy;
  exact V4 carrier at positive dyadic depth;
  first-nontrivial-clock-only V4 corner;
  pointed padding-line C2 stabilizer;
  no entrywise-nonnegative integral projective order-three lift.

NOT CLAIMED:
  a physical symmetry of the THM-2966 scheme;
  preservation of cells, capacities, owners, endpoints, or phase;
  a new extractor construction or lower bound;
  an LRC, Keller, quartic, or tournament consequence.      (24)
```

**QED.**
