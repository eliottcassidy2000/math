---
id: THM-2212
title: "Cofactor, anchor, and ruler pruning of rank-eleven states"
status: >
  PROVED, relative to LRCUpTo13 through THM-763 and THM-2208. Put
  Q=91^6. For each rank-eleven bounded relation matrix B, two coordinate
  anchors identify ker(B) with Q^2, and every primitive row with
  M(v)<=1/14 has a positive anchor pair of sum at most Q^2. Hence a fixed
  B supports fewer than Q^4/2 candidates, independently of the enormous
  transverse-relation height. For a supplied rank-twelve state, maximal
  cofactors reduce exact weak/strict LRC testing to 91 fixed-dimensional
  integer ruler systems, preceded by exact gcd, divisibility, annulus,
  order-chamber, and sparse-circuit sieves. A weak-only survivor has a
  balanced signed mod-14 equality skeleton. This is a finite exact pruning
  architecture, not an executed enumeration or a proof of LRC(14).
source: codex-2026-07-24-projective-state-pruning
depends_on:
  - THM-763-strict-finite-height-for-tight-lrc-instances
  - THM-2052-finite-height-forces-high-rank-bounded-relation-code
  - THM-2208-rank-eleven-positive-cone-projective-slope-collapse
related:
  - THM-668-pair-sum-ruler-witness-structure
  - THM-1002-quantitative-top-gap-margin-deathstar-S56
  - THM-2186-exact-octagon-needle-profile-and-mobius-h-drift
  - THM-2188-finite-phase-bank-and-pairwise-overlap-no-go
---

# THM-2212 -- cofactor, anchor, and ruler pruning

Put

```text
Q=91^6,                    Q^2=91^12.                 (1)
```

Let `B in Z^(11 x 13)` have rank eleven and let

```text
S=ker_Q(B).                                             (2)
```

The motivating states are the support-at-most-three, height-`Q` matrices
supplied by THM-2052, but the first reduction below uses only the rank.

## 1. An `H`-independent per-plane candidate bound

Choose coordinate labels `p,q` such that

```text
pi_(p,q):S -> Q^2,       x |-> (x_p,x_q)              (3)
```

is an isomorphism. Such a pair exists because `S` has dimension two.
THM-763, applied to thirteen speeds and the cited LRCUpTo13, says that
every primitive positive row with

```text
M(v)=max_t min_i ||v_i t||<=1/14
```

satisfies

```text
sum_i v_i<=91^12=Q^2.                                 (4)
```

In particular its positive anchor pair obeys

```text
v_p+v_q<=Q^2.                                         (5)
```

There are exactly

```text
Q^2(Q^2-1)/2 < Q^4/2                                 (6)
```

ordered positive integer pairs satisfying (5), and (3) reconstructs at
most one rational vector in `S` from each pair. Thus each fixed `B`
supports fewer than `Q^4/2` possible primitive non-strict rows. For each
pair, integrality, positivity, distinctness, primitivity, and (4) are exact
rational/integer tests.

This removes the height

```text
H=78*182^13
```

from the per-plane loop. It is strictly sharper than enumerating all
height-`H` transverse relations `a` from THM-2208.

There is a useful order preprocessor. In the anchor coordinates, write

```text
v_i=ell_i(x,y)                                        (7)
```

as rational linear forms. If some `ell_i-ell_j` vanishes identically on
`S`, the plane has no repetition-free candidate. Otherwise the at most
`binom(13,2)=78` equality slopes cut the projectivized positive cone into
at most `79` open order chambers. On each chamber the minimum, maximum,
and second-largest labels are fixed. The necessary zero-Haar annulus

```text
13 v_min<=v_max<=13 v_second-max                      (8)
```

then becomes two rational linear inequalities before any lattice point is
enumerated. The left inequality follows from the `v_min+v_max`, `p=1`
ruler: strict reverse inequality gives an open safe interval. The right
inequality is the contrapositive of THM-1002's quantitative top-gap bound.
The order chamber is only a routing object; the linear forms in (7) retain
the indispensable magnitudes.

## 2. Cofactors give the exact projective candidate

For a supplied THM-2208 state, put

```text
C=[B;a] in Z^(12 x 13),          rank(C)=12,           (9)
```

and let `z` be its signed maximal-cofactor kernel vector. Reject the state
unless all `z_i` are nonzero with one common sign and their absolute values
are pairwise distinct. After the global sign change, write

```text
g=gcd(z_1,...,z_13),            v=z/g.                (10)
```

Then `v` is the unique primitive positive candidate of THM-2208. Raw
cofactors are suitable for determinant hashing, but the divisibility and
mod-14 sieves below must use `v`: an artificial common factor can conceal
the phase obstruction.

## 3. The exact 91-ruler terminal

For every `1<=i<=j<=13`, put

```text
q=v_i+v_j.                                             (11)
```

THM-668 gives the exact equivalence

```text
M(v)>=1/14
 iff some ruler q and some p mod q satisfy
     14|p v_l|_q>=q             for every l,           (12)
```

where `|r|_q=min(r mod q,-r mod q)`. The same test can be run on raw
cofactors with `q=z_i+z_j`, by dilation invariance, but (10) makes the
arithmetic smaller.

For a fixed ruler, weak feasibility is exactly the integer system

```text
1<=p<=q-1,
L_w<=v_l p-q k_l<=q-L_w             (l=1,...,13),
0<=k_l<=v_l-1,
L_w=ceil(q/14).                                        (13)
```

Strict feasibility replaces the lower band by

```text
L_s=floor(q/14)+1.                                    (14)
```

Thus exact LRC testing costs at most `91` integer feasibility systems in
the fixed fourteen variables `(p,k_1,...,k_13)`. It admits a
fixed-dimensional ILP algorithm (for example, Lenstra's), rather than an
explicit `O(q)` scan of every multiplier.

## 4. A gcd/live-ruler sieve before feasibility

For the weak test put

```text
m_w=ceil(q/14)-1,
B_l={p in Z/q:|v_l p|_q<=m_w},
g_l=gcd(v_l,q).                                        (15)
```

Then

```text
|B_l|=g_l(2 floor(m_w/g_l)+1).                        (16)
```

For the strict test use `m_s=floor(q/14)`. Merge labels with

```text
v_l=+/-v_k mod q,                                     (17)
```

because they have the same bad-multiplier set. Every remaining bad class
contains `p=0`, so

```text
sum_classes (|B_class|-1)<q-1                         (18)
```

certifies a nonzero live multiplier. If `q|v_l` for some `l`, the ruler is
dead immediately. Equations (15)--(18) require only gcd and congruence
arithmetic and are exactly THM-668's C1 ledger specialized to the cofactor
candidate.

## 5. Zero-Haar and equality-skeleton filters

Every positive distinct zero-Haar candidate must satisfy all of:

1. the annulus (8);
2. at least one speed is divisible by each `d=2,...,13`, since otherwise
   `t=1/d` is strictly safe;
3. every live ruler with `14` not dividing `q` is impossible, because its
   weak integer band equals its strict band and gives an open safe interval.

For an actual counterexample `M(v)<1/14`, the divisibility list extends
through `d=14` and every one of the 91 rulers must be weak-dead.

Suppose instead a zero-Haar row has a weak ruler witness but no strict one.
Then necessarily

```text
14|q.                                                  (19)
```

At `t=p/q`, define the active boundary sets

```text
A_+={l:v_l p= q/14 mod q},
A_-={l:v_l p=-q/14 mod q}.                            (20)
```

Both sets are nonempty. If, for example, all active runners had only the
`+` sign, then moving a short distance in the direction in which those
sawteeth increase would preserve every active inequality; the inactive
runners have positive clearance. This produces a one-sided safe interval,
contrary to zero Haar. The same argument handles the other sign.

The exact residual record

```text
(q,p,A_+,A_-, all nonactive residue gaps)             (21)
```

is therefore a balanced signed equality skeleton, the one-dimensional
tangent carrier required by THM-2186. Here "balanced" means that both sign
sets are nonempty, not that their cardinalities are equal.

The sparse rows of `B` prune it further. If an integer relation `b.v=0`
has support wholly contained in `A_+ union A_-`, assign
`sigma_l=+1` on `A_+` and `-1` on `A_-`. Multiplying the relation by
`t=p/q` gives the necessary congruence

```text
sum_l b_l sigma_l=0 mod 14.                            (22)
```

Thus every support-two or support-three row supplies a tiny exact sign
table before any tangent continuation is attempted.

## 6. Complete pruning order and no-go

A faithful exact implementation may proceed:

```text
canonical rational rowspace B
 -> at most 79 order chambers
 -> positive anchor pairs satisfying (5) and (8)
 -> reconstruct and test integrality/distinctness/gcd/divisibility
 -> 91 gcd/live-ruler ledgers
 -> weak/strict integer systems only on survivors
 -> balanced equality skeleton plus sparse sign congruences. (23)
```

Alternatively, for an already supplied `(B,a)`, start from the cofactor
vector (10). This is an exact decision architecture, not an executed empty
census.

An oriented-matroid or tournament shadow cannot replace the cofactor
sidecar. At rank twelve every positive one-dimensional kernel has the same
all-positive circuit sign. More decisively, THM-2188 has a common
rank-eleven carrier and arbitrarily matching finite residue banks and
pair-overlap matrices, yet opposite zero-/positive-Haar outcomes. The
destroyed coordinate is precisely the Archimedean projective slope retained
by (3), (7), and (10).

An independent hostile audit checked the anchor count, the scope of
THM-763, the 79 order chambers, both inequalities in (8), raw versus
primitive cofactors, weak/strict integer bands, the gcd count (16), the
balanced-sign necessity, and the sparse congruence (22). QED.
