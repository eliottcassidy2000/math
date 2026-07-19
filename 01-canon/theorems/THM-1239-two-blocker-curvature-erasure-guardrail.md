---
id: THM-1239
title: ONE RESONANT BLOCKER CAN ERASE THE ENTIRE SEVEN-EDGE CURVATURE PATH — exact toothpick self-similarity and the failure of arbitrary-gap active-edge survival
status: PROVED (one-blocker family for every m>=42; two-blocker family for every m>=8; uniform non-BAD deletion quartet; exact global lonely reroute; dependency-free referee; sorry-free Lean arithmetic core)
source: codex-2026-07-19-S78 continuation with active-blocker-scan agent
depends_on: [THM-1156, THM-1219]
related: [THM-1218, THM-1237, THM-1238, THM-1240, HYP-7870]
script: 04-computation/lrc14_two_blocker_curvature_erasure_guardrail_thm1239.py
output: 05-knowledge/results/lrc14_two_blocker_curvature_erasure_guardrail_thm1239.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCCurvatureErasureGuardrail.lean
script_sha256: 5b40596ddb8443d542d5b4de23821e099451d447058f55aa9e41175a267db03d
output_sha256: 043632dbea33bd970b67eb4551c9acc4a56707fe534f7411a7a2654fe8ec5a3e
formalization_sha256: bc37023953410855e2218657c16602322aa8bc0415975c23b54f0a7943282275
---

# THM-1239 — one blocker can erase the curvature path

## 1. The curvature cracks

Let `m>=1`, put

```text
a=7m+1,
P={a,a+1,...,a+6},
G=G_m(a)=[(14m+1)/(14a),(14m+13)/(14a)].              (1)
```

THM-1219 shows that the six teeth of `a+1,...,a+6` leave seven positive
components in `G`:

```text
H_L=[G_L,L_6],
H_r=[R_(r+1),L_r],                  r=5,...,1,
H_R=[R_1,G_R].                                         (2)
```

These are exactly `G intersect Safe_P`.

## 2. A two-blocker family from m=8

Set

```text
c=a+7,                 d=7a+4.
```

For every `m>=8`, the `c`-tooth at address `m+1` strictly contains `H_L`, the
`c`-tooth at `m+2` strictly contains `H_R`, and the `d`-tooth at

```text
n_r=7m+7-r                                               (3)
```

strictly contains `H_r` for `r=1,...,5`.  The exact cleared boundary margins
for `c` are

```text
(14m-5,1),                (66,14m-75),                  (4)
```

and the two margins for `d` on `H_r` are

```text
14m+14r^2-69r-29,
(11-2r)(7r-4).                                           (5)
```

Across `r=1,...,5`, these are

```text
14m-(84,111,110,81,24),
(27,70,85,72,31).                                       (6)
```

They are all positive from `m=8`.  Thus the two blocker incidence words on
the ordered cracks are

```text
c: 1000001,                 d: 0111110.                (7)
```

The smallest row is `a=57,c=64,d=403`.

## 3. The one-blocker toothpick ladder

The construction is much more rigid and self-similar than (7) suggests.  Put

```text
u=14a(t-1/7).                                           (8)
```

Then `G` becomes exactly `[-1,11]`, while the `(a+r)` tooth becomes

```text
( a(11-2r)/(a+r), a(13-2r)/(a+r) ).                    (9)
```

As `a` tends to infinity, the six teeth in (9) tile

```text
(-1,1),(1,3),(3,5),(5,7),(7,9),(9,11),                 (10)
```

and the seven curvature cracks collapse onto the odd joints

```text
J={-1,1,3,5,7,9,11}.                                   (11)
```

Their exact `u`-intervals are

```text
H_L=[-1,-a/(a+6)],
H_r=[q_r a/(a+r+1),q_r a/(a+r)], q_r=11-2r,
H_R=[11a/(a+1),11].                                    (12)
```

Now take the single blocker

```text
z=14a.                                                  (13)
```

Equation (8) gives the exact affine law

```text
zt=2a+u.                                                (14)
```

Thus the `z`-danger comb in `u`-space is the radius-`1/14` needle lattice
around every integer.  Each crack lies strictly inside the needle at its odd
joint once

```text
6/(a+6)<1/14,
q_r(r+1)/(a+r+1)<1/14,             r=1,...,5,
11/(a+1)<1/14.                                         (15)
```

The worst row is `r=2`:

```text
21/(a+3)<1/14
iff a>291
iff m>=42.                                              (16)
```

Therefore, for every `m>=42`, one speed `z=14a` strictly covers all seven
components of `G intersect Safe_P`.  Its tooth addresses form the exact
odd ladder

```text
14m+1,14m+3,...,14m+13.                                (17)
```

This is literal toothpick self-similarity: the cracks shrink like `a^-2`,
but multiplication by the resonant slope `14a` reproduces them at seven
integer needles.

## 4. A deletion quartet does not rescue local coherence

Add the quartet

```text
Q={1,2,3,4}.                                            (18)
```

For every `t in G`, translating its phases by `t` gives

```text
{0,t,2t,3t}.
```

The last cyclic gap is `1-3t`, and uniformly

```text
1-3t-2/7 >=(28m-29)/(14a)>0.                           (19)
```

Hence the entire selected gap is outside the quartet's THM-1203 BAD event,
yet the blocker in (13) consumes every point of `Safe_P intersect G`.
Therefore any arbitrary-gap six-killer coherence lemma using only “outside
one BAD quartet” is false.

The construction is local, not an LRC counterexample.  At the first
one-blocker row, the complete thirteen-speed packet

```text
V={1,2,3,4,295,296,297,298,299,300,301,302,4130}       (20)
```

has the exact lonely phase

```text
t=44/199.                                               (21)
```

Its least-residue numerators, in the order (20), are

```text
44,88,67,23,45,89,66,22,22,66,89,45,33,               (22)
```

so the minimum is `22/199>1/14`.  The proof must move to another gap/address
cell; it cannot demand survival in an arbitrary selected gap.

## 5. The address-word invariant

For a closed crack `H=[L,R]`, a speed `z` strictly contains it in one danger
tooth exactly when the address window

```text
(zR-1/14,zL+1/14)                                      (23)
```

contains an integer.  Its length is

```text
1/7-z|H|<1,                                            (24)
```

so the containing tooth address is unique.  If one blocker contains two
cracks with midpoint separation `Delta u` and address jump `h`, then in the
scaled coordinates

```text
|(z/(14a)) Delta u-h|<1/7.                              (25)
```

Since the limiting crack separations are even integers, a high-degree
blocker is forced near a rational slope channel

```text
z/a approximately 7h/k.                                (26)
```

The exact affine needle is

```text
zt=z/7+(z/(14a))u.                                     (27)
```

At a limiting odd joint `q`, the hit test is therefore

```text
dist((z mod 7)/7+(z/a)q/14,Z)<=1/14.                   (28)
```

The three observed channels are `(residue,slope)=(1,1)`, hitting the two
boundary joints; `(4,7)`, hitting the five internal joints; and `(0,14)`,
hitting all seven.  A future closing theorem must classify these finite
rational-slope/address words and dispatch their resonant rays by a different
gap, beat, or divisor argument.

## 6. Tournament and assumption audit

The runner-order tournament is transitive.  The chronological address order
is also transitive.  Neither sees that one blocker realizes a full `K_(1,7)`
incidence star.  The chi-seven/Fano seam obstruction does not apply: these are
strict third-tooth containments, not exact two-owner abutments.

We challenged runners, gaps, endpoints, active pair edges, crack components,
blocker labels, tooth addresses, rational slopes, and proof obligations as
vertices.  The faithful current vertices are the seven crack/address
obligations, with blocker address word and global gap cell as sidecars.  This
carrier preserves local cover truth and exposes the resonant self-similarity;
it destroys the alternate-gap witness (21), which must remain a global
transport coordinate.

## 7. Verification and scope

The exact referee checks `9,993` two-blocker rows (`m=8,...,10000`), `9,959`
one-blocker rows (`m=42,...,10000`), `3,000` exact scaled-tooth identities,
the sharp failure at `m=41`, the non-BAD margin, and all residues in (22).
Normal and optimized outputs are byte-identical.

The Lean module kernel-checks (8)--(9), every inequality in (15)--(16), the
two-blocker cleared margin ledger, (19), the explicit witness residues and
margin, and the uniqueness scale (24).  Identification of the paper crack
intervals with the safe components remains explicit; there are no proof
placeholders.

Frozen hashes are

```text
source         5b40596ddb8443d542d5b4de23821e099451d447058f55aa9e41175a267db03d
output         043632dbea33bd970b67eb4551c9acc4a56707fe534f7411a7a2654fe8ec5a3e
formalization  bc37023953410855e2218657c16602322aa8bc0415975c23b54f0a7943282275
```

THM-1239 refutes local active-edge survival, blocker-count, and bare
deletion-incidence closure routes.  It does not refute the positioned transfer
of THM-1237, the carrier-spoke cycle of THM-1240, or LRC(14).
