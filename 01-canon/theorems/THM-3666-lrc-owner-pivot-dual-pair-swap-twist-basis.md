---
id: THM-3666
title: "LRC owner-pivot dual pair-swap twist basis"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The dual of every owner-aligned THM-2309 target quotient has an explicit
  sparse basis alpha=e_a-e_ka and beta=e_b-e_kb modulo the scalar gauge.
  Thus THM-2334's abstract target twists act on the present interval product
  by opposite 1/13 shifts of one target blocker and one graft unit, while the
  transported word is exactly neutral.  Combining this with THM-3665 turns
  nonzero target survival into the local pair-swap test
  H(x,y)+H(x-1,y)-2H(x,y-1)!=0.  No such nonvanishing is proved on a
  hypothetical covering row, and no ancestry-digit identification is made.
source: kps-s192 / owner-pivot dualization after THM-3665, 2026-08-21
audit: >
  PASS AFTER ONE SCOPE REPAIR -- agent Schrodinger independently reproduced
  the rank, nullspace, dipole pairings and action, the hostile R=1 control,
  all phase counts, rank 168 versus the rank-156 two-site/dependent controls,
  all 248832 normalized scalar types, all 120 packet choices, and every pinned
  hash.  The audit caught that zeroing one shifted endpoint coefficient is a
  sufficient strategy rather than an equivalent reformulation: the zero-free
  nonconstant profile H(x,y)=2+zeta_13^x is a minimal hostile witness.  Section
  7 now distinguishes the equivalent recurrence/strictness obligations from
  that optional zero-forcing strategy.  THM-2350 is also credited explicitly
  for the earlier dipole normal form and action.
depends_on:
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-3665-lrc-support-minimal-three-twist-target-detector
related:
  - THM-2327-two-colour-marked-unit-c3-triangle
  - THM-2350-owner-pivot-dual-dipole-normal-form
  - THM-3664-lrc-sparse-eight-twist-target-detector-frame
script: 04-computation/lrc_owner_pivot_dual_pair_swap_twists_thm3666.py
output: 05-knowledge/results/lrc_owner_pivot_dual_pair_swap_twists_thm3666.out
script_sha256: 8871c9d8a6eff3d1be33d7944197cffbb8a11c60f6d61132be27a89ce22ff96a
output_sha256: 1015cff6d84f2571cabaaf3e66fad92bee2cba4ba9d616769643dc4809e5c9b9
semantic_sha256: b40649847f46b03b54cc48cb0c6968c861936e4c81441eaaf33adc2f6636ffa0
hash_basis: raw LF bytes
---

# THM-3666 -- target twists are local blocker/unit pair swaps

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
THM-3665 reduced the target-current question to three values of a
function on an abstract plane.  The owner-aligned packet now supplies the
missing concrete coordinates on that plane.

THM-2350 already proved the owner-pivot dipole normal form and its local
action.  The new content here is the adaptation to THM-3665's support-minimal
three-twist recurrence together with the exhaustive 120-choice and scalar-type
ledger; no novelty is claimed for the dipole coordinates themselves.

## 1. Owner-pivot notation

Use THM-2309's labels

```text
U={u_0,...,u_5},
B={j,a,b},                                           (1)
```

where `j` is the selected source blocker, `a,b` are its two target blockers,
and modulo thirteen

```text
w_u!=0 for u in U,
w_j=w_a=w_b=0.                                      (2)
```

Choose an omitted unit `u_0` and distinct graft units

```text
k_a,k_b in U minus {u_0}.                           (3)
```

Let `R_(j,u_0)` be the six owner-aligned rows of THM-2309.  Modulo thirteen,
the ordinary row at pivot label `k` is

```text
r_k=w_k e_(u_0)-w_(u_0)e_k.                         (4)
```

The two graft rows have the additional terms

```text
r_(k_a) <- r_(k_a)-w_(u_0)e_a,
r_(k_b) <- r_(k_b)-w_(u_0)e_b.                     (5)
```

The source row is simply `-w_(u_0)e_j`.  On the pivot columns

```text
P={j} union (U minus {u_0}),                        (6)
```

the matrix is diagonal with nonzero diagonal `-w_(u_0)`.  Hence its row
space `L` has dimension six inside the eight-dimensional scalar kernel
`K=w^perp`.

## 2. The exact dual basis

Define

```text
alpha=e_a-e_(k_a),
beta =e_b-e_(k_b).                                  (7)
```

Equations (4)--(5) give immediately

```text
r.alpha=r.beta=0 for every row r of R_(j,u_0).      (8)
```

The scalar vector `w` is also in `L^perp`.  These three vectors are
independent: `w` has nonzero `u_0` coordinate, while `alpha,beta` do not,
and the latter are separated by their `a,b` coordinates.  Since
`dim L^perp=3`,

```text
L^perp=span(w,alpha,beta).                          (9)
```

THM-2334's target-twist group is

```text
T=L^perp/<w>.                                       (10)
```

Therefore the classes of `alpha,beta` are a basis of `T`.  They are exactly
dual to THM-2309's target quotient basis `e_a,e_b`:

```text
<alpha,e_a>=1, <alpha,e_b>=0,
<beta,e_a>=0,  <beta,e_b>=1.                        (11)
```

This is a lawful relation-quotient duality.  It does not use or construct a
map to the two-current ancestry chart of THM-3657.

## 3. Gauge normalization makes every twist four-sparse

Every class in `T` has the unique representative modulo `<w>`

```text
ell(x,y)=x alpha+y beta,
(x,y) in F13^2.                                     (12)
```

Its only possibly nonzero coordinates are

```text
ell_a=x,       ell_(k_a)=-x,
ell_b=y,       ell_(k_b)=-y.                        (13)
```

All other coordinates, including `j` and `u_0`, vanish.  The phrase
"pair swap" refers to these two opposite coordinate pairs; it is an action
on factor phases, not a permutation of physical speeds.

## 4. Exact action on the THM-2334 current

THM-2334 (31)--(32) realizes a twist `ell` by translating the present
coordinate factors as

```text
I_i(w_i t) -> I_i(w_i t+ell_i/13).                  (14)
```

For (12), only four factors move:

```text
I_a      : +x/13,       I_(k_a): -x/13,
I_b      : +y/13,       I_(k_b): -y/13.             (15)
```

The transported-word factors would shift by

```text
R ell_i/13,
R=13^(lambda_j+1).                                  (16)
```

Every quantity in (16) is an integer.  Hence all word factors are exactly
unchanged on the circle.  This gives a concrete version of THM-2334's
target-neutral/mod-seven-active asymmetry:

```text
target twist = four local present-factor shifts;
transported word = neutral.                         (17)
```

The deep-leg phase remains explicit rather than hidden:

```text
e_13(m ell_c)
=e_13(m[x(1_(c=a)-1_(c=k_a))
       +y(1_(c=b)-1_(c=k_b))]).                     (18)
```

Thus even when it is nontrivial, it depends on at most one of the two pair
coordinates at a time.

## 5. Phase-neutral packet choices

For a fixed source `j`, there are

```text
6*5*4=120                                           (19)
```

choices of `(u_0,k_a,k_b)`.  If the deepest coordinate `c` is a unit, exactly
80 of these choices avoid it in both graft positions, making (18) trivial.
If `c=j`, all 120 choices are phase-neutral because both basis vectors vanish
at the source.  If `c=a` or `c=b`, the corresponding target-coordinate phase
is intrinsic and must be retained.

This freedom is a genuine analytic design parameter: graft units may be
chosen to simplify interval geometry and, except at a target-deep leg, to
remove the extra root-of-unity phase.

## 6. The three-twist test is now a local recurrence

Write

```text
H(x,y)=H(x alpha+y beta).                            (20)
```

THM-3665 turns THM-2334's nonzero-target criterion into

```text
H(x,y)+H(x-1,y)-2H(x,y-1)!=0                       (21)
```

for at least one `(x,y)`.  Each of the three terms in (21) differs from the
others by one of the local pair shifts (15).  This is substantially more
concrete than searching a dense 169-term variance.

If every test in (21) vanished, `H` would satisfy

```text
H(x,y-1)=1/2[H(x,y)+H(x-1,y)]                       (22)
```

on the finite torus.  Equation (22) is harmonicity for an irreducible
directed averaging walk.  The finite maximum-modulus principle forces `H`
to be constant: at a point of maximal modulus, equality in the triangle
inequality propagates the same value along both generating steps.  This is
also an elementary proof of the THM-3665 detector implication in the present
coordinates.

## 7. The sharpened analytic frontier

The target obligation on a hypothetical covering row is equivalent to items
1, 2, and the existential strictness in item 4 below.  Item 3 is instead a
sufficient zero-forcing strategy, not a necessary reformulation: for example,
the nonconstant profile `H(x,y)=2+zeta_13^x` is zero-free.

1. find a packet choice and one `(x,y)` where the pair-swap recurrence (21)
   fails;
2. prove the 169 pair-shifted marked currents cannot be harmonic under (22);
3. force one shifted present rectangle to have a zero endpoint coefficient
   while `H(0,0)` remains the nonzero THM-2327 current;
4. obtain a strict triangle inequality at some maximal-modulus twist.

The most concrete next experiment is to enumerate the 120 packet choices on
the canon-recorded typed rows and record whether one of the three shifted
present products becomes empty, has a unique boundary jump, or acquires a
distinct first nonzero Fourier jet.

No such strictness is proved here.  The all-`91`-unit mask still couples the
septimal action, visible height and terminal phase remain separate sidecars,
and LRC(14) remains open.

## 8. Exact companion

Reproduce with

```bash
python3 -B 04-computation/lrc_owner_pivot_dual_pair_swap_twists_thm3666.py
python3 -B -O 04-computation/lrc_owner_pivot_dual_pair_swap_twists_thm3666.py
```

Both streams match the stored transcript.  The assertion-free companion
source-pins THM-2309, THM-2334, and THM-3665; exhausts all `12^5=248832`
normalized unit scalar types for a fixed label chart; checks all 120 packet
choices, all 169 pair-swap twists, word neutrality through eight depth
exponents, and the exact phase-avoidance counts.  **QED.**
