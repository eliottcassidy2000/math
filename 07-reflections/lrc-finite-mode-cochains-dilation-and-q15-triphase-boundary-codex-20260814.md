# LRC sheet covers: finite modes, dilation, and the first effective triphase

**Date:** 2026-08-14  
**Status:** FINITE-EXACT literal `q=8,...,14` rank-at-most-five reconstruction,
exact `q=15` rank-at-most-six boundary model, and independent full physical
clutter audit through `q=15`; PROVISIONAL ANALYTIC finite-mode theorem
candidate.  No theorem is promoted here.  The `q=15` model is outside the
literal THM-3387 body/divisor atlas, gives no refined decrement, and does not
prove LRC(14).

## 1. Inheritance and the corrected question

THM-3387 identifies the obstruction as a cover of the cyclic sheet fibre.
THM-3395 proves that through seven sheets every firing owner contributes one
kernel coset and that a typed cover is physical exactly when its affine gap
fibres contain a closed complete cochain.  The canonical hostile is the
`q=6` `4+2` partition whose pair fibres are all nonempty but whose star fibre
is empty.  The corrected near miss is therefore “pairwise overlap implies a
common phase.”  The least-used sidecar is the owner mode: its sheet block,
centre lattice, and interval radius.

The live board is:

1. cyclic sheet blocks;
2. affine integral `1`-cochains;
3. quotient/dilation maps between sheet degrees;
4. Boolean covers versus forced tournaments;
5. multiplicative ancestry and harmonic realizations;
6. target-category effectivity, with THM-3383 as the JC hostile.

The new question is not whether the cochain dies at `q=8`.  It does not.  The
question is how to type one owner after a danger arc can meet several phase
classes.

## 2. The finite witness-mode theorem candidate

Fix `q>=2`, a transverse speed `u`, and write

```text
g=gcd(u,q),                         m=q/g.              (1)
```

The map from sheets to distinct phase classes is

```text
phi_u(ell)=(u/g)ell mod m.                              (2)
```

Every nonempty set of phase classes simultaneously lying in the open danger
arc is a proper consecutive cyclic block.  A useful witness bank may contain
every subblock

```text
B(r,s)={r,r+1,...,r+s-1} mod m,
1<=s<=ceil(m/7).                                        (3)
```

This bank is intentionally redundant.  The *actual* dangerous block has
`ceil(m/7)-1` or `ceil(m/7)` classes (discarding the empty case), whereas a
smaller block in `(3)` merely certifies a subset of the sheets that are
dangerous.  Confusing “actual block” with “witness subblock” is the first
failure mode beyond `q=14`.

The sheet block, cleared centre, and cleared width of `(3)` are

```text
K_u(r,s)=phi_u^(-1)(B(r,s)),
h=-g(2r+s-1) mod 2q,
w=g(m-7(s-1)).                                          (4)
```

All sheets in `K_u(r,s)` are dangerous precisely on translates of the open
source interval with

```text
centre lattice:  (1/u)Z+h/(2qu),
radius:          w/(14qu).                              (5)
```

The width is positive exactly in the range `(3)`.  Because the danger arc is
shorter than a half-circle and `u` is transverse, the cyclic start `r` gives
a canonical wraparound convention; changing the integer lift changes only
the `(1/u)Z` translate.

For selected modes `(K_i,h_i,w_i)` and lifted centres `x_i`, put

```text
p_ij=2q u_i u_j(x_i-x_j).                               (6)
```

The exact pair fibre is

```text
p_ij == h_i u_j-h_j u_i  (mod 2q gcd(u_i,u_j)),
7|p_ij| < w_i u_j+w_j u_i.                              (7)
```

It becomes one physical source phase only when the complete antisymmetric
array closes:

```text
u_k p_ij+u_i p_jk+u_j p_ki=0                            (8)
```

for every triple.  Equivalently one compatible star determines every pair.

The proposed exact theorem is

```text
B_q(U) is nonempty
iff selected owner modes cover Z/qZ and admit a cochain (7)-(8). (9)
```

Necessity retains the actual consecutive dangerous block at a covering
phase.  For sufficiency, `(8)` integrates the normalized gaps to rational
potentials; the congruences in `(7)` make their centre-lattice cosets
pairwise compatible, so generalized CRT supplies simultaneous legal centre
lifts.  The strict inequalities in `(7)` make their real intervals pairwise
intersect, and one-dimensional Helly supplies a common source time.  At that
time every selected sheet block is dangerous and their union covers the
fibre.

This proof is the same closure/effectivity mechanism as THM-3395, with the
single-coset type replaced by a finite mode bank.  It is still awaiting an
independent full analytic audit before canon promotion.

## 3. Exact `q=8,...,14` reconstruction

The exact event route and the mode-cochain route agree on every minimal
obstruction through rank five, the complete rank range needed by the literal
body atlas.  The frozen low-rank census is:

| `q` | minimal edges through rank 5 | ranks through 5 | independent five-sets `I_5` |
|---:|---:|---:|---:|
| 8 | 32 | `15` of rank 4, `17` of rank 5 | 1,152 |
| 9 | 22 | `9` of rank 4, `13` of rank 5 | 1,205 |
| 10 | 18 | rank 5 | 1,269 |
| 11 | 0 | none | 1,287 |
| 12 | 8 | `1` of rank 4, `7` of rank 5 | 1,271 |
| 13 | 0 | none | 1,287 |
| 14 | 0 | none | 1,287 |

Every unsafe five-set at these degrees fires outside core clock `1`.  Hence
there are no new core rescues: the exact row vector is

```text
(1152,1205,1269,1287,1271,1287,1287).                  (10)
```

This reconstructs the corresponding THM-3387 slices structurally but does
not add or subtract any ledger row.

It does **not** classify the full physical clutters.  An independent route
using every rational boundary and intervening open-cell midpoint, with no
mode/cochain import, gives the true full rank profiles

```text
q=8:  (4:15,5:17,6:6)       q=12: (4:1,5:7,6:17,7:14)
q=9:  (4:9,5:13,6:54,7:2)   q=13: (7:22,8:29,9:2)
q=10: (5:18,6:70,7:4)       q=14: (7:3,8:60,9:9)
q=11: (6:23,7:91,8:9)       q=15: (6:157,7:16,8:6).   (10a)
```

The first omitted q=8 edge is `(1,3,5,11,13,14)`.  In particular, the zeros
at q=11,13,14 in the first table mean “no edge through rank five,” not “empty
physical clutter.”  The full audit leaves every `I_5` and no-rescue statement
in `(10)` unchanged.

The sharp local boundary is now exact:

```text
q<=7:   one phase class per witness owner;
8<=q<=14: singleton or adjacent-domino witness modes suffice;
q=15:   a three-phase witness mode first becomes available. (11)
```

## 4. Pure dilation is a functor, not three coincidences

Let `d>=1`.  Simultaneously replace

```text
(q,u,t) by (dq,du,t/d).                                 (12)
```

Then

```text
(du)(t/d+k/(dq))=u(t+k/q).                              (13)
```

Thus a sheet block at degree `dq` is exactly the inverse image of the old
block under `Z/dqZ -> Z/qZ`.  The phase count `m` is unchanged, while

```text
g,h,w       -> d(g,h,w),
centre/radius -> (centre/radius)/d,
p_ij        -> d^2 p_ij.                                (14)
```

Both the congruence modulus and the overlap bound in `(7)` also multiply by
`d^2`, so cover, closure, and physical realization correspond in both
directions on the pure `d`-multiple stratum.  This one identity explains the
three frozen lifts

```text
q4 -> q8:   (1,3,5,7)       -> (2,6,10,14),
q5 -> q10:  (1,2,3,4,7)     -> (2,4,6,8,14),
q6 -> q12:  (2,3,5,7)       -> (4,6,10,14).             (15)
```

Mixed edges need not descend because different owners may live in different
gcd strata.  The mode is the sidecar destroyed by forgetting that typing.

## 5. `q=15`: local availability versus global effectivity

The literal speed pool `{1,...,14}` has no cover of fifteen sheets with at
most five owners.  At rank six there are exactly `157` minimal edges.  Of
these, `155` admit no cover using only singleton/domino modes; only

```text
(1,6,7,10,11,13),       (1,7,10,11,12,13)              (16)
```

remain domino-sufficient.  Thus the first local trimode boundary is also a
genuine rank-six effectivity boundary.

The first edge is especially transparent:

```text
E={1,2,3,4,5,7}.                                        (17)
```

At source centre zero choose

```text
speed 5: {0,3,6,9,12},
speed 3: {0,5,10},
speed 1: {0,1,14},
speed 2: {0,7,8},
speed 4: {0,4,11},
speed 7: {0,2,13}.                                      (18)
```

Every cleared centre is `h=0` and all fifteen pair gaps are `p_ij=0`.
The first two blocks cover exactly the zero divisors of `Z/15Z`.  Removing
zero from the last four gives the unit sign-pairs

```text
{+-1}, {+-7}, {+-4}, {+-2}.                             (19)
```

Each pair consists of private sheets, so replacing its displayed three-phase
block by either adjacent domino loses one private sheet.  An independent
minimization over **all** legal mode witnesses for the fixed edge `(17)`
confirms that its minimum triphase count is four.  This is edge-specific: the
two rank-six edges in `(16)` need zero trimodes, while the other 155 need at
least one.  It is not a statement that every q=15 cover needs four.

The four unit modes are indexed by

```text
(Z/15Z)^x/{+-1}={1,2,4,7},
[x] -> [2x]: 1 -> 2 -> 4 -> 7 -> 1.                    (20)
```

This is a natural four-state directed cycle with missing diagonals.  It is
not a tournament certificate, and the physical cochain is in fact a
complete tie.  The Boolean realization is the load-bearing object:

```text
two kernel blocks cover the zero divisors;
four trimodes cover the four unit sign-classes.          (21)
```

## 6. Trees, recurrences, and harmonic subsets

Equation `(20)` is a multiplicative recurrence.  It should not be renamed a
Fibonacci recurrence: the Fibonacci/Stern-Brocot carriers retain ordered
`L/R` ancestry, whereas multiplication in the sheet quotient has already
identified commuting words.  A lawful bridge must therefore retain the word
or exponent vector as a sidecar.

The dilation identity does give a genuine degree-graded ancestry.  Starting
from any finite edge and a free multiplicative monoid generated by selected
primes, exponent vectors index descendants, fibre degrees multiply, and the
corresponding integers form a subset of the natural numbers.  Its harmonic
weight factors as an Euler product whenever root collisions have been
excluded:

```text
sum_(a in monoid) sum_(u in root) 1/(au)
= (sum_(u in root)1/u) product_(p)(1-1/p)^(-1).          (22)
```

This is an exact realization of such an ancestry as a subset of the harmonic
series.  It records mass, not order.  A Fibonacci or ternary branch tree can
map into it only after specifying which ordered branches collide under
commutative multiplication.

## 7. Cross-frontier guardrail

The finite-mode cochain fits the repo's compatibility/closure/effectivity
grammar:

```text
pair gap fibres -> zero-circulation star -> one physical source phase. (23)
```

THM-3383 shows why this is only a grammar shared with the Jacobian lane.
There a rational torsor may close while a boundary principal part remains
nonpolynomial.  No map from the LRC cochain to a JC `H^1` class is supplied
here, and no JC consequence follows.  Likewise the unit four-cycle in `(20)`
does not license a tournament transfer without an intrinsic oriented pair
observable.

## 8. Reproduction and scope

Reproduce the literal reconstruction with

```text
python 04-computation/lrc14_q8_q14_finite_mode_clutter_probe_20260814.py
python -O 04-computation/lrc14_q8_q14_finite_mode_clutter_probe_20260814.py
```

and the boundary model with

```text
python 04-computation/lrc15_first_effective_triphase_mode_probe_20260814.py
python -O 04-computation/lrc15_first_effective_triphase_mode_probe_20260814.py
```

Reproduce the independent full physical clutter audit with

```text
python 04-computation/lrc_q8_q15_full_physical_clutter_audit_20260815.py
python -O 04-computation/lrc_q8_q15_full_physical_clutter_audit_20260815.py
```

Their semantic digests are respectively
`66a69a30c49b72ff8ecbf7de94f495025518e04b73969f2d970debeb6f113023`
and
`bae0a83a2ca7f19906309b1bab951bd9e4dac241e3a3103726da3bb8099e1eb8`.
The full-clutter audit semantic digest is
`2baecf15b57ac7adcfaf41342d74cf12eba363a6d0726566e67eabfcf78b0d22`.

The `q=15` census is a sharp hostile/positive boundary for the proposed
general theorem.  Because no speed in `{1,...,14}` is divisible by `15`, it
is not a THM-3387 literal body row and cannot change the LRC(14) ledger.
