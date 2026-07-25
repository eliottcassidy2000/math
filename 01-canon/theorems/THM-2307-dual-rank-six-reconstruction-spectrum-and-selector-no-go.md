---
id: THM-2307
title: "Dual rank-six reconstruction spectrum and bounded-selector no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Let a
  nine-coordinate word over F_13 have exactly three zero blocker
  coordinates and six nonzero guard/unit coordinates, and let R be a
  six-dimensional relation subspace orthogonal to it. Either R contains
  a literal blocker coordinate vector e_b, or the dual code R-perp has
  at least 77 projective all-unit words, equivalently 924 nonzero
  weight-nine words. The bound is sharp. A bank of selectors whose
  restricted dehomogenized product is nonzero and has total degree at
  most five still leaves at least 12 all-unit impostors. At rank seven
  the corresponding sharp floors are 7 and 2; rank eight is the first
  linear reconstruction point. THM-2299's exact height-13 anchored
  rank-six phase-cancellation row lies in the literal-blocker branch.
  Thus neither branch alone supplies the missing terminal phase, no
  scalar profile is excluded, and LRC(14) remains open.
source: codex-2026-07-25-dual-rank-six-reconstruction-spectrum
depends_on:
  - THM-2299-rooted-current-service-energy-and-base-phase-no-go
related:
  - THM-2069-k-deletion-code-cogirth-crt-wheel
  - THM-2082-rank-one-code-wheel-blindness-and-translated-prime-grid
  - THM-2280-centered-polynomial-grid-avoidance-and-bounded-generic-keller-fibre
  - THM-2285-centered-grid-footprint-and-generic-keller-lines
  - THM-2287-anchored-scalar-rank-six-plucker-flag-and-finite-label-atlas
  - THM-2301-essential-affine-arrangement-and-visible-rank-six-address-bank
  - THM-2303-terminal-component-phase-current-and-defect-rank
script: 04-computation/lrc14_dual_rank_six_reconstruction_spectrum_thm2307.py
output: 05-knowledge/results/lrc14_dual_rank_six_reconstruction_spectrum_thm2307.out
script_sha256: 29af7bec855a51a6ff798ff9a217d267d587443e932124a626b4ef2b39709af7
output_sha256: 3b65e8ce9fb3101a3741559a70bebc10cb3333d87567265ad113abb15fe8ff46
hash_basis: working-tree bytes (LF)
---

# THM-2307 -- the dual rank-six reconstruction spectrum

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2301 and THM-2287 construct six-dimensional mod-thirteen relation
packets. Their coefficient spaces contain many relation addresses. The
orthogonal question is different:

```text
which scalar residue words satisfy every relation in the packet?
```

The physical scalar word has six unit coordinates and three blocker zeros.
At rank six, the relation packet either sees one blocker as a literal
coordinate residue or admits a sharp bank of `77` all-unit scalar
impostors. Even five degrees of additional algebraic selection cannot
remove them all. This is a dual-code obstruction, not another primal
Kakeya count.

## 1. Relation and scalar codes

Work over

```text
k=F_13.
```

Partition the nine labelled coordinates as

```text
U={six guard/unit labels},
B={three blocker labels}.                           (1)
```

Let

```text
wbar in k^9
```

satisfy

```text
wbar_u!=0 for u in U,
wbar_b=0  for b in B.                               (2)
```

Let

```text
R subset wbar^perp,       dim_k R=r,
D=R^perp.                                           (3)
```

Thus `wbar in D` and

```text
dim_k D=9-r.                                        (4)
```

Call a blocker `b` **literally dark** if

```text
e_b in R.                                           (5)
```

By nondegeneracy of the coordinate pairing, (5) is equivalent to

```text
x_b=0 for every x in D.                             (6)
```

This literal coordinate condition is essential. Because the whole
three-dimensional blocker coordinate space lies in `wbar^perp`, a
higher-support blocker relation may occur for dimensional reasons. Such a
relation does not make any one blocker coordinate vanish on `D`.

No guard/unit coordinate can be literally dark: if `e_u in R`, then

```text
0=e_u.wbar=wbar_u,
```

contrary to (2).

An **all-unit dual word** is a vector

```text
x in D,             x_i!=0 for all i.               (7)
```

Its projective direction is unchanged by the twelve nonzero scalar
multiples.

## 2. The sharp rank-six alternative

Assume

```text
r=6,               dim D=3.                         (8)
```

Then exactly one of the following alternatives holds:

1. `e_b in R` for some blocker `b`;
2. `P(D)` contains at least `77` all-unit projective directions.

Equivalently, in the second branch the weight enumerator of the dual code
has

```text
A_9(D)>=12*77=924.                                  (9)
```

### Proof

Assume the first alternative fails. In the projective plane

```text
P(D)=P^2(k),
```

let

```text
W=[wbar].
```

For each blocker `b`, the equation

```text
x_b=0
```

cuts a genuine projective line through `W`: it contains `W` by (2), and it
is not the whole plane by the failure of (5). Let `t`, with

```text
1<=t<=3,
```

be the number of distinct blocker lines.

There are

```text
13+1=14                                             (10)
```

projective lines through `W`. Remove the `t` blocker lines. Every remaining
line `L` contains exactly thirteen points other than `W`. None of those
points lies on a blocker line, since two distinct projective lines through
`W` meet only at `W`.

For every guard/unit label `u`, the coordinate-zero line

```text
x_u=0
```

does not contain `W`, because `wbar_u!=0`. It therefore meets `L` in exactly
one point. The six guard/unit lines remove at most six of the thirteen
points on `L`. At least seven all-unit points remain on every surviving
pencil line. Hence

```text
# all-unit projective directions
 >=(14-t)(13-6)
 >=11*7
 =77.                                                (11)
```

Coincident coordinate lines only improve the bound. This proves the
alternative and (9). QED.

The same count can be written in the affine chart of any unit coordinate.
If the `t` blocker lines are distinct, their union has `1+12t` points.
Each of the other five unit-coordinate lines is parallel to at most one
blocker line and adds at most `14-t` new points. Thus the bad set has size
at most

```text
1+12t+5(14-t)=71+7t<=92,                            (12)
```

again leaving `169-92=77`.

## 3. Sharpness

Take abstract dual coordinates

```text
(x,y,z) in k^3
```

and use the six guard/unit coordinate forms

```text
z, x-z, x-2z, x-3z, x-4z, x-5z                    (13)
```

together with the three blocker forms

```text
x, y, x-y.                                          (14)
```

The coordinate map from `k^3` to `k^9` is injective. Put

```text
wbar=(x,y,z)=(0,0,1).
```

All six forms in (13) are nonzero on `wbar`, and all three forms in (14)
vanish there. None of the blocker forms vanishes identically, so there is
no literal dark blocker.

Every all-unit direction has `z!=0`, hence a unique representative

```text
(x,y,1).
```

Its coordinates are all nonzero exactly when

```text
x notin {0,1,2,3,4,5},
y notin {0,x}.                                      (15)
```

There are therefore exactly

```text
7*11=77                                             (16)
```

all-unit directions. The bad affine union has exactly `92` points. Thus
(11) is sharp.

If `D` is the image of this coordinate map and `R=D^perp`, then

```text
dim R=6,
R subset wbar^perp,
```

so this is a literal sharp code model, not only a line-arrangement
picture.

## 4. Bounded algebraic selectors

Fix any guard/unit anchor `a` and normalize

```text
A_a={x in D:x_a=1}.                                 (17)
```

At rank six, `A_a` is an affine plane over `k`, and every all-unit
projective direction has exactly one representative in it.

Let

```text
F_1,...,F_s in k[X_0,...,X_8]
```

be a finite selector bank. Restrict and dehomogenize the selectors on
`A_a`, and assume that their product

```text
f=product_j F_j|_(A_a)                              (18)
```

is a **nonzero polynomial** of total degree

```text
d=sum_j deg(F_j)<=5.                                (19)
```

Schwartz--Zippel on `k^2` gives

```text
#{x in A_a:f(x)=0}<=13d.                            (20)
```

Combining (20) with (11), at least

```text
77-13d                                              (21)
```

all-unit impostors survive every selector in the bank. For `d=0,...,5`,
the floors are

```text
77,64,51,38,25,12.                                  (22)
```

The nonzero-restriction hypothesis in (18) is load-bearing. An ambient
nonzero linear form coming from any row of `R` vanishes identically on
`D`; calling it a selector would make the statement false and add no
information. Because `d<13`, a nonzero restricted polynomial in (18)
cannot represent the zero function merely through the finite-field
relations.

## 5. Ranks seven and eight

Suppose first that

```text
r=7,               dim D=2,                         (23)
```

and no blocker is literally dark. In the projective line `P(D)`, every
blocker coordinate form has a one-dimensional kernel containing `wbar`.
Hence all three blocker zero sets are the same projective point `W`.
Each of the six guard/unit forms removes at most one other point. Since

```text
#P^1(F_13)=14,
```

at least

```text
14-1-6=7                                            (24)
```

all-unit directions remain.

This is sharp. On dual coordinates `(x,z)`, use the six unit forms

```text
z,x-z,x-2z,x-3z,x-4z,x-5z
```

and take all three blocker forms equal to `x`, with physical point
`(0,1)`. The all-unit directions are exactly

```text
(x,1),       x notin {0,1,2,3,4,5},
```

so there are seven.

On the normalized affine line, a nonzero restricted selector product of
degree `d<=5` has at most `d` roots. Therefore at least

```text
7-d>=2                                              (25)
```

all-unit directions survive.

Finally, if

```text
r=8,
```

then `R` and `wbar^perp` are both eight-dimensional and one contains the
other. Consequently

```text
R=wbar^perp,
D=span(wbar).                                       (26)
```

Rank eight is therefore the first point at which the linear relation code
reconstructs the physical scalar word up to scale. In particular, all
three blocker coordinate vectors lie in `R`.

The exact reconstruction spectrum is

```text
rank 6: literal dark blocker, or at least 77 all-unit impostors;
rank 7: literal dark blocker, or at least 7 all-unit impostors;
rank 8: unique physical projective word.             (27)
```

## 6. The THM-2299 dark-branch obstruction

The literal-dark alternative is not a phase certificate. THM-2299's exact
local profile-`(1,3,5)` witness has scalar row

```text
(H,q_1,...,q_5,c_1,c_2,c_3)
=(1,4,2,3,6,10,13,2197,742586).                    (28)
```

Its five tautological unit relations

```text
e_(q_i)-q_i e_H,            1<=i<=5,
```

together with

```text
p=13e_(q_1)-4e_(c_1)                               (29)
```

form a rank-six packet of height at most `13`. On the pivot columns

```text
(q_1,q_2,q_3,q_4,q_5,c_1)
```

the minor is `-4`, a thirteen-unit. But modulo thirteen,

```text
p=-4e_(c_1),                                       (30)
```

so this packet lies in the literal-dark-blocker branch.

Nevertheless, THM-2299 proves on the same local labelled carrier that

```text
F_hat(4)=E_hat(52)=(W_4)_hat(0)=0                  (31)
```

by antipodal cancellation between two terminal components. The owner,
prescribed clock, target blocker, one-sheet rooted word, and anchored
rank-six packet all survive.

This is a local-interface obstruction, not a global scalar cover or an
LRC counterexample. It proves exactly that the first branch of (27) does
not recover terminal base phase. In the other branch, the `77`-word bank
proves that the mod-thirteen linear code does not even recover the physical
blocker-zero pattern. THM-2303's terminal-component current and relative
`U(1)` phase transport identify one faithful missing sidecar.

## 7. Connection and loss ledger

```text
source:
  the six-row relation packets of THM-2301/THM-2287, the physical
  six-unit/three-zero scalar residue word, and THM-2299's phase witness;

map:
  replace the primal relation-address space by its orthogonal scalar code
  D=R^perp, project from the physical word, and count the surviving
  weight-nine directions;

preserved:
  every mod-thirteen relation in R, all nine coordinate labels, the
  guard/unit versus blocker type of the physical word, and polynomial
  selectors that remain nonzero after restriction to D;

new output:
  the sharp 77/7/1 reconstruction spectrum, 924 nonzero rank-six
  weight-nine words, and the degree-five survivor floors 12 and 2;

destroyed:
  integer height beyond the chosen packet, higher thirteen-adic digits,
  valuation depth, owner ancestry, endpoint address, Fourier amplitude,
  and terminal-component phase;

sharp hostile controls:
  the affine model (13)--(16), ambient selectors vanishing identically on
  D, generic blocker-supported relations mistaken for literal e_b, and
  THM-2299's height-13 dark-blocker phase cancellation;

needed sidecar:
  higher digits or owner/endpoint data to identify the physical dual word,
  together with a component-phase transport when a blocker is dark.     (32)
```

No scalar profile is excluded. LRC(14) remains open.

## 8. Exact reproduction

Run

```bash
python3 04-computation/lrc14_dual_rank_six_reconstruction_spectrum_thm2307.py
python3 -O 04-computation/lrc14_dual_rank_six_reconstruction_spectrum_thm2307.py
```

Both executions must match

```text
05-knowledge/results/lrc14_dual_rank_six_reconstruction_spectrum_thm2307.out
```

byte-for-byte after LF normalization. The companion enumerates the sharp
rank-six and rank-seven models, checks the complete count and selector
invoices, and independently reconstructs THM-2299's height-13 pure-blocker
packet. QED.
