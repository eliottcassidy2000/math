# Every two-ray LRC14 carrier configuration has three strict network certificates

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The `nc2_seed` referee checked the determinant and reciprocal arguments,
complete direction counts, strict cutoffs, and the full finite-head replay
using congruence rows, independent relation boxes, and physical sheet graphs.
No scarce theorem namespace or shared navigation is changed by this note.

Let `w=(a,b,c)` be primitive, distinct, positive, odd, ternary-unit, and sorted.
Use the complete raw carrier support and exact three network projections
`E_i` of
[THM-4414 — six-separated contact capacity collapse](../../01-canon/theorems/THM-4414-lrc14-six-separated-contact-capacity-collapse.md).
Suppose the support consists of exactly **two primitive directions**, where
a direction identifies `v` and `-v` and retains all positive multiples as
separate raw carriers. Then

```text
E_i<6/77 for all three i, at every height.              (1)
```

Together with the complete one-ray theorem in
[the companion note](lrc14_one_ray_overnight_hexagon_sep05.md), this leaves
at least three primitive live directions in every hypothetical failure of
the universal network target. This does not imply chart entry,
synchronization, or LRC(14).

## Inheritance and the repaired local-witness operation

The one-ray note bounds all multipliers on a direction using its largest
primitive coefficient, but explicitly refutes discarding a second direction:
the first hostile is `(17,23,25)`. Here both directions are retained. Their
interaction is an integer cross-product invariant. The least-used sidecar is
the factor of three forced into that cross product by the owner residues.

The closest proved mechanism is
[THM-4386 — canonical component relation and zero-defect incidence,
Lemma 2.1](../../01-canon/theorems/THM-4386-lrc14-canonical-component-relation-and-zero-defect-incidence.md).
It prevents any short owner-live relation from coexisting with an independent
live direction. THM-4422's first dense multi-ray `(19,23,29)` has **three**
directions and remains outside this theorem.

The concept board is: primitive ray; raw multiplicity; mod-three owner coset;
integer determinant; reciprocal coefficient bound; exact network consumer.
The map keeps both complete ordinal multiplier lists and their determinant.
It preserves their sum of network loads. Unlike a geometric extremal
replacement, it discards no live point; the local witness only supplies a
constraint on the other direction. No point-set realizability or SAT
assumption is imported from the empty-hexagon paper.

## 1. Every multi-ray direction has a coefficient of magnitude at least seven

Write the two primitive directions as `u,v`. They are independent and each
has all three coordinates nonzero modulo three. If either had `l1` norm at
most fourteen, THM-4386 Lemma 2.1 would make its relation defect zero on
every live component. Every carrier would then be parallel to that relation,
contradicting the second direction. Thus

```text
||u||_1>=16, ||v||_1>=16.                               (2)
```

The norms are even because all speeds are odd. A ternary-unit coordinate
cannot have magnitude six. If every magnitude were at most five, the norm
would be at most fifteen, hence at most fourteen by parity. Therefore

```text
M_u=max_i|u_i|>=7, M_v=max_i|v_i|>=7.                   (3)
```

This argument applies to every direction whenever the support is multi-ray;
it is not special to having exactly two.

## 2. The owner gate forces a determinant multiple of three

Both directions annihilate `w`, so `u cross v=k w` for a nonzero rational
number `k`. Primitivity of `w` makes `k` integral: apply any integer Bezout
row `z.w=1` to the cross product.

Modulo three, the nonzero elements of the kernel with full coordinate
support are precisely two opposite vectors. Indeed, the three nonzero
values `w_i u_i` sum to zero and so are equal; the same holds for `v`.
Thus `u,v` are parallel modulo three and their cross product vanishes.
Since every coordinate of `w` is a ternary unit, `3|k`, whence `|k|>=3`.
The third component gives the load-bearing height inequality

```text
3c <= |u_1 v_2-u_2 v_1| <=2 M_u M_v,
M_u M_v>=3c/2.                                         (4)
```

No orientation sign is suppressed in the identity `u cross v=k w`; only
its absolute value enters the bound. Both raw sign orbits remain counted.

## 3. Sum the two complete ordinal lists

For either direction `d`, set

```text
B_d=min_i 3(sum(w)-w_i)/(14|d_i|),
K_d=strict_floor(B_d),
N_d=2(K_d-floor(K_d/3)).                                (5)
```

Every multiplier `+/-k d` with `1<=k<B_d` and `3` not dividing `k` is live,
and these are all carriers on that direction. The two lists are disjoint.
For `M_d=max_i|d_i|`, the narrow-coordinate and three-block bounds give

```text
B_d<3c/(7M_d),
N_d<4B_d/3+4/3.
```

Each summand in every network projection is at most `q=3/(7c)`. Therefore

```text
E_i <(12/49)(1/M_u+1/M_v)+8/(7c),            every i.   (6)
```

The reciprocal estimate is elementary and exact in its logical direction:
`(M_u-7)(M_v-7)>=0` implies

```text
1/M_u+1/M_v <=1/7+7/(M_u M_v)
             <=1/7+14/(3c).                            (7)
```

Combining (6)--(7) yields the all-height envelope

```text
E_i <12/343+16/(7c).                                   (8)
```

For `c>=55`, its right side is at most

```text
1444/18865 <1470/18865=6/77.                            (9)
```

Thus the complete remaining proof obligation is the finite two-ray head
`c<55`. This is a proved reduction, not an extrapolation of the census.

## 4. Exact finite head and hostile boundaries

The full universe of primitive sorted distinct positive odd ternary-unit
triples with `c<55` contains `814` rows. Exactly `192` have two primitive
live directions. Each of their three projection sums is strictly below
`6/77`. The maximum across all three coordinates and all these rows is

```text
w=(11,23,29),
directions (8,5,-7), (11,-4,-1),
(E_1,E_2,E_3)=(2/203,114/2233,1634/51359),
max E_i=114/2233 <6/77.                                (10)
```

This maximum is only asserted for the finite head; no sharp global maximum
is inferred from the loose analytic tail. All larger heights are strict by
(8)--(9), so the global theorem has no equality case.

The first multi-ray hostile `(17,23,25)` now satisfies the theorem because
its two directions are retained. The dense `A_2` carrier circuit at
`(19,23,29)` has three primitive directions, so applying (6) to only two
would omit raw mass. The theorem must not be read as a two-dimensional-span
result: every raw carrier lattice lies in a two-dimensional plane, while
the hypothesis counts its distinct unoriented primitive directions.

## 5. Reproducible verification

The primary script is
[lrc14_two_ray_overnight_hexagon_sep05.py](../../04-computation/lrc14_two_ray_overnight_hexagon_sep05.py).
It solves one modular congruence per coordinate row to enumerate the full
support. Before applying the two-ray filter, it compares every dictionary
with the independent integer-box implementation in the companion verifier.
For each of the 192 proof rows it checks the primitive coefficients, complete
strict multiplier lists, exact signed determinant multiplier, reciprocal
bound, and every actual projection.

The separately written literal-sheet engine from the companion verifier is
also used to reconstruct the three network projections on every proof row.
This is an explicit shared verification dependency, not a claim that the
two scripts were written in independent sessions. The kernel-row and kernel-
box implementations themselves are different, and the physical-sheet path
uses neither. Tail controls include `(1,11,55)`, `(1,5,101)`, and
`(5,49,251)`.

```bash
python3 -B 04-computation/lrc14_two_ray_overnight_hexagon_sep05.py
python3 -B -O 04-computation/lrc14_two_ray_overnight_hexagon_sep05.py
```

The finite head is part of the proof. Tail controls are supplementary tests.
All exact checks remain active under optimization. All results are confined
to the stated local network consumer, with three-or-more-ray support,
entry, synchronization, and LRC(14) remaining open.

Normal and optimized output streams are byte-identical. The primary records
`6,370` explicit checks and the imported literal-sheet control records `195`.
The companion script's raw hash is frozen in its note above. This script and
output have raw-byte SHA256:

```text
source 07cb6c797679335152eb850071456bf64d8fc7f1fc7da176b8f5f00790684d26
output 4b387916df346e1b166aee081d78897bae5c089153a5f068d1b44e0ac817ea1c
```
