# Three through six live directions: a quadratic-height network reduction

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED (root).**
This note reserves no theorem ID.
The target is the local network consumer, not chart entry or LRC(14).

Let `w=(a,b,c)` be primitive, sorted, distinct, positive, odd, and nonzero
modulo three. Let `Lambda(w)` be the **complete raw live carrier support** and
`E_1,E_2,E_3` the projections of
[THM-4414 — six-separated contact capacity collapse](../../01-canon/theorems/THM-4414-lrc14-six-separated-contact-capacity-collapse.md).
A direction identifies a primitive vector with its negative; all its live
positive and negative integer multiples remain separate carriers. Write `r`
for the number of directions.

The analytic result is

```text
r>=2, c>=33  ==>  every E_i < max(A_r(c), B_r(c)),

A_r(c)=(12r/49)sqrt(2/(3c))+4r/(7c),
B_r(c)=12/343+(12r-8)/(7c).                             (1)
```

In particular, at **any** finite number of directions,

```text
r>=2 and c>=7r^2+13r  ==>  every E_i<6/77.              (2)
```

The exact finite heads below complete the all-height conclusion

```text
r in {3,4,5,6}  ==>  every E_i<6/77.                    (3)
```

The complete raw support with seven or more directions is still unrestricted
below its quadratic cutoff. Even a universal network bound would leave chart
entry, synchronization, and LRC(14) open.

## Inheritance, connection, and hostile probe

The closest proved mechanisms are THM-4414's exact raw projection formula and
the short-relation obstruction in
[THM-4386 — canonical component relation and zero-defect incidence, Lemma 2.1](../../01-canon/theorems/THM-4386-lrc14-canonical-component-relation-and-zero-defect-incidence.md).
The [two-ray note](lrc14_two_ray_overnight_hexagon_sep05.md) supplies the
pairwise determinant idea.
[THM-4428](../../01-canon/theorems/THM-4428-lrc14-two-direction-network-closure-and-sharp-one-direction-gap.md)
was **RESERVED at inheritance**; incoming commit `4d1ad2a390` subsequently
promotes it. All needed determinant and multiplier arguments are nevertheless
repeated below, so the proof has no dependency on the reservation's former
status. Incoming
[THM-4431 — colored lattice basis and three-direction LRC network closure](../../01-canon/theorems/THM-4431-colored-lattice-basis-and-three-direction-lrc-network-closure.md)
separately reserves the concurrently obtained three-direction closure and its
two circuit types. The three-direction result
is convergent work; this note additionally proves the four-, five-, and
six-direction cases and the general quadratic-height reduction. A reservation
alone is not used as a proved dependency or a novelty claim.

The canonical dense hostile is `(19,23,29)` from
[THM-4422 — projection deficit and Beatty-row reduction](../../01-canon/theorems/THM-4422-lrc14-projection-deficit-and-beatty-row-reduction.md).
It has three directions and defeats the automatic cardinality gate, while
its actual projections are small. The corrected near miss is to retain a
minimal pair of directions and discard the other raw lists. The least-used
sidecar is the mod-three factor in **every pair's signed determinant**.

The live board is:

| Lane | Object and question | Operation and decisive test |
|---|---|---|
| Anchor | Complete multi-direction network support | Sum all raw ordinal lists, then test every actual projection |
| Niche | Positive coefficient maxima with pairwise product floors | Optimize reciprocal sum using its smallest variable |
| Wildcard | Three direction sets without an additive root circuit | Test the incoming `(41,47,49)` arithmetic-progression triple |
| Sidecar | Owner coset inside the rank-two kernel lattice | Retain oriented determinants before taking absolute bounds |

The source-to-target map sends each complete raw list to its primitive maximum
coefficient `M_d` and retains all direction multiplicities separately. It
preserves an upper bound on each of the three network projections. It
forgets actual slopes, boundary deficits, and directional signs **only after**
the signed determinant identity is established; the original raw list is the
restoration sidecar. The cheapest hostile test is to reconstruct each
direction's entire strict multiplier interval and all literal sheet graphs.

No additive circuit or empty-hexagon theorem is assumed. In particular,
incoming root work supplies `(41,47,49)` with directions

```text
(11,5,-14), (14,-7,-5), (17,-19,4).
```

These form an arithmetic progression, but no choice of signs yields a relation
`u+v=z`. The displayed three vectors are checked in the reproduction script.
Their determinant interaction is sufficient for `(1)` without that circuit.
The concept-board update is that a root-circuit classification is unnecessary
for bounded direction counts; it may still matter in the unbounded remainder.

## 1. Every direction is long, and every pair pays a determinant

For any primitive direction `v`, all three coordinates are units modulo three
and `||v||_1` is even because the speeds are odd. If `||v||_1<=14`, THM-4386's
zero-defect lemma forces every live carrier to be parallel to `v`. This
contradicts `r>=2`. Hence `||v||_1>=16`. No ternary-unit coordinate can have
magnitude six; if every magnitude were at most five, parity would force norm
at most fourteen. Therefore

```text
M_v=max_i|v_i|>=7.                                     (4)
```

For two distinct directions `u,v`, both annihilate `w`, so

```text
u cross v = k w.
```

Since `w` is primitive, any integer Bezout vector `z.w=1` gives
`k=z.(u cross v)` integral. It is nonzero because the directions are distinct.
Modulo three the three quantities `w_i u_i` are nonzero and sum to zero, so
they are all equal; the same applies to `v`. Thus `u,v` are parallel modulo
three and their cross product vanishes modulo three. Since all speeds are
ternary units, `k` is a nonzero multiple of three. Taking the third coordinate,

```text
3c <= |u_1v_2-u_2v_1| <=2 M_u M_v.

M_u M_v >= P:=3c/2, for every pair of directions.       (5)
```

This is the step that prevents several small independent directions from
simultaneously carrying large raw multiplier lists.

## 2. A sharp elementary reciprocal lemma for the relaxed variables

Let `r>=2`, `L>0`, `P>=L^2`, and let positive reals `M_1,...,M_r` satisfy
`M_j>=L` and `M_j M_k>=P` whenever `j!=k`. Then

```text
sum_j 1/M_j <= max(r/sqrt(P), 1/L+(r-1)L/P).            (6)
```

To prove this, let `m=min_j M_j`. If `m>=sqrt(P)`, the first branch applies.
Otherwise each other `M_j>=P/m`, whence

```text
sum_j 1/M_j <= 1/m+(r-1)m/P.
```

The right side is convex on `[L,sqrt(P)]`, so its maximum occurs at an
endpoint. The endpoint values are precisely the two branches of `(6)`.
Both branches are attained by admissible **relaxed real variables**:
all `M_j=sqrt(P)`, or one `M_j=L` and the remaining `M_j=P/L`.
This sharpness claim is for the abstract product-constrained optimization;
neither vector realization nor a sharp network constant is claimed.

Apply `(6)` with `L=7` and `P=3c/2`. The condition `P>=49` is supplied by
`c>=33`.

## 3. Restore and sum every raw multiplier

For a direction `v`, put

```text
B_v=min_i 3(sum(w)-w_i)/(14|v_i|),
K_v=(numerator(B_v)-1)//denominator(B_v).
```

The complete list on this direction is exactly

```text
{+/-k v: 1<=k<=K_v, k not divisible by 3}.
```

The strict endpoint formula includes all positive roof conditions and excludes
zero-length boundary carriers. Its cardinality is

```text
N_v=2(K_v-floor(K_v/3)) <4(B_v+1)/3.
```

The coordinate with magnitude `M_v` has a speed-pair sum strictly below `2c`,
so `B_v<3c/(7M_v)`. Every term of every projection is at most `3/(7c)`.
Different primitive direction lists are disjoint. Consequently

```text
E_i <=(3/(7c))sum_v N_v
    <(12/49)sum_v 1/M_v+4r/(7c).                       (7)
```

Combining `(6)` and `(7)` gives `(1)`. The strict inequality comes from the
strict roofs and ordinal count, irrespective of possible equality in `(6)`.

## 4. Exact cutoffs and an explicit general threshold

Both branches in `(1)` decrease with `c`. At an integer `c>=33`, their values
are strictly below `6/77` exactly when

```text
6c>44r,
49(6c-44r)^2 >11616r^2 c,
162c >539(12r-8).                                      (8)
```

The first line preserves the sign before squaring the square-root branch.
There is no floating-point comparison in the verifier. The first integer
cutoffs are

| Directions `r` | First integer `c` passing `(8)` |
|---:|---:|
| 3 | 99 |
| 4 | 159 |
| 5 | 233 |
| 6 | 319 |

To obtain a simple formula valid for every `r>=2`, set `c=r(7r+13)`. Then

```text
49(6c-44r)^2-11616r^2c
 =4r^2(1281r^2-2766r+14161),

162c-539(12r-8)
 =2(567r^2-2181r+2156).
```

Writing `r=s+2`, the two parenthesized polynomials become

```text
1281s^2+2358s+13753,
567s^2+87s+62.
```

Both are positive for `s>=0`; also `6c>44r`. This proves `(2)` by
monotonicity. The implication is sufficient, not necessary. It does not bound
the number of directions of an arbitrary triple.

## 5. Finite proof heads, controls, and reproduction

The exact universe is every primitive sorted distinct positive odd ternary-unit
triple below height `319`. For a row with exactly `r in {3,4,5,6}` directions,
the proof includes it precisely when its height is below the corresponding
cutoff in the table. No pairwise-coprimality or nonresonance filter is added.
Every other row is excluded only because it lies outside the statement or in
its proved analytic tail.

The verifier checks complete support by two separately coded coordinate-row
enumerations **before** applying the direction-count filter. At every height
below `99`, an independent full integer-box loop additionally checks that
support. On **every selected proof row**, it computes the actual three network
projections both from complete raw carriers and by a literal six-sheet interval
graph engine. That literal engine and the low integer-box loop are imported
from the prior [one-ray verifier](../../04-computation/lrc14_one_ray_overnight_hexagon_sep05.py);
this is a declared code dependency, not a claim of independent authorship.

The source is
[lrc14_empty_core_three_ray_sep06.py](../../04-computation/lrc14_empty_core_three_ray_sep06.py),
with matching [output](lrc14_empty_core_three_ray_sep06.out).

```bash
python3 -B 04-computation/lrc14_empty_core_three_ray_sep06.py
python3 -B -O 04-computation/lrc14_empty_core_three_ray_sep06.py
```

The complete support comparison covers `190,951` eligible triples. The exact
finite heads are:

| Directions | Eligible triples below cutoff | Selected proof rows | Largest of all three projections in the head | Speed triple |
|---:|---:|---:|---:|---|
| 3 | 5,409 | 1,791 | `18/301` | `(5,37,43)` |
| 4 | 23,210 | 5,004 | `30/581` | `(1,67,83)` |
| 5 | 72,458 | 13,656 | `36/623` | `(7,77,89)` |
| 6 | 190,951 | 22,971 | `48/1001` | `(13,121,143)` |

All `43,422` selected proof rows have every projection strictly below `6/77`.
The sole failure of `N<=2c/11` within these four declared heads is the inherited
`(19,23,29)`. The finite leaders are not claimed as sharp all-height maxima;
the analytic tail envelope is looser than their exact values. Equality in the
stated strict theorem is impossible.

The three-ray control `(107,127,149)` has **fourteen** raw carriers, showing
that three primitive directions must not be confused with six raw carriers.
Controls also include an actual tail row for each of the four direction
counts, with the six-direction boundary row `(1,5,319)` giving every projection
`36/2233`.

Both normal and optimized runs pass and their complete output streams are
byte-identical: `4,682,652` explicit checks plus `43,428` literal-sheet checks.
The root agent independently audited the signed determinant integrality and
mod-three factor, `M>=7`, the convex reciprocal relaxation, strict raw
multiplier count, radical sign conditions, both polynomial identities, the
finite universes and enumerators, and the literal consumers. Its separate
height-499 integer census supports the finite leaders; that larger census is
supplementary evidence rather than a dependency of the finite-head proof.

Raw LF-byte SHA256:

```text
source       dc0314ecac5733c4854674e5f036fcfec9ff55a4901ac6b85c930ec9c7d9e652
output       8dd3b428ad50c819d5d03aca2b0a99d65164709847a0143a5290f971589f045b
literal dependency
             6b41a879700632aa934650f27dafe9d076c051ddcee3262fabc07556a7aaf875
```

The finite heads are load-bearing proof obligations; high-tail examples and
abstract-polynomial samples are controls. The analytic polynomial identities
above, not the samples, prove the general threshold.

## 6. Audit connection to the concurrent relation-slice compiler

The independent audit script
[lrc14_relation_slice_audit_empty_core_three_ray_sep06.py](../../04-computation/lrc14_relation_slice_audit_empty_core_three_ray_sep06.py)
checks the concurrent
[short-relation slice certificate](lrc14_empty_core_certificate_sep06.md).
Its different representation uses rational two-dimensional polygon clipping,
the third cross-product coordinate, every signed coordinate permutation,
complete finite-field owner/defect fibers, direct enumeration of all eligible
speed triples through the largest coefficient cutoff, and full integer-box
carriers. It reconstructs all `73` coefficient patterns, `4,083` head
memberships, and `1,944` unique head rows, including exact physical mass and
every network projection. Its complete semantic digest agrees with the
producer. See its [audit output](lrc14_relation_slice_audit_empty_core_three_ray_sep06.out).

That connection changes the board: bounded relation length can close families
with arbitrarily many directions, while the reciprocal lemma closes bounded
direction counts without any prescribed coefficient shell. They are
complementary sufficient hypotheses. Neither supplies the other hypothesis
for an arbitrary triple.

The subsequent [actual-zero-coordinate addendum](lrc14_pair_relation_empty_core_certificate_sep06.md)
was also independently audited with
[a separate polygon and pair-ratio verifier](../../04-computation/lrc14_pair_slice_audit_empty_core_three_ray_sep06.py).
All eleven new patterns, `209` memberships, and `182` exact heads reproduce
the producer's complete consumer digest; normal and optimized outputs match.
This removes the full-support restriction from the norm-twenty relation
certificate. The new [audit output](lrc14_pair_slice_audit_empty_core_three_ray_sep06.out)
records the bounded result. No longer-height census was used for that audit.
