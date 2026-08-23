# Vertical packets, cover degree, and rational target operations

**Status:** **PROVED COROLLARY** for the THM-3774 quadratic row over `C`;
**FINITE-EXACT** for the 59-check companion; **SYNTHESIS** for the invariant
hierarchy.  `JC(2)` remains **OPEN**.

The incoming audit separates two invoices that recent rational-exact
Jacobian families had begun to blur:

| operation class | exact invariant | consequence |
|---|---|---|
| keep `Q` fixed and add `H(Q)` | vertical principal-part packet modulo the scalar diagonal | THM-3770 is necessary and sufficient |
| birationally change both targets | `[k(X,T):k(P,Q)]` | field degree is preserved |
| naturally oriented short shear word | completed-DVR leading coefficients | THM-3776 excludes only its typed length-two/three words |

These invariants are independent.  A residue mismatch can survive every
fixed-target correction while an unrestricted rational target inverse still
exists; conversely, target-cover degree two blocks every birational target
escape over `C` even when THM-3776's two-nonzero-residue hypothesis fails.

## 1. The quadratic row of THM-3774

For nonzero `a,b,c,d`, put

```text
A=a+bxy,          B=c+d x^3 A,
U=xAB,             P=(2B-c)/x.                       (1)
```

Exact differentiation and a cleared Bezout identity give

```text
J(P,U)=-bc^2,
(x(2B-c)_x-(2B-c))U_y-x(2B-c)_yU_x=-bc^2x^2.       (2)
```

Hence `U` is smooth.  Its zero fibre consists of the three disjoint reduced
components `x=0`, `A=0`, and `B=0`; the `U`-principal coefficient packet of
`P` is

```text
(ac^2,0,0).                                          (3)
```

The conic identities recover the whole source field and give

```text
(P^2-4dU)x^2=c^2,          [k(x,y):k(P,U)]=2.        (4)
```

Thus THM-3770 forbids a polynomial mate in `P+k(U)`.  Over `C`, the stronger
degree argument applies: any birational rational symplectic target change
preserves the index two in (4).  If both transformed outputs were polynomial,
they would be a degree-two planar Keller map, contrary to THM-3773/THM-1330.
This proves the new THM-3774 `m=1` corollary for all birational target maps
and hence all finite alternating target-shear words, in either orientation.
The `C` and `m=1` restrictions are load-bearing.

The exponent three is also exact, not sampled.  With
`B_n=c+d x^nA`, `U_n=xAB_n`, and `P_n=(2B_n-c)/x`,

```text
J(P_n,U_n)=-b[c^2+2(3-n)B_n(B_n-c)].                 (5)
```

The coefficient of `y^2` in `B_n(B_n-c)` is
`b^2d^2x^(2n+2)`, so for the genuine profile `d!=0` the response is a
nonzero constant exactly at `n=3`.  This is the `m=1` specialization of
THM-3774's all-degree exponent `2m+1`.

## 2. Unequal residues are not an all-map invariant

The smallest nonlinear THM-3771 control has

```text
U=x,       W=x+3xt+1,       q=U(W-2),       p=-W/(3q). (6)
```

Its two zero-fibre coefficients are the distinct nonzero values
`-1/3` and `-2/3`, so it lies inside THM-3776's natural-orientation
hypotheses.  Nevertheless its exact rational symplectic inverse is

```text
x=-q/(3pq+2),
t=(3pq+1)(3pq+2)/(3q)-1/3.                           (7)
```

Therefore unequal vertical coefficients cannot obstruct arbitrary rational
two-output target changes.  Equation (7) is a hostile boundary against any
all-map strengthening of THM-3776; it does not produce, or claim, a
three-shear factorization.

## 3. Two quadratic involutions that must not be identified

THM-3775's generic-fibre sheet involution fixes `(Q,z)`, flips its sheet
coordinate, and anti-commutes with `J(-,Q)`.  The THM-3774 target-cover deck
fixes both `(P,U)` and commutes with the uniquely extended derivation
`J(-,U)`.  The common word quadratic is not an object-level bridge.  The
packet (3) is governed by THM-3770's component quotient, not by importing
THM-3775's anti-invariant sheet representation.

Likewise THM-3776 remains correctly typed: its `q-p` and `q-p-q` conclusions
assume that `q` is the common uniformizer.  The variable-exchanged dual
requires a separate hypothesis with `p` as common uniformizer; merely
starting the original data in the opposite orientation is invalid.

## 4. Exact replay and scope

The independent companion checks 59 exact identities: response and
smoothness, all three components and packet entries, conic recovery and
quadratic deck, exponent uniqueness through hostile values, sharp
degenerations, bounded polynomial-mate ranks, and the unrestricted inverse
(7).  Normal, optimized, and frozen outputs byte-match.

```powershell
python3 -B 04-computation/jc2_vertical_spectrum_cover_degree_target_word_triage_20260823.py
python3 -B -O 04-computation/jc2_vertical_spectrum_cover_degree_target_word_triage_20260823.py
```

This synthesis neither extends the quadratic Galois exclusion to the
non-Galois rows `m>=2` nor supplies a synchronized nonbirational
log-canonical cover.  Those construction boundaries, arbitrary polynomial
Keller pairs, and `JC(2)` remain open.
