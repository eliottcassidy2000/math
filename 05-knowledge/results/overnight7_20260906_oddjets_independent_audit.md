# Independent audit of the odd-prime three-node three-jet Smith law

**Status: INDEPENDENT PROOF AUDIT + COMPLETE EXACT CERTIFICATE PASS.**
This is an independent audit
of the primary `overnight7_20260906_oddjets` report, producer and polynomial
manifest. It does not claim a second discovery. All files remain outside
the repository; no Git mutation was made.

## 1. Exact target, normalization and dependency

The target is the full nine-entry **p-primary Smith exponent list** of the
complete Hasse observer of orders zero, one and two at three distinct
integer nodes, acting on polynomials of degree below nine. The theorem
requires an odd prime. The dyadic unit-sensitive law is not extended.

Translation of a node to zero is an integral unimodular coefficient change.
Its three Hasse rows split off `I_3`. Over `Z_p`, division by a unit in the
node differences is also a unit diagonal source/target change. The residual
six-by-six matrix therefore has rows `(node,order)=(1,0),(1,1),(1,2),
(a,0),(a,1),(a,2)` and columns `j=3,...,8`, with entries

```text
binom(j,r) p^(e(j-r)) a^(j-r)
```

on the `a` rows, and the same expression without `a` on the unit rows.
The close shape is `a=p^d u`, `d>=1`, `u` a p-adic unit; the pairwise
depths are `(e,e,e+d)`. The equilateral shape has both `a` and `1-a` units,
with all pair depths equal to `e`. These exhaust three-node odd-prime
metrics. Unit-normalized coordinates may be rational p-adic integers;
every coefficient argument below applies to those as well as integers.

The already proved dependency
[THM-4443, arbitrary-jet precision and dyadic unit boundary](../../01-canon/theorems/THM-4443-arbitrary-jet-precision-and-dyadic-unit-boundary.md)
supplies the attained largest inverse denominator. It is used only for the
largest exponent, not as a substitute for the intermediate determinantal
ideals. The primary respects [MISTAKE-547](../../01-canon/MISTAKES.md)'s
observer boundary and keeps the dyadic witnesses `(0,4,8)` and `(0,4,12)`.

## 2. Complete residual minors and the close-case formulas

For residual row set `I` and degree set `J`, factor the minor as
`p^(eW) F_(I,J)(a)`, where

```text
W=sum J-sum_(i in I) derivative_order(i).
```

All residual minors of ranks one through four were reconstructed
independently: `36+225+400+225=886` polynomials with 2,623 nonzero monomials.
This is the entire proof universe, not a selected minor bank. A coefficient
`c a^b` in the close case has valuation at least `eW+bd+v_p(c)`.
Writing `d=1+D` gives the three nonnegative-coordinate costs
`(W,b,b+v_p(c))`. The complete lower Pareto sets match the primary
certificate exactly at `p=3`, `p=5`, and with the coefficient-valuation
floor zero used for every `p>=7`.

Let `Delta_r` denote the least valuation of a rank-`r` residual minor,
and put `Delta_0=0`. The audited close formulas are:

| rank | p=3 | p=5 | p>=7 |
|---|---|---|---|
| 1 | min(3e,e+1) | e | e |
| 2 | min(6e,4e+1,3e+d+2) | min(4e,3e+d) | min(4e,3e+d) |
| 3 | min(9e,7e+d+1) | min(9e,7e+d+1) | min(9e,7e+d) |
| 4 | min(13e+d+1,12e+4d+2) | min(13e+d,12e+4d+1) | min(13e+d,12e+4d) |
| 5 | 19e+4d+1 | 19e+4d | 19e+4d |
| 6 | 27e+9d | 27e+9d | 27e+9d |

The full Smith list is three zeros followed by the six successive
differences of these residual determinantal valuations.

Lower coefficient bounds alone would not prove these equalities. The
attaining minors are retained explicitly. The necessary low-weight factors
are `1`, `3`, `6`, `18a(a-1)`,
`30a(a-1)(2a-1)`, `3a(a-1)(5a^2-5a+1)`, and
`90a^4(a-1)^4`. When `a=p^d u`, all displayed factors other than powers
of `a` and the fixed integer coefficients are units. Thus their valuations
are exact at every close parameter, including arbitrarily deep unit lifts.
The residual determinant gives rank six; the largest exponent gives rank
five as audited in Section 4.

The additional quinary rank-three cost `8e+d` is redundant on the entire
real cone `e>=0,d>=1`, not merely at integer depths. Indeed

```text
8e+d >= [(9e)+(7e+d+1)]/2 >= min(9e,7e+d+1).
```

Dropping `d>=1` invalidates the simplification: `e=1/2,d=1/4` makes
`8e+d` strictly smaller than both remaining costs. The primary corrected
its preliminary integer-only description before freezing the theorem.

## 3. Equilateral lower bounds and complete attainment

At `e=0`, the three nodes have distinct residues, and the Hasse CRT
observer has unit determinant. Every Smith exponent is zero. This case
must remain separate from the positive-depth formulas.

For `e>=1`, the first four residual determinantal valuations are

```text
(Delta1,Delta2,Delta3,Delta4)
 =(e,3e,7e,12e)+beta,
beta=(1,2,2,2) at p=3;
     (0,0,1,1) at p=5;
     (0,0,0,0) at p>=7.
```

For `p=3`, equilateral normalization forces `a=2+3z`, `z in Z_3`.
All 8,051 nonzero coefficients after that substitution were independently
verified. If `s=(1,3,7,12)` at rank `r`, every coefficient satisfies
`W>=s_r` and `W+v_p(c)>=s_r+beta_r`. Writing `e=1+E` proves the lower
bound for every `E>=0`, and arbitrary p-integral `z` cannot decrease it.
The analogous original-coefficient inequalities hold at five and with the
zero coefficient-valuation floor for primes at least seven.

The attainment audit needs two rank-three witnesses:

```text
30a(a-1)(2a-1),
30a^4(a-1)(a-2),
```

both of weight seven. Their factors `a,a-1` are units. For `p!=3`,
the two remaining linear factors cannot both vanish modulo `p`, because
`(2a-1)-2(a-2)=3`. Thus the minimum extra valuation is zero. The fixed
coefficient 30 then contributes exactly one at `p=5` and zero at every
`p>=7`. At `p=3`, both linear factors are divisible by three; the same
identity shows that their minimum valuation is exactly one. Along with
`v_3(30)=1`, the exact rank-three intercept is two.

The weight-twelve minor `90a^4(a-1)^4` attains the rank-four intercept.
The rank-one and rank-two coefficients `3` and `18a(a-1)` attain theirs.
Thus there is no unproved generic-unit assumption at a tie.

For `e>=1`, the full exponent lists are therefore:

```text
p=3: (0,0,0,e+1,2e+1,4e,5e,7e-1,8e-1),
p=5: (0,0,0,e,2e,4e+1,5e,7e-1,8e),
p>=7:(0,0,0,e,2e,4e,5e,7e,8e).
```

## 4. Independent largest-loss derivation

At a node with differences `u,v`, the three reciprocal coefficients are

```text
a0=1/(uv)^3,
a1=-3(u+v)/(uv)^4,
a2=3(2u^2+3uv+2v^2)/(uv)^5.
```

At a close-pair node, the two valuations are `e,e+d`. For every odd prime
the quadratic numerator has the unique least-valuation term `2u^2`
after naming `u` the shallower difference. Its valuation is exactly `2e`.
The order-two candidate is `8e+5d-[p=3]`; it exceeds the first two
orders and every outsider candidate. Hence the close largest loss is
exactly that quantity.

In the equilateral case, after removing the common scale the three
quadratic numerators are

```text
N0=2a^2+3a+2,
N1=2a^2-7a+7,
Na=7a^2-7a+2.
```

Some numerator is a unit at every odd prime. For `p!=5`, simultaneous
zeros of the first two imply `a=1/2`, because their difference is
`5(2a-1)`. But `N0(1/2)=4`, impossible at an odd prime. At `p=5`, the
discriminant of `N0` is `-7=3 mod5`, a nonsquare, so `N0` has no residue
root. Thus the order-two maximum is `8e-[p=3]`. Including the attained
order-zero baseline gives

```text
L_equilateral=max(6e,8e-[p=3]).
```

The total determinant valuation is `9` times the sum of pair depths,
namely `27e+9d` in the close case and `27e` in the equilateral case.
Subtracting the largest exponent gives residual rank five. This is a
valid use of the largest factor and determinant after all lower ranks
have been independently controlled; it does not reconstruct a partition
from the determinant alone.

## 5. Independent exact engine and final review

The independent companion uses only Python's standard library. Each minor
polynomial is reconstructed by exact Bareiss integer determinants at
degree-bounded interpolation points. The degree bound is the sum of the
largest columns assigned to variable-node rows minus those rows' derivative
orders, at most eighteen. Finite differences recover its exact integer
coefficients. This is independent of the primary's permutation expansion.
Every reconstructed coefficient matches the primary JSON manifest.

Full matrix controls use p-local pivot elimination on literal nine-by-nine
Hasse matrices. The least-valuation pivot makes all elimination multipliers
p-integral, so the row and column operations are unimodular over `Z_p`.
The method uses neither SymPy nor an integer Smith routine. It matches the
full formulas on **177 fresh matrices**, including prime 13, signed unit
lifts, close and equilateral shapes, and outer depth zero. Two additional
dyadic controls recover largest losses 18 and 19 and stay outside the odd
theorem.

The companion has **24,396 explicit gates** when run against the primary
polynomial certificate. Its polynomial and literal semantic digests are

```text
polynomials f26413b91936dedafca11d18a67372cd43dd9ca92a031ccf17b91d74f6e344e1
literal 25d880af06691156da6a17553972d94d835ab29de4b481c362a03bc0e85603d0
audit source 6ba392a11651dfd36e0b336ff3763be87531ceb8bf12772fc66e54f0eb35ba77
audit output 19cfa766d56ebf7229d01306258d9e6e51aca1feb19f029cee77b737ce066a97
```

```text
python -B 04-computation/overnight7_20260906_oddjets_independent_audit.py --primary-certificate 05-knowledge/results/overnight7_20260906_oddjets.json
python -B -O 04-computation/overnight7_20260906_oddjets_independent_audit.py --primary-certificate 05-knowledge/results/overnight7_20260906_oddjets.json
```

The final primary report was read in full, including its local unit
normalization, complete coefficient universe, lower-bound proof, attaining
witnesses, largest-loss formulas, real-cone simplification and all scope
boundaries. I also read the incoming proved dyadic full-partition result;
it supports the combined conclusion that unit dependence in this precise
three-node uniform-three-jet problem is confined to the final two dyadic
factors. No correction or unresolved mathematical gap remains.

Normal and optimized independent audit outputs are byte-identical. The
frozen primary artifacts audited have hashes

```text
source b894f4860c1bdd824d4b2ce5ef2b443bdff2b73a684a06d5ac699c3afafa45b1
output 1c2b6d20a5a58d8d5a2c7d7565c8a058c3bac54b90d9c0a407785bbc351f5e85
JSON 05fc22845ec10078a2c9d568dea59d00bfdbdc8588d9e314e4134fcbe13e7d11
```
