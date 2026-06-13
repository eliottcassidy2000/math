# Variable: c_k — Fourier Coefficients of W-polynomial

**Symbol:** `c_k` (also `D_k` for degree-2k component)
**Type:** rational or integer, depending on normalization
**Defined in:** THM-055, opus-S11b series

## Definition
W(T,r) = sum_{k=0}^{floor((n-1)/2)} c_k * r^{n-1-2k}

At odd n, only even powers of r appear, so k ranges 0...(n-1)/2.

D_k = degree-2k homogeneous component of W in edge variables s_e.

## Hierarchy (THM-055)

| Level | Depends on | Status |
|-------|-----------|--------|
| D_0 | nothing (universal constant) | PROVED all n |
| D_1 | nothing (universal constant) | PROVED all n |
| D_2 | t_3 only | PROVED all n |
| D_3 | t_3 only | PROVED all n |
| D_4 | t_3, t_5, bc33 at n=7 | PROVED n=7; FAILS to generalize n=9 |
| D_k (k>4) | higher cycle invariants | OPEN |

## Exact formulas at n=7
- c_2 = 24*bc33 - 60*t_3 + 12*t_5 + 231
- c_0 = 2*t_3 - t_5 + 2*t_7 - 2*bc33 + 253/4

## Exact formulas at n=9
- c_2 = 462*t_3 - 60*t_5 + 12*t_7 - 120*bc33 + 24*bc35_w + 48*a3 - 2640

## Perpendicularity
Var(c_0)/Var(H) decays exponentially with n (see [W-polynomial.md](W-polynomial.md)).
At even n, c_0 = 0 identically.

## Key insight (CRITICAL)
Fourier degree-4 at n=9 has dimension >> 200, not 2 like at n=7.
The "two-type" decomposition does NOT generalize. See HYP-113 in hypotheses.

## Tags
#Fourier #W-polynomial #hierarchy #coefficients #perpendicularity
