# Variable: W(T,r) — W-polynomial

**Symbol:** `W(T,r)`
**Type:** polynomial in r with integer coefficients (only even powers at odd n)
**Defined in:** opus-S35c11 skeleton

## Definition
W(T,r) = sum over Hamiltonian paths P of prod_{edges e in P} (r + s_e)

where s_e = A_e - 1/2 is the centered edge variable (+1/2 for forward, -1/2 for backward).

## Key properties

| Property | Status | Source |
|----------|--------|--------|
| Only even powers of r at odd n | PROVED | THM-L (complement invariance) |
| W(T,-r) = (-1)^{n-1} W(T^op,r) | PROVED | THM-L |
| W(T) = W(T^op) for odd n | COROLLARY | of THM-L |
| W(T,r) = (r-1/2)^{n-1} * F(T, (2r+1)/(2r-1)) | PROVED | THM-K (Mobius transform) |
| W(T, 1/2) = H(T) / 2^{n-1} | by definition | evaluation |

## Fourier decomposition
W(T,r) = sum_k D_k where D_k is homogeneous degree 2k in edge variables s_e.

The D_k hierarchy:
- D_0, D_1: universal (same for all tournaments with same n)
- D_2, D_3: linear in t_3
- D_4+: depend on t_5, bc33, higher cycle invariants

## Relationship to F-polynomial (THM-K)
W(T,r) = (r - 1/2)^{n-1} * F(T, (2r+1)/(2r-1))

This Mobius transform connects W (which is even in r at odd n) to F.

**CORRECTION (opus-S37):** F(T,x) is NOT palindromic in general. The correct relation is F_k(T) = F_{n-1-k}(T^op) (complement duality). F is palindromic only for self-complementary tournaments.

## Values
W(T, 1/2) = H/2^{n-1} for any tournament.
W(T, 0) at odd n encodes the "constant term" c_0.

## Perpendicularity phenomenon
Var(c_0)/Var(H) decays exponentially:
| n | Var(c_0)/Var(H) |
|---|-----------------|
| 3 | 1.0 |
| 5 | 0.057 |
| 7 | 0.001 |

c_0 is identically 0 at even n.

## Tags
#W-polynomial #Fourier #complement-invariance #Mobius #perpendicularity
