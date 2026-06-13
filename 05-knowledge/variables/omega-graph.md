# Variable: Omega(T) — Odd Cycle Conflict Graph

**Symbol:** `Omega(T)` or `Omega_k(T)` for k-cycle only version
**Type:** undirected graph
**Defined in:** `01-canon/definitions.md`; THM-002

## Definition
Omega(T) has:
- Vertices: all directed odd cycles in T
- Edges: two cycles share at least one vertex

Omega_3(T) restricts to 3-cycles only.

## Key equations
- **OCF (THM-002):** H(T) = I(Omega(T), 2)
- I(G,x) = independence polynomial = sum over independent sets S of x^|S|

## Structural properties

| Property | Status | Range verified | Source |
|----------|--------|---------------|--------|
| Quasi-regular | TRUE | n=5-20 | INV-041 |
| Claw-free | TRUE n<=8, FALSE n>=9 | exhaustive n<=8 | INV-032 |
| Perfect | TRUE n<=7, FALSE n>=8 | exhaustive | INV-032 |
| S_{1,1,1}-free | TRUE n<=11, FALSE n>=12 | exhaustive | INV-032 |
| Real-rooted I(Omega,x) | TRUE n<=8, FALSE n>=9 | exhaustive n<=8 | THM-025 |
| Line graph | FALSE n>=6 (K5-e found) | | INV-032 |

## Size growth
| n | avg |V(Omega_3)| | max |V(Omega_3)| |
|---|-------------------|-------------------|
| 3 | 0.5 | 1 |
| 5 | ~3 | 8 |
| 7 | ~28 | 35 |

## Related
- [independence-poly.md](independence-poly.md) — I(Omega(T), 2) = H(T)
- [hamiltonian-paths.md](hamiltonian-paths.md) — connected via OCF
- [mu-C.md](mu-C.md) — mu(C) involves restricted Omega

## Tags
#Omega #conflict-graph #OCF #independence-polynomial #claw-free #real-roots
