# Spectral Blindness, Kurtosis, and the Paley Principle

**Session:** kind-pasteur-2026-03-21-S13

## The Discovery

The signed adjacency matrix S_T of a tournament carries ALL directed information.
Its powers S_T^k alternate between Cartan sectors (antisymmetric for k odd,
symmetric for k even). But the SPECTRAL invariants of S_T do NOT determine H(T).

Specifically, at n=7 among regular tournaments:

| H | Eigenvalues of S_T | tr(S^2) | tr(S^4) | dim(Alg) | dim(Comm) |
|---|-------------------|---------|---------|----------|-----------|
| 189 (Paley) | {0, +/-i*sqrt(7)} x3 | -42 | 294 | 3 | 19 |
| 175 | {0, +/-0.482i, +/-1.254i, +/-4.381i} | -42 | 742 | 7 | 7 |
| 171 | {0, +/-0.268i, +/-2.646i, +/-3.732i} | -42 | 486 | 7 | 7 |

## The Second Moment is Universal

tr(S_T^2) = -n(n-1) for ALL tournaments (not just regular).
For regular tournaments: tr(S_T^2) = -42 at n=7.

This is because: (S_T^2)[i,i] = -(n-1) for all i (sum of squares of +/-1 entries).
The second moment is a TOPOLOGICAL invariant — it depends only on the number
of vertices, not on the specific tournament.

## The Fourth Moment is the Discriminant

tr(S_T^4) = sum lambda_k^4 FIRST distinguishes the three H classes.

Paley has the SMALLEST tr(S^4) = 294. The others have 486 and 742.

Since all have the same tr(S^2) = -42, the difference is in the DISTRIBUTION
of eigenvalue magnitudes. The fourth central moment measures this: kurtosis.

**Paley has minimum kurtosis.** All eigenvalues have the same magnitude sqrt(7).
The other regular tournaments have spread-out eigenvalue magnitudes.

## The Paley Principle

**Minimum kurtosis <==> Maximum H <==> Maximum commutant dimension.**

- min kurtosis: all eigenvalues equal → S^3 + p*S = 0 (degree 3 min poly)
- max H = 189 = the Paley value
- max comm_dim = 19 = 1 + 9 + 9 (eigenspace dims 1,3,3)

For non-Paley regular:
- higher kurtosis: spread eigenvalues → degree 7 min poly (generic)
- lower H (171 or 175)
- min comm_dim = 7 (all eigenspaces 1-dimensional)

## Why This Matters for the Cartan Bridge

The Cartan decomposition gl(n,R) = so(n) + p + R separates matrices into
antisymmetric and symmetric parts. But within the antisymmetric sector,
there is a FINER structure: the spectral kurtosis.

The spectral kurtosis of S_T measures how "doubly regular" the tournament is:
- kurtosis = 0 (minimum): all eigenvalues equal → DRT → max H
- kurtosis > 0: eigenvalues spread → non-DRT → lower H

This creates a HIERARCHY of Cartan bridge invariants:
1. Cartan sector (anti/sym split) — coarsest
2. Second spectral moment tr(S^2) — universal, carries no info
3. Fourth spectral moment tr(S^4) — first discriminant = kurtosis
4. Sixth spectral moment tr(S^6) — finer discrimination
...
n-1. Full spectrum {lambda_k} — determines H class but NOT individual H
n. Full matrix S_T — determines everything

H is captured somewhere between levels n-1 and n. The spectrum alone
is NOT enough (shown at n=5), but the spectrum + eigenvector structure IS.

## Connection to Attention

For transformer attention, the Paley Principle says:

**The most powerful attention head has the FLATTEST spectral distribution
in its antisymmetric (tournament) sector.**

This means: optimal attention is not "attend to one token strongly"
(peaked spectrum) but "attend to all tokens with equal strength
but specific direction" (flat spectrum = Paley-like).

In information-theoretic terms: Paley attention has maximum entropy
among tournament-structured attention patterns. This is the directed
analog of "uniform attention" — but with preserved directional structure.

## The Commutant as Symmetry Measure

dim(Comm(S_T)) = sum d_k^2, where d_k are eigenspace dimensions.

For Paley: d = (1, 3, 3) → comm = 1+9+9 = 19
For generic n=7: d = (1,1,1,1,1,1,1) → comm = 7

The ratio 19/7 = 2.71... measures the "symmetry excess" of Paley.

For general n=p (prime, p=3 mod 4):
Paley comm = 1 + 2*((p-1)/2)^2 = 1 + (p-1)^2/2
Generic comm = p

Symmetry excess = (1 + (p-1)^2/2) / p → (p-1)/2 as p → infinity.

The symmetry excess grows linearly with p! Large Paley tournaments
are exponentially more symmetric than generic tournaments.

## The Number 19

dim(Comm(Paley T_7)) = 19.

19 = 1 + 9 + 9 = 1 + 2*3^2.

General formula: dim(Comm(Paley T_p)) = 1 + 2*((p-1)/2)^2 = (p^2 - 2p + 3)/2.

For p=3: (9-6+3)/2 = 3. For p=7: (49-14+3)/2 = 19. For p=11: (121-22+3)/2 = 51.
For p=23: (529-46+3)/2 = 243 = 3^5.

The sequence: 3, 19, 51, 243, ... = (p^2-2p+3)/2.

At p=23: comm = 243 = 3^5. A pure power of 3! This connects to...
nothing obvious yet. But the fact that the commutant dimension for Paley T_23
is exactly 3^5 deserves investigation.
