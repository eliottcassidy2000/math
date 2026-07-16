---
id: THM-908
title: THE SEVEN-SHIFT RESIDUE-SIX SIEVE — averaging x+j/7 turns the exact signed K6 kernel into a finite affine-line ceiling on F7^5; 906 of 2801 projective residue directions close immediately, including every five-tuple of distinct nonzero residues with -F6 <= 25/343 < 0.097
status: PROVED identity and finite projective census; dependency-free verifier being banked
source: codex-2026-07-16-S18
depends_on: [THM-891, THM-903-reflection-frame-residue6]
related: [THM-904, THM-905, THM-906, THM-907, HYP-7064]
verification: exact F7^5 verifier in progress
---

# THM-908 — the seven-shift residue-six sieve

Let `e=(e_1,...,e_5)` be the five moving speeds, let
`s_i(x)=floor(7{e_i x})`, and let `L_6(s)=-K_6(M(s))`, extended by zero when
the inner missed set has size other than one or two as in `THM-891`.  Put
`r_i=e_i mod 7`.  For `u in [0,1)` and `j in F_7`, define

```text
z_i(u)=floor(7 {e_i u/7}).
```

Away from the null wall set,

```text
s_i((u+j)/7) = z_i(u) + j r_i  (mod 7).
```

Therefore the exact signed kernel average satisfies

```text
49(-F_6(e))
 = (1/7) integral_0^1 sum_{j in F_7} L_6(z(u)+j r) du
 <= C(r)/7,                                                     (1)
```

where

```text
C(r)=max_{z in F_7^5} sum_{j in F_7} L_6(z+jr).                 (2)
```

Multiplying `r` by a nonzero scalar only permutes `j`, so `C` lives on the
`(7^5-1)/6=2801` nonzero projective directions.  Exact enumeration gives:

| `C(r)` | projective directions |
|---:|---:|
| 62 | 5 |
| 50 | 60 |
| 48 | 175 |
| 36 | 1075 |
| 34 | 580 |
| 32 | 65 |
| 26 | 60 |
| 25 | 270 |
| 24 | 280 |
| 23 | 120 |
| 22 | 30 |
| 21 | 60 |
| 20 | 10 |
| 18 | 10 |
| 8 | 1 |

Since `33/343<0.097`, every direction with `C(r)<=33` closes.  There are

```text
906 / 2801
```

such projective directions.  In particular, if the five moving speeds have
pairwise distinct nonzero residues modulo seven, all `120` projective directions
have `C(r)=25`, and (1) proves the diameter-free bound

```text
-F_6(e) <= 25/343 = 0.072886... < 0.097.                       (3)
```

The sieve is complementary to `THM-906/907`: it uses exact signed cancellation
before any Fourier absolute value, while their B3/B4 formulas retain the relation
lattice inside the `1895` directions whose raw affine-line ceiling is `34` or more.

## Honest boundary

- [x] prove the seven-shift affine-line identity;
- [x] enumerate all `2801` projective residue directions exactly;
- [x] close `906` directions and the distinct-nonzero residue subcase;
- [ ] bank the dependency-free verifier and tournament audit;
- [ ] use quotient-speed / relation data inside the `1895` bad directions;
- [ ] compose the finite-`t` wall remainder after the limiting sign closes.

## Assumption challenge

The vertices are projective residue directions, not runners or arcs.  This quotient
preserves the exact seven-shift signed-kernel ceiling, but destroys the within-cell path
`z(u)`, quotient speeds `floor(e_i/7)`, additive relations, and wall chronology.  The
failure set is therefore not a counterexample set: it is precisely the set where those
discarded sidecars must be restored.
