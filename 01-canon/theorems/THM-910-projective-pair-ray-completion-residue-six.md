---
id: THM-910
renumber_note: was THM-909 at first claim; renumbered after pulling mac-mini's independently claimed BV-Fourier theorem
title: THE PROJECTIVE PAIR-RAY COMPLETION OF RESIDUE SIX — after the seven-shift signed quotient, exact zero-coordinate marginals and affine-invariant labeled pair-ray certificates close every one of the 2,801 projective speed-residue directions, proving the universal limiting bound -F6 <= 32/343 < 0.097
status: PROVED and solver-free VERIFIED — all 7^5 states, 2,801 projective directions, 1,505 zero-direction marginal rows, three sparse unordered-pair certificates, twenty nonzero invariant-pair certificates, 36 ordered ray endpoints, and 3,540 ordered-pair referees pass exactly
source: codex-2026-07-16-S18
depends_on: [THM-891, THM-908]
related: [THM-903-reflection-frame-residue6, THM-904, THM-905, THM-906, THM-907, THM-909, HYP-7081]
verification: 04-computation/lrc14_residue6_projective_pair_completion_codex_S18.py -> 05-knowledge/results/lrc14_residue6_projective_pair_completion_codex_S18.out
---

# THM-910 — the projective pair-ray completion

This theorem closes the **limiting negative residue-six sign** left by `THM-891`.
It does not yet compose the finite-`t` wall remainder and therefore does not by itself
prove the full lonely runner conjecture for fourteen runners.

Let

```text
S_r(z)=sum_{j in F_7} L_6(z+jr),
```

with `L_6=-K_6` as in `THM-908`.  That theorem proves

```text
49(-F_6) <= (1/7) E_u S_r(z(u)).                               (1)
```

Its raw affine-line ceiling already closes `906/2801` projective directions.  The
remaining directions close after restoring exactly the pair information that the raw
residue quotient discarded.

## 1. Directions with a zero coordinate

If `r_i=0`, primitivity gives `e_i=7q_i` and

```text
z_i(u)=sec(q_i u),
```

so `z_i` is exactly uniform on `F_7`.  Consequently

```text
U(r)=min_{i:r_i=0} sum_{a in F_7} max_{z:z_i=a} S_r(z)
```

satisfies

```text
-F_6 <= U(r)/2401.                                             (2)
```

There are `1505` projective directions containing a zero coordinate.  Exact enumeration
gives `U(r)<=204` for `1470` of them, hence

```text
-F_6 <= 204/2401 = 0.08496... < 0.097.                         (3)
```

The remaining `35` directions are precisely three projective permutation orbits:

| representative | directions | zero-coordinate pairs | pair-ray cap `t` | bound for `-F_6` |
|---|---:|---:|---:|---:|
| `(1,0,0,0,0)` | 5 | 6 | `9626/2527` | `57756/866761` |
| `(1,6,0,0,0)` | 10 | 3 | `97/14` | `291/4802` |
| `(1,2,4,0,0)` | 20 | 1 | `912/35` | `912/12005` |

For each row there is a literal rational weight `alpha` on the 28 unordered pair-sector
types such that

```text
S_r(z) <= sum_{zero pairs i<j} alpha({z_i,z_j})                 (4)
```

pointwise on all affine cosets.  The 22 exact pair vertices from `THM-891` give
`E alpha<=t`; the displayed bounds follow from (1).  Every zero-containing direction
therefore closes.

## 2. Directions with no zero coordinate

Among the `1296` projective directions in `(F_7^*)^5`, the raw ceiling closes `621`.
The remaining `675` directions form only `20` orbits under coordinate permutation and
common nonzero scaling.

For a representative `r`, define the affine-line invariant on coordinate pair `i<j`

```text
d_ij(z)=r_j z_i-r_i z_j mod 7.
```

The verifier stores, for each of the twenty representatives, ten seven-entry rational
functions `phi_ij`.  Exact enumeration of the `2401` affine cosets proves

```text
S_r(z) <= sum_{i<j} phi_ij(d_ij(z)).                            (5)
```

Because `d_ij` is unchanged by `z -> z+jr`, its expectation on `z(u)` equals its
expectation on the original pair of runner sectors.  After dividing that speed pair by
its gcd, its ordered residue pair is `lambda(r_i,r_j)` for some `lambda in F_7^*`.
The labeled form of `THM-891`'s ray law places the distribution on the segment from
uniform mass to the product-minimal endpoint for that ordered residue pair.  Checking
those seven vertices for every `phi_ij` proves a rational cap on the right side of (5).

All twenty exact caps clear the propagation threshold.  The largest is at orbit
representative `(1,1,1,1,4)`:

```text
sum t_ij = 3240009/140000,
-F_6 <= 3240009/48020000 = 0.06747... .                         (6)
```

Thus every no-zero direction also closes.

## 3. Universal limiting conclusion

The raw sieve, zero-marginal rows, three sparse pair certificates, and twenty
affine-invariant pair certificates cover all `2801` projective directions.  The largest
bound among the four branches is the raw ceiling `32/343`, so

```text
-F_6(E) <= 32/343 = 0.093294... < 0.097                        (7)
```

for every primitive five-speed slow core.  This proves the sole limiting sign left in
`THM-891`.

## Assumption challenge

Pair marginals failed in the original sector-count frame, and boxeph's diagonal ray mass
is exactly neutral.  The successful pair observable is different: it is the affine-line
determinant `r_j z_i-r_i z_j` **after signed seven-shift averaging**.  That quotient
preserves the kernel cancellation and the exact labeled pair rays while destroying
within-cell chronology and higher relation lattices.  The result therefore does not say
that raw runner pairs suffice; it says the correct signed residue quotient makes a
specific pair sidecar sufficient.

## Honest boundary

- [x] close all zero-coordinate directions;
- [x] reduce the nonzero collision locus to twenty symmetry orbits;
- [x] obtain rational pair-ray certificates on all twenty;
- [x] bank literal weights and an independent exact verifier;
- [ ] compose the limiting closure with the finite-`t` wall remainder.

Tournament Analysis on the twenty nonzero collision orbits uses rigorous obstruction
bounds as the pairwise observable.  Switching the gauge from raw affine ceiling to the
completed pair-ray cap flips `57` edges.  Both tournaments remain transitive with
singleton SCCs and unique tie Hamiltonian paths, so the scalar order is telemetry; the
determinant-class pair sidecar is the proof-bearing object.
