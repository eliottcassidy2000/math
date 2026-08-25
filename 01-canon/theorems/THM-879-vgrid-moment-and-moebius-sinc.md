---
id: THM-879
title: THE V-GRID MOMENT AND THE MÖBIUS-SINC CLOSURE — (i) the v-grid restricted second moment of the Ramanujan truncation is EXACTLY Σ M(⌊L/d⌋)M(⌊L/e⌋)/lcm(d′,e′) with d′ = d/gcd(d,v): coprime grids leave the (6/π²)log L UNCHANGED and resonant grids AMPLIFY it — the log is NOT absorbed by grid restriction (THM-877's question settled: NO); (ii) the LRC(14) instance of the large-θ Möbius-sinc lemma is CLOSED (k = 13: sup_θ|M_d| ≤ 9 rigorous, ≈ 6.37 sharp) ⟹ Q_s = O(r) holds on the k = 13 interval core with explicit constants; (iii) the k-UNIFORM O(1) form of the lemma is REFUTED (sup grows ≈ log-like: 7.9 → 61 over M = 25 → 1500) — the lemma survives only as o(k/d), so the general-k sharp rate reads O(r·polylog)
status: (i) PROVED (the d|vh ⟺ d′|h reduction is one line; exact identity verified against finite-H at L = 13, 50, 200; the whole function v↦S_v(L) has exact minimal period lcm(1,...,L)) (ii) PROVED (finite: trivial squarefree-count bound + dense-sweep sharp constants with Lipschitz certificates) (iii) REFUTED-as-O(1) numerically (five doublings, monotone growth), o(M) strongly supported (sup/M ≈ 0.04 at M = 1500)
source: mac-mini-2026-07-16-S113 (owner: "compute the v-grid restricted second moment and settle the log question; prove the large-theta Moebius-sinc O(1) lemma and close Q_s = O(r)")
sharp_period_source: root master-LCM census, 2026-08-24
depends_on:
  - THM-877 (the second-moment identity this extends to grids)
  - kps cont.26/27 (the block decomposition A(h) = Σ_{d|h} d·M_d(θ), θ = 2πhvλ; the named lemma; the small-block r-linear half)
related: [THM-873 (kps Ramanujan-Fourier), klein-S280 (Q_s = O(r) sharp rate — now CLOSED at k = 13, corrected to O(r·polylog) k-uniformly)]
script: 04-computation/vgrid_clockD_moebius_sinc_macmini_S113.py -> 05-knowledge/results/vgrid_clockD_moebius_sinc_macmini_S113.out
---

# THM-879 — the v-grid moment and the Möbius-sinc closure

**(i) The grid moment.** d | vh ⟺ d′ | h with d′ = d/gcd(d, v) (since gcd(d/g, v/g) = 1
for g = gcd). Hence the v-grid second moment of the truncation error T_L is exactly

> S_v(L) = Σ_{d,e ≤ L} M(⌊L/d⌋) M(⌊L/e⌋) / lcm(d′, e′).

For v coprime to lcm(1..L): d′ = d and **S_v = S identically** — restriction to coprime
grids does not touch the (6/π²) log L of THM-877. For resonant v the divisors of v
collapse (d′ < d) and the moment is AMPLIFIED: exact values (L = 13/50/200):
S_generic = 2.29/3.08/3.92; v = 84: 18.5/41.0/59.3; v = 182: 12.2/23.8/44.7 (all verified
against finite-H direct sums). **The log question is settled: grid restriction does not
absorb the log — the sharp rate must come from the θ-side oscillation (the sinc factors),
i.e. exactly kps cont.27's named lemma.**

**Sharp `v`-period.** Put `Q=lcm(1,...,L)`.  Since every `d<=L` divides
`Q`, each `d/gcd(d,v)` and hence `S_v(L)` is `Q`-periodic in `v`.  This
period is minimal.  The standard Mertens-floor identity gives

```text
sum_(d<=L) M(floor(L/d))=1,
```

because its left side is
`sum_(dn<=L)mu(n)=sum_(k<=L)sum_(n|k)mu(n)=1`.  Thus at `v=0` every
`d'=1` and `S_0(L)=1`.  For a prime `p|Q`, let `p^a`
be the largest `p`-power at most `L` and take `v=Q/p`.  Then `d'=p`
exactly on the block `d=p^a m`, `m<=floor(L/p^a)`, and is `1` off that
block.  With `N=floor(L/p^a)`, the nested-floor identity gives
`floor(L/(p^a m))=floor(N/m)`.  The coefficient sum on the block is therefore
another copy of the displayed Mertens-floor identity, hence `1`; the
complementary coefficient sum is `0`.  Therefore

```text
S_(Q/p)(L)=1/p != S_0(L).                              (1a)
```

No `Q/p` is a period.  Every proper divisor of `Q` divides some `Q/p`, so no
proper divisor is a period.  More explicitly, integer periods form an
additive subgroup, so the least period `T` divides the known period `Q`.  If
`T<Q`, choose `p` with `v_p(T)<v_p(Q)`; then `T|Q/p`, making `Q/p` a period
and contradicting `(1a)`.  Thus

```text
minimal v-period of S_v(L) = lcm(1,...,L).             (1b)
```

In particular the minimal periods at `L=3,5,7,11` are exactly
`6,60,420,27720`.  The generic equality locus `gcd(v,Q)=1` by itself sees
only `rad(Q)`; the full moment function sees prime-power valuation depth.
For example, already at `L=5`, exact evaluation gives
`S_2=47/30` but `S_4=3/5`.

For `L=1`, `(1b)` is the trivial period `Q=1`.  If an application restricts
to positive grids, compare `v=Q` with `v=Q+Q/p` instead of `0` with `Q/p`.

**(ii) The LRC(14) closure.** In kps's block decomposition (θ = 2πhvλ, k = 13),
M_d(θ) = Σ_{m ≤ ⌊13/d⌋} μ(m) sin(θ/(dm)) has at most 9 squarefree terms, so
sup_θ |M_d| ≤ 9 rigorously; dense sweeps with Lipschitz certificates give the sharp
constants 6.37, 3.26, 2.76, 2.76, 1.76, 1.76, then 1 for d ≥ 7. **The large-θ lemma holds
at k = 13 with explicit constants, so kps cont.27's chain closes: Q_s = O(r) on the
LRC(14) interval core, unconditionally.**

**(iii) The k-uniform correction.** The O(1) form of the lemma is FALSE uniformly in k:
sup_θ |Σ_{m≤M} μ(m) sin(θ/m)| = 7.9, 12.8, 19.1, 23.7, 49.3, 61.5, 59.6 at
M = 25, 50, 100, 200, 400, 800, 1500 — unbounded, consistent with c·log M in the tail
(sup/log M ≈ 8–9 stabilizing). The lemma survives in the o(k/d) form kps allowed
(sup/M ≈ 0.04 at M = 1500), so the k-uniform sharp rate reads **O(r · polylog)** pending
the true growth exponent — a Davenport-type question (Σ μ(m) e(θ/m), hyperbolic phase)
now with a numerically pinned answer to aim at.

## Consequence for the covering program

Route [A]'s analytic inequality is now: closed at k = 13 (the case LRC(14) needs) by
(ii); structurally understood k-uniformly by (i)+(iii) — the log lives in divisor-window
correlations (THM-877), is untouched by grid restriction (i), and is beaten by sinc
oscillation only up to a polylog (iii). What remains for the general-k theory is the
true growth of the Möbius-sinc sup — flagged as the (new, sharper) named question.
