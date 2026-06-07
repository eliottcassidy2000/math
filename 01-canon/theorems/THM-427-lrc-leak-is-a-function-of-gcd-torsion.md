# THM-427 — The composite-LRC cell-leak of a single-coordinate defect is a function of gcd(r,n) (its torsion order); closed form leak = N_i·n − g·W_i(g)

**Status:** PROVED (exact, elementary). Closed form verified EXACT for **all** (coordinate i,
residue r) at n = 6, 10, 12, 14, 15, 18, 20, 21
(`lrc_torsion_leakage_proof_monad_s3.py`). Reproduces the S377 numbers
(n=14 → 56 at coord 6 / res 7; n=15 → 120 at coords 6,14 / res 5,10).
**Source:** monad-explorer-2026-06-06-S3. Makes the S377/HYP-1832 "torsion/CRT leak" picture
precise and proves its core.

## Setup (the S367 full-cell model, verbatim)

`n` runners, `k = n−1` mover speeds `v = (v_1,…,v_k) ∈ (ℤ/n)^k`. A **cell** is a distinct
floor-bin pattern `bins_i(α) = ⌊n·{i·α}⌋`, `i=1..k`, as `α` ranges over `(0,1)` (breakpoints
`a/(ni)`). A **candidate** is a pair `(s,p)` with shift `s ∈ ℤ/n` and cell `p`; there are
`N_cand = n·#cells`. Coordinate `i` with residue `r` **blocks** `(s,p)` iff
```
        (s·r + bins_i(p)) mod n  ∈  {0, n−1}
```
(the speed-`i` runner is within `1/n` of the origin → not a lonely time). The **leak** of `v` is
```
        leak(v) = #{ candidates blocked by NO coordinate }.
```
The scalar ramp `v_i = m·i` blocks every candidate (`leak = 0`); it is a gauge (adding `m·i`
reparametrizes), normalized to `v_1 = 0`. The target finite lemma (HYP-1823) is "every nonzero
gauge-normalized vector has a leak ≥ 1." THM-427 resolves the **single-coordinate** (support-1)
slice exactly.

## The theorem

Write `e_i` for the unit defect (coordinate `i` only). Call a cell `p` **i-exposed** if
`bins_j(p) ∉ {0, n−1}` for every `j ≠ i` (only coordinate `i` can avert the residue-0 baseline
block). Let
```
   N_i      = #{ i-exposed cells },
   W_i(g)   = Σ_{i-exposed p} ( [ g | bins_i(p) ] + [ g | bins_i(p)+1 ] ).
```

> **THM-427.** For every coordinate `i` and residue `r ≠ 0`, with `g = gcd(r,n)`:
> ```
>        leak(r·e_i) = N_i·n − g·W_i(g).
> ```
> In particular **leak(r·e_i) depends on r only through g = gcd(r,n)** — equivalently only through
> the **order** `ord(r) = n/g` of `r` in `ℤ/n`.

### Proof (elementary)

A candidate `(s,p)` survives the `k−1` zero coordinates `j≠i` iff none of them blocks it. Coordinate
`j` (residue 0) blocks `(s,p)` iff `(s·0 + bins_j(p)) mod n ∈ {0,n−1}`, i.e. iff
`bins_j(p) ∈ {0,n−1}` — independent of `s`. So `(s,p)` survives all `j≠i` **iff `p` is i-exposed**
(for every `s`). Restrict to i-exposed cells.

Fix an i-exposed cell `p`, `b = bins_i(p)`. The candidate `(s,p)` is a leak iff coordinate `i` also
fails to block, i.e. iff `s·r ∉ {−b, −b−1} (mod n)`. The map `s ↦ s·r (mod n)` has image
`gcd(r,n)·ℤ/n` and hits each value of that image exactly `g = gcd(r,n)` times; the congruence
`s·r ≡ c (mod n)` has `g` solutions if `g | c` and none otherwise. Hence the number of shifts `s`
that *do* block is
```
   g·[ g | b ]  +  g·[ g | (b+1) ]
```
(the two targets `c = −b` and `c = −(b+1)`). The number of leak shifts for this cell is therefore
`n − g([g|b]+[g|b+1])`. Summing over the `N_i` i-exposed cells,
```
   leak(r·e_i) = N_i·n − g·Σ_{exp p}([g|b_p]+[g|b_p+1]) = N_i·n − g·W_i(g).  ∎
```
The right-hand side involves `r` only through `g`. Since `g = gcd(r,n) = n/ord(r)`, it is equivalently
a function of `ord(r)`. ∎

## Corollaries (all PROVED)

> **C1 (torsion-constancy).** leak is constant on the nonzero elements of each cyclic subgroup of a
> fixed order. In particular it is constant on the nonzero **order-p torsion**
> `(n/p)ℤ/n \ {0} = { j·(n/p) : 1≤j<p }` (every such element has `g = n/p`). These residues are
> **exactly** the nonzero residues that project to 0 under `ℤ/n ↠ ℤ/(n/p)` — invisible in the
> `(n/p)`-runner base. This is the seed's statement, now a theorem (the n=14 half-turn `r=7=n/2`,
> the n=15 order-3 pair `r∈{5,10}=n/3·{1,2}`).

> **C2 (g=1 ≡ g=2 identity).** `leak(coprime r) = leak(g=2 r) = N_i·(n−2)` for every coordinate `i`.
> Proof: `W_i(1) = 2N_i`; and `W_i(2) = #{b even}+#{b odd} = N_i`, so `g·W = 2N_i` in both cases.
> So a coprime (full-order) residue and an order-`n/2` residue leak identically and maximally; the
> maximal-leak (least-blocking) defects are the high-order ones.

> **C3 (gcd partition).** The `n−1` nonzero residues split into `gcd`-classes; leak is one value per
> class. Over a coordinate, the leak values are indexed by the proper divisors `g | n`, `g<n`.

## Scope / what is NOT claimed here

- THM-427 is the **support-1** law. Multi-coordinate defects are governed by unions of masks and are
  not a function of coordinatewise gcd alone (HYP-1832's "support ≤ 4 scans"); the support-1 result
  is what pins the extremal single defects S377 found.
- The statement "the GLOBAL minimum leak sits at `g = n/p*` (smallest prime `p*`)" — i.e.
  `g·W_i(g)` is maximized at the largest proper divisor — is **verified n ≤ 21** but is a separate
  combinatorial fact about the bin distribution `W_i`; it is recorded as **HYP-2294**, not proved
  here.

## Why it matters

- Turns "the leak is locked into torsion" (the seed) from a numerical observation into the identity
  `leak = N_i·n − g·W_i(g)`, with the torsion appearing **exactly** as `g = gcd(r,n) = n/ord(r)`.
- Cleanly separates the two LRC moduli: this leak lives on the **n-grid / mod-n floor side**
  (THM-369), where the relevant invariant is `ord(r)` in `ℤ/n`; it is the cyclotomic-torsion face,
  **complementary** to the signed `2n−1` shell side (S709/HYP-2286 showed the signed apparatus is
  blind to the n-grid). Same "torsion = roots of unity = cyclotomic floor" theme as THM-403/S699o.

## Sources
- `04-computation/lrc_torsion_leakage_proof_monad_s3.py` (+ `.out`) — (A) gcd-class law, (B) closed
  form exact for all (i,r), n≤21.
- `04-computation/lrc_torsion_leakage_census_monad_s3.py` (+ `.out`) — full support-1 census.
- `07-reflections/lrc-leakage-is-cyclotomic-torsion-of-the-base-projection-s3.md`.
- Builds on S367 (cell model), S377/HYP-1832 (torsion/CRT leak observation), THM-369 (n-grid floor),
  HYP-2286 (signed/2n−1 blind to n-grid). Related: HYP-2294 (min-at-smallest-prime), HYP-1823.
