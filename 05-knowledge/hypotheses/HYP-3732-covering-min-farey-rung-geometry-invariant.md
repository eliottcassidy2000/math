---
id: HYP-3732
title: The covering-min's GEOMETRY INVARIANT is the STERN-BROCOT / FAREY RUNG k(n) -- M(n) = [0; n-1, k] = k/(k(n-1)+1), the unique Farey fraction ADJACENT to the ceiling 1/(n-1) (|det|=1) at Farey-HEIGHT k above the floor 1/n (|det with 1/n| = k-1 = the margin numerator). The construction n/Phi_6(n) is rung k=n; the floor 1/n is rung 1; the ceiling 1/(n-1) is rung inf. The rung k ~ D/n is the binding-modulus STRETCH per runner. New sequences: rung k(n)=2,2,4,4,3 (n=7..11, irregular), binding D(n)=k(n-1)+1, defect delta=(k-1)/k, construction-deficit n-k(n). HONEST CORRECTION: the user/earlier-search value n=13 = 1/12 is WRONG -- the construction {1..11,156} gives 13/157=[0;12,13]<1/12 (under-resourced: small speeds can't host the killer 156)
status: VERIFIED for n=7..11 (covering-min spreads 2/13,2/15,4/33,4/37,3/31 = rungs 2,2,4,4,3; each a Farey neighbor of 1/(n-1), det=1, re-checked by exact gap). HONEST: n=13=1/12 REFUTED (construction 13/157<1/12 verified, binding D=157). For n>=12 the best-known covering-min is the construction (rung n); whether a larger-speed spread gives a lower rung (the spread/construction transition) is being tested (n=12, Smax=80). The rung sequence has no known closed form (genuinely irregular).
source: klein-2026-06-30-S37
depends_on:
  - HYP-3731   # the covering-min IP (which pins the values)
  - HYP-3718   # the semiconvergent ladder q_j=q_{j-1}+(n-1) (the rungs)
related:
  - HYP-3724   # Phi_6 = the construction denominator (rung n)
  - HYP-3551   # 14/183 = the n=14 construction (rung 14)
  - HYP-3717   # three-gap / torus lift (the defect as rotation number)
results:
  - 04-computation/covering_min_rung_klein.py
  - 05-knowledge/results/covering_min_rung_klein.out
---

# HYP-3732 — the covering-min's Farey-rung geometry invariant

## The invariant: the Stern-Brocot rung k(n)
Writing the covering-min `M(n) = j/D`, its continued fraction is **`M(n) = [0; n-1, k]`**, so
`M(n) = k/(k(n-1)+1)` -- the rung-`k` fraction on the Farey ladder from the floor `1/n` (rung 1) to the
ceiling `1/(n-1)` (rung `inf`). Two clean Farey adjacencies (verified `n=7..11`):
- **`|det(M, 1/(n-1))| = |k(n-1) - (k(n-1)+1)| = 1`** -- `M(n)` is ALWAYS a Farey NEIGHBOR of the ceiling
  `1/(n-1)`.
- **`|det(M, 1/n)| = k-1`** -- `M(n)` is at Farey-HEIGHT `k` above the floor `1/n`; and `k-1` is exactly the
  **margin numerator** (`M - 1/n = (k-1)/(nD)`).

So the geometry invariant is the **Stern-Brocot DEPTH `k(n)`** (geometric: the Farey tessellation; arithmetic:
the continued fraction). Equivalently the **defect `delta(n) = n - 1/M(n) = 1 - 1/k = (k-1)/k`** (a
rotation-number-like fractional deficit, ties to three-gap, HYP-3717). And `k ~ D/n` is the **binding-modulus
stretch per runner**: tight (`k=2`, `D~2n`) for spreads, loose (`k=n`, `D=Phi_6(n)`) for the construction.

```
n   M(n)    rung k   delta=(k-1)/k   D=k(n-1)+1   construction n/Phi6 = rung n
7   2/13      2        1/2            13           7/43   (rung 7)
8   2/15      2        1/2            15           8/57   (rung 8)
9   4/33      4        3/4            33           9/73   (rung 9)
10  4/37      4        3/4            37           10/91  (rung 10)
11  3/31      3        2/3            31           11/111 (rung 11)
```

## New sequences (extensions)
- **rung `k(n)`** = `2,2,4,4,3` (`n=7..11`) -- the covering-min's Farey depth. Irregular, no closed form.
- **binding `D(n) = k(n-1)+1`** = `13,15,33,37,31`. `D = 2n-1` only at `n=7,8` (the Paley cases, HYP-3731).
- **defect `delta(n) = (k-1)/k`** = `1/2,1/2,3/4,3/4,2/3`.
- **construction-deficit `n - k(n)`** = `5,6,5,6,8` -- how far the spread pulls the covering-min below the
  construction rung `n`.
- **margin numerator `eps = k-1`** = `1,1,3,3,2` (`= |det(M,1/n)|`).

## HONEST correction: n=13 is NOT 1/12
The value `n=13 -> 1/12` (rung `inf`, the ceiling) is REFUTED: the construction `{1,...,11, 156}`
(`156 = 13.12`) is a primitive covering with `M = 13/157 = [0;12,13]` (rung 13) `< 1/12` -- verified exactly
(binding `D=157`). The `1/12` came from a small-`Smax` search that could not host the killer `156`. So the
covering-min at `n=13` is `<= 13/157` (rung `<= 13`), NOT `1/12`. (Same failure mode as my own under-resourced
S36 run; cf. MISTAKE-087 on the construction being non-extremal.)

## The irregularity
The rung `k(n) = 2,2,4,4,3` for `n=7..11` is genuinely irregular (no closed form): consecutive pairs `(7,8)`
and `(9,10)` share a rung (2 and 4), `n=11` drops to 3. For `n>=12` the best-known covering-min is the
construction (rung `n`), and the spread->construction TRANSITION near `n=11->12` is **weakly confirmed**: the
`Smax=48` IP returns the construction, the `Smax=80` IP timed out, and a randomized probe of 40000 primitive
coverings (speeds `<=60`, medium killers `18..60`) found NO spread beating `12/133`. So the rungs go
`2,2,4,4,3` (spreads, `n<=11`) then jump to `n` (construction, `n>=12`) -- a genuine regime change, not yet
proved. The irregularity is ARITHMETIC -- it depends on how the resonances `<= n` can be killed while spreading
(a prime `n` whose only small killer is `n` itself behaves differently from a composite).

## Connections (structural, not numerical)
The covering-min rung shares its character with the repo's other no-closed-form sequences -- `width(G_n) =
1,2,3,6,15,49` (breaks `C(n-2,floor)` at `n=7`), `chi(E_n) = 2,3,5,10,28`, `A000568`, `NS-merged =
0,1,2,22,184` -- but there is NO numerical coincidence between them. The genuine link is STRUCTURAL: all are
arithmetic invariants of the SAME staircase `delta_{n-2}` / `Phi_6` object (CLAUDE.md "everything is the
triangle"). The covering-min's construction rung is literally `n/Phi_6(n)`; the binding `D = 2n-1` (Paley) at
`n=7,8` is the tournament side (HYP-3731). They are different spectra of one circulant/Farey machine, irregular
for the same reason: they encode genuine number theory, not a closed form.

## Net
The covering-min's geometry invariant is the **Stern-Brocot / Farey rung `k(n)`**: `M(n)=[0;n-1,k]`, a Farey
neighbor of `1/(n-1)`, height `k` above `1/n`, `k ~ D/n` the binding stretch. Verified rungs `2,2,4,4,3`
(`n=7..11`); the construction is rung `n`; `n=13 = 1/12` is REFUTED (`13/157 < 1/12`). New sequences (rung,
binding, defect, deficit) extend it. The irregularity is arithmetic and shares only a structural (`Phi_6` /
circulant) anchor -- not a numerical one -- with the tournament sequences.
