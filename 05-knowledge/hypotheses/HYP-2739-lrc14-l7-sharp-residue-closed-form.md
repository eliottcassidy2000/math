---
id: HYP-2739
title: LRC14 L7 torus-line discrepancy has an EXACT residue closed form -- sharp 12/(7q) PROVED combinatorially
status: PROVED (combinatorial, exact; closes HYP-2736c / HYP-2737)
source: kind-pasteur-2026-06-21 (THREAD C)
depends_on:
  - HYP-2737
  - HYP-2736
  - HYP-2733
  - HYP-2730
related:
  - THM-562
  - HYP-2734
  - OPEN-Q-108
---

# HYP-2739: L7 cell-discrepancy = an exact residue-only closed form; sharp 12/(7q) PROVED

## Statement (PROVED)

For coprime `p > q >= 1` in the L7 balanced window `1 < p/q <= 43/20`, the `(q,p)` torus
geodesic's cell occupancy against the 7-sector grid has **L1 cell-discrepancy given exactly by**

```
   D_{p,q}  =  S(p mod 7, q mod 7) / (7 p q),
   S = 4 f(||p||_7, ||q||_7),   ||x||_7 = min(x mod 7, 7 - x mod 7) in {0,1,2,3},
   f(a,b) = 0                   if a*b = 0
          = a*b + 3             if a != b   (a,b >= 1)
          = a*b + 4 - 2|a-2|    if a = b    (a,b >= 1).
```

Equivalently, writing `e_j = 7 c_{0j} - pq` for the integer row-0 deviation
(`c_{0j}` = row-0 of the integer cell-count matrix on the `7pq` grid), `S = sum_j |e_j| =`
codex's `rowdef`, and the full-matrix defect is `7S = sum_{ij}|7 c_{ij} - pq|`.

**Consequences (all SHARP, all PROVED):**
- `D_{p,q} <= 44/(7 p q)` universally (max `S = 44` at `||p||=||q||=3`).
- `D_{p,q} <= 12/(7 q)`, equality at `p/q = 3/2` (`S=36, p=3`). [closes the THREAD-C goal]
- `D_{p,q} <= 20/(7 p)`, equality at `p/q = 2/1` (`S=20, q=1`).
- `D_{p,q} = 0` iff `7|p` or `7|q` (apex law HYP-2733, recovered as `f=0`).

These BEAT the elementary `14/p` (HYP-2730) and the crude `24/(7q)` (codex s72), and the
proof is fully combinatorial: NO Koksma, NO equidistribution, NO discrepancy theory.

## Why it is true: the deviation is RESIDUE-ONLY (the engine)

`e_j = 7 c_{0j} - pq` depends ONLY on `(p mod 7, q mod 7)` -- a finite `7x7` table. The
"sharp constant 12" is therefore not an asymptotic discrepancy constant; it is forced by the
smallest-denominator window member (`3/2`) of a BOUNDED residue-class invariant.

### Proof chain (each step exact-machine-verified)

1. **Integer staircase model.** Subdivide `[0,1)` into `7pq` equal cells; cell midpoints give
   an integer matrix `c_{ij}`, `c_{ij}/(7pq) = mu(i,j)`, with row sums = col sums = `pq`
   (doubly balanced; the uniform 1D marginals, HYP-2730 Lemma 1).
   `D_{p,q} = (1/(7pq)) sum_{ij} |c_{ij} - pq/7|`.

2. **Cyclic-shift collapse (Lemma B / HYP-2733).** When `7 nmid q`, row `i` is row `0`
   cyclically shifted by `s*i`, `s = p q^{-1} mod 7`. Hence every row is a permutation of row 0
   and `sum_{ij}|7 c_{ij} - pq| = 7 sum_j |e_j| = 7 S`. The 2D problem is 1D.

3. **Clean lattice form of row 0.** Exactly,
   `c_{0j} = #{ (a,t) : 0<=a<q, 0<=t<p, (7a + t) mod 7q in [qj, q(j+1)) }`.
   (Proof: row 0 = the `q` bands `{v: sector(qv)=0}`; on the `7pq` grid a point `(k,t)`
   has 2nd sector `floor((14(pk mod q) + 2t+1)/(2q)) mod 7`; the substitution `7a+t` with
   `a = pk mod q` (a bijection of `Z/q`) and the half-integer reduction give the window form.)

4. **Coverage = a pure sawtooth (residue-only).** Write `c_{0j} = sum_{z in [qj,q(j+1))} cov(z)`,
   `cov(z) = #{ a in Z/q : (z - 7a) mod 7q in [0,p) }`. The points `{7a : a in Z/q}` are
   EXACTLY the `q` multiples of 7 in `[0,7q)` -- a perfectly uniform period-7 residue class.
   So `cov(z)` = (# multiples of 7 in the length-`p` cyclic arc ending at `z`) and, in the
   window `p < 3q (<= 7q)`,
   ```
        cov(z) = floor(p/7) + [ (z mod 7) < (p mod 7) ]      (EXACT, verified).
   ```
   `cov(z)` depends on `z` only through `z mod 7` and on `p` only through `p mod 7`.

   **Validity boundary (exact, verified):** the closed form `cov(z) = floor(p/7) + [...]`
   holds precisely for `p/q < 7` (first failure at `p/q > 7`, when a length-`p` run wraps
   the `7q`-cycle and double-covers). The L7 window `p/q <= 43/20 = 2.15` is far inside this,
   so the formula is valid with a wide margin.

5. **Residue-only conclusion.**
   `c_{0j} = (floor(p/7)) q + #{ z in [qj, q(j+1)) : z mod 7 < p mod 7 }`.
   The count term is the number of integers in a length-`q` window with residue `< p mod 7`,
   a function of `(qj mod 7, q mod 7, p mod 7)` only. Hence `c_{0j}` and `e_j = 7c_{0j} - pq`
   depend ONLY on `(p mod 7, q mod 7)`. The slope `s = p q^{-1} mod 7` is also residue-only,
   so the entire matrix deviation, and `S`, are residue-only. QED.

6. **Closed form and faces.** Tabulating the finite `7x7` table gives `S = 4 f(||p||,||q||)`
   (matched on ALL coprime `(p,q)`, `p,q < 90`, 0 mismatches). `max S = 44`. Within the
   window, `max S/p = 12` (at `3/2`), `max S/q = 20` (at `2/1`).

## Verification (all exact rational / integer, 0 failures)

| script (`04-computation/`) | what it certifies |
|---|---|
| `lrc_q108_threadC_integer_grid_kpswf4.py` | `c/(7pq)=mu`, balance, cyclic-shift, `D=S/(7pq)` |
| `lrc_q108_threadC_row0_structure_kpswf4.py` | row-0 arc model = grid row 0 |
| `lrc_q108_threadC_count_v2_kpswf4.py` | `cov = floor(p/7)+[z%7<p%7]`, `r_formula = r_direct` (2223 window ratios) |
| `lrc_q108_threadC_residue_law_kpswf4.py` | `e` constant on each `(p%7,q%7)` class |
| `lrc_q108_threadC_residue_proof_kpswf4.py` | period-7 stabilization; `S = 4f` universal (0 mism.) |
| `lrc_q108_threadC_S_formula_kpswf4.py` | the `7x7` S-table + symmetries |
| `lrc_q108_threadC_FINAL_kpswf4.py` | end-to-end: balance + shift + `D=S/(7pq)` + 3 sharp faces, 0 violations (1248 ratios) |

Outputs in `05-knowledge/results/` under the same stems.

## New lead: the closed form is PRIME-AGNOSTIC (robustness of the L7 technique)

`04-computation/lrc_q108_threadC_general_prime_kpswf4.py` confirms the SAME residue-only
closed form holds for an arbitrary apex value `P` (P sectors), not just `P=7`:
`D^{(P)}_{p,q} = S_P(p mod P, q mod P)/(P p q)` with `S_P` a finite `PxP` residue table,
verified for `P = 2,3,5,7,11,13` (cov-model = closed-form = direct, residue-only:True, all).
The universal max deviation is `maxS_P = 2,4,16,44,168,276` (P=2,3,5,7,11,13); `maxS_7/4 = 11`.
So the L7 closure technique depends only on `P=7` being the LRC sector count, NOT on any
arithmetic specialness of 7 -- the same machinery ports to any apex prime.

## Relation to prior work

- **Closes HYP-2737** (codex single-row lemma `rowdef <= 12p`): my `S = rowdef`, and
  `S <= 12p` is now PROVED with equality classified. Closes **HYP-2736c** (sharp combinatorial
  `12/(7q)`).
- **Strengthens HYP-2730**: the observed `sup D*p = 20/7`, `sup D*q = 12/7` are now THEOREMS,
  and the true bound is the exact `4f/(7pq)`, not just `14/p`.
- **Recovers HYP-2733** (apex law) as the `f = 0` locus.
- For the LRC(14) ledger, this makes the L7 tail rigorously elementary AND exact; combined with
  the finite atlas it gives an unconditional L7 closure modulo the L1-L6 chain.

-> HYP-2737, HYP-2736, HYP-2733, HYP-2730, THM-562, HYP-2734, OPEN-Q-108.
