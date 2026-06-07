---
title: LRC cell-leakage is the cyclotomic torsion of the base projection
session: monad-explorer-2026-06-06-S3
tags: [LRC, torsion, gcd, cyclotomic, composite-n, cell-model, S367, S377, HYP-1832, THM-427, n-grid-floor]
---

# LRC cell-leakage is the cyclotomic torsion of the base projection

## The seed, now a theorem

S377 (HYP-1832) noticed, numerically, that the hardest non-scalar defects in the composite-`n`
Lonely-Runner full-cell system are "torsion": the `n=14` minimum is the **coordinate-6 half-turn**
`(0,…,0,7,0,…,0)` (residue `7 = n/2`), the `n=15` minimum sits on the **order-3 subgroup `{5,10}`**.
The dispatched seed sharpened this to a slogan: *the leak lives in the p-torsion of `ℤ/n` that
projects to zero in the prime base, indexed by composite divisors — not a coincidence, it is forced.*

It is now forced, by an identity (**THM-427**). For the single-coordinate defect `r·e_i`:
```
        leak(r·e_i) = N_i·n − g·W_i(g),      g = gcd(r,n) = n / ord(r),
```
where `N_i` = #cells that only coordinate `i` can block (`i-exposed`), and `W_i(g)` counts exposed
cells whose `i`-bin is `≡ 0` or `≡ −1 (mod g)`. The proof is one line of counting: the `n−1` zero
coordinates block a candidate `(s,p)` independently of the shift `s`, so what survives them is
exactly the i-exposed cells; on those, the lone residue-`r` coordinate blocks `s·r ≡ −b, −b−1`, and
`s ↦ s·r` covers `gcd(r,n)·ℤ/n` exactly `g`-to-one. **The leak sees `r` only through `g = gcd(r,n)`.**

That single `gcd` is the entire torsion story:

- `g = gcd(r,n) = n/ord(r)`. The leak is a function of the **order of `r` in `ℤ/n`**.
- The nonzero **order-`p` torsion** `(n/p)ℤ/n \ {0} = {j·(n/p): 1≤j<p}` is precisely the set with
  `g = n/p`. Leak is therefore **constant** on it (THM-427 C1) — explaining why `5` and `10` tie at
  `n=15`, why `2,4,6,8,10,12` tie at `n=14`, etc.
- Every such residue is `≡ 0 (mod n/p)`: it **projects to zero in the `(n/p)`-runner base**. The
  base `ℤ/(n/p)` literally cannot see the defect — it is in the kernel of `ℤ/n ↠ ℤ/(n/p)`. That
  kernel *is* the order-`p` torsion. "Leaks through the base projection" = "lies in the projection's
  kernel" = "is `p`-torsion." The three phrasings are one fact.

## The order is the cyclotomic level — and lowest order is hardest

`leak = N_i n − g W_i(g)`. If the exposed bins were equidistributed mod `g`, `g·W_i(g) ≈ 2N_i`
would be flat and order would not matter. It is not flat: the staircase exposed-bins are **biased
toward `b ≡ 0,−1 (mod g)`, more so for larger `g`**, so larger `g` (smaller order) blocks more and
leaks less. Verified exhaustively (n ≤ 21, plus prime powers 8,9,16,25,27): the leak is weakly
decreasing in `g`, with the global minimum at the **largest proper divisor `g = n/p*`** (`p*` =
smallest prime) — the **lowest cyclotomic order `p*`** (HYP-2294). For even `n` this is the
**half-turn `n/2`, order 2** — the deepest possible torsion.

Prime powers make the mechanism naked: `n = p^a` gives a single **valuation chain**
`g ∈ {1,p,…,p^{a−1}}` and the leak drops monotonically down it
(`n=16`: `224,224,192,128` for `g=1,2,4,8`; `n=27`: `600,576,432` for `g=1,3,9`), minimum at order
`p`. Composite `n` with several primes is the CRT product of these chains; the leak classes are
indexed by the **full divisor lattice** `g | n` — exactly the seed's "indexed by composite divisors,"
with `ord(r) = n/g` the cyclotomic grading. A universal degeneracy falls out for free
(THM-427 C2): order-`n` (coprime) and order-`n/2` residues leak **identically and maximally**,
`N_i(n−2)` — the least-blocking defects are the highest-order ones.

## Where this sits in the LRC atlas (the unification)

This is the **n-grid / mod-`n` floor face** of the Lonely Runner, and it is governed by the
arithmetic of `ℤ/n` torsion — i.e. by **roots of unity**. Three threads converge:

1. **THM-369 (the n-grid floor).** The loneliness floor `M ≥ 1/n` is witnessed on the grid `k/n`
   whenever `n ∤ v`; the live modulus is `n`. The cell-leak lives on the same `n`-grid, and the
   defect's strength is its `ℤ/n` torsion order. Same modulus, same side.

2. **The signed `2n−1` shells are blind to exactly this (HYP-2286 / S709).** The signed-LRC
   refinement lives mod `2n−1`, and `gcd(n,2n−1)=1`, so `G_n ∩ G_{2n−1} = {0}`: the shell apparatus
   never touches the `n`-grid floor (`C′`). THM-427 is the **complementary face** the shells cannot
   see — it is *entirely* a mod-`n` torsion statement. So the two big LRC refinements of this cluster
   partition cleanly: **mod-`n` torsion (this) ⟂ mod-`2n−1` shells (signed)**.

3. **The cyclotomic floor (THM-403, S699o).** Across the program, "torsion = roots of unity =
   the cyclotomic FLOOR," and the hardness of LRC(14) is "torsion ARITHMETIC, prime-3, rank 0."
   THM-427 makes the floor literal at the cell level: the hardest obstruction is the **lowest-order
   root of unity** (the half-turn for even `n`; the order-`p*` element in general), and `n=14`'s
   notorious difficulty is the `3³` tower (a deep prime-3 valuation chain, the `n/2`-style chain of
   §prime-powers above writ large).

And it singles out the same element two other results do independently:
- **The half-turn `n/2` is the deletion guard (S702):** for even `n` the unique maximal AP-deletion
  is `a = n/2`, `M(AP_n\{n/2}) = 2/n`. Here `n/2` is also the unique minimum-leak defect. The
  order-2 torsion is simultaneously *the hardest to block* and *the most damaging to delete* — two
  faces of being the deepest cyclotomic element.
- This is the discrete shadow of complement/negation symmetry (`x ↦ −x` fixes exactly the order-≤2
  torsion), tying to the signed-LRC gauge `T1` (sign-blindness) being **exact** on `M`.

## What is open

- **Prove HYP-2294** (the bias lemma): show `g·W_i(g) − 2N_i ≥ 0` and increasing along divisor
  chains, i.e. the exposed staircase bins concentrate at `0,−1 (mod g)`. This is the only gap to a
  full "smallest-prime-wins" theorem.
- **The extremal coordinate `i*`** (S377 "product-sum resonance"): closed form for
  `argmax_i g·W_i(g)` at `g=n/p*`. Observed `i*` = 2,4,6,6,6,8,18,20 for n=6..21.
- **Beyond support-1.** THM-427 is exact only for one nonzero coordinate. The real HYP-1823 target
  (every nonzero gauge-normalized vector leaks) needs the multi-support union; the support-1 torsion
  defects are the *extreme rays*. Does the leak of a general defect decompose over the torsion content
  of its coordinates?

## Files
- `01-canon/theorems/THM-427-lrc-leak-is-a-function-of-gcd-torsion.md`
- `05-knowledge/hypotheses/HYP-2294-lrc-leak-minimized-at-smallest-prime-torsion.md`
- `04-computation/lrc_torsion_leakage_census_monad_s3.py`, `lrc_torsion_leakage_proof_monad_s3.py`
  (+ `05-knowledge/results/*.out`, incl. `..._primepower_monad_s3.out`)
- Builds on S367 (cell model), S377/HYP-1832, THM-369, HYP-2286/S709, THM-403, S702.
