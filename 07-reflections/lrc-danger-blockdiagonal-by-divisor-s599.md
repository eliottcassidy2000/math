---
source: claudebox-2026-06-03-S599
status: REFLECTION — the danger structure block-diagonalizes by divisor (= Dirichlet level = orbit
  type); H3–H6 are one identity; the n=2q defect is a rank-1 two-block.
tags: [rigidity, divisor-block, Dirichlet-character, Ramanujan, gcd, isostatic, Maxwell, CRT,
  v2, 2-adic, apex, LRC]
---

# One identity blocks the danger: gcd(j,n) = Σ φ(d)[d|j]

**Prompt (human):** four refinements — H3 spectral/character (block-diagonalize by Dirichlet
characters; rigidity leaks through gcd>1), H4 isostatic (worry-set critically rigid, over-rigid
forbidden, #coincidences = φ(n)), H5 CRT (n=2q defect in the 2-block, q-block = prime lifted), H6
rigidity-height = v₂(n).

All four are the same sentence in four dialects, and the sentence is a hundred-year-old identity.
The danger depth of the arithmetic progression on the n-clock is `d(j) = gcd(j,n) − 1`, and

```
gcd(j,n) = Σ_{d | n} φ(d) · [d | j].
```

Read it as block-diagonalization: the danger is a direct sum over divisors `d|n`, and a divisor is
exactly a Dirichlet-character *level* (the characters of period `n/d`) and exactly an *orbit type*
of the `(Z/n)*` action (the gcd-`d` clock points). So the user's "symmetry-adapted orbit rigidity
matrix" is literally this decomposition; the spectral, the arithmetic, and the orbit pictures are
one object seen in three bases.

Once you have the identity the four hypotheses fall out as readings of it:

- **H3.** The danger lives on the `gcd>1` clocks (`d(j)=0` iff `gcd(j,n)=1`), and in the character
  spectrum its nontrivial mass sits in the imprimitive (`gcd(a,n)>1`) blocks — at `n=14`, the
  2-block and the 7-block. (The one honest caveat: the top block `d=n` is the trivial `j=0`
  collapse and it does spray onto the primitive characters; the *proper*-divisor defect is what is
  imprimitive. Stated precisely, H3 is exact.)
- **H5.** For `n=2q` the divisors are `{1,2,q,2q}`, so the proper-divisor defect is the 2-block plus
  the q-block. The q-block is a single clock — the apex `j=q` — of rank `φ(q)=q−1`: the solved
  odd-prime case lifted. The 2-block is the even-`j` defect, of rank `φ(2)=1`. **The whole `n=2q`
  obstruction is one-dimensional.** That is why a *pair* sieve closes it: there is a single rank-1
  block to clear, sitting over the `±`-orbit, and a pair-arity tool is exactly its size.
- **H4.** The unit block (`gcd=1`) carries the witnesses: the AP is lonely at exactly the `φ(n)` unit
  clocks, two binders apiece. That is the isostatic count — one time DOF pinned by a `±`-pair, at
  `φ(n)` symmetry copies. A counterexample would be over-rigid (no witness), but the unit orbit
  always supplies `φ(n) ≥ 1`. Maxwell's bound from below, in orbit language: you cannot have fewer
  coincidences than the orbit floor.
- **H6.** The 2-adic part is a *tower*. The doubling map is nilpotent on `ℤ/2^a` with index `a`, so
  the rigidity-height is `v₂(n)`: `n=14` is height 1, one collapse above the prime; `n=4q` is height
  2; powers of two are the deep water. Each level is one more 2-block to clear.

## What the assembly buys

The architecture writes itself. LRC at `n` is the problem of clearing the proper-divisor blocks of
`d(j)`. The odd-divisor blocks are lifted prime cases — solved by the polynomial method on each
prime. The obstruction is the 2-adic tower of height `v₂(n)`, whose bottom rung is a rank-1 block.
The sieve needed has depth `v₂(n)` and, on each block, arity equal to the orbit it sits over. The
`n=2q` frontier is the minimal nontrivial instance: height 1, rank 1, cleared by the pair-sum
sieve — which is exactly what the repo found empirically (HYP-2075, "multi-sieving has no apex,"
complete at n=14..22) without yet seeing that it was clearing a single rank-1 two-block.

## The transcending pattern

The deepest theorems often turn out to be a known identity read in the right basis. Here the entire
"rigidity zoo" — spectral, isostatic, CRT, height — is `gcd(j,n) = Σ φ(d)[d|j]`, the statement that
the clock decomposes into its divisor towers. The Lonely Runner is hard exactly on the 2-adic tower,
trivial on the unit block, and lifted-prime on the odd towers; the difficulty index is the height of
the 2-tower and the cyclotomic profile of the odd ones. To prove the conjecture is to clear a tower
of blocks whose total height is `v₂(n)` and whose hardest rung is, at the frontier, one-dimensional.

**Artifacts:** `04-computation/lrc_rigidity_spectral_crt_isostatic_s599.py` (+`.out`); new
**HYP-2145**. Develops the human's H3–H6 on HYP-2140 / HYP-2124, with HYP-2063, HYP-2075, HYP-2135.
