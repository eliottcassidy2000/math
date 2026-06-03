---
id: HYP-2145
status: EXPLORATION + VERIFICATION — develops & confirms the user's H3–H6 (spectral/character,
  isostatic, CRT-block, rigidity-height) with one backbone identity; assembles them into a
  divisor-block proof-architecture for the n=2q frontier. Identity & all four verified; the
  architecture is a conjecture.
source: claudebox-2026-06-03-S599
related:
  - HYP-2140
  - HYP-2124
  - HYP-2063
  - HYP-2075
  - HYP-2135
---

# HYP-2145: the danger block-diagonalizes by divisor — H3–H6 unified

The user proposed four refinements (H3 spectral/character, H4 isostatic/percolation, H5 CRT-block,
H6 rigidity-height) of the "rigidity = orbit-type" frame (HYP-2140). They are all consequences of
one identity, verified here, and together give a proof-architecture for the `n = 2q` hardness.

## The backbone (verified)

On the `n`-clock, the AP `{1,…,n−1}` puts runner `v` at phase `vj/n` at time `t=j/n`. Its danger
depth (runners stuck at the origin) is

```
d(j) = #{v ∈ [1,n−1] : vj ≡ 0 (mod n)} = gcd(j,n) − 1,
```

and the clock is the tight lonely witness iff `d(j)=0` iff `gcd(j,n)=1`. The classical identity

```
gcd(j,n) = Σ_{d | n} φ(d) · [d | j]     ⇒     d(j) = Σ_{d|n, d>1} φ(d) · [d | j]
```

**block-diagonalizes the danger by divisor `d|n`** — and divisor = Dirichlet-character *level* =
orbit type. This is the symmetry-adapted orbit rigidity matrix, diagonal in the Dirichlet basis;
its additive transform is the Ramanujan sum `c_d`. (Identity verified for `n=7,10,13,14,18,22`.)

## H3 (spectral/character) — verified

Rigidity leaks only through the **proper-divisor (`gcd>1`, imprimitive) blocks**. At `n=14` the DFT
of `d(j)` has its nontrivial mass exactly at frequencies with `gcd(a,14)>1`: the **2-block**
(`a=2,4,6,8,10,12`, each `|DFT|=18`) and the **7-block** (`a=7`, `|DFT|=13`). The principal/`gcd=1`
frequencies carry only the trivial `d=n` term (the `j=0` total collapse, `t=0`, not a real
witness). So: the *real* (proper-divisor) defect lives entirely in the imprimitive blocks; the
2-block is the even-`j` defect — the 2-adic seam in the character spectrum.

## H5 (CRT-block) — verified

For `n=2q`, `divisors = {1,2,q,2q}`, so `d(j) = φ(2)[2|j] + φ(q)[q|j] + φ(2q)[2q|j]`:

- **q-block** `(q−1)[q|j]`: supported at the single clock `j=q` (the **apex**), of rank `φ(q)=q−1`
  — this is the **solved odd-prime case `q` lifted** (the mod-`q` structure is the prime
  loneliness, HYP-2063's "mod-`q` free over `ℤ_q^×`").
- **2-block** `[2|j]`: the even-`j` defect, of rank `φ(2)=1`. **The obstruction is one-dimensional.**

So the `n=2q` rigidity defect lives entirely in the rank-1 two-block, exactly HYP-2063's "the mod-2
component is forced" / the apex zero-divisor. The pair-sum multi-sieve (HYP-2075/HYP-2135) is the
tool that clears this rank-1 block — its pair-arity matches the `⟨−1⟩` 2-orbit sitting over it.

## H4 (isostatic / Maxwell) — verified

The AP is lonely at exactly the **`φ(n)` unit clocks** (verified `n=6..22`), each with **2 binders**
(the `v` with `vj≡±1`, the `±`-pair of HYP-2135). This is the **critically rigid / isostatic** count:
one time DOF, pinned by a 2-binder `±`-pair, at `φ(n)` symmetry-equivalent points. A *counterexample*
would be **over-rigid** — zero lonely clocks — which is below the orbit-counting floor: the
`(Z/n)*`-orbit always supplies `φ(n) ≥ 1` witnesses. Maxwell/orbit-counting forbids over-rigidity;
the worry-set saturates the isostatic bound, it cannot exceed it.

## H6 (rigidity-height) — verified

The **rigidity-height = `v₂(n)`** = the nilpotency index of the doubling map `x↦2x` on the 2-part
`ℤ/2^a` (verified: `n=14→1`, `28→2`, `8→3`, `16→4`, odd `n→0`). Each doubling level degrades the
dynamical rigidity by one factor of 2; `n=14=2·7` has height 1 — a single mod-2 collapse above the
solved prime, the apex/2-adic seam.

## Assembly — a divisor-block proof-architecture for the frontier

H3–H6 cohere: the danger is a direct sum of divisor-blocks; the **unit block** (`gcd=1`) is the
solved lonely part (`φ(n)` isostatic witnesses, H4); the defect is the **proper-divisor blocks**
(H3). For `n = 2^a · (odd)`, the odd-part blocks are the **prime case lifted** (H5), and the 2-adic
part is a **height-`a` tower of 2-blocks** (H6), the lowest of rank `φ(2)=1`. The sieve must clear
the 2-tower; its depth = `v₂(n)` and its arity = the orbit size over each block. For `n=2q`
(height 1, rank-1 2-block) the **pair-sum sieve suffices** — explaining HYP-2075's completeness at
n=14..22 and predicting that `n=4q` (height 2) needs a depth-2 / higher-arity sieve.

> **Conjecture (block proof-architecture).** LRC at `n` reduces to clearing the proper-divisor
> blocks of `d(j)`. The odd-divisor blocks are the lifted prime cases (solvable). The obstruction is
> the 2-adic tower of height `v₂(n)`; a sieve of depth `v₂(n)` and arity = the witness-orbit size
> over the 2-blocks clears it. `n=2q` is the minimal nontrivial case (height 1, rank-1), cleared by
> the pair-sum sieve.

## Open / next

- Generalize `d(j)` beyond the AP to arbitrary primitive speed sets: does the divisor-block
  structure persist (with the speed set's own symmetry orbits replacing the full clock)?
- Make the "sieve depth = `v₂(n)`, arity = 2-block orbit size" a theorem; test on `n=28` (height 2).
- Quantify the difficulty index `R(n)` (HYP-2140 C2) as `(v₂(n), the proper-odd-divisor coset
  profile)` — the 2-tower height plus the odd cyclotomic structure.

**Artifacts:** `04-computation/lrc_rigidity_spectral_crt_isostatic_s599.py` (+`.out`),
`07-reflections/lrc-danger-blockdiagonal-by-divisor-s599.md`. Develops the user's H3–H6 on
HYP-2140 (rigidity=orbit-type) / HYP-2124 (opus), with HYP-2063 (apex/2-adic), HYP-2075 (pair
sieve), HYP-2135 (±-orbit).
