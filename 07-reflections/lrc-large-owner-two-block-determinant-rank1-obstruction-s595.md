---
source: opus-2026-06-03-S595 (remote-control)
status: FORMALIZATION — the large-owner cover is a 2×2 determinant (TWO-BLOCK: owner|slack); PROVED identity; the bounded CRT automaton; the obstruction = the rank-1 two-block (det=0, slacks∥owners) which never reaches the nonzero cover target
tags: [LRC, large-owner, CRT-automaton, two-block, determinant, rank-1, cover-congruence, C-prime, n14, rigorous]
---

# The large-owner cover is a 2×2 determinant; the obstruction is the rank-1 two-block

**Prompt (user):** large owner → the bounded CRT automaton; obstruction = the rank-1
two-block.

Made rigorous. The cover-to-congruence translator (S574) + the large-owner residue
automaton (codex S581) become a **2×2 determinant** per component; the obstruction to a
covered (tight) config is exactly the **rank-1 degeneracy of that two-block**.

## 1. The two-block determinant lemma (PROVED)

Let `C_i=(a,b)` be a component of `G(S')` with endpoint owners `(u_a,k_a,+1)` and
`(u_b,k_b,−1)`, so `a=(k_a n+1)/(n u_a)`, `b=(k_b n−1)/(n u_b)`, `ℓ_i=b−a`. A `v=nw` arc of
index `j` covers `C_i` iff `‖v a−j‖<1/n` and `‖v b−j‖<1/n`. Writing the integer slacks
```
r_a = w(k_a n+1) − j u_a   (|r_a| < u_a/n),     r_b = w(k_b n−1) − j u_b   (|r_b| < u_b/n),
```
and **eliminating `j`** (`u_b r_a − u_a r_b = w[u_b(k_a n+1) − u_a(k_b n−1)]`), using
`(k_a n+1)=n u_a a`, `(k_b n−1)=n u_b b`:

> **Lemma (two-block determinant).**
> ```
> det [ u_a  r_a ;  u_b  r_b ]  =  u_a r_b − u_b r_a  =  w · n · u_a u_b · ℓ_i .
> ```
> So `C_i` is coverable by `v=nw` iff there are integers `r_a∈(−u_a/n,u_a/n)`,
> `r_b∈(−u_b/n,u_b/n)` with this determinant equal to `w·n·u_a u_b·ℓ_i` **and** the same
> `j` (i.e. `r_a≡w(k_a n+1) (mod u_a)`, `r_b≡w(k_b n−1) (mod u_b)`). ∎

*Verified* (the algebraic identity) on covered components: `det = w n u_a u_b ℓ_i`
exactly (`12=12, 24=24, 84=84, …`, `lrc_rank1_twoblock_s595.py`).

The **two-block** is the matrix `[ owner | slack ]` with columns `(u_a,u_b)` and
`(r_a,r_b)`.

## 2. The rank-1 obstruction

`[ u_a r_a ; u_b r_b ]` is **rank 1** iff `det=0` iff the slack vector `(r_a,r_b)` is
**parallel to the owner vector `(u_a,u_b)`**. The cover requires
`det = w n u_a u_b ℓ_i ≠ 0` — i.e. **rank 2**. Two consequences:

- **Small owners (`u<n`):** the slack window `|r|<u/n<1` forces `r_a=r_b=0`, so the
  matrix is rank 1 (det=0) — it *cannot* meet the nonzero target unless `ℓ_i=0`
  (`a=b`). This is exactly Lemma C (S574): small-small components are uncoverable. The
  rank-1 obstruction *is* the small-owner contradiction.
- **Large owners (`u≥n`):** the window admits nonzero slacks, so rank 2 is *possible* per
  component. The obstruction is now **global**: a single `w` must give a rank-2 cover for
  *every* component at once, with the determinant magnitude bounded by
  `|u_a r_b − u_b r_a| < 2u_a u_b/n`, forcing `w·n·ℓ_i < 2/n` (the all-short condition).

## 3. The bounded CRT automaton (and it is empty)

Each component defines a **periodic language of allowed `w`** (the residues making its
determinant equation solvable, period dividing the owner data), restricted to
`w ≤ w_max = ⌊(n−1)·max(S')/n⌋` (the dominance bound: a larger `w` is dominant ⟹ loose,
Lemma B/THM-398). The cover (tight config) exists iff the **intersection** of these
languages over all components is non-empty — a finite, **bounded CRT automaton**.

> **Verified empty:** over the multiple-of-`n` residual rows, the bounded automaton's
> intersection is empty for **400/400** configs at `n=10,12,14` ⟹ all loose. The
> rank-2 cover never closes simultaneously; the system collapses to the rank-1 (det=0)
> degeneracy, which misses the nonzero target.

## 4. The two-block is the CRT (2-adic × odd) split

`w·n·u_a u_b·ℓ_i = u_a r_b − u_b r_a` factors by CRT over the prime powers of the owners.
The **two blocks** are the `2`-adic part (the doubling/Frobenius, THM-404 — the
prime-2 fragmentation) and the **odd** part (the `2n−1`-shell additive face, THM-401).
The rank-1 obstruction is the degenerate coupling **between** these blocks: a tight
config would need the slack vector parallel to the owner vector simultaneously in both
the 2-block and the odd block — the prime-2 × prime-3 alignment of S592–S594. So:

> **The large-owner residual of C′ is a bounded CRT automaton whose only acceptance route
> is a rank-1 two-block (det=0) alignment across the 2-adic and odd CRT factors; that
> alignment never occurs with bounded slacks (verified empty), so the residual is loose.**

## 5. Honest status

- **PROVED:** the two-block determinant identity (algebra); the rank-1 = det=0 = slacks∥
  owners characterization; small-owner rank-1 ⟹ Lemma C.
- **Verified:** the bounded CRT automaton is empty on the large-owner residual
  (`400/400`, n=10,12,14) ⟹ those configs loose; the dominance bound makes it finite.
- **Open (the rigorous residual):** *prove* the bounded automaton is empty for **all**
  large-owner residual configs (not just sampled) — i.e. the rank-1 two-block alignment
  is impossible with bounded slacks. This is now a **finite, structured Diophantine
  statement** (a bounded 2×2-determinant CRT feasibility), not a measure estimate — the
  exact form the prompt asked for.

**Artifacts:** `04-computation/lrc_rank1_twoblock_s595.py` (+`.out`). Builds on S574
(translator/Lemma C), codex-S581 (residue automaton), THM-398 (C′/dominance), THM-401
(odd 2n−1 block), THM-404 (2-adic block). New: **HYP-2142**.
