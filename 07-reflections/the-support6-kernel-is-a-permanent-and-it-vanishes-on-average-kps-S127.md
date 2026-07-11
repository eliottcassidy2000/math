# The support-6 kernel is a permanent, and it vanishes on average

*kind-pasteur-2026-07-11-S127. Owner: "consider topics related to the support-6 relation-lattice Minkowski
count, search heavily, and prioritize closing mathematics before formalizing." I searched the frontier,
found the one thing in it that was genuinely closable, and closed it. This note is that closure and why it
is the right handle on the open lemma.*

---

## The frontier, and the closable thing inside it

The support-6 relation-lattice count is the project's **single named open lemma** for the LRC(14)
wide-spread half (MISTAKE-078, HYP-2613/2614/2644). Its shape: `meas(S7(E)) = M7(k) + Σ_{n∈Λ°(E)} K(n)`,
`K(n) = D7(n mod 7)/∏ n_j` (HYP-2646), and one must bound `corr(E) = Σ_c D7(c)·S_c(E)` uniformly over wide
`E`, where `S_c(E) = Σ_{n≡c, supp≥6} 1/∏n_j` are conditionally-convergent reciprocal sums. The naive absolute
Minkowski/successive-minima count is *required for convergence* but provably too lossy (5×–32×). So the whole
lemma is genuinely open and hard — the reciprocal sums `S_c` are the wall.

But the *weights* `D7(c)` are a finite object — 46656 numbers over `(F_7^*)^6` — and HYP-2613 had already
called them "a sector root-of-unity permanent" without proof, HYP-2646 had measured `|Re D7| ≤ 0.1431`
without a closed form. That is the closable piece: not the wall, but the exact algebraic weight the wall's
eventual summation-by-parts will run on. Prioritizing "close mathematics," I closed it.

## Three facts, all proved

**(1) It literally is a permanent.** Writing `h_T(c_j) = −A(c_j)·Σ_{i∈T} ζ^{−c_j i}`, the 64-term subset sum
collapses by inclusion–exclusion:
`Σ_{T⊆[6]} (−1)^{|T|} ∏_j (Σ_{i∈T} g_{ij}) = perm(g)`, because expanding the product over maps `φ:[6]→[6]`
with `im φ ⊆ T`, the alternating `Σ_{T⊇im φ}(−1)^{|T|}` kills every non-surjective `φ` and leaves exactly the
`6!` bijections. Hence `D7(c) = ∏_j A(c_j)·perm(ζ^{−c_j·i})_{i,j∈[6]}` — a 6×6 root-of-unity permanent. (The
same argument gives a *surjection* sum for supports `s > 6`.)

**(2) The sharp bound is closed-form, and its extremum has a mechanism.** `|D7(c)| = ∏|A(c_j)|·|perm| ≤
(max_r sin(πr/7)/π)^6 · 6! = 720·(sin(3π/7)/π)^6 = 0.643`, with equality **iff `c` is a constant coset**
`c≡3` or `c≡4`. The mechanism is a coincidence of extremes: the constant coset simultaneously maximises every
`|A(c_j)|` *and* makes the permanent matrix **rank-1** (all columns equal), so `perm = 6!·ζ^{−21c} = 6!`
(`21≡0 mod 7`). The maximally-coherent relation — all six coordinates in one residue class — is the peak.
This turns HYP-2646's measured `0.1431` into a finite-check theorem with a reason.

**(3) It vanishes on average — and *that* is the useful part.** Because the permanent is column-multilinear
and column `j` depends only on `c_j`, `Σ_c D7(c) = 6!·∏_{m=1}^6 B(m)`, `B(m) = Σ_{r=1}^6 A(r)ζ^{−rm}`. A
one-line character telescoping gives `B(m) = 0` for `m = 1,…,5` (only `B(6) = −7/(2πi)` survives), so
**`Σ_{c} D7(c) = 0`.** The kernel has zero mean over the coset space.

## Why the zero-mean fact is the right handle

`Σ_c D7(c) = 0` means the correction is *centered*: `corr(E) = Σ_c D7(c)·S_c(E) = Σ_c D7(c)·(S_c(E) − S̄)`
for **any** constant `S̄`. So the wide-spread correction never sees the *average* size of the reciprocal
sums — only their **variation across cosets**. This is exactly why the signed series converges where the
absolute envelope diverges (MISTAKE-078): the divergent, coset-uniform part of `S_c` is annihilated by the
kernel before it can contribute. The open lemma's summation-by-parts (HYP-2614, THM-504-D) is precisely a
tool for turning "the mean cancels" into "the tail is `O(1/w)`"; this corollary is the algebraic identity that
licenses it. It also explains the guardrail HYP-2614 flagged — the cancellation is *global* (`Σ_c`, over the
whole coset space), not *per-coordinate* (`Σ_{c_1}`, which the same computation shows does **not** vanish,
`B`-telescoping only closes when all six columns are summed).

## The honest scope

This closes the **weight** structure of the support-6 kernel — a complete, verifiable piece (THM-699). It does
**not** close the open lemma: the reciprocal sums `S_c(E)` and their conditional convergence over wide `E`
remain the wall, and that is a multi-session Diophantine problem (Erdős–Turán / cotangent–Dedekind summation
on the relation hyperplane after finite low-height wall deletion). But the frontier is now better armed: the
weight is a bounded, zero-mean, closed-form permanent, and the correction is provably a centered sum. When
the `S_c` bound is written, these are the exact facts it will multiply against.

The lesson I keep relearning: at a hard open frontier, *look for the finite object inside it*. The lemma is
infinite and open; the weight `D7` is 46656 numbers with a permanent structure, and that structure was
sitting there, named but unproved. Closing it is real progress even when the wall still stands — it is the
part of the wall you can carve into a keystone.

*Files: `04-computation/lrc14_support6_permanent_kernel_kps_S127.py` (+`.out`), canon `THM-699`. Search
synthesis over THM-538/532/537/533, THM-501/503/504/515/518, HYP-2613/2614/2616/2644/2645/2646, MISTAKE-078,
THM-681/684. The open lemma is the signed reciprocal hyperplane tail bound; this arms its weight side.*
