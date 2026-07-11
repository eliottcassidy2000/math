---
id: THM-699
title: "The support-6 seven-sector kernel IS a 6×6 root-of-unity permanent — D7(c) = ∏_j A(c_j)·perm(ζ^{−c_j·i})_{i,j∈[6]} (inclusion–exclusion collapses the 64-subset T-sum to the 720-permutation permanent, only surjections survive); hence the SHARP closed-form bound |D7(c)| ≤ 720·(sin(3π/7)/π)^6 = 0.64308 for ALL cosets, with equality iff c is a constant coset c≡3 or c≡4 mod 7 (where the permanent matrix is rank-1, perm = 6!), proving/sharpening HYP-2646's empirical |Re D7| ≤ 0.1431"
status: PROVED (permanent identity by inclusion–exclusion, exact; verified max err 1.68e-19 over 300 random cosets, 1.68e-16 here; the sharp |D7| bound by |perm|≤6! + max_r|A(r)|; the |Re D7|=0.14310 ruler value by exhaustive finite check over all 6^6 = 46656 cosets). Resolves the "permanent" layer of HYP-2613; proves HYP-2646's ruler bound.
source: kind-pasteur-2026-07-11-S127 (cont.21) — searching the support-6 relation-lattice / Minkowski frontier for closable structure
depends_on:
  - HYP-2646   # K(n) = D7(n mod 7)/∏n_j, the coset/reciprocal factorization
  - THM-538    # the support-6 kernel setup (active-support ≥6)
related:
  - HYP-2613   # relative support-6 permanent count (this proves its "root-of-unity permanent" layer)
  - HYP-2614   # support-6 residue-cusp tail (this supplies the exact algebraic weight for the summation-by-parts)
  - HYP-2606   # the signed relation-lattice Fourier identity
external: inclusion–exclusion; the permanent of a matrix; Minkowski geometry-of-numbers (the count this weight feeds).
---

# THM-699 — the support-6 kernel is a 6×6 root-of-unity permanent

## Statement

For a support-6 offset relation `n` (six coordinates nonzero and ≢0 mod 7), write `c = n mod 7 ∈ (F_7^*)^6`.
The seven-sector signed Fourier kernel `D7(c) = Σ_{T⊆{1,…,6}} (−1)^{|T|} ∏_{j=1}^6 h_T(c_j)` (HYP-2646),
`h_T(r) = −A(r)·Σ_{i∈T} ζ^{−r i}`, `A(r) = (1−ζ^{−r})/(2πi)`, `ζ = e^{−2πi/7}`, equals

> **`D7(c) = ∏_{j=1}^6 A(c_j) · perm( ζ^{−c_j·i} )_{i,j∈[6]}`**  —  the 6×6 **root-of-unity permanent**.

Consequently:

> **(SHARP BOUND)** `|D7(c)| ≤ 720·(sin(3π/7)/π)^6 = 0.64308` for **every** coset `c`, with equality iff
> `c` is a constant coset `c ≡ (3,…,3)` or `c ≡ (4,…,4)` mod 7.
>
> **(RULER)** `max_c |Re D7(c)| = 0.14310`, attained at `c ≡ (4,…,4)` — proving HYP-2646's empirical
> `|Re D7| ≤ 0.1431` (`Re D7` is the ruler, since `K(n)+K(−n) = 2·Re D7(c)/∏n_j`).

## Proof

### The permanent identity
Factor `h_T(c_j) = (−A(c_j))·Σ_{i∈T} ζ^{−c_j i}`, so `∏_{j=1}^6 h_T(c_j) = (−1)^6 ∏_j A(c_j) · ∏_{j=1}^6
(Σ_{i∈T} g_{ij})`, `g_{ij} := ζ^{−c_j i}`. Then, with `(−1)^6 = 1`,
`D7(c) = ∏_j A(c_j) · Σ_{T⊆[6]} (−1)^{|T|} ∏_{j=1}^6 (Σ_{i∈T} g_{ij})`. Expand the product as a sum over maps
`φ:[6]→[6]` with `im(φ)⊆T`:
`Σ_{T⊆[6]} (−1)^{|T|} ∏_j (Σ_{i∈T} g_{ij}) = Σ_φ (∏_j g_{φ(j),j})·Σ_{T ⊇ im(φ)} (−1)^{|T|}`.
The inner alternating sum is `Σ_{S⊆[6]∖im(φ)} (−1)^{|im φ|+|S|} = (−1)^{|im φ|}(1−1)^{6−|im φ|}`, which is `(−1)^6`
if `im(φ)=[6]` and `0` otherwise. Surjections `[6]→[6]` are exactly bijections, so only permutations `σ∈S_6`
survive, each contributing `∏_j g_{σ(j),j}`. Hence the inner sum is `perm(g)`, giving the identity. ∎
*(For support `s>6` the same argument gives a sum over **surjections** `[s]↠[6]` in place of the permanent.)*

### The sharp bound
`|D7(c)| = ∏_{j=1}^6 |A(c_j)| · |perm(ζ^{−c_j i})|`. The permanent is a sum of `6! = 720` products, each of
modulus `∏_j |ζ^{−c_j σ(j)}| = 1`, so `|perm| ≤ 720`. And `|A(r)| = |1−ζ^{−r}|/(2π) = sin(πr/7)/π`, maximised
at `r∈{3,4}` (`sin(3π/7)/π`). Therefore `|D7(c)| ≤ 720·(sin(3π/7)/π)^6`.

**Equality** requires *both* `|perm| = 720` and every `|A(c_j)|` maximal. `|perm| = 720` (all 720 unit terms
in phase) forces the matrix rank-1, i.e. all six columns equal, i.e. `c` constant; then
`perm = 6!·∏_{i=1}^6 ζ^{−c·i} = 6!·ζ^{−21c} = 6!` (since `21 ≡ 0 mod 7`), so
`D7((c,…,c)) = 720·A(c)^6` exactly. And `|A(c_j)|` maximal forces `c_j∈{3,4}`. Together: `c≡3` or `c≡4`. ∎

### The ruler value
`|Re D7| ≤ |D7| ≤ 0.643` does **not** give the tighter `0.1431` (which is a phase fact). The finite check over
all `6^6 = 46656` cosets gives `max|Re D7| = 0.14310` at `c≡4` (and `c≡3`, its conjugate), equal to
`|Re(720·A(4)^6)|` — a finite, exact certificate of HYP-2646's ruler bound.

### Corollary — the kernel annihilates coset-constants: `Σ_c D7(c) = 0`
The permanent is multilinear in its columns and column `j` of `ζ^{−c_j i}` depends only on `c_j`, so summing
the identity over `c∈(F_7^*)^6` decouples:
`Σ_c D7(c) = Σ_{σ∈S_6} ∏_{j=1}^6 (Σ_{r∈F_7^*} A(r)ζ^{−r σ(j)}) = 6!·∏_{m=1}^6 B(m)`, `B(m):=Σ_{r=1}^6 A(r)ζ^{−rm}`.
Now `B(m) = (1/2πi)[Σ_{r=1}^6 ζ^{−rm} − Σ_{r=1}^6 ζ^{−r(m+1)}]`, and `Σ_{r=1}^6 ζ^{−rm} = 7·[7∣m] − 1 = −1`
for `m∈{1,…,6}`. So for `m∈{1,…,5}` (both `m, m+1 ∈ {1,…,6}`) the bracket is `(−1)−(−1)=0`, i.e. **`B(m)=0`**;
only `B(6) = (1/2πi)[−1 − 6] = −7/(2πi)` survives. Since `B(1)=0`, the product vanishes:

> **`Σ_{c∈(F_7^*)^6} D7(c) = 0`** — the support-6 kernel has **zero mean** over the coset space.

**Consequence for the open tail.** Because `Σ_c D7(c) = 0`, the correction is *centered*:
`corr(E) = Σ_c D7(c)·S_c(E) = Σ_c D7(c)·(S_c(E) − S̄)` for any constant `S̄`. The wide-spread correction
therefore sees **only the coset-variation** of the reciprocal sums `S_c(E)`, not their average — precisely
the cancellation the summation-by-parts / cotangent–Dedekind route (HYP-2614, THM-504-D) is built to exploit.
This is the structural reason the *signed* series converges where the absolute envelope diverges (MISTAKE-078).

## Why it matters

- It makes HYP-2613's "the exact support-6 layer is a six-sector root-of-unity permanent" **precise and
  proved**, and HYP-2646's `|Re D7| ≤ 0.1431` a closed-form/finite-check theorem rather than a measurement.
- The extremal structure is the mechanism: the **constant coset** (all six coordinates in one residue class
  mod 7 — the maximally-coherent support-6 relation) is where the permanent degenerates to rank-1 and the
  kernel peaks at `720·A(3)^6`. This is the exact algebraic weight the summation-by-parts / cotangent–Dedekind
  route (HYP-2614, THM-504-D) needs for `corr(E) = Σ_c D7(c)·S_c(E)`.
- **Scope (honest).** This closes the *algebraic-weight* structure of the support-6 kernel. It does **not**
  close the open wide-spread lemma itself — the *signed reciprocal hyperplane tail bound* on
  `Σ_c D7(c)·S_c(E)` uniform over wide `E` after finite low-height wall deletion (the project's single named
  open lemma; MISTAKE-078, HYP-2614, HYP-2644) — which requires bounding the conditionally-convergent
  reciprocal sums `S_c`, not the weights. Guardrail (HYP-2614): the one-coordinate residue marginal does NOT
  vanish (the permanent is column-multilinear but `Σ_{c_1} A(c_1)·[col 1]` need not be 0), so the cancellation
  must run on the integer relation hyperplane, not per residue coordinate.

## Files

`04-computation/lrc14_support6_permanent_kernel_kps_S127.py` (+ `.out`): the identity verification (max err
1.68e-16), the exhaustive sharp-bound / ruler check (46656 cosets), and the constant-coset closed form.
