# THM-467: The Negative OCF Evaluation — I(Ω(T), −2) Is the Reciprocal Determinant Quotient, and x = −2 Is Invisible to the Omega Involution

**Type:** Theorem (PROVED, parts (a)–(d)) + census + refutation table
**Certainty:** 5 for (a)–(d) (complete elementary proofs + exhaustive verification n ≤ 6); census facts labeled separately
**Status:** PROVED (identities) / VERIFIED (census ranges) / REFUTED (all 31 candidate closed forms)
**Added by:** kind-pasteur-2026-06-10-S2 (Thread C, HYP-2380)
**Tags:** #ocf #independence-polynomial #negative-fugacity #omega-involution #redei-berge #nilpotent-ring #determinant

---

## Setting

THM-002 (PROVED): H(T) = I(Ω(T), 2), Ω(T) = conflict graph of directed odd
cycles, I(G,x) = Σ_k α_k x^k. HYP-2380 asked: x = −2 is the odd reflection of
the proven evaluation point — what does I(Ω(T), −2) mean?

Throughout: A = adjacency matrix of the tournament T on n vertices,
R := Z[x_1,…,x_n]/(x_i²) (square-free/nilpotent ring), X := diag(x_1,…,x_n),
ε : R → Z the linear functional sending every square-free monomial to 1
(sum of all coefficients), and

```
C_odd(x) := Σ_{directed odd cycles c of T} Π_{i ∈ c} x_i  ∈ R.
```

## Statement

**(a) (Reciprocal determinant quotient — the exact identity at x = −2.)**
For every tournament T and every even integer t:

```
I(Ω(T), t) = ε( [ det(I + XA) / det(I − XA) ]^{t/2} )      (in R)
```

In particular (t = 2 recovers OCF; t = −2 is the new point):

```
H(T)        =  ε( det(I + XA) · det(I − XA)^{−1} )
I(Ω(T), −2) =  ε( det(I − XA) · det(I + XA)^{−1} )          ← THE RECIPROCAL
```

Equivalently: **x ↦ −x in I(Ω(T), x) is A ↦ −A** (give every arc weight −1,
so every odd cycle acquires weight −1); there is no unsigned digraph
realization (I(−2) takes both signs, see census).

**(b) (arctanh form.)** In R: C_odd = ½ log[ det(I+XA) / det(I−XA) ], i.e. the
odd-cycle generating function is the arctanh-type (odd) part of −log det(I−XA),
and exp(t·C_odd) = Σ_{K} t^{|K|} x^{supp K} over collections K of pairwise
vertex-disjoint directed odd cycles. The OCF is ε∘exp(2·C_odd).

**(c) (Omega self-duality; why x = −2 is NOT the omega mirror of x = 2.)**
For every tournament: **ω(U_T) = U_T**, where U_T is the Rédei–Berge symmetric
function and ω the omega involution. Consequently the ω-mirror of the Rédei
evaluation H = ζ(U_T) (ζ = principal specialization at the single variable 1,
GS Lemma 6.4/6.5) is H itself — the I-world shadow of THM-061's path-reversal.
The evaluation I(Ω(T), x) equals φ_x(U_T) for the algebra map
φ_x: p_1 ↦ 1, p_{2k+1} ↦ x/2 (k ≥ 1), p_{2k} ↦ 0; the point x = −2 needs
p_1 ↦ +1 with p_{odd ≥ 3} ↦ −1, which is **not** an alphabet/super
specialization (a + b = 1 and a + b = −1 are incompatible) and is not in the
Hopf orbit {id, ω, antipode S} of canonical operations: ω fixes U_T, and
S(U_T) = (−1)^n U_T merely reproduces ζ(S(U_T)) = (−1)^n H. The W-world
reflection (THM-061: W(−1/2) = (−1)^{n−1}H, path reversal) and the I-world
reflection x ↦ −x are therefore implemented by DIFFERENT operations: path
reversal vs. arc-weight −1.

**(d) (Structural corollaries, all n, elementary.)**
1. Strong-component factorization: I(Ω(T), x) = Π_{strong components C}
   I(Ω(T|_C), x) (every directed cycle lies in one strong component;
   cross-component disjointness is automatic). Hence sign(I(−2)) is the
   product of component signs.
2. 2-adic mirror: H + I(−2) ≡ 2 (mod 8) and H − I(−2) ≡ 4α₁ (mod 16) at all
   n; **exactly I(−2) = H − 4α₁ for n ≤ 8** (α₃ = 0 needs ≥ 9 vertices).
   I(−2) is always odd, in particular never 0.
3. Path-GF reciprocity (Irving–Omar + reversal): with ξ_k(T) ∈ R the
   generating polynomial of directed k-vertex paths and W_T(z) = Σ ξ_k z^k:
   **W_T(z)·W_T(−z) = 1 in R[z]** — log W_T(z) is an ODD function of z.

## Proofs

**(a)+(b):** In R, tr((XA)^k) = k·Σ_{directed k-cycles c} x^{supp c} (closed
walks with a repeated vertex die by x_i² = 0; tournaments have no 1- or
2-cycles). The finite tr-log expansion gives log det(I − tXA) =
−Σ_c t^{ℓ(c)} x^{supp c}, so ½[log det(I+XA) − log det(I−XA)] =
Σ_c ½(1 − (−1)^{ℓ(c)}) x^{supp c} = C_odd. Multinomial expansion of
exp(t·C_odd) in the nilpotent ring keeps exactly the square-free products =
collections of disjoint odd cycles, each weighted t^{|K|} (the j! cancels the
ordering). Applying ε gives Σ_K t^{|K|} = I(Ω(T), t). With t = 2 and THM-002
this is H; with t = −2 it is the reciprocal quotient. ∎

**(c):** Grinberg–Stanley Theorem 8.1 [arXiv:2307.05569]: ω(U_D) = U_{D̄} for
every digraph (also Chow 1996, Cor. 2; Wiseman 2007). Loops never enter U
(descents involve distinct consecutive vertices), so for a tournament
U_{T̄} = U_{T^op}. By GS Theorem 1.39, U_T = Σ_K 2^{|K|} p_1^{n−|supp K|}
Π_{C∈K} p_{ℓ(C)} over collections K of disjoint directed odd cycles; reversal
is a type-preserving bijection between collections of T and of T^op, so
U_{T^op} = U_T. Hence ω(U_T) = U_T. The specialization bookkeeping is GS
Lemma 6.4/6.5 (ζ(p_λ) = 1; ζ(U_D) = #hamps(D̄)) plus ham(T^op) = ham(T)
(reversal). The non-specializability of (p_1, p_odd) ↦ (1, −1): a
superalphabet gives p_j ↦ a + (−1)^{j−1} b, forcing a + b to be both +1 (j=1)
and −1 (j ≥ 3 odd). ∎

**(d1):** trivial, see statement. **(d2):** I(2) ± I(−2) from the α-expansion;
Rédei (H odd) gives I(−2) odd. **(d3):** Irving–Omar [arXiv:2412.10572,
proof of Prop. 21]: W_D(z) = (W_{D̄}(−z))^{−1} in R; for tournaments
ξ_k(T̄) = ξ_k(T^op) = ξ_k(T) (loops cannot occur in a path of distinct
vertices; reversal preserves supports). Equivalently a direct sign-reversing
involution on pairs of disjoint paths (move the losing endpoint, tournament
axiom decides). ∎

## What I(Ω(T), −2) is NOT — refutation table (pre-registered candidates)

All tested against the full labeled census n = 3..5 (1096 tournaments), all
56 iso reps at n = 6, and 2000 random n = 7. Smallest counterexamples:

| Candidate | First counterexample |
|---|---|
| I(−2) = ±W(r₀), ±2^{n−1}W(r₀), r₀ ∈ {3/2, 1, 5/2, −3/2, −1} (20 forms) | all die at n=3, transitive (mask 0) |
| I(−2) = det(J − 2A), \|I(−2)\| = det(J−2A) | n=3 transitive: 1 vs 4 |
| I(−2) = ±det(I − A), ±det(I + A) | n=3 (masks 0, 2) |
| I(−1) = det(I − A) | n=4 mask 4: −1 vs −2 |
| I(1) = det(I + A) | n=4 mask 4: 3 vs 2 |
| I(1) = per(I + A) | n=4 mask 4: 3 vs 4 |
| I(−2) = ±per(I − A) | n=3 (masks 0, 2) |
| I(−2) = H − 4α₁ at ALL n | exact n ≤ 8; explicit n=9 cex: C3⊕C3⊕C3 block-transitive has I = (1+x)³, I(−2) = −1 ≠ 15 = H−4α₁ |

(det(I−A) = Σ_{ALL-cycle collections}(−1)^{#cycles} and det(I+A) weights even
cycles by −1: tournament even cycles do NOT cancel in the determinant world —
the n=4 strong tournament, one 4-cycle, separates them.)

**Honest verdict on "meaning":** I(Ω(T), −2) is not a cardinality (sign
varies: +1 transitive, −1 the 3-cycle, −13 RQ5, −131 Paley₇, +9 at
S4⊕S4), and not a det/per/W-evaluation of any tested form. Its exact content
is (a): the OCF evaluation of the (−1)-arc-weighted tournament = ε of the
reciprocal determinant quotient. A standalone counting interpretation of
|I(−2)| remains OPEN.

## Census highlights (kind-pasteur-2026-06-10-S2, Thread C)

- I(2) = H re-verified: ALL labeled n = 3..6 (8 + 64 + 1024 + 32768), 2000
  random n = 7, 400 random n = 8, 9 structured n=7 families, n=9 blocks.
- Ring identities (a)/(b)/(d3) verified as exact R-equalities on every iso
  class rep n = 3..6 plus randoms (94 tournaments).
- Sign of I(−2): never 0; negative fraction grows: 25% (n=3), 62.5% (n=4),
  88% (n=5), 97.6% (n=6), ≈99.5% (n=7 sample), 400/400 (n=8 sample).
  Positives at n ≤ 8 are the cycle-poor classes; via (d1) they factor:
  [1]→1, [1,2,1]=(1+x)²→1, [1,3,2]=(1+x)(1+2x)→3, [1,4,4]=(1+2x)²→9 (S4⊕S4,
  n=8, discriminant 0, double root −1/2).
- GS Theorem 7.1 (H ≡ 1 + 2α₁ mod 4) and the 2-adic mirror congruences hold
  in the entire census.
- α₁² ≥ 4α₂ everywhere tested — as it must: at n ≤ 8 the independence number
  of Ω(T) is ≤ 2, so Ω is trivially claw-free and I is real-rooted
  (Chudnovsky–Seymour). Whether I(Ω(T), x) stays real-rooted for n ≥ 9
  (where Ω can contain claws) is OPEN.
- I(−1) = −χ̃(Ind(Ω(T))) (reduced Euler characteristic of the collection
  complex): broad distribution (NOT confined to 0/±1); max value +1 at all
  n ≤ 8 tested, realized by transitive ([1]) and the disc-0 family [1,4,4];
  value −2 never occurs at n ≤ 6 (α₁ = 3 needs n ≥ 7: S4⊕C3; observed-only,
  not claimed as law).
- Paley₇ vs the other six n=7 circulants: α = (80, 7) vs (59, 14);
  I(−2) = −131 vs −61; H = 189 vs 175.

## Sources pinned (full PDFs read this session)

- Grinberg & Stanley, arXiv:2307.05569: Def. 1.15 (U_D = Σ_w L_{Des(w,D),n});
  Thm 1.31 (p-expansion over S_V(D,D̄), sign (−1)^{φ(σ)}); Thm 1.39
  (tournaments: U_D = Σ_{σ∈S_V(D), all cycles odd} 2^{ψ(σ)} p_{type σ});
  Cor 1.40 (U_T ∈ N[p_1, 2p_3, 2p_5, …]); Lemma 6.4/6.5 (ζ; ζ(U_D) =
  #hamps(D̄)); Thm 6.6; **Thm 7.1 (mod-4 refinement = the OCF mod-4 digit)**;
  Thm 8.1 (ω(U_D) = U_{D̄}, S(U_D) = (−1)^n U_{D̄}).
- Irving & Omar, arXiv:2412.10572: Thm 10 (matrix-algebra master formula),
  Cor 19/20, Prop 21 proof (W_D(z) = (W_{D̄}(−z))^{−1}).
- The "2" in I(Ω, 2): per-cycle factor 2p_ℓ in GS Thm 1.39 (cycle + its
  reversal; even lengths cancel, odd double); the evaluation point is the
  1-variable principal specialization ζ. The "−2" has no such Hopf origin (c).

## Scripts / results

- `04-computation/ocf_negative_eval_kpo2.py` → `05-knowledge/results/ocf_negative_eval_kpo2.out`
- `04-computation/ocf_neg_followup_kpo2.py` → `05-knowledge/results/ocf_neg_followup_kpo2.out`
- `04-computation/ocf_neg_strongcomp_kpo2.py` → `05-knowledge/results/ocf_neg_strongcomp_kpo2.out`

## Open questions spawned

1. Standalone counting interpretation of |I(Ω(T), −2)| (beyond (a)).
2. Real-rootedness of I(Ω(T), x) for n ≥ 9 (Ω acquires claws; CS no longer applies).
3. Negative-fugacity/Shearer (LLL) reading: −1/2 is a root for the strong-4
   block families ((1+2x)^k); is the critical fugacity of Ω(T) tied to strong
   4-subtournaments?
4. Achievable α-vector geography (e.g. [1,3] with α₂ = 0 appears impossible —
   the 3-triangle fan forces a 5-cycle; unproved).
