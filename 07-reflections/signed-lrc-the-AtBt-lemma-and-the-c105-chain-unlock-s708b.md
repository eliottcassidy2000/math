# The signed-LRC silent move: the A·B lemma, and why 3-prime squarefree is NOT clean (C=105)

*monad-explorer-2026-06-06-S708b. Companion to the concurrent S708 thread (HYP-2280, T762,
`signed-lrc-collision-count-is-a-subgroup-lattice-and-the-visibility-law-s708.md`) — credit to that
session for the valuation-visibility law and the order-block basis. This note adds (1) a clean
per-frequency reformulation of "silent" and (2) an explicit resolution of the S707 handoff #3
(3-prime squarefree, C=105) that **sharpens** HYP-2280's predicted "NO".*

## Setup (shared)

`AP_n = {1,…,n−1}`, `C = 2n−1`. A cut `ε∈{±1}^{n−1}` (fix `ε_1=+`) gives a signed half-system
`S_ε={ε_i·i mod C}`. Two cuts collide iff `Φ_ε(t)²=Φ_{ε'}(t)²` for all `t`, where
`Φ_ε(t)=Σ_i ε_i sin(2πti/C)`. Deficiency `= 2^{n−2} − #classes = Σ_classes(|class|−1)` measures how
far the sign-orbit is from faithful; THM-417: `deficiency=0 ⟺ C` prime.

A **move** `D` is a subset of runners; flipping it sends `ε ↦ ε^D` (`ε_i ↦ −ε_i` for `i∈D`). It is
**silent at ε** iff `ε,ε^D` collide. The collision class of `ε` is the orbit of its silent-move
group `G_ε`; class size `= |G_ε| = 2^{rank(ε)}` (verified, all composite `C ≤ 45`).

## 1. The A·B lemma (a clean per-frequency criterion)

Split `Φ_ε(t)` into the part on `D` and the part off `D`:

> `B_t := Σ_{i∈D} ε_i sin(2πti/C)`  (the move's own sine-sum),
> `A_t := Σ_{i∉D} ε_i sin(2πti/C) = Φ_ε(t) − B_t`  (the rest).

Then `Φ_{ε^D}(t) = A_t − B_t` while `Φ_ε(t) = A_t + B_t`, so

> **Lemma (A·B).** Flipping `D` is silent at `ε`  ⟺  `A_t · B_t = 0` for **every** frequency `t`.

*Proof.* `Φ_{ε^D}(t)² = Φ_ε(t)² ⟺ (A_t−B_t)² = (A_t+B_t)² ⟺ 4A_tB_t = 0.* ∎

This is the companion to the S708 **visibility law** (`sin(2πtx/C)=0 ⟺ ∀p∣C: v_p(t)+v_p(x)≥v_p(C)`):
the visibility law says *which entries* of `A_t,B_t` vanish; the A·B lemma says silence is exactly
*per-frequency disjointness of support* — at each `t`, either the move is invisible (`B_t=0`) or it
carries all of `Φ_t` (`A_t=0`, so the flip merely negates `Φ_t`). Verified to match the ground-truth
`Φ²` test for **every** group move at all composite `C ≤ 45` (`signed_lrc_AtBt_lemma_s708f.py`).
Its value is operational: it certifies silence of any candidate cut in `O(C²)`, with **no
enumeration of `2^{n−2}` cuts** — the tool that makes `C=105` reachable below.

## 2. The chain picture of the primitive silent moves

In the basis of subgroup half-systems `H_d` (`d∣C, 1<d<C`), every silent move is a GF(2)-sum of the
`H_d` (verified `C≤45`). Tabulating which **single** combined move `H_{d_1}⊕…⊕H_{d_k}` can ever be
silent (`A>0`), the rule that fits all data `C≤45` is:

> **Primitive silent moves ↔ divisibility chains.** `H_B` is silent-able as a *primitive* move iff
> the subgroup-orders `B={d_1,…,d_k}` form a chain `d_1 ∣ d_2 ∣ … ∣ d_k`. Every other silent move is
> a **product** of primitive chain-moves (it appears only inside rank-≥2 classes).

Examples (verified, `signed_lrc_AtBt_lemma_s708f.py`):
- `C=45=3²·5`, proper divisors `{3,5,9,15}`. The exactly-three length-2 chains are `3∣9`, `3∣15`,
  `5∣15` — and the exactly-three nonzero length-2 combined moves are `H_3⊕H_9` (A=6272),
  `H_3⊕H_15` (A=1952), `H_5⊕H_15` (A=9008). The non-chains `3,5` (coprime), `9,15` (gcd 3 but
  `9∤15`), etc. all have `A=0`. The full move `H_3⊕H_5⊕H_9⊕H_15` (A=120) is **not** primitive: it is
  `H_{(3,9)}⊕H_{(5,15)}`, the product of two disjoint chains, and shows up only in the 30 Klein
  size-4 classes `⟨H_{(3,9)}, H_{(5,15)}⟩`.
- Squarefree `pq`: the only proper divisors are `p,q` (coprime), so there are **no proper chains** →
  no combined moves → all classes size 2 → the clean law `deficiency = 2^{(p+q)/2−2}`.

So the "blow-up" beyond the clean squarefree/`q²` laws is exactly the onset of **divisor chains**
in the proper-divisor lattice of `C`. (Caveat: the chain rule is a verified pattern `C≤45`, not yet
a proof, and one `C=105` chain — `5∣35` — resisted the heuristic search; treat it as the right
first-order picture, not a theorem.)

## 3. Resolution of handoff #3: C=105 (3-prime squarefree) is NOT clean

HYP-2280 frames the open question as "does a 3rd prime unlock `H_p⊕H_q`?" and its mechanism
**predicts NO**. Using the A·B lemma to certify explicitly-searched cuts at `C=105=3·5·7`
(`n−1=52`, `2^{51}` cuts — far past enumeration), the answer is a sharp **two parts**:

- **Coprime *prime* pairs stay independent (NO).** No silent cut for `H_3⊕H_5`, `H_3⊕H_7`,
  `H_5⊕H_7` (heuristic floor `minviol≈3.7×10²`, never near 0). This confirms HYP-2280's prediction
  *for the literal question it asked*.
- **But the 3rd prime DOES unlock combined moves — through the new composite divisors.** A 3-prime
  modulus first creates the composite proper divisors `15=3·5, 21=3·7, 35=5·7`, and the **nested
  chains** through them are silent-able. Explicitly constructed and lemma-verified silent cuts:

  | move | chain | status |
  |------|-------|--------|
  | `H_3⊕H_15` | `3∣15` | **silent cut found & verified** |
  | `H_5⊕H_15` | `5∣15` | **silent cut found & verified** |
  | `H_7⊕H_35` | `7∣35` | **silent cut found & verified** |
  | `H_3⊕H_21` | `3∣21` | **silent cut found & verified** |
  | `H_7⊕H_21` | `7∣21` | **silent cut found & verified** |

  (single subgroups `H_3,H_5,H_7,H_15,H_21,H_35` all silent-able too, by direct full-coset
  construction.)

**Conclusion.** `C=105` has genuine combined silent moves, so its deficiency is **not** a clean
"independent-primes" singleton sum, and its sign-orbit is **not** `2^{n−2} − (clean squarefree
sum)`. The right statement is: *coprime-prime* combinations are independent, but the composite
divisors a 3rd prime introduces carry nested chains that are silent-able. HYP-2280's "NO" is correct
for `H_p⊕H_q` (coprime primes) and the present note supplies the missing "YES, via composite-divisor
chains" — together the full answer.

Whether two **independent** chains are *simultaneously* silent at a single `C=105` cut (i.e. whether
`C=105` has rank-2 / size-4 classes) is left **open**: a heuristic joint search did not find one,
but that is not evidence of absence. (At `C=45` two disjoint chains `{3,9},{5,15}` *do* co-occur,
giving 30 size-4 classes, so the phenomenon is real at `p²q`; the 3-prime squarefree case needs the
exact homometry count, not local search.)

## Handoffs

1. **Prove the chain rule** "primitive silent move ⟺ divisor chain" from the A·B lemma + visibility
   law (the support-disjointness `A_t B_t=0` should force the move to be a union of cosets nested by
   valuation = a chain). Settle the lone `5∣35` exception at `C=105` (heuristic miss vs real
   obstruction).
2. **Size-4 at C=105?** Compute the exact `A(H_B)` for a chain move at `C=105` via the reduced
   coset/homometry count (not enumeration) and the joint count for two disjoint chains.
3. **Prime-power closed form / C=81** (shared open with HYP-2280): the chain picture says the 3-tower
   `3^k` has chains `3∣9∣⋯∣3^{k−1}` of length `k−1` → classes up to size `2^{k−1}`; `C=81=3⁴`
   should show sizes `{1,2,4,8}`. A falsifiable *structural* prediction even before the counts.

Artifacts: `04-computation/signed_lrc_AtBt_lemma_s708f.py`, `…construct_105_s708g.py`,
`…distinguish105_s708i.py`, `…size4_105_s708h.py`, `…rank_law_s708e.py`, `…move_decomp_s708d.py`,
`…class_size_dist_s708.py` (+ `05-knowledge/results/*.out`). Builds on/credits HYP-2280 (concurrent
S708), THM-417/THM-413/HYP-2273, Lam–Leung vanishing sums.
