# Independent audit of THM-4864 (reversal_cube_mod4_anchor_20261104)

Referee: independent Claude referee subagent with a fresh context. I did not write the packet. Date 2026-11-04.

Objects reviewed, pinned at commit `97632252db`:
- `PROOF.md`, sha256 `020f999e…b02fef`;
- `01-canon/theorems/THM-4864-…md`, sha256 `73bb733d…ce9744`.

The cited inputs were read only to check scope:
- parity_profile_transfer §§2–3;
- prime_seven_frontier §3;
- highest_quotient_reduction §§2–3;
- S189 collision §§4–5;
- two_port_parity_separation E1;
- odd_pair_quotient §2;
- even_endpoint_obstruction §1.

No packet code was imported or copied. I grepped `confirm.py` only to check its bookkeeping.

## Verdict: SOUND WITH CORRECTIONS

- **No MAJOR issue.** Lemmas 1–3, Proposition 4, Theorem A, Theorem B and the Corollary are correct.
- **Required changes.** Three MINOR corrections (a missing hypothesis, the evidence bookkeeping, and one inaccurate
  sentence about evidence) and some NOTE-level wording.
- **Independent computations.** Every claimed identity and functional value was recomputed. That covers two
  from-scratch C++ engine modes plus Python path enumeration, and about 25,000 individual checks. There were
  **0 failures**.

## 1. Analytic check, step by step

**(1) and the odd-cycle expansion.**
- The port cycles split into four disjoint classes: x-only, y-only, x→y→P→x, and x→P→y→P′→x. No other class
  exists because y→x is impossible.
- Each literal is right, and each parity condition is right once the ports are given weight 1.
- Singleton terms such as (1+a_s)a_s vanish identically.
- Modulo 4 the expansion only needs cycles of length ≥ 3, and a tournament has no others.

**Lemma 2.**
- I re-derived both coefficients. The only monomials that could alias are excluded by disjointness, namely the
  cases s=e′ and s′=e.
- The cubic coefficient arises only from P={v} in the double sum.

**Converse duality.** Reversing the arcs and swapping x with y preserves port weights, because both ports weigh 1.

**Lemma 1.**
- *Segment step.* The S_r action is free, because the cycle is linearised from a distinguished port and the
  segments are disjoint and nonempty.
- *Several modules.* The runs between segments stay fixed. So taking the least-indexed module with at least 2
  segments gives a well-defined free action.
- *One segment.* Rédei gives the parities: even subsets ≡ 1 and odd subsets ≡ 0 when |A_i| ≥ 2. A singleton
  has only odd subsets.
- *Interiors.* Interiors affect only the constant term, through the port-free cycles. This is why arbitrary
  same-size target interiors are enough: every functional used omits the constant.
- *Anchored use in Theorem B.* beta is a module of the augmented W and of the augmented W′, because u is an end of
  both A_j and A′_j. The other children are modules on both sides, and the weights agree.

**Lemma 3 and its consequence.** Both are correct. The consequence telescopes as used.

**Proposition 4.** All of the following are correct:
- the N-terms give O_v;
- the Π-terms reduce to Σ eps_v(S)f(S), with Rédei applied to P′;
- the two cubic bijections (prepend v; complement then append v);
- g(v) ≡ w(v) for n ≥ 3;
- the table f(S) = w(S) + ω(S)[|V∖S| = 1], checked in all 4 cases;
- the identity Σ_{S∋v,x; S≠V} eps_v(S) ≡ eps_v(V∖x);
- the cancellation of the correction sum.

At n = 2 one gets L_v = g(v) = w(u)+w(v), so the stated boundary is real.

**Theorem A.**
- n = m ≥ 3 as required.
- L_z = 1 forces a nonzero nonconstant coefficient. Möbius inversion then puts a residue-2 query inside the L_z
  battery, which has at most 4m queries.

**Theorem B.**
- **Reversed type.** The target kernel is Z^op, and Proposition 4 at u applies with n = m+1.
- **Unchanged type.**
  - Y^op = Z_u, so g_{Z_u} = g_Y∘J.
  - The τ-correspondence keeps the vertex set, hence the weight. It permutes only the port literals of u and beta,
    and τJ = Jτ. So the orbit sums agree.
  - The monomial set of L_z is τ-closed, because u and beta have the same arc to z.
  - For u = b_1, mono(L_u) ∪ mono(L_beta) ∪ {a_u a_β c_β} is τ-closed.
- **Hard case, h_ub = p_j(Q).** Every path of Z−beta with u at an end has weight 1, and Z−beta ≅ Q.
- **The p_j = 1 branch.** With zero weights, (1) gives [a_j]δ_Q = (S_j+1)+(E_j+1) = p_j.
- **k = 2.** A port cycle using the u–β arc has weight ≡ number of ports. So the odd ones have exactly one port and
  pure literals, while L_u is mixed.
- **Case split.** The cases are exhaustive: reversed | unchanged non-hard | hard k ≥ 3 (anchored if p=0, coarse if
  p=1) | hard k = 2.

**Corollary.**
- *Lift.* It is used within scope. The q_2 and q_3 residues depend only on the augmented order, and c_1 is odd. No
  equality of child profiles is needed at M.
- *Restoration.* prime_seven_frontier §3 never looks inside M. The ports copy M's exterior incidences, so every
  other highest node N stays a module, and N has equal profiles by hereditary equality. Anchored queries are
  therefore admissible.
- *Transitivity.* By module_decomposition, H(U[A]) = 1, so U[A] is transitive.
- *Ends.* A TT2 or TT3 always has b_1 or b_k at an end of its target order. For TT4 the only residue is target ends
  {b_2,b_3}. The ≥ 5 case follows from S189 §5 (U[A] ∈ {A, A^op}).

## 2. Issues

**MINOR-1. Missing hypothesis m ≥ 2 in Theorem B.**
- The hypothesis is absent from PROOF §4 and from THM-4864 Statement item 4 and title (B).
- The proof uses n = m+1 ≥ 3 in Proposition 4.
- At m = 1 the "unchanged" case is W = W′, so it cannot be separated. Check A7 recomputes this.
- The Corollary is unaffected, since m ≥ 7 there.

**MINOR-2. The FINITE-EXACT layer is not recorded as claimed (state at audit time).**
- PROOF Status and README say the outcomes are "appended" to PREREGISTRATION.md, but its Outcomes section is empty.
- `confirm.out` is untracked and stops at "P4 done", with no P5 line.
- `confirm.json` does not exist, yet README and the THM-4864 `scripts:` field list it.
- The THM-4864 status says "FINITE-EXACT pre-registered confirmations (P1-P5) are in the packet".

**MINOR-3. One evidence sentence is inaccurate.**
- PROOF §6 says: "The FINITE-EXACT evidence for those [non-transitive children, TT4 crossed] is in
  PREREGISTRATION.md, with 0 blind cases". THM-4864 Scope says the same.
- The 2,272-comparison sweep is labelled exploratory ("NOT pre-registered") and uses the weighted model.
- It covers transitive children and their reorderings only. HYP-9456 itself says "Non-transitive children have not
  yet been tested".
- The packet's own `residual_c3` probe exhibits a **modulo-4-blind** non-transitive configuration: Paley-7 with C3
  children, unchanged.
- Related: the title of HYP-9456 still conjectures a one-child **modulo-four** separator for every residual child,
  which that same probe contradicts. Its body already revises the conjecture to modulo 8.

**NOTE-1. The mirror case u = b_k in Theorem B(b).**
- For u = b_k one has beta→u. Then a_u a_β c_β already lies in L_β, so the literal "union" is not τ-closed.
- The mirror reduction is valid: apply the u = b_1 case to (W^op, W′^op) and pull back by J, which maps queries of
  the stated form to queries of the stated form.
- The THM-4864 proof summary states L_u+L_β+[a_u a_β c_β] = 1+p_j without restricting to u = b_1. That remains true
  for u = b_k as a sum mod 2, because the symmetric difference is τ-closed. The τ-completion, which uses
  [a_u a_β c_u] instead, gives the same value. Both forms are verified in M6 and A3.
- Suggest one sentence.

**NOTE-2. The symbol [alpha_j] in §4(b) is undefined.**
- It means the coarse coefficient [a_j] in the zero-weight model.
- The cross-reference "also parity_profile_transfer (1)" is a consistency remark only, since that result assumes
  equal child profiles. The self-contained route through Lemma 1 and (1) at zero weights is what is valid here.
- Suggest writing that one extraction line.

**NOTE-3. Lemma 1 step 1 with several modules.** Say "least-indexed module with ≥ 2 segments" (see §1).

**NOTE-4. Theorem B never uses transitivity.**
- The proof only needs u to be a source or a sink of both A_j and A′_j.
- Check A6 confirms the widened statement on non-transitive A_j.
- The widened statement would exclude, at a highest changed node of a collision, any child with a vertex that is an
  end (source or sink) of both T[A] and U[A]. This is optional.

**NOTE-5. The §6 "Sharpness" remark.**
- "Blind to all coarse queries modulo 4" is true for the nonconstant coarse part, or in S188's equal-interior
  setting.
- With this packet's arbitrary same-size interiors it fails: take Q8 with a TT3 child in the source and a C3 child
  in the target. Then H(W) − H(W′) ≡ 2 mod 4, so the (0,0) query separates (`aud_sharpness_constant.out`).
- Suggest qualifying the remark.

**NOTE-6. Minor wording.**
- README says the anchored query is "named only at the source", but Theorem B also anchors at the sink b_k. The TT3
  case of the Corollary needs that.
- In §1, "directed odd cycles" should say "nontrivial".
- The battery bounds hold, with room: at most 4m queries for L_z and at most 7(m+1)−2 for Λ′.
- The bound m ≥ 3 in Theorem A is sufficient but not sharp. At m = 2 with one singleton, L_z = 1 still holds.

## 3. Independent computations (all 0 failures)

The engine `aud_engine.cpp` has two modes:
- (H) exact Hamilton counts of the augmented tournament by subset DP;
- (C) the brute-force parity of odd-weight port cycles, by a subset DP from each port. It never uses formula (1).

The driver modules are `aud_lib.py`, `aud_model_checks.py` and `aud_actual_checks.py`.

| check | scope | cases |
|---|---|---|
| M1: (1) equals the brute-force g_Z on **all** queries | every tournament and every weight vector for n ≤ 4; 40 random tournaments at n = 5 | 1,138 instances |
| M2: Lemma 2, both coefficients | n = 3..6 | 4,120 |
| M3: converse duality, all queries | n = 2..5 | 120 |
| M4: Lemma 3 and its consequence | n = 1..8 | 2,228 |
| M5: Proposition 4, L_v = w(v) | exhaustive at n = 3, 4 (192 + 4,096); random at n = 5..9 (7,700) | 11,988 |
| M5: n = 2 boundary | L_v = w(u)+w(v) in 16/16 cases; L_v ≠ w(v) in 8/16 | 16 |
| M6: weighted kernels (m = 2..5) | τ-orbit equality on **every** monomial orbit (500 kernels); k = 2 correction has no mixed monomial (185); reversed L_u = 1 (500); hard k ≥ 3 Λ′ = 1+p_j for u = b_1 / b_k in both forms (84 / 91); p_j = 1 gives coarse [a_j] = 1 (124); hard k = 2 L_u = 1 (185); non-hard L_z = 1 (223 × 2) | 2,290 |
| A1: actual blow-ups (weight-0 vertex → random 2- or 3-vertex child; independent random target interiors), Proposition 4 through actual F | n = 3..8, ≤ 13 names | 656 L_v; 132 Lemma 1 |
| Lemma 1 on **all** queries, actual F vs the model | 24 coarse instances; 32 anchored instances, source and target | 56 |
| A2: Theorem A on actual graphs | m = 3..6, independent random target interiors | 109 L_z = 1; 80/80 residue-2 witnesses |
| A3: Theorem B on actual graphs | m = 2..6, k = 2,3,4, u = b_1/b_k, 4 target types; see breakdown below | 680/680 witnesses |
| A4 / A5: ambient lift / lift with restoration of a second changed node N → N^op | half parity equals the bare core at every battery query | 28 / 28 |
| A6: non-transitive A_j with a common source/sink (generalisation) | every prediction holds | 48 contrasts |
| Q8 (m = 8, coarse-blind quotient): actual graphs, u = b_1 | k = 2 at all j, k = 3 at j ∈ {0,3,5,7}; unchanged and reversed; every predicted value equals 1, with a witness | 24 |

The four target types in A3 are: unchanged, reversed, a random transitive order with u at the same end, and a random
transitive order with u at the other end. Its results were:
- reversed type: L_u = 1 in 340 cases;
- unchanged type, hard, k ≥ 3: Λ′ = 1+p_j in 108 cases (p = 0: 66; p = 1: 42), and coarse [a_j] = 1 in all 42
  p = 1 cases;
- unchanged type, hard, k = 2: L_u = 1 in 62 cases;
- non-hard: L_z = 1 in 236 + 236 cases.

The outputs are in `aud_model_checks.out` (seed 20261104), `aud_actual_checks.out` (seed 4864),
`aud_lemma1_allqueries.out`, `aud_q8_anchor.out` and `aud_sharpness_constant.out`.

## 4. Required corrections

1. **m ≥ 2 (MINOR-1).** Add "m ≥ 2" to Theorem B in PROOF §4 and to THM-4864, both Statement item 4 and the title
   clause (B).
2. **Evidence bookkeeping (MINOR-2).**
   - Finish the confirmatory run and commit `confirm.out` and `confirm.json`.
   - Append the P1–P5 outcomes to `PREREGISTRATION.md`.
   - Otherwise, change the PROOF, README and THM statements to "pending".
3. **Evidence sentence (MINOR-3).**
   - Restrict the "0 blind cases" sentence in PROOF §6 and in the THM-4864 Scope to transitive reorderings, labelled
     as an exploratory weighted-model result.
   - State that non-transitive children are untested by that sweep, and that `residual_c3` contains a mod-4-blind C3
     configuration.
   - Align the title of HYP-9456 with its modulo-8 revision.
4. **Recommended wording (NOTE-1, NOTE-2, NOTE-5).**
   - Add the u = b_k sentence (NOTE-1).
   - Define [alpha_j] (NOTE-2).
   - Qualify the sharpness remark (NOTE-5).

## 5. Reproduction

    g++ -O2 -o aud_engine.exe aud_engine.cpp        # GCC 15.2.0 used; Python 3.13.14
    python3 aud_model_checks.py 20261104 > aud_model_checks.out      # ~25 s
    python3 aud_actual_checks.py 4864   > aud_actual_checks.out      # ~100 s
    python3 aud_lemma1_allqueries.py    > aud_lemma1_allqueries.out
    python3 aud_q8_anchor.py            > aud_q8_anchor.out          # ~100 s
    python3 aud_sharpness_constant.py   > aud_sharpness_constant.out

The engine binary is a build artifact. It was deleted after the runs.
