# THM-470: the algebra ladder above (sign, v₂) — interval structure, the coarsening collapse, and the t=7 wall

**Status:** parts A–B PROVED (short proofs below); part C COMPUTED/VERIFIED
(per-witness independent re-verification; scripts
`04-computation/erdos592_uniform_rungs_kp0611.py`,
`04-computation/erdos592_invariant_wall_kp0611.py`, outputs in
05-knowledge/results/). Part C's master verdict: see inline status.
**Source:** kind-pasteur-2026-06-11-S1 (HYP-2392/2393/2396; continues THM-465 B,
POKE Task 2).

**Setting.** A *feature map* F assigns each lex-ordered pair x <lex y in ℕ³ a
class, independently of any ambient size ("t-independent"); F is *gap-determined*
if F(x,y) = φ(y−x). A *rule* is a set T of classes (blue iff F ∈ T); a rule *wins
at t* if it is triangle-free and hits every binary subgrid on [t]³. An
*F-measurable strong witness* is an infinite rule winning on all of ℕ³ (kills
every binary subgrid of the full grid, hence every full-type set, THM-453 C).

## A. Interval structure and the König equivalence (PROVED)

**A1 (antitonicity).** For any t-independent F: a table winning at t+1 wins at t.
(Restriction: subgrids and triangles of [t]³ are subgrids/triangles of [t+1]³,
and the hitting edge of a subgrid lies inside it.) So the SAT region of every
algebra is an interval [3, t_dead(F)).

**A2 (König equivalence).** An F-measurable strong witness on ℕ³ exists ⟺ the
F-instance is SAT at EVERY finite t (t_dead = ∞).
Proof. (⟹) restriction. (⟸) For each t the set of winning tables, restricted to
the (finite) set of classes realizable at t, is finite, nonempty, and coherent
under restriction (A1's argument). König's lemma yields a thread, i.e. one table
defined on all classes whose restriction wins at every t; its rule graph on ℕ³
is triangle-free (every triangle lives in some [t]³) and hits every binary
subgrid (it lives in some [t]³). ∎
Consequence: the "conjoined-t" SAT of THM-465 B is equivalent to the single
largest-t instance; per-size death anywhere kills the whole rung. The
infinite-witness question for ANY t-independent algebra is a question about
one threshold t_dead(F).

**A3 (coarsening collapse).** If F' refines F (F-classes are unions of
F'-classes), every F-rule is an F'-rule, so t_dead(F) ≤ t_dead(F'); in
particular any algebra FACTORING THROUGH (sign, v₂) — comparisons, tail
conditions, unions — inherits THM-465 B's refutation. Climbing the ladder
requires refinement. The finest gap-determined algebra is Finv: φ = identity on
the gap vector (= the translation-invariant game of THM-453 F at n=3);
t_dead(F) ≤ t_dead(Finv) for every gap-determined F.

## B. The rungs computed (VERIFIED witnesses at every SAT cell)

* **F2J (dyadic 1-jet)**: (sign, v₂(|d|), next binary digit) per coordinate.
  Per-size: SAT (3,4) 474 edges, (3,5) 1740, (3,6) 4008 — all independently
  re-verified — and **feature-UNSAT at (3,7)** (665 classes, 516 s).
  HYP-2392 REFUTED: t_dead(F2J) = 7.
* **Corollary (new even for THM-465): t_dead(F2) = 7.** F2J refines F2, so the
  bi-dyadic algebra dies PER-SIZE at t=7 — the conjoined refutation of THM-465 B
  was the shadow of a per-size wall one size past the (3,6) frontier.
* **F2X (cross-gap / Larson-flavored)**: sgn_v2 of (d₁,d₂,d₃,d₁−d₂,d₂−d₃,d₁−d₃).
  Per-size: SAT (3,4) 441, (3,5) 1303, (3,6) 4062 (verified); (3,7): see .out
  (master verdict below subsumes it).
* **Control:** F2 conjoined t=4..7 re-derived feature-UNSAT (0.3 s) — environment
  reproducibility vs mac-mini's run.

## C. The t=7 wall (master experiment)

Question: t_dead(Finv) — the full translation-invariant algebra at (3,7);
by A3 this bounds EVERY gap-determined rung (jet towers of any depth, cross-gap,
leading-digit, all of them at once).

**VERDICT: [PENDING — run in flight at write-up; see
`05-knowledge/results/erdos592_invariant_wall_kp0611.out`; the session-close
message will carry the result. If feature-UNSAT: no translation-invariant strong
witness exists on [7]³ and the entire gap-determined ladder is decided (alive
t ≤ 6, dead t ≥ 7); if SAT: the wall is jet-specific and the invariant witness's
structure is the next rung.]**

**HYP-2396 (the 2n+1 conjecture), claimed on this evidence:** R(n,2) = 2n+1.
Data: R(1,2) = 3, R(2,2) = 5 (THM-459), and at n=3 the gap-determined wall sits
at 7 with matching lower bound Q(3,6) SAT (THM-465 A); Q(n, 2n) SAT also at
n = 1, 2. If true: strong witnesses NEVER exist (any n), HYP-2363 is false,
Specker's witnesses are necessarily non-strong, and the constructive program
must either go non-invariant (value-dependent features — not gap-determined, so
outside A3's collapse) or weaken "strong" (kill shapes between binary subgrids
and full type — the THM-460 tower ladder becomes the natural home).

## Honesty

- The wall (even if confirmed) is about FEATURE algebras/translation-invariant
  rules; free Q(3,7) remains undecided (raw instance timed out in S2). The
  invariant = free cutoff equality is proved only at n = 2 (THM-453 G).
- A2's König needs t-independence of the feature map; value-dependent algebras
  satisfy it too (classes of a pair don't change with ambient), so the framework
  covers the non-invariant route — only A3's Finv bound stops applying.
- All feature-UNSAT verdicts: THM-465 semantics (algebra too coarse), no negative
  claim about Q(3,t) itself.

**Cross-refs:** THM-465 (the rung below), THM-469 (which refinements are not
instantly dead: sum-free gradings), THM-453 F/G (invariant game, n=2 equality),
HYP-2392/2393 (resolved here), HYP-2396 (opened here).
