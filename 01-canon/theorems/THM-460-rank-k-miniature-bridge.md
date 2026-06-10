# THM-460: the tower miniature for ω^(ω^m) — full-type = stacked towers, and the König bridge to the Erdős 592 ladder

**Status:** PARTIAL — parts A–C PROVED (proofs below); part D = computations (statuses
inline). The open problem (m=3) is NOT resolved; this builds its finite probe.
**Source:** mac-mini-2026-06-09-S2 (T768, HYP-2366; continues THM-453)

## Context

Erdős 592's surviving frontier is α_m := ω^(ω^m) for finite m: α_1 = ω^ω is positive
(Chang), α_2 = ω^(ω²) positive (Schipperus, Darby), α_4 = ω^(ω⁴) negative
(Schipperus), and **α_3 = ω^(ω³) is the smallest open case**. THM-453's grid
miniature covered ω^n (n < ω). This file lifts it across the first limit exponent.

## A. Presentation and block structure (standard, fixed here for the code)

Ordinals below ω^(ω^m) ↔ finitely supported functions f: N^m → N (the exponent
δ < ω^m in CNF ↔ an m-tuple; f(δ) = its coefficient), ordered by largest
disagreement: ξ < η iff f_ξ(δ*) < f_η(δ*) at the lex-largest δ* where they differ.
For m-tuples δ, lex order on N^m ↔ the order on exponents. Nesting direction
(ordinal arithmetic, fixed here because it is easy to get backwards): for
exponents x + y (x the CNF-leading part), ω^(x+y) = ω^x · ω^y = "ω^y copies of
ω^x" — the TRAILING (smaller) CNF terms give the OUTERMOST/most-significant
nesting. E.g. ω^(ω+1) = ω copies of ω^ω. So a type-ω^δ set, δ = (d_1,…,d_m),
peels from its LAST CNF term: an ω^(d_m)-indexed nest of ω^(ω·d_{m-1} + …)-blocks,
and so on inward, the innermost scale being the ω^(ω^(m-1))-power structure.

## B. Full-type characterization: STACKED TOWERS (PROVED)

Call a **tower of shape (δ_1 < … < δ_M)** (δ_j ∈ N^m, lex increasing) a subset
T = T_1 ∪ … ∪ T_M of ω^(ω^m) with T_1 < T_2 < … < T_M (as sets, in the ordinal
order) and otp(T_j) = ω^(δ_j). A **full tower** has M = ω and {δ_j} lex-cofinal
in N^m.

**Lemma B1.** otp(Σ_j ω^(δ_j)) over an increasing ω-sequence δ_1 < δ_2 < …
lex-cofinal in N^m equals ω^(ω^m).
Proof. The ordered sum Σ_j ω^(δ_j) is ≥ ω^(δ_j) for every j, hence ≥
sup_j ω^(δ_j)·1; since the δ_j are cofinal, sup_j ω^(δ_j) = ω^(ω^m). For ≤:
Σ_{j≤J} ω^(δ_j) ≤ ω^(δ_J)·J < ω^(ω^m), and the total sum is the sup of its
partial sums, each < ω^(ω^m), with the sum a countable increasing union, so
≤ ω^(ω^m). (Both bounds use that ω^(ω^m) is additively indecomposable.) ∎

**Lemma B2 (the tower characterization).** X ⊆ ω^(ω^m) has order type ω^(ω^m)
⟺ X contains a full tower.
Proof. (⟸) By B1 the tower has type ω^(ω^m); X ⊇ tower and X ⊆ ω^(ω^m) gives
equality. (⟹) ω^(ω^m) is additively indecomposable, so every final segment of X
has full type ω^(ω^m), in particular type ≥ ω^δ for every δ < ω^m. Build the
tower recursively: given δ_1 < … < δ_j chosen and T_1 < … < T_j built inside X,
pick δ_{j+1} > δ_j with {δ_i} marching cofinally (e.g. δ_j = (j,j,…,j)); the final
segment of X above sup T_j has type ω^(ω^m) > ω^(δ_{j+1}), and any set of type
≥ ω^(δ_{j+1}) contains a subset of EXACTLY type ω^(δ_{j+1}) (initial segment). ∎

**Lemma B3 (recursive binary miniatures — the shape grammar).** Define, by
recursion on the exponent ε ≤ ω^m (CNF-peeling the LAST/smallest term) and a
uniform tower-length parameter M:
* BinM(ω^0) = a single element;
* ε = η + ω^0 (last term finite part): BinM(ω^ε) = two lex-ascending
  BinM(ω^η)-blocks, all elements of the two blocks agreeing above the single
  cross-split position (every cross pair then splits at that same position —
  forced, see below);
* ε = η + ω^e with e ≥ 1 (last term a LIMIT power): BinM(ω^ε) = an M-TOWER of
  BinM(ω^η)-blocks stacked along shapes marching lex-cofinally in ω^e
  (recursively: the index structure is itself BinM at the lower rank).
Then (i) every set of type ω^ε contains a BinM(ω^ε) for every M; (ii) at m=1
this grammar yields: pairs at e=0 and stacked binary grids of heights 1..M at
the top — split-position hierarchies: a binary height-h grid is exactly 2^h
elements whose internal split positions lie strictly below their cross-split
positions level by level (cross pairs at a level automatically split at one
common position, since both sides agree above it).
Proof. (i) Induction along the CNF peel using two facts: a set of type σ·τ is an
ordered τ-indexed sum of type-σ blocks with cofinally many full blocks (the
ordered-sum argument of THM-453 Lemma C with σ in place of ω^(n-1)); and for
τ = ω^(ω^e) limit, Lemma B2 at rank e gives the tower form of the index
structure. Pruning each ω-branching to its first two usable blocks and each
tower to length M gives the finite object. (ii) is the e=0 chain of the
grammar; the split-position characterization: elements of distinct top blocks
agree above the top split (construction), hence every cross pair splits exactly
there; internal splits are below by recursion. ∎

**Definition (binary M-tower at rank m).** BT(m, M) := BinM(ω^(ω^m)) — by the
grammar, an M-tower whose parts are the recursive binary shapes. At m=1: parts
are binary grids of heights 1,…,M (with free tails), stacked lex-ascending.
By B2+B3(i), every full-type X ⊆ ω^(ω^m) contains a BT(m, M) for every M.

## C. The König bridge (PROVED, same skeleton as THM-453 D1)

Finite ambient: A(m; s, c) := functions [s]^m → [c] (a truncation of the
function presentation, order = largest disagreement; a finite linear order).
**Q^(m)(M; s, c):** there is a triangle-free graph on A(m; s, c) such that every
binary M-tower inside A(m; s, c) contains an edge.

Antitone in (M ≤, s, c): a bigger-parameter witness restricts (hitting edges lie
inside their towers; towers of the smaller ambient are towers of the larger;
binary (M+1)-towers contain binary M-towers... NOTE the M-monotonicity runs the
OTHER way — killing all M-towers is HARDER for small M; we track the SAT region
in the full parameter cube and only use restriction in (s, c)).

**Theorem C1.** If for some fixed M, Q^(m)(M; s, c) holds for all s, c, then
ω^(ω^m) ↛ (ω^(ω^m), 3)², witnessed by a strong witness (no independent binary
M-tower).
Proof. König on the (s, c)-directed system of witnesses (restriction-coherent as
in THM-453 D1) gives a triangle-free G on all of ω^(ω^m) (= ∪ A(m; s, c)) with
no independent binary M-tower anywhere; by B2/B3 a full-type independent set
would contain one. ∎

**Corollary C2.** ω^(ω^m) → (ω^(ω^m), 3)² (TRUE for m = 1 by Chang, m = 2 by
Schipperus) ⟹ for every M there are finite (s, c) where Q^(m)(M; s, c) fails.
**Chang's theorem forces a family of finite tower-Ramsey cutoffs at m = 1** —
to our knowledge never computed; same at m = 2. At m = 4, Schipperus's negative
makes persistence possible (if his witness is strong); **at m = 3 the behaviour
of Q^(3) for growing (s,c) is a computational probe of the open $1000 case**:
persistence + a uniform rule would PROVE ω^(ω³) ↛ (ω^(ω³),3)² via C1.

## D. Computations (exact; structural tower verifier cross-validated against brute
force, 0 disagreements; `erdos592_chang_towers{,_v2}_macmini_s2.py`)

m=1, M=2 tower game on ambient [C]^s (kill all binary 2-towers = independent pair
below an independent binary 2-grid). Equivalent form used by the fast verifier:
since models are triangle-free, ANY 3 below-candidates contain an independent
pair, so a witness must leave every independent 2-grid with ≤ 2 surviving
elements below it.

| (s,C) | vertices | result |
|-------|----------|--------|
| (2,3) | 9 | SAT, 5-edge witness |
| (3,2) | 8 | SAT, 5-edge witness |
| (2,4) | 16 | SAT, 34-edge witness |
| (3,3) | 27 | grinding at session close (tens of thousands of CEGAR rounds — UNSAT-convergence-shaped; result lands in `erdos592_chang_towers_v2_macmini_s2.out`) |

* Chang's theorem GUARANTEES the sweep hits UNSAT (C2); the first UNSAT (s,C) is
  the first computed "Chang number".
* Minimal grid-witness sizes for the THM-453 game (new sequence, certified
  optimal at t ≤ 3): **f_grid(2)=1, f_grid(3)=7, f_grid(4)≤30** (35-edge and
  30-edge witnesses known; optimality run timed out below 30).
* m ≥ 2 instances require the general-shape recursive enumerator (the B3
  grammar); not yet implemented — the m=3 instance is the open case's probe and
  the designated next-session centerpiece.

## Honesty

C1/C2's "strong witness" caveat is inherited from THM-453: positive theorems force
the finite cutoffs (C2 unconditional), but known NEGATIVE theorems only predict
persistence if the classical witnesses kill binary towers — not just full-type
sets. The m=3 probe is evidence-generating, not (yet) theorem-generating in the
negative direction; in the positive direction (if m=3 cutoffs appear and a
pattern emerges) it generates the right finite conjectures to attack the open
case from the Chang/Schipperus side.
