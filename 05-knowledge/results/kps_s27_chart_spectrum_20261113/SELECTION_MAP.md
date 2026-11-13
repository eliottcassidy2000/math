# Where target selection stands, and what each construction points at

kind-pasteur-2026-11-13-S27. This is a reading map built from three read-only sweeps of the repository (2026-11-13).
Session 313 works the same owner prompt in parallel; its packet
([session313_zeckendorf_target_selection_20261113](../session313_zeckendorf_target_selection_20261113/README.md)) holds
the Zeckendorf selector, the independent physical switches and the mod-3 bypass, and is not repeated here.

Tags: **P** proved, **FE** finite-exact with range, **C** conjecture, **S** speculation. CP is Chan–Pak, which remains
**OPEN**.

## 1. The reduction everyone works under

- **The prime/7p reduction** is THM-4094 Theorem 6.1, not THM-5108 (P). CP holds iff every odd prime other than 7
  and every 7p with p not 3 is a strong count. THM-5108 closes only the smooth sector 3^a 5^b 7^c 31^d 41^e.
- **One strong component.** A prime count has exactly one nontrivial strong component (Session 309 HYPERGRAPH_SEAMS,
  P). Every selector must therefore build one new strong core per prime target.

## 2. Near-complete families and what completeness would need

| family | status | the missing all-order statement |
|---|---|---|
| THM-4875 (D) F2: one-apex word v joined to a bank of rooted triples of order <= 7, H = g(v) h − e(v) e + (1+v) R | FE: every odd N <= 2^32 except 7, 21 (audited). 2^35 in one unreplicated run | a fixed bank; HYP-9480 (C) predicts sporadic misses of size exp(−capacity) |
| O/Z/J grammar, uniquely Hamiltonian (HYP-9494 (D)) | FE to 10^9: 29 misses, largest 2647. With two root changes, 140 misses, largest 3,780,967 | an effective two-root-change threshold, or the four-block lemma (session199 hidden-structure) |
| 63 block matrices with runs <= 5 (THM-5148 (E), HYP-9677) | FE: none missed in (10391, 2^19] | first coordinates controlled mod 2^k and in short intervals |
| 2-adic lines 1 + c 2^s and anti-seed lines c 2^s − 1 | P for c <= 129 (THM-5126); FE seeds and anti-seeds for every c <= 4097 (THM-5130) | a seed and an anti-seed for every c with a uniform threshold (HYP-9658, HYP-9662). Any N ≡ 1 mod 4 lies on its own line s = v_2(N−1), and N ≡ 3 mod 4 on s = v_2(N+1), so the target's low bits select the line with no search |
| types (c, r) (THM-5145, HYP-9675) | FE: 92.8 percent of small pairs at order 10 | HYP-9675 in quantitative form, which is equivalent to CP plus (R) |
| one-arc closure of a point (HYP-9656) | FE: N(4095) = 145,433 | N(r0) -> infinity |

## 3. Selector mechanisms and the obstruction each met

- **Fixed component banks.** Session 309 FIXED_BANKS uses the ineffective S-unit theorem; prime_terminal T1 gives
  O((1 + log X)^D) counts for a fixed frame and bank.
- **Determinant repairs.** Session 310 shows every e = ±1 output is ±1 mod 96, an infinite residue obstruction.
- **Neutral refills.** THM-5151 (9) shows a clock fed only by neutral refills stops after cap(Q_N0) transfers.
- **Fixed word languages.** symmetry_port_return gives at most C (Y/K_0)^(1/4) outputs below Y.
- **Absorption.** prime_terminal absorption needs a repair radius of order X/(log X)^(D+1) to cover every prime.
- **Collatz descent.** For every depth, infinitely many primes have only composite bounded descendants (Session 309
  PRIME_COLLISION).
- **Anchors.** THM-5114 gives a multiplicatively closed countermodel with every even gap. THM-5115 leaves exact
  anchor selection OPEN.
- **The 3-adic chart itself.** THM-5152 (this packet) shows the chart class omits 31, 35, 53, 113, 127, 149, 173,
  181, 203, 217 at every order. A chart-only selector cannot be complete.

The common pointer: every route needs **new strong cores per target**, built by budget-raising (non-neutral) moves
through composite intermediate states, with a bank or repair radius that grows with the target.

## 4. The Zeckendorf tricolour and why it has not yet transferred

The owner's construction (2026-09-26, crossroads_crossing_20260926_colour.md, collatz_crossings_20260926_potential_and_seeds.md):

- Write 1, 2, 3, ... row by row in a triangle; row r holds labels r(r−1)/2 + 1, ..., r(r+1)/2.
- Colour each cell by its diagonal d = r − c + 1, according to the lowest digit of d's Zeckendorf representation:
  - red if it is F_2 = 1;
  - black if its index is odd and at least 3;
  - blue if its index is even and at least 4.
- This is the Wythoff AA/B/AB split.
- The extension with three typed units 1_K, 1_B, 1_R gives every n exactly 2 or 3 representations
  (reset_20260926_colours.md, P).
- The zero-slot alphabet gives 3 or 4 (odd_cell_gluing_20261023/ZERO-SLOT.md, P).

There is no proved link to path counts. The one direct attempt, inverse65_zeckendorf_20261103, proves that a
compiler which reads Zeckendorf digits, only appends, and is correct after every prefix is impossible (Z2). A
finite-state compiler that is right only at the end needs a unimodular emitted cycle with characteristic
polynomial (X − 1)(X^2 − L_q X + (−1)^q).

The sweep also checked THM-5148's O/Z/J blocks. Every block with a nonempty word has determinant 0 or 2^n, and the
pure root change has characteristic polynomial t^3 − t^2 − t − 1 (tribonacci). So the gadgets' native numeration is
binary or tribonacci, not golden. Session 313 builds its Zeckendorf decoder from typed units with this in mind.

**What the tricolour still offers (S).** Its multiplicity (2-4 representations per n, with typed units that look
the same but carry different hidden charges) is the feature target selection lacks. Each target should have a few
representations, and a carry (merge) step should change hidden data without changing the count. Session 276 has a
pure +2 step; the missing piece is a count-neutral merge at every scale with bounded defect.

## 5. A route not yet tried here: escalation and finite truth sets

The sweep found no Bhargava-style escalation and no finite-truth-set argument in the repository. The analogy (S):

- the 15/290 theorems decide universality of a quadratic form from a finite truth set;
- the reason is that escalator lattices, built by adjoining a vector for the smallest missed value, run out after
  finitely many steps;
- the tournament analogue escalates a selector family by adjoining a smallest-order tournament for the smallest
  missed odd value (the truant), with fixed composition rules, for example THM-4875's one-arc join F2.

Two facts bound what this can give.

- **The Poisson heuristic** (HYP-9480 (C)) predicts that every fixed finite bank eventually misses. Escalation should
  therefore never stop, and the useful quantity is how slowly it grows (number of escalators against X).
- **Partial finite-truth phenomena already exist.** The two-step family's six all-order omissions are decided on
  [3, 225) (THM-5129 (C)), and this packet's chart-class omissions are decided below Busch's bound. In each case a
  minimum-H bound (Busch, THM-4744) converts finite data into all-order statements for prime or unfactorable values.

That last mechanism is the closest thing to a truth set this repository has. Busch's bound n <= 1 + 3 log_5 H for
strong tournaments makes every statement of the form "the class C, closed under deleting singleton components,
contains/omits the prime p" decidable by a finite search, provided C is hereditary.

## 6. The single-lemma reductions closest to a proof

1. **Ternary gate.** Session193 prime_induction and THM-4094 §6: for attained h >= 3, h not 7 or 21, d in {2, 4}, if
   3h + d is prime or 7p then 3h + d is attained. The parent is forced by N mod 6, so selection is automatic.
   prime_terminal calls this "exactly the source obligation in Theorem B".
2. **Bounded root changes** (HYP-9494): two root changes suffice above 3,780,967, which is FE to 10^9.
3. **A finitely generated matrix orbit** (HYP-9677, HYP-9640): 63 or 56 explicit nonnegative matrices; FE to 2^19
   and 2 · 10^5 respectively.
4. **Two-step line supply** (HYP-9661): r(N) >= 1 suffices; FE to 382,687.

## 7. Bypass ideas on record

| idea | state |
|---|---|
| circle method over k-apex families (s90 README) | k >= 6 apexes pass the size test; the multi-apex saving and Mersenne-type major arcs are missing |
| local limit law (HYP-9659 (A)) | implies CP by THM-5127; unproved |
| thin-orbit density (Bourgain–Kontorovich) | gives density one only; Haag–Kertzer–Rickards–Stange warn that local-to-global principles can fail |
| automata (Session 297) | needs two automaticity theorems for the same set |
| sumsets of one-apex values | impossible by THM-5124 (2) |
| residue or covering arguments | repeatedly proved insufficient (THM-5114) |
