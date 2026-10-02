---
id: THM-4524
title: "Selfie tournaments. The N optional loops are the switching gauge, and the N-1 base-path arcs are their discrete derivative (an exact 2-to-1 map). In the loop spins, H has Walsh degree at most 2 floor(N/3), so the alternating loop sum of H vanishes. An arc on an odd number of Hamiltonian paths for EVERY arc needs N = 1, 2 (mod 4); this first happens at N = 6 (QR_7 minus a vertex, unique). Cayley tournaments of odd abelian groups have every arc count even. A049313 counts the Euler graphs whose automorphisms all reverse an even number of edges (answers OPEN-Q-060)"
status: >
  PROVED + INDEPENDENTLY AUDITED:
  - A1, the gauge theorem;
  - A2, the oriented-path class sums;
  - A3, the loop-Walsh expansion and the degree bound with its corollaries;
  - B1, the fixed-point OCF, and selfie Redei;
  - the arc-parity constraints;
  - the dead-arc formula;
  - no tournament on n >= 3 vertices has every arc on some HP and some arc
    on every HP;
  - C1, Cayley all-even;
  - the antipodal half of C2;
  - the shaving lemmas (rho = 2 exactly for all-even; beta <= Hall);
  - E1, OPEN-Q-060.
  FINITE-EXACT: the census of arc-HP counts for N <= 10 (N = 10 lane-only);
  beta = Hall for N <= 9 (N <= 7 audited); Redei shaving for N <= 7;
  QR_q minus a vertex is all-odd for q = 7, 11, 19, 23, 27 (q <= 23
  audited).
  REFUTED:
  - the owner's conjecture, in both natural readings;
  - "even N implies at least N - 1 odd arcs" (minimum 5 at N = 10).
  OPEN: HYP-9167 (Paley minus a vertex has every arc odd; all-odd
  tournaments exist iff N = 2 mod 4); HYP-9168 (beta = Hall bound); the
  bijective form of E1.
  CITED: THM-474; THM-479; THM-1467; Camion 1959; Alspach 1967; Brauer's
  permutation lemma; Babai-Cameron 2000; OEIS A049313, A002854, A038375.
  Notation: T a tournament on N vertices; H(T) = #Hamiltonian paths (HPs);
  c(e) = #HPs through the arc e; switch_L reverses every arc with exactly
  one end in L.
  (A1) The map (tiling t, loop set L) -> switch_L(T_t) is exactly 2-to-1
  onto the labeled tournaments, with fibres {L, V \ L}. The base-path arc
  (i+1 -> i) is reversed iff exactly one of i, i+1 is looped. So
  C(N-1,2) + N = C(N,2) + 1, and the +1 is the global complement.
  (A3) H(switch_L T) = sum over even A, |A| <= 2N/3, of h_T(A) times the
  product of sigma_v over v in A. Hence sum_L (-1)^|L| H(switch_L T) = 0
  for every T. The constant term is N!/2^(N-1), and the class sum is
  2 N!.
  (C) Parity constraints:
  - all arcs odd forces N = 1, 2 (mod 4);
  - all arcs even forces N odd;
  - #odd arcs = N - 1 (mod 2).
  No tournament with 3 <= N <= 5 or N = 9 has all arcs odd. N = 6 has
  exactly one class (QR_7 minus a vertex, 240 labeled tournaments).
  N = 10 has exactly two.
  (C1) For every Cayley tournament of an abelian group of odd order, every
  arc lies on an even number of HPs.
  (E1) A049313(n) = the number of unlabeled Euler graphs F on n vertices
  such that every automorphism g of F reverses an even number of edges of
  a fixed orientation. Equivalently: no automorphism has an odd number of
  even cycles whose antipodal pairs are edges of F.
source: collatz-procgen-20260922 session, selfie lane (2026-10-01), answering the owner's selfie-tournament / Hamiltonian-path-parity / shaving prompt of 2026-10-01; audited and promoted by the session orchestrator 2026-10-01
depends_on:
  - 01-canon/theorems/THM-474-tilings-are-switching-classes.md (the gauge: one tiling tournament per switching class)
  - 01-canon/theorems/THM-479-level-holonomy-and-branch-split.md (the affine S_n action on switching classes; the A049313 Burnside setup)
related:
  - 01-canon/theorems/THM-1467-oriented-spanning-tree-switching-sum.md (switching-class sum of H = N!)
  - 01-canon/theorems/THM-1470-even-tournaments-are-the-tournament-two-graph-theorem.md (the n = 1 mod 4 incarnation of A049313)
  - 01-canon/theorems/THM-062-forward-edge-distribution.md (deformed Eulerian numbers; their class sums are binomial)
  - 05-knowledge/hypotheses/HYP-9167-paley-minus-a-vertex-is-redei-rigid.md
  - 05-knowledge/hypotheses/HYP-9168-hp-blocking-number-equals-hall-bound.md
  - 00-navigation/OPEN-QUESTIONS.md (OPEN-Q-060)
note: 05-knowledge/results/procgen_selfie_20261001_selfie_tournaments.md
scripts: 04-computation/experiments/procgen_selfie_20261001_{run,lib}.py, procgen_selfie_20261001_{arcs,euler,paley2}.c
script_audit: 04-computation/experiments/procgen_selfie_20261001_orchestrator_check.{py,c}
output: 05-knowledge/results/procgen_selfie_20261001.out
output_sha256: 3529e8a0a03f19407779ed1e4258ea5ace0ce1b7c285003e3932bb5ab8ef204a
output_n10: 05-knowledge/results/procgen_selfie_20261001_n10_census.out (sha256 f0bfccd69c1d52b46510f3ca35821d1d81e02cfb4b2beea7e42c84a7a9f13db0), procgen_selfie_20261001_n10_oddfilter.out (sha256 dda760efd761e28ef452debcb894a00de445c8dd4bb6b14deb2b0eacff44464f)
output_audit: 05-knowledge/results/procgen_selfie_20261001_orchestrator_check.out
output_audit_sha256: 25c584f9dca62286c1b70fb14fbac8a4e3a22abe5bfc1f625e02e2656c913cb0
hash_basis: raw bytes
audit: >
  The orchestrator read the proofs of A1, A3(b) (the sign-reversing block
  reversal of the first odd segment), the parity constraints, the dead-arc
  formula, C1 (P -> reverse(phi(P)) with phi(x) = u + v - x), the
  antipodal half of C2 (orbit sizes q - 1 and (q - 1)/2 under F_q^*),
  rho = 2 (an odd end(v) forces an odd consecutive pair) and E1 (Brauer's
  lemma on the torsor F_2^E / cut space; the cocycle makes eps_F a
  homomorphism on Aut F). All were found sound.
  Independent code (the orchestrator's own C engine and Python; the lane's
  code was not read) confirms 58 checks plus E1 for every n <= 10:
  - all labeled tournaments N <= 7: the HP-covered counts 2, 40, 664,
    26048, 1934528; all-odd 0, 0, 0, 240, 0; all-even 2, 0, 184, 0, 96608;
    the min/max number of odd arcs; the dead-arc formula; "no tournament
    is HP-covered with an arc on every HP"; sum_L H = 2 N!;
    sum_L (-1)^|L| H = 0; the Walsh degree bound on every tiling, with top
    degree attained as stated; the two constant-H classes at N = 4; the
    x = 0 fixed-point identity; selfie Redei (N <= 6);
  - all classes N <= 9 (gentourng): the census table (N = 9: 191536
    classes, 189992 covered, 0 all-odd, 899 all-even, 167 strong with a
    dead arc);
  - for every class N <= 7: beta = Hall bound; Redei shaving to one HP
    fails only for the all-odd class; beta = N - 1 for every regular
    class;
  - the N = 7 counterexample to the sigma bound (beta = hall = 3 < 4);
  - QR_q minus a vertex all-odd for q = 7, 11, 19 (exact) and 23 (mod 2),
    with QR_q itself all-even;
  - every circulant on Z_3..Z_13 and every Cayley tournament of
    Z_3 x Z_3 all-even;
  - both all-odd N = 10 classes (H = 3929 and H = 15745) and all four
    5-odd-arc N = 10 witnesses;
  - E1: the per-permutation identity #Fix(g) = sum_F eps_F(g) with eps
    from the orientation definition (every permutation n <= 5, every cycle
    type n = 6, 7); orientation and cycle forms agree; Burnside = A049313
    for n <= 7; the direct orbit count of untwisted Euler graphs equals
    A049313 for n <= 6; cycle-form Burnside by F_2 linear algebra equals
    A049313 and A002854 for n = 8, 9, 10.
  The lane's runner was re-run (220.8 s, 133060 checks, ALL CHECKS
  PASSED). It is identical up to timing fields.
  Scope notes:
  - The exhaustive N = 10 census (59 min) and the q = 27 mod-2 run were
    not repeated. The lane's two N = 10 passes (exact and an independent
    mod-2 DP) agree.
  - beta = Hall for N = 8, 9 rests on the lane's branch and bound.
  - Literature (2026-10-01, THM-4531): E1 is the switching-class analogue of Royle, Praeger, Glasby, Freedman and
    Devillers, "Tournaments and even graphs are equinumerous" (J. Algebraic Combin. 57 (2023), arXiv:2204.01947). That
    paper proves von Broemssen's conjecture with the same sign homomorphism and Cauchy-Frobenius mechanism. No published
    statement of E1 itself was found. Both are the U = 0 and U = cuts instances of THM-4531's Theorem G.
---

# THM-4524 — selfie tournaments, arc parity of Hamiltonian paths, and the odd Mallows–Sloane partner

**PROVED + INDEPENDENTLY AUDITED.** Full note: [procgen_selfie_20261001_selfie_tournaments](../../05-knowledge/results/procgen_selfie_20261001_selfie_tournaments.md).

## 1. Loops versus the Hamiltonian path

The owner's selfie tournament adds `N` optional loops to a tournament. In the tiling model, the `N - 1` arcs of the base path `N -> N-1 -> ... -> 1` are fixed and the `C(N-1,2)` non-path arcs are the binary tiles.

The two kinds of extra structure fit together exactly:
- **The loops are the switching (gauge) bits.** Switching at the looped vertices sends a tiling tournament onto every labeled tournament, exactly twice (`L` and its complement).
- **The path arcs are the loops' discrete derivative.** The path arc `i+1 -> i` is reversed iff exactly one of `i, i+1` is looped.
- **The counts match up to one bit.** `C(N-1,2) + N = C(N,2) + 1`. The one extra bit is the global complement (`H^0` of the path), which acts trivially on tournaments.

So the path that "must exist" is the gauge-fixing condition, the loops that "may exist" are the gauge frame, and the tiles are the gauge-invariant data. Letting the loops vary turns the directed base path into an oriented Hamiltonian path of every type; each type occurs equally often (A2).

In the loop spins, `H` is a polynomial of degree at most `2 floor(N/3)`. So the alternating sum `sum_L (-1)^|L| H(switch_L T)` is 0 for every `T`, and for `N <= 5`, `H` on a switching class is an Ising Hamiltonian. Loops can also be read as odd 1-cycles in the odd-cycle formula (B1). The resulting "selfie Rédei" count is odd and is `H + 2|L| (mod 4)`.

## 2. The owner's 5/6 threshold

The conjecture was that every arc lies on a Hamiltonian path up to `N = 5`, and that from `N = 6` on every class has arcs on none. It is **false in both readings.**

**(i) Arcs on no HP.** HP-covered tournaments exist at every `N >= 3`: regular ones by Alspach's theorem, and for even `N` a regular one followed by a sink. Their share rises to 1 (0.25, 0.63, 0.65, 0.80, 0.92, 0.97, 0.99, 0.998 for `N = 3..10`). The dead arcs are exactly those that skip a strong component, plus the dead arcs inside the components.

**(ii) Parity.** Here the threshold is real but runs the other way:
- **No `N` from 3 to 5 has a tournament in which every arc lies on an odd number of HPs.**
- `N = 6` has exactly one such class, Paley `QR_7` minus a vertex (`H = 45`, the maximum).
- `N = 10` has exactly two such classes. `N = 9` has none. Parity forbids `N = 0, 3 (mod 4)`.
- Cayley tournaments of odd abelian groups (circulants, Paley) sit at the opposite extreme: every arc count is even.
- HYP-9167 conjectures that Paley minus a vertex is always all-odd, and that all-odd tournaments exist exactly when `N = 2 (mod 4)`.

## 3. Shaving

- **Single arcs.** Deleting `e` changes `H` by `c(e)`. So a single arc can be shaved keeping `H` odd iff `c(e)` is even, and the all-odd tournaments are exactly the ones that cannot be shaved at all.
- **Down to one HP.** For every other class with `N <= 7`, arcs can be shaved one at a time, with `H` odd at every step, down to a single Hamiltonian path.
- **Breaking parity.** All-even tournaments need exactly two deletions to make `H` even.
- **Killing every HP.** The fewest deletions that kill every HP equal a Hall-type deficiency bound in every class with `N <= 9` (HYP-9168). Regular tournaments survive any `N - 2` deletions.

## 4. OPEN-Q-060

The question was: what does A049313 (switching classes of tournaments) count, in the way that A002854 counts Euler graphs?

**Answer (E1).** A049313(n) counts the unlabeled Euler graphs whose automorphisms all reverse an even number of edges of a fixed orientation. This is Mallows–Sloane twisted by the cocycle of the tournament torsor:
- Labeled switching classes form an affine space over `F_2^E / (cut space)`.
- The dual of that space is the cycle space, i.e. the Euler graphs.
- Brauer's lemma, applied to the affine action, turns each fixed-point count into a sum of a sign character over the Euler graphs that `g` fixes.

The bijective form remains open.

**FORMALIZED 2026-10-01 (Lean 4.30 core, no Mathlib; package [`04-computation/lean/ProcgenSelfieEdim/`](../../04-computation/lean/ProcgenSelfieEdim/README.md); `verify.py` PASS, 503 theorems, axioms within {propext, Quot.sound}; orchestrator re-built from scratch and re-audited the axioms).**
- A1, the gauge theorem for all N: `path_arc`, `selfie_fibre`, `loops_ne_compl`, `selfie_onto`.
- A2: `switch_sum`.
- A3(e): `constant_switching_class`.
- C1 for cyclic groups: `circulant_arcCount_even`.
- The double counting and the parity constraints: `sum_arcCount`, `allOdd_mod_four_of_tournament`.
- Rédei's theorem itself, proved from scratch: `redei`.
- rho = 2: `parity_break_two`.
- Finite facts: QR_7 / QR_7 - v (`qr7_counts`, `qr7del_counts`); no all-odd tournament for N <= 5 (`no_allOdd_three_to_five`); `no_covered_universal_arc`.
- Not formalized: A3(a)-(d), B1, C2, E1, the dead-arc formula, the censuses.

**UPDATE 2026-10-01 (THM-4531): the bijective form of E1 is answered negatively.** No relabelling-invariant construction gives a bijection on isomorphism types between switching classes and even (untwisted) Euler graphs. This is PROVED for 5 <= n <= 100 and is impossible in the reverse direction for every n >= 3. For odd n, the canonical object on the tournament side is the unique member of each switching class with out-degree = in-degree (mod 4) at every vertex.

**Later work (2026-10-02; independent audit owed).** [`natural_matchings_switching_classes_20261002.md`](../../05-knowledge/results/natural_matchings_switching_classes_20261002.md)
was written in parallel with THM-4531 and extends its range: no natural bijection from switching classes to untwisted
Euler graphs exists for any 5 <= n <= 10^6 (PROVED, using rigid blocks of every odd order with no prime factor = 1
mod 8 and a computed additive condition; CONDITIONAL beyond 10^6). Its odd-n canonical member is THM-4531's P1, found
independently.
