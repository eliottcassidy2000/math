# Build required_edits.json and validate every old string against the current worktree files (read-only).
import json, sys, os
ROOT = "/Users/e/Documents/GitHub/math-wt-pyth-cuboids"
T = "01-canon/theorems/THM-4811-arc-strong-connectivity-is-the-exact-complement-distance-to-the-nonstrong-locus-black-refinement-and-the-dipole-law.md"
P = "04-computation/experiments/black_terminal_distance_20261027/classification/PROOF.md"
E = []
def add(f, old, new, reason):
    E.append({"file": f, "old": old, "new": new, "reason": reason})

# ---- Finding 1 (MAJOR): extremal-predecessor overclaim ----
add(T,
"maximum-score one), by blue-cut cosets, the prime-orbit block normal form and two vertex arguments; this also proves the\n"
"open-path form of the extremal-predecessor statement left open in black_descent_20261020c/cycles.md section 7. If",
"maximum-score one), by blue-cut cosets, the prime-orbit block normal form and two vertex arguments. The same case\n"
"analysis does not use `lambda = delta^0` and shows that every nonregular strong self-converse tournament has a black\n"
"Hamiltonian path starting at a minimum-score vertex; this is weaker than the extremal-predecessor statement of\n"
"black_descent_20261020c/cycles.md section 7 (closed form, or its open-path form: a Hamiltonian path from a minimum-score\n"
"vertex to a vertex of less than maximum score), which stays OPEN outside the special cases listed in\n"
"classification/PROOF.md (Remark after Theorem C). If",
"Finding 1 (major, status): Theorem C's case analysis proves only C' (a black Hamiltonian path from a minimum-score vertex). "
"A black path from m may end in M, so this does not give a path from m to a non-maximum vertex, i.e. the open-path form; "
"Steps 1 and 3 use blueness, not only endpoint scores.")
add(T,
"(HYP-9348 (b)); the closed extremal-predecessor form of black_descent_20261020c/cycles.md section 7 in the case of two\n"
"  minimum-score vertices.",
"(HYP-9348 (b)); the extremal-predecessor statement of black_descent_20261020c/cycles.md section 7, closed and open-path\n"
"  forms (Theorem C gives only the weaker black Hamiltonian path from a minimum-score vertex; the closed form is settled only\n"
"  when a minimum-score vertex has a non-strong deletion or the minimum score is 1 and attained once, the open-path form\n"
"  also when the minimum score is attained once; FINITE-EXACT through order 10: both forms hold for every nonregular strong\n"
"  self-converse class, independent classification audit).",
"Finding 1: the OPEN entry implies that the closed form is proved whenever the number of minimum-score vertices is not 2. "
"In fact it is proved only in case (c1) (F10) and when mu = 1 and delta = 1; for mu = 1 and delta >= 2, and for every "
"mu >= 2, the argument starts from the blueness hypothesis (H), not from the negation of the statement.")
add(T,
"(E)); the closed extremal-predecessor statement when two vertices have minimum score.",
"(E)); the extremal-predecessor statement of black_descent_20261020c/cycles.md section 7 (closed and open-path forms)\n"
"outside the special cases listed in classification/PROOF.md (Remark after Theorem C).",
"Finding 1: same correction in section (F).")
add(P,
"Remark: this is the open-path form of the \"extremal predecessor\" statement left OPEN in cycles.md Section 7; the\n"
"closed form is not needed and is not claimed (Step 3 is non-constructive).",
"Remark (what (c1)-(c2) prove without lambda = delta). The case analysis never uses lambda = delta (only the last\n"
"sentence does), so it proves: **every nonregular strong self-converse T (n >= 4) has a black Hamiltonian path starting at a\n"
"minimum-score vertex** (a black path ending in M becomes one starting in m by sigma and reversal). This is weaker than\n"
"the \"extremal predecessor\" statement left OPEN in cycles.md Section 7, both in its closed form and in its open-path form\n"
"(a Hamiltonian path from a minimum-score vertex to a vertex of less than maximum score): a black path from m may end in M,\n"
"and Steps 1 and 3 use the blueness in (H) (the coset lemma, F2 at inner positions), not only \"every path from m ends in\n"
"M\". The argument settles the open-path form only in case (c1) and when mu = 1 (Step 2 uses (H) only through \"every\n"
"Hamiltonian path from v ends in M\"), and the closed form only in case (c1) and when mu = 1 and delta = 1 (then every\n"
"Hamiltonian path from v is closed). Both forms remain OPEN; an independent check finds the closed form true for every\n"
"nonregular strong self-converse class of order <= 10.",
"Finding 1: the Remark claimed Theorem C is the open-path form of the extremal-predecessor statement; it is not.")
add(P,
"- Theorem C, Step 3 is by contradiction: it shows a singleton-type black candidate exists but does not exhibit it,\n"
"  so the closed \"extremal predecessor\" statement of cycles.md Section 7 stays open in the case mu = 2.",
"- Theorem C (c2) is by contradiction from (H); Steps 1 and 3 use the blueness of the paths, not only their endpoint\n"
"  scores. It shows that a singleton-type black candidate exists but does not exhibit it, and it does not prove the\n"
"  \"extremal predecessor\" statement of cycles.md Section 7 (closed or open-path form) except in the cases listed in the\n"
"  Remark after Theorem C; both forms stay OPEN, not only for mu = 2.",
"Finding 1: the honest note implied the closed form is open only for mu = 2.")

# ---- Finding 2 (minor): Benhocine-Wojda is CITED ----
add(T,
"  (C) Theorems A and B, all orders, given the PROVED repo inputs cited in classification/PROOF.md (THM-2477 via\n"
"  Benhocine-Wojda, LEM-003,",
"  (C) Theorems A, B and C, all orders, given the CITED Benhocine-Wojda theorem (the n >= 7 step of THM-2477's all-order\n"
"  rigidity, taken there from the publisher abstract) and the PROVED repo inputs cited in classification/PROOF.md (THM-2477,\n"
"  LEM-003,",
"Finding 2 (minor, status): F6 rests on a CITED external theorem. merged_metagraph_finite_pure_atoms_20261020 types it CITED "
"and gives it at abstract level, so it is not a PROVED repo input. Theorem C is also proved at all orders, and (D1) uses it.")
add(T,
"  CITED / CREDITED: Camion 1959; THM-4729 (lines, blue/black),",
"  CITED / CREDITED: Camion 1959; Benhocine-Wojda 1983 (Hamiltonian bypasses in tournaments of order >= 7, through THM-2477\n"
"  and merged_metagraph_finite_pure_atoms_20261020 section 3); THM-4729 (lines, blue/black),",
"Finding 2: list the cited external theorem.")
add(P,
"Camion. All are marked PROVED in the repo and were read for this note.",
"Camion. All are marked PROVED in the repo and were read for this note; the n >= 7 step of THM-2477's all-order\n"
"  rigidity is the CITED Benhocine-Wojda theorem (publisher abstract), so Theorems A-C are PROVED given that citation.",
"Finding 2.")

# ---- Finding 3 (minor): which classes the censuses test ----
add(T,
"At order 11 the census of all 884,464,590 strong classes\n"
"(classification/census_sc.cpp) finds exactly one exception.",
"At order 11 the census of all 884,464,590 strong classes\n"
"(classification/census_sc.cpp; the descent test runs directly on the 276,971 nonregular self-converse classes, the others\n"
"being covered by Part I of the proof) finds exactly one exception, `E_2(S5)`.",
"Finding 3 (minor, wording): census_sc tests only the nonregular self-converse classes; the rest follow from Part I (a proof).")
add(T,
"  (classification/): the order-11 census of all 884,464,590 strong classes (exactly one exception, E_2(S5)), bottom-layer",
"  (classification/): the order-11 census of all 884,464,590 strong classes (descent test on the 276,971 nonregular\n"
"  self-converse ones, the rest by Part I; exactly one exception, E_2(S5)), bottom-layer",
"Finding 3.")
add(P,
"- **FINITE-EXACT** controls: complete gentourng censuses n <= 11, exhaustive structured searches n <= 16",
"- **FINITE-EXACT** controls: complete gentourng censuses n <= 11 (the descent test runs on every strong class for n <= 9 and\n"
"  on the nonregular self-converse classes for n = 10, 11, the others being covered by Part I), exhaustive structured\n"
"  searches n <= 16",
"Finding 3.")

# ---- Finding 4 (minor): 'stays in that regime' ----
add(T,
"For Theorem B, a singleton-type descent keeps `lambda = delta^0`, a descent from\n"
"`2 <= lambda < delta^0` stays in that regime and cannot enter `E`,",
"For Theorem B, a singleton-type descent keeps `lambda = delta^0`, a descent from\n"
"`2 <= lambda < delta^0` stays in the regime `lambda < delta^0` and cannot enter `E` (for `lambda = 2` at most one vertex of\n"
"the target has outdegree below 3, while members of `E` have three vertices of outdegree 2) or `R5`,",
"Finding 4 (minor, wording): from lambda = 2 the target has lambda = 1, outside '2 <= lambda < delta^0'.")

# ---- Finding 5 (minor): title, n >= 4 ----
add(T,
"the black distance to the non-strong locus is lambda(T) for every strong T other than R5 (no black move) and the members of E (black distance 2 = lambda + 1).",
"the black distance to the non-strong locus is lambda(T) for every strong T on n >= 4 vertices other than R5 (no black move) and the members of E (black distance 2 = lambda + 1).",
"Finding 5 (minor, wording): C3 (n = 3) is strong and pure blue, so it has no black move either; Theorem B is stated for n >= 4.")

# ---- Finding 6 (minor): P_t ----
add(T,
"`(P_0, c_1, ..., c_(t-1), P_t)` with `P_0` a Hamiltonian path of `S` into `x0`, and the stabiliser of `x0` (independent",
"`(P_0, c_1, ..., c_(t-1), P_t)` with `P_0` a Hamiltonian path of `S` into `x0` and `P_t` one of `sigma(S)` from\n"
"`sigma(x0)`, and the stabiliser of `x0` (independent",
"Finding 6 (minor, completeness).")

# ---- Finding 7 (minor): audit status ----
add(T,
"  The classification proof was written in this session by a dedicated proof agent and reviewed by a second agent with its\n"
"  own code (no logical error found); not yet independently audited at the repository level.",
"  The classification proof was written in this session by a dedicated proof agent, reviewed by a second agent with its\n"
"  own code, and independently audited (2026-10-27, PASS WITH CORRECTIONS: every step of Parts I-VIII and Theorems A, B, C\n"
"  confirmed; an overstated side claim about the extremal-predecessor statement withdrawn).",
"Finding 7 (minor, status): record this audit.")
add(T,
"E(M; Z) properties, with the census through order 10 reproduced; Theorems A, B, C are under a separate independent audit.",
"E(M; Z) properties, with the census through order 10 reproduced; Theorems A, B, C audited separately (PASS WITH\n"
"  CORRECTIONS; own code: Theorems A and C on every strong class through order 9 and every nonregular self-converse class of\n"
"  order 10, the family E from its definition through order 18, the last-block lemma for blocks of order <= 9, Lemma BL for\n"
"  |S| <= 9, exhaustive lambda = 1 < delta^0 searches through order 17 and lambda = 2 < delta^0 searches at orders 14-17, and\n"
"  exact class-level black distances through order 9).",
"Finding 7.")

# ---- Finding 8 (minor, optional): exhaustiveness at n = 15, record the audit's searches ----
add(P,
"| `family3.cpp` (n = 15) | lambda = 2 < delta with a fixed middle vertex (t = 2) | 1131 configurations, 0 failures |",
"| `family3.cpp` (n = 15) | lambda = 2 < delta with a fixed middle vertex (t = 2); exhaustive in that regime, since L_0 and L_t have 7 vertices each and the single middle vertex is sigma-fixed | 1131 configurations, 0 failures |",
"Finding 8 (minor, optional): the n = 15 search is exhaustive for the regime 2 <= lambda < delta.")
add(P,
"| `audit/` (independent auditor, own code, no `tour.hpp`) | last-block lemma for all strong R with m <= 9 and all A; Lemma BL hypotheses for a <= 9; Theorem C on 78,387 random self-converse nonregular strong T with lambda = delta, n = 10..13; the 19 members of E with n <= 15 | 0 counterexamples; BL pairs 1, 2, 5, 12, 28; all 78,387 have a singleton-type black descent; 19/19 |",
"| `audit/` (independent auditor, own code, no `tour.hpp`) | last-block lemma for all strong R with m <= 9 and all A; Lemma BL hypotheses for a <= 9; Theorem C on 78,387 random self-converse nonregular strong T with lambda = delta, n = 10..13; the 19 members of E with n <= 15 | 0 counterexamples; BL pairs 1, 2, 5, 12, 28; all 78,387 have a singleton-type black descent; 19/19 |\n"
"| second independent audit (2026-10-27, own C + nauty code) | Theorems A and C on all strong classes n <= 9 and all nonregular self-converse classes n = 10; E from its definition (112 members, n <= 18); last-block lemma m <= 9 with (H1) tested directly; Lemma BL a <= 9; exhaustive lambda = 1 < delta searches n <= 17 and lambda = 2 < delta searches n = 14..17; all F/G shapes with a <= 8; exact class-level black BFS n <= 9 | no counterexample; every exception found is R5 or a member of E; d_black = lambda for every strong class n <= 9 except R5 |",
"Finding 8 (optional): record the extended searches.")

# validate
ok = True
for i, e in enumerate(E, 1):
    txt = open(os.path.join(ROOT, e["file"]), encoding="utf-8").read()
    c = txt.count(e["old"])
    if c != 1:
        ok = False
        print("EDIT", i, "old occurs", c, "times in", e["file"])
print("edits:", len(E), "all unique:", ok)
out = sys.argv[1]
json.dump(E, open(out, "w"), indent=2, ensure_ascii=False)
