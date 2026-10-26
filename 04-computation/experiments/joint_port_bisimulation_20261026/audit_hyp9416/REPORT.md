# Independent audit of HYP-9416 and `fresh_remarking_quotient_and_divisor_arrangements_20261026.md`

2026-10-26. Separate-agent audit, run in this directory. All code is the auditor's own and imports nothing from the
lane. The 57 edits in `required_edits.json` (built by `make_required_edits.py`) were applied by the lane session.

## Verdict

Every finite number in both documents reproduces exactly with independent code. The interpretation does not hold up:

- **Item (1) is mis-labelled.** It is filed as an OPEN hypothesis, but it is a one-line theorem for every n ≥ 3.
- **Item (2) misdescribes its mechanism.** The separation of the collision pairs is a one-step joint-port statistic
  built into V1, not a closure phenomenon.
- **Three sentences in the note are false.** The lane's own V1/sand column and new runs contradict them.
- **The link to C5 is wrong.** The construction is not the quotient that vitali-response C5 left open.
- **Smaller problems:** several misattributions and novelty overclaims, and one arithmetic wording error.

## What was run

- **`tcanon.py`: canonical labelling** by individualisation and refinement, written from scratch. Self-tests:
  - all labelled tournaments for n ≤ 6 give 1, 1, 2, 4, 12, 56 classes, with Σ n!/|Aut| = 2^C(n,2);
  - the rooted partition equals the brute-force one at n = 5 (176 classes);
  - n = 7 gives 456 classes, with Σ 7!/|Aut| = 2^21;
  - codes are invariant under relabelling at n = 9, 10, 12.
- **Check A, `burnside_counts.py`.** It uses no enumeration: Burnside over S_n, plus the strong-component identities
  T = 1/(1−A), V = T²B, R = T²Cc + 2T³B².
- **Checks B and C:** `quotients.py`, `collision_anatomy.py`, `single_context.py`, `extra_observers.py`,
  `sandwich_reconstruction.py`. Contents:
  - the V1 and V2 multiset quotients;
  - a normalized-distribution variant;
  - extra observers: the number of s→t paths, h*e only, the scalar H(h*e), and both one-sided classes;
  - each two-sided context tested alone.
- **Check D, `atoms.py`:** strong-atom histograms, coarser invariants, and a coefficient-by-coefficient comparison
  with `strong-atlas.json`.
- **Check E, `control_audit.py`:** the 81/27 control and the Δc3 flux table.
- **Check F, `arith_checks.py`:** the divisor arrangement dictionary and the arithmetic threads.

## Numbers

All lane values were confirmed except where marked.

**A.**
- Rooted two-port states, n = 3..7: 8, 32, 176, 1408, 17824. Of these, 2, 12, 92, 946, 14162 are strong.
- New: 376576 / 331032 at n = 8, and 13694400 / 12750980 at n = 9.

**B.**
- Confirmed table entries:
  - V1/free, V1/sand, V1/glue1 and V2/free;
  - the V1/V2 transition counts;
  - V2/sand = 8685.
- glue2 already separates every strong state at n = 3 and 4 (2 and 12), so the table's "–" entries were wrong.
- One-sided collisions: 7, 20, 70 pairs in 2, 12, 39 fibres. The lane omitted n = 3 (1 pair, 1 fibre) and n = 4 (0).
- How the pairs separate:
  - glue1 alone separates 0;
  - the V1/glue1 closure separates all;
  - glue2 separates all;
  - sand alone separates 6, 17, 64.
- New results:
  - The bare class observation `free`, closed under V1, already separates every pair. In every case this happens at
    the **first** refinement round.
  - The number of s→t Hamilton paths alone separates 1/1, 4/7, 16/20, 65/70.
  - V2 closure separates none of the pairs.
  - The normalized-distribution version of V1/free has 11031 strong cells instead of 11033 at n = 7 (11603 vs 11605
    cells overall).
  - With only the free class of h*e, V1 leaves 12371 of the 14162 strong cells; with (H(h*e), H(e*h)) it leaves 13543.
  - 3381 of the 14162 strong states have no s→t path.

**C.** glue1 depends only on the one-sided classes: an isomorphism of (C, t) extends by the identity to h*k, and one
of (C, s) to k*h. This was also verified on all states for n ≤ 7.

**D.**
- Confirmed counts at n = 3..7:
  - atoms: 1, 1, 6, 35, 353;
  - distinct joint targets: 1, 5, 44, 547, 9985, of which 0, 2, 32, 484, 9516 are strong;
  - atoms with blue presentations: 1, 1, 6, 9, 79;
  - presentations: 1, 5, 50, 903, 30773;
  - largest fibres: 1, 1, 2, 3, 8.
- F_A is pairwise distinct, also with the ports forgotten (vacuously at n = 3, 4).
- The atlas agrees on all 396 sources, coefficient by coefficient.
- Much less information already classifies the atoms:
  - the bare **set** of partner classes separates all 353 at n = 7;
  - the multiset of partner score sequences gives 351 classes.

**E.**
- Confirmed:
  - H = 462915 for both, with 67 labelled cyclic triples;
  - |Aut| = 81 / 27, and tc = 5715 / 17145;
  - the partner map is injective in both, with 74 partner classes shared;
  - T−'s three entrance phases carry 154305 each;
  - the phase families have 5715 classes each and are disjoint;
  - the vitali Δc3 flux table reproduces exactly.
- New:
  - All 74 shared classes lie in phase family 0.
  - T+'s blocks also carry 154305 each.
  - The joint (entrance, exit) phase mass depends only on the phase difference (42930, 45036, 66339).
  - Equal phase mass is therefore forced, because A− = A+^op has the same ordered path-cover counts.

**F.**
- Confirmed:
  - the region and vertex counts, and the weights s, σ, μ (847 squarefree n);
  - s(170) = s(182) = 154;
  - s⁻¹(42) = {30, 1681}, and this is complete;
  - the Gersonides contacts {1+1=2, 1+2=3, 1+3=4, 1+8=9}, with radicals 2, 6, 6, 6 and quality 1.2263;
  - the twin-prime double-hit set for v > 3 (v = 3 is an extra member, via 4 = 2²), with 6 the only self-hit.
- Goormaghtigh: 31 is the only value shared by two prime bases. This is complete for every value below 10¹², which
  means every prime base below 10⁶. With arbitrary bases up to 10¹⁰, the only coincidences are 31 and 8191.

## The lemma that proves (1) for every n ≥ 3

- In e*h*e, the new entrance vertex a has out-degree 1 and the new exit vertex b has in-degree 1.
- Because C is strong, every other vertex has out-degree and in-degree at least 2, with two possible exceptions:
  - s may have out-degree 1, when d⁺_C(s) = 1; its out-neighbour then has out-degree at least 2;
  - t may have in-degree 1; its in-neighbour has in-degree at least 2.
- Since a → s and t → b, the vertices a and b are determined by the unlabelled tournament.
- Then s = N⁺(a), t = N⁻(b), and C = e*h*e − {a, b}.

So e*h*e alone is a complete invariant of strong two-port states.

This was checked constructively. [C; s, t] was recovered from randomly relabelled sandwiches for all 6 / 12 / 120 /
1050 / 14826 labelled strong cases at n = 3..7. Each of the nine two-sided contexts is individually injective on
strong states for n ≤ 7. The context A2(s→t)·h·A2(s→t) is injective on **all** 17824 states; whether that holds in
every order is open.

## Verdict per claim

| Claim | Verdict |
|---|---|
| (1) glue2 separates all strong joint classes, n ≥ 5 | Numbers correct; holds from n = 3; PROVED trivially (lemma); "no compression" is a tautology for a complete invariant |
| (2) glue1 is one-sided and separates no collision pair | Correct; one-line proof |
| (2) V1 closure of glue1 separates all collision pairs | Correct but misleading: closure and glue1 play no role; one V1 step from `free` does it; the s→t path count does 65/70; V2 separates none |
| (3) F_A classifies atoms, also with ports forgotten | Correct (n ≤ 7); the histogram is not needed (the support set suffices); the data are the vitali atlas |
| (4) V1/free compresses; V2/free is the class partition | Correct; compression already at n = 4; the V2 part is trivial by construction |
| "This is the stochastic quotient left open by C5" | Inaccurate. Different carrier (two-port classes, not presentations X_n); counts, not distributions (the distribution quotient differs at n = 7); V1 is restricted to the current ports. Session70 also calls it a different operation menu |
| 81/27 numbers | All correct. The 74 shared classes are already in HYP-9414 (3), uncredited. The phase masses are forced (uncredited). "Cannot predict a single partner" and "the datum must keep …" are overstated or unproved |
| §5 dictionary | Correct; elementary and classical, not a new result |
| §6 arithmetic | Correct except the "6 is the first n" wording; novelty overclaims |

## Errors by severity

### Major

1. **Three false generalizations in note §1.**
   - "With fresh re-marking and any two-sided attachment observer, it is the joint two-port class": false. sand is a
     two-sided observer, and V1/sand has 5889 of 14162 strong cells.
   - "Only for observables that ignore attachment": false. sand compresses, h*e alone compresses (12371), and
     V2/glue1 stays at 14111.
   - "With only one-sided attachment observers, closure forces the joint class anyway": false for one-sided observers
     that see only one side, and false under V2.
   - The HYP's "for two-sided attachment observations there is no smaller quotient" has the same flaw.
2. **Item (1) is a trivial theorem presented as an OPEN hypothesis.** Its listed open question, "why order-2 contexts
   suffice", is answered: order 1 suffices.
3. **Item (2) and note §1 misdescribe the mechanism.** The phrases at fault are "closure turns one-sided observers
   into joint ones" and "a dynamic version of 289/253 … propagates through the dynamics".
   - V1 successor multisets are built from s→t paths, so they are joint statistics by definition.
   - The 289/253 pair itself has 4 vs 2 entrance-to-exit paths.
   - The HYP title says "closure under the dynamics" without saying V1.
4. **The identification with C5's stochastic quotient is wrong**; see the C5 row of the table above.

### Minor

5. **Thresholds and table entries.** The n ≥ 5 thresholds should be n ≥ 3 for (1) and n ≥ 4 for (4). The n = 3
   collision pair is omitted, and the table had "–" where 2 and 12 belong.
6. **Labels.** The trivially provable parts (glue1 one-sidedness, V2/free) are labelled OPEN or FINITE-EXACT.
7. **§4 misattributions.**
   - The masks 82468 / 84707 with H = 141 / 143 come from THM-168 and `collatz_effective_prefix_mass_20261005.md` §9.
     They are not in THM-169/170, and not in HISTORY.md, which uses 1354195 / 1351956 with H = 109 / 111.
   - The recovery-principle quote is from `vitali_atom_recovery_20261011.md`, not HISTORY.md.
   - "15/11" cannot be traced; vitali C2's control is 15 vs 23.
   - "The one formal object behind all of these" contradicts HISTORY.md's own caveat.
8. **A wording error in §6.** "6 is the first n whose square lies in a later doubling row" is wrong as worded: 3 is
   first, and 6 is the first even n. It also contradicts the general condition stated on the next line at a = 0. The
   cited "chains §2" should be §3.
9. **Novelty overclaims in §6.**
   - Goormaghtigh's 31 was marked "new here" but is already in HYP-9362 and the Leech note. Its classical status was
     missing:
     - Ratat–Goormaghtigh;
     - Bateman–Stemmler found 31 the only case below 1.275·10¹⁰ with x, y and N all prime;
     - Grantham (2024) found no new Goormaghtigh primes below 10⁷⁰⁰.
   - The radical-6 triples were marked a "new remark" but are already in THM-3833 (30), CORE-PAPERS-ABC-IUT (2) and
     abc_difference_families.
   - The Gersonides adjacencies {2, 4, 8} are in THM-4640.
10. **C4 is mis-scoped** by "beyond one run": its family is a single global Hamilton path through k blocks (Session70
    flags this too).
11. **Hygiene.**
    - `validate_289_253.out` has no producing script.
    - `aliquot_preimages.out` printed an unexplained failing check ("False", caused by v = 3).
    - bisim3's glue observation silently includes the free class.
    - "THM-1370" is ambiguous: the number is used by two files.
    - The HYP INDEX merges HYP-9416 under the hourglass entry's "independently audited" header.

A new MISTAKE entry is warranted for items 1, 3, 4 and 8.

## Sources for the Goormaghtigh literature status

- Goormaghtigh conjecture: https://en.wikipedia.org/wiki/Goormaghtigh_conjecture
- Grantham, "No new Goormaghtigh primes up to 10^700", arXiv:2410.03677
