# Critical expansions: the fifty-vertex tournament, the sixty-four hourglass twists, and how the repository's results fit together

klein, 2026-10-25, worktree `codex/session-tournament-twist-unify`. Reflection: idea provenance and a
proposed unifying reading. Truth lives in the linked canon and packets, not here.

## 1. What the owner's prompt turned out to be

The prompt named a 50-vertex tournament with one symmetry-respecting toggle, cancellation depth 5 over
`Q`, 21 modulo 5, 25 modulo 2, and separately an average over 64 quadratic twists equal to
`q^2 - 6q + 12 = (q-1)^2 - 4(q-1) + 7`. Reading it against the repository's own vocabulary:

- "cancellation depth" is Session62's order of vanishing at `t = -1` of the subset-orbit polynomial of the
  clone fibre (`session62_agent_theme_bridges_20261025/clone/README.md`), whose 24-vertex example realizes the
  regular `C3^2` action on nine gaps; "one independent symmetry-respecting toggle" is its invariant-toggle rank
  `r = 1` (a transitive gap action);
- the 50-vertex object is the next case: inflate a rigid Cayley tournament on `Z5^2` to `TT2` blocks, 25 gaps,
  regular `Z5^2` action — `Cay(Z5^2,S)[TT2]` with an explicit half-system `S`;
- "depth 21 modulo 5 / 25 modulo 2" are reductions of the integer polynomial, never computed before in the
  repository; they obey a law: `p^e - p + 1` for non-cyclic `p`-groups, `p^(e-1)` for cyclic ones
  (THM-4764), and `|A|/o` jumps at primes dividing class multiplicities;
- "64 quadratic twists" are the `2^7/2` twist classes of the seven-root hourglass surface of Sessions 58–60, and
  the twist average is the branch complement `P^2 \ 7 lines` (THM-4765), whose count polynomial is the
  deconed characteristic polynomial of the line arrangement;
- "the negative coefficient blocks a decomposition into split tori" is Session59's torification test applied
  to the complement; "the useful structure is the multiplicative root construction together with its additive
  branch equations" is exactly right in the following sense: the Kummer construction determines the twist
  average and nothing else; the 63 other Walsh coefficients of the twisted counts are character sums that see
  the linear relations among the seven forms (and show the `1+6` deletion split of Session60's K3 ranks).

Everything was verified from literal arcs and literal field arithmetic, in two modes, with an independent
`2^25` census in C. The packet is `05-knowledge/results/clone_depth_modular_twist_average_20261025/`.

## 2. The themes the repository's novel results live under

A bounded survey (navigation layer plus the October packets; two read-only subagents) suggests that the
repository's proved mechanisms cluster under a small number of recurring operations. I record them as themes
with their strongest instances, because the owner asked for the fit, not only the facts.

**Theme 1 — Average over a finite group, then read the expansion at the critical point.**
Burnside/Pólya averages over odd automorphism groups (A000568, THM-587's per-level signed cycle index with
`P_n(-1) = SC(n)`, THM-4720/4721/4725 twisted and signed Burnside, Session62's `P_A(t)`, THM-305's 2-adic
valuation `(n-1)/2`), twist/character averages (THM-4765; Session59 §6's `(C+3T_2+2T_3)/6` Frobenius
correction), and `q → 1` specializations of counting polynomials (Session59 §2: `D_X(q) = (q-1)^j R(q)` with
`j` the minimum torus rank; the torification tests). The new modular depth law is the `p`-adic layer of this
theme; its lesson is that one lemma ("the depth jumps modulo `ℓ` iff `ℓ` divides the leading coefficient")
governs both the tournament and the surface, while the objects themselves have no map between them.

**Theme 2 — Odd order is a structural fact, not a curiosity.** `|Aut(T)|` odd (THM-048), solvable
(THM-868 §5), inside `A_n` (THM-1420); regular orbits carry invariant tournaments iff the group is odd
(THM-3447); every transitive group of degree `p^2` contains a regular `C_{p^2}` or `C_p^2`
(`center_transitivity_20261021c`), which is why the 50-vertex case has exactly two twins; Session62's
exact root order and this packet's `p`-group law both use odd cycle lengths (`1 + t^m` has a simple zero at
`-1` only for odd `m`) and `p` odd (so `m - v_p(m) ≥ 2` for `m ≥ 2`).

**Theme 3 — Sign-reversing involutions and their toggles.** The clone-cut involution (CLONE_INVOLUTION),
the invariant-toggle group `C_2^r` acting freely on cut orbits (Session62 §2), the parity halves of the
eight-frame cube (Session62/HN8), THM-3181's `Q(-t) = (-1)^n Q(t)`. The repair pattern is always the same: a
labelled toggle does not descend; an orbit-sized invariant toggle does; extra cancellation beyond the toggle
rank is cancellation *between* toggle orbits, which is what the modular depth measures prime by prime.

**Theme 4 — Multiplicative carrier, additive relations.** Session59's graph tori and magic torus (monomial,
`F1`-descendable) with additive tangent/cycle equations, the flat family `G_t` interpolating `G_a` and
`G_m`, the square-root covers as finite-free pullbacks that are not monomial; Session60/Gale: `K = (n-2d-2)H`,
six balanced, seven the first surplus. THM-4765 adds the cleanest finite-field statement of the split:
twist-average = multiplicative data; Walsh spectrum = additive data.

**Theme 5 — Quotients owe a loss ledger; numbers are not bridges.** The repository's hardest-won discipline
(MISTAKE-214/222/227–235/260/261, the Oct-25 "28 ≠ 28", "16 ≠ 16", "3,1,3 ≠ 3+1+3" non-bridges). This session
adds two entries to the ledger of coincidences: `χ(U) = p_g = 7` holds only because the hourglass has three
triple points (18 of 36 seven-subsets fail), and `p^2 - p + 1 = Φ_6(p)` is an artefact of exponent two.

**Theme 6 — Finite exactness with hostile controls as the proof's companion.** Every October packet pairs a
proved mechanism with a complete finite universe and at least one structured adversary; this packet follows
the pattern (22 regular actions, 16 fields × 3 arrangements, 36 subsets, two twins, a mixed-order witness).

Outside tournaments and `F1`, the frontier's other lanes (LRC(14) relation codes and finite banks; GMC(2)
via Frobenius amplification of a balanced face; the planar Jacobian program; Collatz/Mahler block pairs;
Leech deep-hole theta identities) share Themes 1, 5 and 6 and, in GMC(2) and the Leech lane, Theme 2's
`p`-adic/Frobenius amplification. I did not find a way to make Theme 1 bear on LRC(14) that passes the
connection contract, and I record that as a stopping reason rather than a tangent.

## 3. The proposed unifying reading, stated as a card

**Expand at the critical point.** When a count is a projection onto invariants (Burnside average, character
average, twist average) of a product of elementary factors each vanishing at one point (`1 + t^m` at
`t = -1`, `(q-1)^r` at `q = 1`), read three things from its expansion there: (i) the order of vanishing is a
*minimum rank* over the summands (minimum cycle count; minimum torus rank; minimum fixed-space dimension),
with a positive leading coefficient when summands cannot cancel (odd cycles; real representations);
(ii) the quotient by the critical factor has nonnegative coefficients iff there is a *global* free factor
(exterior algebra on the common fixed space; a split-torus stratification), so a negative coefficient is an
obstruction; (iii) modulo a prime `ℓ` the order jumps iff `ℓ` divides the leading coefficient, and for
regular `p`-group actions the jump is `p^e - p^J + p^(J-1)` with `J` read from cyclic-subgroup counts.
Counterindications: even actions, complex scalar actions, mixed-order collisions, and equal critical values
of unrelated expansions.

This is a method, not a theorem about the world; its value is that it told me, before any computation, what
the prompt's numbers had to be, and then what to verify.

## 4. What is open and what I would do next

1. HYP-9374: equality of the candidate-order lower bound for mixed-order regular actions (collisions of
   candidate orders are the only possible failure). A proof would classify when two `p'`-parts can produce
   the same `s`-order and show the leading coefficients cannot cancel.
2. Non-regular gap actions: real tournaments whose automorphism group acts on gaps with several orbits of
   different sizes (Session62's 33-vertex example has four orbits). The `p`-group proof needs element-wise
   cycle types rather than orders; a computation on that example would already be informative.
3. The Walsh spectrum of the twisted hourglass counts as a function of `q`: the six-subset sums are traces
   of Frobenius on the transcendental parts of the seven K3 deletions, so the `1+6` split should be provable
   from Session60's Picard ranks via the Tate conjecture for K3 surfaces (known over finite fields); this
   packet proves only the symmetry half.
4. The period-12 quasi-polynomial of the parallel-pair sums, which Session59's branch addendum would predict
   from the 2-, 3-, 4-torsion of the multiplicative group; a direct inclusion–exclusion proof should be short.

## 5. Cards used, and stopping reasons

Used: "Search the statement before the method" (the prompt's numbers were found in Session62 and Session59
before any derivation), "Respect symmetries by searching orbit representatives" (the half-systems were
searched exhaustively, the gap action verified from arcs), "Type every analogy" (the two sevens), "Quotient
by the actual repair, then test the native survivor" (toggle rank versus depth). Candidate card with evidence
from distinct threads: "Expand at the critical point" (§3 above). It is NOT promoted to `META-PATTERNS.md`:
the bounded startup surface had 13 bytes of slack at the inherited HEAD and `agents/check_docs.py` rejects
any addition; a maintainer must compress an existing entry first. The same budget kept this packet off
`CURRENT-FRONTIER.md`; routing is through the results index, the hypotheses digest and THM-4764/4765.
Stopped: no LRC(14) transfer (no map); no claim that the depth mod `p` is visible in A000568 congruences
(THM-307's `T(p)` mod `p` is a different average).
