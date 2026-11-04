---
id: THM-4867
title: "Every uniquely Hamiltonian tournament is a one-arc join. KNOWN (Shelor 1977; Douglas 1970; Garey 1972), recorded with its proof: Shelor's complete maximum-outdegree insertion construction gives a vertex v of in-degree 1, and by converse symmetry a vertex of out-degree 1. Hence ({v}, V - v) is a single-vertex one-arc cut, T = J({v}, v; T - v, p) with e_p(T - v) = 1, and T has arc-strong connectivity 1. PROVED: the equivalences one-arc cut <=> cut arc on the cycle <=> arc-strong connectivity 1, and the completeness for every n of the extension generator (Moon deletion; cyclic-interval out-neighbourhoods, which is klein's same-day block lemma; hc(T' + v) = sum P_T'(O -> I), with THM-4862 (b) the one-vertex case). FINITE-EXACT: the generator reproduces Garey's F(2n - 6) through order 16 (klein's generators through 14); three natural proof routes fail (strong tournament-minus-arc Hamiltonicity at order 4; the 2-arc-strong arc-avoiding route at order 9; degree 1 at non-critical vertices from order 5); PROVED by the S74 audit (given Shelor's grammar): the number of cut arcs on the unique cycle is 3 + #{b, e} - #{c} for a Shelor word, so T(n, k) = sum_d C(n-4, k-3+2d) C(k-2+d, d), with generating function z^2 (x + (1 - x) z)/(1 - 2x - x(1 - x) z) (HYP-9458)."
status: >
  KNOWN: the theorem (Shelor 1977 sections 2-5, as cited in the repository by session103,
  operation_closed_synthesis_20261029 section II, two_entry_funnels_20261103 and klein THM-4862; Douglas 1970; Garey
  1972). The S73 audit identified S73's HYP-9457 (b) as this fact (MISTAKE-705).
  PROVED (elementary): the equivalences and the generator's completeness. Credit: klein's block lemma (two_branch_rearming_sigma_sweep_20261104/PROOFS.md section 6, Lemma 5); THM-4862 (b); the S73 audit's deletion run.
  PROVED by the S74 audit: the cut-arc law and the generating function (machine-checked on 9,227,464 Shelor words of orders 4-20).
  FINITE-EXACT: the independent generation through order 16, the counterexamples, and T(n, k).
  PRE-REGISTERED: HYP-9458 (the cut-arc formulas, confirmed at orders 14-16, then proved).
  AUDIT: PASS WITH CORRECTIONS, all applied (MISTAKE-706); report 04-computation/experiments/uh_one_arc_cut_20261104_audit_S74/AUDIT.md.
session: opus-2026-11-04-S74
source: 05-knowledge/results/uniquely_hamiltonian_one_arc_cut_shelor_20261104.md
related:
  - THM-4866
  - THM-4862
  - THM-4846
  - THM-4744
  - HYP-9457
  - HYP-9458
  - HYP-9459
  - HYP-9460
  - THM-4868
verification:
  - 04-computation/experiments/uh_one_arc_cut_20261104_audit_S74/AUDIT.md
  - 04-computation/experiments/uh_one_arc_cut_20261104.c
  - 04-computation/experiments/uh_one_arc_cut_20261104.out
---

# THM-4867 - every uniquely Hamiltonian tournament is a one-arc join

## Theorem (KNOWN: Shelor 1977)

Let T be a tournament with exactly one directed Hamiltonian cycle C. Then T has a vertex of in-degree 1 and a vertex of
out-degree 1.

**Proof.**
- Orders 3 and 4: the tournaments are C_3 and S_4, where these vertices are visible.
- Order n >= 5: Shelor's construction (bijective and complete) obtains T from a uniquely Hamiltonian tournament of order
  n - 1 by inserting a vertex v that beats every old vertex except its parent p. So v has in-degree 1.
- For out-degree, apply the same to the converse of T, which is uniquely Hamiltonian.

**Corollary.** The set A = {v} has exactly one entering arc, p -> v. Hence:
- T = J({v}, v; T - v, p), the one-arc join of THM-4866 (B);
- the cycle formula gives e_p(T - v) = hc(T) = 1;
- p -> v is a cut arc of C, so the arc-strong connectivity of T is 1.

## Equivalences (PROVED)

For a tournament with a Hamiltonian cycle C, the following are equivalent:
- a one-arc cut exists;
- some arc of C is a cut arc;
- the arc-strong connectivity is 1.

The single backward arc of a one-arc cut must lie on C, and then both sides of the cut are intervals of C.

## Extension lemma (PROVED) and the independent check (FINITE-EXACT)

**Lemma.** Let T be uniquely Hamiltonian of order n >= 4, and let v be non-critical (it exists by Moon 1966). The interval part is klein's block lemma.
- T - v is uniquely Hamiltonian.
- N^+(v) is a proper cyclic interval of its cycle C'.
- hc = sum over w in O and u in I of P_(T - v)(w -> u).

**Check.** Generating by this lemma reproduces Garey's counts through order 16, without using Shelor's grammar. Every
class has vertices of in-degree 1 and of out-degree 1.

## Routes that fail (FINITE-EXACT)

- **Strong tournament minus an arc.** Such a digraph can be non-Hamiltonian. At order 4 this happens in 24 of the 72
  labelled strong cases.
- **2-arc-strong tournaments.** "T - e is Hamiltonian for every arc e" holds through order 8 but fails in 4 (class, arc)
  pairs at order 9.
- **Non-critical vertices.** A non-critical vertex of a uniquely Hamiltonian tournament can have in- and out-degree
  >= 2, from order 5 on.

## Cut arcs (PROVED by the S74 audit, given Shelor's grammar; HYP-9458)

**Count law.** For a Shelor word, the number of cut arcs on C is k = 3 + #{b, e} - #{c}. One insertion keeps exactly
the old cut arcs whose largest one-arc cut contains the parent, and every such cut contains the newest vertex.

**Closed form.** With m = n - 4,

    T(n, k) = sum_d C(m, k - 3 + 2d) C(k - 2 + d, d),
    sum T(n, k) x^(n-4) z^k = z^2 (x + (1 - x) z) / (1 - 2x - x(1 - x) z).

**Consequences.**
- T(n, 2) = 2^(n-5);
- T(n, 3) = n 2^(n-7), which is A045891;
- T(n, 4) = (n^2 + n - 8) 2^(n-10), which is A158920;
- the diagonals are 1, n - 4, C(n - 3, 2), and then n(n - 4)(n - 5)/6 (A005581);
- k >= 2 for every n.

**Checks.** Pre-registered and confirmed at orders 14-16; machine-checked through order 20.

**Concurrent.** kind-pasteur-S10's HYP-9460 (same day, 12:28) conjectured the equivalent recurrence
T(n, k) = T(n - 1, k - 1) + 2 T(n - 1, k) - T(n - 2, k - 1), the denominator of the generating function above. It also
announces a proof of this theorem that does not use Shelor's grammar (THM-4868, pending).
