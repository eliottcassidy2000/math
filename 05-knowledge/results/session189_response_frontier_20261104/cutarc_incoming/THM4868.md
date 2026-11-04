---
id: THM-4868
title: "The non-critical run of a uniquely Hamiltonian tournament. For T of order n >= 4 with exactly one Hamiltonian cycle C: (A) every non-critical vertex v (T - v strong) beats exactly the d+(v) vertices that follow it on C, and d+(v) >= 2 makes v+ non-critical (dually d-(v) >= 2 makes v- non-critical); (B) hence T has a sole exit and a sole entry (KNOWN: Shelor 1977, THM-4867; new Shelor-free proof by roundness); equivalently, minimum out-degree >= 2 forces at least two Hamiltonian cycles, sharp at RT5; (C) the non-critical vertices form ONE run of C from a sole entry y to a sole exit x, 2 <= r <= n - 1, and y is the canonical Shelor vertex; (D) an arc of C is a cut arc iff it is not internal to the run, so #cut arcs + #non-critical = n + 1 and HYP-9458's T(n, k) = N(n, n + 1 - k); (E) T is asymmetric; (F) F(2n - 8) and F(2n - 7) classes with one and two sole exits, and F(2n - 10) + 1, F(2n - 9) - 1, F(2n - 9) - 1, F(2n - 8) + 1 by (sole exits, sole entries) = (1,1), (1,2), (2,1), (2,2); (G) the admissible Shelor parents (e_p(T) = 1) are exactly y, y- and, when y- is a sole entry, y--, so T -> (T - y, y-) is a bijection onto admissible pairs one order down and |UH(n)| = F(2n - 6) (Garey's theorem re-derived); (H) the child of an admissible parent y, y-, y-- has run length r, r + 1, r + 2, so N(n, r) = N(n - 1, r) + 2 N(n - 1, r - 1) - N(n - 2, r - 1) and every formula of HYP-9458 holds"
status: >
  PROVED (elementary; kind-pasteur-2026-11-04-S10): (A), (C), (D), (E), (F), (G), (H);
  - (B) is a new proof of a KNOWN theorem (Shelor 1977), recorded first in this repository as THM-4867 by
    opus-2026-11-04-S74, which this theorem defers to for the statement;
  - (G) re-derives KNOWN results: Garey 1972's F(2n - 6) (also Moon 1982) and the completeness of Shelor 1977's
    insertion;
  - (H) gives a second, Shelor-free proof of opus-S74's HYP-9458 and of this session's pre-registered HYP-9460 P7.
    HYP-9458 was conjectured and confirmed at orders 14-16 by opus-S74, and first PROVED by its independent audit
    (a0518f53cd), given Shelor's grammar: the count law k = 3 + #{b, e} - #{c}, whose parent letters j = 0, 1, 2 are
    (G)'s y, y-, y--. That audit noted that P7 follows from it with (C).
  CITED: Moon 1966 (a strong tournament of order n >= 4 has a strong subtournament of order n - 1), Camion 1959,
  Garey 1972 (F(2n - 6); re-derived in (G), used nowhere), klein THM-4862 (b) (entry attachment C(U + a_entry) =
  e_r(U), re-proved in the note as Lemma 3), Las Vergnas 1975 via Burzio-Demaria 1988 (the bineutral A_n is the only Hamiltonian tournament
  with exactly two neutral vertices; used only for the top diagonal in (D)).
  FINITE-EXACT: every statement on all F(2n - 6) classes at orders 4-18 (orders 14-18 pre-registered as HYP-9460 in
  a8db6218f9 before the run).
  FINITE-EXACT also: (G) and (H) on every parent of order 4-14 (shelor_parents_20261104.c).
  NOVELTY: not claimed against Douglas 1970, Shelor 1977 or Moon 1982, whose bodies could not be read (paywall, CAPTCHA,
  Cloudflare challenge); new to the repository.
  AUDIT: see the note, section 9.
session: kind-pasteur-2026-11-04-S10
source: 05-knowledge/results/uniquely_hamiltonian_noncritical_run_20261104.md
related:
  - THM-4867
  - THM-4866
  - THM-4862
  - THM-4846
  - THM-4744
  - HYP-9460
  - HYP-9458
  - HYP-9459
  - HYP-9457
verification:
  - 04-computation/experiments/uh_noncritical_run_20261104.c
  - 04-computation/experiments/uh_noncritical_run_20261104_explore.out
  - 04-computation/experiments/uh_noncritical_run_20261104_confirm.out
  - 04-computation/experiments/shelor_parents_20261104.c
  - 04-computation/experiments/shelor_parents_20261104.out
---

# THM-4868 - the non-critical run of a uniquely Hamiltonian tournament

Let T be a tournament of order n >= 4 with exactly one Hamiltonian cycle C. Notation:
- v+ and v- are the successor and predecessor of v on C;
- d+ and d- are out- and in-degree;
- NC is the set of non-critical vertices v, those with T - v strong;
- a sole exit has d+ = 1, and a sole entry has d- = 1.

## Statements

- **(A) Local roundness and propagation.**
  - If v is non-critical, then N+(v) = {v+, v++, ..., v^(+d+(v))}.
  - If moreover d+(v) >= 2, then v+ is non-critical. Dually, d-(v) >= 2 makes v- non-critical.
- **(B) Degree-1 theorem** (KNOWN: Shelor 1977; THM-4867). T has a sole exit and a sole entry. Equivalently, a strong
  tournament with minimum out-degree >= 2 has at least two Hamiltonian cycles. Two is attained by RT_5 (n = 5). At
  orders 6-10 the minimum is 3 (klein's exhaustive HYP-9459 table, FINITE-EXACT); raising the bound to 3 for all
  n >= 6 is HYP-9459 (1)'s open two-cycle half.
- **(C) One run.**
  - NC is one run y -> ... -> x of C, of length 2 <= r <= n - 1, where y is a sole entry and x is a sole exit.
  - For every sole exit z, z+ is critical; for every sole entry z, z- is critical.
  - y is the only non-critical sole entry. So T is Shelor's insertion of y into the uniquely Hamiltonian tournament
    T - y with parent y-, canonically.
- **(D) Cut arcs.**
  - An arc of C is a cut arc iff it is not internal to the run.
  - So #cut arcs + |NC| = n + 1, and T(n, k) (HYP-9458) = N(n, n + 1 - k).
  - With Las Vergnas, T(n, n - 1) = 1.
- **(E) Asymmetry.** Aut(T) is trivial.
- **(F) Refined counts** (n >= 5).
  - Exactly F(2n - 7) classes have two sole exits. They are in bijection with rooted pairs (T', sole exit) of order
    n - 1.
  - By (number of sole exits, number of sole entries):

    | pair | (1, 1) | (1, 2) | (2, 1) | (2, 2) |
    |---|---|---|---|---|
    | number of classes | F(2n - 10) + 1 | F(2n - 9) - 1 | F(2n - 9) - 1 | F(2n - 8) + 1 |

  - There are F(2n - 5) rooted pairs (T, sole exit).
- **(G) Admissible parents.** For p in T, let T^p add a new vertex beating everything except p, with p -> new. Then
  T^p is uniquely Hamiltonian iff e_p(T) = 1 (Hamiltonian paths ending at p).
  - These p are exactly y and y-, together with y-- iff y- is a sole entry.
  - So T -> (T - y, y-) is a bijection from UH(n + 1) onto the admissible pairs of UH(n).
  - Hence |UH(n)| = F(2n - 6). This is Garey's theorem, re-derived without Douglas's characterization.
- **(H) Run transitions.** The children T^y, T^(y-) and T^(y--) have run lengths r, r + 1 and r + 2. Therefore
  N(n, r) = N(n - 1, r) + 2 N(n - 1, r - 1) - N(n - 2, r - 1) for n >= 6, and HYP-9458's T(n, k) satisfies
  T(n, k) = T(n - 1, k - 1) + 2 T(n - 1, k) - T(n - 2, k - 1). The latter gives 2^(n-5), n 2^(n-7),
  (n^2 + n - 8) 2^(n-10) and the diagonals 1, n - 4, C(n - 3, 2).

## Proof (outline; full proofs in the note, section 2)

**Insertion.** Let v be non-critical.
- Every Hamiltonian cycle of T - v takes v at each I -> O switch of v's in/out word along it.
- hc(T) = 1 forces a unique cycle C' = C - v + (v- -> v+) and a single switch, located at v. This gives (A)'s first
  claim.
- If d+(v) >= 2, then v -> v++, and the cycle v -> v++ -> ... -> v- -> v shows that T - v+ is strong.

**(B).** Suppose minimum out-degree >= 2.
- Moon's theorem gives a non-critical vertex. Propagation then makes every vertex non-critical, so T is round with all
  steps >= 2.
- At odd n the constant +2 walk is a second Hamiltonian cycle.
- At even n >= 6, some vertex w has out-degree >= 3 (the degrees sum to n(n - 1)/2 > 2n). The walk with steps +2, one
  +1 and one +3 from w is a second Hamiltonian cycle.

**(C).** A sole exit z makes z+ critical, since z becomes a sink of T - z+.
- Each maximal run ends at a sole exit (by propagation).
- There are at most two sole exits, and if two, they are z and z+.
- So there is exactly one run.

**(D).**
- An internal run arc u -> u+ has u -> u++ and u- -> u+ (by (A) and (C)), so deleting it keeps T strong.
- At a critical w, write T - w = Q_1 => ... => Q_m. Any second arc from w into Q_1 would give a second Hamiltonian cycle,
  so w -> w+ is the only arc entering Q_1. Dually, w- -> w is the only arc leaving Q_m.

**(E).** An automorphism is a rotation of C fixing the set of sole exits.

**(F).** Map a class with two sole exits x -> y to (T - x, y), and back by attaching a new sole exit before y. This
gives b(n) = a(n - 1) + 2 b(n - 1); then use (G)'s total. The bijection preserves the number of sole entries.

**(G).** Induction on the canonical deletion T = U^(p_0), with U = T - y and p_0 = y-. For p in U, classify the
Hamiltonian paths of T ending at p by the position of y:

    e_p(T) = [e_p(U) - hp_U(p_0, p)] + W(p_0, p),

where W counts two-path covers of U ending at p_0 and at p.
- The cycle path of U ending at p contributes to both terms unless it starts at p_0, i.e. p = p_0-. So every p outside
  {y, p_0, p_0-} has e_p(T) >= 2.
- For p = p_0-: if p_0 is a sole entry, e_p(T) = hc(U) = 1. Otherwise p_0- is not admissible in U (induction), and
  W >= max(1, hp) gives e_p(T) >= 2.

**(H).** In T = U^p, every non-critical vertex of U other than p stays non-critical. p becomes critical, and x_U+
stays critical (x_U is still a sole exit). When p = y_U--, the vertex y_U- becomes non-critical, because p lies in the
last strong component of U - y_U-. Eliminating the count of two-sole-entry classes between the two resulting recursions
gives the recurrence. QED.

## Verification (FINITE-EXACT)

`uh_noncritical_run_20261104.c` generates every class by interval extension (Moon plus insertion; Shelor-free).
- It reproduces F(2n - 6) at orders 4-18.
- It checks (A), (C), (D), the degree-1 counts, (E) and (F) on every class.
- It confirms hc = 1 by an independent subset DP, and one-arc cuts by brute force over all subsets, at orders <= 14.
- Orders 14-18 were pre-registered (HYP-9460).
- `shelor_parents_20261104.c` checks (G) and (H) on every parent of order 4-14.
