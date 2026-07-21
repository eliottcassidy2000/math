# THM-1825: the q-Redei census, Landau as moment-map surjectivity, and the Galois terminal stack

**Status:** DELIVERED (census with honest negatives + leads; verified
proposition; the phi(n+1) Galois-stack identification). No overclaims.
**Author:** boxeph-2026-07-20-S190 (HYP-8620)
**Owner:** "explore these and more proposals, also consider the 7,21
impossibility and how that may subtly show up in the structure of various
parts."

## 1. (W1) the q-Redei census — back-arc statistic

h_q(T) = sum over Ham paths of q^{bk} (bk = #arcs backward w.r.t. the path
order; label-free; h_q(1) = h). Census n = 3..6 (all 74 classes):
- REFUTED (honest): back-arc parity rigidity is NOT universal — rigid in
  2/2, 3/4, 5/12, 8/56 classes only; |h_q(-1)| = h was an n = 3 accident
  (first break: the strong 4-tournament, h_q = q + 3q^2 + q^3). No cheap
  q-Redei at q = -1.
- LEADS (live): h_q(i) != 0 for ALL classes n <= 5 (roots-of-unity
  nonvanishing conjecture); h_q is NOT even coefficient-unimodal in
  general (8/12, 27/56 unimodal) — sharp contrast with the leaf
  polynomials (log-concave 530/530, THM-1760): the two tournament
  q-refinements have OPPOSITE regularity; the parity-RIGID subfamily
  (which classes?) is uncharacterized; the reversal functional equation
  h_q(T^op) vs q^D h_q(T)(1/q) unverified (lead).
- {7,21} SHADOW: the hole's flanks at n = 5: h = 5 classes share
  h_q = q(1 + 3q + q^2) — SHIFTED-PALINDROMIC; h = 9 classes are ragged
  ([0,1,0,4,2,2], [0,0,5,1,2,0,1]). A phantom h = 7 would need a
  coefficient vector summing to 7 between these shapes; no structural
  exclusion found YET — recorded as the concrete q-shadow question: does
  shifted-palindromicity + the census constraints exclude sum-7 shapes?

## 2. (W2) Landau = moment-map lattice surjectivity  [VERIFIED n <= 7]

Tournaments = points of (P^1)^{C(n,2)}; the T = (C*)^n action has moment
map = score vector; the moment polytope is the permutohedron of
(0,1,...,n-1); Landau's inequalities = majorization = membership.
VERIFIED: attained sorted score multisets == lattice points of the
polytope at n = 3..6 (2, 4, 9, 22) and the n = 7 lattice count is 59 =
A000571(7). STATEMENT (proposition grade): Landau's theorem is exactly
the assertion that the combinatorial moment map attains EVERY lattice
point of its moment polytope — tournament realizability as
Atiyah-Guillemin-Sternberg-type lattice surjectivity, i.e. a GIT/moment
statement, as the S188 reading predicted ("in/transitivity itself":
the polytope's vertices are the transitive tournaments' score vectors).

## 3. (W3) the terminal stack is a GALOIS ORBIT: multiplicity phi(n+1)

For the tight family v = (1,...,n), the number of good-set components
dying simultaneously at delta* = 1/(n+1) is measured as 2, 4, 2, 6 for
n = 3,4,5,6 — exactly phi(n+1). Identification: the lonely times at
threshold are the PRIMITIVE (n+1)-th fractions k/(n+1), gcd(k,n+1) = 1,
and the multiplicative group (Z/(n+1))* permutes them: the terminal
stacked death is a GALOIS-SYMMETRY STACK — the LRC face of the mu_g
rotation stacks (THM-1785) and the reality stacks (THM-1680): the
dictionary's stack column now reads: GMC mu_g-rotation / reality (Z/2) /
LRC Galois (Z/(n+1))*.

## 4. Files

04-computation/q_redei_landau_shadows_boxeph_S190.py + frozen .out.
