# Exact correction in a soluble recursive boundary

**Scoped mathematical results and executable certificates; Chan–Pak OPEN.**

This lane turns a fixed tournament quotient with transitive child modules into
an exact rational walk series after a Stirling change of coordinates. The
underlying quotient substitution law is inherited; the new junction is the
recurrence-class dichotomy, geometric-moment correction certificate, and
complete finite target inverse for this family.

- A fixed quotient's multinomial-weighted Hamiltonian aggregate is P-recursive
  **iff the quotient is transitive**. Cycles create infinitely many complex
  poles after the exponential coordinate change, despite a finite run recurrence.
- Every regular quotient of the same order shares this aggregate. This exact
  collapse is not a classification of those graphs or their individual targets.
- Conjugate pole corrections with directed intervals recover every integer
  digit of the regular-quotient aggregate. Relative accuracy alone fails badly.
- Every nontrivial strong quotient Q has a complete target inverse over the box
  `1 <= module_size <= floor(log2(c-1))`. It yields literal matrices when found.
  This is a complete fixed-family decision procedure, not global coverage.
- The C3 face `H=2+(2^a-1)(2^b-1)` gives a productive Mersenne-factor inverse;
  source witnesses compose multiplicatively under ordinal sum. For example
  23 and47 have sizes(2,3,1) and(2,4,1); 45 has(2,2,2). The complete C3
  transitive-inflation family cannot produce25,27,29,31 at any size.

For the C3 weighted aggregate (not A000568), the signed leading-pole errors
and the number of conjugate pairs needed by the certificate are:

| n | leading approximation minus exact, rounded display | certified pole pairs |
|---:|---:|---:|
| 10 | 0.01581312603326517 | 1 |
| 20 | -82.41454270212493 | 2 |
| 30 | -31124666.24226621 | 2 |
| 50 | -8.267449170263519E+23 | 3 |
| 100 | 3.659404525283840E+77 | 6 |
| 250 | 6.183785750894861E+291 | 15 |

[Proof and exact scope](PROOF.md). [Full certificates](RESULTS.json).
[Literal target witnesses and complete fixed-family empty results](TARGETS.json).
Run `python3 recurrence_bridge.py` with Python3, no external packages.
Tests:5,357 independent subset-DP inflation controls,8 aggregate controls,
16 Mersenne-face controls,18 directed integer recoveries, and10 complete
C3 target inverse queries. No A000568 numerical claim is made by this lane.

The remaining constructive gap is explicit: a finite bank of fixed strong
quotients with transitive modules has only polylogarithmically many targets
below X. Progress beyond this soluble boundary requires growing quotients,
nontransitive children, or nonuniform reusable interfaces. The exact run
compiler can absorb the inherited full cover response; its aggregate cannot
replace that response without losing the prescribed target.
