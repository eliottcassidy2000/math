---
id: HYP-9172
title: "(a) PROVED: for every n >= 5 no S_n-equivariant map from tournaments to sign-even graphs, or switching classes to sign-even Euler graphs, induces a type bijection. (b) Branch integrality PROVED for every n >= 3; an explicit decorated combinatorial pairing remains OPEN"
status: >
  (a) PROVED for every n >= 5, in both the tournament and switching-class formulations, by the linked
  uniform 3/5 obstruction (2026-10-20; independent audit PASS). The threshold is sharp: both containment
  matching problems have solutions for n <= 4. This concerns S_n-equivariant labelled maps inducing
  type bijections, not arbitrary or differently natural bijections between unlabelled types.
  (b) NUMERICAL INTEGRALITY PROVED for every n >= 3, by Frobenius root-count divisibility and termwise
  2-integrality; n = 2 has the sharp 1/2 + 1/2 exception. An explicit structural decorated pairing for
  n = 0 mod 4 remains OPEN. A fixed-point-free pairing of all unadorned even-Aut types fails at n = 8:
  there are 27 such types, whereas N_lev(8) = 15.
  HISTORICAL FINITE EVIDENCE retained: block certificates for 5 <= n <= 100 (switching classes),
  5 <= n <= 15 plus 16 sporadic n <= 86 (tournaments); exhaustive matchings for n <= 10
  (deficits 1, 1, 2, 3, 4, 10) and n <= 9 (deficits 1, 2, 5, 12, 44). Branch values through n = 16
  remain valid. The old Goldbach/Lemoine-type route is no longer needed for the all-order theorem.
source: collatz-procgen-20260922 session, tbij lane (2026-10-01), conjecture C1 and the open n = 0 mod 4 case; promoted with THM-4531
related:
  - 01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md
  - 01-canon/theorems/THM-479-level-holonomy-and-branch-split.md


---

# HYP-9172 — precise equivariance closure and the remaining pairing question

## (a) Uniform nonexistence, with its exact notion of naturality

The [all-order proof](../results/uniform_equivariance_obstruction_20261020/UNIFORM_EQUIVARIANCE_OBSTRUCTION.md) combines a signed Walsh stabilizer-product identity with containment forced by equivariance. Every odd-prime valuation must then agree at each matched pair. The empty graph would require a source stabilizer with full ambient Sylow-3 and Sylow-5 subgroups. Overlapping single 3-cycles would create a forbidden double transposition; the resulting disjoint-triple partition has insufficient 5-part. This proves both stated cases for every n>=5, without a prime-additive conjecture or group classification, and gives an explicit Hall-deficient family.

**Concurrent contribution preserved.** Incoming commit `b9aac7f88921f2bf18f060fccc8fc5ed80db56a5` from session `mac-mini-2026-10-19-hyp9172` supplies [THM-4722's separate pre-audit proof](../../01-canon/theorems/THM-4722-no-equivariant-bijection-tournaments-even-graphs-for-every-n-sylow-hall-violator.md): the tournament case at every n>=5 via Sylow refinement and Jordan, and the switching case at odd n>=5 via an Euler-parity decomposition. Its theorem, note and exact certificates remain unchanged. Its contemporaneous status left even-n switching and part(b) open; those limitations are superseded here by the independently audited all-n obstruction and numerical-integrality proofs. The separate incoming argument is recorded as provenance, not used as an audited dependency of those proofs.

The earlier finite certificates remain in THM-4531 and its original packet. Arbitrary bijections of unlabelled types still exist by equinumerosity. The broader natural-bijection request in Royle et al. did not impose this precise equivariant-map definition; it is not resolved in every possible sense by this theorem.

## (b) Numerical integrality is settled; a decorated pairing remains open

The [branch proof](../results/switching_branch_integrality_20261020/BRANCH_INTEGRALITY.md) puts the odd-order restricted Burnside average in Z[1/2], using Frobenius's root-count theorem on each stabilizer. The explicit odd-partition formula independently puts it in Z_(2), apart from n=2. Their intersection is Z, proving both branches integral for every n>=3 (also n=1). The established n=2 mod4 pairing is retained in its stated range n>=6.

The [pairing boundary](../results/switching_branch_integrality_20261020/PAIRING_BOUNDARY.md) shows why this arithmetic proof does not produce the requested structural pairing. At n=8 the 27 even-Aut switching-class types cannot support a fixed-point-free involution, and N_lev=15 is not half their number. The old count of 19 types in the Sylow-C2 stratum remains valid. A different, suitably decorated carrier or weighted grouping is still OPEN; failure of these naive carriers does not refute every pairing formulation.
