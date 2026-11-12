---
id: HYP-9669
title: "The 2-adic derivative law of the Hamiltonian path count on a switching class. For a tournament T, a vertex set U and c_U(F) = sum_{S in U} (-1)^(|U|-|S|) F(T^S): v_2(c_U(H)) >= ceil(|U|/2) and v_2(c_U(hc)) >= ceil(|U|/2) - 1 for every non-empty U (hc = number of Hamiltonian cycles). Consequence: the binary digit k of H has F_2-degree at most 2^k on every switching class"
status: >
  OPEN (CONJECTURE). PROVED for |U| <= 4 (THM-5139: |U| = 1, 2 is Redei; |U| = 3, 4 is the quadratic law).
  PROVED reductions (source note, section 2.3): the statement for hc on all odd subtournaments implies the
  statement for H, and for odd order the two are equivalent; the digit lemma (the valuations imply degree at
  most 2^k for digit k); equivalent form: for every k the residue 2^(-k) c_U(H) mod 2, |U| = 2k, is the same at
  every member of the class.
  FINITE-EXACT evidence: every U in random tournaments of orders 6-10 and in one circulant and the transitive
  tournament of order 9 (in sample); pre-registered run R1 (registration pushed at 9f260a24f7 before the runs):
  orders 11 (three random, one circulant), 12 (two random, transitive), 13 (one random, one circulant), 14 (one
  random) - no violation, digits 1, 2, 3 of H of degree exactly 2, 4, 8. For odd order the bound is attained at
  every |U| <= 6; for even order the least valuation at |U| = 5, 6 was 5 (orders 10, 12, 14).
  KNOWN, credited: c_U(H) = 0 for |U| > 2 floor(n/3) (THM-4524); the digit-degree law of a vertex cube
  (THM-5127 (C) (8)), where H is a quadratic polynomial with even coefficients.
  NOT TRUE on the full arc cube (FINITE-EXACT, all labelled tournaments of orders 6 and 7): the mixed arc-flip
  differences of H are only even in general (two disjoint arcs: valuation 1), so the law is special to
  switching cubes and vertex cubes.
  Falsifier: one U with a smaller valuation.
  AUDIT: AUDIT-PENDING (pre-audit checkpoint)
session: mac-mini-2026-11-11-gauge
source: 05-knowledge/results/switching_quadratic_census_gauge_20261111/README.md
related:
  - THM-5139
  - THM-4524
  - THM-5127
  - HYP-9667
---

# HYP-9669 - the 2-adic derivative law on a switching class

Source: [the note](../../05-knowledge/results/switching_quadratic_census_gauge_20261111/README.md), section 2.3;
registration and outcome in `PREREGISTRATION.md` of the packet (R1).

**Conjecture.** For every tournament `T` and every non-empty vertex set `U`, with
`c_U(F) = sum_{S in U} (-1)^(|U|-|S|) F(T^S)`: `v_2(c_U(H)) >= ceil(|U|/2)` and `v_2(c_U(hc)) >= ceil(|U|/2) - 1`.

**Next step.** The case `|U| = 5, 6` for Hamiltonian cycles: `c_U(hc) = 0 (mod 4)`; equivalently
`(1/2) c_{abcd}(hc) mod 2` is a class invariant (a 4-uniform "second endpoint hypergraph"). The proof of THM-5139
is a recursion modulo 2 on the pair (endpoint vector, endpoint graph); a recursion modulo 4 would start from
the exact system, for `e = s + t`, `d = s - t` (integers):
`e_T(u) = <iota_u, d_{T-u}> + H(T-u)`, `d_T(u) = <iota_u, e_{T-u}> - H(T-u)`, `hc(T) - hc(T^a) = <iota_a, d_{T-a}>`,
and `H(T^u) - H(T) = - sum_{i in N^+(u)} w_i(T-u)` with `w_i(R) = 2 d_R(i) + K_R(.,i) - K_R(i,.)` (`K(y,z)` = number of
ordered covers of `R` by two non-empty paths, the first ending at `y`, the second starting at `z`; from the cut
quadratic of THM-5127; checked on 360 cases, orders 4-8). In this form the law for `H` at order `n` is a law with
offset for the vector `w` at order `n - 1`, and the missing step is a symmetry of the next tensor
(`sum_{v in U} (1/2) c_{U - v}(d(v)) = 0 (mod 2)` for `|U| = 4`), the analogue of the symmetry of `Lambda`.
