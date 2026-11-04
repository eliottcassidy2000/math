---
id: THM-4866
title: "The single-cycle floor, the bottom band and one-arc joins. (A) CITED by combination, with an independent cut proof: every interval chain IC(p) has exactly one Hamiltonian cycle, so every Busch-extremal strong tournament is uniquely Hamiltonian. (B) PROVED: the one-arc join J = A => B with a -> b reversed has H(J) = x_h y_h + x_R y_R - x_e y_e (a Lorentzian inner product of rooted triples x = (h, R_a, e_a), y = x(B^op, b)) and hc(J) = s_a(A) e_b(B). (C) FINITE-EXACT through order 10: the cycle ceiling H >= f(n) + 2(hc - 1) with equality only at the floor; second floor f_2(n) = 13, 23, 41, 69, 115, 205 (n = 5..10) below f(n + 1); the bottom band [f(n), f_2(n)) is uniquely Hamiltonian, and all 610 uniquely Hamiltonian classes of order <= 10 are one-arc joins; the crevasses 7 and 21 are holes of the uniquely Hamiltonian stratum (7 between U(4) = {5} and f(5) = 9; 21 the hole of U(6) = {15, 17, 19, 23}). (D) Lorentzian joins with both sides of order <= 7 represent every odd number in [23, 84705], with a 3-adic singular series. (E) PROVED: a self-converse tournament has H = H(core) x s^2, so H = 5 never occurs at odd order and H = 3 never at even order."
status: >
  CITED: (A)'s conclusion, by combination of three sources:
  - Moon-Yang 2007 (extremal = interval chains; THM-4744 (4), THM-4745);
  - THM-4744 (2) (interval chains have n - 2 cyclic triples);
  - Burzio-Demaria, J. Graph Theory 14 (1990) 663-672, as stated by Savchenko, arXiv:2311.13377 (strong tournaments
    with n - 2 cyclic triples have exactly one Hamiltonian circuit).
  PROVED (elementary): the cut proof in (A), (B) and (E).
  FINITE-EXACT: (C) and (D), with program scf_crevasses_20261104.c and its census, chains and joins outputs.
  PRE-REGISTERED: HYP-9457, with P1-P5 confirmed and P6 refuted (singular-series ratio 1.37, not in [1.6, 3.0]).
  CITED: THM-1370 (the impossibility of 7 and 21), THM-4846 (Douglas 1970; Garey F(2n - 6)), and klein's THM-4865 /
  HYP-9453 (root response B_t; near-join; (B) contains the tower step and the near-join as special cases).
  AUDIT: see the note, section 10.
session: opus-2026-11-04-S73
source: 05-knowledge/results/single_cycle_floor_crevasses_lorentzian_joins_20261104.md
related:
  - THM-1370
  - THM-4744
  - THM-4745
  - THM-4846
  - THM-4862
  - THM-4865
  - THM-4863
  - HYP-9453
  - HYP-9454
  - HYP-9457
verification:
  - 04-computation/experiments/scf_crevasses_20261104.c
  - 04-computation/experiments/scf_crevasses_20261104_census.out
  - 04-computation/experiments/scf_crevasses_20261104_chains.out
  - 04-computation/experiments/scf_crevasses_20261104_joins.out
---

# THM-4866 - the single-cycle floor, the bottom band and one-arc joins

## (A) Interval chains are uniquely Hamiltonian

**Statement.** Let IC(p) be the interval chain (THM-4744). Its only backward arcs are the block arcs
x_(l_(i+1)) -> x_(l_i), and the blocks [l_i, l_(i+1)) partition the cuts. Then IC(p) has exactly one Hamiltonian
cycle.

**Proof.**
1. Each cut is spanned by exactly one backward arc.
2. A Hamiltonian cycle crosses each cut exactly once in each direction, so it uses every back arc.
3. Its forward arcs tile the line, with tails forced: every vertex except the back-arc tails.
4. So at most one Hamiltonian cycle exists. Camion gives one.

**Single-cycle floor.** With Moon-Yang, every strong tournament with H = f(n) has exactly one Hamiltonian cycle. This
conclusion also follows from Burzio-Demaria (1990) and THM-4744 (2), as recorded in the status. It is checked on every
floor class through order 10, and on all 1,601 interval chains of order <= 17.

## (B) The one-arc join law

For J = A => B with a -> b reversed, and p(empty) = 1:

    x(A, a) = (h(A), R_a(A), e_a(A)),    R_a(A) = sum over X subset V(A) - a of p(X) s_a(A - X),
    y(B, b) = x(B^op, b),

    H(J)  = x_h y_h + x_R y_R - x_e y_e,
    hc(J) = s_a(A) e_b(B)    (|A| + |B| >= 3).

**Proof.**
1. A Hamiltonian path crosses back from B to A at most once, through b -> a.
2. Paths that never cross back contribute h_A h_B - e_a s_b.
3. Paths that cross back once, A_1 B_1 (b -> a) A_2 B_2, contribute R_a(A) R'_b(B).
4. A Hamiltonian cycle uses b -> a once. What remains is an A-path from a followed by a B-path to b.

**Special cases.** klein's T_1 = h - e + B (B a point) and the near-join h + 2^k B (B a transitive chain), from
THM-4865 / HYP-9453; Moon's NT_s.

**Checks.** The law is checked on 9,000 random joins.

## (C) Bottom band, cycle ceiling, crevasses (FINITE-EXACT through order 10)

**Definitions.** f_k(n) = min{H : hc = k}, and U(n) is the uniquely Hamiltonian stratum.

**Cycle ceiling.** H(T) >= f(n) + 2(hc(T) - 1) for every strong T of order <= 10, with equality only at the floor
classes. At odd order this is the deck form G = H - 2hc >= f(n) - 2 of THM-4863.

**Second floor.** f_2(n) = 13, 23, 41, 69, 115, 205 for n = 5..10, and f_2(n) < f(n + 1) in each case.

**Bottom band.** The band [f(n), f_2(n)) consists of uniquely Hamiltonian classes. The number of such classes at orders
3-10 is 1, 1, 3, 8, 21, 55, 144, 377 (Garey), and every one of them has a cut crossed backward by exactly one arc. By
(B) each is therefore a one-arc join with s_a = e_b = 1.

**The crevasses.**
- 7 lies between U(4) = {5} and f(5) = 9.
- 21 is the only hole of U(6) = {15, 17, 19, 23} below f(7) = 25.
- The other band holes (47, 55, 67 at order 8; 77, 79, 83 at order 9; 13 values at order 10) are filled at other
  orders or as products.

The impossibility of 7 and 21 is THM-1370's.

## (D) Coverage by joins (FINITE-EXACT)

Both sides of the join are drawn from the 2,732 distinct rooted triples of order <= 7.
- **Coverage.** Every odd N in [23, 84705] is a join value; 7 and 21 are not.
- **Comet.** The minimum representation count per dyadic window rises to 198 on [8192, 16384).
- **Singular series.** The mean count on N ≡ 0 mod 3 is 1.37 times that on the other residues (window [25000, 50000)).
  Residues mod 5 and mod 7 are flat.

## (E) Real points

A self-converse T has a palindromic condensation, so H(T) = H(M) (prod over pairs of H(C_i))^2 with |T| ≡ |M| (mod 2),
where M is the self-converse middle component (or empty). Hence:
- H = 5 is impossible at odd order;
- H = 3 is impossible at even order.

## Not proved here

- That every uniquely Hamiltonian tournament has a one-arc cut (HYP-9457 (b)).
- The cycle ceiling for every n (HYP-9457 (a)).
- Chan-Pak completeness (OPEN; THM-4094, HYP-9453, HYP-9454).
