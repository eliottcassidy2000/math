---
id: HYP-9170
title: "For every odd m, some anti-circulant tournament on 2m vertices has every arc on an odd number of Hamiltonian paths (P1); moreover, for every odd m there are infinitely many primes p = 3 mod 4 with QR_p restricted to mu_2m all-odd (P2). P1 would prove the existence half of HYP-9167(b)"
status: >
  OPEN.
  PROVED (THM-4529): anti-circulants exist only for N = 2 mod 4, and their antipodal arcs are always odd. P2 holds for m = 3
  (exactly the p = 7 mod 8) and for m = 5 (exactly the p = 3 mod 8), by Legendre-sign algebra plus the finite
  classification at N = 6, 10.
  FINITE-EXACT (lane; audited for N <= 14): all-odd anti-circulant classes number 1, 1, 1, 2, 4, 2 for N = 6, 10, 14, 18, 22, 26.
  Examples: P7 - v, QR11 - v, QR127[mu_14], QR19 - v, QR_p[mu_22] for p = 23, 199, 727, 1783, QR27 - 0, QR2003[mu_26].
  First open cases: N = 30 (QR31 - 0, i.e. HYP-9167(a) at q = 31) and N = 34 (no Paley tournament minus a vertex exists),
  both beyond the 2^N subset DP.
  By Chebotarev, P2 for a fixed m is a finite question about Frobenius classes of cyclotomic units.
source: collatz-procgen-20260922 session, petersen lane (2026-10-01), proposals (P1), (P2) of section 6.5; promoted with THM-4529
related:
  - 01-canon/theorems/THM-4529-petersen-and-heawood-families-in-paley-coordinates-anti-circulant-all-odd-tournaments.md
  - 05-knowledge/hypotheses/HYP-9167-paley-minus-a-vertex-is-redei-rigid.md
---

# HYP-9170 — all-odd anti-circulant tournaments

**Definition.** A tournament on `N` vertices is *anti-circulant* if it has an anti-automorphism that is a single `N`-cycle; that is, it maps the tournament onto its reverse. Equivalently it is `T_s` on `Z/2m`: `x -> y` iff `s(y - x) = (-1)^x`, for a sign function with `s(-d) = -(-1)^d s(d)`.

**Why they are natural.**
- They are exactly the tournaments in which the all-odd question has a built-in partial answer: the antipodal arcs are always odd (THM-4529, Theorem 6.1).
- Every Paley tournament minus a vertex is anti-circulant.
- So is every `QR_p` restricted to a cyclic subgroup of order `2m`, when `p = 3 mod 4`.

**(P1)** For every odd `m >= 1` some anti-circulant `T_s` on `2m` vertices is all-odd. This would give HYP-9167(b)'s existence half: all-odd tournaments exist for every `N = 2 mod 4`.

**(P2)** For every odd `m` there are infinitely many primes `p = 3 mod 4` with `2m | p - 1` and `QR_p[mu_2m]` all-odd.

**Evidence.**
- Exhaustive censuses up to `N = 26`.
- Proofs of (P2) for `m = 3, 5`.
- The 14-vertex example `QR127[mu_14]`, built from the Mersenne prime `2^7 - 1` and the `x2`-codes of length-7 parity words.

**UPDATE 2026-10-01 (THM-4532): an explicit candidate family for (P1).**
The two-sheet clocks `D_r` are anti-circulant: the sheet swap composed with a shift is a single `2r`-cycle anti-automorphism. They are all-odd for every odd `r <= 19` (`N <= 38`), and `D_9` is the second all-odd class at `N = 18` in THM-4529's census.

**Conjecture C1 (the explicit form of P1):** `D_r` is all-odd for every odd `r`.
- Only the vertical arcs are proved odd so far (THM-4532 T1).
- The non-vertical arcs lie in orbits of size `2r` and need a new mechanism.
- Their parities are not affine in the cross set over GF(2).
