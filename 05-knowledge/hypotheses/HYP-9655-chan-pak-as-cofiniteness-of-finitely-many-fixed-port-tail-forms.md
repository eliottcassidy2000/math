---
id: HYP-9655
title: "Chan-Pak as cofiniteness of finitely many fixed-port tail forms. (A) Fixed-port tails over cores of order <= 7 realise every odd value except 7 and 21 (base 2) and every value = 1 (mod 3) except 10 and 31 (base 3); each implies the corresponding Chan-Pak conjecture. (B) Cores of order <= 6 leave infinitely many misses, so the core order needed is at least 7. (C) Question: does the one-apex share of odd numbers per dyadic window converge (observed 0.563 at 2^6, 0.5006 at 2^16, 0.497 at 2^19)?"
status: >
  OPEN.
  (A): FINITE-EXACT below 10^6 in both bases (THM-5085 (4), pre-registered Q1/Q2 CONFIRMED).
  (B): evidence only. There are 18 sparse misses in (30000, 10^6) at order <= 6, and no closed-form reason is known.
  (C): an open question, with data only (W < 2^22).
session: kind-pasteur-2026-11-09-S19
source: 05-knowledge/results/kps_s19_fixed_port_tails_20261109/README.md
related:
  - THM-5085
  - HYP-9640
  - HYP-9530
  - THM-5070
  - THM-1370 (h-spectrum-omits-7-21-all-n)
---

# HYP-9655 - Chan-Pak as cofiniteness of finitely many fixed-port tail forms

## (A) Finitely many forms suffice

**Statement.** Let C_7 be the finite set of core states (h0, R0, Delta0, h_u) of all tournament classes of order <= 7
with every port. There are 2,831 distinct states in base 2 and 2,901 in base 3. By THM-5085 (2), each value of a
fixed-port tail over a core in C_7 is

    h0 + (x/2) [ R0 N + Delta0 [N odd] - h_u Q(N) ],

where N ranges over the positive integers and Q(N) is the pair sum of THM-5085 over the 1-bits of N. (The values also
include h0 itself, with no bridges.) The conjecture: these values are
- in base 2, every odd number except 7 and 21;
- in base 3, every number = 1 (mod 3) except 10 and 31.

**Why it implies Chan-Pak.** Each value is H (base 2) or I(Omega, 3) (base 3) of an explicit tournament.

**Evidence.** Both statements hold below 10^6.

**Falsifier.** A value v < 10^9, say, that no form reaches. A computation can test this.

**Route.** The forms are quadratic in the bits of N. A proof might start from a handful of cores whose drifts h_u and
offsets cover all residues, then control Q(N) along explicit bit patterns. The minimum number of cores needed on
[1001, 2*10^5] is between 17 and 81 (THM-5085 (4)).

## (B) Order 6 is not enough

**Statement.** Infinitely many odd values are missed by fixed-port tails over cores of order <= 6.

**Evidence.** 18 misses in (30000, 10^6), with the smallest 34,869, and no sign of thinning to zero.

**Falsifier.** A proof that the order-<= 6 forms are cofinite.

## (C) The one-apex share (question, not a conjecture)

The one-apex tournaments (the core {u}) attain H = 1 + 2a(W). The share of odd numbers in [2^k, 2^(k+1)) that they
attain, computed from all W < 2^22, is:

| k | 6 | 8 | 10 | 12 | 14 | 16 | 17 | 18 | 19 |
|---|---|---|---|---|---|---|---|---|---|
| share | 0.5625 | 0.5547 | 0.5391 | 0.5127 | 0.5023 | 0.5001 | 0.4976 | 0.4992 | 0.4966 |

It is open whether the share has a positive limit, tends to 0, or oscillates. The overall share in [1001, 200000] is
49,728 / 99,500.
