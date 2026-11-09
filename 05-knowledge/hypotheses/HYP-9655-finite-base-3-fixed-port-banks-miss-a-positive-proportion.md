---
id: HYP-9655
title: "Base-3 fixed-port banks: the ternary counterpart of HYP-9480 (C). (A) Every finite bank of base-3 fixed-port tail forms h0 + (3/2)(R0 N + Delta0[N odd] - h_u Q(N)) misses a positive proportion of the integers = 1 (mod 3), so no finite bank proves ternary Chan-Pak; the order-<= 7 bank nevertheless misses only 10, 31 below 10^7. (B) Data for HYP-9480 (C) in base 2: the order-<= 6 bank misses 169 odd values below 10^7 at a stable density of about 3.3e-5, and the order-<= 7 bank misses none besides 7, 21. (C) The one-apex share question is HYP-9480 (A)'s root case (conjectured limit 1/2)."
status: >
  OPEN.
  (A): conjecture, by the Poisson heuristic of HYP-9480 (C) transferred to base 3. Data: the order-<= 6 base-3 bank
  misses 15 values below 10^6 (10, 31, 5401, 40597, 97282, ..., 904876); the order-<= 7 bank misses only 10, 31 below
  10^7 (THM-5085 (4)).
  (B): data, not a separate conjecture; it is evidence for HYP-9480 (C). The capacity brackets are 8.74 and 31.2, from
  the S19 audit's evaluation of Session206's bracket.
  (C): replaced by HYP-9480 (A). Its pre-registered P2 measured 0.498538 on [2^27, 2^28), and Session206 proves the
  root image has at most (log 2) X/2 + O(sqrt X) elements below X.
  The pre-audit version of this hypothesis conjectured the opposite of (A): cofiniteness of the order-7 bank. That
  ignored HYP-9480 (C) and was withdrawn after audit M2; see MISTAKE-749.
session: kind-pasteur-2026-11-09-S19
source: 05-knowledge/results/kps_s19_fixed_port_tails_20261109/README.md
related:
  - THM-5085
  - HYP-9480
  - THM-4888
  - HYP-9640
  - HYP-9530
  - THM-4940
  - THM-1370 (h-spectrum-omits-7-21-all-n)
---

# HYP-9655 - base-3 fixed-port banks

## (A) Finite base-3 banks miss a positive proportion

**Statement.** Let C be any finite set of core states (h0, R0, Delta0, h_u) of base-3 fixed-port tails (THM-5085).
Then the union of the forms

    h0 + (3/2) [ R0 N + Delta0 [N odd] - h_u Q(N) ],   N >= 0,

misses a positive proportion of the integers v = 1 (mod 3).

**Consequence.** No finite base-3 bank of fixed-port tails proves ternary Chan-Pak (HYP-9530 (A)). A proof needs
unboundedly many cores, or another mechanism. This parallels HYP-9480 (C) and its root-change program.

**Evidence.**
- The heuristic is HYP-9480 (C)'s Poisson argument, which carries over verbatim: the tail forms have the same shape,
  and only the factor x/2 changes.
- The order-<= 6 base-3 bank already misses 15 values below 10^6.
- The order-<= 7 bank misses only the permanent gaps 10 and 31 below 10^7. This is consistent with an exponentially
  small miss proportion.

**Falsifier.** A proof that some finite base-3 bank is cofinite.

**Note.** THM-5085 (4) gives the census: below 10^7 the order-<= 7 bank attains every v = 1 (mod 3) except 10 and 31.

## (B) Base-2 bank data for HYP-9480 (C)

- **Order <= 6** (273 core states, capacity bracket 8.74): 169 odd values below 10^7 are missed, other than 7 and 21.
  Their density per odd integer is stable at about 3.3*10^-5 on [10^5, 10^7] (S19 audit `out/b2_le6_density_1e7.out`).
- **Order <= 7** (2,831 core states, capacity bracket 31.2): no miss below 10^7 other than 7 and 21.

Both are consistent with HYP-9480 (C)'s positive miss proportion of Poisson order exp(-capacity).

## (C) The one-apex share

This is HYP-9480 (A)'s root case. The S19 table of the share per dyadic window (exact, computed from all W < 2^22) is:

| k | 6 | 8 | 10 | 12 | 14 | 16 | 17 | 18 | 19 |
|---|---|---|---|---|---|---|---|---|---|
| share | 0.5625 | 0.5547 | 0.5391 | 0.5127 | 0.5023 | 0.5001 | 0.4976 | 0.4992 | 0.4966 |

This is consistent with HYP-9480 (A)'s conjectured limit of 1/2, which was measured at 0.498538 on [2^27, 2^28).
