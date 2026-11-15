---
id: HYP-9692
title: "Band fullness: for every order n >= 12, every admissible point (h odd, f = 0 mod 4 (mod 8 for even n), f = (-1)^n h mod 3) of the ridge band of the (H, F_2) profile cloud of n-vertex tournaments is an actual profile. With the census below 2^32 and cyclic joins with small readers above it, this would close Chan-Pak; occupancy 0.30 at order 9, 0.71 at order 10, predicted 0.89 at order 11 and full at order 12"
status: >
  OPEN (CONJECTURE), mac-mini-2026-11-14-hosts, wave 25 of chanpak_waves_20261113.
  FINITE-EXACT:
  - the profile census of every class of orders <= 10 (gentourng; 9,733,056 classes and 1,550,357 distinct profiles
    at order 10);
  - band(0.5, 0.5) occupancy 0.236, 0.337, 0.302, 0.714 at orders 7, 8, 9, 10;
  - coverage: cyclic joins of two cores of order <= 10 attain every odd N in (36,833, 2,398,737]; with order <= 9
    alone they miss 10,330 odd N <= 2*10^5;
  - certificate: a 12-vertex cyclic join TT_2 -> B -> c -> TT_2 has exactly 100,003 Hamiltonian paths.
  NUMERICAL (model): a Poisson occupancy model fitted on order 10 (collision factor 4.60) reproduces order 10
  from samples (0.713) and predicts 0.892 at order 11 and 1.0000 at order 12 (lambda = 112).
  ARGUMENT, not a proof: the reduction below assumes the band geometry continues its measured trend beyond order 12.
  Pre-registered P66, P67, P69 REFUTED; P68, P70, P71, P72 CONFIRMED; P73 half confirmed (order 12) and half
  refuted (order 11: 0.892 < 0.95).
session: mac-mini-2026-11-14-hosts
date: 2026-11-14
source: 05-knowledge/results/chanpak_waves_20261113/README.md#wave-25-hosts-and-cores--no-heisenberg-orbit-the-cyclic-join-in-zsqrt-2-and-band-fullness
scripts:
  - 05-knowledge/results/chanpak_waves_20261113/scripts/prof.c, profhf.c (census), band.py, occupancy.py,
    fullness_model.py, sample.c, windows.py, cover.py, cert_100003.py
related:
  - THM-5172 (congruences, Z[sqrt(-2)] form, no Heisenberg orbit), chanpak_absolute_reflection_20261114 (the host),
    HYP-9682 (the invertible orbit), THM-4094 (census to 2^32)
---

# HYP-9692 - band fullness closes the two-core host

**Setting.** A tournament T of order n has the profile (h, f) = (H(T), F_2(T)), the numbers of ordered 1-slot and
2-slot path covers. The cyclic join A -> B -> c -> A has count

    N = (f_A - h_A)(f_B - h_B) + 2 h_A h_B,     equivalently     3N = Re(alpha_A conj(alpha_B))

in Z[sqrt(-2)], with alpha = (f - 3h) + f sqrt(-2) (THM-5172 (C)). Fixing N and the first core A puts the second
core's profile on the bounded line a f_B - c h_B = N, with a = f_A - h_A and c = f_A - 3h_A.

**Ridge band.** The (h, f) cloud of order n is a thin band:

| order | ridge | residual sd | sd of f given h / sd of f |
|---|---|---:|---:|
| 9 | f = 6115 + 13.68 h | 766 | 0.085 |
| 10 | f = 30574 + 15.69 h | 3545 | 0.072 |

B_n(kappa) denotes the admissible points (THM-5172 (B)) with |h - m_n| <= kappa sd_n and
|f - ridge_n(h)| <= kappa s_n.

**Conjecture.** For every n >= 12, every point of B_n(1) is the profile of an n-vertex tournament.

**Why it would close Chan-Pak (argument).**
1. A reader of scale a, with (a, c), c/a in [1/3, 1), crosses the band in an h-interval of length at least
   2 s_n/beta_n. Admissible points on its line recur at most every 24a in h: one condition modulo a for f to be an
   integer, then parity, f mod 4 or 8, and f mod 3. So for a <= a*(n) = s_n/(12 beta_n), fullness gives a profile on
   the line.
2. Readers TT_k (a = 2^k - 1) step by about a factor 2. Each order's window spans a factor U_n/L_n of about 2.1.
   So the readers alone, and certainly together with a few small cores, cover the window [L_n, a*(n) U_n].
   Measured (scripts/windows.py):

   | order | a*(n) | window covered |
   |---|---:|---|
   | 10 | 18.8 | [8.9·10^4, 3.6·10^6] |
   | 11 | 85.2 | [5.4·10^5, 9.8·10^7] |
   | 12 | 428.8 | [3.6·10^6, 3.2·10^9] |

   Consecutive windows overlap with growing margin: a*(n) grows about fivefold per order, while L_n grows about
   sixfold. So each N lies in two consecutive windows, one of even order, which handles the mod-3 condition of
   THM-5172 (C).
3. Every odd N < 2^32 is already a count (THM-4094 census).

So band fullness from order 12 on, together with the persistence of the measured geometry, would prove Chan–Pak.
A one-dimensional local limit law for H alone cannot do this: the bulks of consecutive orders are separated by
factors of about (n + 1)/2. The readers' scales bridge those gaps; that is what the second count buys.

**Data.**
- *Occupancy* of B_n(0.5): 0.236, 0.337, 0.302, 0.714 at orders 7–10. The Poisson model predicts 0.892 at
  order 11 and 1.0000 at order 12, with about 0.005 empty points expected in all of B_12(0.5).
- *Coverage* actually achieved with the census of orders <= 10: every odd N in (36,833, 2,398,737] is a
  cyclic-join count. Primes in [10^5, 2*10^5] have at least 24 representations, mean 48. The predicted order-10
  window is [8.9·10^4, 3.6·10^6] at full occupancy.
- *Certificate:* the prime 100,003 picks the reader TT_2. The line 3 f_B - h_B = 100,003 then contains the
  order-9 profile (1997, 34000), and the 12-vertex join has exactly 100,003 Hamiltonian paths (outgoing masks
  [2046, 2044, 3192, 2800, 4064, 3520, 2944, 3844, 2572, 3108, 2376, 3]).

**What a proof needs.** A two-dimensional local limit theorem for (H, F_2) of random tournaments, at scale 1 on the
admissible sublattice (index 24 in Z^2, 48 at even order), in the ridge band. That is the "remaining task" of the owner's
formulation, stated as a property of random tournaments.
