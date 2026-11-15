# Wave 25 pre-registration (mac-mini, registered before the coverage and order-10 programs reported)

Registered: 2026-11-15T01:02:30Z

Setting (chanpak_absolute_reflection_20261114): cyclic join T = A -> B -> c -> A, with
H(T) = (f_A - h_A)(f_B - h_B) + 2 h_A h_B, where h = P(1) = H and f = P(2) = F_2. Cores: every tournament class
of order <= 9 (gentourng census, 52,095 distinct profiles). Padding = ordinal sum with a transitive tournament,
which sends (h, f) to (h, 2^k f).

Already measured at registration time:
- the congruences f = (-1)^n h mod 3, f = 0 mod 4, and f = 0 mod 8 at even order hold for all 52,095 profiles;
- band occupancy (admissible points within half a residual sd of the ridge and half an sd of the centre) is 0.24,
  0.34 and 0.30 at orders 7, 8 and 9.

- **P66.** Cyclic joins of two cores of order <= 9 (padding allowed) miss no odd N in (1000, 2*10^5].
- **P67.** Every prime in [10^5, 2*10^5] has at least 5 such representations (ordered pairs (A, B)).
- **P68.** At order 10 (the full gentourng census, 9,733,056 classes), the band(0.5, 0.5) occupancy of admissible
  points lies in [0.40, 0.75], above order 9's 0.302, and no profile violates the congruences.
- **P69.** Every odd N in (1000, 2*10^5] has a representation in which the line for B crosses the ridge band of
  B's order (|f_B - ridge(h_B)| <= 2 residual sd).

## Addendum (registered 2026-11-15T01:05:58Z, after the order-<=9 coverage run REFUTED P66, P67 and P69, before the order-10
## census finished)

Order-<=9 result being tested against: 10,330 odd N <= 2*10^5 missed. By window, the misses are 1507, 1269,
206, 1775, 3676 on [2^13, 2^14), [2^14, 2^15), [2^15, 2^16), [2^16, 2^17), [2^17, 2*10^5]. N = 0 mod 3 is
never missed from 2^13 on. Diagnosis: N != 0 mod 3 needs a core of even order, and order 8 is poor.

- **P70.** With every order-10 profile added (padding allowed), the misses in [2^16, 2*10^5] fall by at least a
  factor 10 (from 5451). The singleton reader then puts N = f + h through the order-10 ridge band at about 30
  admissible hits per line.
- **P71.** The misses below 2^16 change by less than a factor 2: order 10 does not reach them through the bulk.

## Addendum 2 (registered 2026-11-15T01:28:17Z, after P68 (0.714), P70 and P71 confirmed with order 10; before any run beyond 2*10^5)

- **P72.** Cyclic joins with cores of order <= 10 miss no odd N in (36,833, 10^6], and miss fewer than 1% of the
  odd numbers of [10^6, 2*10^6].

## Addendum 3 (registered 2026-11-15T01:30:58Z, before the sampling program was written)

Model: occupancy = 1 - exp(-lambda/k), where lambda = (classes whose profile lies in the band)/(admissible band
points), and k is a collision factor fitted on the exact order-10 census (band(0.5, 0.5) occupancy 0.714).

- **P73.** With lambda estimated from uniform random tournaments of orders 11 and 12, and with the class counts
  903,753,248 and 1,523,198,812,064, the predicted band(0.5, 0.5) occupancy is at least 0.95 at order 11 and at
  least 0.999 at order 12.

Correction (2026-11-15T01:33:32Z, before the fit was computed): the order-12 class count in P73 was misquoted. A000568(12) is
154,108,311,168 (the repository's copy of A000568), not 1,523,198,812,064. P73 is evaluated with the correct
count; its thresholds are unchanged.
