---
id: HYP-9172
title: "(a) For every n >= 5 there is no S_n-equivariant map from switching classes of tournaments to even (untwisted) Euler graphs that is bijective on isomorphism types, and none from tournaments to even graphs. (b) The THM-479 branch counts N_odd(n), N_lev(n) are integers for n = 0 (mod 4) by a combinatorial pairing, as Theorem B does for n = 2 (mod 4)"
status: >
  OPEN.
  (a) PROVED for 5 <= n <= 100 (switching classes) and for 5 <= n <= 15 plus 16 sporadic n <= 86 (tournaments) by
  block-lemma certificates (THM-4531). FINITE-EXACT exhaustive matchings for n <= 10 (deficits 1, 1, 2, 3, 4, 10) and
  n <= 9 (deficits 1, 2, 5, 12, 44). The certificate method needs additive representations n = b_1 + b_2 (n even) or
  n = 2b + b' (n odd) over "twist-rigid" block sizes (1, prime powers q = 3 mod 4, q = 5 mod 8, 9). That is a
  Goldbach/Lemoine-type question, so (a) for all n follows from a weak form of those conjectures, or needs a different
  certificate. The growing deficits suggest it is true.
  (b) Integrality itself holds for n <= 16 (THM-479 closed form). At n = 8 the even-Aut classes mix Sylow types
  C_2 (19), C_4, V_4, C_8, D_8, so no pairing inside one Sylow type works. OPEN.
source: collatz-procgen-20260922 session, tbij lane (2026-10-01), conjecture C1 and the open n = 0 mod 4 case; promoted with THM-4531
related:
  - 01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md
  - 01-canon/theorems/THM-479-level-holonomy-and-branch-split.md
---

# HYP-9172 — no natural bijection for the twisted Mallows–Sloane identities at any n >= 5

**Progress on (a) (2026-10-02; independent audit owed).**
[`natural_matchings_switching_classes_20261002.md`](../results/natural_matchings_switching_classes_20261002.md),
written in parallel with THM-4531, proves the switching-class half of (a) for every 5 <= n <= 10^6. Lexicographic
products of rigid blocks are rigid (Lemma 5.3(d) there), so rigid blocks exist for every odd order with no prime
factor = 1 mod 8 (the set Sigma), which contains every twist-rigid size listed above. The forced classes are told apart
by hand, and the number A(n) of available constructions is >= 2 for every 5 <= n <= 10^6 except n = 8, 16, 40, which
have computer certificates. So the switching half of (a) reduces to A(n) >= 2: binary additive problems over Sigma,
a set of density about x/(log x)^(1/4), with one summand a prime = 3 mod 4 when n = 0 mod 4. The tournament half
gains two infinite families: n in Sigma with at least two distinct prime factors, and twice such n (21, 30, 35, 39,
42, 45, ... are new). The exhaustive deficits there agree with the ones above for n <= 9 and n <= 8.
