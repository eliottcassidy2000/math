---
id: HYP-9294
title: "Singular support of Vir_{p,q} beyond the first correction: (a) in polynomial degree N = (p-1)(q-1)/2 the minimal differential-ideal generators of grI sit exactly at Li degrees 0, 3, 5, 7, ..., 2p-3 (one each); (b) ker(J_infty R_V -> gr_F Vir_{p,q}) is finitely generated as a differential ideal for every minimal model, e.g. for Vir_{3,5} by b_4 (weight 11) and one degree-5 relation of weight 19"
status: >
  OPEN (conjecture).
  Evidence FINITE-EXACT (mod-p ranks, one prime unless stated). Every dim L checked against Rocha-Caridi.
  (a) degree-N generators observed:
    - (3,4) (3,5) (3,7) (3,8) (3,10) (3,11): Li {0, 3};
    - (4,5) (4,7): Li {0, 3, 5};
    - (5,6): Li {0, 3, 5, 7};
    - (2,5) (2,7) (2,9): Li {0}.
  (b) Generators (weight, degree) found:
    - Vir_{3,4}: (6,3), (9,3), then none to weight 38 (AEH 2020 prove this case);
    - Vir_{3,5}: (8,4), (11,4), (19,5), then none to weight 38;
    - Vir_{3,7}: (12,6), (15,6), (21,7), (23,7), (31,8) to weight 32;
    - Vir_{3,8}: (14,7), (17,7), (23,8), (25,8), (31,9) to weight 32;
    - Vir_{4,5}: (12,6), (15,6), (17,6), (23,7), (25,7) to weight 32;
    - Vir_{5,6}: (20,10), (23,10), (25,10), (27,10), (31,11), (32,11) to weight 32.
  PRE-REGISTERED TESTS of (a), all PASSED; predictions written in the source note before the runs:
    - (5,7): weights 24 27 29 31;
    - (5,8): 28 31 33 35;
    - (4,9): 24 27 29;
    - (6,7): 30 33 35 37 39;
    - (5,9): 32 35 37 39.
  The degree-N generators appear exactly there, with no other degree-N generator up to weight 34-40.
  Next layers observed: (5,7) degree 13 at 35 (x2) and 36; (4,9) degree 13 at 33.
  The weight-(2N+3) generator is THM-4671 (PROVED). Part (b) for Vir_{3,4} is AEH's theorem; Vir_{3,5} is open in the
  literature (AEH remark the weight-11 relation alone fails from weight 19, which the engine reproduces).
session: kind-pasteur-2026-10-18-S1
source: 05-knowledge/results/singular_support_arc_spaces_classical_freeness_20261018.md
scripts:
  - 04-computation/voa_singular_support/vir_ss.py (outputs in 04-computation/voa_singular_support/out/)
related: [THM-4671, HYP-9293, Andrews-van Ekeren-Heluani 2020, Feigin-Jimbo-Miwa-Mukhin-Takeyama 2005 (monomial bases)]
---

# HYP-9294 — the degree-`N` layer and finite generation for `Vir_{p,q}`

## (a) The degree-`N` layer

The Li-degree-`l` part of degree `N` in weight `2N + l` is spanned by monomials `x_0^(N-s) x_(j_1) ... x_(j_s)` with
`j_1 + ... + j_s = l`.

**Conjecture.** The minimal generators of `grI` of polynomial degree exactly `N` are:
- `x_0^N` in Li degree 0;
- one generator in each odd Li degree `3, 5, ..., 2p-3`.

So there are `p - 1` of them, and the last sits at weight `2N + 2p - 3`.

**Heuristic.**
- The weight `2N + 2p - 1` is where the image of `M(c, 2N)` in `V^c` acquires its first singular vector (the
  Rocha-Caridi term `-q^(2N+2p-1)`), and the layer stops just below it.
- For `p = 2` the layer is `{0}`, which is classical freeness.
- THM-4671 is the first rung (Li 3), and its coefficient `A` carries the factor `(p - 2)`.
- One expects the rung at Li `2j+1` to carry a coefficient with factors `(p - 2)(p - 3)...(p - j - 1)`. This prediction
  is untested.

## (b) Finite generation

The Arakawa-Linshaw question asks whether `ker pi` is always finitely generated as a differential ideal (for `W_3`
and in general). AEH answered it for Ising. The data show:
- a finite list of generators up to the computed weight in every case;
- generators in degree `N + 1` at Li degrees 7 and 9 once `q >= 7` (for `p = 3`);
- a single degree-`N+1` generator (Li 9) for `(3,5)`.

**Cheapest tests.**
- Push `Vir_{3,5}` and `Vir_{3,4}` to weight 38-40 (running).
- Compute the bigraded Hilbert series of `gr_F Vir_{3,5}` and fit an AEH-type fermionic form.
