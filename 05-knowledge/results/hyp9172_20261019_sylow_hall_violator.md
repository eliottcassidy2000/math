# hyp9172_20261019 — HYP-9172 for every n: the Sylow Hall violator (mac-mini-2026-10-19-hyp9172)

Owner's prompt: prove HYP-9340 beyond Sylow (first open case Z_3 at n = 6), find a uniform rule or a definability obstruction for the prime-circulant bijection, and prove or disprove HYP-9172, thinking in odd primes and their powers.

Sibling note (same day): [even_graphs_20261019](even_graphs_20261019_sylow_refinement_components_walsh_duality.md), which lists the earlier 2026-10-19 notes.

Canon: [THM-4722](../../01-canon/theorems/THM-4722-no-equivariant-bijection-tournaments-even-graphs-for-every-n-sylow-hall-violator.md) (HYP-9172 (a): tournament half for every n >= 5, switching half for every odd n >= 5), [THM-4723](../../01-canon/theorems/THM-4723-sylow-profile-form-of-hyp-9340-and-the-defect-reduction-for-the-three-cycle-type.md), [HYP-9341](../hypotheses/HYP-9341-what-remains-open-after-thm-4722.md). Scripts in `04-computation/experiments/even_graphs_20261019/`.

## 1. The idea in one paragraph

THM-4531 reduces "natural bijection" (= S_n-equivariant, bijective on types) to Hall's condition for the stabiliser-containment graph between tournament types and even-graph types. THM-4721's Sylow refinement says the two sides have exactly as many types with a Sylow-3 symmetry, and exactly as many with a Sylow-5 symmetry. The symmetry groups of tournaments have odd order, and no odd-order group can carry both a Sylow 3- and a Sylow 5-subgroup of S_n (a 3-cycle and a 5-cycle through a common point generate A_5, A_6 or A_7: Jordan). The empty graph carries both. So the tournaments with either symmetry are one more than the even graphs they could be mapped to. This is the lane's n = 5 hand proof (Z_5 and the four Z_3's against ∅, K_{1,3}+K_1, K_{1,4}, K_{2,3}) turned into a uniform argument: "odd primes and their powers" enter exactly through the Sylow subgroups, which are the only p-subgroups for which the refined identity is proved.

| n | tau_T(P_3) = tau_E(P_3) | tau_T(P_5) = tau_E(P_5) | |S| | |N(S)| <= | THM-4531 exhaustive deficit |
|---|---|---|---|---|---|
| 5 | 4 | 1 | 5 | 4 | 1 |
| 6 | 1 | 2 | 3 | 2 | 2 |
| 7 | 4 | 4 | 8 | 7 | 5 |
| 8 | 16 | 12 | 28 | 27 | 12 |
| 9 | 1 | 48 | 49 | 48 | 44 |

## 2. Switching classes

For odd n, THM-4531 P1 identifies switching classes with the tournaments whose scores are all ≡ (n−1)/2 (mod 2), equivariantly, so class stabilisers have odd order. The restricted principle transfers because (i) for odd n, V = Z ⊕ B (Euler graphs ⊕ cuts) as S_n-modules, (ii) an invariant cut has no self-paired h-orbit, so the h-odd invariant graph of THM-4721's Lemma P has an h-odd Euler component (Lemma E'), and (iii) invariant cuts realise every invariant even-weight parity vector (Lemma Q). Then Sylow + Frattini give tau_C(P) = tau_EE(P) and the same Hall violator works: the switching half of HYP-9172 holds for every odd n >= 5. Numerically the refined identity also holds at n = 6 for every odd-order class (`switching_even_n.out`), but Lemma E' fails at n = 6 and n = 8 (36 pairs each; `euler_parity_lemma.out`), so even n needs another proof (THM-4531's certificates reach n = 100).

## 3. What this does not do

It does not produce a bijection (none equivariant exists), and it does not touch HYP-9340 beyond Sylow: the argument uses only Sylow subgroups precisely because those are the cases where multiplicity one and Frattini are available.
