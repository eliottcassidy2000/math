# Pre-registered predictions for the exact order-11 census (committed before the census finished)

**Session:** klein-2026-11-06. The order-11 census (`gentourng -q 11 r/8 | hspec 11`, eight parts, all 903,753,248
isomorphism classes) was launched at about 14:07 local time on 2026-11-06 and is expected to take about two hours. The
predictions below are read off the exact censuses of orders 3–10 in [`data/`](data/) (`census_<n>.txt`, produced by
[`scripts/hspec.c`](scripts/hspec.c)) and the recorded canon; each is checked by `spectrum_tables.py` / `level_tables.py`
on `census_11.txt` once it exists. Notation: `S(n)` = set of Hamiltonian-path counts of `n`-tournaments; `Ss(n)` = those of
strong `n`-tournaments; `n(h) = min{n : h ∈ S(n)}`; `m(h) = min{n : max H(n) ≥ h}`; level `k = c3` = number of cyclic
triangles; `Λ(m)` = the largest number of cyclic triangles of an `m`-tournament (`(m³ − m)/24` for odd `m`, `(m³ − 4m)/24`
for even `m`); `φ_n(k)` = least `h` at level `k` of order `n`; `f(n)` = Busch's strong floor; `B(h) = max{n : f(n) ≤ h}`.

| # | prediction | source of the expectation |
|---|---|---|
| P1 | `max H(11) = 95,095` and the census reads exactly 903,753,248 classes | OEIS A038375 (Royle 2022), A000568 |
| P2 | every one of the 208 odd values below 15,745 unattained at order 10 ([`data/unattained_10.txt`](data/unattained_10.txt)) is attained at order 11; hence `n(h) − m(h) ∈ {0, 1}` for every odd `h ≤ 15,745` | true at orders 6→7, 7→8, 8→9, 9→10 (all 2, 16, 9, 157 values fill at the next order) |
| P3 | every value first attained at order 11 is attained by a strong 11-tournament (`N(11) ⊆ Ss(11)`) | true at every order 3–10 |
| P4 | zone law: `φ_11(41) > 15,745`, i.e. every 11-tournament with more than `Λ(10) = 40` cyclic triangles has more Hamiltonian paths than any 10-tournament; equivalently no value of `S(10)` occurs at order 11 above level 40 | true at orders 5–10 with margins `φ_n(Λ(n−1)+1) / max H(n−1)` = 1.80, 1.53, 1.13, 1.39, 1.35, 1.47 |
| P5 | the level floor `φ_11(k)` is nondecreasing in `k`; `φ_11(k) = φ_10(k)` for `k ≤ 9` (exact: by Moon's bound a strong class at level `k` has order `≤ k + 2`); the weaker guess `φ_11(k) = φ_10(k)` for `10 ≤ k ≤ 17` (values 51, 69, 97, 107, 141, 171, 213, 227, 321) is recorded with the caveat that order 10 lowered the floors of levels 18–22 (e.g. 621 → 467 at level 20) | orders 5–10 |
| P6 | THM-4744's multiplicities `N(11, h)` = 157, 8, 121, 170, 85, 67, 278, 148, 117 for `h` = 45, 63, 75, 81, 105, 125, 135, 189, 223; exactly 3,536 strong 11-classes have `H ≤ 625`; the least strong value is `f(11) = 225`; the `z^11` coefficients of the recorded polynomials `P_d` (`d = 225..373`, `strong_ears_375.out`) equal the census counts | THM-4744 (orders 11, 12 computed by one-vertex ears) |
| P7 | every odd `h ≤ 80,830 = 0.85 · 95,095` other than 7 and 21 is attained at order 11 (HYP-9452 (1)); THM-4104's strong solid interval `[429, 80,265]` is confirmed | HYP-9452, THM-4104 |
| P8 | the unattained values below 95,095 lie at levels `≥ 53` (the top three of 0..55), and there are more than 1,000 of them: the top level consists of the 1,223 regular 11-tournaments (A096368) spread over thousands of odd values | the band lives in the top two or three levels at orders 7–10; the top level at odd orders is sparse (15 classes at order 9 gave 157 gaps) |
| P9 | for odd `h < 375` the strong window `W(h) = {n : h ∈ Ss(n)}` is an interval `[n(h), e(h)]` with `e(h) ∈ {B(h) − 1, B(h)}`, and `n_s(h) = n(h)` for every `h ∈ S(11)` | true for `h < 225` from orders ≤ 10 (90 windows end at `B(h)`, 19 one short) |
| P10 | the number of distinct values at order 11 lies in `[40,000, 48,000]`; the solid prefix `s(11)` (every odd `≤ s` except 7, 21 attained) is at least `0.85 · 95,095` | the sampled `h_spectrum_n11.out` saw 39,951 distinct values in 10⁶ samples; `s/max` = 0.92, 0.86, 0.93 at orders 8–10 |
| P11 | strong classes first appear at level `9 = 11 − 2` (Moon's bound), and the level widths (classes per level) sum to 903,753,248 | Moon; orders 3–10 |

**Addendum (14:32, census still running).** Two further laws read off orders 6–10 after the first commit:

| # | prediction | source of the expectation |
|---|---|---|
| P12 | top-level separation: `φ_11(55) > ψ_11(54)`, i.e. every regular 11-tournament (level 55) has more Hamiltonian paths than every non-regular 11-tournament; equivalently the values of the top level occur at no lower level | true at orders 6–10: `(ψ(Λ−1), φ(Λ))` = (37, 41), (159, 171), (613, 621), (3109, 3159), (14733, 14877); the margin at order 9 is only 1.6 %, and the 1,223 regular 11-tournaments are more varied than the 15 regular 9-tournaments, so this test is genuinely uncertain |
| P13 | the vertical gradient `h(T) / max_v h(T − v)` over all 11-tournaments stays below 11 (trivial, since `h(T) ≤ Σ_v h(T − v)`) and above 1 (insertion lemma), with maximum at least 5.35 (the order-9 maximum); the least ratio among classes whose value is new at order 11 lies in `[1.6, 1.8]` | 1.727, 1.690, 1.646, 1.681 at orders 6–9 — not tested at order 11 by this session's census (it records no deletions); recorded for a later session |
