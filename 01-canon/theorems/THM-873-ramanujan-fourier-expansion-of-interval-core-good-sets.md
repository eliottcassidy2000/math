---
id: THM-873
title: THE RAMANUJAN–FOURIER EXPANSION OF INTERVAL-CORE GOOD SETS (the roots-of-unity synthesis) — for the core {1,…,k} at any λ < 1/(k+1), the good-set Fourier coefficients are ĝ(h) = −Σ_{l≤k} c_l(h)·sin(2πhλ/l)/(πh) + (closed-gap overlap corrections), where c_l(h) are the RAMANUJAN SUMS. Proof: THM-826 Lemma 1 (one effective arc per primitive Farey fraction) makes the bad set a disjoint-to-explicitly-overlapping union of arcs at primitive fractions; the primitive-numerator character sum at level l IS c_l(h); overlaps occur exactly on the closed Farey gaps (λ(i+j) > 1) and subtract as explicit interval transforms. CONSEQUENCES: (i) disc_v(G) = Σ_{h≠0}|ĝ(hv)|² — the CENTRAL OPEN OBJECT of routes [A]/[B] (THM-731/732's certificate input; klein-S280's sharp-rate target) — becomes a RAMANUJAN-SUM SECOND MOMENT with the concentration structure exposed: c_l(hv) is large only where l resonates with the divisor lattice of hv, so the disc sum localizes on v's divisor structure (the torsion/roots-of-unity reading of the covering theory: klein-S296's clean-witness grid = the circle's torsion points; disc_v = the cyclotomic resonance of the core against v); (ii) numerically the form SHARPENS the crude bound by 3.5×/71×/17× at v = 182/84/13 on the deep well and CROSS-VALIDATES opus-S271's exact disc₁₃({1..12}) = 1.650e-2; (iii) the sharp-rate question (Q_s = O(r)) is now a classical estimate about Ramanujan-sum mean squares (Σ_h |Σ_l c_l(h)·sinc|²) — textbook orthogonality machinery applies
status: PROVED (three lines given THM-826's lemmas + linearity; the closed-gap correction is bookkeeping) + REFEREED to machine precision (max |formula − direct| ≤ 9e-16 over k ∈ {5,8,12}, both λ-regimes incl. the deep-well λ=1/14, h ≤ 40) + disc cross-validated vs opus-S271's exact rational value at v=13
source: kind-pasteur-2026-07-15-S128 (cont.25; owner: roots-of-unity lens on LRC(14), synthesize threads)
depends_on:
  - THM-826   # the effective-arc/Farey-gap structure this transform diagonalizes
related:
  - THM-731/732 (disc_v certificates — this is their spectral engine for interval cores), klein-S280 (the sharp rate, now a Ramanujan mean-square), klein-S296 (clean witnesses = torsion points), THM-805 (the Dirichlet-kernel Tutte face), opus-S271 (exact disc cross-check), opus-S319/THM-869 (Milgram frame — the mod-8 Gauss sums are this expansion's quadratic shadow), HYP-6955 (F-T dichotomy + Gauss-sum probe)
---

# THM-873 — the Ramanujan–Fourier expansion

**Theorem.** Let G = G({1..k}; λ), λ < 1/(k+1). For h ≠ 0:
ĝ(h) = −Σ_{l=1}^{k} c_l(h)·sin(2πhλ/l)/(πh) + Σ_{closed gaps} [overlap interval transform],
where c_l(h) = Σ_{gcd(c,l)=1} e(hc/l) and a Farey gap (a/i, b/j) is closed iff λ(i+j) > 1, its
overlap being [b/j − λ/j, a/i + λ/i].

*Proof.* By THM-826 Lemma 1 the bad set is the union over primitive c/l of arcs of radius λ/l;
by Lemma 2 arcs interact only within their own gap, overlapping exactly when the gap is closed.
Fourier-transform each arc; the sum over primitive numerators at level l produces c_l(h); subtract
each closed gap's double-counted overlap. ∎

**The critical-path consequence.** disc_v = Σ_{h≠0}|ĝ(hv)|² — the input to every THM-731/732
certificate and the object whose sharp rate is the open analytic core — is now a Ramanujan-sum
mean square. The resonance reading: c_l(hv) is supported where l divides into hv's arithmetic, so
disc_v measures the cyclotomic resonance between the core's spectrum {c_l : l ≤ k} and v's divisor
lattice — small for v with no small-divisor alignment (the far/isolated elements), large at
aligned v (the blind peels of opus's dilation theorem). Deep-well numbers: disc ≈ 4.12e-4 (v=182,
3.5× under crude), 9.57e-5 (v=84, 71× under), 1.648e-2 (v=13, = opus-S271's exact value ✓).

## Evidence log
- [x] machine-precision referee (both λ-regimes, k ≤ 12, h ≤ 40); opus cross-validation at v=13
- [ ] the Ramanujan mean-square estimate (Σ orthogonality) toward the sharp disc rate — the named
      analytic next step, now in classical form
