# The mod-2 collapse of odd-symmetry cancellation polynomials is exact and lossy, and a Collatz strategy audit against the 2026 trusted advances

**kind-pasteur-2026-10-27-S6.** Owner dispatch: "for any odd-order symmetry action on `n` gaps, the cancellation
polynomial becomes `(1+t)^n` modulo two; that does not mean `n` independent toggles exist"; organize around integral
coordinates, operations that must be preserved, and positivity; look for ways to improve the Collatz strategies in the
light of the recent trusted advances. Companion ledger: [TRUSTED-ADVANCES-2026-10-27](../reference/TRUSTED-ADVANCES-2026-10-27.md).

**Status.**
- **PROVED (Lemma 1, with proof):** the owner's statement, in session 65's pair-orbit language, with the exact value
  of `n` and the exact loss. **FINITE-EXACT:** [`mod2_collapse_odd_symmetry_kps_20261027.py`](../../04-computation/experiments/mod2_collapse_odd_symmetry_kps_20261027.py),
  5 checks (session 65's product formula re-derived by direct enumeration on all permutations of `n <= 6` vertices;
  the collapse and its integral failure on all 1,852 odd-order non-identity permutations of `n <= 7` vertices).
- **AUDIT (section 3):** six strategy shapes from the T1/T2 advances mapped to the repository's named Collatz
  obstacles, each with preserved predicate, loss, and the cheapest decisive test. No Collatz theorem is claimed.
- **ANALOGY (typed):** the sparse-polynomial counting shape of the ETH manuscript against Collatz word counting.

## 0. Inheritance

- [Session 65, general signed pair-orbit defect polynomial](session65_cryptid_tournament_interfaces_20261026/session65-fifteens/GENERAL_PERMUTATION.md):
  `F_σ(z) = Π_O [(1+z)^{L_O} + ε_O (1-z)^{L_O}]`, minimum defect = number of even vertex cycles, invariant tournament iff
  odd order. Session 63 (`cancellation/RESULT.md`): orbit polynomials, coefficient characteristic and legal toggles are
  separate; a root multiplicity is not a symmetry-kernel dimension.
- [mac-mini, why the advances work II](../../07-reflections/why-the-advances-work-ii-potentials-symmetry-budgets-positive-coordinates-macmini-20261027.md):
  potential plus descent certificate; symmetry budgets; positive coordinates and coboundaries; the open Collatz
  boundary term. The owner's three organizing things are this reflection's sections 2, 1 and 2 respectively, and
  Lemma 1 is the smallest exact instance of all three at once.
- Collatz obstacles named in canon and used below: the Periodicity Conjecture's HARD class (HYP-9127, HYP-9134); THM-4469
  (no-divergence on an adjacent supercritical block pair is a Mahler Z-number statement); THM-4474 (bounded-lookahead
  provability window); THM-4476 (thin divergence); THM-4517 (no root-uniform basin density); THM-4581 (Haar
  coalescence); THM-4627 (Chowla for Collatz); the renewal/first-exit ledger and its "universal entry OPEN".

## 1. Lemma 1 (mod-2 collapse)

Let `σ` be a permutation of odd order on a vertex set, with pair-orbits of lengths `L_O` (all odd, all `ε_O = +1`). Write
`E_L(z) = Σ_{j even} C(L, j) z^j`, so that `F_σ(z) = 2^{q} Π_O E_{L_O}(z)` with `q` the number of pair-orbits.

1. For odd `L`, `E_L(z) ≡ (1+z)^{L-1} (mod 2)`.
2. Hence `F_σ(z) / 2^q ≡ (1+z)^n (mod 2)` with **`n = Σ_O (L_O - 1) = #pairs - #orbits`**: modulo two, the normalized
   cancellation polynomial is the polynomial of `n` independent toggles.
3. Integrally it never is (for `σ ≠ id`): `E_3 = 1 + 3z^2`, `E_5 = 1 + 10z^2 + 5z^4`, and the product `Π_O E_{L_O}` is
   an even polynomial, so it has no `z^1` term, while `(1+z)^n` has coefficient `n >= 1` at `z`. (The `z^2` coefficient
   alone does not separate them: for two 3-orbits it is `3 + 3 = 6 = C(4,2)`; the separation is the vanishing of every
   odd coefficient, i.e. the defect is always even, session 65's parity law.) The integral polynomial determines the
   multiset `{L_O}` (checked through `n = 7`); its mod-2 image determines only `n`.

*Proof.* (1) `(1+z)^L = (1+z)(1+z)^{L-1}` and `(1+z)^{L-1} ≡ (1+z^2)^{(L-1)/2}` is a polynomial in `z^2` modulo two, so the
even part of `(1+z)^L` modulo two is `(1+z)^{L-1}`. (2) Multiply over orbits. (3) Each `E_L` is even, so the product has
only even powers, whereas `(1+z)^n` has the odd term `n z`; they differ integrally for every `n >= 1`, i.e. for every
odd-order `σ ≠ id`. The orbit multiset is read off the integral polynomial by the check S5 (an exhaustive comparison;
a general proof would use the distinct factorizations of the `E_L`). `QED`

**The three coordinates.** *Integral coordinates:* the counts `C(L_O, j)`, which are positive and orbit-indexed; the
mod-2 image forgets them. *Preserved operation:* the orbit product, which survives the reduction (the mod-2 image is
still a product, of the factors `(1+z)^{L_O - 1}`). *Positivity:* the only thing that recovers the genuine model is the
positive integral polynomial, which still carries the orbit multiset. The "toggles" of the mod-2 image are not legal
operations on tournaments: a single pair cannot be flipped inside an orbit without breaking `σ`-invariance, which is
exactly session 63's separation of coefficient characteristic from legal toggles. The owner's statement is therefore
exact with `n = #pairs - #orbits`, and the warning is exact too.

## 2. What the 2026 advances have in common

The trusted items of the ledger share one shape with the mac-mini reflection:

| advance | the positive object | the finite compression | the certificate |
|---|---|---|---|
| zeta zeros `>= 2/3` (Alpöge–Furman / Lamzouri) | Weil's Hermitian form | a finite-rank compression; one Hilbert-space inequality | rank–trace inequality; Montgomery pair correlation (unconditional BGST) |
| local smoothing manuscript (T3) | a conserved current; a determinant energy | one frequency annulus; a thin cylinder | boundary flux of the current |
| forced Navier–Stokes blow-up (T2) | a self-similar profile's stability | a finite-dimensional spectral problem | computer-assisted plus Lean |
| unit-distance disproof (T1) | a number field with many units | an explicit finite point set | counting unit vectors algebraically |
| Erdős problems in Lean (T1) | — | — | a kernel-checked certificate |

In every case the content lives in a **positive object on a finite compression**, and the proof is a certificate on that
compression plus a lemma that the compression loses nothing that matters. That is the owner's triple.

## 3. Collatz strategy audit

For each shape: the Collatz obstacle it addresses (repository file), the map, what is preserved, what is lost, the
needed sidecar, and the cheapest decisive test. Every item is a proposal; nothing below is a result.

1. **Positivity on a finite compression (zeta shape).** Obstacle: no uniform potential exists on the integers (THM-4474's
   provability window; Conway's undecidability of generalized Collatz maps). Map: replace a potential on `Z` by a
   positive form on the finite transfer operator of the Syracuse map on residues modulo `3^J` (Tao's Syracuse random
   variables) or `2^K` (the repository's level 17–20 spectrum, `collatz_spectrum_equality_and_stitched_sheets_20261005`).
   Preserved: exact residue dynamics below the stopping time (THM-4627's periodicity). Lost: the integer size, i.e. the
   drift. Sidecar: the "compression loses nothing" lemma, which for Collatz is precisely the open boundary term of
   mac-mini II section 2. Cheapest test: certify with exact rational arithmetic that the second eigenvalue of the
   finite-section operator at levels `K = 10..20` is bounded away from 1 uniformly in `K`, using the inherited level data;
   if the gap decays with `K`, the shape is dead for Collatz and the note should say so.
2. **Kernel-checked exact lemmas (Erdős/zeta shape).** Obstacle: the Collatz canon has hundreds of exact elementary
   lemmas (carry cylinders (C3), the run law of THM-4793, the rank law, the cycle identities of THM-4515/4516) whose
   trust is FINITE-EXACT plus prose. Map: Lean 4 statements in the repository's existing lane (GMC2Main, TournamentH7).
   Preserved: everything. Lost: nothing but time. Cheapest test: formalize (C3) `w` is the first `n` carry bits of `z`
   iff `Q_z(w) ≡ 0 (mod 2^n)`, a one-induction statement; it is the lemma every Antihydra and root-27 argument rests on.
3. **Arithmetic extremizers (unit-distance shape).** Obstacle: Collatz hostiles found by random sampling miss structured
   extremizers (MISTAKE-101/102); the repository already found its best lines algebraically (the Mersenne and Fermat lines,
   THM-4556, THM-4793). Map: enumerate resonant families `x = (2^a - 3^b)/c` and `x ≡ ±1 (mod 2^k)` as the unit-distance
   disproof enumerated number-field points. Preserved: exact run structure. Lost: typicality. Cheapest test: for every
   `(a, b)` with `2^a - 3^b` small relative to `3^b`, compute the valuation word of the odd part and compare its first-exit
   charge with the paid-band records of `collatz_band_budget_hostile_20261012`; a resonant family that beats the random
   record `c* = 4.64` would be a new hostile, a family that cannot is evidence the records are typical.
4. **Conserved current with boundary flux (local-smoothing shape, T3).** Obstacle: the ROOT forcing cocycle equals the
   drift cocycle plus a coboundary (THM-4813 (A)), and the boundary term is uncontrolled. Map: the manuscript's current is
   an exact identity whose boundary flux is then estimated on a thin cylinder; the Collatz analogue is a current on the
   cylinder of fixed valuation words whose flux is the first-exit charge. Preserved: the identity. Lost: the estimate,
   which is the whole problem. Cheapest test: write the first-exit charge of `collatz_renewal_source_20261011` as the
   flux of an explicit current on cylinder boundaries and check numerically whether the flux is sign-definite on the
   paid bands; sign-definiteness would be the first positivity statement in this lane.
5. **Sparse-polynomial counting (ETH manuscript, T3) — ANALOGY.** The Collatz cycle equation `x(2^A - 3^j) = C_w` and
   Terras's coefficient-stopping sets are sparse Diophantine systems in the word `w`; the manuscript counts common zeros
   of sparse polynomial systems over `F_2`. Map: none exact (the Collatz systems are over `Z`, with the carry). Lost: the
   field. Cheapest test: none decisive; recorded so that the analogy is not re-proposed.
6. **Finite-compression blow-up certificates (Navier–Stokes shape, T2).** Obstacle: the repository's thin-divergence
   theorem (THM-4476) and Mahler reformulation (THM-4469) say a divergent orbit is a transcendence-type object. Map: the
   blow-up proofs certify a self-similar profile's stability on a finite spectral problem; a Collatz divergence would need
   a self-similar valuation word, which THM-4469 shows is a Mahler Z-number. Preserved: the self-similar ansatz. Lost:
   everything else; the analogy predicts no divergence, which is the conjecture. Recorded as a direction, not a plan.

**Ranking.** Items 2 and 1 are actionable this month (one is pure formalization, one is an exact computation on
inherited data with a yes/no outcome). Item 3 is a cheap hostile sweep. Items 4–6 are typed directions.

## 4. Reproduction

```bash
python3 04-computation/experiments/mod2_collapse_odd_symmetry_kps_20261027.py
```

Standard library; a few seconds; output beside the script.
