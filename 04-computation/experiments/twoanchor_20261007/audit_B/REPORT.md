# Audit B — tournaments as chirotopes, Legendre friezes, partner product law (THM-4602, HYP-9241; results note sections 8-9)

Independent adversarial audit, 2026-10-07, for session mac-mini-2026-10-07-twoanchor. The results note was audited at sha1 702c3d47. Saved by the main session from the auditor's final message.

## Verdicts

| Item | Verdict |
|---|---|
| THM-4602 (1): Pf ∈ {±1, ±3}, \|Pf\| = 3 iff exactly one 3-cycle | CONFIRMED on all 64 labelled 4-tournaments; proof sketch correct |
| THM-4602 (2): (a) ⇔ (b) ⇔ (c) | CONFIRMED exhaustively for n ≤ 9, with integer-vector realisation certificates |
| (2) "the locus of positive SL2 friezes" | CORRECTION |
| (3)(a)–(c) | CONFIRMED for the 17 primes p ≡ 3 mod 4, p ≤ 131, plus 360 random strips |
| (3)(d) | CONFIRMED but trivial; CORRECTION: the Singer strip closes as a frieze |
| (3)(e) | CONFIRMED as an upper bound; CORRECTION: "at most"; quiddity −2 gives a closed frieze through p points |
| "Paley = finite-field totally positive tournament" | OVERCLAIM |
| §8.3 "covering all p+1 points needs the PGL2 Singer cycle" | FALSE |
| §8.4 Catalan dictionary | CORRECTION: repeats the mechanism retracted in MISTAKE-060/061 |
| HYP-9241 numbers | CONFIRMED; 0 unequal-time merges |
| HYP-9241 "Karlin–McGregor / total positivity" | OVERCLAIM: the data show independent pair events |
| Prior art | (1), (2), (3)(a)–(b) are themselves KNOWN |
| §8.6 | CONFIRMED |
| §8.1 | algebra CONFIRMED; minor wording |

## Corrections (applied in the canon files)

1. **Positive friezes are a slice.** The totally positive part is the transitive pattern. Positive real friezes are its slice `p_(i,i+1) = p_(1n) = 1`; Conway–Coxeter friezes are the integer points of that slice. For odd n the slice ≅ Gr⁺(2,n)/T; for even n it meets only special torus orbits.
2. **(3)(d).** Any ordering of the p+1 points gives the class. The PGL2 Singer ordering moreover closes as a frieze of width p−2: `v_(p+1) = χ(det g)v_0 = −v_0`, with 2-periodic quiddity.
3. **(3)(e).** Elliptic constant quiddities cover at most (p+1)/2 points, and hyperbolic ones at most (p−1)/2. The parabolic quiddity −2 gives a closed frieze through p points. No constant quiddity covers all p+1 points.
4. **The Paley analogy.**
   * The real chart gives the transitive pattern, while the F_p chart gives Paley + sink: this is an ANALOGY.
   * A χ-positive configuration is, after an SL2 change of coordinates, a transitive subtournament of Paley + sink, so it has ≤ tt(QR_p) + 1 points. This bound is attained for p ≤ 43.
   * The proved bound is tt ≤ 2√p + 1. Logarithmic growth is unproved.
5. **Catalan.** THM-438's C_k is a signed Möbius sum over even-series patterns (MISTAKE-060/061), not a count of plane-tree tours. The match with frieze counts is a NUMERICAL COINCIDENCE, not a dictionary.
6. **HYP-9241: independence, not repulsion.**
   * `q_2 = q_01 q_12 q_02` within 3% for T ≥ 64; Karlin–McGregor would give π/4 for Brownian walkers.
   * `c_2 ≈ q_02/q_01 ≈ 1.2` is a lag effect; `q_03/q_01 ≈ 1.38`.
   * The R(R+1)/4 exponents are the vicious-walker ones, but no determinantal structure is exhibited.
   * `q_2/q_1^3` = 1.11–1.19 for T = 16..512; the 1.30 at T = 1024 is a fluctuation.
7. **KNOWN.**
   * (1) is the Babai–Cameron 2000 §3 / Knuth 1992 / Gunderson–Semeraro 2017 (Fact 18, Lemma 21) criterion.
   * (2) is Babai–Cameron 2000 Lemma 3.3 (with Brouwer 1980 and Lachlan 1984).
   * (3)(a)–(b) are Gunderson–Semeraro 2017 Def. 17 / Thm 19 and arXiv:2204.10775 Lemma 3.4, Prop. 3.5, going back to Paley 1933.
   * New to the repo: the Pfaffian / frieze-strip phrasing, the closure of Singer strips, the parabolic closed frieze, and the χ-positive bound.
8. **§8.1 wording.**
   * The tori share the identity.
   * The "mutation path" is an analogy.
   * The affine-group "Borel > torus" is an analogue of THM-4553's ladder, not the same group. Audit B notes that the ladder follows from (3)(b) together with Babai–Cameron Lemma 3.1.
9. **§8.6.** Quiver mutation at a source or sink is switching, so it preserves the tournament class.

## Confirmed numbers

* **Labelled counts.** LT = rank-2 = Pfaffian-±1 tournaments number 48, 384, 3840, 46080, 645120, 10321920 for n = 4..9, i.e. (n−1)!·2^(n−1).
* **Orbit sizes.** The maximum elliptic orbit is (p+1)/2; parabolic orbits have size p; hyperbolic orbits ≤ (p−1)/2. Checked for p = 7, 11, 19, 23.
* **tt(QR_p)** for the primes p ≡ 3 mod 4 from 3 to 199: 2, 3, 4, 5, 5, 7, 7, 7, 9, 8, 9, 9, 8, 11, 9, 11, 11, 11, 11, 11, 11, 12, 11, 11.
* **HYP-9241 runs.**
  * Run 1 (40,000 samples, 2048-bit y): `q_2/q_1^3` = 1.14–1.22.
  * Run 2 (20,000 samples, 4096-bit y): `q_2/q_1^3` = 1.11–1.18 and `q_3/q_1^6` = 1.66 → 1.91 over T = 16..512.
  * Slope ratios 1 : 2.91 : 5.6–5.8.
  * Unequal-time merges are impossible here: |m ln 3 − n ln 2| ≥ 4.4e-5 ≫ 2^(−1021).

## Sources

* Babai–Cameron, Electron. J. Combin. 7 (2000) R38.
* Gunderson–Semeraro, JCTB 126 (2017), arXiv:1509.03268.
* Gunderson–Semeraro, arXiv:2204.10775.

## Scripts

All in this directory:
* `chirotope_audit.py`
* `lt_enum.c`
* `legendre_audit.py`
* `km_audit.py`
* `km_audit_N40000_B2048_T1024_s20261007.out/.json`
* `km_audit_N20000_B4096_T2048_s777.out/.json`
* `km_audit_summary.py`
