# Audit2: the required edits of sierpinski_audit2_20261022.md as data, with a check that every old string occurs exactly
# once in its file (nothing is modified).  Run from the repo root or from this directory.
import os, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../..'))
THM = '01-canon/theorems/THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md'
HYP = '05-knowledge/hypotheses/HYP-9365-two-loop-law-for-hamiltonian-paths-in-regular-tournaments.md'
RES = '05-knowledge/results/sierpinski_steiner_collatz_tournaments_20261022.md'
EDITS = [
 ('T1', THM,
  "  PROVED, pre-audit (mac-mini-2026-10-22-sierp): (1)-(3) are elementary identities, (3) as a formal power series in z,",
  "  PROVED (mac-mini-2026-10-22-sierp; independently audited 2026-10-22, 05-knowledge/results/sierpinski_audit2_20261022.md:\n"
  "  (1)-(4) confirmed by exact rational brute force; one side remark corrected): (1)-(3) are elementary identities, (3) as a formal power series in z,"),
 ('T2', THM,
  "  the closed formula's z^4 and z^6 coefficients equal the brute-force excess-1 cluster sums for R_15 and P_19, and the\n"
  "  brute-force excess-0 sums equal the log rho coefficients for L <= 6 (cluster_check.py, two_loop_series_check.py).",
  "  the closed formula's z^4 and z^6 coefficients agree with the brute-force excess-1 cluster sums for R_15 and P_19 to the\n"
  "  4 digits printed (vertex-transitive cases; the dumbbell has no term below z^14 for any regular T), and the brute-force\n"
  "  excess-0 sums equal the log rho coefficients for L <= 6 (cluster_check.py, two_loop_series_check.py). AUDIT\n"
  "  (sierpinski_ocf_20261022/audit2/): exact rational equality of the excess-0 and excess-1 cluster sums with (3) and (4),\n"
  "  by set-partition enumeration and by the blow-up scaling log Psi(M (x) J_m/m) = sum_e m^(-e) X_e(M), for all regular\n"
  "  tournaments with n = 7 (to z^12) and n = 9 (to z^10), R_11, P_11 and five non-vertex-transitive n = 11 classes (to\n"
  "  z^8), R_15 and P_19 (to z^8), the non-circulant regular 7-tournament to z^16 (dumbbell nonzero at z^14, z^16) and\n"
  "  Eulerian skew matrices (one to z^20; dumbbell from z^6); the kernel table equals the brute-force coefficient of every core class\n"
  "  of total length <= 14."),
 ('T3', THM,
  "  walks; Wormald's generating functions use (I - AD)^-1 for blow-ups.",
  "  walks; Wormald's generating functions use (I - AD)^-1 for blow-ups. Wormald's Lemma 3 (a fixed linear forest with a arcs\n"
  "  and s nontrivial paths lies in a uniformly random regular tournament with probability ~ 2^(-a) e^((a-s)/n), by McKay's\n"
  "  method) is the averaged counterpart of the forest identity (1)."),
 ('T4', THM,
  "      log Φ_T(z) = Σ_(e ≥ 0) X_e(z),",
  "      log Φ_T(z) = Σ_(e ≥ −1) X_e(z),"),
 ('T5', THM,
  "of **excess** `e = #arcs − #blocks`, i.e. cyclomatic number `e + 1`.",
  "of **excess** `e = #arcs − #blocks`, i.e. cyclomatic number `e + 1`. Trees (`e = −1`) occur exactly when `T` is not regular (`audit2/b_validate.out`)."),
 ('T6', THM,
  "* The `z⁴` and `z⁶` coefficients match brute force: −0.19911 against −0.19910 and −0.060250 against −0.060250 for `R_15`; −0.16271 against −0.16270 and −0.0068139 against −0.0068140 for `P_19` (`two_loop_series_check.py`).",
  "* The `z⁴`, `z⁶`, `z⁸` coefficients equal the brute-force cluster sums exactly: −0.1991111, −0.0602495, −0.0526505 for `R_15`; −0.1627059, −0.0068139, −0.0012067 for `P_19` (`audit2/c_session_cases.out`; `two_loop_series_check.py` compared 4 printed digits). Non-vertex-transitive tournaments and the dumbbell: `audit2/c_two_loop_check.out`, `audit2/c_dumbbell_n7_L16.out`."),
 ('T7', THM,
  "* For vertex-transitive `T`, `d` is constant and `G·1 = 0`, so the dumbbell vanishes. The figure-eight becomes",
  "* For vertex-transitive `T`, `d` is constant and `G·1 = 0`, so the dumbbell vanishes. For every regular `T` the first term of `d` is the constant `2 diag(Ŝ²)`, so the dumbbell has no `z^L` term with `L < 14` (first term: loops `a = b = 6`, bridge `c = 2`, i.e. `(a, b; c) = (6, 6; 2)`). The figure-eight becomes"),
 ('T8', THM,
  "a function of the spectrum alone. The theta becomes `(1/n) Σ_(k_1+k_2+k_3 ≡ 0)`",
  "a function of the spectrum alone. For circulant `T` the theta becomes `(1/n) Σ_(k_1+k_2+k_3 ≡ 0)`"),
 ('T9', THM,
  "The top term is `N_(n−1) = Σ_σ (−1)^(#backward steps) = (−1)^(n−1) F(T, −1)`, the alternating value of the forward-step polynomial of THM-4756.",
  "The top term is `N_(n−1) = Σ_σ (−1)^(#backward steps) = (−1)^(n−1) F(T, −1)`, the alternating value of the forward-step polynomial of THM-4756.\n\n"
  "Reversing one odd path (say the one through the least vertex lying on an odd path) is a sign-reversing involution, so only forests whose paths all have even length contribute to `N_L`; in particular `N_L = 0` for odd `L` (audit, `audit2/a_forest_identity.out`). The same involution on a single run kills every cluster that contains an odd run."),
 ('T10', THM,
  "That is where, for example, Paley character sums stop collapsing to polynomials in `q`.",
  "For doubly regular tournaments, Paley included, the three-loop terms (`K_4` among them) still reduce to functions of `n`, because they are fixed by the 4-vertex census (HYP-9365 (d)); Paley character sums do not stop collapsing at three loops."),
 ('H1', HYP,
  "  OPEN, pre-audit (mac-mini-2026-10-22-sierp). The expansion is derived (THM-4757: exact forest identity, proved one- and",
  "  OPEN (mac-mini-2026-10-22-sierp; independently audited 2026-10-22, 05-knowledge/results/sierpinski_audit2_20261022.md:\n"
  "  every exact H, y and y_pred recomputed and confirmed; wording and credits corrected). The expansion is derived (THM-4757: exact forest identity, proved one- and"),
 ('H2', HYP,
  "* Since `Σ_k x_k² = (n−1)/n` for every regular tournament, `ρ = e(1 − 1/n + O(n^(−2)))`. Hence",
  "* Since `Σ_k x_k² = (n−1)/n` for every tournament, and here `Σ_k x_k^(2m) = (n−1)/n^m`, `ρ = ((n+1)/(n−1))^((n−1)/2) = e(1 − 1/n + O(n^(−2)))`. Hence"),
 ('H3', HYP,
  "* This sharpens HYP-9274 (A) with an explicit rate. The `1/n` coefficient",
  "* This sharpens HYP-9274 (A) with an explicit rate. The constant `e` itself is the classical random-regular value: a uniformly random regular tournament has `E[H] ∼ e·n!/2^(n−1)` (Wormald, arXiv:2610.16245: the case `d = 1` of Theorems 1 and 5, via Lemma 3 and McKay's method; Alon's observation that the Adler–Alon–Ross bound holds for regular tournaments, recorded by Yuster, arXiv:1511.07784). New here are the single-tournament statement and the rate `3/(2n)`. The `1/n` coefficient"),
 ('H4', HYP,
  "* Data: `y(P_q) = 0.398, 0.417, 0.431, 0.440` at `q = 19, 23, 27, 31`; the two-point extrapolation gives 0.502.",
  "* Data: `y(P_q) = 0.398, 0.417, 0.431, 0.440` at `q = 19, 23, 27, 31`. Polynomial extrapolation in `1/q` gives 0.502, 0.479, 0.425 from the last 2, 3, 4 points, so it fixes the limit only loosely. For every doubly regular `T` the two-loop value is exactly `y_pred = (3n² − 8n + 1)/(6(n − 1)²) = 1/2 − 1/(3n) − 1/n² + O(n^(−3))`, and the residual `q(y − y_pred) = −1.55, −1.53, −1.50, −1.49` is stable, as an `O(1/q)` three-loop term should be; `κ(P_q)/(e(1 − 3/(2q))) − 1 = 0.0094, 0.0063, 0.0045, 0.0034`, about `3.3/q²` (`audit2/f_families.out`; `H(P_27)` recomputed in `audit2/d_paley27.out`)."),
 ('H5', HYP,
  "* Exact data: `y(R_33) = 1.2106`, with extrapolations 1.3975 (3 points) and 1.406 (4 points).",
  "* Exact data: `y(R_33) = 1.2106`. Mixed polynomial fits in `1/n` are unstable (1.381, 1.3975, 1.406, 1.475, 1.968 from the last 2–6 points) because of the `n mod 4` term of (e); fitted within each residue class mod 4, the 4- to 7-point fits give 1.3971–1.4002 (`audit2/f_families.out`)."),
 ('H6', HYP,
  "* Together with `ρ(R_n) → cosh 1/cos 1`, this settles the shape of the approach to Wormald's constant: `κ(R_n) = (cosh 1/cos 1)(1 − c/n + …)`.",
  "* Since `ρ(R_n) = (cosh 1/cos 1)(1 − 1/n + O(n^(−2)))` (Bhattacharjee's expansion of the regular-block gain `e^(1/n) ρ(R_n)` has no `1/n` term), the conjecture predicts `κ(R_n) = (cosh 1/cos 1)(1 − (1 + y_∞)/n + O(n^(−2)))`, with `1 + y_∞ = 2.398872`. This is conditional on HYP-9365; it does not yet settle the approach to Wormald's constant."),
 ('H7', HYP,
  "Here the theta term survives, `n·Θ → −0.23703`, because",
  "Here the theta term survives, `n·Θ → −0.237054`, because"),
 ('H8', HYP,
  "* Hence non-isomorphic doubly regular tournaments of one order should agree to relative `O(n^(−3))`.",
  "* Hence `X_0`, `X_1`, `X_2` are the same formal power series for all doubly regular tournaments of one order (a theorem, given THM-4757 (2)). That non-isomorphic ones then agree to relative `O(n^(−3))` in `H` is the conjectural part."),
 ('H9', HYP,
  "They have the same size at both orders, as an `O(n^(−3))` (four-loop) difference requires.",
  "They have the same size at both orders, which is consistent with an `O(n^(−3))` (four-loop) difference; two orders, with a sign change between them, cannot separate `n^(−3)` from `n^(−2)` or `n^(−4)`. At order 23, `D_23.2` and `D_23.3` have identical 5-vertex censuses (their 6-vertex censuses differ), and their `H` differ by only `2.8·10^(−8)` relative (`audit2/e_drt_noIE.out`, `audit2/e_drt_census6.out`)."),
 ('H10', HYP,
  "and `w_n κ_n / (F(R_n,−1)/n!) ≈ 1.62–1.72` for `n = 13–21` (`fneg1.out`).",
  "and `w_n κ_n / (F(R_n,−1)/n!) ≈ 1.62–1.69` for `n = 13–21` with the cubic smooth fit of `rot_limit_fit.py` (`fneg1.out`). The ratio depends on the fit: 1.60–1.70 with a quartic, 0.40–2.76 with a quadratic (`audit2/g_wobble.out`). The values and signs of `F(R_n, −1)` are FINITE-EXACT (`n ≤ 21`); the ratio is NUMEROLOGY."),
 ('H11', HYP,
  "| `R_n`, `n ≤ 33` | 1.2106 (33) | 1.2415 | 1.398872 | 1.3975 / 1.406 |",
  "| `R_n`, `n ≤ 33` | 1.2106 (33) | 1.2415 | 1.398872 | 1.3971–1.4002 (within residue classes mod 4) |"),
 ('H12', HYP,
  "**Conjecture.** `κ(T) = H(T) 2^(n−1)/n! = ρ(T) (1 − y(T)/n + O(n^(−2)))`, uniformly over regular tournaments.",
  "**Conjecture.** `κ(T) = H(T) 2^(n−1)/n! = ρ(T) (1 − y(T)/n + O(n^(−2)))`, uniformly over regular tournaments.\n\n"
  "*Notation below:* `y` is the exact residual `n(1 − κ/ρ)` and `y_pred` the two-loop value `y(T)` defined above; the conjecture says `y = y_pred + O(1/n)`."),
 ('R1', RES,
  "| rotational | `y → tan 1 tanh 1 + ½ tan² 1 − 1 = 1.398872`; exact-data extrapolation 1.3975–1.406 |",
  "| rotational | `y → tan 1 tanh 1 + ½ tan² 1 − 1 = 1.398872`; exact-data extrapolation within residue classes mod 4: 1.3971–1.4002 |"),
 ('R2', RES,
  "The scaled gap `n³(H/H(P) − 1)` is ≈ 0.2 at both orders 19 and 23, the `O(n^(−3))` four-loop scaling (`drt_results.out`).",
  "The scaled gap `n³(H/H(P) − 1)` is ≈ 0.2 at both orders 19 and 23, consistent with the `O(n^(−3))` four-loop scaling (two orders only; `drt_results.out`)."),
 ('R3', RES,
  "with ratio ≈ 1.65 for `n = 13–21`.",
  "with ratio ≈ 1.65 for `n = 13–21` (fit-dependent, NUMEROLOGY; `audit2/g_wobble.out`)."),
 ('R4', RES,
  "There the arithmetic of the tournament (Paley character sums, subtournament censuses) appears for the first time.",
  "There the 4-vertex subtournament census enters as a genuine 4-index sum; for doubly regular tournaments, Paley included, that census is fixed by `n`, so their individual arithmetic (5- and 6-vertex censuses) enters only at four loops."),
]

if __name__ == '__main__':
    ok = True
    for tag, f, old, new in EDITS:
        txt = open(os.path.join(ROOT, f), encoding='utf-8').read()
        c = txt.count(old)
        status = 'ok (occurs once)' if c == 1 else f'PROBLEM: occurs {c} times'
        if new in txt: status += '; NOTE: new text already present'
        ok &= (c == 1)
        print(f"{tag:4s} {f.split('/')[-1][:60]:60s} {status}")
    print('ALL OLD STRINGS UNIQUE' if ok else 'SOME OLD STRINGS NOT UNIQUE')
    if len(sys.argv) > 1 and sys.argv[1] == '--markdown':
        for tag, f, old, new in EDITS:
            print(f"\n**{tag}** `{f}`\n\nold:\n```\n{old}\n```\nnew:\n```\n{new}\n```")
