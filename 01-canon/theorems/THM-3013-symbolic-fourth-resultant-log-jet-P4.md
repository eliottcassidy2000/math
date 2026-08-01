---
id: THM-3013
title: "The symbolic fourth resultant log-jet P_4(M,U,V) of the first-gap family"
status: VERIFIED-EXACT (interpolated under a verified degree ansatz; independently reproduced by a distinct algorithm) / AWAITING INDEPENDENT HOSTILE AUDIT
source: klein-S428 (two-agent build + adversarial rebuild; shape, digest and top-column checked here)
depends_on:
  - THM-3011-fourth-resultant-jet-and-the-third-edge-invoice
  - THM-3006-first-gap-wall-is-a-four-band-charge-density-with-all-order-multipole-law
related:
  - THM-2997-first-gap-wall-stripped-all-width-second-edge-circuit-positivity
  - THM-3000-fixed-edge-cumulant-curvature-universality-and-bounded-jet-transfer
  - THM-3003-antipodal-circuit-rigidity-and-the-multipole-spread-criterion
table: 05-knowledge/results/gmc_first_gap_fourth_resultant_jet_P4_table_thm3013.json
table_family: 05-knowledge/results/gmc_first_gap_fourth_resultant_jet_P4_table_family_thm3013.json
script: 04-computation/gmc_first_gap_fourth_resultant_jet_verify_thm3013.py
output: 05-knowledge/results/gmc_first_gap_fourth_resultant_jet_verify_thm3013.out
table_sha256: 200fae225af6c381733a386f31e107b81e6d33104699e5f69e8d7cd7e445d163
table_family_sha256: c9db12f0e92fa134102fbdb86c0d6dc2c66e925839d4930d872d897a40a7caa3
script_sha256: 2caec38a9104088b60ae4e42d23dd31690d3b2fb6a9ca338a4618e4c788f5678
output_sha256: 15d8653be373631c553320747d00988ea6904118bd3a05435c341ae559902730
hash_basis: LF-normalized bytes
---

# THM-3013 -- the symbolic P_4(M,U,V)

**VERIFIED-EXACT.  NOT symbolically derived: the object is an exact interpolation
inside a degree ansatz that was measured, then verified far outside its own box.**

This supplies the obligation left open by THM-3011 section 5 and by THM-3000
section 7: the fourth reduced resultant log-jet of the first-gap family, as an
explicit polynomial in `(M,U,V)` in the frozen THM-2997 table format.

## 1. The object

    P_4 = c_4 * D^8 * L_4,    D = U^2+3U-3V-1,    L_4 = ell_4(R_M)/4,

**717 terms**, degrees `(M,U,V) = (8,16,8)`, support exactly `{b+2c <= 16}`
crossed with `M^0..M^8` (all `81` `(b,c)` pairs occur, nothing outside), minus
`12` absences, all at the three corners `(16,0)`, `(0,8)`, `(0,0)`:

    (2,0,0) (2,0,8) (2,16,0) (6,0,0) (6,0,8) (6,16,0)
    (7,0,0) (7,0,8) (7,16,0) (8,0,0) (8,0,8) (8,16,0)

Tables, in the frozen `[m,u,v,num,den]` row format:

| file | `c_4` | content | sha256 |
|---|---|---|---|
| `..._P4_table_thm3013.json` | `1658880 = 2^12*3^4*5` | `1` (integer) | `200fae22...` |
| `..._P4_table_family_thm3013.json` | `61440 = 2^12*3*5` | `1/27`, denominators `{1,3,9,27}` | `c9db12f0...` |

## 2. `c_4` is CONVENTION-DEPENDENT -- do not canonize a single value

The normalization-free object is `D^8 * L_4`, with

    content(D^8 * L_4) = 1/1658880 = 1/(2^12 * 3^4 * 5).

At least three values are defensible and the frozen data cannot separate them:

- `c_4 = 1658880` -- integer coefficients, content `1`;
- `c_4 = 61440` -- content `1/27`, matching the *observed* pattern
  `content(P_j) = 3^-(j-1)` (the frozen `P_2`, `P_3` are **not** integral:
  contents `1/3`, `1/9`);
- `c_4 = 4096 = 2^12` -- the rule "`c_j` is the 2-part of `1/content`", which fits
  `c_1,c_2,c_3 = 1,16,128` equally well because `G_2, G_3` carry no prime past `3`.

The **sign** is also free: `c_1,c_2,c_3 = +1,+16,-128` shows no alternation.  Both
tables here take `c_4 > 0`.  Whoever extends the THM-2997 companion should pick
the convention deliberately and record it; `5` enters for the first time at `j=4`
because `L_j ~ (-1)^(j+1) 46/(j(j+1)) M^(j+1)` and `1/(4*5)` carries a `5`.

## 3. Controls (the reason this is trustworthy)

- **The pipeline re-derives the frozen table.**  Run at `j=1,2,3` it reproduces
  `P_1,P_2,P_3` **row by row** (counts `27/122/333`, every exact `Fraction`,
  including `-368/3` at `(3,8,0)` and `-4913/9` at `(6,10,1)`) and re-emits the
  frozen `EXPECTED_JET_DIGEST` `cfb36557...` **byte-for-byte**.
- **Two algorithmically distinct builds agree.**  Build: `X_k = m_0^(-1)m_k` plus
  trace formulas.  Audit: Gaussian elimination over `Q[t]/t^5` with unit pivoting,
  own chart-entry builder, own log-jet recursion, own interpolation.  Identical
  coefficient dictionaries.
- **Independent raw-minor anchor.**  Full interpolation of the true `36x36`
  Macaulay minor at `58M-36+1` depths, exact division by
  `q200^6 c300 curvature`, Newton on the top coefficients: `deg R_M = 46M-26` and
  `L_1..L_4` agree exactly at `M=6..12`.
- **Disjoint grids agree.**  Fit on `M=4..12, U=2..18, V` even `2..18` and on a
  disjoint `M=21..29, U=21..37, V=20..36`: identical polynomials, ruling out an
  interpolation artifact.  Plus `135` fully out-of-box points, `0/540` mismatches.
- **Out-of-box evaluation.**  `P_4(M,2^M,3^M) = c_4 D^8 L_4` exactly at
  `M = 6,7,8,9,11,13,19,23,25,31,40,62` -- astronomically outside the fit box --
  and at free points including `M=3`, negative `U`, negative `V`.
- **Checked here** (section 5 script): row count, degrees, support, integrality,
  digest, the `12` corner absences, and section 4 below.

## 4. It confirms THM-3011 in full

The top column gives

    [U^16] P_4 / c_4 = -(23/10)M^5 + 3M^4 - (47/24)M^3 + (707/3840)M - 15797937/128,

**exactly** THM-3011 section 2's `A_4`.  THM-3011 section 2a could predict only
the top three coefficients from the frozen table's laws; `707/3840` and
`-15797937/128` were reachable only by that section's large-width fit, and the
symbolic route now reproduces both.  So THM-3011's reconstruction is confirmed by
a genuinely independent computation, and with it

    p_4(N_M)/M^5 -> 46/5 - 90211/19440 = 88637/19440 = 4.55951646...,
    w = m_4/m_1^4 -> 337994862976/119272468005 = 2.83380455... = O(1),

i.e. the repaired third-edge invoice of THM-3000 `(19)` is cleared by an
unbounded asymptotic margin: `w/d->0>-923/60000`.

## 5. Corrections carried forward

- **Status.**  The build agent labelled this `PROVED-SYMBOLIC`; that is wrong and
  its own residual says so.  The degree bounds `M<=8, U<=16, V<=8` were **measured**
  by Newton divided differences with oversampling, not proved.  The correct label
  is the one in the frontmatter.
- **`m_0` is NOT the identity.**  THM-2997 section 4 previously conflated its
  locally inverse-normalized slice with the raw selected chart; that wording is
  now repaired.  For the raw chart, `det(m_0)` vanishes on
  `{U=1} u {V=1} u {V=2U-1} u {D=0}` and carries a high power of `D`
  (exactly `D^24` at `(M,U,V)=(7,3,2)`).  No symbolic inverse is needed because
  `U,V` stay numeric per evaluation.
- **"Absences are always corners" is FALSE at `j=3`.**  `P_3` omits `(5,0,3)`, and
  `(0,3)` is not a corner.  The statement happens to hold at `j=4`.  Record the
  enumerated lists, not the generalisation.
- **The `M=12` wall trap is real** and was handled: encoded wall degree `281` vs
  atlas `282` (the `(12,5)` Smith sporadic).

## 6. Scope

- The wall for `M>=27` is still the **encoded continuation**, not proved
  (THM-2969 certifies only `M=6..26`), so `p_4(N)>0` for `M=6..62` is
  unconditional only through `M=26` -- the same conditionality THM-2997 carries.
- No symbolic proof of the degree ansatz; no Lean; no proof that `p_4(N_M)>0` for
  `M>62`, though `P_4` now makes that a finite-analysis question rather than a
  computation.
- No structural simplification of the `717` coefficients was attempted.  `det m_0`
  carrying `D^24` suggests a `D`-adic or `(V-2U+1)`-adic factorisation is worth a
  look.
- Nothing here proves no-return, ULC, or GMC(2).
