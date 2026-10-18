# Audit F2: THM-4690(5) and the results note, sections 2 and 3

Independent referee, 2026-10-18, session `mac-mini-2026-10-18-sixcolour`. The subagent could not write this file
itself, so the session saved it verbatim from the referee's final message.

**Scope.**
- THM-4690 statement 5, its proof and its remarks, including "Why the denominator is 7".
- `05-knowledge/results/first_arrival_one_tick_coefficient_copies_gibonacci_20261018.md`, sections 2 (six-colour) and
  3 (Fibonacci from (5,1) and (7,1), Solfeggio merge), with rows 8-12 of its table.
- Statements 1-4 and note section 1 belong to audit F1 and are not covered here.

**Versions.**
- The audit started at commit `044e5bdf47`.
- The session kept committing in the same worktree, so the text was re-checked at `f6e0e18450`. Line numbers below
  refer to that commit.
- Changes since the start that touch this scope:
  - the K-direction connector probe (`kdirs.py`, `copy_connectors_kdirs*`), audited below;
  - the new title of THM-4690;
  - audit F1's rewrite of the "denominator 7" passage, which now reads correctly.

**Method.**
- Nothing from the session was run. Its scripts import `pysat`, which the referee kept outside its own toolset.
- All checks use code written for this audit, in this folder:
  - Python with python-flint for one linear solve and one determinant;
  - its own arithmetic in `F` and its own reduction to `F_361`;
  - its own degree-1 places;
  - a separate PARI/GP cross-check with exact `nfeltval` valuations.
- OEIS entries were read as HTML pages. Nothing was downloaded or installed.

## Verdicts

| Item | Verdict |
|---|---|
| A. THM-4690(5): statement, proof, remarks | **CORRECT WITH FIXES.** The 19-adic mathematics reproduces exactly. The "sheets" claim needs scoping, and the remarks draw a wrong inference from the census. |
| B. Note section 2 | **CORRECT WITH FIXES.** Numbers right. One premise is overstated, one heading overclaims, and the census is wrongly said to be explained by 2.2. |
| B. Note section 3 | **CORRECT WITH FIXES.** Every number is right. Receipts are missing for the squares and the prime terms, and one coincidence needs a label. |
| (1) Maximality at 19 | **CORRECT.** The proof should state the one-line reason. |
| (2) Valuations at the 24 primes | **CORRECT** |
| (3) Sheet-move counts | **CORRECT** |
| (4) Logic of 5 and its sheets extension | **CORRECT WITH FIXES:** scoping of the "recipe sum" and wording. |
| (5) Two-step criterion and the screens | **CORRECT WITH FIXES:** survivor count, and the interpretation of the census. |
| (6) Fibonacci facts | **CORRECT WITH FIXES:** receipts and a COINCIDENCE label. |
| (7) Solfeggio merge | **CORRECT** in the current text, after audit F1's fix; minor typing notes. |
| (8) Typing and wording | See the fixes. No PROVED claim is wrong as mathematics. |

## Evidence

### (1) Maximality at 19 (`f2_disc.py`)

**The discriminant.**
- The trace form of the power basis `ζ^a j^b w^c`, with traces taken from multiplication matrices and no tensor
  assumption, has determinant `disc(O) = 2^112 · 3^24 · 5^36 · 11^24 · 37^16`.
- This is exactly `disc(Φ15)^6 · disc(x^2+11)^24 · disc(g)^16 = (3^4 5^6)^6 (-44)^24 592^16`, with
  `disc(g) = 592 = 2^4·37`.
- 19 does not divide it. So 19 does not divide the index, `O` is 19-maximal, and 19 is unramified.

**Why `F` has degree 48.**
- `Q(√-11)` is not in `Q(ζ15)`, whose quadratic subfields are `Q(√-3)`, `Q(√5)` and `Q(√-15)`.
- `Q(w)` is a non-Galois cubic (592 is not a square), so it meets the abelian field `K` only in `Q`.

**`F` is CM.** `g` has three real roots (`disc > 0`), and `F = K·Q(w)`. So `conj` is complex conjugation in every
embedding.

**The equivalence.** "No 19 in the coordinate denominators ⟺ integral at every `𝔭 | 19`" holds in both directions,
since `O ⊗ Z_(19) = O_F ⊗ Z_(19)`.

**Roots of unity.** `μ(F) = μ_30`. Any quadratic subfield of `F` already lies in `K`, which rules out `i`. The other
candidates fail by degree or Galois type.

### (2) The 24 primes and their valuations

Computed by `f2_valuations.py`, with `F_361 = F_19[s]/(s^2 - 2)`, and by `f2_pari_check.gp`.

**The primes.**
- There are 48 ring homomorphisms `Z[ζ, j, w] -> F_361`, forming 24 Frobenius orbits, all of size 2.
- The types `h ≡ +j`, `-j`, `7` occur 8 times each.
- PARI: `K` has 8 primes above 19 with `(e, f) = (1, 2)`, and `Q(w, j)` has 3.

**Valuations.**

| Element | Valuation | Evidence |
|---|---|---|
| `b`, `c` | -1 at all 24 | coordinate denominator exactly 19; `19b`, `19c` nonzero at every prime; PARI agrees |
| `b/c`, `c/b` | 0 (19-units) | both 19-integral |
| `u` | -1 at `h ≡ -j` | `19u` nonzero there |
| `u` | **+1 at `h ≡ +j`** | `19/u` is 19-integral and nonzero there; the single-reduction rule only gives `>= 0` |
| `u` | 0 at `h ≡ 7` | `h - j` and `h + j` both nonzero there |
| `u_2`, `a = conj(d)`, the 30 roots of unity | 0 at all 24 | |
| `A` (54 points), `C` (10), `M` (7) | 19-integral | |
| the 34 nonzero spindle differences | 0 at all 24 | |

PARI gives the `u` row exactly: `N_{L/Q}(h ± j) = 5776 = 76^2`; `v(u) = 0, +1, -1` at the types `7`, `+j`, `-j`; and
`v(u_2) = 0` at all 6 primes of `Q(w, j, √5)` above 19.

**The construction is reproduced independently.** `f2_field.py` checks `u^3 = m_2 u^2 + m_1 u + 1`,
`u^2 = a + b u + c Y`, `|d|^2 - 1 = φ`, and `b, c ∈ K` (in fact in `Q(ζ5, √-11)`).

### (3) Sheet moves

- `|M - M| = 35`, counting 0. The 42 ordered nonzero differences collapse to 34, because the spindle has parallel edges.
  So there are 1224 nonzero pairs.
- The number of `h ≡ 7` primes at which `b Δm_2 + c Δm_3` is integral is distributed `{0: 1166, 1: 16, 2: 42}`.
- No pair is 19-integral. Every nonzero `b Δm` and `c Δm` has a pole at all 24 primes.
- The same distribution holds at the `+j` and the `-j` primes, and over the 8 primes of `K` in PARI. All these elements
  lie in `K`, and each prime of `K` above 19 has exactly one prime of each type above it.

### (4) Logic

**Correct:**
- The ultrametric argument: steps integral at `𝔭` give a displacement integral at `𝔭`.
- `A`, `C` and `a` are 19-integral.
- For the sheets, integrality at **all eight** `h ≡ 7` primes is the right hypothesis: 42 moves are integral at two of
  them. For one factor, integrality at a single prime above 19 already suffices.
- `ζ^l u^k u_2^m q`, with `q` 19-integral, has `v >= 0` at the `h ≡ 7` primes for all `k, m ∈ Z`.
- The `conj(γ)/γ` remark. PARI finds that none of the 8 primes of `K` above 19 is fixed by `conj`, and builds an
  explicit `δ` with `δ·conj(δ) = 1` and `v = -1` at one prime.

**The problem is scoping (fix 1).** The "recipe sum" is the plain Minkowski sum in `F`. In clock_duality's recipe
graph, the vertices are `x_0 + u x_1 + v x_2 + uv x_3 + Y x_4`.
- There, the `u`-square contact family `q(a, b, 0, 0, c)` joins `(m_1, m_2, m_3)` to `(m_1 + q, m_2 + q, m_3 + q)`
  whenever `q` is a spindle unit difference. The planar step is `q(a + bu + cY) = q u^2`, a unit.
- `golden_contact_circuit_20261017i.py` enumerates exactly these `q` (unit differences of `(C + bM)/b`; family count
  194600).
- So the sheets of the recipe graph are joined, by steps integral at the `h ≡ 7` primes. The claim is true only for the
  auxiliary sum.

### (5) Two-step criterion and the screens

**The derivation is correct.**
- Write `z = (t/2)(1 + s)`. Then `|z| = |t - z| = 1` forces `conj(s) = -s` and `s^2 = 1 - 4/N(t)`.
- Conversely, suppose `D = 1 - 4/N(t)` is a square `s^2` in `F` and `|t| < 2` in one embedding. Then `conj(s) = ±s`.
  The sign `+` would put `s` in the totally real subfield, so `D >= 0` in every embedding, contradicting `D < 0` where
  `|t| < 2`. Hence `conj(s) = -s`.
- At `|t| = 2`, `D = 0` and the middle point is the midpoint.
- A pre-filter in a single embedding is a valid necessary condition.

**The certificate is applied correctly.**
- The screen certifies only when `red(X) ≠ 0` (the session's code skips residue 0) and `red(X)` is a non-residue.
- The places are `p ≡ 1 (mod 15)` with `p > 1000`, which do not divide `disc(O)`. So the coordinate reduction is the
  residue map of a degree-1 prime, and `X` is a unit there.
- If `red(X) = 0` at every place, nothing is certified and the test is reported as a survivor. That never happens.
- The newer rule, skipping places where a coordinate denominator is divisible by `p`, is conservative and valid.

**Full independent recount** (`f2_connectors.py`, `f2_kdirs.py`). The places are all embeddings at 12 primes from
1021 to 3391, with a fallback to 40 primes.

| Family | Two-step tests | Three-step tests | Uncertified |
|---|---|---|---|
| (none) | 172 + 118 | – | 0 |
| 22 unit differences of `C` | – | 4377 + 3017 | 0 |
| 90 integral directions | – | 17854 + 12340 | 0 |
| 120 directions `ζ^l u^k` | – | 23729 + 16508 | 0 |
| 240 K-directions | – | 47515 + 33031 | 0 |

**Every recorded survivor is reproduced and certified.**
- The 15 survivors of the 22-, 90- and 120-direction runs are certified at the session's primes and again at
  `p = 10321`. The 37 K-direction survivors are certified at `p = 1021`.
- `|t|^2` agrees with the recorded value to 6 decimals in every case.

**Controls** (`f2_controls.py`). The reduction is a ring homomorphism, and conjugate places behave as expected.
Negative control: actual squares were certified as non-squares 0 times in 1440 tests.

**K-directions.**
- Each of the four `γ = j - c_0 - c_1 ζ` vanishes at exactly one prime of `K`, not at its conjugate. Together they hit
  all 4 conjugate pairs.
- Each `δ` has `v = -1` at the 3 primes of `F` above `𝔮` (one of each type), `+1` above `conj(𝔮)`, and 0 elsewhere.

**Duplicates.** The 2 survivors of the 22-direction run are the same tests as 2 of the 3 survivors of the
90-direction run: `U0[8] = D90[39]` and `U0[20] = D90[88]`. So "5 survivors" is wrong; the correct count is 3 (fix 6).

### (6) Fibonacci (`f2_gibonacci.py`, `f2_windows.py`)

**All confirmed:**
- the 24 listed terms, and "5 / 7 followed by A022096 / A022098" (30 OEIS terms each);
- `G_n = aF_(n-1) + bF_n = (aL_n + (2b-a)F_n)/2`;
- `(7,1) = 4L_n - L_(n+1) = L_(n+2) - L_(n-3)` and `(5,1) = F_(n+3) - F_(n-4)`, with the 7- and 5-window sums;
- `θ(5,1) = 5φ - 4 = φ^-4(φ^7 - 1)` and `θ(7,1) = 7φ - 6 = φ^-3 √5 (φ^5 - 1)`;
- `θ φ^n = G_n φ + G_(n-1)`, and `N(φ^r - 1) = -L_r` for odd `r`;
- Cassini `G_n^2 - G_(n-1)G_(n+1) = (-1)^n · 29`, resp. `(-1)^n · 55`, for `-30 <= n <= 2000`;
- the growth constants 4.0902 and 5.3262;
- `α = (25 - 3√5)/10` and `(7 - √5)/2`, with `α·conj(α) = 29/5` and `11`;
- `round(αφ^n)` is exact from `n = 4`, resp. `5`, and fails at `n = 3`, resp. `4`;
- the eigen-sequences: `λ = 6` of order 14 mod 29; `λ = 3` of order 4 mod 5; `λ = 8` of order 10 mod 11. No term is
  divisible by these primes.

**Squares** (`n < 3000`): (5,1) has `G_1 = 1` and `G_10 = 225`; (7,1) has `G_1 = 1`, `G_3 = 9` and `G_15 = 3249`;
Fibonacci has 0, 1, 1, 144.

**The Siegel identity** `(2G_(n+1) - G_n)^2 = 5G_n^2 + 4(-1)^n N` holds.

**Prime terms.**
- (5,1): 5, 7, 13, 53, 139, 953, 44771, 189653, ... These match OEIS A161611, "Primes that are a sum of 7 consecutive
  Fibonacci numbers", from 53 on.
- (7,1): 7, 17, 43, 181, 293, 13763, 399601, ...

**Periods** (`m < 3000`): (5,1) has exactly the Pisano periods at every modulus, multiples of 29 included. (7,1) has
exactly the Lucas periods, as Purath's comment states; these differ from Pisano's only at multiples of 5.

**OEIS.**
- Titles confirmed; the McRae (A022096, 2014), del Arte (A022098, 2013) and Purath (A022098, 2019) comments are present.
- Aside: A022096's formula "a(n) = Fibonacci(n+3) - Fibonacci(n-4)" (Dresden and Neale) is off by one at OEIS offset 0.
  It returns 5, 1, 6, 7, ..., our `G_n`. McRae's comment gives `F(n+4) - F(n-3)` in OEIS indexing. The note's G-indexed
  formula is right.

### (7) Solfeggio (`f2_solfeggio.py`)

**Mod-9 facts, all confirmed:**
- `(7,1) ≡ (F_9, F_10)`, so `G_n ≡ F_(n+9)`.
- `G_n(5,1) ≡ 4F_(n+14) ≡ 5F_(n+2)`. These two are the only unit-multiple representations.
- The digits 3, 6, 9 occur at `n ≡ 3`, resp. `2 (mod 4)`, and the 24-cycle is as printed.
- The golden classes are `±1` and `±2`. Multiplication by 4 fixes each fibre and cycles `±1 -> ±2 -> ±4 -> ±1`.

**The denominator 7.**
- Recomputing the THM-4690(4) tabulation gives `34/7` and `5/7` exactly. The 63 is `Σ_(3|L) 4^-L = 1/63`.
- The 9 in `63 = 9·7` cancels against the weight factors: `(3/8)·24 = 9` for `E_0` and `(1/8)·18 = 9/4` for `SIB`.
- Only `B(x) ≡ x + 1 (mod 3)` is used. The current text, after audit F1, says exactly this and types the Solfeggio link
  STRUCTURAL.

**Numerology.** The primes are `{5, 11, 13, 19, 29, 71, 107, 139}`: 5 ramified; 11, 19, 29, 71, 139 split; 13, 107
inert. The hypergeometric tail is `P(>= 3 | 35, 8, 5) = 0.0665`, so "p ≈ 0.07" is right.

## Required fixes

THM = `01-canon/theorems/THM-4690-first-arrival-ladder-port-and-integer-clock-decompositions-and-the-nineteen-adic-copy-obstruction.md`;
NOTE = `05-knowledge/results/first_arrival_one_tick_coefficient_copies_gibonacci_20261018.md`.
Line numbers are at `f6e0e18450`.

### Substantive

These change what a reader concludes.

**Fix 1. Scope the sheets claim to the auxiliary sum.**
- Where: the THM title (line 3), "(or two sheets) of clock_duality's cubic golden recipe"; THM statement 5, "Sheets of
  the whole recipe sum" (lines 161-166); NOTE 2.2, "For the whole recipe sum ..." and "So the 49 sheets ..." (lines
  291-294); NOTE results-table row 8 (line 36).
- Change: in the title, write "(or two sheets of the Minkowski sum of its three translated factors)". In the THM and the
  NOTE, add: "Here the recipe sum is the auxiliary Minkowski sum `(A + aM) + (C + bM) + (C + cM)` in `F`. A sheet is its
  subset with fixed `(m_2, m_3)`; sheets are disjoint because no nonzero move is 19-integral. This is not the vertex set
  `x_0 + u x_1 + v x_2 + uv x_3 + Y x_4` of the recipe graph. There the `u`-square family `q(a, b, 0, 0, c)`, with `q`
  a spindle unit difference, joins `(m_1, m_2, m_3)` to `(m_1+q, m_2+q, m_3+q)` by the planar unit step `q u^2`; for
  example `q = 1` joins sheet `(0,0)` to `(1,1)`. So the recipe graph's sheets are not separated."
- Why: as written, the claim reads as a property of clock_duality's recipe, and there it is false.

**Fix 2. Retract "2.2 explains the census".**
- Where: NOTE lines 313-314, "2.2 explains the second family completely: `u` cannot touch the poles at the `h ≡ 7`
  primes. The third family has the right kind of pole, but one `δ` covers only 3 of the 24 primes."
- Replace with: "2.2 does not explain these negative results. Only one step of each three-step test is restricted. The
  other two are arbitrary unit vectors of `F`: in a two-step remainder `e_1 + e_2 = t'` the steps are `(t'/2)(1 ± s)`,
  and they carry every pole of `t'`. Valuations allow such connectors in all three families. What fails is the
  squareness of `N(t')(N(t') - 4)`."
- In the THM remark title (line 297), "Connector census, inside `F`, agrees with 5." becomes "Connector census inside
  `F` (consistent with 5, not implied by it)."
- Why: 5 only constrains paths all of whose steps are integral at the `h ≡ 7` primes. The integral family, where
  valuations give no obstruction at all, also found nothing, for the same squareness reason.

**Fix 3. Valuations do not force the length.**
- Where: NOTE 2.4 (lines 321-323), "The next search should combine several `δ`'s so that the step poles can cover all
  24 poles of `b(m' - m)`. That means four or more steps, or `δ`'s attached to products of primes."; and the THM remark
  (lines 316-318), "A connector must clear the poles at all 24 primes, so `δ`-steps covering several pairs ... are the
  next family to try."
- Change: keep the pole-covering requirement, but say that free steps can supply the poles (the two-step remainder
  already does). So this does not force four or more steps, and covering by named `δ`'s is a search heuristic, not a
  consequence of 5.
- Why: the same faulty inference as fix 2, and it steers the next search.

**Fix 4. Correct the independence premise.**
- Where: NOTE 2.1 (line 271), "Each copy is coloured independently, so forcing cannot propagate between copies."
- Replace with: "No axis unit edge joins two copies, so each copy is coloured independently within the factor's own unit
  graph. The contact families do couple choices across copies (golden_contact_circuit section 4; the `u`-square family
  shifts all three spindle indices together). The saved witnesses show that this coupling blocks none of the 570
  patterns."
- Why: the sentence contradicts the cited source.

**Fix 5. Retitle NOTE 2.2.**
- Where: NOTE line 274, "### 2.2 Why short connectors cannot exist: 19-adic separation (THM-4690(5))".
- Replace with: "### 2.2 What 19-adic valuations exclude (THM-4690(5))".
- Why: 2.2 excludes paths, of any length, whose steps are all integral at the `h ≡ 7` primes. It excludes no "short
  connector" in general; only the finite census of 2.3 does that, for named families.

### Minor

These are wrong numbers, proof wording and receipts.

**Fix 6. Survivor count.**
- Where: THM line 299, "(5 survivors of the 12-place screen, ..."
- Replace with: "(3 distinct survivors of the 12-place screen; the 2 survivors of the 22-direction run are 2 of these 3,
  ..."
- Why: double count (`U0[8] = D90[39]`, `U0[20] = D90[88]`). NOTE 2.3 already says 3.

**Fix 7. Garbled sentence in the proof.**
- Where: THM line 274, "A path with steps in a set that is integral at every prime above `𝔭` has an integral
  displacement."
- Replace with: "A path whose steps are all integral at `𝔭` has a displacement integral at `𝔭` (ultrametric inequality)."
- Why: `𝔭` is already a prime of `F`, so "every prime above `𝔭`" means nothing.

**Fix 8. Notation clash.**
- Where: THM line 272, "`v_𝔭((c' + b m') - (c + b m)) < 0`".
- Replace with "`v_𝔭((γ' + b m') - (γ + b m)) < 0` for `γ, γ' ∈ C`".
- Why: `c` is the circuit coefficient.

**Fix 9. Give the reason for maximality.**
- Where: THM line 265, "Maximality at 19".
- Add: "`F = Q(ζ15) ⊗ Q(j) ⊗ Q(w)` (degree 48: `√-11 ∉ Q(ζ15)`, and the non-Galois cubic meets the abelian `K` only in
  `Q`), and the order is the tensor product of the three monogenic orders. So
  `disc = disc(Φ15)^6 disc(x^2+11)^24 disc(g)^16 = 2^112 3^24 5^36 11^24 37^16` (audit F2, `f2_disc.out`). Conversely,
  a 19-free denominator means integral at every `𝔭 | 19`."
- Why: a PROVED theorem asserts the support without a reason.

**Fix 10. What the reduction rule decides.**
- Where: THM line 290, "With this, the listed valuations and the counts 1166 / 16 / 42 are exact finite computations."
- Add: "A vanishing reduction only gives `v >= 1 - k`. `v(u) = +1` at `h ≡ j` and `v(u) = 0` at `h ≡ 7` come from the
  relative norms, and `nineteen_places.out` prints `>=0` there." In the contact-unit proof, add the symmetric line:
  "Symmetrically, `v(h - j) = 1` exactly at the prime with `h ≡ j`."
- Why: the receipts do not show these two values; the norm argument does.

**Fix 11. Make the "Consequence" consistent.**
- Where: THM lines 168-169, "To change sheets it needs a pole at an `h ≡ 7` prime, which `u` does not have. The contact
  unit `u = ...` is one."
- Replace with: "To change sheets it needs a pole at an `h ≡ 7` prime. The contact unit `u = (h - j)/(h + j)` has poles
  above 19, but only at the eight `h ≡ -j` primes; it is a unit at the `h ≡ 7` primes."
- Why: as written it reads as self-contradictory.

**Fix 12. Reproduction details.**
- Where: NOTE line 434, "`ROOT_OF_UNITY` / `U_DIRECTIONS`".
- Change to "`ROOTS_OF_UNITY` / `U_DIRECTIONS` / `K_DIRECTIONS`". Add `kdirs.py` and `copy_connectors_kdirs*` to NOTE
  section 4 and to the THM `verification:` list (line 34).
- Why: the scripts read `ROOTS_OF_UNITY`; the documented name silently runs the 22-direction default. The K-probe is
  cited but not listed.

**Fix 13. Receipts for squares and prime terms.**
- Where: NOTE lines 379-386, "Prime terms ..." and "These lists are FINITE-EXACT for `n < 1500`".
- Change: add the computation to `gibonacci_51_71.py` and its `.out`, or cite `audit_F2/f2_gibonacci.out` (squares for
  `n < 3000`, prime terms for `n < 120`). Cite OEIS A161611 for the (5,1) primes.
- Why: neither list is computed in the session's script or receipt.

**Fix 14. Label the window coincidence.**
- Where: NOTE line 358, "So the start 5 pairs with a window of 7, and the start 7 pairs with a window of 5."
- Add: "(COINCIDENCE of small numbers: `5^2 + 5 - 1 = 29 = L_7` and `7^2 + 7 - 1 = 55 = 5 L_5`. Among `(a, 1)` with
  `a <= 40`, only `a = 1, 2, 5, 7` are window sums of Fibonacci or Lucas numbers. `(3, 1)`, of norm `-11 = -L_5`, is
  not; see `f2_windows.out`.)"
- Why: the wording suggests a duality.

## Optional improvements

**Stronger or clearer statements.**
- `b, c ∈ Q(ζ5, √-11)`, and `A, C, M ⊂ K`. Their valuations are therefore equal at the three primes of `F` above each
  prime of `K`; the `h`-type matters only for `u`-type directions. The sheet statement could be phrased with the 8
  primes of `K`.
- The single-factor statement already holds when the steps are integral at **one** prime above 19.
- The proof's "spindle" bullet uses that `m - m'` is 19-integral. Without that, "norm prime to 19" does not give a
  `𝔭`-unit. Say so.
- Write the two-step derivation into the NOTE, including the CM sign argument.
- "Example": the 30 roots of unity are `±ζ15^l`, since `ζ = ζ15` has only 15 powers. `μ(F) = μ_30`.

**Periods** (NOTE lines 376-377). State them as:
- (5,1) has the Pisano periods at every `m < 3000`, multiples of 29 included, because `ord_29(6) = 14 = π(29)`.
- (7,1) has the Lucas periods; this is Purath's OEIS comment, checked for `m < 3000`.
- PROVED for `gcd(m, N) = 1`: `G = (aI + Q)·F` with `det(aI + Q) = a^2 + a - 1 = -N`. Likewise `(7,1) = (4I - Q)L`,
  with determinant 11.

**Citations.**
- Fibonacci squares: list 0, and cite J. H. E. Cohn, J. London Math. Soc. 39 (1964) 537-540.
- Credit the Dresden–Neale and Dresden–Nyc formulas on OEIS (2022), noting the A022096 index shift above.
- Context: Guyer and Mbirika, "GCD of sums of k consecutive Fibonacci, Lucas, and generalized Fibonacci numbers",
  J. Integer Seq. 24 (2021).

**Numerology** (NOTE lines 415-418).
- The overlap of 3 counts the seed `G_0 = 5`. Without it the overlap is {13, 139} out of 4 prime terms, with
  `p ≈ 0.22`.
- With the population restricted to primes 5..150 (2 and 3 are excluded by construction), `p = 0.078`.
- "Six of eight split or ramify": 5 split out of 8 has binomial tail 0.36.
- The NUMEROLOGY label stays right.

**Smaller points.**
- "Pushed through the linear part of `B`" (NOTE line 399): the multiplier `±4` is forced by the norm class `±2`
  (`G ≡ 4F_(n+14) ≡ 5F_(n+2)`). Type the `B`-reading as ANALOGY.
- "`diagonal_19adic_check.py` confirms independently" (THM line 293): change to "by a second criterion (same field
  code)". Audit F2 supplies independent code.
- Script docstrings: `copy_connectors_screen.py` promises a gp `nffactor` certification it never performs;
  `copy_connectors.py` has a garbled sentence about survivors; the session's scripts need `pysat` at import (through
  `golden_forcing_search_20261017.py`), which NOTE section 4 does not mention.
- "Why the denominator is 7": add that the 9 in 63 cancels against the weights.

## Reproduce

Run from `04-computation/experiments/first_arrival_20261018/audit_F2`. Each job takes under 3 s and under 50 MB. Only
`python3` (numpy, sympy, python-flint) and `gp` are used.

```bash
python3 f2_field.py        > f2_field.out        # construction: u^3, u^2 = a + b u + c Y, denominators 19
python3 f2_disc.py         > f2_disc.out         # disc(O) = 2^112 3^24 5^36 11^24 37^16
python3 f2_valuations.py   > f2_valuations.out   # 24 primes, valuations, 1166/16/42, no 19-integral move
gp -q f2_pari_check.gp     > f2_pari_check.out   # exact nfeltval cross-check; conj-splitting; explicit delta
python3 f2_connectors.py   > f2_connectors.out   # full recount + survivor certificates (writes f2_connectors.json)
python3 f2_controls.py     > f2_controls.out     # homomorphism / negative controls; extra prime 10321
python3 f2_kdirs.py        > f2_kdirs.out        # K-direction probe: gammas, delta valuations, 47515 + 33031 tests
python3 f2_gibonacci.py    > f2_gibonacci.out    # section 3.1-3.3
python3 f2_solfeggio.py    > f2_solfeggio.out    # section 3.4, 34/7 and 5/7, hypergeometric p
python3 f2_windows.py      > f2_windows.out      # which (a,1) starts are window sums
```

Sources read: oeis.org/A022096, oeis.org/A022098, oeis.org/A161611, and arxiv.org/abs/2104.12262 (abstract page only).
