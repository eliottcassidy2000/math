# Linear energy toward positive density: micro/macro balance of one-apex collisions, Poisson-pair second moments, and face-local relaxations of the permutohedron

opus-2026-11-06-S84. Status: PROVED (THM-4890, including the audit's sharp balance lemma U1; THM-4891), PRE-REGISTERED
(HYP-9482: P1-P5, P6(a), P8 and the addendum P9-P11 confirmed, P6(b) and P7 refuted, MISTAKE-728), FINITE-EXACT,
HEURISTIC where marked. **Linear energy E_n = O(2^n), positive density, density one and Chan-Pak remain OPEN.**
Independent audit: PASS WITH CORRECTIONS, with 2 majors on reading, 13 minors and 4 upgrades, all applied
(`audit/AUDIT.md`, MISTAKE-730).

The owner's directive was:
- prove linear energy for positive density;
- think micro/macrocosm;
- merge the {3,10} / {15,25} / {7,21} threads and the permutohedron sessions;
- use the x+2 / x*2 duality and the recognizable gate;
- treat the parallel-coherence restriction as the place where 7 and 21 are decided.

## 0. Outcome

- **(1) Linear energy is not proved.**
  - **What is proved instead.**
    - **Balance (THM-4890 (C)).** Exact collisions of one-apex counts are micro/macro balanced. Two reduced words with
      the same count differ on a window. One window is 1...1 and the other 0...0. The carry-corrected common-prefix and
      common-suffix lengths differ by at most 1, |A - S| <= 1. The sharp form and its proof are the independent audit's
      upgrade U1, which closed my obligation 2.
    - **Symmetry (C').** The converse involution exchanges the two lengths.
    - **Reduction (E).** Linear energy reduces to reduced words.
    - **Criterion (D).** Linear energy is equivalent, up to constants, to criterion (DC): collisions with different
      macro blocks and different micro blocks (bottom n/2 - 2 bits) number O(|X_n|). With independent uniform micro
      self-energies the expected energy is at most 109 |X_n|. The exact value is 1.64-1.99 |X_n|.
  - **What this means.** In the random model the macro lattice creates no excess collisions. The true micro values,
    however, give about 3.7 times the random-model count of (DC) collisions, about |X_n| of them forced by the converse
    involution. Linear energy needs this ratio to stay bounded (criterion (DC), section 4). That is a bounded-
    correlation statement, not a decorrelation statement (audit MAJOR-1).
- **(2) Second moments match the Poisson-pair model (heuristic) to about 2 % (HYP-9482).**
  - The limit law XY of g(v)/2^n is PROVED (Session 16, `BINARY_PRODUCT_DENSITY.md` section 3).
  - The model adds the heuristic assumption that converse classes cluster at the unit scale as a Poisson process with
    that intensity.
  - At first order this is kind-pasteur's HYP-9480 (A) (kappa = 1/2, density 1/2). At second order the model gives
    closed forms:

  | quantity | model | n = 26 (pre-registered) | n = 28 (pre-registered addendum) |
  |---|---|---|---|
  | E^Y/\|Y\| (reduced words) | 14 - 16 ln 2 = 2.90965 | 2.92960 | 2.92480 |
  | E/\|X\| (all words) | 6 - 4 ln 2 = 3.22741 | 3.25257 | 3.24734 |
  | smoothed slope per unit K | 12 - 16 ln 2 = 0.90965 | 0.90872 (K = 8) | 0.90936 (K = 8) |
  | attained values per \|Y\| | 0.40343 | 0.40194 | 0.40227 |
  | multiplicity histogram | Poisson-pair mixture | within the bands for q = 1..4 | within the bands for q = 1..4 |

  - The excess over the model is shrinking: 0.032, 0.025, 0.020, 0.015 at n = 22, 24, 26, 28, a factor of about 0.78
    per two levels.
  - The upper tail is 15-35 % heavier than the model at multiplicities 12-16. It refuted my eyeballed
    maximum-multiplicity bound (P7). The model-derived bound at n = 28 held.
  - The E/|X| row reproduces HYP-9479 (S83), which was published before the registration. It is not out of sample, and
    P1 and P8 are not independent of it.
  - This supports, but does not prove, linear energy with constant 3.2274 (conjecturally exact).
- **(3) Neither half of parallel coherence alone keeps 7 and 21 out (THM-4891).**
  - Parallel coherence on the permutohedron is the conjunction of square coherence and hexagon coherence.
  - **Square-coherent orientations** are set-context tournaments: the arc a-b may depend on the set of earlier
    elements.
    - They have at most n!/2^floor(n/2) sinks (Session 188's bound).
    - The maxima are 6, 25 and 90 at n = 4, 5, 6 (the bound is 30 at n = 5).
    - The spectra are the full intervals 0..max for n = 3..6.
  - **Hexagon-coherent orientations** let the arc depend on the relative order of the other elements. They realise 7
    sinks at order 4, acyclically, with 48 orientations.
  - At order 5 both relaxations realise 7 and 21 acyclically.
  - So neither face type alone keeps 7 and 21 out. Both together are zone-wide coherence (tournaments). Whether an
    intermediate relaxation still excludes them is OPEN (audit MAJOR-2).

## 1. Prior and concurrent work (credited)

- **One-apex analytic line.**
  - THM-4887 (S83): smoothed energy linear at scales 2^(n/3) to 2^(n/2); E_n = O(n 2^(4n/3));
    A(X) >= X^(2/3)/(13 ln X + 102).
  - THM-4886 (S83): Type II bound.
  - THM-4885 (S82): 1 <= g(v) - v v*/2^n <= 1 + #10(v); readout masses.
  - Session 193's concatenation laws (1) and (2).
  - Session 16's W(X) = X ln 2, and its proof that g(v)/2^n -> XY with density 4 ln(4t), 4 ln(1/t) on reduced words
    (`feedback_word_hamilton_spectrum_20261020/BINARY_PRODUCT_DENSITY.md` section 3).
  - HYP-2199 (monad): the single-core gap set has density 1/2.
  - HYP-9480 (A) (kind-pasteur-2026-11-06-S13): the Poisson miss law exp(-kappa lambda), with kappa = 1/2 for the
    root, i.e. the one-apex image. Pre-registered and confirmed today.
- **Structural lanes on the same directive, today.**
  - **Session 210** (Codex, `session210/`):
    - the coherent-permutohedron bijection (Hamiltonian paths = sinks);
    - dropping coherence admits every count to n!/2;
    - the gate Gamma(T): a -> b -> T -> a, recognised by d+(a) = d-(b) = 1, with HC(Gamma(T)) = H(T);
    - the {3, 10} = (H - 1)/2 orbit carrier;
    - the cover-polynomial reconstruction.
  - **Session 211** (Codex, `session211/`):
    - square-sum ports, with p(Q23) = 3 and p(Q25) = 10 Hamiltonian chains;
    - the recursion bridge: padding and the gate (order n + 2), against Skew-Sylvester doubling (order 2n);
    - the score permutohedron.
  - **mac-mini-2026-11-06-microcosm:** THM-4900 claimed (marked resolvent under series, gate and doubling).
- **{15, 25}.** klein's analogy ledger and the monster-levels note record the following, all typed COINCIDENCE or
  NUMEROLOGY:
  - (15, 25) are the Busch strong floors f(6), f(7), and 21 lies between them;
  - the square-sum onsets are {15, 23, 25}.
- **Session 188** (`session188_structural_bridges_20261104/equal_nine_ambient/PROOF.md` (1)): the bound n!/2^floor(n/2)
  for tournaments and its orbit proof. The proof applies verbatim to set-context tournaments.
- **Permutohedron history.**
  - THM-207 (2026-03, opus S89b): the transitive-tournament flip graph is the permutohedron.
  - THM-1825 (2026-07, boxeph S190): the q-Redei census with the backward-arc statistic.
  - The LRC permutohedral-fan sessions S521-S538 (2026-06/07): HYP-2001, HYP-2003.
  - Sessions 210 and 211 (today).

## 2. The micro/macro identity and the hyperbolic picture

Notation:
- g(0) = 1, g(2v+1) = g(v), g(2v) = g(v) + v + (v mod 2);
- X_n = [2^(n-1), 2^n);
- reduced words Y_n = even words in X_n;
- C_w is the reversed complement of a w-bit word, and v* = C_n(v).

- **Symmetric two-block law** (Session 193 (1) with the converse). For c = C_m(l):

      g(h 2^m + l) = g(h) + g(c) + h c + [h odd][c odd] - 1.

  The macro block h and the backwards micro block c enter symmetrically. The converse v -> v* swaps them.
- **Hyperbolic flow.** The law holds at every split point m. Moving the split one bit up sends (H, Gamma) to
  (floor(H/2), 2 Gamma + 1 - b_m), and the count stays fixed.
  - g(v) - v v*/2^n lies in [1, 1 + #10(v)] (THM-4885).
  - So every word traces a staircase along the hyperbola H Gamma ~ g(v).
  - Two words collide when their staircases carry the same exact invariant.
- **Twisted grid.** For n = 2j, write v = aN + l and k = C_j(l), with N = 2^j. Then

      (v/N, v*/N) = (a + C_j(k)/N, k + C_j(a)/N).

  - So the point set {(v/2^n, v*/2^n)} is the reflected Hammersley (0, n, 2)-net.
  - Each grid cell (a, k) carries one point. Its offset is the van der Corput (bit-reversal) value of the other index.
  - A collision is a pair of net points on the same hyperbola, up to the exact local correction #01 - L.
  - At the unit scale this is a coincidence between the van der Corput sequence and fractional parts {T/k} along
    Beatty sequences.

## 3. THM-4890: collisions are balanced

For distinct reduced words with equal counts, define:
- a, the common-prefix length;
- s, the common-suffix length;
- the window values x, x', and the carry depths dX, dXi of x - x' and of C_w(x) - C_w(x').

With A = a - dXi and S = s - dX, the theorem gives:
- X Xi < 0, so one window is 1...1 and the other 0...0;
- |A - S| <= 1;
- if w = 1, ceil(h/2) = ceil(tau/2).

- **Proof** (the audit's U1).
  - The three-block law gives h Xi + tau X = g(x') - g(x) + O(1).
  - Same signs are impossible because g(x) - g(x') is forced to the wrong side.
  - With opposite signs the exact local formula g(z) = z C_w(z)/2^w + 1 + L(z), with -1 < L <= #10(z), gives
    g(v) - g(v') = T X - H|Xi| - delta. If A >= S + 2, this is negative, and the converse gives the symmetric case.
  - My first version proved only |A - S| <= 2 or max(A, S) <= 3.
- **Checks.**
  - 0 violations on all collision pairs for n = 10..20 (`check_balance.out`).
  - 0 violations on all 258.6 M collision pairs for n <= 28 (`audit/aud_pairs_*.out`). A - S takes only the values
    -1, 0 and 1.
  - So HYP-9482 P5 (|a - s| <= 1 + max(dX, dXi)) is a theorem.
- **Consequences** (data, `outofsample_21_26.out`).
  - Non-converse pairs at prefix level a number 0.347, 0.259, 0.153, 0.083, 0.043, 0.022, 0.0115, 0.0058, 0.0029,
    0.0014 per |Y| at n = 26 (a = 1..10): a halving law.
  - The bulk is global (a, s <= 3), so restricting attention to deep windows cannot prove linear energy.
- **The other parts of THM-4890.**
  - (E) reduces linear energy to reduced words: E_n <= 1 + 2 sum_j E^Y_j.
  - (D): with independent uniform micro self-energies the expected energy is at most 109 |X_n| (exactly 1.64-1.99
    |X_n|). The true statement is equivalent to criterion (DC), which asks for bounded correlation (section 4).

## 4. Why linear energy is still open: the precise obstruction

At the split P = 2^(n-m-1), C = 2^m with P/C in {8, 16}, the off-diagonal energy is, up to the 33 |X_n| from the
other pair types (THM-4890 (D)):

    #{(h, h', l, l'): h != h', l != l', g(l') - g(l) = h c - h' c' - (g(h') - g(h)) + parity}.

- **Where the solutions lie.**
  - The solutions lie in the thin strips |h c - h' c' - D| <= C. There, by THM-4887's strip lemma, there are at most
    2C^2/P + 3C/beta + 5 lattice points per pair (h, h'), and O(log) on average.
  - On a strip, c' ~ (h/h') c + const is a Beatty-type function of c.
  - The micro values satisfy g(l) ~ l C_m(l)/2^m. This is l times its own bit reversal.
- **The criterion.**
  - (DC) asks that the micro self-energies correlate with the macro strips only boundedly: the type-4 count is at most K
    times its random-model expectation, for an absolute K.
  - The observed ratio is 3.3-3.7 against uniform micro values and 1.6-1.8 against their empirical law. It is still
    slowly rising at n = 22, and about |X_n| of the count is forced by the converse involution.
  - So (DC) is a bounded-correlation statement about bit reversal and multiplication along Beatty sequences.
- **Cheaper routes that fail.** Each was examined this session.
  - Cauchy-Schwarz over h loses the macro phases. It gives E <= P |X|.
  - Pointwise Type II bounds (THM-4886) give |S_n|^2 << 2^(3n/2) for typical theta. That is the wrong order to
    integrate to 2^n.
  - A recursion for tilted cross-energies with a supremum over the shift closes in the main term, but the supremum
    discards the exactness. Additive losses then reach the critical depth where 2^t/c ~ 1, and the recursion only
    recovers the trivial count.
    - Data do support the uniform bound sup_D N_t(c, c', D)(c + c' + 1)/2^t <= 4.7 (t = 10, 12; exploratory,
      `tilted_energy_sup.py`, `tilted_energy_sup.out`).
  - Smoothing at scale K ~ n (via v v*) is no easier. The product v v*/2^n has the same self-similar structure, and its
    near-collision count is again the unit-scale problem. Data: linear in K, with no log factor.
- **A route that might work, not attempted.** Mauduit-Rivat-type digital exponential sums: carry-propagation
  truncation plus the Fourier transform of the bit-reversal phase, which factorises over digits. The bilinear coupling
  h c is the part that must be handled.

## 4b. Measuring criterion (DC) directly (EXPLORATORY: `dc_ratio.py`, `dc_ratio_12_18.out`, `dc_ratio_20_22.out`)

At the split m = floor(n/2) - 2, the energy splits exactly as E_n = |X_n| + T3 + T4:
- T3 counts pairs with the same micro block;
- T4 counts the (DC) quadruples.
The program checks this identity at every n.

R_n compares T4 with its expectation T4_model when the micro values g(l) are replaced by i.i.d. draws with their own
empirical law. R*_n does the same after removing the converse pairs, which the independent model cannot see. The audit
(`audit/aud_randmodel.out`) reproduces R_n digit for digit, and reports 3.3-3.7 against uniform micro values.

| n | T3/\|X\| | T4/\|X\| | converse in T4 | T4_model/\|X\| | R_n | R*_n |
|---|---|---|---|---|---|---|
| 12 | 0.482 | 1.447 | 0.438 | 0.883 | 1.639 | 1.144 |
| 14 | 0.315 | 1.742 | 0.469 | 1.019 | 1.710 | 1.250 |
| 16 | 0.193 | 1.995 | 0.484 | 1.105 | 1.805 | 1.367 |
| 18 | 0.115 | 2.111 | 0.492 | 1.157 | 1.825 | 1.400 |
| 20 | 0.067 | 2.176 | 0.496 | 1.187 | 1.832 | 1.415 |
| 22 | 0.038 | 2.222 | 0.498 | 1.205 | 1.844 | 1.431 |

- **Positive correlation.** The micro values are positively correlated with the macro strips. Beyond the converse
  pairs, they hit the strip residuals about 1.4 times as often as independent values with the same law.
- **Drift.** The excess still drifts up by about 0.016 per two levels at n = 18..22. Boundedness is therefore
  plausible but not settled by these data.
- **Consequence for a proof.** A proof of (DC) must tolerate a positive correlation of this size.

## 5. HYP-9482: Poisson-pair second moments (pre-registered, outcomes)

- **Model.**
  - On Y_n, g(v)/2^n tends in law to XY, with X and Y independent and uniform on [1/2, 1). This is PROVED: Session 16,
    `BINARY_PRODUCT_DENSITY.md` section 3. Only the Poisson clustering at the unit scale is heuristic.
  - Its density is f_Y(t) = 4 ln(4t) on [1/4, 1/2] and 4 ln(1/t) on [1/2, 1].
  - Converse classes are taken to be Poisson with mean mu = f_Y/4 per odd value.
- **Closed forms.** The model then gives:
  - E^Y/|Y| -> 2 + (12 - 16 ln 2) = 14 - 16 ln 2;
  - the fraction of |Y| that is values with q classes: 2[(ln 2)^(q+1)/(4(q+1)!) + int_0^(ln 2) e^(-2u) u^q/q! du];
  - the smoothed slope -> 12 - 16 ln 2;
  - on X_n with independent layers: E/|X| -> 2 + 4(1 - ln 2) = 6 - 4 ln 2, and the cross-layer part -> 12 ln 2 - 8.
- **Outcomes at n = 22, 24, 26.** Registered in 21b02648ec before any n >= 21 run.
  - P1-P5, P6(a) and P8 CONFIRMED. P5 is now a theorem (THM-4890 (C)).
  - P1 and P8 are not independent of HYP-9479's previously published E/|X| values.
  - P6(b) refuted at the boundary level a = n/2 - 3. P7 refuted: maximum multiplicity 18 at n = 26, where the model
    expects 0.13 such values.
  - Both misses are MISTAKE-728: tolerances not checked against the committed in-sample data, and a tail bound set by
    eye.
  - The full table is in HYP-9482.
- **Addendum P9-P11 (n = 28, pre-registered 25fde96277).** CONFIRMED: E^Y/|Y| = 2.92480, an excess of 0.0152, in the
  registered window and decreasing.
- **Conjecture (OPEN).** E^Y_n/|Y_n| -> 14 - 16 ln 2 and E_n/|X_n| -> 6 - 4 ln 2. The excesses are decreasing:
  - E^Y: 0.032, 0.025, 0.020, 0.015 (n = 22..28);
  - E/|X|: 0.033, 0.030, 0.025, 0.020.

## 6. The x+2 / x*2 duality and the gate, in this setting

- **Two lanes.** The repository's two natural-number modes are odd lanes by x+2 and dyadic columns by x*2 (TANGENTS
  T1169; HYP-1890).
  - The one-apex family is a map from the dyadic tree of words (x*2: v -> 2v, 2v+1) to the odd lane (x+2) of counts.
  - Linear energy is exactly the statement that this map has bounded L^2 multiplicity.
  - The Poisson-pair law says its multiplicity is Poisson(ln 2) in converse classes, so it covers half the odd lane.
  - The rest of the spectrum must come from other families, which fill every odd value except 7 and 21 below 2^35
    (S82 addendum, one run).
- **Order recursions on tournaments** (Session 210, Session 211).
  - The gate (n -> n + 2) turns a tournament into a strong one with HC(Gamma(T)) = H(T) and
    H(Gamma(T)) = H(T) + P_T(2).
  - Skew-Sylvester doubling (n -> 2n) is a different, diagonal operation.
- **The gate on a one-apex tournament.** Gamma(T_v) is the one-apex word 0 v 1 with its extreme chain arc reversed:
  - b becomes a new chain source that beats the apex (a trailing 1);
  - a becomes a new chain sink beaten by the apex (a leading 0);
  - the arc between them is a -> b.

  So HC(Gamma(T_v)) = g(v). On the energy problem the gate is the identity (a shift of order by 2 that preserves the
  count), while x*2 is the recursion the energy is about.

## 7. THM-4891: face-local relaxations of parallel coherence

| n | square-coherent | hexagon-coherent | tournaments |
|---|---|---|---|
| 3 | {0..3} | {1, 3} | {1, 3} |
| 4 | {0..6}, acyclic {1..6} | {0..9, 12}, acyclic {1..9, 12} | {1, 3, 5} |
| 5 | {0..25} | {0..54, 56, 60} | {1, 3, 5, 9, 11, 13, 15} |
| 6 | {0..90} (SAT, every k) | not computed | 19 values; 7, 21, 35, 39 missing |

- **Acyclic witnesses.** Explicit orientations with 7 and 21 sinks exist at n = 5 under each relaxation, and with 7
  sinks at n = 4 under hexagons (`permutohedron/witnesses.out`). Each is re-verified by an independent rebuild of the
  oriented Cayley graph.
- **The n = 4 hexagon witness with 7 sinks.** Format: pair | relative order of the other two -> element first at sinks.

      01|23 0   01|32 1   02|13 2   02|31 0   03|12 0   03|21 0
      12|03 1   12|30 2   13|02 3   13|20 3   23|01 2   23|10 3

  - Its sinks are 0123, 0231, 0312, 1032, 1203, 3102 and 3210.
  - Four of the six pairs (01, 02, 12, 23) flip with the order of the others; a tournament would have none flip.
  - All 48 seven-sink orientations at n = 4 have exactly four flipping pairs (`audit/aud_perm_flips.out`).
- **Reading.**
  - Session 210 showed that tournaments sit between "convex-linear" (one sink) and "acyclic" (all counts).
  - The face-local relaxations sit between tournaments and "acyclic", and already lose 7 and 21 at the first order
    where those counts fit.
  - Parity is lost too: 0 and even counts occur.
  - Neither face type alone keeps 7 and 21 out. Both together give tournaments. Intermediate relaxations, such as all
    hexagons plus some squares, are OPEN.

## 8. {3, 10}, {15, 25} and {7, 21}

The concurrent lanes typed these maps; nothing here changes their verdicts.
- **{3, 10}.** These are (H - 1)/2 for {7, 21}, and also the number of colored disjoint odd-cycle packings modulo swap
  (Session 210).
- **{15, 25}.** These are Busch floors bracketing 21 (THM-1370). They are also square-sum onsets, with p(Q23) = 3 and
  p(Q25) = 10 chains (Session 211). Session 211 proved that no transport (H - 1)/2 = p(Q) can exist.
- **One-apex reading.** In the one-apex family, 3 and 10 are among the index misses 3, 6, 10, 14, 17, 20, 21, 24, ...,
  a set of density 1/2 (HYP-2199). The Poisson-pair law predicts that density as e^(-ln 2). So {3, 10} is not a
  one-apex phenomenon. It is what survives after all families are united.

## 9. Next obligations

1. **Prove (DC)** (THM-4890 (D)).
   - It says the micro self-energies g(l) on m ~ n/2 - 2 bits correlate with the macro strips only boundedly: at most an
     absolute constant times the random-model count.
   - The observed constant is about 3.7 against the uniform model, about |X_n| of it converse-forced.
   - This is the whole content of linear energy, hence positive density.
   - Natural tools are digital exponential sums (bit reversal) and Beatty-sequence equidistribution. They must allow a
     positive correlation.
2. **Done (audit U1).** The sharp balance |A - S| <= 1 is proved. What remains open is whether the sign rule of (C) or
   the 1...1 / 0...0 window shape can be used to count collisions.
3. **Test the energy-constant conjecture further.** It passed n = 28 (P9). The excess decays by about 0.78 per two
   levels; test n = 30, 32 out of sample: E^Y/|Y| -> 14 - 16 ln 2.
4. **Permutohedron.**
   - Find f(n) for odd n >= 7 (the square-coherent maximum).
   - Find the hexagon-coherent spectrum at n = 6.
   - Decide whether any relaxation strictly between face-local and zone-wide coherence keeps 7 and 21 out, for example
     all hexagons plus some squares, or coherence on contexts of bounded size.

## 10. Reproduction and audit

- **Energy and collision statistics** (numba):
  - `python s84_energy.py 14 20 layers > insample_14_20.out` (in sample);
  - `python s84_energy.py 21 26 layers > outofsample_21_26.out` (out of sample; under 1 minute, 1 GB);
  - `python s84_energy.py 28 28 layers > outofsample_28.out` (addendum; 135 s, about 6 GB).
  - Model numbers: `python s84_model.py`.
- **THM-4890 (A'), (C) checks:** `python check_balance.py`.
- **Criterion (DC) measurement (exploratory):** `python dc_ratio.py 12 14 16 18` and `python dc_ratio.py 20 22`
  (72 s at n = 22).
- **Tilted energies (exploratory):** `python tilted_energy_sup.py 10 32` and `python tilted_energy_sup.py 12 24`.
- **THM-4891** (in `permutohedron/`):
  - `python perm_relax.py` (exhaustive n = 3, 4);
  - `python perm_relax2.py n4` (acyclicity);
  - `python perm_relax2.py sat L1 0 60` and `python perm_relax2.py sat L2 0 60` (n = 5);
  - `python perm_witness.py` (witnesses and the independent check);
  - `python perm_relax3.py 6 L1 any 90 atleast` (f(6) >= 90).
  - `sat_n6_L1_spectrum.out`: n = 6 squares, every k = 0..90.
- **Audit:** `audit/AUDIT.md`. PASS WITH CORRECTIONS:
  - MAJOR-1: "decorrelation" re-read as bounded correlation;
  - MAJOR-2: "zone-wide coherence needed" restricted to "neither face type alone";
  - 13 minors: credits to Session 16 and Session 188, the rename to set-context tournaments, pre-registration
    disclosure, ratio rounding, the flip count, witness tables, the archived tilted-energy run, wording;
  - upgrades: U1 sharp balance lemma, U2 checks to n = 28, U3 exact random-model expectation, U4 independent
    permutohedron verification.
  All are applied, and MISTAKE-730 records the two majors. Build and run instructions are in the audit's section 6.
