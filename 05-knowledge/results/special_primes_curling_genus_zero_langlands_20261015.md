# Special primes, curling numbers and genus zero: the owner's list read through the repository's objects

2026-10-15, session opus-2026-10-15-S37. Companion to [the energy-budget note](collatz_octave_memory_bias_energy_budgets_20261015.md). Audited by the independent S37 referee (`04-computation/experiments/s37_audit_20261015/REPORT.md`, scripts `a7_special_primes.py`, `a7b_orbit40487.py`, `a8_curling_check.py`); its corrections (I12) are applied below.

**Status.**
- **PROVED, elementary:**
  - non-generous primes are exactly the primes whose least primitive root is a Wieferich base, with the corresponding orbit explosion of circulant tournaments;
  - the first odd ancestors of 1 under `3x − 1` are the numbers `(2^m + 1)/3` with `m` odd, and the Wagstaff numbers are the subfamily with `m` prime;
  - Woodall and Cullen numbers shadow the 2-adic points −1 and 1;
  - `2^y + y²` prime forces `y = 1` or `y ≡ 3 (mod 6)`;
  - the Paley half-imbalance reading of supersingularity for `p ≡ 3 (mod 4)` (a reformulation of classical formulas).
- **FINITE-EXACT:**
  - the non-generous primes in `(2, 10^7)`;
  - the Leyland primes below `10^12`;
  - the supersingular characterization for `p < 2000`;
  - the Stern primes below `10^6`;
  - the curling records for all starts of length `<= 25`, reproducing OEIS A094004.
- **CITED:** Ogg, Conway–Norton, Borcherds; Deuring–Eichler; Dirichlet; Chaffin–Linderman–Sloane–Wilks; Kapustin–Witten; Gaitsgory–Raskin et al.; Feigin–Frenkel–Reshetikhin; Frenkel.
- **ANALOGY / NUMEROLOGY:** typed where they occur.

## 1. Non-generous primes are the canonical-base orbit explosions

**Fact.** A prime `p` is non-generous (OEIS A055578) when its least primitive root `g` is not a primitive root of `p²`. Since `g` has order `p − 1` modulo `p`, its order modulo `p²` is `p − 1` or `p(p − 1)`. So `p` is non-generous iff `g^(p−1) ≡ 1 (mod p²)`: `g` is a Wieferich base at `p`.

**Corollary (the THM-4644 / THM-4656 orbit law with base `g`; PROVED by the same proof).** Let `d'` be the odd part of `p − 1` and `g' = g^((p−1)/d')`. The circulant tournaments on `Z_(p²)` invariant under `g'` number `2^E`, where:
- `E = (p − 1)/d'` if `p` is generous;
- `E = (p + 1)(p − 1)/(2d')` if `p` is non-generous.

*Proof.* `g'` has odd order `d' p^ε` on units, with `ε = 0` exactly in the non-generous case, and order `d'` on the nonzero non-units. Negation pairs the orbits freely; count orbit pairs as in THM-4644 part 2. ∎

**Data (FINITE-EXACT).**
- The only non-generous prime in `(2, 10^7)` is `40487`, with `g = 5`.
- The next known one, `6692367337`, is confirmed (`g = 5`, `5^(p−1) ≡ 1 (mod p²)`).
- At `p = 40487` the invariant circulant count jumps from `2^2` to `2^40488`. At `6692367337` it jumps from `2^8` to `2^26769469352`.
- Both known odd non-generous primes are base-5 Wieferich primes (OEIS A123692) whose least primitive root is 5.
- `p = 2` is non-generous by convention.

**Credit.** THM-4644 (opus S33) and THM-4656 (mac-mini) for the orbit law with bases 2 and `φ`.

## 2. Wagstaff, Woodall, Cullen: first ancestors and 2-adic shadows

- **Wagstaff (PROVED).** An odd `n` has `3n − 1 = 2^m` iff `m` is odd and `n = (2^m + 1)/3`.
  - So the first odd ancestors of 1 under the `3x − 1` map (`(3n − 1)/2 = 2^(m−1) → … → 1`) are exactly the numbers `(2^m + 1)/3` with `m` odd, just as `(4^k − 1)/3` are under `3x + 1`.
  - The Wagstaff numbers `(2^p + 1)/3`, `p` prime, form a strict subfamily: `m = 1, 9, 15, 21, 25, 27` give the non-Wagstaff ancestors `1, 171, 10923, 699051, 11184811, 44739243` (correction by the S37 referee; the draft said "exactly").
  - Wagstaff primes `(2^p + 1)/3` for prime `p < 1300`: `p = 3, 5, 7, 11, 13, 17, 19, 23, 31, 43, 61, 79, 101, 127, 167, 191, 199, 313, 347, 701` (probable-prime tests).
  - The "new Mersenne conjecture" of Bateman–Selfridge–Wagstaff links them to Mersenne primes (CITED).
  - Repository context: the `3x − 1` versus `3x + 1` lanes of S32 and S36.
- **Woodall and Cullen (PROVED; checked `2 <= n <= 80` by the referee).** Under `T(x) = x/2, (3x+1)/2`:
  - `W_n = n·2^n − 1` begins with exactly `v_2(W_n + 1) = n + v_2(n)` consecutive odd steps (the draft said `n`), and after `n` of them it sits at `n·3^n − 1`. It follows the itinerary of the 2-adic fixed point −1 of `x ↦ (3x+1)/2`, like the Mersenne numbers `2^n − 1`.
  - `C_n = n·2^n + 1` follows the trivial-cycle word `(OE)^j` through `3^j n 2^(n−2j) + 1` for `2j <= n − 1`. That is the 2-adic shadow of the fixed point 1 of `x ↦ (3x+1)/4`.
  - So Woodall and Cullen numbers are the simplest "shadow" families of S6's shadow theorem, sitting over the two simplest 2-adic cycle points.

## 3. Leyland primes on the alphabet {2, y}

- **PROVED.** `2^y + y²` can be prime only for `y = 1` or `y ≡ 3 (mod 6)`. For even `y` the value is even. For odd `y` with `3 ∤ y` it is `2 + 1 ≡ 0 (mod 3)`.
- **FINITE-EXACT (probable-prime tests).** For `y < 400` the primes occur at `y = 1, 3, 9, 15, 21, 33`.
- Every Leyland prime `x^y + y^x` below `10^12` is of this form: `17 = 2^3 + 3^2`, `593 = 2^9 + 9^2`, `32993`, `2097593`, `8589935681 = 2^33 + 33^2`.
- `17` is the unique Leyland number on the Collatz pair `{2, 3}`.

## 4. Supersingular primes, Paley tournaments and the fifteen rational function fields

**Classical (CITED).**
- Ogg (1974): `X_0(p)^+` has genus 0 exactly for `p ∈ {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71}`, which are the primes dividing the order of the Monster (moonshine: Conway–Norton 1979; Borcherds 1992).
- For `p >= 5` this is equivalent to: every supersingular `j`-invariant in characteristic `p` lies in `F_p`.
- The total number of supersingular `j` is `⌊p/12⌋ + ε_p` (`ε_p = 0, 1, 1, 2` for `p ≡ 1, 5, 7, 11 (mod 12)`).
- The `F_p`-rational ones number `(h(−p) + h(−4p))/2` for `p ≡ 3 (mod 4)` and `h(−4p)/2` for `p ≡ 1 (mod 4)` (Deuring; Eichler).

**Check (FINITE-EXACT, `p < 2000`, class numbers by reduced forms).** The primes `5 <= p < 2000` with all supersingular `j` rational are exactly `5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71`.

**Paley reading (exact reformulation).** For `p ≡ 3 (mod 4)`, `p > 3`, Dirichlet's class number formula gives `h(−p) = I(p)/(2 − (2/p))`, where

`I(p) = Σ_(0<a<p/2) (a/p) = 2·#{out-neighbours of 0 among 1, …, (p−1)/2 in P_p} − (p−1)/2`

is the half-interval score imbalance of a vertex of the Paley tournament.
- Hence the number of `F_p`-rational supersingular `j` is `I(p)` for `p ≡ 7 (mod 8)` and `(2/3) I(p)` for `p ≡ 3 (mod 8)`.
- The eight supersingular primes `≡ 3 (mod 4)` (`7, 11, 19, 23, 31, 47, 59, 71`) are exactly the primes `p ≡ 3 (mod 4)` at which this quantity saturates the total count `⌊p/12⌋ + ε_p`. For example `I(71) = 7 = #ss(71)` and `(2/3)I(59) = 6 = #ss(59)`.
- Non-examples: `43` gives 2 against 4, and `79` gives 5 against 7.
- Checked for every `p ≡ 3 (mod 4)` below 200.

**The rational function field.** Genus 0 means the function field of `X_0(p)^+` is a rational function field `C(T_p)`, generated by a Hauptmodul: the McKay–Thompson series of the Monster class `pA`. So the fifteen supersingular primes are the fifteen levels at which this modular curve is a Riemann sphere and its function field is literally the field of rational functions on it.

**Coincidences (NUMEROLOGY).**
- `17` is supersingular, a Stern prime, a Leyland prime and a Fermat prime.
- The golden Paley primes up to 71 (`11, 19, 31, 59, 71`; THM-4656 part 4) are all supersingular, while the next, `79`, is not.
- Small primes are in many lists; no structure is claimed.

## 5. Stern primes

- **Data (FINITE-EXACT, below `10^6`).** The primes not of the form `q + 2a²` (`q` prime, `a >= 1`) are `2, 3, 17, 137, 227, 977, 1187, 1493` (OEIS A042978). The odd composites below 20001 not of the form `q + 2a²` (`a >= 0`) are `5777` and `5993` (Stern 1856, against Goldbach's "other conjecture").
- **Connection.** The quadratic form `2a²` is half of the Pell norm `x² − 2y²` behind the Q–Pell unit `3 + 2 sqrt(2)`. No structural link to the repository objects is claimed (NUMEROLOGY).

## 6. Curling numbers, Gijswijt's sequence and a 2-3 termination conjecture

**Definitions.**
- The curling number of a word is the largest `k` with `word = X Y^k`, `Y` nonempty.
- Gijswijt's sequence (A090822) starts from 1 and repeatedly appends the curling number of the word so far.
- The **Curling Number Conjecture** (Chaffin–Linderman–Sloane–Wilks, J. Integer Seq. 16 (2013)) states that every finite starting sequence, extended this way, eventually produces a 1. The case of starts over `{2, 3}` is the essential one, and it remains open.

**Results (FINITE-EXACT; `curling_number_stopping_times_20261015.py`).**
- Incremental curling uses `O(1)` work per period per appended term, and is checked against brute force on random starts.
- Gijswijt's first `4` is at position 220.
- For every one of the `2^n` starts over `{2, 3}` with `n <= 25`, the extension reaches 1.
- The longest extensions, in appended terms up to and including the 1, are `1, 3, 3, 5, 5, 9, 9, 59, 60, 61, 113, 113, 113, 119, 119, 119, 119, 119, 120, 120, 120, 121, 121, 121, 121`.
- In OEIS A094004's convention (length before the terminating 1) these are `1, 4, 5, 8, 9, 14, 15, 66, 68, 70, 123, …, 145`, which match the published values term by term.
- The mean extension stabilizes at 3.74.

**Collatz analogy (ANALOGY).** The curling map is a deterministic process over the same two letters as the Collatz parity alphabet. Its termination conjecture has the same shape: every start reaches the absorbing "1". Its stopping-time records grow extremely slowly, consistent with Gijswijt's tower-type first appearances (`4` at 220; `5` near `10^(10^23)`). No transfer of methods is claimed.

## 7. Geometric Langlands, S-duality and the Riemann sphere

**The established picture (CITED).**
- Kapustin and Witten (Commun. Number Theory Phys. 1 (2007)) derive geometric Langlands from S-duality of four-dimensional `N = 4` super-Yang–Mills. The S-transformation exchanges the gauge group `G` with its Langlands dual `G^∨`.
- After topological twisting and compactification on a Riemann surface `C`, S-duality becomes mirror symmetry of Hitchin moduli spaces (T-duality on the Hitchin fibres).
- 't Hooft line operators act as Hecke operators. Wilson operators give the eigenvalues, the dual local systems, and eigenbranes give Hecke eigensheaves.
- The unramified de Rham geometric Langlands conjecture was proved in 2024 by Gaitsgory, Raskin and collaborators.
- In genus 0, `C = P¹` has no nontrivial local systems, so the content comes from ramification. For tame ramification at marked points, the Gaudin model realizes the correspondence (Feigin–Frenkel–Reshetikhin 1994; Frenkel 2004): Gaudin eigenvalues correspond to opers on `P¹` with regular singularities. Those are differential operators whose coefficients lie in the rational function field `C(z)`.

So the field of rational functions on the Riemann sphere is where the dual side lives in genus 0. In the function-field (arithmetic) Langlands programme, `F_q(T)` plays the same role.

**Bridges to the repository (ANALOGY; each states what matches and what does not).**
1. **Hecke-type relations.** The exact class-sum recursion `C(A; N) = −C(A/2; ⌊N/2⌋) − C((3A+1)/2; 3⌈N/2⌉)` (THM-4658) is a correspondence between the 2-adic and 3-adic "places" acting on residue-class functions. Haar measure is an eigen-distribution precisely at the zeros of `D(s)`, so `D` plays the role of a Hecke polynomial.
   - *What does not match:* there is no group, no automorphic side and no dual group.
2. **Critical lines.** The 3-adic transfer operator is balanced in `L²(Z_3)` exactly on `Re s = 1/2` and `Re s = 1` (THM-4658, statement 5). That is formally like the unitarity line of tempered spectra (Ramanujan bounds). The square-root residue biases proposed in HYP-9282 did not survive their out-of-sample test, so nothing is known to sit on that line.
   - *What does not match:* no spectral theorem or Hilbert-space unitarity is proved; the operator norms on `L²(Z_3)` are about 1.35-1.40, not 1 (S37 referee).
3. **Rational functions.**
   - S29's rational approximants of the Collatz transfer law are elements of `C(w)`.
   - The Q covariance update is a Möbius map, an element of `PGL_2 = Aut(C(w))`.
   - The supersingular genus-0 curves of section 4 have rational function fields.

   These are three distinct appearances of `C(z)`, not one structure.
   - *Added after the out-of-sample census (same session).* The parity sum `S(N)` is fitted to 2% three octaves ahead by an explicit formula over the zeros of `D(s) = 1 + 2^(−s) + (3/2)^s/3` (HYP-9285), the analogue of the Riemann–Weil explicit formula with `D` in the role of the zeta function and `H` as the remainder. Along a convergent `p/q` of `log_2 3`, `D` is approximated by the Laurent polynomial `1 + z^q + z^(q−p)/3` in `z = 2^(−s/q)`, with error controlled by `|s| ln 2 |log_2 3 − p/q| ≪ 1` (0.0095 for `s_0` and the convergent 84/53). The modal expansion is then a limit of finite, rational-function-field-like explicit formulas, which is the S29 rational-approximant picture (HYP-9260) read on the zeros instead of the poles (ANALOGY).
4. **Electric–magnetic exchange.** The repository's mirror pairs, such as `3x + 1 ↔ 3x − 1` and tournament reversal `T ↔ T^R` (HYP-9280's in- and out-degree duality), are involutive dualities. S-duality acts on the complexified coupling `τ` through `SL(2, Z)`: `S: τ ↦ −1/τ` exchanges electric and magnetic charges (and `G` with `G^∨`), and `S²` is charge conjugation. Only the involutive exchange shape matches (ANALOGY, weak).

## 8. Reproduction

```sh
python 04-computation/experiments/special_primes_collatz_tournament_bridges_20261015.py
python 04-computation/experiments/curling_number_stopping_times_20261015.py 25 3000
```

Outputs: `05-knowledge/results/special_primes_collatz_tournament_bridges_20261015.out`, `05-knowledge/results/curling_number_stopping_times_20261015.out`.

**References.**
- Number theory:
  - A. Ogg, Automorphismes de courbes modulaires, Sém. Delange–Pisot–Poitou 16 (1974–75);
  - J. H. Conway and S. P. Norton, Monstrous moonshine, Bull. LMS 11 (1979);
  - R. Borcherds, Invent. Math. 109 (1992);
  - M. Deuring (1941); M. Eichler (1955) on supersingular counts;
  - OEIS A055578, A123692, A042978, A094004, A090822, A094133, A064539.
- Curling numbers: B. Chaffin, J. P. Linderman, N. J. A. Sloane and A. R. Wilks, On curling numbers of integer sequences, J. Integer Seq. 16 (2013), Article 13.4.3; F. J. van de Bult, D. C. Gijswijt, J. P. Linderman, N. J. A. Sloane and A. R. Wilks, A slow-growing sequence defined by an unusual recurrence, J. Integer Seq. 10 (2007).
- Langlands and S-duality:
  - A. Kapustin and E. Witten, Electric-magnetic duality and the geometric Langlands program, Commun. Number Theory Phys. 1 (2007);
  - D. Gaitsgory, S. Raskin et al., proof of the geometric Langlands conjecture (2024);
  - B. Feigin, E. Frenkel and N. Reshetikhin, Comm. Math. Phys. 166 (1994);
  - E. Frenkel, Gaudin model and opers (2004).
- Repository: THM-4644, THM-4656, THM-4658, HYP-9280, HYP-9282, S6 shadow theorem, S29 (HYP-9260).
