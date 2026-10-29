# 101, 19, 79, 169: the 3x+101 cycle is the rotational tournament on seven vertices — a circulant-tournament ↔ rational-Collatz dictionary, the locking theorem behind 35 and 175, and multiplier descent for Paley run observers

**Session:** mac-mini-2026-10-29-cycle101, 2026-10-29.

**Owner's directive (verbatim):** "think of prior work regarding 101, 19, 79, 169 and really try to find many examples
of sessions on each and merge in related ideas they discuss as you investigate those too and look for higher order
connections, keep focusing on tournament and collatz structure and symmetry. two cycles sharing denominator \(101\) and
identical clocks acquire denominators \(35\) and \(175\) after the same appended step."

**Status.**
- *PROVED (elementary):*
  - THM-4849: dictionary, two-base-point lemma and composite-q structure, attachment law in fixed-point form, locking
    theorem, interval corollary, 5x+1 criterion;
  - THM-4848: multiplier descent, multiplicity law, converse conjugation.
- *FINITE-EXACT:* every table here (scripts in `04-computation/experiments/cycle101_tournament_collatz_20261029/`).
- *CITED:* Böhm–Sontacchi, Lagarias, Turner/Babai, Dvoretzky–Motzkin, Ailon–Rudnick, Bugeaud–Corvaja–Zannier.
- *OPEN:* HYP-9446. One pre-registered prediction FAILED (HYP-9446 (1)) and one passed (HYP-9446 (3)).
- *Audit.* PRE-AUDIT checkpoint; independent audit pending (section 8).

---

## 0. The answer in one screen

1. **The four numbers are one object.** 19, 79, 169 are the three points of the rational 3x+1 cycle with word
   (1,1,5) and denominator `101 = 2^7 − 3^3`, which is the integer cycle 19 → 79 → 169 of `n ↦ (3n+101)/2^v`.
   - Read as a 7-periodic parity word, that cycle is the out-neighbourhood `{1,2,3}` of the **rotational tournament
     R7**: THM-4849 (A), the dictionary.
   - 19 is the marking at the first letter, `19 = 3^3 − 2^3`.
2. **Every 101-cycle is a tournament.** Every rational cycle of shape (3,7) comes from a circulant tournament on
   `Z_7`. The five 3x+101 cycles are exactly R7 (three achiral cycles) and Paley P7 (a chiral pair swapped by
   reversal = the converse). Seven is the last prime where this happens: tournament cycles number `2^k − k`, all
   cycles `Cat_k`.
3. **The owner's sentence.**
   - The two cycles are the same tournament R7 in two labellings: the interval labelling (1,1,5) and its ×3
     decimation (1,3,3).
   - Appending the letter (1) cancels a prime `p` of the new denominator iff the marked point is `p`-adically close
     to `−1`, the fixed point of `x ↦ (3x+1)/2`.
   - The interval cycle is **locked** at 5: all its letters are ≡ 1 mod `ord_5(2) = 4`, so all three points are
     ≡ −1 (mod 5). Hence 35 at every marking. The ×3 labelling is unlocked, and 31/101 ≡ 1 (mod 5): 175.
   - The locking theorem is an iff (THM-4849 (E)). In general the interval cycle of `R_q` is locked under (1)
     exactly at the primes of `gcd(2^m − 1, 3^m − 1)`, `m = (q+1)/2`; at `q = 7` this is `gcd(15, 80) = 5`.
4. **The same pair is the pair of positive 5x+1 three-cycles.** The positive 5x+1 cycles {13,33,83} and {17,43,27}
   have exactly these two words (collatz_fermat_reframes C1). On the (3,7) shape, a word gives an integral 5x+1
   cycle iff all its letters are odd.
5. **Symmetry on the run-observer side (THM-4848).**
   - The multiplier group of a circulant tournament, which acts by decimation on the Collatz side, acts on run
     observers by Galois descent.
   - For Paley, `det M_k = det B_0 · N_(Q(√−q))(β)^((q−1)/2)`, and Smith multiplicities are divisible by
     `|QR·⟨p⟩|`. In P11 the primes 89 and 62987 occur as `(Z/p)^5`, while in R11 such primes occur once.
   - The converse conjugates `β`, so the "one cubic, zero of its dual" pattern of the overnight packet is a
     chirality signal.

---

## 1. The four numbers in prior work (merged sweeps)

The digests below were made by four read-only sweeps of the repository. They have been merged and cut down; only
items used or typed here are kept.

### 1a. 101

- **The origin of the cycle (opus, 2026-09-21; `collatz_guards_20260921_squarefree.md`).**
  - The actual 3x+1 path 23 → 35 → 53 → 5 has word (1,1,5) and reset map `R(n) = (27n+19)/128`, with slope gap
    `128 − 27 = 101`. Its fixed point 19/101 has the phases 19/101, 79/101, 169/101.
  - Repetition with all `3r+1` nodes squarefree is possible iff `r ≤ 4985`. 101 is the only local obstruction, and
    13 is the one exceptional prime because `169 = 13²`.
  - **This is where the owner's four numbers first appear together, and 35 is on the path.**
- **Free shapes (THM-4484, procgen sporadic 2026-09-26).** The shift criterion forces all five (3,7) necklaces to be
  3x+101 cycles, but 101 is never named there. THM-4849 (C) is that statement, with the tournament reading added.
- **The 5x+1 transport (Session100 and collatz_fermat_reframes, today).** The 5x+1 cycle 13 → 33 → 83 decodes to
  19/101, and denominator 101 excludes its basin. The C1 table adds (1,3,3) ↔ 17, 43, 27 and (1,4) ↔ 1, 3.
- **Two roles of 101 in one algebra.** 101 is the denominator `2^7 − 3^3` and also the carry `C(1,2,2,·)`, with reset
  `(81x+101)/512` and centre 101/431 (October 11–19 notes). The cycle numerator ignores the last letter.
- **On the tournament side.**
  - THM-4746 (Camion pseudoprimes, mac-mini fanopsp): `(2/101) = −1`, sporadic zero class 19 (`101 | s(19) − 2^9`),
    and `12019 = 7·17·101`.
  - ORBITS(8) = 101 in the metagraph atlas (THM-853/854).
  - 101 of the 456 order-7 classes reach zero ear chi (HYP-9080).
  - P19, P23 and P31 Hamilton ranks certified mod 101 (session54).
- **Other.** `g(101) = ⌊101²/3⌋ + 2` (THM-4820, MISTAKE-687); 101 negative discriminants of exponent ≤ 2
  (THM-4566); `ord_101(2) = 100`, today's carry clock.

### 1b. 19

- **The carry 19 arises in exactly two shapes** (19-sweep observation, verified):
  - `19 = 3 + 2^4 = C(4, ·)`, giving the cycle {5, 19} of 3x+23. This is T5's interval cycle (THM-4849 (C), q = 5).
  - `19 = 9 + 6 + 4 = C(1,1,·) = 3^3 − 2^3`, giving the ladder `19/(2^q − 27)`: −1, −19/11, 19/5, 19/37, 19/101,
    19/229 for `q = 3, …, 8`. Earlier sessions met single rungs of this ladder (codex 2026-10-04, posets, run
    classification).
- **The strongest earlier Collatz–tournament bridge** (`checked_switch_phase19_20261004.md`, difference families,
  2026-10-04).
  - The branches `C_a(x) = 2^(−a)(3x+1)` mod 19 act on `P_19` and preserve orientation iff `a` is odd, since 2 and 3
    are non-squares mod 19.
  - `19 = 3^3 − 2^3 = N(3 − 2ω)`.
  - This acts on residues mod `p`. THM-4849 reads words as connection sets; the two are different mechanisms.
- **The P_19 tournament line.**
  - P19 is not H-maximal (THM-135; the interval circulant wins). The Betti formula is THM-130. `|Aut P19| = 171`.
  - P19 is 3 flips from R19 because `h(−19) = 1` (THM-4751). `tc(P19)` fails Korselt (THM-4747). `399 = 3·7·19` is a
    Camion pseudoprime (THM-4746).
  - THM-4848 here adds `β_2(P19) ∈ Q(√−19)`, the ninth-power law, and a one-sided 5.
- **The 18s.** `ord_19(2) = 18`, so `(Z/27)^* ≅ F_19^*` (mac-mini mod1819, THM-4554). The golden notes typed most
  "18" coincidences as NUMEROLOGY.

### 1c. 79

- **On the Collatz side.**
  - `79 = C(1,5,·) = 3² + 3·2 + 2^6` is the numerator of every 3-letter word starting (1,5): 79/101 for (1,5,1) and
    79/229 for (1,5,2) (79-sweep observation, verified).
  - The plain path 79 → 119 → 179 → 269 → 101 → 19 (difference families, golden; typed NUMEROLOGY there).
  - 79 is the least element of one 9-cycle of 7x+169 (THM-4516): `{79, 361 = 19², 337}`, since
    `7·79 + 169 = 2·19²`.
- **The doubling tournament.** `ord_79(2) = 39`, so `⟨2⟩ = QR_79`: the doubling clock mod 79 *is* `P_79` (THM-4631,
  opus S30). This is the bridge note's reading at a Mersenne-free prime; see THM-4849 (F) for how the two readings
  meet.
- **The 79 in Paley7's three-run group** (THREE_RUN_GLUE R12, CLOCK_MINIMAL_PAIRING C5–C8). `79³` divides the order;
  the module is one cubic `X³+13X²+12X−1` and none of its dual.
  - The sweep found that the prior work explains *behaviour* (Frobenius orbits), not *which* primes occur, and that
    no `Q(√−7)` norm formula existed. THM-4848 supplies it: `N_(Q(√−7))(β) = 2^22·5^4·11^2·79·5861^2·392809`.
  - The sweep's untested question, whether reversing P7 swaps the occupied cubic, is answered **yes** by THM-4848 (v).
  - **That this 79 is the cycle's 79 is NUMEROLOGY.** It is a prime splitting in `Q(√−7)`, and nothing links it to
    `C(1,5,·)`.
- **Other.** Golden Paley primes 11, 19, 31, 59, 71, 79 (THM-4656). `h(−79) = 5` and the beaten natural order of
  `P_79` (HYP-9287). Waring: 79 needs 19 fourth powers (THM-4824 note, typed numerology).

### 1d. 169

- **On the Collatz side.**
  - `169 = C(5,1,1) = 3² + 3·2^5 + 2^6` (THM-4849 (G)). It is the third point of the interval cycle.
  - It is the only prime square obstructing the 19/101 orbit (2026-09-21).
  - `169 = 2^9 − 7^3` is a Fermat–Catalan clock (THM-4515/4516), with the Eisenstein square `8 − 7ζ = (3 − ζ)²`.
  - 3x+169 has fixed points 13 and 169 and a 9-cycle through 101 (today, `cycles_3xb.out`).
- **On the tournament side.**
  - 169 passes the strong Camion test: 13 is a Camion–Wieferich prime, `13² ∥ D(13)` (THM-4746/4750, HYP-9358).
    Rechecked today (`camion_cycle.out`); the products `19·79·169`, `19·101` and so on fail.
  - `det(C + J) = 169` for the regular cyclic tournament on `Z_13` (THM-2532).
  - `H = 169` is a gap at `n = 7` and first appears at `n = 8` (THM-4810).
  - No Paley tournament exists on `F_169` (`169 ≡ 1 mod 4`).
- **That 169 is a square is NUMEROLOGY for the cycle.** The square matters only for the 13-adic squarefree
  obstruction of 2026-09-21.

### 1e. 35 and 175 (the owner's sentence)

- **Sources.** collatz_fermat_reframes C3 gives the test 185/175 = 37/35 against 221/175. session101 (concurrent,
  read before writing) gives the local centre law, the fibre `{7:1, 25:1, 35:5, 175:8}`, the 175/175/7 rotation
  witness and the R7/Paley7 tournament controls.
- **Coincidences.** `H(R7) = 175 = 2^8 − 3^4` was already flagged by session101 as not identified with the appended
  denominator. 35 is a coefficient of R7's characteristic polynomial (session101 controls), also not identified.
- **The 3x+35 cycle.** 13/35 lies on the same 3x+35 cycle {13, 37, 73, 127} as 37/35 (101-sweep): the shape (4,8),
  denominator `2^8 − 3^4 = 175` reduced. GENERIC.

---

## 2. The dictionary (THM-4849)

The Terras map `T` acts on rationals with odd denominator. A circulant tournament `T(S)` on `Z_q` (`q = 2k+1` prime,
`S ⊔ −S = Z_q^*`) gives the `q`-periodic parity word `1_S`, hence a rational `T`-cycle of length `q` with `k` odd
points and denominator dividing `2^q − 3^k`.

| tournament side | Collatz side |
|---|---|
| vertex `v` | cycle point `z_v` with parity vector `1_(v+S)` |
| translation `v ↦ v − 1` | the Collatz map `T` |
| the in-neighbours of a vertex | the odd points (the `U`-cycle) |
| gaps of `S` | the valuation word |
| multiplier `x ↦ mx` (isomorphism) | decimation of the word (moves between cycles with the same denominator) |
| converse `T(−S)` | word reversal (preserves the trace, five-mirrors Prop. 6) |
| connection set of `T(S)` | tournament base point: backwards = complement of forwards |
| rotational `R_q` | the only cycles with two base points = the achiral ones |

**The two readings.**
- *The bridge note's reading* (2026-10-01) reads a word as a binary fraction. There `×2` is the Collatz shift and
  translations die.
- *Here* translations are the shift and multipliers move between cycles.
- *Where they meet.* At the Paley base points the binary values satisfy the Gauss-period equation
  `c² + c + (q+1)/4 ≡ 0 (mod 2^q − 1)`.

**Counts** (`circulant_collatz.out`).

| `q` | `D_q = 2^q − 3^k` | tournament cycles `2^k − k` | all cycles `Cat_k` | classes (decimation orbits) |
|---:|---|---:|---:|---|
| 5 | 23 | 2 | 2 | T5 (the regular 5-tournament) |
| 7 | 101 | 5 | 5 | R7 (3 cycles), P7 (2) |
| 11 | 1805 = 5·19² | 27 | 42 | R11 (5), two classes of 10, P11 (2) |
| 13 | 7463 = 17·439 | 58 | 132 | R13 (6), four of 12, the ⟨3⟩-cyclotomic class (4) |

**At composite `q`.** For `q = 9, …, 27`, doubly-based cycles are lexicographic products `C[R_(q/g)]` (THM-4849
(B')). No cycle with three base points was found (HYP-9446 (4)).

## 3. The 35 and the 175 (THM-4849 (D), (E))

**Appending `v` to the marked point `x = C/D`.** The new numerator is

    C(wv) = 2^H (1 − (2^v − 3)x) + x·D(wv),

so `p^e` cancels iff `x ≡ ξ_v = 1/(2^v − 3) (mod p^e)`.

**The centres are the one-letter cycles.** `ξ_1 = −1`, `ξ_2 = 1` (the trivial cycle), `ξ_3 = 1/5`, `ξ_4 = 1/13`.

**For (1): `D(w·1) = 175 = 5²·7`.**
- The six (3,7) markings ≡ −1 (mod 5) are 19, 29, 49, 79, 89, 169; 49 is also ≡ −1 (mod 25), giving 7.
- The one marking ≡ −1 (mod 7) is 53, giving 25.
- This reproduces session101's fibre exactly.

**Locking theorem** (iff; 26665 tests, all agree). Every marking is ≡ `ξ_v (mod p^e)` iff every letter is
≡ `v` (mod `ord_(p^e)(2)`).

**The interval cycle.** For `R_q`'s interval cycle under (1) the locked primes are those of `gcd(2^m − 1, 3^m − 1)`,
equivalently `lcm(ord_p 2, ord_p 3) | m`.
- Locked at 5 iff `q ≡ 7 (mod 8)`, at 7 iff `q ≡ 11 (mod 12)`, and so on.
- Table to `q < 100` in `locking_iff.out`: 5 at `q = 7`; 7 at 11; 11 at 19; 5, 7, 13 at 23; …

**Locking does not see the tournament class.** The pre-registered prediction "locked tournament cycles are
rotational" held for `q ≤ 19` and FAILED at `q = 23`, where 14 of 18 locked cycles are not rotational
(HYP-9446 (1)).

## 4. Multiplier descent and chirality (THM-4848)

**The decomposition.**
- The blocks `B_j` of `M_k(T(S))` over the characters of `Z_q` satisfy `σ_m(B_1) ~ B_1` (monomial conjugation) for
  every multiplier `m ∈ H_S`.
- So `β = det B_1 ∈ K_H` exactly, `det M_k = det B_0·N_(K_H/Q)(β)^|H|`, and each `Z/p^e` of the nontrivial cokernel
  has multiplicity divisible by `|H⟨p⟩|`.

**Data.**

| case | quadratic norm of `β` | the decisive check |
|---|---|---|
| P7, k = 3 | `2^22·5^4·11^2·79·5861^2·392809` | the overnight packet's primes, now with a norm formula |
| P11, k = 2 | `2^12·5^2·11^2·31·53·89·62987·…` | **89, 62987 ≡ 1 (mod 11) occur as (Z/p)^5** |
| R7, R11 controls (`H = 1`) | — | primes ≡ 1 (mod q) occur **once** |
| P19, k = 2 | 310 digits: `2^8·5^2·11·47·61·263·c299` | invariance under `σ_4`, ninth-power law |
| `Z_13`, `H = ⟨3⟩` | quartic field | invariance; the law holds |

**Chirality.**
- *The converse conjugates.* `β(T(NQR)) = β̄(T(QR))` (THM-4848 (v)), so the occupied prime of `Q(√−q)` flips under
  the converse. The overnight packet's "one cubic, zero of its dual" (C5) is therefore an orientation invariant.
- *Odd primes look one-sided.* At every split odd prime tested, `β` is one-sided (HYP-9446 (3), pre-registered,
  passed four non-automatic instances).
- *The prime 2 is two-sided.* 2 is two-sided in P7.

## 5. Higher-order connections, typed

| connection | type | why |
|---|---|---|
| {19, 79, 169}/101 = the interval cycle of R7 | STRUCTURAL | THM-4849 (A), (G) |
| all five 101-cycles = all circulant tournaments on Z_7; 7 is the last such prime | STRUCTURAL | (B): `2^k − k = Cat_k` iff `k ≤ 3` |
| C3 pair = two labellings of R7 | STRUCTURAL | decimation ×3 |
| 35 at every marking of {19,79,169} | STRUCTURAL | locking at 5, `gcd(2^4 − 1, 3^4 − 1) = 5` |
| C3 pair = the two positive 5x+1 three-cycles | STRUCTURAL (cited transport plus a new criterion) | reframes C1, plus all letters odd ⟺ `3 \| C_5` |
| T5's interval cycle {5, 19} of 3x+23 = the 5x+1 cycle {1, 3} | STRUCTURAL (cited) | reframes C1 table |
| Paley pair = the chiral pair, equal traces 197/101 | STRUCTURAL | reversal = converse; five-mirrors Prop. 6 |
| "one cubic, zero of its dual" flips under the converse | STRUCTURAL | THM-4848 (v) |
| Paley binary values = Gauss periods mod `2^q − 1` | STRUCTURAL (classical in substance) | THM-4849 (F) |
| 19 = 3³ − 2³ and 19 = 3 + 2⁴ (on R7 and T5) | NUMEROLOGY | `27 = 8 + 16 + 3` |
| the 79 of Paley's three-run group vs the cycle's 79 | NUMEROLOGY | a split prime of `Q(√−7)` vs `C(1,5,·)` |
| 169 = 13², a Camion pseudoprime | NUMEROLOGY for the cycle | a 13-adic fact (THM-4746) |
| `H(R7) = 175` = appended denominator | NUMEROLOGY | flagged by session101 |
| minima 19, 23, 29, 31, 37 consecutive primes | NUMEROLOGY | — |
| `19² \| D_11`; P11's cycle has reduced denominator 361; the cycle {5,7,17,31,35} of 3x+19 is a Z_11 tournament cycle | GENERIC | divisibility of `C` by `19²` only |
| Paley pair's 5-integrality at q = 11 (QR) and 19 (NQR) | GENERIC | consistent with a 1/5 chance per cycle |
| interlocks: 101 periodic for 3x+79 and 3x+169; 175's cycles through 101 and 169; 7x+169 through `19²`; the plain orbit of 79 through 101, 19, 13 | GENERIC / NUMEROLOGY | small integers on small cycles (`cycles_3xb.out`) |
| inversion pairs `p/q`, `q/p` both periodic | GENERIC | 192 mutual pairs among 2805 periodic reduced rationals with `p, q ≤ 401`, every odd `b ≤ 401`, starts `≤ 2·10^6` (`periodic_rationals.txt`) |

## 6. Cycle spectra (FINITE-EXACT, every odd start ≤ 2·10^7, `cycles_3xb.out`)

- **3x+101.** Eight cycles:
  - the five (3,7) tournament cycles;
  - two (6,14) cycles, {7,61,71,157,143,265} and {11,67,151,277,233,25}, with `D = 2^14 − 3^6 = 5·31·101` reduced
    to 101;
  - the fixed point 101.
- **3x+23.** The two T5 cycles {5, 19} and {7, 11}, a 26-cycle, and the fixed point 23.
- **3x+19.** The Z_11 tournament cycle {5, 7, 17, 31, 35} and the fixed point 19.
- **3x+175.** 18 cycles, including {101, 239, 223, 211}, an 8-cycle through 169, and the fixed point 35.
- **3x+35.** {13, 37, 73, 127}, {17, 43, 41, 79}, {25, 55}, and six more.
- **3x+169.** 12 cycles: fixed points 13 and 169, and a 9-cycle through 101.
- **3x+79.** Five cycles, the cycle through 1 containing 101 and 13.
- **Camion** (`camion_cycle.out`). Only 169 passes the strong test. Every product of 19, 79, 101 and 169 fails, even
  though `19·101` satisfies the 101-local condition (`r = 19`, THM-4746's sporadic class).

## 7. Open (HYP-9446)

1. REFUTED: locked tournament cycles are rotational.
2. Ailon–Rudnick type: infinitely many primes `q` with no interval-locked prime.
3. Paley one-sidedness at split odd primes (pre-registered, passed so far; P23 next).
4. At most two base points for every odd `q`.

## 8. Reproduction

| what | command | output |
|---|---|---|
| dictionary | `python3 circulant_collatz.py` | `out/circulant_collatz.out`, `.json` |
| (3,7) family and futures | `python3 family37.py` | `out/family37.out`, `.json` |
| fixed-point law (8055 tests) | `python3 attachment_lemma.py` | `out/attachment_lemma.out` |
| locking theorem, interval table | `python3 locking_iff.py` | `out/locking_iff.out` |
| HYP-9446 (1), (2) and Paley 5-integrality | `python3 locked_tournaments.py` | `out/locked_tournaments.out` |
| composite base points | `python3 composite_basepoints.py` | `out/composite_basepoints.out` |
| Paley Fourier, P7 k = 3 | `python3 p7_fourier.py` | `out/p7_fourier.json` |
| Paley Fourier, P7 and P11 k = 2 | `python3 paley_fourier.py` | `out/paley_fourier.out`, `.json` |
| P19 | `python3 paley_fourier2.py 19; python3 p19_norm.py` | `out/paley_fourier2*.json/.out`, `out/p19_norm.out` |
| Smith multiplicity law | `python3 galois_smith.py` | `out/galois_smith.out`, `.json`, `.log` |
| chirality, converse | `python3 chirality.py; python3 converse_conjugate.py` | `out/chirality.out`, `.json`; `out/converse_conjugate.out` |
| Camion tests | `python3 camion_cycle.py` | `out/camion_cycle.out` |
| cycle spectra | `cc -O2 -o cycles_3xb cycles_3xb.c; ./cycles_3xb b 20000000` | `out/cycles_3xb.out` |
| inversion pairs | `cc -O2 -o inversion_pairs inversion_pairs.c; ./inversion_pairs 401 2000000 401` | `out/periodic_rationals.txt` |

The binaries are not committed. The scripts import helpers from `../growth_gf_runs_fermat_20261029/` and call PARI/gp.

## 9. Audit

PRE-AUDIT checkpoint. An independent audit follows.
