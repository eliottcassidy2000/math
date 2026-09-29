# Three preprints on twisted Bernoulli zeros as structural mirrors of the 3-adic Syracuse law: the character spectrum of the law (two exact identities tying the seed-1 mass and the negative-cycle spike to sums of multiplicative moments), the full-period recursion (the primitive Fourier maximum over all units to level 19), the Poisson census of the 3-adic zero lines behind the ridges, and the saddle chain of 1729 = 1 + 12^3 with its basin

**Session:** opus, `collatz-poset-dag-20260927` (S23), 2026-09-29.
**Owner's directive:** "dissect the attached preprints for captivating ideas you
can apply creatively as they merge with our work and consider extending
previous work related to the fact there are no integers `x, y, z` with `xyz ≠
0`, `gcd(x, y) = 1` and `x^4 + y^3 = z^17`, and how 17 relates creatively with
other past work in the repo, also 1729 and the proof that at least `X^0.84`
numbers return to 1729 in Collatz; be free." The three preprints are Peter
Chocian, *Explicit twisted Hilbert class components beyond classical
irregularity* (arXiv:2607.23177, 25 Jul 2026), *Character Fourier spectra of
circular units and twisted Bernoulli class components* (arXiv:2607.27503, 29
Jul 2026), *Twisted Bernoulli zeros in quasi-linear time* (arXiv:2608.08724,
9 Aug 2026). Per the standing pattern for fresh papers (memory: structural
analogy, the shape of the hard direction), each is read for the shape of its
argument and matched to the thread's frontier.
**Parallel session, complemented not duplicated.** The owner's `(4,3,17)`
statement, the Krasikov–Lagarias `X^0.84` bound applied to any root, the
`3x-1` sheet and the necklace picture were worked the same day by the
mac-mini session `collatz-necklace-20260929`
([note](collatz_necklace_20260929_fair_splits_power_clocks_basins.md); THM-4515,
THM-4516, THM-4517, HYP-9165; the statement recorded UNVERIFIED there). This
note cites that work, adds an independent audit of its three theorems
(section 6), and works the parts it left: the preprints, the character
spectrum, the census, the saddle chain and the basin of `1729`.
**Inherits (cited):** the S19 Theorems A–C of
[`mazur_positive_density_20260928.md`](mazur_positive_density_20260928.md)
(harmonic mass `H_n(y) = 3^n mu_n(y)`, reference density `rho_n = (2/3) H_n`,
the negative-cycle spike `rho_n(-1) ≈ 0.975 (3/2)^n`, the seed-1 test); the
S20–S22 Fourier profile of
[`collatz_five_mirrors_20260929.md`](collatz_five_mirrors_20260929.md)
(the closed family `m_n(k) = mu_hat_n(2^k mod 3^n)`, its recursion, the
primitive maxima `M(h)` to `h = 18` by FFT, the ridges of the negative family,
the multiplier families, Theorem C and its hypothesis H); Proposition 4 of
[`collatz_hedgehog_family_20260927.md`](collatz_hedgehog_family_20260927.md)
(`T^(2j)(1 + 4^j t) = 1 + 3^j t`); Theorem 2.3 of
[`collatz_mod6_20260917_extended_collatz_scc.md`](collatz_mod6_20260917_extended_collatz_scc.md)
(the 3-adic hostile family, whose `j = 6` row lists `730, 973, 1297, 1729,
2305, 3073, 1024`); THM-4484 (the sporadic `-17` cycle, clock `2^11 - 3^7 =
-139`); Krasikov–Lagarias 2003 (CITED, as in
[`arithmetic_braids2_20260917_inverse_completion.md`](arithmetic_braids2_20260917_inverse_completion.md)).

**Status: PROVED (Theorem 1, the character-spectrum duality and the two
increment identities; Proposition 2, the full-period one-pole recursion is
exact; the saddle-chain identities; the critical-rate identity `e^(-I)
3^(θ*/ln 2) = log_2 3 - 1`; the lower bound `X_0(1729) > 2^30` for the
Krasikov–Lagarias threshold) + VERIFIED (the character spectrum to level 16
against S19's independent values, five digits for `rho_n(1)` at every level `2..16` and the printed three to four digits for the parity differences; the
primitive maxima over all units to level 18 against the S20 FFT to seven
digits, extended to level 19) + FINITE-EXACT (the zero-line census on `16.7
million` pairs; the basins of the saddle chain, of `27` and of `137` to
`2^30`, confirmed by a second code at `2^21`) + OBSERVED (the multiplicative
spectrum is non-Rayleigh with a heavy tail and maxima growing like `n/2`
Parseval units; the basin densities) + NUMEROLOGY (typed as such: the `139`
at `973`, the `137` at the merge of `1729` and `27`, the exponents `4, 3`
of the saddle) + UNVERIFIED (the owner's `(4,3,17)` statement: a web search
on 2026-09-29 found no source; the same search found a September 2026
preprint claiming the smallest open Beal signature `(3,5,7)`). No Collatz
proof step. Audit of mac-mini's THM-4515–4517: section 6. Own audit (a second
subagent, `collatz_three_mirrors_20260929_audit.py/.out/.md`): SOUND WITH
CORRECTIONS, eighteen applied, record in section 8.**

Scripts and outputs: `04-computation/experiments/collatz_three_mirrors_character_spectrum_20260929.py`
→ `05-knowledge/results/collatz_three_mirrors_character_spectrum_20260929.out`
(the full-period recursion with FFT, the character spectrum, the identities,
levels `<= 16`); `collatz_three_mirrors_fullperiod_max_20260929.py` →
`..._fullperiod_max_...out` (the lean recursion, maxima over all units to
level 19); `collatz_three_mirrors_zero_lines_20260929.py` and
`..._zero_lines_wide_...py` → `.out` (the 3-adic zero-line census);
`collatz_three_mirrors_saddle_chain_20260929.py` → `.out` (1729, 17, the
chain, the merge with 27); `collatz_three_mirrors_basins_20260929.c` →
`..._basins_...out` (the basins to `2^30`); `collatz_three_mirrors_basins_mask_20260929.c` → `..._basins_mask_...out` (the `1 + 12^m` family, `13, 17, 137, 9232`); `collatz_three_mirrors_jacobi_transfer_20260929.py` → `..._jacobi_transfer_...out` (the transfer recursion and its structure to level 7).

---

## 0. What is new, in one screen

1. **The character spectrum of the Syracuse law is the Gauss-sum dual of
   the exponent profile, and two of its sums are old friends (Theorem 1,
   PROVED; VERIFIED).** For a character `psi` of `(Z/3^n)^×` and the
   multiplicative moment `S_n(psi) := E[psi(Y_n)]`,
   `sum_(k mod L_n) mu_hat_n(2^k) conj(psi(2^k)) = tau(conj psi) S_n(psi)`
   (`L_n = 2·3^(n-1)`, `psi` primitive), and, summing over the primitive
   characters, `sum_psi S_n(psi) = rho_n(1) - rho_(n-1)(1)` and `sum_psi
   psi(-1) S_n(psi) = rho_n(-1) - rho_(n-1)(-1)`: the seed-1 reference
   density of S19 Theorem C is the accumulated sum of the multiplicative
   moments, and the negative-cycle spike of Theorem B is their accumulated
   parity-twisted sum. The seed-1 identity reproduces S19's independently computed `rho_n(1)` to five digits at every level `2..16` (`rho_16(1) = 0.29011` both ways); the parity identity reproduces S19's printed differences to the three or four digits S19 prints (`rho_16(-1) - rho_15(-1) = 213.4` there, `213.43` here). The
   multiplicative spectrum has the Parseval rms `0.84·3^(-n/2)` but is not
   Rayleigh: its largest moments grow like `n/2` Parseval units (`9.1` at `n
   = 16`, where random phases would give `3.5`), a heavy tail (`P(z > 3) =
   0.074` against `0.050`), and at low levels sit on the characters
   `psi_(±2), psi_(±4), psi_(±8)` — the 3-adic-logarithm characters.
2. **The full-period recursion (Proposition 2, PROVED; VERIFIED).** Because
   `2` generates the units, the closed recursion of S20 can be run on the
   whole cycle `Z/L_n` as a one-pole filter `m_n(k) = (g(k-1) + m_n(k-1))/2`
   with no valuation truncation, in `O(L_n)` time and memory: it returns
   every Fourier coefficient of `mu_n` at every unit, organised by the
   discrete logarithm. It reproduces the S20 FFT maxima to seven digits at
   every level `1..18` (`0.0111873` at `18`) and extends them: `M(19) = 0.0098157` at `k = 24` (`= 19 log_2 3 - 6.1`) and its mirror — the maximum
   over all units is still at `±2^s`, `s = h log_2 3 - 6 ± 1` (`s = 20, 21,
   23, 24, 26` at `h = 16..20`; `floor(h log_2 3) - 6` is exact at `h = 19`
   only).
3. **The 3-adic zero lines behind the ridges are Poisson (FINITE-EXACT).**
   Chocian's survey design (a `p`-adic coincidence per pair, counted against
   a Poisson null model, with a depth) applied to the ridge seeds of S22: the
   pairs `(u, Q)` with `u·2^Q ≡ ∓1 mod 3^d`, depth `d = v_3(u 2^Q ∓ 1)`. On
   `16.7·10^6` pairs (`u <= 5000` odd prime to 3, `Q <= 10^4`) the depth
   histogram is geometric to `1%` through depth 10 (forced by the
   equidistribution of the powers of two for `d <= 8`, where `L_d <= W`) and
   the tail is Poisson-consistent at every depth `>= 10` (`14, 5, 2, 1` lines of depth `>= 14,
   15, 16, 17` against `10.5, 3.5, 1.2, 0.4`). The S22 seed `(55, 423)` is
   one of five depth-15 lines; the deepest is `1187·2^5031 ≡ 1 mod 3^17`. A
   small box (`666,000` pairs) had shown a `1%`-level excess of deep lines,
   which the wide box dissolves. The identity `e^(-I) 3^(θ*/ln 2) = log_2 3 - 1` (the definition of `I`
   rewritten) reads the critical rate of Theorem C as the Chernoff amplitude
   law of the multiplier families extended to `u ≈ 3^n` under `M(n) ≈
   e^(-nI)` (CONJECTURAL); on that reading Parseval would force the
   multiplier exponent to steepen beyond `-0.438` for large `u`, and H's
   margin `0.577` against `0.585` is the margin between the Parseval scale
   and that extrapolation (DIRECTION).
4. **1729 is the balanced point of a saddle chain (EXACT), and its basin
   (FINITE-EXACT).** `1729 = 1 + 12^3 = 1 + 3^3 4^3` is the midpoint `j = 3`
   of the chain `1 + 3^j 4^(6-j)` (`4097, 3073, 2305, 1729, 1297, 973, 730`)
   that the orbit of `1 + 4^6` runs under `T^2` (Proposition 4 of the
   hedgehog note); `1 + 12^m` is the midpoint of the chain of `1 + 4^(2m)`
   for every `m`; `17 = 1 + 4^2` heads the chain `17, 13, 10`. Along the
   chain `1297 = 6^4 + 1` is prime, `973 = 7·139` carries the clock of the
   `-17` cycle (NUMEROLOGY), the chain bottom `730 = 3^6 + 1` gives the
   second taxicab representation `1729 = 3^6 + 10^3`, and the orbits of
   `1729` and `27` merge at `137` (NUMEROLOGY). Basins to `2^30`:
   `dens B(1729) = 0.004395`, nested along the chain from `0.003339`
   (`4097`) to `0.014459` (`730`); `dens B(137) = 0.2992` (three in ten
   integers pass through `137`); `B(27)` is its doubling ray. The
   Krasikov–Lagarias count `|B(1729) ∩ [1, X]|` is `12%` of `X^0.84` at `X =
   2^30`, so the theorem's threshold `X_0(1729)` exceeds `2^30` (PROVED by
   the count) and is near `2^49` if the density persists (OBSERVED).
5. **The `(4,3,17)` statement** stays UNVERIFIED: a web search found no
   source for `x^4 + y^3 = z^17`; it found arXiv:2609.26996 (September 2026,
   "The primitive generalized Fermat equation `x^3 + y^5 = z^7`: a
   computer-assisted proof"), which, if correct, removes the smallest open
   Beal signature named in mac-mini's note (AUTHOR-CLAIMED). mac-mini's sentence that the clock shadow of `x^4 + y^3 = z^17` is empty was refuted as unsupported by the audit of section 6 (the shadow — the solutions whose even variable is a power of two — is as open as the statement), as was the Beal-shadow consequence of THM-4516(5); both are corrected in their files and logged in MISTAKE-551.

---

## 1. The three preprints and their shape

**What they do.** For an odd primitive Dirichlet character `chi` of conductor
`f` and a prime `p` with `p ∤ f`, the divided generalized Bernoulli values
`b_(chi,j) = f B_(1, chi ω^(-j)) mod p` control, through the characterwise
Main Conjecture, the odd isotypic components of the `p`-class group of
`Q(zeta_(fp))`. The first paper (conductor `3`, `p ≡ 1 mod 6`) finds that
the local expansion of the reflected circular unit `mu_p = (1 + z zeta_p)/(1
+ conj(z) zeta_p)`, `z = -zeta_3^2`, selects exactly the zeros of the
`chi_(-3)`-twisted Bernoulli values: with universal polynomials `P_m` (`sum_m
P_m(X) Y^m = -log(1 - X(1 - e^(-Y)))`), `P_(p-j)(h) - P_(p-j)(1-h) =
-(2h-1)(j-1)! b_j`; the zeros ("blind lines") below `500` are twelve, at `p
= 67, 103, 139, 199, 241, 271, 331, 337 (twice), 409, 421, 457`, seven of
them at classically regular primes, and each carries an explicit Kummer
radical generating an order-`p` Hilbert class component, certified by a
finite split-prime computation. The second paper isolates the mechanism as a
character Fourier transform, `sum_t conj(chi)(t) P_m(h_t) = tau(conj chi)
B_(m,chi)/(m·m!)` (odd `m`), `h_t = zeta_f^t/(zeta_f^t - 1)` — a Gauss sum performs the
transform — and applies it at conductor `5` (eleven lines). The third
computes the whole spectrum `(b_(chi,j))_j` in quasi-linear time by a
residue-class weight formula and a Bluestein chirp factorisation, surveys
`27,508` zero lines over `55,121` character–prime pairs, and finds the
counts consistent with a Poisson model (`chi^2 = 0.05` and `3.93` at
conductors `3` and `5`), the divided digits uniform, conjugate quartic
characters never vanishing at a common index, no association with classical
irregularity, and eight non-simple zeros including one of depth three at
`(f, p, j) = (19, 37, 16)`.

**Their shape, and what it mirrors here.**

* *A universal object's local expansion, Fourier-transformed against the
  characters of a cyclic group by a Gauss sum, is a special value.* The
  thread's universal object is the closed family `m_n(k) = mu_hat_n(2^k mod
  3^n)`, a function on the cyclic unit group through the discrete logarithm;
  its character Fourier transform is, by the same Gauss-sum step, the
  multiplicative moment `E[psi(Y_n)]` of the Syracuse law (Theorem 1). The
  "special values" that appear are the reference densities at `±1` — the
  seed-1 mass of Theorem C and the negative-cycle spike of Theorem B — as
  telescoping sums over the levels (section 2).
* *A twist reveals degeneracies invisible to the untwisted invariant.* Seven
  of the twelve `chi_(-3)`-twisted lines sit at classically regular primes.
  The S22 ridges are the same phenomenon in the thread: the multiplier
  families `u·2^j` (the "twists" of the pure family) resonate at 3-adic
  coincidences `u 2^Q ≡ ∓1 mod 3^(n_0)` that the pure family does not see.
* *A `p`-adic coincidence per pair, counted against a Poisson null model,
  with a depth and a second-order digit.* This is exactly the design the
  ridge seeds needed; section 3 runs it.
* *A reflection `psi = ω θ^(-1)` pairing character lines, with an
  antisymmetric functional equation whose antisymmetric part is the special
  value.* The thread's reflection is word reversal with the carry reciprocity
  `C_(rev w)(u,v) = u^(d-1) v^A C'_w(1/u, 1/v)` (S20, Proposition 5) and the
  reversal-invariant trace (Proposition 6); the antisymmetric part of the
  carry under reversal, evaluated at `(3,2)`, is the difference of the cycle
  points of a word and its reverse (`-17` against `-13801/139` for the
  sporadic cycle). Typed ANALOGY; not pursued further here.
* *A finite certificate for a global statement (a split-prime Artin
  calculation proving a radical is not a `p`-th power).* The thread's
  finite certificates are the orbit checks of the Fermat–Catalan clocks
  (THM-4516) and the basin counts of section 4: finite computations that
  decide a stated instance, never the general statement.

---

## 2. The character spectrum of the Syracuse law

**Setting.** `Y_n` is Tao's Syracuse random variable at level `n` (law `mu_n`
on `Z/3^n`, supported on the units); `L_n = 2·3^(n-1)`; the family `m_n(k) =
mu_hat_n(2^k mod 3^n)` is periodic in `k` with period `L_n` since `2`
generates `(Z/3^n)^×`. The characters of `(Z/3^n)^×` are `psi_j(2^k) =
e(jk/L_n)`, `j mod L_n`; `psi_j` is primitive (conductor exactly `3^n`) iff
`3 ∤ j` for `n >= 2`, and `psi_j(-1) = (-1)^j`. `S_n(psi) := E[psi(Y_n)] =
sum_y mu_n(y) psi(y)`.

**Theorem 1 (PROVED).** (i) For every `j`,
`F_n(j) := sum_(k mod L_n) m_n(k) e(-jk/L_n) = sum_y mu_n(y) g(conj psi_j, y)`,
`g(conj psi, y) = sum_(u unit) conj psi(u) e(uy/3^n)`; for `psi_j` primitive
and `y` a unit, `g(conj psi_j, y) = psi_j(y) tau(conj psi_j)`, so
`F_n(j) = tau(conj psi_j) S_n(psi_j)`, `|tau| = 3^(n/2)`. (ii) Summing over
the primitive characters (orthogonality of characters on the two levels),
`sum_(psi prim) S_n(psi) = L_n mu_n(1) - (L_n/3) mu_(n-1)(1) = rho_n(1) -
rho_(n-1)(1)` and `sum_(psi prim) psi(-1) S_n(psi) = rho_n(-1) - rho_(n-1)(-1)`
for `n >= 2`, with `rho_n(y) = (2/3) 3^n mu_n(y)` the reference density of
S19 (for `n = 1` the sums over the non-trivial character are `-1/3` and
`+1/3`). Hence `rho_n(±1) = 1 ∓ 1/3 + sum_(m=2)^n sum_(psi prim mod 3^m)
(±1-twisted) S_m(psi)`. Proof: (i) substitute `u = 2^k` in the sum over `k`
and use the Gauss-sum identity `sum_u conj psi(u) e(uy/3^n) = psi(y)
tau(conj psi)` for units `y` (`Y_n` is always a unit); (ii) `sum_(psi prim)
psi(y) = L_n 1_(y ≡ 1 mod 3^n) - (L_n/3) 1_(y ≡ 1 mod 3^(n-1))`, and `mu_n(y ≡
1 mod 3^(n-1)) = mu_(n-1)(1)` by consistency; `L_n mu_n(1) = (2/3) 3^n
mu_n(1)` since `L_n = (2/3) 3^n`. ∎ (Corollary: `E[chi_(-3)(Y_n)] = -1/3`
exactly, `Y_n ≡ 2 mod 3` having probability `2/3`.)

**Proposition 2 (the full-period recursion; PROVED).** With `omega_n(k) =
e((2^k mod 3^n)/3^n)` and `g_n(k) = omega_n(k) m_(n-1)(k mod L_(n-1))`, the
S20 recursion `m_n(k) = sum_(a>=1) 2^(-a) g_n(k-a)` on `Z/L_n` is the
one-pole filter `m_n(k) = (g_n(k-1) + m_n(k-1))/2` around the cycle, whose
unique periodic solution is obtained by running the filter from any state for
`>= 80` steps and then once around (the transient is `2^(-steps)`). No
valuation is truncated. The residues `2^k mod 3^n` are generated in blocks by
int64 arithmetic for `n <= 19`. Cost `O(L_n)` per level; memory two vectors.
The first implementation ran two passes for every cycle and inherited a
`2^(-2)` transient from level 1 (all maxima `12%` low) — the check against
the S20 FFT caught it; the corrected values agree with the FFT at every
level `1..18`: `0.5773503, 0.3779236, ..., 0.0144095 (16), 0.0125107 (17),
0.0111873 (18)` (fullperiod output; seven digits), and at level 19 `M(19) = 0.0098157` at `k = 24` (`19 log_2 3 - 6.1`; ratio `M(19)/M(18) = 0.877`; independently reproduced by the audit, mass `0.7089`) with the mirror `-2^24` equal to seven digits and `k = 25` (`0.009769`), `k = 23` next: the maximum over all `774,840,978` units is on `±2^s`
with `s = floor(19 log_2 3) - 6`, extending the S20 observation by one level; and in complex64 (`37 GB`, `_fullperiod_max20_` output) `M(20) = 0.0088846` at `k = 26` (`20 log_2 3 - 5.7`; ratio `0.905`) and its mirror, `k = 25` next (`0.008380`): the maximum over all `2,324,522,934` units at level 20 is on `±2^26`. Fourier mass per level `0.709` at
`n = 19` (this is `3^n/L_n = 3/2` times S20's `0.462 -> 0.472`: the same sum
`sum_(units) |mu_hat_n|^2`, normalised per unit here and as a total there;
it equals the rms² of the character spectrum).

**The spectrum (VERIFIED to level 16, character-spectrum output).** Per
level, over the primitive characters: the rms of `|S| 3^(n/2)` is `0.8452,
0.8321, 0.8345, 0.8356, 0.8362, 0.8356, 0.8360, 0.8367, 0.8373, 0.8380,
0.8387, 0.8393, 0.8399, 0.8404, 0.8408` (`n = 2..16`) — the Parseval scale,
`rms^2 = 0.70` = the typical `|mu_hat|^2 3^n` of S20, as it must be (the
Gauss sum has modulus `3^(n/2)`). The distribution of `z = |S|^2/mean` is
not exponential: `P(z > 1, 2, 3) = 0.27, 0.13, 0.074` against `0.37, 0.14,
0.05` for random phases — fewer moments above the mean, a heavier tail. The
largest `|S| 3^(n/2)` grows steadily: `1.00, 1.41, 1.59, 2.02, 2.21, 2.82,
2.95, 3.60, 4.09, 4.35, 5.75, 6.04, 7.41, 7.79, 9.11` for `n = 2..16`
(roughly `n/2`), where the maximum of `L_n` independent Rayleigh variables
would be `0.84 sqrt(ln L_n) = 3.5` at `n = 16`; the top four conjugate pairs `psi_(±j)`
at `n = 16` are all `8.9–9.1`, a cluster, not an outlier. Which characters: at `n = 3,
4, 6` the largest is `psi_(±2)`, at `n = 5` `psi_(±8)`, at `n = 4` `psi_(±4)`
next, at `n = 7` `psi_(±80)`, at `n = 8` `psi_(±278) = psi_(∓2^12)` (`278 = L_8
- 4096 = 2(3^7 - 2^11)`); from `n = 9` on the indices are not small and not
powers of two (`±521, ±1193, ±1127, ±4219, ±40162, ±71333, ±154112 =
±2^9·301, ±1516613`). The character `psi_(±2)(y) = e(±2 log_2(y)/L_n) =
e(±log_2(y)/3^(n-1))` is the 3-adic-logarithm character of the 1-unit
component of `y` (`log_2(y) = 2 log_4(±y)` on `(Z/3^n)^× = {±1} × (1 +
3Z)/(1 + 3^n Z)`); since `Y_n = 2^(-a_n)(1 + 3z)` with `z = sum_(i>=1)
3^(i-1) 2^(-(a_(n-1) + ... + a_(n-i)))`, `psi_(±2)(Y_n) = e(∓a_n/3^(n-1))
e(±2 log_4(1 + 3z)/3^(n-1))`, a phase product over the suffix sums of the walk
read 3-adically — the exponent walk of S22 with the roles of `2` and `3`
exchanged (DIRECTION; the dominance of these characters at low levels is
OBSERVED, and it does not persist past `n = 8`). The parity split: `mean_odd
S = -mean_even S` to three digits at every level (`∓1.116·10^(-5)` at `n =
16`), which is Theorem 1(ii): the difference `sum_even - sum_odd = rho_n(-1)
- rho_(n-1)(-1)` grows like `0.325 (3/2)^n` while the sum `rho_n(1) -
rho_(n-1)(1)` is `O(0.05)`.

**The Jacobi transfer (exact recursion PROVED; its structure VERIFIED to level 7, CONJECTURED beyond; `_jacobi_transfer_` output).** Since `Y_n = 2^(-a)(3 Y_(n-1) + 1)` and `3y + 1 mod 3^n` depends on `y mod 3^(n-1)`,
`E[psi(Y_n)] = G_psi · sum_(psi' mod 3^(n-1)) c(psi, psi') E[psi'(Y_(n-1))]`, `G_psi = sum_(a>=1) 2^(-a) psi(2)^(-a) = (psi(2)^(-1)/2)/(1 - psi(2)^(-1)/2)`, `c(psi, psi') = (1/L_(n-1)) sum_(y unit mod 3^(n-1)) psi(3y + 1) conj psi'(y)`
— an exact linear recursion on the multiplicative spectrum (checked to `10^(-13)` at `n <= 7` against the spectra of the exact law). Its structure, found by computation: for `psi` primitive mod `3^n` and `psi'` primitive mod `3^(n-1)` (`n >= 3`) the Jacobi-type coefficients have **constant modulus** `|c(psi, psi')| = (#prim')^(-1/2)` (`0.5, 0.289, 0.167, 0.096, 0.056` at `n = 3..7`, equal to `1.000000` after scaling for every one of the pairs), and `c = 0` for imprimitive `psi'` (`< 2·10^(-13)`); the geometric factor satisfies `1/3 <= |G_psi| <= 1`, `rms |G_psi| = 3^(-1/2)` exactly over the primitive characters (`n >= 4`), with `|G_psi| -> 1` iff `psi(2) -> 1`, i.e. for the characters of small index `j` (`|G| = 0.950, 0.994, 0.9993, 0.9999` at `j = ±2`, `n = 4..7`). So each level is a *flat unitary-like mixing* of the previous spectrum (gain `3` in `ell^2`, the column sums of `|c|^2`) followed by the *diagonal contraction* `G_psi` (mean square `1/3`): the Parseval rms is preserved, and the shape of the spectrum is that of a Rayleigh variable multiplied by a factor ranging over `[1/3, 1]` — which is the observed non-Rayleigh form (fewer moments above the mean, a heavier tail) — with the largest moments on the characters where `|G_psi| ≈ 1`: `corr(|S_n|, |G_psi|) = 0.77, 0.74, 0.68, 0.65, 0.62` at `n = 3..7`, and the mean of `|S_n| 3^(n/2)` is `1.19–1.35` where `|G| > 0.9` against `0.48–0.50` where `|G| < 0.45`. This is Proposition 3 of S20 (no uniform one-step gap of the geometric Gauss sums) in the multiplicative picture: the contraction `|G_psi| = |sum_a 2^(-a) psi(2)^(-a)|` has no gap at the characters with `psi(2) ≈ 1`, the 3-adic-logarithm characters, which is why `psi_(±2), psi_(±8)` lead at low levels and why the sum of the spectrum (the seed-1 mass) is carried by the near-trivial characters. The constant-modulus law is the prime-power analogue of the classical `|J(chi_1, chi_2)| = p^(n/2)` for Jacobi sums of primitive characters; a proof for this restricted sum (over `y` a unit, with `psi'` evaluated at `y` and `psi` at `1 + 3y`) is left as an obligation.

**The identities, checked (VERIFIED).** Against the exact law by a forward
DP for `n <= 9`: `sum_(prim) S_n = (2/3)(H_n(1) - H_(n-1)(1))` and the parity
sum `= (2/3)(H_n(-1) - H_(n-1)(-1))` agree to six digits at every level
(e.g. `n = 5`: `+0.024449` and `+2.472826`; `H_2(1) = 8/7`, `H_2(-1) = 22/7`).
Against S19's independent tree computation (`mazur_harmonic_mass_deep18`
output, `rho_n(±1)` for `n = 10..18`): the accumulated primitive sums give
`rho_n(1) = 0.42497, 0.39428, 0.35824, 0.33343, 0.31549, 0.30591, 0.29011`
for `n = 10..16` against S19's `0.42497, 0.39428, 0.35824, 0.33343, 0.31549,
0.30591, 0.29011`; the parity sums give `28.108, 42.169, 63.247, 94.861,
142.291, 213.434` for `n = 11..16` against S19's differences `28.11, 42.2,
63.2, 94.9, 142.3, 213.4` (and `0.325 (3/2)^n = 28.11, 42.17, 63.25, 94.88,
142.32, 213.47`). Two different algorithms (a tree-layer count; a Fourier
recursion on the unit group) agree to five digits: each certifies the other.

**What this changes (DIRECTION).** Theorem C's question `liminf H_n(1) > 0`
is the question whether the partial sums `H_n(1) = 1 + sum_(2<=m<=n) T_m`, `T_m = (3/2)
sum_(prim mod 3^m) S_m(psi)`, stay away from zero (the `m = 1` term `T_1 = -1/2` is already inside `1 = H_1(1)`); the `T_m` are `-1/2,
+0.143, -0.151, -0.064, +0.037, -0.009, -0.095, -0.086, -0.078, -0.059,
-0.046, -0.054, -0.037, -0.027, -0.014, -0.024` for `m = 1..16`, each a sum
of `(2/3) L_m` moments of size `0.84·3^(-m/2)` whose random-phase size would be
`0.56` — so the moments' phases cancel almost completely in the sum at `y =
1` (the sum *is* the mass at `1`); the identity is a bridge, not a bound.
The negative-cycle spike, by contrast, is the parity asymmetry of the
spectrum, and its growth `(3/2)^n` is the statement that odd and even
characters see the law differently by a margin `≈ 0.33 (3/2)^n / L_n ≈ 0.5
(1/2)^n` per character, below the Parseval scale `3^(-n/2)` by `(0.866)^n`.

---

## 3. The 3-adic zero lines behind the ridges (Chocian's census design)

**Objects.** A ridge of the negative family (S22, section 2d of the
five-mirrors note) is born at a coincidence `u·2^Q ≡ ∓1 mod 3^d` with `u`
small: the multiplier family `(∓u) 2^j` then sits at the negative exponents
`-Q + j` for every level `n <= d`. Define the *zero line* of the pair `(u,
Q)` (`u` odd, `3 ∤ u`, `Q >= 1`) as the sign for which `u 2^Q ∓ 1 ≡ 0 mod 3`
(exactly one sign works) and its *depth* `d(u, Q) = v_3(u 2^Q ∓ 1) >= 1`. Under
the uniform model (the residue `u 2^Q mod 3^d` uniform on its class) `P(d >=
k) = 3^(-(k-1))`; the number of lines of depth `>= k` in a box of `N` pairs is
approximately Poisson with mean `N 3^(-(k-1))`.

**Census (FINITE-EXACT; zero-lines outputs).** Box A: `u <= 1000` (`333`
multipliers), `Q <= 2000`, `N = 666,000`. Box B: `u <= 5000` (`1667`), `Q <=
10^4`, `N = 16,670,000`, depth capped at `20` (residues mod `3^20` in int64).
Box B's depth histogram against `N (2/3) 3^(-(d-1))`: ratios `1.000, 1.000,
1.000, 1.000, 1.000, 1.000, 0.999, 1.004, 0.998, 1.011, 0.951, 0.956, 1.100,
1.291, 1.291, 1.291, 3.87` for `d = 1..17` (the ratios at `d <= 8` are
forced: `L_d = 2·3^(d-1) <= W`, so each multiplier contributes `2W/L_d +
O(1)` lines of depth `>= d` deterministically by the equidistribution of the
powers of two — the audit measured a per-multiplier variance `0.13` against a
Poisson `41` at `d = 6`; the Poisson test begins at `d >= 10` here and at
`d >= 8` in box A); cumulative counts of depth `>= d`
against the Poisson mean, with the tail probability `P(X >= observed)`: `d =
11`: `276` vs `282.3` (`0.65`); `12`: `97` vs `94.1` (`0.40`); `13`: `37` vs
`31.4` (`0.18`); `14`: `14` vs `10.5` (`0.17`); `15`: `5` vs `3.5` (`0.27`);
`16`: `2` vs `1.16` (`0.32`); `17`: `1` vs `0.39` (`0.32`); none deeper. Box A
had shown `3` lines of depth `>= 14` against `0.42` expected (`P ≈ 0.009`)
and `4` of depth `>= 13` against `1.25` (`0.04`); the twenty-five-fold larger
box B contains those lines and dissolves the excess. The deepest lines:
`1187·2^5031 ≡ 1 mod 3^17`; `2441·2^8384 ≡ -1 mod 3^16`; depth 15 at `(3025,
846, -)`, `(1685, 8663, -)`, `(55, 423, +)` — the S22 seed `2^(-423) ≡ -55 mod
3^15` is one of five depth-15 lines in the box (expected `3.5`), not an
anomaly; depth 14 at `(4639, 4171)`, `(4573, 2226)`, `(4541, 7132)`, `(3853,
9384)`, `(3037, 5995)`, `(1997, 873)`, `(917, 1557)`, `(521, 1365)`, `(235,
8097)`. Also checked exactly: `v_3(13·2^154 - 1) = 7`, `v_3(2^486 - 1) = 6`,
`v_3(2^480 - 1) = 2` (the S22 seeds). Lines of depth `>= 5` with `Q <= 60`
(the window entering `Ñ_n`) exist at every `Q`, the deepest `(497, 41)` of
depth 11 and `(205, 5)`, `(133, 17)`, `(647, 27)`, `(901, 30)`, `(997, 35)`,
`(379, 45)`, `(491, 58)` of depth 8: at levels `n <= 8` the negative window
already carries several multiplier families, all of them with `u` in the
hundreds and hence of amplitude `0.05–0.35 M(n)` (the S22 trend), i.e.
`1–4` Parseval units — the size of the fluctuations seen in `Ñ_n 3^(n/2)`
at those levels.

**Reading.** The ridge seeds are Poisson-distributed 3-adic coincidences, as
Chocian's twisted-Bernoulli zeros are (his `chi^2 = 0.05, 3.93`; here every
tail probability is above `0.17`). The deepest line a box of `N` pairs
contains has depth `≈ log_3 N + O(1)` (the maximum of geometric variables), so
the largest ridge seeded by multipliers `u <= U` and exponents `Q <= W` has
depth `≈ log_3(UW)`, which is the quantitative form of S22's "polynomial in
the window" heuristic; a ridge's amplitude at its birth level `d` is `0.05–
0.35` of `M(d)`, `1.64^d` Parseval units up to the same prefactor, and
afterwards it decays as a remnant. What the census does not decide is the
mechanism by which one remnant (the level-5 ridge of S22) grew for forty
levels; that stays OPEN.

**The critical rate as a multiplier statement (PROVED identity; reading
DIRECTION).** The Chernoff mass law of S22 shifts by `log_2 u` for the family
`u·2^j`, predicting an amplitude `M(n) u^(θ*/ln 2) = M(n) u^(-0.438)`. At the
largest multipliers the family can have, `u ≈ 3^n`, this gives `e^(-nI)
3^(nθ*/ln 2) = (log_2 3 - 1)^n` **exactly**: `e^(-I) 3^(θ*/ln 2) = e^(-θ* m +
ln(m-1)) e^(θ* m) = m - 1` with `m = log_2 3` (the definitions of S22, `I =
θ* m - ln(m - 1)`). So the critical rate `log_2 3 - 1 = 0.585` of Theorem C
is the rate at which the resonant multiplier families would fill the unit
group if the Chernoff law held to `u ≈ 3^n`; Parseval (`sum_t |mu_hat|^2 =
O(1)` over `L_n` units) forbids a rate above `3^(-1/2) = 0.577` for the
generic unit, so the multiplier exponent must steepen below `-1/2` for large
`u` — which the S22 measurement (`-0.47..-0.55` at `u <= 127`) already shows.
Hypothesis H (rate of `Ñ_n` below `0.585`) is therefore the statement that
the negative family sits nearer the Parseval scale than the Chernoff
extrapolation: its `1.3%` margin is the gap between `3^(-1/2)` and `log_2 3
- 1`, two constants of the problem, not a fitted number.

---

## 4. 17, 1729, the saddle chain, and the basin of 1729

**4.1 The saddle chain (EXACT; saddle-chain output).** Proposition 4 of the
hedgehog note: for `x ≡ 1 mod 4^j`, `T^(2j)(x) - 1 = (3/4)^j (x - 1)`. Hence
the orbit of `1 + 4^m` under `T^2` runs down the chain `1 + 3^j 4^(m-j)`, `j =
0..m`, ending at `1 + 3^m` — `m` consecutive steps `x -> (3x+1)/4`, the
exactly linear part of the trivial cycle's saddle (real contraction `3/4`,
2-adic expansion `4`, 3-adic contraction `1/3`). For `m = 6` the chain is
`4097, 3073, 2305, 1729, 1297, 973, 730`, and

* `1729 = 1 + 12^3 = 1 + 3^3 4^3` is its balanced point `j = m/2`, the only
  chain point where the 3-adic and 4-adic exponents are equal; in general
  `1 + 12^m` is the midpoint of the chain of `1 + 4^(2m)` (`13, 145, 1729,
  20737, 248833, 2985985` for `m = 1..6`, checked by iteration);
* the chain bottom is `730 = 3^6 + 1`, so Ramanujan's second representation
  reads `1729 = 3^6 + 10^3 = (730 - 1) + 10^3`: the two taxicab
  representations are the midpoint (`12^3 + 1^3`) and the bottom (`3^6 +
  10^3`) of the same chain (EXACT as arithmetic; that the second is "of the
  chain" is a reading);
* `1729 = 7·13·19` — all three primes `≡ 1 mod 6`, the primes Chocian's first
  paper is about (`x^3 + 1 = (x+1)(x^2 - x + 1)` at `x = 12`: `13 · Φ_6(12)`,
  and the prime divisors of `Φ_6(x)` are `≡ 1 mod 6` or `3`); `1297 = 6^4 +
  1` is prime; `973 = 7·139` with `139 = 3^7 - 2^11` the clock of the
  sporadic `-17` cycle (THM-4484), by the accidental relation `4·3^5 + 1 =
  7(3^7 - 2^11)` — NUMEROLOGY;
* `17 = 1 + 4^2` heads the chain `17, 13, 10 = 3^2 + 1`; `-17` is the seed
  of the seven-odd-step cycle of `3x+1` on the negative integers (`-17, -25,
  -37, -55, -82, -41, -61, -91, -136, -68, -34`), whose clock is `139` and
  whose 3-adic resonance is S19 Theorem B's `(2187/2048)^(n/7)`; the
  exponents `4, 3` of `x^4 + y^3 = z^17` are the saddle's `4` and `3` —
  NUMEROLOGY, recorded because the owner asked how 17 relates; nothing
  structural connects the Fermat–Catalan signature to the cycle;
* the orbits of `1729` and `27` merge at `137` (`1729 -> 2594 -> 1297 -> 1946
  -> 973 -> 1460 -> 730 -> 365 -> 548 -> 274 -> 137`, ten `T`-steps; `27`
  reaches `137` in twelve), after which they share the well-known highway
  through `9232 = 2·4616` down to `1`; the coincidence with the fine-structure
  number is NUMEROLOGY.

The chain is the `j = 6` row of Theorem 2.3 of the extended-Collatz note
(`730, 973, 1297, 1729, 2305, 3073, 1024`: the greedy backward word of `3^6 +
1`, whose growth block lands on `4^5`), read forward; that `1729` sits on it
was not remarked there.

**4.2 The basins (FINITE-EXACT counts; densities OBSERVED).** `B(a) = {n >=
1 : the orbit of n passes through a}`. A memoised sieve (`_basins_...c`; the
code of each `n` is the first target hit, or that of the first orbit element
below `n`) on `[1, 2^30]`, checked against an independent Python sieve at
`2^21` (identical counts):

| `a` | `4097` | `3073` | `2305` | `1729` | `1297` | `973` | `730` | `27` | `137` |
|---|---|---|---|---|---|---|---|---|---|
| `dens B(a)` on `[1, 2^30]` | `0.003339` | `0.003524` | `0.004163` | `0.004395` | `0.004881` | `0.005795` | `0.014459` | `0` (`26` numbers) | `0.299245` |

Stable to four decimals from `2^24` on (`B(1729)`: `0.004388, 0.004384,
0.004386, 0.004389, 0.004391, 0.004393, 0.004395` at `2^24..2^30`). The chain
basins are nested (`B(4097) ⊂ ... ⊂ B(730)`); the increments are the side
entries `(2^k a - 1)/3`, `k >= 4`: tiny at `1729` (`0.00023`, through `9221,
36885, ...`), large at `730` (`0.0087`, through `3893 = 1 + 4·973`, the
saddle point above `973`). `B(27)` is the doubling ray of `27` (`27 ≡ 0 mod
3` has no odd preimage): the famous starting value has an empty tree above
it. `B(137)` holds three integers in ten: the path `137 -> 206 -> 103 -> 155 ->
233 -> 350 -> ... -> 577 -> ... -> 5` carries `30%` of all integers into the trunk-entry point `5`, whose basin is `0.938` (mac-mini's `e_2`).

**The family `1 + 12^m`, and the highway `11 -> 17 -> 13 -> 10 -> 5` (FINITE-EXACT counts to `2^30`; `_basins_mask_...c`, a bitmask sieve that assumes no orbit relation between its targets).** `dens B(1 + 12^m) = 0.476517, 0.035740, 0.004395, 0.000385, 0.000020` for `m = 1..5` (`13, 145, 1729, 20737, 248833`): about a decade per unit of `m`, i.e. per two extra links `(3x+1)/4` that the tree above must funnel through (ratios `13.3, 8.1, 11.4, 19`). The first two are not small numbers' basins but highways: `dens B(13) = 0.4765` and `dens B(17) = 0.4614` — nearly half of all integers pass through `17 -> 26 -> 13 -> 20 -> 10 -> 5` (and `11 -> 34 -> 17` above it), while the other odd preimages of `13` (`277, 1109, ...`) carry `1.5%` together; against `dens B(5) = 0.938` (mac-mini's `e_2`), the trunk entry `5` is fed `49%` through `13` and `49%` otherwise. `dens B(9232) = 0.000066`: the celebrated peak of the orbit of `27` has a basin of `6.6·10^(-5)` (it is `16·577`; its odd preimages `12309, 49237, ...` feed it). `dens B(137) = 0.299245` confirms the first sieve.


**4.3 The Krasikov–Lagarias bound at the root 1729.** Krasikov–Lagarias
(CITED; Theorem 6.1 of arXiv math/0205002, Acta Arith. 109 (2003), read by
the audit; the count is over the map `T`, the same as the sieve's): for every `a ≢ 0 mod 3` there is `X_0(a)` with `|B(a) ∩ [1, X]| >=
X^0.84` for `X >= X_0(a)`. Here `|B(1729) ∩ [1, 2^30]| = 4,719,191` while
`(2^30)^0.84 = 3.85·10^7`: the count is `12.2%` of the bound, rising by
`2^0.16 = 1.117` per doubling (`0.046` at `2^21`, `0.098` at `2^28`, `0.122`
at `2^30`). Consequences: (i) **`X_0(1729) > 2^30`** (PROVED by the count:
the inequality fails at `X = 2^30`); (ii) if the density `0.0044` persists,
the inequality first holds near `X = 2^30 · 1.117^(-log(0.122)/log(1.117)) =
2^(30 + 18.9) ≈ 2^48.9 ≈ 5·10^14` (equivalently `0.0044^(-1/0.16)`; OBSERVED
extrapolation). The theorem is asymptotic with an ineffective `X_0(a)`; for the root `1729`
the extrapolated threshold `≈ 2^49` lies inside the range checked
exhaustively before the theorem (`3·2^53`, Oliveira e Silva 1999, CITED via
the audit) and far below Barina's `2^68`, so for this root the sublinear
bound is weaker than the truth throughout the verified range (the sentence
first written here placed the threshold "astronomically beyond" that range:
MISTAKE-552).
mac-mini's THM-4517 shows why no root-uniform *positive proportion* can
exist; the numbers here show, for one root, how far the sublinear bound sits
below the truth. The owner's phrase "return to 1729" is, for `n > 1729`, the
event that the orbit of `n` passes through `1729` — `0.44%` of the integers.

---

## 5. The `(4,3,17)` statement

A web search on 2026-09-29 for `x^4 + y^3 = z^17` and the signature `(3,4,17)`
found no paper, survey entry or announcement (the survey arXiv:2412.11933v2
of solved signatures, checked by mac-mini, lists among `(3,4,n)` only `n = 4,
5`). The statement stays **UNVERIFIED**; if the owner has a source, it should
be cited by the next session. The same search found arXiv:2609.26996
(September 2026), *The primitive generalized Fermat equation `x^3 + y^5 =
z^7`: a computer-assisted proof* — the signature mac-mini's note names as
the smallest open Beal signature; recorded as AUTHOR-CLAIMED, not read, and
added as an addendum to that note's section 2.5. mac-mini's further sentence that the equation's clock shadow is empty ("no term can be a pure power of two while the other two are powers of one odd base") was refuted as unsupported by the audit below: a coprime solution has exactly one even variable, and the shadow is precisely the sub-case where it is a power of two (`x = 2^a` gives `2^{4a} = z^{17} - y^3`, a clock of shape `(4a, 17)`), which is as open as the equation itself. Nothing in this session bears on the truth of the statement.

---

## 6. Audit of mac-mini's THM-4515 / THM-4516 / THM-4517 (their first obligation)

An independent auditor (own code, `collatz_necklace_20260929_audit.py/.out/.md`)
was launched on the necklace note and its three theorem files with the brief
to re-derive before reading, recompute the censuses and densities, and check
statuses. **Verdict (report `collatz_necklace_20260929_audit.md`, 38 claims, own sieves to `2^30` and `2^32`, own perfect-power sieve over every exponent, own DFT checks on 1206 fair splits): SOUND WITH CORRECTIONS.** Holds: the discrete IVT, the circulant share equation, the DFT diagonalisation and the clock factorisation (exact in `Z[zeta_j]`), the CRT form, the existence criteria, the `K <= 18` census to the last row, `E[N_j]` and the Stirling form, the fair Eliahou identity, the 30 necklaces of shape `(11,7)`, the dictionary of THM-4516, all four perfect-power censuses (now covered for every exponent, where the session's `iroot` covered `r <= 64`), every free cycle by direct iteration, the Eisenstein square, the disjoint trunk-entry basins, the superadditivity argument, the sheet cap, every FINITE-EXACT density table to the last digit, the Krasikov–Lagarias citation, HYP-9165's OPEN status. Refuted: (1) THM-4516(5)'s consequence "Beal implies no map `py ± m^r` (`r >= 3`) has a free mixed shape with `K, sX >= 3`, `X >= 2`" — freeness needs only `(2^K - p^X) | m^r` (THM-4484), and `3y + 125` has the free shape `(5,3)` with clock `2^5 - 3^3 = 5 | 125` (two integer cycles, no Beal solution); the correct consequence is that no perfect-power *clock* has all three exponents `>= 3`; (2) "the clock shadow of `x^4 + y^3 = z^17` is empty" — unsupported (section 5). Corrections of statement: `B(5)` and `B(32)` partition `B(1) \ {1,2,4,8,16}`, not the integers minus the powers of two (that would be Collatz); `B(2^{2i-1})` includes the trunk above; the trunk-entry densities summing to `1` implies almost every orbit reaches `1` but the converse is unproved; `e_3, e_6, e_9, e_12` vanish to seven decimals (doubling rays), not exactly; the third basin is flat to three decimals, not four; the sheet cap is the limit of the still-rising `{5,7,10}` basin, `>= 0.3250`, not `0.3248`; the `1/log m` lower bound is PROVED only at the computed `m` (its order EMPIRICAL); minor wording (the eigenvalues up to units; the least-absolute-value elements of the `(11,7)` necklaces; `21 = 7·3`; the ordering sense of "smallest open Beal signature"). All eighteen textual corrections were applied to mac-mini's note, theorem files and HYP-9165 with attribution tags; the two refuted consequences are MISTAKE-551. Nothing changes the Collatz status.

---

## 7. Typing, what changes, obligations

| item | status |
|---|---|
| Theorem 1 (character-spectrum duality; `sum_prim S_n = rho_n(1) - rho_(n-1)(1)`; parity sum `= rho_n(-1) - rho_(n-1)(-1)`) | PROVED; VERIFIED against the exact law (`n <= 9`, six digits) and against S19's independent values (`n = 10..16`, five digits) |
| Proposition 2 (full-period one-pole recursion, exact) | PROVED; `M(n)` agrees with the S20 FFT at `n = 1..18` to seven digits; level 19 reproduced independently by the audit |
| `M(19) = 0.0098157` at `±2^24`, `24 = floor(19 log_2 3) - 6`; the maximum over all `774,840,978` units on `±2^s` | VERIFIED (float64, no truncation; `_fullperiod_max19_` output) |
| `M(20) = 0.0088846` at `±2^26` (`20 log_2 3 - 5.7`), over all `2,324,522,934` units | VERIFIED (complex64 at level 20, no truncation; `_fullperiod_max20_` output) |
| the multiplicative spectrum: rms at the Parseval scale, non-Rayleigh, largest moments `≈ n/2` Parseval units, on `psi_(±2), psi_(±8)` at low levels | VERIFIED (`n <= 16`); the reading of `psi_(±2)` as the 3-adic-logarithm character DIRECTION |
| the Jacobi transfer `E[psi(Y_n)] = G_psi sum_(psi') c(psi, psi') E[psi'(Y_(n-1))]`; `|c| = (#prim')^(-1/2)` on primitive pairs, `c = 0` off them; `rms |G_psi| = 3^(-1/2)`; `corr(|S_n|, |G_psi|) = 0.6–0.8` | recursion PROVED; the constant modulus VERIFIED to `n = 7` (exact to `10^(-6)`), CONJECTURED for all `n`; the shape reading OBSERVED |
| the 3-adic zero-line census: geometric depth law (deterministic for `d <= 8`, Poisson-testable from `d >= 10`), Poisson tail, the S22 seed one of five depth-15 lines | FINITE-EXACT (`16.7·10^6` pairs) |
| `e^(-I) 3^(θ*/ln 2) = log_2 3 - 1` | PROVED (a one-line identity); the reading of H's margin DIRECTION |
| the saddle chain `1 + 3^j 4^(6-j)`; `1729 = 1 + 12^3` its midpoint; `1 + 12^m` the midpoint for all `m`; `17 = 1 + 4^2` | EXACT (PROVED from the hedgehog note's Proposition 4; checked) |
| `139 | 973`, the merge at `137`, the exponents `4, 3` | NUMEROLOGY |
| basins of the chain, of `27`, of `137` to `2^30`; `dens B(1 + 12^m) = 0.48, 0.036, 0.0044, 3.9·10^(-4), 2.0·10^(-5)` (`m = 1..5`); `dens B(13) = 0.477`, `dens B(17) = 0.461`, `dens B(9232) = 6.6·10^(-5)` | FINITE-EXACT counts; densities OBSERVED (two codes agree at `2^21`; the bitmask sieve agrees with the chain sieve at `2^30`) |
| `X_0(1729) > 2^30` for the Krasikov–Lagarias threshold | PROVED (by the count) |
| `X_0(1729) ≈ 2^49`, inside the range verified before the theorem (`3·2^53`, CITED) | OBSERVED extrapolation |
| `x^4 + y^3 = z^17` has no primitive solution | UNVERIFIED (no source found); its "clock shadow" is not empty by any known argument (audit of THM-4516) |
| arXiv:2609.26996 on `(3,5,7)` | AUTHOR-CLAIMED (found, not read) |

**What changes for the repo.** (1) The seed-1 mass and the negative-cycle
spike have a second exact expression each, as accumulated primitive
character sums; the two S19 computations are now cross-validated by a
different algorithm. (2) The primitive Fourier maximum over all units is
known two levels further (19 and 20) by a method whose cost is `O(3^n)` in memory only, with no truncation. (3)
The ridge seeds are Poisson; the S22 obligation "the ridge inventory as a
function of `v_3(u 2^Q ∓ 1)`" is done, and H's margin is identified as the
gap between two constants. (4) `1729` enters the Collatz thread with an exact
place (the balanced point of a saddle chain) and a measured basin, and the
Krasikov–Lagarias threshold for it is bounded below. (5) The `(4,3,17)`
statement is still without a source.

**Obligations.** (a) The Jacobi transfer is now written and its constant-modulus structure verified to level 7 (section 2); to prove `|c(psi, psi')| = (#prim')^(-1/2)` for all `n` (the prime-power Jacobi-sum evaluation) and to derive the spectrum's tail law from `|G_psi|` and the mixing. (b) Level 20 done in complex64 (`M(20) = 0.0088846` at `±2^26`); the audit remarks that streaming the last level needs only `m_19` in memory (`12.4 GB`) with residues by doubling, so level 21 (`m_20`, `37 GB`) is within the machine and level 22 is not. (c) The growth mechanism of the level-5
remnant (S22), untouched here. (d) `dens B(1 + 12^m)` falls by about a decade per unit of `m` (measured to `m = 5`); a law, and the basins of the other chain points `1 + 3^j 4^(m-j)`, remain to be found.
(e) A source for the `(4,3,17)` statement.

---

## 8. Reproduction and audit record

```
python 04-computation/experiments/collatz_three_mirrors_character_spectrum_20260929.py 16 > 05-knowledge/results/collatz_three_mirrors_character_spectrum_20260929.out   # 1 min, 5 GB
python 04-computation/experiments/collatz_three_mirrors_fullperiod_max_20260929.py 18 2   > 05-knowledge/results/collatz_three_mirrors_fullperiod_max_20260929.out   # 2 min, 6 GB
python 04-computation/experiments/collatz_three_mirrors_fullperiod_max_20260929.py 19 19  > 05-knowledge/results/collatz_three_mirrors_fullperiod_max19_20260929.out   # 10 min, 17 GB
python 04-computation/experiments/collatz_three_mirrors_fullperiod_max_20260929.py 20 19 20 > 05-knowledge/results/collatz_three_mirrors_fullperiod_max20_20260929.out   # 30 min, 37 GB (complex64 at level 20)
python 04-computation/experiments/collatz_three_mirrors_jacobi_transfer_20260929.py 7      > 05-knowledge/results/collatz_three_mirrors_jacobi_transfer_20260929.out
python 04-computation/experiments/collatz_three_mirrors_zero_lines_20260929.py 1000 2000  > 05-knowledge/results/collatz_three_mirrors_zero_lines_20260929.out
python 04-computation/experiments/collatz_three_mirrors_zero_lines_wide_20260929.py 5000 10000 > 05-knowledge/results/collatz_three_mirrors_zero_lines_wide_20260929.out
python 04-computation/experiments/collatz_three_mirrors_saddle_chain_20260929.py         > 05-knowledge/results/collatz_three_mirrors_saddle_chain_20260929.out
gcc -O3 -o basins3 04-computation/experiments/collatz_three_mirrors_basins_20260929.c && ./basins3 1073741824 > 05-knowledge/results/collatz_three_mirrors_basins_20260929.out   # 15 s, 1 GB
gcc -O3 -o basins_mask 04-computation/experiments/collatz_three_mirrors_basins_mask_20260929.c && ./basins_mask 1073741824 13 145 1729 20737 248833 17 137 9232 > 05-knowledge/results/collatz_three_mirrors_basins_mask_20260929.out   # 15 s
```

**Own audit (2026-09-29, a second subagent, own code: forward DP of the law to `n = 9`, the full-period recursion re-implemented in C to level 19, an own bitmask sieve to `2^30`, own censuses in Python and C, the preprint texts checked; report `collatz_three_mirrors_20260929_audit.md`, 40 claims): SOUND WITH CORRECTIONS, eighteen applied.** Holds: Theorem 1(i), (ii) with own proofs, the corollary `-1/3`; Proposition 2 and every level `1..18` to seven digits, level 19 reproduced (`0.0098157`, mass `0.7089`); the spectrum statistics, the outlier characters, the parity split; the census, the deepest lines, the four exact valuations; the rate identity; the saddle chain and every factorisation; the basins to the last digit (own bitmask sieve), the nestedness, the empty tree above `27`, `dens B(137) = 0.2992`; the Krasikov–Lagarias quotation and the logic of `X_0(1729) > 2^30`; the attributions to the preprints. Corrections: (1) "five digits" holds for `rho_n(1)` at `2..16`; the parity differences are known to the three or four digits S19 prints; `rho_12(1) = 0.35824`; (2) the Fourier-mass explanation (the factor `3/2` is `3^n/L_n`, the same sum over the units, not extra characters); (3) the summary's `s = floor(h log_2 3) - 6` holds at `h = 19` only — the rule is `h log_2 3 - 6 ± 1`; (4) the top four conjugate pairs, not "top eight"; (5) a factor `2` in the `psi_(±2)` formula; (6) `H_n(1) = 1 + sum_(2<=m<=n) T_m` (the note's `1 + sum_(m<=n)` double-counted `T_1`); (7) the random-phase size `0.56` for `(2/3) L_m` primitive characters; (8) the histogram is geometric through depth 10, and for `d <= 8` its ratios are forced by the equidistribution of the powers of two (each multiplier contributes `2W/L_d + O(1)` lines deterministically; per-multiplier variance `0.13` against Poisson `41` at `d = 6`) — the Poisson test starts at `d >= 10`; (9) the critical-rate reading is typed CONJECTURAL/DIRECTION in the summary; (10) `137 -> 206 -> 103`; (11) `2^48.9 ≈ 5·10^14`, not `2^49.3 ≈ 7·10^14`; (12) the sentence placing the threshold "astronomically beyond any range checked exhaustively" was false — `2^49` lies inside the `3·2^53` of Oliveira e Silva 1999 and far below Barina's `2^68` (MISTAKE-552); (13) at audit time section 6 was a placeholder cited as a result ("confirmed by the audit") — the placeholder had been filled and the sentence corrected before the report arrived, and the record of that slip is MISTAKE-552; (14) "(odd `m`)" in the Gauss-sum identity of the second preprint; (15) the level-20 memory remark. The audit did not cover section 6, the history remarks, the truth of the `(4,3,17)` statement, or arXiv:2609.26996.
