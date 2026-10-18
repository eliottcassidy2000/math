---
id: THM-4690
title: "First arrivals of the signed ROOT transport into a stopping set: the backward ladder, port-diagonal energy, integer-clock Gram and clock spectrum (E_V is the clock spectrum at the Nyquist frequency), and the exact source pair ledger E0 = 34N/7, sibling pairs 5N/7, so E_V = 24N/7 + (non-sibling pair sums) + O(|V| + log N); plus the 19-adic obstruction: no path of unit steps integral at the eight h = 7 primes above 19 joins two translated coefficient copies (or two sheets) of clock_duality's cubic golden recipe"
status: >
  PROVED. Statements 1-3 are elementary bookkeeping. Statement 4 is an elementary residue computation along B-chains. Statement 5 is a
  valuation argument with exact field checks.
  FINITE-EXACT companions:
  - the identities of 1-4 hold at every N < 2^27 (fap.c, fap_levels.c), and the pair ledger of 4 at every N < 2^28 (fap_lean.c),
    with the invariant R(N) = H(4N) - H(N) checked at every N;
  - the constants of 4 are derived symbolically and checked by brute force (source_ledger_constants.py);
  - the ingredients of 5 are exact arithmetic in the degree-48 field and its residue fields F_361
    (nineteen_obstruction_check.py, nineteen_places.py, diagonal_19adic_check.py).
  The inequalities suggested by the data are HYP-9325 (OPEN).
  No priority is claimed for the bookkeeping. Statement 1 is clock_duality_oct15's R6 setting enlarged by backward closure.
  The constant 34/7 of statement 4 is KNOWN in the repo: it is the d = 0 value E_0 = 4.857143 of THM-4667's ideal energy
  ladder (statement 6), and T_0/X^2 -> 51/7 = (3/2)(34/7) is THM-4667(7)'s observed octave-trace limit. What is new in 4
  is the proof for the actual sources with error O(log N), the exact sibling constant 5/7, and the pair ledger.
  Statement 3 carries out clock_duality's instruction to keep the literal integer clocks in the first-arrival Gram
  (root_sparse_boundary_module_20261017i.md, end of section 8).
  Renumbered from THM-4670, which is ceded to opus-2026-10-18-S43.
source: mac-mini-2026-10-18-sixcolour
depends_on: [THM-4667 (signed ROOT transport, energy ledger), HYP-9290 (frontier energy cone; Lemma F; basin form), clock_duality_oct15 R5/R6 (root_sparse_boundary_module_20261017i.md), golden_contact_circuit_20261017i.md section 4 (cubic recipe)]
verification:
  - 04-computation/experiments/first_arrival_20261018/fap.c, fap_lean.c (+ fap_*.out): port, sibling-pair, depth, integer-clock, spectrum and pair-ledger quantities at every N, with the invariant checked
  - 04-computation/experiments/first_arrival_20261018/fap_levels.c (+ fap_levels_*.out): the backward ladder E^(k)(N)
  - 04-computation/experiments/first_arrival_20261018/source_ledger_constants.py (+ .out): derivation of 34/7 and 5/7 and brute-force check
  - 04-computation/experiments/first_arrival_20261018/nineteen_obstruction_check.py, nineteen_places.py, diagonal_19adic_check.py, copy_connectors.py, copy_connectors_screen.py (+ .out, .json)
note: 05-knowledge/results/first_arrival_one_tick_coefficient_copies_gibonacci_20261018.md
related: [HYP-9325, HYP-9290, HYP-9289, THM-4667, THM-4668]
---

# THM-4690 — First arrivals: ladder, ports, integer clocks, the source pair ledger; the 19-adic copy obstruction

## Setting (ROOT part)

This uses the notation of THM-4667 and HYP-9290.

- **Seeds and transport.** Seeds are odd `s` with `s ≢ 5 (mod 8)`. The transport is `P(s) = seed((3s+1)/2^a)`, where
  `a = ν_2(3s+1) ∈ {1, 2}` and the sign is `σ = (-1)^a`. Here `seed` strips `B(x) = 4x + 1`.
- **Literal clock.** One transport step has literal clock `a + 2·(number of strips)`. This is the difference of the
  shortcut stopping times, since a strip changes the stopping time by 2 and the step by `a`.
- **Sources.** Every odd `m <= 4N` carries `K_N(m) = w(m) + thr(m)[m < 2N]` and sits on its seed. Here `w(m) = 2` if
  `m ≡ 2 (mod 3)` and `-1` otherwise, and `thr = 1, -2, 1` for `m ≡ 0, 1, 2 (mod 3)`. The seed coefficient is `c_0(x)`.
- **Stopping set.** `V` is a finite set of seeds. Each source is stopped at its first seed in `V`; every source considered
  reaches `V`. When `V` is forward closed with a rank descending to ROOT (clock_duality's R6), it gives the certified
  frontiers `V_100` (74 seeds) and `V_1000` (627 seeds, largest 63665).
- **Entry data.** For a source `x`:
  - `u(x)` is its entry seed;
  - `t(x)` is the literal clock of its path to `u(x)`;
  - `dep(x)` is the number of `P`-steps to entry;
  - its **entry port** is the last seed of the path outside `V` (none if `x ∈ V`).
- **Energies.** `c_V(u) = Σ_(u(x) = u) σ(x → u) c_0(x)` and `E_V = Σ_u c_V(u)^2`.

## Statements

1. **Backward ladder.** Let `V` be forward closed, and for `k >= 0` let `V^(k) = {y : dep(y) <= k}`.
   - `V^(k)` is forward closed.
   - Stopping on first entry into `V^(k)` freezes each source at the point of its path with `dep = k`, or at itself if
     `dep(x) < k`. Hence

     `E^(k)(N) = Σ_(dep(y) = k) J(y)^2 + Σ_(dep(x) < k) c_0(x)^2`,

     where `J(y) = Σ_(x in tree(y)) σ(x → y) c_0(x)` is the subtree current. For `y` in the ROOT basin,
     `J(y) = f(y) R_y(N)`, with `R_y` the ROOT current restricted to the sources through `y` (HYP-9290 basin form).
   - `E^(0) = E_V`, and `E^(K) = E_0 := Σ_x c_0(x)^2` for `K >= max dep`.
   - Every step is a merge identity:

     `E^(k) - E^(k+1) = 2 Σ_(dep(y) = k) [ Σ_(z < z' in P^(-1)(y)) σ(z)σ(z') J(z)J(z') + c_0(y) Σ_(z in P^(-1)(y)) σ(z) J(z) ]`.

2. **Ports.** `E^(1)` is the **port-diagonal energy**

   `E^(1) = Σ_(u in V) c_0(u)^2 + Σ_p J(p)^2`,

   where `p` runs over the ports `P^(-1)(V) ∖ V`. The first-arrival cross term `X(N) := E_V - E^(1)` is the level-0 merge
   term of statement 1.

   The preimages of a seed `y` are described by the residues `B^j(y) ≡ y + j (mod 3)`:
   - if `B^j(y) ≡ 1`, there is an `a = 2` preimage `s_j = (4B^j(y) - 1)/3` with step clock `2j + 2`;
   - if `B^j(y) ≡ 2`, there is an `a = 1` preimage `(2B^j(y) - 1)/3` with step clock `2j + 1`;
   - if `B^j(y) ≡ 0`, there is none.

   So the preimages are sibling pairs `(s_j, 2s_j + 1)` with step clocks `(2j + 2, 2j + 3)` and signs `(+, -)`. There is
   one unpaired preimage `(2y - 1)/3`, with step clock 1, when `y ≡ 2 (mod 3)`.

3. **Integer clocks and the clock spectrum.**
   - Along every path `σ(x → u) = (-1)^(t(x))`.
   - Define `h_t(u) = Σ_(u(x) = u, t(x) = t) σ(x → u) c_0(x)` and `G(D) = Σ_u Σ_t h_t(u) h_(t+D)(u)`. Then

     `E_V = G(0) + 2 Σ_(D >= 1) G(D)`.

   - Equivalently, take the unsigned clock polynomial `Q_u(z) = Σ_(u(x) = u) c_0(x) z^(t(x))` and the clock spectrum
     `S(θ) = Σ_u |Q_u(e^(iθ))|^2`. Then
     - `E_V = S(π)`: the signed energy is the clock spectrum at the Nyquist frequency;
     - `E^+ := Σ_u (Σ_(u(x)=u) c_0(x))^2 = S(0)`: the unsigned arrival energy;
     - `G(0) = (1/2π) ∫ S(θ) dθ`;
     - `S(θ) = G(0) + 2 Σ_(D>=1) (-1)^D G(D) cos(Dθ)`.

     So first-arrival cancellation, `E_V` small against `G(0)` and `E^+`, is a dip of the clock spectrum at `θ = π`.
     A negative one-tick term `G(1)` produces exactly such a dip.

4. **Source pair ledger.** Let `SIB(N) = Σ c_0(x) c_0(2x+1)`, summed over `x ≡ 1 (mod 8)` with `2x + 1 <= 4N`. Then:
   - `|c_0| <= 3`;
   - `E_0(N) = (34/7) N + O(log N)` and `SIB(N) = (5/7) N + O(log N)`.

   A sibling source pair `x`, `2x + 1` has a common image: `P(x) = P(2x+1)`, because `3(2x+1) + 1 = 2·B((3x+1)/4)`. The
   two step clocks are `2 + 2r` and `3 + 2r`, where `r` is the number of strips of `(3x+1)/4`. So the pair enters every
   stopping set at a common seed with clock lag exactly 1 and opposite signs, unless `x` or `2x + 1` lies in `V`. Hence,
   for every finite stopping set `V`:

   `E_V = (24/7) N + (G(0) - E_0) + 2 Γ_1 + 2 Σ_(D >= 2) G(D) + O(|V| + log N)`,

   where `Γ_1 := G(1) + SIB_V` collects the lag-1 pairs other than sibling sources. `SIB_V` is `SIB` restricted to pairs
   with a common entry, and `|SIB - SIB_V| <= 18|V|`.
   - `G(0) - E_0 = 2 Σ` (same-clock pairs of distinct sources with a common entry);
   - `G(D) = Σ` (pairs at clock lag `D`).

   In words: the source energy is `34N/7`; the sibling pairs remove exactly `10N/7` of it in every stopping set; the rest
   of the cancellation is carried by same-clock pairs, by other lag-1 pairs and by longer lags. For `V = {1}`, the ROOT
   itself, `E_V = R(N)^2` whenever every source `<= 4N` reaches 1 (true throughout the verified Collatz range). Then
   `R(N)^2 - 24N/7` is the sum of the non-sibling pair correlations of the full transport.

5. **The 19-adic obstruction to connecting coefficient copies.**
   - **Setting.** `F = Q(ζ15)(√-11)(w)`, with `w^3 + 2w^2 - 4w - 4 = 0`, is the degree-48 field of the golden contact
     construction. `b, c ∈ F` are the coefficients of clock_duality's `u`-square circuit `u^2 = conj(d) + b u + c Y`.
     `C` is the 10-point coefficient set and `M` the 7-point Moser spindle.
   - **Statement.** No path of unit steps whose directions all lie in a 19-integral set joins a point of `C + b m` to a
     point of `C + b m'` for `m ≠ m'`; the same holds for `c`. Here 19-integral means integral at every prime above 19 of
     `F` (or of any finite extension of `F` containing the steps).
   - **Example.** This applies to the 90 directions `{ζ^l q}`, where `ζ^l` runs over the 30 roots of unity of `F` and `q`
     over the unit differences of `C` together with 1.
   - **Prime by prime.** The 19-maximal tower order has exactly 24 primes above 19, all with residue field `F_361`. They
     fall into three types of eight: `h ≡ j`, `h ≡ -j` and `h ≡ 7 (mod 𝔭)`, where `h = w + 3`.
     - `b` and `c` have simple poles at all 24 primes, and `b/c` is a 19-unit.
     - The spindle differences are units at all 24.
     - At the eight `h ≡ 7` primes the following are all units: the contact unit `u = (h - j)/(h + j)`, the mixed-golden
       unit `u_2 = (h + φ - j)/(h + φ + j)` and the 30 roots of unity.
     - `u` has valuation -1 / +1 at the primes of type `h ≡ -j` / `h ≡ j`.
   - **Sheets of the whole recipe sum.** Take the sum `(A + aM) + (C + bM) + (C + cM)` with `a = conj(d)`; `A`, `C` and
     `a` are 19-integral. For every nonzero `(Δm_2, Δm_3) ∈ (M - M)^2`, the element `b Δm_2 + c Δm_3` has a pole at some
     `h ≡ 7` prime: of the 1224 pairs, 1166 have poles at all eight, 16 at seven and 42 at six. Hence no path whose
     steps are integral at all eight `h ≡ 7` primes joins two different sheets `(m_2, m_3)`. In particular this holds at
     every length for paths with steps in `{roots of unity}·u^Z·u_2^Z·{19-integral elements}`. It also holds for two
     copies inside one coefficient factor.
   - **Consequence.** A coefficient context connecting the seven copies must use unit directions with a pole above 19.
     To change sheets it needs a pole at an `h ≡ 7` prime, which `u` does not have.
     The contact unit `u = (h - j)/(h + j)`, with `h = w + 3` and `j = √-11`, is one. Its relative norms are
     `N_(Q(w,j)/Q(j))(h ± j) = -p(∓j) = -76 = -4·19`, where `p(x) = x^3 - 7x^2 + 11x - 1` is the minimal polynomial of
     `h`. So `u` has valuation `-1` at the prime above 19 where `h ≡ -j`.

## Proofs

**1.** If `dep(y) = k + 1 > 0`, then `y ∉ V` and `dep(P(y)) = k`, so `P(V^(k)) ⊆ V^(k)`.

A source with `dep(x) >= k` meets `dep = k` exactly once on its path, because `dep` drops by one per step. No earlier
point of the path lies in `V^(k)`. A seed with `dep < k` receives no stopped source other than itself, since any path
through it meets `dep = k` first.

For the merge identity:
- `J(y) = c_0(y) + Σ_(z in P^(-1)(y)) σ(z) J(z)`.
- As `y` runs over `{dep = k}`, the preimage sets `P^(-1)(y)` partition `{dep = k+1}`.
- Square, and subtract `E^(k+1)`.

Basin form: `σ(x → y) = f(x) f(y)`, because `f(P(s)) = (-1)^a f(s)` and strips preserve `f` (HYP-9290). ∎

**2.** `E^(1)` is the case `k = 1` of 1.

For the preimages: `(3s + 1)/2^a = B^j(y)` with `s` integral forces `B^j(y) ≡ 2^(-a) (mod 3)`. Both kinds are
automatically seeds: `a ∈ {1, 2}` exactly, because `B^j(y)` is odd, and `a ∈ {1, 2}` excludes `s ≡ 5 (mod 8)`. Its sibling is
`(2B^(j+1)(y) - 1)/3 = (8B^j(y) + 1)/3 = 2 s_j + 1`. The step clocks are `2 + 2j` and `1 + 2(j+1)`. ∎

**3.** Along a path, `σ = Π (-1)^(a_i) = Π (-1)^(a_i + 2·strips_i) = (-1)^t`.
- `c_V(u) = Σ_t h_t(u)`, and expanding the square by lags gives the first identity.
- `h_t = (-1)^t [z^t] Q_u`, so `c_V(u) = Q_u(-1)`.
- Parseval gives `G(0) = ∫ S/2π`.
- The cosine series is the Fourier expansion of `|Q_u|^2`. ∎

**4.** Write the B-chain of a seed `x` as `x, B(x), B^2(x), ...`; then `B^j(x) ≡ x + j (mod 3)`.

Only the first chain element can be a seed, since `B(x) ≡ 5 (mod 8)`. So every odd `m <= 4N` lies on exactly one chain.

The weights are:
- below `2N`: `K = κ(m mod 3)`, with `κ = (0, -3, 3)`;
- in `[2N, 4N]`: `K = λ(m mod 3)`, with `λ = (-1, -1, 2)`.

Since `B` multiplies by about 4, at most one chain element lies in `[2N, 4N]`. Let `L` be the number of chain elements
below `2N`. Then

`c_0(x) = P(x mod 3, L mod 3) + [chain ends in [2N, 4N]]·λ((x + L) mod 3)`,

where `P(r, ℓ) = Σ_(j<ℓ) κ(r + j)` and `κ` sums to zero over a period.

Writing `y = 3x + 1` and `Y = 6N + 1`:
- level `L` means `4^(L-1) y < Y <= 4^L y`;
- the chain ends in `[2N, 4N]` iff `4^L y <= 2Y - 1`.

So, up to `O(1)` at the endpoints, the `x`-ranges have length `2N/4^L` (chain ends in `[2N, 4N]`) and `4N/4^L` (it does
not, `L >= 1`).

Seeds have density `3/8`, with `x mod 3` equidistributed (CRT). Tabulating:
- If the chain ends in `[2N, 4N]`, the possible sums are `-1, -1, 2` for each `L mod 3`, so the average of `c_0^2` is 2.
- Otherwise the average of `c_0^2` is `0, 6, 6` for `L ≡ 0, 1, 2`.

Hence `|c_0| <= 3`, and

`E_0/N -> (3/8) [ 2·2·Σ_(L>=0) 4^(-L) + 4·6·Σ_(L>=1, 3∤L) 4^(-L) ] = (3/8) [16/3 + 24·(20/63)] = 34/7`.

For `SIB`, use `x ≡ 1 (mod 8)`, which has density `1/8`, and `2x + 1 ≡ 2r + 1 (mod 3)`. The sibling has
`3(2x + 1) + 1 = 2y + 2`. So:
- (A) if `x`'s chain ends in `[2N, 4N]` at level `L`, the sibling's chain does not end there, at the same level;
- (B) otherwise the sibling's chain ends there, at level `L - 1`.

Both products average 3 when `L ≢ 0 (mod 3)`, and 0 when `L ≡ 0`. Hence

`SIB/N -> (1/8)·3·(2 + 4)·Σ_(L>=1, 3∤L) 4^(-L) = (18/8)(20/63) = 5/7`.

The errors are `O(1)` per level, and there are `O(log N)` levels.

The common image: `(3x + 1)/4 = z` gives `3(2x + 1) + 1 = 2(4z + 1) = 2B(z)`, and `seed(B(z)) = seed(z)`. The step clocks
are `2 + 2·strips(z)` and `1 + 2(1 + strips(z))`.

If neither sibling lies in `V`, they share the rest of the path. Their contribution to `G(1)` is
`σσ' c_0 c_0' = -c_0(x) c_0(2x+1)`. Each seed of `V` belongs to at most two sibling pairs, so pairs split by `V` number at
most `2|V|`, each with `|c_0 c_0'| <= 9`.

The displayed identity is `E_V = E_0 - 2 SIB_V + (G(0) - E_0) + 2Γ_1 + 2 Σ_(D>=2) G(D)`, rewritten with the asymptotics.

The table computation is carried out symbolically in `source_ledger_constants.py`. The brute force gives
`|E_0 - 34N/7| <= 11.2` and `|SIB - 5N/7| <= 4.9` for the listed `N <= 4,000,037`. ∎

**5.** Four facts, then the conclusion.
- **Maximality at 19.** The tower order `Z[ζ15][j][w]` has discriminant supported on `{2, 3, 5, 11, 37}`; the cubic has
  discriminant `592 = 2^4·37`. So the order is maximal at 19, and an element with 19 in a coordinate denominator is
  non-integral at some `𝔭 | 19`.
- **The coefficients.** By exact computation, `b` and `c` have coordinate denominator exactly 19, and
  `v_19(N_(F/Q)(b)) = v_19(N_(F/Q)(c)) = -48`. So `v_𝔭(b) < 0` for some `𝔭 | 19`.
- **The spindle.** For all 21 pairs `m ≠ m'` of the spindle, `N_(F/Q)(m - m')` is prime to 19. So `m - m'` is a
  `𝔭`-unit, and `v_𝔭(b(m - m')) = v_𝔭(b) < 0`.
- **C is integral.** The points of `C` are 19-integral, so `v_𝔭((c' + b m') - (c + b m)) < 0`.

A path with steps in a set that is integral at every prime above `𝔭` has an integral displacement. So no such path
exists.

The 90 named directions are 19-integral by exact check (`nineteen_obstruction_check.out`).

For the contact unit:
- With `j^2 = -11` one gets `p(j) = p(-j) = 76`, so `Π_i (h_i ± j) = -p(∓j) = -76`.
- 19 is inert in `Q(j)`, since `-11` is a non-residue mod 19.
- Mod 19, `p(x) ≡ (x^2 + 11)(x - 7)`. So over `Q(j)` the prime 19 splits into three primes of `Q(w, j)` of relative
  degree 1, with `h ≡ j`, `h ≡ -j` and `h ≡ 7` respectively.
- The relative norm `-76` then gives `v(h + j) = 1` exactly at the prime with `h ≡ -j`, and `v(h - j) = 0` there.
- 19 is unramified in `Q(ζ15)`, so these valuations persist in `F`.

Prime by prime (`nineteen_places.py`):
- The 24 primes are the Frobenius orbits of the 48 homomorphisms `Z[ζ15, j, w] -> F_361` (8 × 2 × 3 choices).
- `v_𝔭(x) = -k` when 19^k exactly divides the coordinate denominators and `19^k x` does not vanish at `𝔭`.
- With this, the listed valuations and the counts 1166 / 16 / 42 are exact finite computations.
- The sheet statement follows as in the single-factor case, applied at the `h ≡ 7` prime where `b Δm_2 + c Δm_3` has its
  pole.
- `diagonal_19adic_check.py` confirms independently that no nonzero `b Δm_2 + c Δm_3` is 19-integral. ∎

## Remarks

- **Connector census, inside `F`, agrees with 5.** `copy_connectors.py` and `copy_connectors_screen.py` found:
  - no two-step connector between copies (172 + 118 point pairs at distance `<= 2`, all certified non-squares);
  - no three-step connector whose first step lies in the 90 integral directions (5 survivors of the 12-place screen,
    each certified non-square at a degree-1 place, `p = 4051` or `4261`);
  - none whose first step lies in the 120 19-adic directions `ζ^l u^k`, `k = ±1, ±2` (10 survivors, each certified
    non-square at `p = 3631`, `4051` or `4261`);
  - none whose first step lies in the 240 K-directions `ζ^l δ_𝔮^(±1)`. Here `δ_𝔮 = conj(γ_𝔮)/γ_𝔮`,
    `γ_𝔮 = j - c_0 - c_1 ζ` vanishes at one prime `𝔮` of `K` above 19 from each of the four conjugate pairs, and
    `δ_𝔮` has a pole at the three primes of `F` above `𝔮` (37 survivors, each certified non-square at an 80-place
    screen; places where an element is not integral are skipped).

  A single quadratic non-residue at a degree-1 place where the element is a unit proves that it is not a square in `F`.
  Connectors with intermediate points outside `F` are not covered.
- **Where connectors must come from.** Unit directions with a pole at an `h ≡ 7` prime exist already in
  `K = Q(ζ15, √-11)`:
  - every prime of `K` above 19 is split from its complex conjugate, since `ζ̄^(-1) = ζ̄^4` would force `ζ̄^5 = 1`;
  - so `δ = conj(γ)/γ` with `γ ∈ 𝔮 ∖ conj(𝔮)` is a unit vector with a pole at the three primes of `F` above `𝔮`, one
    of which has type `h ≡ 7`.

  The census above includes a first probe with one such `δ` per conjugate pair. It finds no three-step connector. A
  connector must clear the poles at all 24 primes, so `δ`-steps covering several pairs (two or more `δ`'s, or `δ`'s
  attached to products of primes) are the next family to try.
- **Witnesses.** Inside a coefficient factor's own unit graph, a new connector point with at most five neighbours never
  blocks a 6-colouring. A connector can invalidate saved additive witnesses only through the contact families it creates.
- **Ports versus HYP-9290.** `X(N)` (statement 2) plays the role of HYP-9290's frontier terms (`MF` and the same-step
  merges), re-indexed by backward level instead of synchronous depth. The two decompositions are different bookkeeping
  and are not claimed equal.
- **Relation to THM-4667's ideal ladder.** In the synchronous transport, every source moved one step,
  `E_1 - E_0 = 2 Σ` over depth-one meet pairs (pairs with `P(x) = P(x')`). THM-4667(6) lists this cross term as
  -1.4687. By clock lag it splits as follows (`pair_series_depth1.py`, `N` up to `4·10^6`):
  - lag 1 is exactly the siblings, `-10/7 = -1.4286`, about 97% of the term;
  - lags 5 and 7 give `+2·0.008929` each, matching `1/112` to six digits;
  - lag 6 gives `2·(-0.037946)`, matching `-17/448`;
  - lags 11-13 give about `10^-3` in total.

  These matches are observed, not derived (FINITE-EXACT numerics).
- **Why the denominator is 7.** The chain residues have period 3 because `B^3(x) = 64x + 21 ≡ x + 3 (mod 9)`. These are
  the Solfeggio three-lift fibres (S42 L4). So the geometric sums run over `4^(-3k)`, and
  `Σ_k 4^(-3k) = 64/63 = 64/(9·7)`.
