# Two preprints, the Solfeggio numbers, Kathara and golden clocks, and digit patterns: an analogy atlas around ROOT energy and six-colour forcing

2026-10-17, opus-2026-10-17-S42.

**Owner prompt.**
- Merge in results related to, or creatively adjacent to, `{396, 417, 528, 639, 741, 852, 963, 174, 285}`.
- Keep working on uniform ROOT energy and cross-fibre six-colour forcing.
- Reroute creatively using Kathara and golden clocks.
- Read two attached preprints very carefully and pattern-match liberally, even where logic says there should be no
  connection:
  - OpenAI, *The Quasi-Riemann Hypothesis*, 5 Oct 2026: a zero-free half-plane `Re s > 11/12`;
  - OpenAI, *Conformal flow and the Riemannian Penrose inequality with minimizing frontiers*, 27 Sep 2026.
- Engage the pasted questions on transcendence, digit patterns, normality and elementary approximation.

The concrete mathematical result of the session is the Penrose frontier for the ROOT transport:
[root_frontier_penrose_20261017.md](root_frontier_penrose_20261017.md) and HYP-9290, with ten pre-registered tests
passed to `N < 2^28`. This note is the wider atlas.

**Typing.** Every item carries one of these labels:
- PROVED, with a proof here (all short);
- CITED, with a source;
- CONDITIONAL, on the unrefereed preprint;
- ANALOGY, a structural pattern-match that claims no implication;
- NUMEROLOGY, a coincidence of numbers with no mechanism found.

**Proved here (elementary; repaired after the independent S42 audit, `04-computation/experiments/s42_audit_20261017/REPORT.md`, MISTAKE-628).**
- L1, Lipschitz amplification.
- L2, the transport sign is a Legendre symbol mod 3, and the stopping-time parity is a product of such symbols.
- L3, the Kathara/Fibonacci factorisation as the seed map.
- L4, the mod-9 dictionary (Solfeggio = three-lift fibres, golden classes, CRT, `K_6 = K_(3,3) ⊔ 2K_3`).
- L5, digit corollaries for divergent orbits (chains of cited theorems through THM-4528; item 3 rests on a citation
  not yet checked against the paper).
- L6, the CRT form of six-colourings and the odd-reflection remark.

## 1. The quasi-Riemann preprint, read in full

**What it proves.** Every finite-order Hecke `L`-function over `K = Q(sqrt(-3))`, hence every Dirichlet `L`-function and
`ζ`, has no zeros in `Re s > 11/12`. It is a shorter companion of the `7/8` preprint ("openai/math #003", already
cited by S31 and THM-4566).

**The four steps** (section 2 of the preprint):
1. **Littlewood.** A power saving `A_1(D) ≪ D^(1 - δ + ε)` for smoothed twisted Möbius sums gives a zero-free region.
2. **Embedding in a sextic family.** `A_u(D) = Σ_n μ(n) ν(n) (u/n)_6 W(N n / D)`, with the mean square
   `Σ_(N u <= H) |A_u(D)|^2 ≪ D^(1+ε) H` (rows `u`, columns `n`, `H = D^(1+θ)`).
   - **Amplification.** For primes `p` of norm about `Y`, `A_(p^6)(D) = A_1(D) + O(D/Y)`. The single sum is replicated
     in about `Y/log Y` rows, giving `|A_1|^2 ≪ D^(1+ε) H^(5/6) + D^2 H^(-1/3)`, i.e. exponent `11/12`.
3. **Poisson summation in `u`.** Sextic Gauss sums absorb `μ` (Hasse, Heath-Brown) and become cubic Gauss sums, which
   are coefficients of Kubota's cubic theta function (Patterson). The dual row scale is `H* = D^2/H`.
4. **Cube completion and transfer.** Cube completion `n b^3` (Dunn-Radziwiłł), with Möbius inversion over `b`, is
   followed by a transfer estimate with *power contraction*. A finite recursion closes. Quadratic twists `χ_(3h)` and
   the Goldmakher-Louvel quadratic large sieve over `K` enter at this stage.

**Consequences claimed.**
- Effective `h(D) ≫ sqrt|D| / log log |D|` (Littlewood).
- Least quadratic nonresidue `<= (log p)^C`.
- No Landau-Siegel zeros.
- Completeness of Euler's 65 idoneal numbers (with Elsenhans-Klüners-Nicolae and Grube).
- Corollary 1.2: `π(x; q, a) = li(x)/φ(q) + O(x^(11/12) log x)` uniformly for `q <= x`.

### 1a. Conditional merges into the repository (CONDITIONAL on the preprint)

| repository result | what the preprint would add |
|---|---|
| THM-4666 (S40): the natural order of `P_p` has excess `p h(-p)` | **effective** `p^(3/2)/log log p ≪ excess ≪ p^(3/2) log log p` |
| S34 section 5: `{0, ..., n(p)-1}` is the maximal initial transitive run of `P_p` | that run has length `<= (log p)^C`, at a fixed power, by a zero-free half-plane |
| THM-4666 statement 5: Kalmynin's `L+` primes have no real zero (Fekete) | no prime has a zero in `(11/12, 1)`. `L+` still says more: no real zero in `(0, 1)` |
| THM-4566 and S40 idoneal layers | the 19 idoneal planes and Euler's 65 numbers are complete |
| THM-4631 (Artin multiplicity one) | **nothing new.** Hooley needs square-root Chebotarev errors for about `x^(1/2)` Kummer fields; a fixed half-plane does not supply them (S31, MISTAKE-604) |

**Details.**
- **Paley excess.** The upper bound uses `L(1, χ) ≪ log log q` from the half-plane (Littlewood's short Euler product).
  The lower bound uses the claimed effective `h`. Unconditionally the upper bound is `p^(3/2) log p`, and the lower
  bound is only Siegel's ineffective `p^(3/2 - ε)`; the effective unconditional lower bound (Goldfeld-Gross-Zagier) is
  of logarithmic size in `h`.
- **Initial transitive run.** Under GRH the bound would be `2 log^2 p` (Bach). S34 called this "the GRH channel"; the
  preprint would replace GRH there by its half-plane.
- **Kalmynin's `L+`.** The preprint excludes zeros only in `(11/12, 1)`, while `L+` excludes real zeros in all of
  `(0, 1)`. So the preprint does not reach THM-4666's `L+` statement, which remains the sharper (and finite) assertion.
- **Idoneal completeness.** THM-4566 already records "complete modulo openai/math #003"; the 11/12 paper is a shorter
  route to the same input.

### 1b. Structural pattern-matches (ANALOGY unless marked)

1. **Amplification ↔ Lipschitz amplification (L1, PROVED; constant sharpened after the S42 audit, F8-F9).** The
   current `R(N) = H(4N) - H(N)` is 4-Lipschitz in `N`.
   - The defining difference is
     ```
     R(N+1) - R(N) = Σ_(k=4N+1)^(4N+4) w(k) f(k) - w(N+1) f(N+1).
     ```
   - The last terms cancel: `4N + 4 ≡ N + 1 (mod 3)` gives equal weights, and `f(4m) = f(m)` (two halvings). So
     `ΔR = Σ_(k=4N+1)^(4N+3) w(k) f(k)`.
   - Among three consecutive `k` exactly one has `w = 2`, so `|ΔR| <= 2 + 1 + 1 = 4`. The bound is attained (max 4 over
     `N < 2^26`).
   - So a value `|R(N_0)| = A` is *replicated*: `|R(N)| >= A/2` for `|N - N_0| <= A/8`, exactly as the preprint
     replicates `A_1(D)` in the rows `u = p^6`.

   **Lemma L1.** Let `D(X) = Σ_(X <= N < 2X) R(N)^2`. For `N_0` in `[X, 2X)`,
   ```
   |R(N_0)| <= max( (32 D(X))^(1/3), 2 (D(X)/X)^(1/2) ),     with (32)^(1/3) < 3.18.
   ```
   *Proof.*
   - Let `A = |R(N_0)|` and `h = ⌊A/8⌋`.
   - The integers of `[X, 2X)` within distance `h` of `N_0` number at least `min(h + 1, X)`: either one side has room
     for `h` more, or the whole octave is within distance `h`.
   - On them `|R| >= A - 4h >= A/2`, so `D(X) >= min(h + 1, X) A^2/4`.
   - If `h + 1 <= X`, then `D >= (A/8)(A^2/4)`, so `A <= (32 D)^(1/3)`.
   - Otherwise `D >= X A^2/4`, so `A <= 2 (D/X)^(1/2)`. ∎

   The audit checked the lemma in every octave `2^6-2^25`. The trivial bound is `|R(N)| <= 4N + 2`; the observed
   maximum of `|R|/N` is 2.

   **Consequence.** HYP-9288 (B) bounds the octave traces `sup_d T_d(X) >= D(X)`. So an averaged octave bound
   `D(X) <= X^2 (a + b log X)` already gives the pointwise power saving `|R(N)| ≪ N^(2/3) (log N)^(1/3)`. The second
   term of L1 is then only `O((X log X)^(1/2))`. The averaged statement is far weaker than the pointwise `sqrt(N)`
   target, yet it yields a fixed exponent below 1. The preprint's `11/12` has exactly this shape: replication count
   times `|A|^2` is at most the family mean square. The frontier lemma of the companion note is the `l2` alternative
   that reaches exponent `1/2` pointwise, given the frontier cone.
2. **The family.** The preprint's mean square `Σ_u |A_u|^2 ≪ D^(1+ε) H` says "square-root cancellation on average over
   the rows". The octave trace `D(X) ≍ X^2 = X · X` says the same for the ROOT current over `N` in `[X, 2X)`. The
   rows of the preprint are character twists; ours are the scales `N`.
3. **Duality and the missing contraction.**
   - The preprint's Poisson step sends rows of size `H` to dual rows of size `D^2/H`. Its transfer step contracts the
     ratio `H/X` by a power, so a finite recursion closes.
   - The ROOT transport's Fourier form (THM-4667 statement 10; clock_duality's R8) is the analogous duality. There each
     odd transport *triples* the ternary modulus and the band period (`q_d = 3^(d+1)`, `P_d = 2·3^d`). Nothing
     contracts.
   - This is the depth obstruction, stated in the preprint's language: we have the duality but not the power
     contraction. The frontier route avoids needing it by stopping at a fixed scale.
4. **Cube completion ↔ seeds and B-chains.**
   - The preprint completes squarefree `n` to `n b^3` and removes `b` by Möbius inversion.
   - The ROOT transport reduces every source to its seed: it strips the chain `s -> 4s + 1 = B(s)`, whose members share
     an odd image, `3B(s) + 1 = 4(3s + 1)`.
   - clock_duality's R8.2 shows that the three last-digit lifts of every residue sum to zero. This is a Möbius-type
     cancellation on the completed chain (`B^(3^d)` is a translation by a multiple of `3^d` modulo `3^(d+1)`).
   - ANALOGY: the seed is the "cube-free part", and the B-chain is the completion.
5. **Sextic symbol ↔ transport sign and six colours (L2, PROVED).**
   - The preprint's twist is `(u/n)_6 ∈ μ_6 = μ_2 × μ_3`. Its quadratic part is absorbed by Poisson (Gauss-Jacobi),
     leaving cubic Gauss sums; quadratic twists return in the transfer.
   - The ROOT transport sign is itself a quadratic character: see L2 in section 4. Six colours are `Z/6 = Z/2 × Z/3`
     (L6, section 7).
6. **Collatz has no quasi-RH for its parity sum, and that is the point** (object corrected after the S42 audit, F10).
   - THM-4659: the zeros of `D(s) = 1 + 2^(-s) + (3/2)^s/3` approach `Re s = 1` at every parity-matched near-return of
     `log2 3`, with `1 - Re ρ ≈ 199.165 ε^2`. Such near-returns with `ε -> 0` exist because `(2r+1) log2 3 mod 2` is
     dense.
   - So `D` has no zero-free strip of any width to the left of `Re s = 1`. Under THM-4660's hypothesis, its zeros are the
     *poles* of the parity Dirichlet series `F = G/D`.
   - THM-4660: if the mod-3 discrepancy `H` has square-root energy, then the second Riesz mean of the parity sum is an
     explicit sum over these zeros plus `O(x^(1/2 + δ))`.
   - So the analogue of the preprint's theorem is not a half-plane for `D` or `F`. It is a square-root ("RH-strength", not
     merely quasi) bound for the numerator side `H`, whose two-octave differences `H(4N) - H(N)` are the ROOT current.
     HYP-9288 and HYP-9290 are exactly such statements.
   - The preprint shows that a fixed-width zero-free half-plane is reachable for `ζ` and all Dirichlet `L`. For Collatz
     the corresponding object `D` provably has no zero-free strip, and the Diophantine near-returns are the reason.
7. **The Eisenstein base field and THM-4285.** The preprint lives in `O = Z[ω]`: primary `n ≡ 1 (mod 3)`, and six units
   representing `(O/3)^×`.
   - THM-4285 (codex, the JC(2) programme) computes the theta series of `A_2 ⊕ A_2`, two Eisenstein norms, with
     coefficient `12(σ(n) - 3σ(n/3))`, and splits its contacts into norm and non-norm quotients (`5Q` against `3Q`).
     That is the split/inert dichotomy of `χ_(-3)`, the character of the preprint's base field.
   - The golden six-colour field `K = Q(ζ15, sqrt(-11))` contains `Q(ζ3) = Q(sqrt(-3))`. So the preprint, THM-4285 and
     the golden contact share one base field. ANALOGY; no transfer claimed.
8. **NUMEROLOGY.**
   - `11/12` is also the error exponent in THM-4028 (Sun's 2-4-6-8 count, `V X^(25/24) + O(X^(11/12))`).
   - `7/8` is the half-plane used by the owner's Ramanujan-inequality manuscript (S31).
   - No mechanism is found in either case.

## 2. The Penrose preprint: the refined dictionary

The frontier note gives the transfer. Re-reading section 1.2 of the preprint ("four quantities that control the
proof": area `A`, mass `m`, deficit `μ = m - c`, outer radius `R_t`) sharpens one row of that dictionary:

| conformal flow | stopped ROOT transport |
|---|---|
| area `A`, fixed by the flow | `R(N) = Σ f(u) c(u)`, fixed by every transport step |
| mass `m(t)`, decreasing | energy `E_d`, decreasing in the ideal regime |
| capacity `c_t` = the minimal Dirichlet energy carried by the frontier | `R(N)^2 / F` = the minimal energy carrying the invariant on `F` seeds (the Cauchy-Schwarz minimiser) |
| mass-capacity inequality `μ = m - c >= 0` | **Lemma F**: `E_Y - R^2/F_Y >= 0` |
| `m' = -2μ`, the deficit drives the mass down | `E_(d+1) - E_d = 2 Cross_d`; measured: merges above `Y` remove 54-76% of `E_0` |
| outer radius `R_t`; enclosing spheres | the frontier height `Y`, i.e. `F_Y`; the final Cauchy-Schwarz |
| odd reflection `ψ = ±φ` across the doubled seam | sibling pairs `{s, 2s+1}` with opposite clock parities (THM-4667 statement 8; clock_duality's R5), and, in six colours, `c_2 = -c_1` (section 7) |
| singular frontier; "regularity at a smooth obstacle" | ROOT and the few small seeds where arrivals are coherent (frontier decomposition: the large climbs are arrival terms below `Y ≈ 30-100`, apart from transient bumps) |

**Lemma F as mass-capacity.** In the preprint, mass-capacity is the hard positive-mass input, proved by doubling. Here
it is one line of Cauchy-Schwarz, and the hard part moves to the energy bound (HYP-9290).
- The equality case `m = c` corresponds to `c_Y(u) = (R/F_Y) f(u)` on every frontier seed, with `E_Y = R^2/F_Y`: the
  **maximally spread** state (corrected after the S42 audit, F11).
- The one-seed state `E = R^2` is only the trivial case `F = 1`, the end of the unstopped flow.
- In basin form (frontier note, section 2), equality means that all basin currents `R_u` are equal.

## 3. Kathara, Fibonacci and golden clocks are one object: the golden-mean shift

**(L3, PROVED; the decomposition itself is THM-4667 statement 10 and clock_duality's even-cost register.)**
- **Every halving exponent splits.** An exponent `A >= 1` is uniquely `A = a + 2j` with `a ∈ {1, 2}` and `j >= 0`. In
  generating functions:
  ```
  z/(1 - z) = (z + z^2) · 1/(1 - z^2).
  ```
  The left side is the Kathara letter alphabet (clock_duality's (9)-(10): `a_n = 2^(n-2)` counts compositions). On the
  right, `z + z^2` is the Fibonacci restriction (parts `{1, 2}`), and `1/(1 - z^2)` is the even register.
- **The seed map realises this split.** Seeds have `a ∈ {1, 2}` only, since `3s + 1 ≡ 4, 2, 6 (mod 8)` for
  `s ≡ 1, 3, 7`. Each strip `u -> (u - 1)/4` of the B-chain spends two halvings, so `j` counts strips. So the seed
  transport runs on Fibonacci words, and the B-chain is the even register.
- **Parity words are golden-mean words.** THM-4528 (S15) proves that parity words never contain `11`. Here `T` is the
  non-shortcut map `n/2, 3n+1`; L2 in section 4 uses the shortcut map `n/2, (3n+1)/2`. Read as blocks `0` and `10`, they are compositions with parts `{1, 2}`: the same Fibonacci restriction again.
- **The golden clock is this shift's transition matrix.** The golden-mean shift has transition matrix
  `Q = [[1,1],[1,0]]`, clock_duality's golden clock. The number of `T`-parity words of length `l` is a Fibonacci
  number. Modulo 9 these counts repeat with Pisano period 24, and `Q^12 = -I (mod 9)` makes them change sign after 12
  steps. (O4's norm classes are the orbits of `Q` on primitive pairs.)

**Consequence for the frontier.** The stopped transport is a walk on Fibonacci words.
- Its path counts grow like `φ^t`, while Haar mass decays like `2^(-t)`. The ratio, `φ/2 = cos 36°`, is also the first
  point of the trivial cycle in THM-4528's golden reading (section 6). That coincidence is NUMEROLOGY.
- ANALOGY: the coalescence that the frontier energy measures is a property of the golden-mean shift with the
  Collatz weights.

## 4. Solfeggio, Kathara and golden classes modulo 9 (L4, PROVED; mostly elementary, labels mark the rest)

**The Solfeggio set.** The nine numbers are `396, 417, 528, 639, 741, 852, 963, 174, 285`. Each lists, in its decimal
digits, one fibre of `Z/9 -> Z/3`, descending by 3 (`0` written as `9`):
```
396: 3, 0, 6   639: 6, 3, 0   963: 0, 6, 3     (fibre {0, 3, 6}: the non-units)
417: 4, 1, 7   741: 7, 4, 1   174: 1, 7, 4     (fibre {1, 4, 7}: the principal units 1 + 3Z/9)
528: 5, 2, 8   852: 8, 5, 2   285: 2, 8, 5     (fibre {2, 5, 8} = 2·{1, 4, 7})
```
- **B³-orbits and three-lifts.** With `B(x) = 4x + 1`, the seed chain of the transport, `B^3(x) = 64x + 21 ≡ x + 3
  (mod 9)`. So each Solfeggio triple is a `B^3`-orbit, i.e. the three lifts of one residue mod 3. These are the triples
  whose coefficients sum to zero in clock_duality's R8.2, at the first level `d = 1`.
- **Arithmetic.** Since `10 ≡ 1 (mod 9)`, every Solfeggio number is `≡ 3a (mod 9)` for its first digit `a`. So all nine
  are divisible by 3, and exactly `396, 639, 963` by 9.
- **Kathara mod 9 is the doubling circuit.** `a_n = 2^(n-2)` runs through `1, 2, 4, 8, 7, 5` (2 is a primitive root
  mod 9). It never meets `3, 6, 9`, the non-unit fibre, and it alternates between the fibres `{1,4,7}` and `{2,5,8}`
  because `2 ≡ -1 (mod 3)`.
- **Golden classes are orthogonal to Solfeggio fibres (CRT).** The three `Q`-orbits mod 9 are the norm classes
  `{±1}, {±2}, {±4}`, the cosets of `{±1}`. The Solfeggio unit fibres are the cosets of `{1, 4, 7}`. Since
  `(Z/9)^× = <2> ≅ Z/6 = Z/2 × Z/3`, the unit `2^k` is determined by the pair (`k mod 2` = Solfeggio fibre,
  `k mod 3` = golden class):

  | `k` | 0 | 1 | 2 | 3 | 4 | 5 |
  |---|---|---|---|---|---|---|
  | `2^k mod 9` | 1 | 2 | 4 | 8 | 7 | 5 |
  | Solfeggio fibre | 417 | 528 | 417 | 528 | 417 | 528 |
  | golden class | ±1 | ±2 | ±4 | ±1 | ±2 | ±4 |

  - The linear part 4 of `B` is `2^2`. It fixes each Solfeggio fibre and cycles the golden classes. This is O4's
    remark that doubling sends a norm to four times itself.
- **The colour graph.** Use the six units as six colours. Then `K_6 = K_(3,3) ⊔ 2K_3`, edge-disjointly:
  - the two triangles are the unit Solfeggio fibres `{1,4,7}` and `{2,5,8}` (same `Z/2` coordinate, different `Z/3`);
  - the `K_(3,3)` joins them (different `Z/2` coordinate).

  This is the Z/2 × Z/3 split of a six-colouring (section 7).
- **Fibonacci mod 9 (NUMEROLOGY made exact).** The 24-term Pisano cycle of digital roots is
  `1,1,2,3,5,8,4,3,7,1,8,9, 8,8,7,6,4,1,5,6,2,8,1,9`.
  - `Q^12 = -I` pairs `F_n` with `F_(n+12) ≡ -F_n`, so the two halves are complementary (`1+8, 2+7, 4+5, 3+6, 9+9`).
    The unit pairs are exactly the golden classes.
  - The non-unit digits `3, 6, 9` occur exactly at `n ≡ 0 (mod 4)`, because `3 | F_n ⟺ 4 | n`.

**L2 (PROVED): the transport sign is a Legendre symbol mod 3.**
- For any odd `m`, write `3m + 1 = 2^a u` with `u` odd. Then `u ≡ 2^(-a) ≡ (-1)^a (mod 3)`, so
  ```
  (-1)^a = (u/3).
  ```
  The ROOT transport sign `σ = (-1)^a` (sign `+` exactly when `a = 2`) is therefore the quadratic character mod 3 of
  the unstripped image.
- Under the shortcut `T`, the step from an odd `m` to the next odd orbit element `u` takes `a` steps. Iterating gives,
  for odd `n`,
  ```
  (-1)^(σ_T(n)) = Π (u/3),   the product over the odd orbit elements u ≠ n of n, down to 1.
  ```
  Examples:
  - `n = 3`: the product is `(5/3)(1/3) = -1`, and `σ_T(3) = 5`;
  - `n = 7`: the odd elements are `11, 17, 13, 5, 1`, the product is `-1`, and `σ_T(7) = 11`.
- So the clock parity `f` is a product of quadratic characters along the orbit. The ROOT current
  `R = Σ w f` is a twisted sum of such products, with `w = 2, -1, -1` on residues `2, 0, 1 mod 3`. This is the
  analogue of the preprint's quadratic twists `χ_(3h)`. The identity is elementary and is probably implicit in the sign
  conventions of THM-4658 and clock_duality's R8; it is recorded here because it makes the mod-3 nature of the parity
  explicit.

## 5. The nine numbers in the repository (merge of related results)

Theorem files carrying these numbers, and the adjacency each has to the present threads. The mathematics in the
second column is the cited theorem's; the third column is new and typed.

| number | repository results | adjacency |
|---|---|---|
| 396 | THM-396: a universal blocker of a small `n = 14` pinch is a sum-multiple shield. THM-4396: LRC(14) finite dual / exact-pair certificate | **ANALOGY.** A *fixed-scale blocker* that settles a whole family, the role the frontier plays for the ROOT climb |
| 417 | THM-417 (with THM-413/415): the signed-LRC sign orbit is full iff `C = 2n - 1` is prime. THM-4417: width-two Laurent first-return parabolic bound | **EXACT link.** At `C = 9` (`n = 5`) the only proper subgroup is `K_3 = {0, 3, 6}`, the Solfeggio 3-6-9 fibre, and the unique collision is the single `H_3` flip (THM-417: deficiency `2^(q-3) = 1`). At `C = 27 = 3^3`, the shell modulus of LRC(14) (THM-407), the nested `K_3 ⊂ K_9` give deficiency 69. This is the three-lift tower of section 4 and of R8.2, in a third thread |
| 528 | THM-528: the LRC four-window lemma (windows at `0, 1/2, 1/3, 2/3`). **THM-4528: Collatz parity words are golden-mean words** | THM-4528 *is* the bridge between the golden clocks (section 3) and the digit-pattern questions (section 6). Its trivial cycle reads `{cos 36°, cos 72°, cos 60°}`, the half-lengths `φ/2, 1/(2φ), 1/2` of the golden pentagon whose diagonals are the `d`-steps of the six-colour ports (NUMEROLOGY for the six-colour link) |
| 639 | THM-639: Hamiltonian-path frame for runner families. The gap law is Haar on the annihilator loop; the AP is the unique wall-count minimiser ("the transitive tournament of the runner world") | **ANALOGY.** Haar on a loop is the ideal regime of the transport. The AP-as-transitive-order matches THM-4666's natural order of `P_p` |
| 741 | THM-741: near-AP four-slot closure, all 2002 = `C(14,9)` bodies | |
| 852 | THM-852: Klein-four structure; `n = 8` refutation of the black self-line SC law | |
| 963 | THM-963: scale-eight Hamming-six `K_(3,3)` owner-nerve obstruction | **ANALOGY.** `K_(3,3)` is the `Z/2`-part of the six-colour graph `K_6 = K_(3,3) ⊔ 2K_3` (section 4) |
| 174 | THM-174: `det(I + 2A) = Pf(S)^2` for tournaments (even `n`), because `1^T S^(-1) 1 = 0`. THM-4174: six-deletion completion of newcomer Haar transfer | **ANALOGY.** The all-ones vector is invisible to the skew part, as R8.2's three-lift sums vanish |
| 285 | THM-285: closure of the `β_2 = 0` gap at `n = 4, 5`. THM-4285: Eisenstein norm shells (`A_2 ⊕ A_2` theta), contact trichotomy | **Same base field.** THM-4285 shares `Q(sqrt(-3))` with the quasi-RH preprint and the golden contact field (section 1b item 7) |

Further notes on the table:
- **396.** THM-396 and THM-4396 are LRC(14) results; no transfer is claimed.
- **639.** No transfer is claimed.
- **741, 852.** No adjacency was found beyond the LRC(14) programme.
- **963.** The `K_(3,3)` is the same graph in both places; nothing more is claimed.

**NUMEROLOGY of the numbers themselves.**
- `417 = 3 · 139`, and 139 is the second cycle-count prime of S32 ("3 vs 1 cycles = Gersonides + 139").
- The pairs `396, 528` (both divisible by 11), `639, 852` (by 71) and `741, 285` (by 19) share a prime.
- `528 - 396 = 12 · 11`, `852 - 639 = 3 · 71` and `741 - 285 = 24 · 19`.
- No mechanism.

## 6. Digit patterns and algebraic properties (the pasted questions)

The repository already owns the key bridge, **THM-4528**. With `F_n(z) = Σ_j b_j(n) z^j` the parity series:
- (G2) the orbit of `n` reaches 1 ⟺ `(1 - z^3) F_n ∈ Z[z]` ⟺ `7 Φ_T(n) ∈ Z` ⟺ `2Θ(n) ∈ Z[φ]`;
- (G3) the orbit is bounded ⟺ `F_n` is rational ⟺ `F_n` is D-finite.

(G3) holds by Szegő's theorem: a power series with finitely many distinct coefficients is rational or has the unit
circle as a natural boundary.

So "no divergent orbit" is a statement about the *algebraic nature of a digit generating function*. This is exactly
the kind of theorem the pasted block asks for.

**L5 (items 1-2: PROVED from THM-4528 and K. Schmidt's theorem; item 3: CITED, with the β-expansion citation
unverified).** For `n >= 1`:
1. **Boundedness is a field condition.** The orbit is bounded ⟺ `Θ(n) ∈ Q(φ)`.
   - The parity word of `n >= 1` is the greedy base-`φ` expansion of `Θ(n)`. It contains no `11` (THM-4528 G1). It
     also does not end in `(10)^∞`: such a tail would force `T^k(n) = -1`, the 2-adic point following `(10)^∞`.
     (Golden-mean words in general need not be greedy: `0.(10)^∞ = 1`.)
   - For the Pisot number `φ`, every element of `Q(φ) ∩ [0, 1)` has an eventually periodic expansion (K. Schmidt
     1980; classical for `φ`). The converse is a geometric sum.
   - An eventually periodic parity word means an eventually periodic orbit, by 2-adic uniqueness of the point
     following a periodic word (Lagarias), which is the same as bounded.
2. **A divergent orbit gives a number outside `Q(φ)`.** If the orbit of `n` diverges, then `Θ(n) ∉ Q(φ)`, and the
   binary reading `x(n) = Σ b_j 2^(-j-1)` is irrational.
3. **Low complexity forces transcendence.** If, moreover, the parity word of a divergent orbit has subword complexity
   with `liminf p(k)/k < ∞` (for example Sturmian, or automatic), then:
   - `x(n)` is transcendental (Adamczewski-Bugeaud, Annals 2007);
   - `Θ(n)` is transcendental, by Adamczewski-Bugeaud, *Dynamics for β-shifts and Diophantine approximation*, ETDS 27
     (2007) 1695-1711, for a Pisot base.
     - The statement used is: for `β` Pisot or Salem and `ξ` algebraic with `ξ ∉ Q(β)`, `liminf p(n)/n = ∞`.
     - It is quoted via the subspace-theorem survey and is UNVERIFIED against the paper; it must be checked before any
       further use.

   A divergent Collatz orbit with a *simple* parity pattern would manufacture transcendental numbers. One with a
   *complex* pattern is unconstrained by these theorems.

**The pasted questions, answered with theorem pointers.**

1. *"One function is transcendental iff a logarithmic perturbation dovetail around every center is normal."*
   - Rigorous relatives exist for series with finitely many coefficients. By Szegő, `F_n` is rational (algebraic) iff
     it continues analytically across *one* arc of the unit circle. It is transcendental iff *every* boundary point
     ("every center") is singular.
   - For entire functions, transcendence ⟺ an essential singularity at ∞. In logarithmic coordinates `w = log z`,
     rescalings `z -> 2^k z` become translations, and non-normality of the translated family on an annulus is a
     Montel/Julia-type signature of the essential singularity.
   - ANALOGY in our setting: THM-4632 shows that the Collatz transfer function `R(s)` has poles on its line iff
     `log2 3` were rational, and otherwise boundary growth exponent `2(μ - 1)`. The algebraic nature of a number is
     read from the boundary behaviour of a family built from one elementary function.
2. *"More theorems connecting digit patterns to algebraic properties."* Here is the list we would use (all CITED):
   - rational ⟺ eventually periodic (any base);
   - Szegő (1922) and Pólya-Carlson (1921): integer coefficients and radius 1 give rational or a natural boundary;
   - Christol: algebraic over `F_p(z)` ⟺ `p`-automatic coefficients;
   - Adamczewski-Bugeaud: automatic irrational ⟹ transcendental, and algebraic irrational ⟹ complexity `p(n)/n -> ∞`;
   - Bailey-Borwein-Crandall-Pomerance (JTNB 16, 2004): an algebraic irrational of degree `D` has at least
     `c N^(1/D)` ones among its first `N` binary digits;
   - Ridout and the subspace theorem (long zero blocks), already used in S15b;
   - Mahler's method (lacunary and automatic numbers);
   - Borel: almost every number is normal; normality of algebraic irrationals is open;
   - Stoneham and Bailey-Crandall: normality in one base but not another, and normality of BBP constants under
     "Hypothesis A".
3. *"If the generating function of the digit sequence has an elementary inverse, does that constrain normality?"*
   Yes, strongly, in the direction **away** from normality:
   - **Characteristic 0.** Suppose the digit generating function `Σ a_n z^n` (bounded integer digits) is D-finite (for
     example algebraic). It then has finitely many singularities and continues across an arc of the unit circle, so by
     Pólya-Carlson the number is rational. Any generating function that continues across an arc suffices.
   - **Characteristic p.** If it is algebraic over `F_p(z)`, then by Christol the digits are automatic. The subword
     complexity is then `O(n)`, so the number is **never disjunctive and never normal**, and by Adamczewski-Bugeaud it
     is rational or transcendental.
   - **Compositional inverse.** A compositionally elementary-invertible `A` blows up near `|z| = 1` like
     `(1 - |z|)^(-1)`. Its behaviour there is decided by the elementary inverse's behaviour at ∞, so the question
     reduces to that function's singularity at infinity. No general theorem is claimed.
   - **The constructive side.** Normal numbers are built by the opposite recipe: concatenation (Champernowne,
     Copeland-Erdős) or lacunary-geometric series with a chaotic shift map (Stoneham). These are generating functions
     with natural boundaries.
4. *"How does the degree of the approximating circuit grow near certain points? Can a number be proved transcendental
   because it is too well approximated by elementary functions?"*
   - **There are two dual families of theorems.**
     - *Too well approximated by rationals ⟹ transcendental*: Liouville, Roth, Mahler's S/T/U classes. Adamczewski-
       Bugeaud is the automaton (finite-circuit) version. The **Hartmanis-Stearns conjecture**, "real-time computable
       ⟹ rational or transcendental", is the exact open statement of the owner's intuition for circuits.
     - *Too well approximated by rational functions ⟹ the function is rational*: Kronecker (Hankel determinants),
       Borel-Dwork, Pólya-Carlson, Szegő. THM-4528 uses this family, and S15b's "Borel-Dwork 1" criticality is the
       same family.
   - **Hermite's proof of the transcendence of `e`.** It is the prototype of "degree grows with precision". The `[n/n]`
     Padé approximant of `e^z` has error `(n!)^2/((2n)!(2n+1)!) z^(2n+1)`.
   - **In our setting.** THM-4632 is a theorem of exactly the asked shape. Near the frequencies of the convergents of
     `log2 3`, the elementary function `R(s)` grows with exponent `2(μ - 1)`. S29's rational approximants along the
     convergents need degree about `q_n` to resolve frequency `≈ q_n`. The irrationality measure of `log2 3` is the
     growth rate of the approximating degree.

## 7. Six colours: the CRT form, odd reflection, and what cross-fibre forcing must defeat

**L6 (PROVED, elementary).**
- **A six-colouring is two colourings.** Write `c = 3c_2 + 2c_3 (mod 6)`; this is a bijection
  `Z/2 × Z/3 -> Z/6`. A map `c` into `Z/6` is a proper colouring iff every edge is properly coloured by `c_2` *or* by
  `c_3`. Equivalently, the edge set is covered by a bipartite spanning subgraph and a 3-colourable one. This is the
  product bound `χ(G_1 ∪ G_2) <= χ(G_1) χ(G_2)`, read as an equivalence for 6 = 2 · 3.
- **The additive problem splits the same way.** The additive colourings of HYP-9289,
  `c = c_A + c_1 + c_2 (mod 6)`, split into an additive `Z/2` problem and an additive `Z/3` problem. A mixed edge is
  satisfied iff its shift is nonzero in at least one component.
- **For forcing, both components must fail on the same edge.** An additive obstruction must make every pair of
  component solutions leave a common unsatisfied edge. It is not enough to block the `Z/2` problem (always solvable on
  bipartite pieces) or the `Z/3` problem separately.
- **Odd reflection kills shifts.** In S41, `c_2 = -c_1` with homomorphism-like coefficient colourings makes every
  contact shift zero. Then the additive problem collapses to plain 6-colouring of the golden two-distance graph, and
  465 of 570 patterns fail. The even choice `c_1 = c_2` frees everything.
  - This is literally Bray's doubling with the odd reflection `ψ = -φ`, which forces the doubled function to vanish on
    the seam. In the colouring it forces the coupling to vanish on the contact.
  - ANALOGY, with an exact computation on the colouring side.
  - Lesson for forcing: a coefficient graph whose *every* proper colouring is "odd-reflection-like" on the contact
    directions, i.e. rigid with forced zero shifts, would reduce forcing to the two-distance graph. That graph already
    blocks 465 patterns.
  - The zero-shift statement is exact for the `v`-contact in `A + uC + vC`. With the reflected and `e^2` contacts the
    shifts involve further layers.
- **Where the cross-fibre edges are.**
  - THM-4668: the degree-192 direction `e` is Cartesian over `K + uK + vK`.
  - S41's product census found exactly two cross-unit families, the reflected `uv` contact and the `e^2` contact. They
    are the only known cross-fibre edges.

**S42 runs (FINITE-EXACT; CaDiCaL 1.9.5; HYP-9289 has the table).**
- **All three families coupled.** The five-coordinate source `A' + uC' + vC' + (uv)C' + e^2 C'` blocks **0 of 570**
  patterns, with shared and with independent coefficient colourings, for:
  - base × base;
  - base × MM;
  - dM × MM (144,716 directed `v`-contact patterns);
  - contact1 × MM;
  - contact2 × base. This is S41's unfinished recipe: `|A'| = 10,597`, a CSP of `1.12 × 10^7` clauses, solved in
    2,203 s (shared) and 1,961 s (independent).
- **Caveat** (broadened after the audit, F16). These graphs are subgraphs of the unit-distance graphs on the same
  addresses. No census of differences with both `uv`- and `e^2`-parts has been made. Beyond the base difference sets,
  even the single-family edge sets are not known to be complete. A block would still need certification on the full
  graph.
- **Reflection-forcing experiment.** Enlarging the `uv`-layer to `W = C' + conj(d)({0} ∪ U)` (211 points), and once more
  (2,113 points, 460-896 `conj(d)q`-pairs per direction), still blocks 0 of 570.
  - In the solver's colouring of the 2,113-point layer, **4 of the 22 directions show all five nonzero differences**.
    So that solution has zero shift forced on those directions.
  - The other 18 directions leave a nonzero difference value missing.
  - The difference-forcing target is therefore partly met on individual directions. It is not met uniformly.
- **Fixed coefficient colourings are not enough** (sampled). Over 400 coefficient colourings produced by the solver,
  none extends more than 120 of the 570 base patterns; the union is also 120. The joint additive solver extends all 570.
  The additive freedom is pattern-adaptive in the coefficient layers.
- **Design requirement recorded.** Forcing needs a layer whose every proper colouring shows all five nonzero
  differences in *every* contact direction (a difference-forcing gadget). Alternatively, it needs a no-go theorem for
  these golden sources.

## 8. What is new here, and what is credited

**New (elementary).**
- L1, Lipschitz amplification. `R` is 4-Lipschitz (proved), and `|R(N_0)| <= max((32 D)^(1/3), 2 (D/X)^(1/2))`, so averaged
  octave bounds give pointwise `N^(2/3)`.
- L2, the transport sign equals `(u/3)`, and the clock parity is a product of Legendre symbols mod 3 along the orbit.
- L3, the identification of the Kathara, Fibonacci, THM-4528 golden-word and golden-clock structures as one
  golden-mean shift.
- L4, the mod-9 CRT dictionary, and the exact `C = 9` / `C = 27` link to the signed LRC (THM-413/417).
- L5, the digit corollaries for divergent orbits. Items 1-2 are proved; item 3 relies on a citation not yet checked
  against the paper.
- L6, the CRT form and the odd-reflection remark.

**Credited.**
- THM-4528 (S15) for the golden reading and Szegő.
- THM-4667 statement 10 and clock_duality's even register for `A = a + 2j`.
- clock_duality's R8.2 and O4.
- THM-413/415/417 (monad-explorer).
- THM-4285 (codex).
- THM-4566 and S31 for the earlier conditional uses of the quasi-RH preprints.

**Conditional.** Section 1a, on the unrefereed preprint.

**Tests run** (`preprint_connections_checks_20261017.py`, output `.out`/`.json`).

1. **L2.** Zero violations among all odd `n < 10^6` (exact).
2. **Fibonacci mod 9.** The digital-root cycle, the complementary halves and the positions of 3, 6, 9 are as stated.
3. **L1 against the data, every `N < 2^26`.**
   - The maximum of `|R(N+1) - R(N)|` is 4. This is now the proved bound (audit F9; the first draft proved 8).
   - The amplified bound `(32 D(X))^(1/3)` exceeds the octave maximum of `|R|` by factors of about 3.7 (`2^6`) to 29
     (`2^25`).
   - That growth is roughly `X^(1/6)`, the predicted looseness of exponent `2/3` against the true `1/2`.
4. **Character structure of the frontier coefficients** (FINITE-EXACT at `N = 7910, 66463, 2^20, 6904411,
   133970036`; `Y = 100, 1000`; numbers corrected after the audit, F13-F14).
   - **Correlation with the Legendre symbol.** The statistic `ρ = Σ (u/3) c_Y(u) / sqrt(F_Y E_Y)` is small:
     `|ρ| <= 0.16` at `Y = 100` and `<= 0.08` at `Y = 1000`. The true cosine with the character vector is `1.23 ρ`,
     i.e. at most 0.19 and 0.09.
   - **The Solfeggio fibres are not equal.** At `Y = 1000` the energy shares against the seed shares are:

     | fibre of `u` mod 3 | energy share | seed share |
     |---|---|---|
     | `{0, 3, 6}` | 0.38-0.45 | 0.336 |
     | `{1, 4, 7}` | 0.37-0.43 | 0.333 |
     | `{2, 5, 8}` | **0.17-0.23** | 0.331 |

     All five values of `N` agree. At `Y = 100` the deficit holds at 4 of the 5 values; at `N = 7910` the 2-5-8 share
     is 0.343 against a seed share of 0.316.
   - **Mass, not coherence** (`root_transport_frontier_fibre_mass_20261018.py`, FINITE-EXACT at `N = 2^20, 6904411,
     2^24`).
     - The energy share of each fibre tracks the share of sources that *first land* there: energy share divided by
       landing share lies in 0.74-1.17.
     - The landing shares are 0.414, 0.361, 0.225 at `Y = 1000` and 0.540, 0.289, 0.171 at `Y = 100`, for the fibres
       0, 1, 2.
     - The diagonal (no-interference) energy is exactly `8N`.
   - **Why the 2-5-8 fibre receives less** (the structural fact is PROVED; that it yields the observed shares is a
     HEURISTIC mechanism).
     - Seeds `u ≡ 2 (mod 3)` have the singleton preimage `(2u - 1)/3 < u` (THM-4667 statement 8). For the other
       classes the least preimage exceeds `u`.
     - In the seed transport the heads of `u` sit at `4^(j+1) u/3` with `j ≡ 1 - u (mod 3)`. So the first head of a
       2-fibre seed is at `64u/3`, of a 0-fibre seed at `16u/3`, of a 1-fibre seed at `4u/3`.
     - The singleton lies below every frontier containing `u`, so that branch never arrives.
     - A harmonic first-entry model (mass `∝ 1/s`) predicts 16-20% landing on the 2-fibre. It does not robustly order
       the 0 and 1 fibres.
   - **Correction of the first draft's branching table.** The first draft read the fibres as the branching types of the
     *classical* odd Collatz tree (3-6-9 = leaves). In the *seed transport* measured here, every seed `u ≡ 0 (mod 3)`
     has preimages through `B(u) = 4u + 1`, the least near `16u/3` (e.g. `P(17) = P(35) = 3`). The 0-3-6 fibre carries
     the largest share. Only the 2-5-8 row of that reading shows up in the frontier energy.
   - **Golden classes mod 9: a weak but consistent pattern.** At `Y = 1000` the `±2` class is above its seed share
     (0.221) at all five `N` (0.234-0.274), and `±1`, `±4` are below at all five.
5. **Six colours.** See section 7 (S42 runs).
