# Odd primes and their powers: the (p, q, 1^k) defect family of HYP-9340, the disjoint-orbit lemma, and the fold reading of the mirror remarks — mac-mini-2026-11-04-oddprimes

**Session:** mac-mini-2026-11-04-oddprimes (worktree `math-wt-oddprimes-20261104`, branch `codex/session-oddprimes-20261104`).
**Owner directive (2026-11-04):** prove HYP-9340 beyond Sylow (first open case `Z_3` at `n = 6`), find a uniform rule or
definability obstruction for the prime-circulant bijection, prove or disprove HYP-9172, think in odd primes and their
powers; the arc-connectivity shortcut and `{-3}`-self-dual decomposable tournaments; the cycle-count table with the
odd-indexed Fibonacci numbers `2, 5, 13, 34, 89`; the unit-ball volume peak and the `8 : 4` split at `n = 5`; the past
work around `-7/4` and `-29/16`; the MRB constant and the microcosm-macrocosm reading of "the same value realised twice
on either side of the peak"; the uploaded PDF *The MRB Multiverse* (two constants `C` and `B`, mirror bases, Lambert W).
**Truth typing:** PROVED (THM-4871 statements 1-4) + FINITE-EXACT (64 census rows) + REFUTED (HYP-9344 (2)) + NUMERICAL
(HYP-9463 (1)-(3)) + OPEN (HYP-9463 (4)) + CITED/typed (the fold section, the PDF).

## 0. Concurrency and inheritance

* The Sylow/`Z_3` part of the directive was already settled on 2026-10-19/20 (same wording, sessions evengraphs, hyp9172,
  oddprimes, S45 and the Codex packets): HYP-9340 (a) is REFUTED in general (C15 of type `(3,5,1,1)` at `n = 10`, `32`
  tournament types against `31`) and PROVED for every odd `p`-group; HYP-9172 (a) is PROVED for every `n >= 5` (THM-4722
  and the uniform obstruction packet); the affine and bounded-degree rule classes for prime circulants are settled
  (Codex); THM-4728 gives the defect identity and the census through `n = 13`; HYP-9344 (2) pre-registered the next C15
  value, `defect(5) = 88`, as "too large" (`2^24` objects) for the Python census. **Nothing of that is re-derived here.**
* The arc-connectivity shortcut and the few-cycle Fibonacci law `u_3(n) = F(2n - 9)` (klein's HYP-9459 (2), exactly the
  owner's `2, 5, 13, 34, 89`) were reserved by **opus-2026-11-04-S75** (THM-4869, HYP-9461) an hour before this session;
  this session claimed its lanes by message (MSG-3303 broadcast) and does not touch them. The `{-3}`-self-duality
  statement the owner quotes ("a decomposable tournament with at least 9 vertices is `{-3}`-self dual iff it is
  transitive or almost transitive") is a literature pointer (the Boudabbous-Lopez line on `(-k)`-self-duality of
  tournaments); it is passed to S75 and not verified here.
* The PDF's mathematical content (mirror bases `r_n`, `r_n = exp(-W_0(-ln n / n))`, the Cayley/tree-function expansion
  `r_n = sum (m+1)^(m-1)/m! (ln n/n)^m`, `B = 3 - length(E)`, `C = sqrt 2 - 1 - int_E f'`, the absolutely convergent
  bridge `B - C`) coincides with what session31 (`session31_mirror_mrb_20261022`) already proved on 2026-10-22
  (`R(z) = -W_0(-z)/z`, the Galton-Watson reading `r_n/n` = extinction probability of Poisson(`ln n`), the Gelfond-Schneider
  classification with `2 <-> 4` the only rational pair). Verified again here numerically (section 5); nothing new claimed.
* **Closest proved mechanism:** THM-4721 Theorem R (normaliser orbits on tournaments = normaliser orbits on `N`-even
  graphs) with THM-4728's defect identity. **Hostile example:** the C15 counterexample (`K_(5,5)` with two fixed points
  joining the 3-orbit). **Corrected near miss:** HYP-9344 (2)'s `88` (it took the witness count `b_(k-2)(4)` for the
  defect and missed the tournament coincidences `E_T(k - 2)`; MISTAKE entry below if the audit confirms the reading).
  **Least-used sidecar:** THM-4723's remark that type counts need Sylow or inclusion-exclusion over disjoint markings —
  here the disjoint-orbit lemma makes the markings disjoint for free in the tournament half.

## 1. Concept board

| item | typed as |
|---|---|
| normaliser-orbit parametrisation `R = 2 t_p t_q b_k(4)`, product criterion for `N`-even | map (PROVED, THM-4871 (1)) |
| disjoint-orbit lemma (odd-order groups; prime-length cycles) | preserved predicate (PROVED, THM-4871 (2)) |
| `E_T = 0` for `k < p`, `E_T(p) = 4 t_p t_q` | PROVED (THM-4871 (3)) |
| `K_(q,q)` witness at `k = q - p`; infinite failure family `p < q < 2p` at `n = 2q` | PROVED (THM-4871 (4)) |
| census, HYP-9344 (2) refuted (`84`, not `88`) | FINITE-EXACT / REFUTED (THM-4871 (5)) |
| first-defect law, `b_(k-d)(4) - E_T^F(k-d)` formula, sign | NUMERICAL (HYP-9463 (1)-(3)) |
| first-order rules | OPEN with a quantifier-free reduction (HYP-9463 (4)) |
| fold / mirror / ball / `-7/4` | typed; small exact facts (section 5) |

## 2. The computation and what it proves

`04-computation/experiments/cpq_defect_20261104/cpq_defect_census.py p q kmax [--brute kb]` builds, for `H = <c_P c_Q>`
of type `(p, q, 1^k)`, one representative of every `N`-orbit of `H`-fixed tournaments (dilation classes of skew sets on
`Z_p` and `Z_q`, the cross orientation, a 4-coloured tournament on the fixed points up to isomorphism) and of every
`N`-even `N`-orbit of `H`-fixed graphs (multiplier-even dilation classes, cross, a colour-sign-even 4-coloured graph),
canonically labels each full object with nauty (pynauty 2.8.8.1) and reads the automorphism generators for the sign.
The `--brute` mode enumerates all `2^c` `H`-fixed objects and the normaliser by closure; it agrees with the orbit method
at `(3,5)` `k <= 3`, `(9,5)` `k <= 1`, `(3,25)` `k <= 1`. The `(3,5)` rows `k <= 4` reproduce THM-4728 / HYP-9344 exactly
(`R = 176, 1408`; `tau_T = 172, 1376`; `tau_E = 168, 1360`; `W = 4, 14`; `E_T = 4, 32`; `E_even = 4, 34`).

One subtlety worth recording: nauty's canonical form under a vertex colouring is canonical only within a fixed ordered
partition shape, so coloured types must be keyed by (colour-class sizes, canonical form); keying by the canonical form
alone silently merged non-isomorphic coloured tournaments (`82` instead of `b_3(4) = 88`) and produced wrong rows.

Running times on the 8 GB mac-mini: `(3,5)` through `k = 6` (`376,576` objects per side on 14 vertices) in 43 s;
`(5,7)` and `(3,7)` through `k = 6` (`753,152` objects on 18 and 16 vertices) in a few minutes each.

### 2.1 The census (THM-4871 statement 5)

Defects by family (`k = 0, 1, 2, ...`): `(3,5)`: `0, 0, 1, 4, 16, 84, 672`; `(3,7)`: `0, 0, 0, 0, 1, 4, 16`; `(5,7)`:
`0, 0, 1, 4, 16, 88, 704`; `(3,11)`, `(3,13)`, `(7,13)`: `0` through `k = 4, 3, 3`; `(9,5)`: `0, 0, 0, 0, 1`; `(9,7)`,
`(9,11)`: `0, 0, 1, 4`; `(9,13)`: `0, 0, 0, 0, 1`; `(3,25)`: `0` at `k = 0`. The full `tau_T, tau_E, W, E_T, E_even, E_wit`
table is in THM-4871 and in the `.out` files. In every row THM-4721 Theorem R (`R = RG`) and THM-4728's identity hold,
and `tau_T >= tau_E`.

### 2.2 HYP-9344 (2) is refuted

The pre-registered `defect(5) = 88 = b_3(4)` (against `64`) for C15 `(3,5,1^5)` at `n = 13` is wrong: the defect is `84`
(`tau_T = 17564`, `tau_E = 17480`, `W = 60`, `E_T = 260`, `E_even = 282`, `E_wit = 2`), and `defect(6) = 672` (`b_4(4) = 704`).
What HYP-9344 (2) got right is the *scale*: the defect is `b_(k-2)(4)` minus the tournament coincidences at level `k - 2`
(`88 - 4`, `704 - 32`), which is HYP-9463 (2).

### 2.3 The three proved pieces behind the law

* **Parametrisation and the product criterion** (THM-4871 (1)): the fixed points carry a 4-coloured object, the two
  orbits carry circulant classes, and `N`-evenness factors over the three blocks because an automorphism supported on
  one block can only reverse edges inside that block (vertex order `P < Q < F`).
* **The disjoint-orbit lemma** (THM-4871 (2)): in an odd-order permutation group the supports of two prime-length cycles
  are equal or disjoint (distinct primes: disjoint). Proof: the group they generate is transitive on the union, primitive
  by a block argument, and then Jordan (union size `>= r + 3`), parity (`r + 1`), or Galois + Feit-Thompson (`r + 2`,
  fixed-point set of a linear map has prime-power size, never `2`) rule out every overlap. Sanity: over all tournaments
  of orders 3-8, the 838 prime-length cycles in automorphism groups never overlap.
* **Consequences** (THM-4871 (3)-(4)): tournament coincidences need a second `p`-cycle living on the fixed points, so
  `E_T = 0` for `k < p` and `E_T(p) = 4 t_p t_q` exactly (the circulant constant-colour fixed-point tournaments, each type
  with exactly two `N`-orbits, never `Aut`-conjugate); `K_(q,q)` on the large orbit against the small orbit plus `q - p`
  fixed points is an `N`-even sign-odd `H`-fixed graph, so for primes `p < q < 2p` the defect at `k = q - p` is `>= 1`:
  HYP-9340 (a) fails for every such pair at degree `2q`, with `K_(q,q)` as the explicit obstruction (the C15 failure at
  `n = 10` is the first member).

### 2.4 The cancellation at `k = p` and the witness anatomy

At `k = 3` for `(3,q)`, `q >= 7`: `W = 2 t_q` (the `K_(3,3)` swaps of the 3-orbit with the three fixed points, `Q` joined to
both or neither), `E_even = 2 t_q` (the sign-even graphs with `Q` joined to exactly one of the two triples, two
non-conjugate `C_3q`'s), `E_T = 4 t_q` (the monochromatic cyclic triangles), and `2 + 2 - 4 = 0`. The lower bounds and
`E_T` are proved; that nothing else occurs is the census (`q = 7, 11, 13`). `witness_examples.py` prints the witness
configurations: at `k = q - p` the unique witness is `K_(q,q)`; at `(3,7)` `k = 4` the 17 witnesses are the one
`K_(7,7)` and sixteen `K_(3,3)`-type configurations (`2` cross x `2` circulant classes x `4` colours of the extra point),
and `17 + 48 - 64 = 1`.

## 3. The empirical laws (HYP-9463)

* First defect exactly at `k = |p - q|`, value `1`, then `4`, `16` (all eleven families).
* `defect(k) = b_(k-d)(4) - E_T^F(k - d)`: `(5,7)` gives the pure sequence `1, 4, 16, 88, 704` (`E_T(j) = 0` for `j < 5`),
  C15 gives `1, 4, 16, 88 - 4, 704 - 32`.
* `tau_T >= tau_E` in all 64 rows.

## 4. Uniform first-order rules

HYP-9463 (4) states the question exactly and records the one easy reduction: a quantifier-free rule built from the
atoms `S(c_i x)` is dilation-equivariant with Boolean degree bounded by the number of atoms, so the Codex degree theorem
excludes all of them uniformly in `p`; non-homogeneous atoms and quantifiers break equivariance, and then the
well-definedness of the induced map on dilation classes is itself a uniform first-order constraint. No obstruction for
the general case was found; the model-theoretic counting of definable sets in pseudofinite fields does not control sets
defined with an arbitrary predicate `S`.

## 5. Folds, mirrors, the ball, and `-7/4`

All statements here are elementary; the PDF's identities were checked at 40 digits (`scratch: mirror_ball.py`).

* **Mirror bases.** `f(x) = x^(1/x)` rises on `(1, e)` and falls on `(e, infty)`; the mirror map `m` with
  `f(m(x)) = f(x)`, `m(x) != x`, is a decreasing involution of `(1, infty) \ {e}`, `m(2) = 4` (the only rational pair,
  Gelfond-Schneider for the rest, session31), `r_n = m(n) = exp(-W_0(-ln n/n)) = e^(T(ln n/n))` with `T` the tree function,
  whence the Cayley coefficients `(m+1)^(m-1)/m!` (convergent for `ln n/n <= 1/e`, slowly at `n = 3`: 60 terms leave
  `0.09`). `m'(x) = f'(x)/f'(m(x))` with `f'(u) = f(u)(1 - ln u)/u^2`, so as `x -> 1+` (`m -> infty`)
  `m'(x) ~ -m^2/(ln m - 1) -> -infty`: at `x = 1.01`, `m = 658.8` and `m' = -76,722`. This is the owner's "small changes in
  the first are big changes in the latter", made quantitative: the steep side near `1` and the flat side near infinity
  are the two preimages of a fold.
* **The two constants.** `C = sum_(n >= 1) (-1)^n (n^(1/n) - 1) = 0.18785964246206712024851793405...` (the MRB constant)
  and `B = sum (-1)^n (r_n - 1) = 2.09679112439152612572503067508...` with `r_1 := 1`, `r_2 = 4`. `B = 3 - length(E)`,
  `E = union_(k >= 1) [r_(2k+2), r_(2k+1)]`, is the pairing of the alternating series (`r_n - 1` decreases to `0`);
  `C = sqrt 2 - 1 - int_E f'(u) du` is the mirror property `int_[r_(2k+2), r_(2k+1)] f' = f(2k+1) - f(2k+2)`; and
  `B - C = 4 - sqrt 2 + sum_(n >= 3) (-1)^n (r_n - n^(1/n)) = 1.90893148192945900547651274103...` converges absolutely
  because `r_n - n^(1/n) = e^(T(z)) - e^z = z^2 + O(z^3)`, `z = ln n/n` (the ratio to `(ln n)^2/n^2` is `1.018, 1.002,
  1.0003` at `n = 10^3, 10^4, 10^5`). All as in the PDF; nothing false was found in its mathematical notes. The
  irrationality of `C` and `B` is open, as the PDF says.
* **The unit ball.** `V(x) = pi^(x/2)/Gamma(x/2 + 1)`: `V(n) = 1, 2, 3.1416, 4.1888, 4.9348, 5.2638, 5.1677, 4.7248, ...`,
  integer maximum at `n = 5`, continuous maximiser `x* = 5.2569` (`V(x*) = 5.2778`), mirror dimensions `V(n') = V(n)`:
  `0 <-> 12.764, 1 <-> 10.814, 2 <-> 9.205, 3 <-> 7.820, 4 <-> 6.603, 5 <-> 5.517`. **No two distinct integer dimensions
  have the same ball volume** (PROVED): for `n = m (mod 2)` or `|n - m| >= 2` the ratio `V(n)/V(m)` is a nonzero integer
  power of `pi` times a rational, so equality would make `pi` algebraic (Lindemann); for `|n - m| = 1` the ratio
  `V(n)/V(n-1) = sqrt pi Gamma((n+1)/2)/Gamma(n/2 + 1)` is strictly decreasing in `n` and crosses `1` between `n = 5`
  (`16/15`) and `n = 6` (`5 pi/16`). **The alternating ball sum has a closed form** (PROVED):
  `sum_(n >= 0) (-1)^n V(n) = e^pi erfc(sqrt pi) = 0.28205917617568264518411963...`, from the generating function
  `sum V(n) t^n = e^(pi t^2)(1 + erf(sqrt pi t))` at `t = -1` (checked at 25 digits); the MRB analogue for the ball is
  elementary where the MRB constant itself is not.
* **`5` and `8 : 4`.** The twelve tournaments on five vertices split `8` self-converse `+ 4` not (two converse pairs),
  `7` asymmetric `+ 4` with `|Aut| = 3` `+ 1` with `|Aut| = 5` (the regular tournament), `6` strong `+ 6` not. The
  self-converse counts `1, 1, 2, 2, 8, 12, 88` match the terms of A002785 (self-complementary oriented graphs) as far as
  computed here; the ball-volume peak at dimension `5` and the `8 : 4` split share only the number `5`: NUMEROLOGY.
* **`-7/4` and `-29/16`.** Prior roles in the repo: the Thue closure `c in {0, -1, -2, -7/4}` of the mod-6 session (the
  dynamical Zsigmondy failures at `n = 3`: `c(c+1)^2 = -63/64` at `c = -7/4`, so `f^3(0) = c/2^6 = -7/256`), the airplane
  residue `51/49` (parabolic scaling, THM-4714), and THM-4139/THM-4146 (`x^2 - 29/16` is the unique centered monic
  quadratic over `Q` with a nondegenerate AP-supported rational 3-cycle, `-7/4 -> 5/4 -> -1/4 -> -7/4`, step `3/2`;
  no rational 6-cycle). New exact facts (gp): at `c = -7/4`, `(f^3(x) - x)/(f(x) - x) = (8x^3 + 4x^2 - 18x - 1)^2/64`, a
  perfect square (the period-3 saddle-node, i.e. the root of the airplane component), and **the parabolic 3-cycle
  `{-1.74698, -0.054958, 1.30194}` generates the real cyclotomic field `Q(zeta_7)^+`** (field discriminant `49`,
  cyclic cubic, `nfisisom` with `x^3 + x^2 - 2x - 1`): the `7` of `-7/4` is the conductor. At `c = -29/16` the rational
  AP 3-cycle has multiplier `35/8` and its companion 3-cycle (`64x^3 + 16x^2 - 164x + 23`) has multiplier `-23/8` and
  generates the cyclic cubic field of discriminant `31^2`. The fold reading: `c = -7/4` is where the attracting and the
  repelling 3-cycles of the real family merge — the dynamical "same value twice on either side", with the cyclotomic
  field at the fold.
* **Microcosm-macrocosm.** THM-4527's exact microcosm (the Syracuse map in the owner's labels) is a conjugacy, not a
  fold; the Galton-Watson duality behind the mirror bases (supercritical mean `ln n` conditioned on extinction becomes
  subcritical mean `ln r_n`, session31) is the fold with a probabilistic meaning. The two are typed as different
  mechanisms; no theorem connects them.

## 6. Sources typed

| source | status |
|---|---|
| *The MRB Multiverse* (PDF, fan fiction with mathematical notes) | notes verified (section 5); identical in substance to session31's packet |
| HYP-9459 (2), the Fibonacci law `u_3(n) = F(2n - 9)` | klein's; S75 lane; A001519 odd-indexed Fibonacci (OEIS 2026-11-04) |
| `{-3}`-self-duality of decomposable tournaments | literature pointer (Boudabbous-Lopez line), unverified, passed to S75 |
| `b_k(4) = 1, 4, 16, 88, 704, 8912, 188288, 6847200, 437069312` | not in OEIS (2026-11-04); `1, 4, 16, 84, 672` not in OEIS |
| THM-4721, THM-4728, HYP-9340, HYP-9344, Codex packets 2026-10-20 | inherited, CITED |
| Jordan 1873, Burnside (prime degree), Galois, Feit-Thompson, Jones 2014, Alspach-Parsons 1979, Palfy 1987, Muzychuk 1997 | CITED (classical) |

## 7. Reproduction

    cd 04-computation/experiments/cpq_defect_20261104
    python3 cpq_defect_census.py 3 5 6 --brute 3      # C15 rows k <= 6 (43 s); brute-force cross-check k <= 3
    python3 cpq_defect_census.py 5 7 6 ; python3 cpq_defect_census.py 3 7 6
    python3 cpq_defect_census.py 9 5 4 --brute 1 ; python3 cpq_defect_census.py 9 7 3 ; python3 cpq_defect_census.py 3 25 2 --brute 1
    python3 witness_examples.py                        # witness configurations at the first defects
    # requires: pip install pynauty (2.8.8.1)

## 8. What was not done

No proof of `defect = 0` below `|p - q|` or of the exact value `1` at `|p - q|`; no proof of the `b_(k-d)(4) - E_T^F` formula;
no C15 `k = 7` row (memory); no definability obstruction beyond the quantifier-free reduction; nothing on HYP-9172 (b)'s
decorated pairing (THM-4726 fixed what a carrier must carry; untouched here); the arc-connectivity and Fibonacci lane
(S75); `{-3}`-self-duality not verified.

## 9. Next obligations

1. Prove the first-defect law: classify the odd automorphisms of `N`-even `H`-fixed graphs for `k < |p - q|` through the
   modular decomposition (`P`, `Q` are modules; the only candidate odd elements are swaps of twin classes of odd size).
2. Prove the `k = p` cancellation in general (bijection between the `K_(p,p)` witnesses plus even-graph coincidences and
   the tournament coincidences), then the formula of HYP-9463 (2) by induction on `k`.
3. The C15 `k = 7` row with hashed canonical forms (predicted `8652`).
4. Model-check the two-atom first-order rules at `p = 13, 17` (HYP-9463 (4)).

## 10. Audit

Pending at this checkpoint; to be filled when the independent audit returns.

## 11. META-PATTERNS

Cards used: *disjoint markings before counting types* (THM-4723's lesson; here the lemma supplies the disjointness),
*compute the pre-registered test before theorising* (HYP-9344 (2) fell to a 2-second run once the orbit parametrisation
replaced brute force), *agreement of numbers is not agreement of mechanism* (`b_(k-2)(4)` was the right scale and the wrong
quantity). Candidate card: *when a census is "too large", look for the normaliser-orbit parametrisation* — the `2^24`
objects were `17,824` orbits.
