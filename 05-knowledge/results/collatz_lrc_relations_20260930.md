# Lonely runner relations and the Collatz cycle equation: Beck–Everett's harmful relations sharpen the repo's LRC(14) relation lane (square norm `195 → 65`, 1-norm `50 → 29`, and the odd coordinate sum proved for the short relation), the window-product Fourier argument and the Collatz character sum have the same shape, the cycle half of Collatz is of lonely-runner type (a finite check at each period) while the divergence half is not, Dirichlet's theorem is the common root, and the Kawasaki fixed-point "proof" is the repo's THM-4471

**Session:** opus, `collatz-posets-zeta5-20260927` (S15, sixteenth note), 2026-09-30.
**Owner's directive:** "both of the major open problems this repo has focused
on are finely tuned to be aligned with the major subtle patterns in the
primes. lonely runner being local while collatz being global means that each
lacks information about the other's type; the former is hard to pin down
whether all possible speed sets behave this way, the latter is difficult to
see if no large cycles form. Understanding both at the same time may unlock
the number theoretical insights to each", with two attached papers:
M. Beck, S. Everett, *Lonely runner relations*, arXiv:2609.06259v1 (5 Sep
2026), and T. Kawasaki, *A proof of the Collatz conjecture*, arXiv:2502.20642v1
(math.GM, 28 Feb 2025).

**Status: CITED (Beck–Everett's Theorems 1, 2, 5, Proposition 4, Corollary 6,
read in full; their Lemma 3 identities re-verified numerically); the import
into the LRC(14) lane is a CITED sharpening of THM-4009 with FINITE-EXACT
support-two counts; FINITE-EXACT checks of Theorem 2 on all tight instances
in small boxes (`k = 2..5`); the per-period bound for Collatz cycles PROVED
(elementary) with a FINITE-EXACT table; the Kawasaki paper is REFUTED by
THM-4471 (cited, not re-derived); the joint reading typed ANALOGY where it
is one; LRC(14) and Collatz OPEN. Independent audit OWED.** Script
`04-computation/experiments/collatz_lrc_relations_20260930.py` (output `.out`).

## 1. Beck–Everett, digested

For `n ∈ Z^k_{>0}` let `lon(n) = max_t min_j ‖t n_j‖`; LRC says `lon(n) ≥ 1/(k+1)`
for distinct entries (Dirichlet's theorem makes `n = (1, ..., k)` tight). A
*relation* is `m ∈ Z^k` with `m · n = 0`; it is *harmful* if `m_1 + ... + m_k` is
odd.

* **Theorem 1.** If `lon(n) ≤ 1/(k+1)` (a counterexample or a tight instance)
  then (a) `n` has a harmful relation with `‖m‖_1 ≤ 2k + 3`, and (b) a relation
  with `‖m‖_1 ≤ (k+1)/(k-1) · flt(k)`, `flt(k)` Khinchin's flatness constant.
  So all counterexamples and tight instances lie on finitely many
  hyperplanes `m · x = 0`.
* **Theorem 2 (the Fourier bound).** If `n` has no strict lonely time
  (`‖t n_j‖ > 1/(k+1)` for all `j` fails for every `t`), then there is a harmful
  relation with `‖m‖_2 ≤ 2(k+1)/√(k-1)` and `‖m‖_1 ≤ 2k + 3`.
  *Mechanism.* With `h = (k-1)/(2(k+1))`, a window `W_h` (the autocorrelation
  of a cosine bump of width `h`, shifted by `1/2`) is positive exactly when
  `‖x - 1/2‖ < h`, i.e. when the runner is in the lonely arc, and has Fourier
  series `Σ_r (-1)^r a_r e(rx)` with `a_r ≥ 0`, `Σ a_r = 1`, `Σ r^2 a_r = 1/(4h^2)`
  (both identities re-verified numerically here). The product
  `H(t) = Π_j W_h(n_j t)` is positive exactly at strict lonely times, so under
  the hypothesis `H ≡ 0` and every Fourier coefficient of `H` vanishes: for
  each frequency `q` the mass of even-sum frequency vectors `c` with `c · n = q`
  equals the mass of odd-sum ones. A second-moment identity,
  `Σ_c w(c) ‖P c‖_2^2 = (k-1)/(4h^2)` for the projection `P` onto `n^⊥`, plus
  averaging, produces one even and one odd vector in a common fibre whose
  difference `m` is a harmful relation in `n^⊥` with `‖m‖_2^2 ≤ (k-1)/h^2`.
* **Proposition 4 / Theorem 5 (the zonotope).** `lon(n, s) ≥ λ` iff the
  zonotope `Z_{λ,γ}(n, s) = s + [λ, 1-λ]^k + [0, γ] n` contains a lattice point;
  equality iff it is hollow. A hollow zonotope has lattice width
  `≤ flt(k)`, and its width is `min_{m ≠ 0} (1 - 2λ)‖m‖_1 + γ |m · n|`, so a
  relation with `‖m‖_1 ≤ flt(k)/(1 - 2λ)` exists. This gives Theorem 1(b) and,
  for random speeds, Czerwiński's theorem: with probability `→ 1` the
  loneliness `1/(k+1)` can be replaced by `1/2 - ε` (Corollary 6). The shifted
  LRC (arbitrary starting points) is false for `k ≥ 5` (`n = (1,2,3,4,5)`).
* **Open problems named by the authors:** `inf lon(n, s)` in small dimensions,
  the statistics of `lon` on a hyperplane `m · x = 0`, re-deriving the lacunary
  cases, flatness constants of lonely-runner zonotopes (the known tight
  instances have width `< 3`), and a possible connection to Schmidt's subspace
  theorem (points of small height lie on finitely many hyperplanes).

**Checks (FINITE-EXACT).** For `k = 2..5` and speeds `≤ 24, 18, 14, 11`, every
tight instance (`12, 6, 5, 3` of them: the dilations of `{1..k}` and, for
`k = 4, 5`, the sets `{1,3,4,7}`, `{2,6,8,14}`, `{1,3,4,5,9}`) has a harmful
relation of 1-norm `3` (`2n_1 = n_2` or `n_1 + n_2 = n_3`), inside both bounds;
no counterexample appears. The converse fails as it must: among the `k = 3`
sets with speeds `≤ 12`, `196` have a harmful relation of 1-norm `≤ 9` and `192`
of them are strictly lonely.

## 2. The import into the repo's LRC(14) relation lane (CITED sharpening)

LRC(14) is `k = 13`. The repo's necessary conditions on a primitive
counterexample are THM-3743 (a relation with `‖a‖_1 ≤ 356`, superseded) and
THM-4009 (a Graver relation with `‖a‖_2 < 14`, `Σ a_i^2 ≤ 195`, `‖a‖_1 ≤ 50`,
`|a_i| ≤ 13`; "the projected half-lattice centre also forces some odd-sum
relation, not proved to be the short one"). Beck–Everett's Theorem 2 gives,
for any counterexample *or tight instance*:

| bound | THM-4009 (2026-08-24) | Beck–Everett (2026-09-05) |
|---|---|---|
| `Σ m_i^2` | `≤ 195` (`‖a‖_2 < 14`) | `≤ 65` (`‖m‖_2 ≤ 2·14/√12 = 8.083`) |
| `‖m‖_1` | `≤ 50` | `≤ 29` |
| odd coordinate sum | some odd-sum relation exists, not proved short | **the short relation is harmful** |
| support-two branch (coprime ratios `a : b`) | `47` with `a^2 + b^2 < 196` | `11` with `a^2 + b^2 ≤ 65` and `a - b` odd: `1:2, 1:4, 1:6, 1:8, 2:3, 2:5, 2:7, 3:4, 4:5, 4:7, 5:6` |

Caveats for the LRC sessions: Beck–Everett's relation is any integer
relation, not necessarily a Graver element; both statements are necessary
conditions and may be imposed simultaneously; THM-972's relation lock
(weight `≤ 14`) does not reach `29`. The relation-ledger coordinate `R` of
the lens map should carry the harmful short relation as a sidecar
(Direction D53). The repo's THM-3743 already used Khinchin's flatness on the
lonely-runner polyhedron, so Theorem 1(b) is a cousin; Theorem 2's Fourier
route and the parity are the new ingredients.

## 3. The joint reading (the owner's thesis, typed)

| | LRC | Collatz cycles | Collatz divergence |
|---|---|---|---|
| statement | `∀` finite speed sets `∃` a time | `∀` words `w ≠ (2)^k`, `x_w ∉ Z_{>0}` | `∀ n`, the orbit descends |
| local parameter | `k` (number of runners) | the period `p` (or the shape `(A, p)`) | none |
| finite check at fixed parameter | Tao: counterexamples have bounded speeds, so LRC(`k`) is decidable | every cycle element of shape `(A, p)` is `≤ S_max(A,p)/(2^A - 3^p)` (table in the script: `≤ 381` for `p ≤ 12` at the least positive clock), so cycles of period `p` are a finite check | none known (Conway: the family is Turing-complete) |
| settled range | `≤ 13` runners | all shapes `A ≤ 22` (eleventh note); no cycles with `≤ 91` runs (Steiner; Simons–de Weger; Hercher, Baker's method) | almost all orbits (Tao, logarithmic density) |
| the exceptional object satisfies a relation | `m · n = 0`, `‖m‖_1 ≤ 2k+3`, odd sum (Beck–Everett) | `Σ_i 3^{p-1-i} 2^{d_i} = x (2^A - 3^p)`: an `S`-unit relation with one unknown integer (THM-4490's affine-word specialization) | — |
| parity layer | harmful = odd coordinate sum | the clock is `(-1)^A mod 6`, the carry is odd (twelfth note) | — |
| Fourier shape | `H ≡ 0 ⟹` all coefficients vanish `⟹` even/odd masses balance in every fibre | no cycle of shape `(A,p) ⟺ (1/D) Σ_{r ≠ 0} Σ_w e(r S_w/D) = -C(A-1, p-1)` (thirteenth note §3.4; S19's resonances) | — |
| Dirichlet root | the constant `1/(k+1)` *is* Dirichlet's theorem; `{1..k}` saturates it | cycle shapes are Dirichlet approximations of `log_2 3` (the convergents; Baker bounds the saturation) | — |
| primes | prime-13 tightness, modulus supply (THM-979); harmful parity is 2-adic | the clock primes form a lattice (eleventh note, Prop 2); `2` and `3` | — |

So the owner's dichotomy refines: **the cycle half of Collatz is of lonely
runner type** (a universal statement over a discrete parameter with a finite
check at each value and an unbounded family of values, where the exceptional
objects satisfy short relations and the Fourier side turns non-existence into
exact cancellation), and **the divergence half is the genuinely global part**,
with no parameter, no finite check, and no relation to satisfy. "Large cycles"
in Collatz and "all speed sets" in LRC are the same worry: the per-parameter
checks are feasible only for small parameters (`k ≤ 13`; few runs), and no
uniform argument exists on either side. What each side can lend the other
today is one-directional and typed ANALOGY: Beck–Everett's window-product
argument needs a *product* structure (the Fourier expansion of `Π_j W(n_j t)`
factorizes over runners), while the Collatz word sum `Σ_w e(r S_w/D)` is a
transfer-matrix product over levels, not a product over independent
factors — the obstacle D54 names.

## 4. The Kawasaki paper (REFUTED, THM-4471)

The repo already refuted arXiv:2502.20642: its fixed-point theorem
(Theorems 2.1(5)–2.3(5), weighted generalized pseudocontractions) is false
— the successor map on `(N, |x - y|)` satisfies every hypothesis with the
paper's constants and has no fixed point; the proof needs an alternative at
the swapped pair that condition (5) does not provide, and the claimed decay
fails at `x = 5`; the paper's coefficient table equally "proves" that `3n - 1`
reaches `1`, whose orbit `5 → 7 → 10 → 5` does not (SHEET control); and no
contraction-type argument in `|x - y|` can prove Collatz, since every odd step
stretches consecutive distances by at least `4/3`. The surviving metric
principles (Caristi; a uniformly discrete contraction metric; Bessaga) are
reformulations of the conjecture. Nothing in the paper touches the joint
structure above.

## 5. Verdicts

| claim | status |
|---|---|
| Beck–Everett Theorems 1, 2, 5; Lemma 3 identities | CITED (paper read in full) + numerically re-verified |
| tight instances for `k ≤ 5` in small boxes satisfy Theorem 2; converse fails | FINITE-EXACT |
| LRC(14): harmful relation with `Σ m_i^2 ≤ 65`, `‖m‖_1 ≤ 29`; support-two ratios `47 → 11` | CITED sharpening of THM-4009 (necessary condition, not a proof) |
| per-period finite check for Collatz cycles; table `p ≤ 12` | PROVED (elementary) + FINITE-EXACT |
| the joint table; cycle half of LRC type, divergence half global | typed ANALOGY / EXACT where marked |
| Kawasaki | REFUTED (THM-4471) |
| LRC(14), Collatz | OPEN |

## 6. Directions

* **D53 (for the LRC sessions).** Enter Beck–Everett's harmful relation into
  the relation ledger: recompute the support-two and support-three branches
  of THM-4009 under `Σ m_i^2 ≤ 65` with odd sum, and test whether the shortest
  harmful relation can be chosen Graver-minimal; check whether THM-972's lock
  extends from weight `14` to the relevant weights.
* **D54.** A product structure for the Collatz word sum: the full-period
  recursion of S23 and the Gauss-twisted circulant of THM-4520 write the
  level sums as products of operators; if the character sum modulo a clock
  can be expressed as a product over levels, Beck–Everett's second-moment
  averaging may transfer.
* **D55.** The tight instances of LRC beyond dilations (`{1,3,4,7}`,
  `{1,3,4,5,9}`; Goddyn–Wong's list, CITED from memory) against the Collatz
  "tight" shapes (the convergents of `log_2 3`): both are the saturating
  cases of Dirichlet's theorem; is there a common parametrization?
* **D56.** Schmidt's subspace theorem as the shared Diophantine root
  (Beck–Everett's last open problem; the cycle equation as an `S`-unit
  relation): what the subspace theorem says about the hyperplanes containing
  the Collatz cycle vectors `(2^{d_0}, ..., 2^{d_{p-1}}, x 2^A, x 3^p)`.
