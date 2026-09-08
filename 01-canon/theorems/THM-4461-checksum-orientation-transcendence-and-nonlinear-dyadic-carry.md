---
id: THM-4461
title: "Checksum orientation transcendence and nonlinear dyadic carry"
status: >
  PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED. The THM-2225 checksum
  extractor has two transcendental orientation head-probability germs,
  given explicitly by the lacunary series L(t)=sum t^(2^r) at t=p(1-p).
  Its finite dyadic prefix probabilities admit a three-component nonlinear
  integer polynomial recurrence. Modulo two the composed lacunary germ
  collapses to p; the integral divided carry retains the lost coefficients.
  The scale index and polynomial degrees are unbounded, and the actual
  checksum decoder remains necessary. No stopping deadline is improved.
source: cross-concepts-20260908
depends_on:
  - THM-2225-dyadic-critical-run-extractors-and-cyclic-checksum-shell-bisection
related:
  - THM-2966-spine-normal-form-for-critical-run-fair-extractors
  - THM-3337-cross-shell-compression-attains-the-T4-floor
  - THM-3340-single-donor-cyclic-rotation-proves-all-pointwise-AMM-floors
  - THM-3342-sublinear-deadline-excess-is-impossible-for-fair-critical-run-extractors
  - THM-3343-shifted-donor-rotation-bisects-exactly-the-dyadic-annuli
  - THM-3344-orientation-splitting-saves-exactly-one-dyadic-donor-bit
  - THM-3577-amm-r512-offset-transition-and-causal-horizon
  - THM-4210-rule30-lossless-dyadic-block-current-cartier-tree
script: 04-computation/cross_concepts_donor_20260908.py
output: 05-knowledge/results/cross_concepts_donor_20260908.out
script_sha256: 107538d57fd064c5dbc2dd17ffe0b89118317d4ce9ffa4c950e6f8008edc362f
output_sha256: a47ea577eb359ad8425f3c86146937fa6495818600fb24c0bd7d768e6ad5e5ba
hash_basis: working-tree bytes; output LF
independent_audit:
  - cross-concepts-20260908-root
  - cross-concepts-20260908-tournament-repair
---

# THM-4461 -- Checksum orientation transcendence and nonlinear dyadic carry

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.** This recovers the older checksum construction as a concrete positive model for the scale-aware state left open by the September 6 rational-transport obstruction. It supplies exact orientation formulas, proves that both germs are transcendental, and exhibits an integer nonlinear dyadic recurrence for their finite prefixes. It changes no stopping deadline and does not solve the open uniform coin envelope, Mahler `3/2`, Rule 30, or LRC(14).

## Inheritance and portfolio

The inheritance baseline is `origin/main` at `47bf7e622`; namespace reservation was pushed at `8b25fa628`. The closest proved mechanism is **THM-2225**, `01-canon/theorems/THM-2225-dyadic-critical-run-extractors-and-cyclic-checksum-shell-bisection.md`: the checksum bisects each strictly interior tail-weight class, and the two constant-half words receive opposite verdicts. The least-used sidecar for the present question is precisely which orientation owns each of those two constant-half words. Their equal total mass hides a nonzero orientation transfer.

The inherited donor results are **THM-3340**, `.../THM-3340-single-donor-cyclic-rotation-proves-all-pointwise-AMM-floors.md`; **THM-3343**, `.../THM-3343-shifted-donor-rotation-bisects-exactly-the-dyadic-annuli.md`; and **THM-3344**, `.../THM-3344-orientation-splitting-saves-exactly-one-dyadic-donor-bit.md`. They separate pointwise attainment from one uniform rule and keep the donor capacity and orientation. **THM-3342**, `.../THM-3342-sublinear-deadline-excess-is-impossible-for-fair-critical-run-extractors.md`, forbids sublinear excess only after fixing an extractor; MISTAKE-368 prevents promoting this to a uniform gap above slope one.

The corrected near miss is compulsory shellwise closure: **THM-3337**, `.../THM-3337-cross-shell-compression-attains-the-T4-floor.md`, has a legal row-three spill past its old shell endpoint. THM-3344's one-bit sharpness keeps all other rows at their floors. Neither statement forbids a general cross-shell transport. **THM-3577**, `.../THM-3577-amm-r512-offset-transition-and-causal-horizon.md`, concerns a fixed policy, not alternative-prefix infeasibility.

The newer inherited obstruction is [September 6 rational orientation transport](../../05-knowledge/results/continuing8_20260906_coin_rational_handoff.md), independently audited there: for an actual exactly fair linear-deadline extractor, both orientation germs are nonrational. Eventual finite-dimensional **linear** polynomial row realizations, including periodic ones, are excluded. The canonical hostile is the ordinary von Neumann extractor, which has rational half-integral orientation germs but no finite deadline at fixed critical value two. Its missing causal row box matters.

The portfolio is: **Anchor**, recover an exact retained state for the live coin transport question; **Niche**, characterize the old checksum rule's orientation arithmetic; **Wildcard**, test what characteristic-two Cartier compression destroys. The five live concepts are:

| Concept | Representation and predicate | Operation / missing coordinate / cheapest test |
|---|---|---|
| Donor capacity | Actual integer Bernstein head-count rows | Compress by orientation sums; retain causal row decoder; replay literal checksum words |
| Rational obstruction | Two integer analytic orientation germs | Replace stationary row shift by dyadic scale; test an explicit nonrational model |
| Checksum shells | Cyclic tail response plus two constant-half words | Sum complete scales; retain which constant-half word is heads |
| Rule 30 Cartier state | Characteristic-two current sections | Test coefficient reduction on the exact coin state before proposing any dynamics map |
| Mahler reset clocks | Arithmetic native residual plus follower state | Compare what a reset erases; retain arithmetic carry rather than only finite reset labels |

No external literature claim or priority claim is needed. Targeted searches for coin orientation sums, lacunary series, dyadic series, and coin/Mahler transport found the September 6 obstruction but not the formulas below. That is a repository search statement, not novelty certification.

## 1. Recover the orientation that total fairness erases

Write `p=P(0)`, `q=1-p`, and let `F(p)` be the probability that the checksum extractor outputs heads and the initial bit is zero. Write `G(u)` for the corresponding initial-one heads probability in its native variable `u=P(1)`. Thus exact fairness is `F(p)+G(1-p)=1/2`.

For `m=2^r`, the shell consists of words of length `2m` whose first `m` bits are constant and whose whole word is not constant. At `m=1`, THM-2225 assigns `01` heads and `10` tails. Its contributions are

```
F_1(p)=pq,                     G_1(q)=0.
```

For every dyadic `m>=2`, the checksum splits every tail-weight class `1<=j<m` exactly in half, separately within each initial orientation. In the two remaining nonconstant words, `0^m1^m` is tails and `1^m0^m` is heads. Therefore the shell contributions are exactly

```
f_m(p) = (p^m/2)(1-p^m-q^m),
g_m(q) = (q^m/2)(1-q^m-p^m) + (pq)^m.                 (1)
```

These are probabilities of actual legal words, rather than a capacity relaxation. The half factors do not introduce half-integral ordinary coefficients: for dyadic `m>=2`, `1-p^m-(1-p)^m` is even coefficientwise.

Define the one-variable lacunary series

```
L(t)=sum_(r>=0) t^(2^r),       |t|<1,
L(t)=t+L(t^2).                                         (2)
```

Summing (1) over complete dyadic scales telescopes the pure powers. With `t=pq`, the exact limiting formulas are

```
F(p) = (3p-2p^2-L(p-p^2))/2,
G(u) = (2u^2-u+L(u-u^2))/2.                            (3)
```

For a finite prefix including shells `m=1,2,...,2^R`, define

```
L_R(t)=sum_(r=0)^R t^(2^r),       N=2^(R+1).
```

The exact formulas, including their undecided constant rays, are

```
F_[R](p) = (3p-2p^2-p^N-L_R(p-p^2))/2,
G_[R](u) = (2u^2-u-u^N+L_R(u-u^2))/2,
F_[R](p)+G_[R](1-p) = (1-p^N-q^N)/2.                 (4)
```

For `0<p<1`, the remaining mass tends to zero. Thus (3) are the actual orientation sums. Their germs at zero also agree coefficientwise with the limit because later shells have increasing minimum degree. The sidecar `(pq)^m` in (1), which disappears in the total head probability, is the source of `L(pq)`.

## 2. Both orientation germs are transcendental

This is stronger for this specific inherited extractor than the general September 6 nonrationality statement. It does not assert transcendence for every fair extractor.

**Lemma.** `L(t)` is transcendental over `Q(t)`.

**Proof.** Put `delta=1-t`, with `1/2<=t<1`, and `K=floor(log_2(1/delta))`. For the first `K+1` terms, Bernoulli's inequality gives

```
sum_(r=0)^K t^(2^r)
 >= K+1-delta(2^(K+1)-1) > K-1.
```

For the remaining terms use `(1-delta)^n <= 1/(1+n delta)`. Since `delta*2^K>1/2`, for `j>=1` the term at `r=K+j` is less than `2^(1-j)`. Consequently

```
K-1 < L(t) < K+3.                                     (5)
```

In particular `L(t)` tends to infinity, but `(1-t)^a L(t)^b` tends to zero for every positive integer `a` and every fixed nonnegative integer `b`.

Suppose a nonzero polynomial relation `sum_j A_j(t)L(t)^j=0` held. Let `v` be the least order of vanishing of any nonzero `A_j` at `t=1`, and let `j_0` be the largest index attaining that order. Divide the relation by `(1-t)^v L(t)^j_0`. Its `j_0` term tends to a nonzero constant. The terms of the same vanishing order and smaller index tend to zero because `L` tends to infinity; every term of larger vanishing order tends to zero by the preceding subpower bound. This is a contradiction. QED.

Now `p` is algebraic of degree two over `Q(t)` under `t=p-p^2`. If `L(p-p^2)` were algebraic over `Q(p)`, transitivity of algebraic extensions would make it algebraic over `Q(t)`. Substitution by `p-p^2` is injective on formal series, and has a formal inverse at zero, so this would make `L(t)` algebraic. Hence `L(p-p^2)`, and therefore both germs in (3), are transcendental.

No natural-boundary theorem, algebraic-function classification, or external analytic input is used. The growth estimate (5) is enough.

## 3. A fixed-dimensional nonlinear integer scale state nevertheless exists

The literal identity

```
L(p-p^2) = p                         in F_2[[p]]        (6)
```

follows from Frobenius: the summands become `p^(2^r)+p^(2^(r+1))` and telescope coefficientwise. This proves that the divided carry

```
C(p)=(L(p-p^2)-p)/2                         in Z[[p]]   (7)
```

is integral. In these coordinates (3) simplify to

```
F(p)=p-p^2-C(p),                 G(u)=u^2+C(u).         (8)
```

The first nontrivial terms are

```
C(p)=-p^3+p^4-2p^5+3p^6-2p^7+p^8+... .               (9)
```

There is an exact autonomous polynomial update over the integers at the **dyadic scale index**. Set

```
a_0=p,       B_0=0,       C_0=0,
a_(r+1)=a_r^2,
B_(r+1)=a_r^4-a_r^3+2a_r(1-a_r)B_r+2B_r^2,
C_(r+1)=C_r+B_(r+1).                                  (10)
```

Then, for every `r>=0`,

```
a_r=p^(2^r),
B_r=((p-p^2)^(2^r)-p^(2^r)+p^(2^(r+1)))/2,
C_r=(L_r(p-p^2)-p+p^(2^(r+1)))/2.                     (11)
```

**Proof.** The first formula is repeated squaring. If `a=a_r` and `B=B_r`, then `(p-p^2)^(2^r)=a-a^2+2B`. Squaring this identity and subtracting `a^2-a^4` gives exactly twice the displayed formula for `B_(r+1)`. Summing the `B_r` telescopes the pure powers and gives `C_r`. The initial identities are literal. QED.

The actual finite orientation probabilities are recovered by

```
F_[R](p)=p-p^2-C_R(p),
G_[R](u)=u^2+C_R(u)-u^(2^(R+1)).                      (12)
```

This is a fixed number of polynomial state components, with unbounded polynomial degree and an unbounded scale index. It is neither a finite set of states nor a stationary linear realization at every critical row. It therefore satisfies the September 6 obstruction rather than contradicting it. The squaring operation retains precisely what a finite matrix resolvent erased.

The causal decoder remains the actual checksum with its shell endpoint, initial orientation, and tail-weight response. Equations (10)--(12) summarize a known legal extractor; they do not turn arbitrary aggregate recurrences into legal row boxes. The non-dyadic hostile `m=6`, tail weight `j=2`, gives checksum counts `7/8`; merely running the same scalar formulas beyond their dyadic scope loses literal realization.

## 4. Connection contracts and what the new object changes

**Actual coin rows to scale state.** The source is THM-2225's causal checksum word coloring. The map first groups by dyadic shell and initial orientation, applies (1), and then applies (10). The preserved predicate is the exact finite heads probability in each orientation, including the undecided mass in (4). The lost information is which individual words and critical rows own that mass. The needed sidecar is the literal checksum decoder with its dyadic endpoint. The decisive test compares literal word probabilities with the independent binomial formula; it passes all 548 words in the four shells through endpoint 16.

**Characteristic-zero state to characteristic-two Cartier series.** The source is the exact lacunary series (2); the map is coefficient reduction modulo two and composition by `p-p^2`. It preserves the parity of every coefficient, and equation (6) is exact. It destroys every nontrivial coefficient of the transcendental germ `L(p-p^2)`: only `p` remains. The needed sidecar is the divided integral carry `C`, not a label saying that the parity state reset. The cheapest decisive witness is the coefficient `-2p^3` of `L(p-p^2)`, which becomes zero while `C` retains `-p^3`. The nonlinear recurrence (10) restores the entire missing coordinate.

This is a concrete control for **THM-4210**, `01-canon/theorems/THM-4210-rule30-lossless-dyadic-block-current-cartier-tree.md`, whose Rule 30 Cartier tree retains an infinite physical current sequence and a transverse quadratic channel. A small algebraic or finite-state parity observer does not imply a small characteristic-zero state; conversely the coin's transcendence does not imply anything about Rule 30's center sequence. No map from coin words to physical Rule 30 temporal evolution has been supplied, so no prize predicate is transported.

**Mahler reset comparison, with a recorded stopping reason.** The [September 6 reset-cost theorem](../../05-knowledge/results/continuing8_20260906_mahler_reset_cost.md) already shows actual launches `180` and `148` with the same length-eight follower/reset quotient and Rule 30 inverse monoid image, but different next reset events. The present source state is an integral coefficient carry, whereas that source is a native arithmetic residual. Both examples provide an exact failure of deleting a named coordinate. There is no identified intertwiner between the polynomial recurrence (10) and the Mahler `3/2` native-address update. This lane stops at that type mismatch; shared squaring, digits, or the word “Mahler” supplies no map. “Mahler functional equation” for (2) and the Mahler `3/2` Z-number problem are distinct uses of the name.

Revisiting the board: the donor lane gains a precise benchmark that a proposed scale-aware state should reproduce; the rational-obstruction lane gains an explicit transcendental survivor; the checksum lane identifies its otherwise invisible constant-half transfer; the Rule 30 lane gains an exact parity-loss hostile but no dynamics bridge; the Mahler lane retains its native residual and declines an unjustified recurrence transfer. No cosmetic tournament is introduced: the native structure here is a cyclic response and a polynomial scale action, not an intrinsic pairwise orientation.

## 5. Exact evidence and remaining frontier

Reproduce from the repository root:

```
python 04-computation/cross_concepts_donor_20260908.py
python -O 04-computation/cross_concepts_donor_20260908.py
```

The source imports no repository implementation and uses standard-library integer polynomials and exact fractions. The universe is all 548 shell words with `m=1,2,4,8`, all scale recurrences through `r=8` (degree 512), four exact interior biases per scale, the non-dyadic `m=6,j=2` hostile, and six rational near-one controls of (5). Independent paths compare literal word probabilities to binomial expansions, and nonlinear updates to closed polynomials. Finite fairness including the undecided mass is the printed consequence object. All **1,241 always-active gates** pass. The proof of transcendence and the all-scale identities are the mathematical arguments above, not finite extrapolations.

Both ordinary and optimized runs give the displayed transcript. Frozen working-tree SHA256 values are `107538d57fd064c5dbc2dd17ffe0b89118317d4ce9ffa4c950e6f8008edc362f` for the source and `a47ea577eb359ad8425f3c86146937fa6495818600fb24c0bd7d768e6ad5e5ba` for the LF output.

This result gives a reusable positive control for nonlinear scale transport. It supplies no new donor profile, no below-two deadline, and no feasible alternative to the fixed-policy failures in THM-3577. The next nontrivial target remains a capacity-preserving cross-shell decoder with an improved profile. It must retain enough orientation and causal information to realize actual row boxes; merely fitting (10) or another aggregate scale recurrence is insufficient.

## Independent promotion audit

The root cross-concepts agent and the tournament-repair agent independently accepted the full proof before promotion. Both checked the actual THM-2225 output conventions, the two orientation formulas, finite telescoping with the undecided constant rays, the logarithmic-growth transcendence proof and algebraic substitution, the integer nonlinear carry recurrence, and the causal decoder versus aggregate-state boundary. Neither found a required mathematical correction. The only proved dependency is the inherited checksum bisection in THM-2225; the transcendence and recurrence arguments are proved here. The September 6 rational obstruction and the other cross-domain results are related context, not proof dependencies.
