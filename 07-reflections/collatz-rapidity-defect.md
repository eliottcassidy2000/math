# The Collatz harmonic defect: what the "+1" costs in rapidity

*Session: collatz-defect-rapidity. Builds on `collatz_rapidity_s116n.py` and the
"arctanh is the universal linearizer" thesis (README, grand-trichotomy.md).*

## The one-line result

The Collatz map, viewed in the formal-group logarithm `rho = (1/2)ln = arctanh`
(rapidity = `log_F` for `F(x,y)=(x+y)/(1+xy)`), obeys an **exact conservation
law**. Earlier work (`collatz_rapidity_s116n.py`) saw the rapidity walk with
drift `ln(3/4)/2` "plus corrections." Those corrections are not slop — they are
a single, sign-definite, **telescoping harmonic defect**:

```
    ln n = K·ln2 − L·ln3 − D(n),     D(n) = Σ_{i=0}^{L-1} ln(1 + 1/(3 a_i)) > 0
```

for the Syracuse trajectory `n = a_0 → a_1 → … → a_L = 1`, with `L` odd steps,
`k_i = v_2(3a_{i-1}+1)`, `K = Σ k_i`. Proof is two lines: each step is
`2^{k_i} a_i = 3 a_{i-1} + 1`; take logs and telescope.

## Why this is the *right* coordinate

This is the repo's central thesis applied to Collatz. `arctanh = log_F`
linearizes the formal group: multiplication-by-3 and halving become **pure
rapidity translations** (`+½ln3`, `−½ln2`). A clean formal-group map would be a
straight additive walk. Collatz is *not* clean — and `D(n)` measures, exactly,
the entire non-formal-group content. The only non-algebraic ingredient of
Collatz is the `+1`, and in rapidity the `+1` is precisely `−D(n)/2`.

## The arithmetic shadow

Exponentiating gives an **exact rational identity** (verified with `Fraction`):

```
    n · 3^L · Π_i (3 a_i + 1)/(3 a_i) = 2^K
```

The pure map `x → 3x` would demand `2^K = 3^L · n`, which is impossible for
`n > 1` because `3` never divides a power of `2`. The fudge factor
`Π (3a_i+1)/(3a_i) = e^{D(n)} > 1` is *exactly* the correction that repairs the
divisibility. So the `+1` is not arbitrary: it is the minimal arithmetic
perturbation that lets the multiplicative `2^K / 3^L` flow land on an integer.

## Two consequences

**(1) Forced inequality (unconditional, per-orbit).** Since `D(n) > 0` and
`ln n ≥ 0`,
```
    K/L = log2(3) + (ln n + D(n))/(L·ln2)  >  log2(3) = 1.5849625…
```
Every `n > 1` that reaches 1 has mean halving exponent *strictly above*
`log2(3)`, and `K − L·log2(3) = log2(n) + D(n)/ln2` exactly. This is not an
"on average" statement — it holds for each individual trajectory. (The only
conjectural input is that `L, K` are finite, i.e. Collatz itself.)

**(2) The defect is bounded (CONJ-defect-bound, HYP-2148).** `ln(1+1/(3a)) ~
1/(3a)`, so `D(n)` is dominated by the handful of *small* odd values in the
terminal cascade to 1. Empirically `D(n) < D* ≈ 0.2257` for all `n`; the
champion `n = 993` (`D = 0.225654`) is unbeaten over all odd `n` up to
5,000,000, and the record progression converges rather than growing. The
heuristic: the Syracuse predecessor tree of 1 is `{1,5,21,85,341,…} =
(4^j−1)/3`; a single trajectory threads one root-path and can only collect a
bounded budget of small odds. Collecting *all* small odds (which would diverge
like `(1/6)ln M`) is forbidden by the dynamics.

## Where it connects in the repo

- **Formal group / rapidity:** sharpens `collatz_rapidity_s116n.py` and the
  `rapidity_numbertheory_s116b.py` thread; `D(n)` is the exact obstruction term
  the "universal linearizer" framing was missing for Collatz.
- **Grand trichotomy** (`07-reflections/grand-trichotomy.md`, Collatz row):
  halving = INERT (prime 2, removes bits), tripling = RAMIFIED (prime 3,
  perturbation), drift = SPLIT. `D(n)` quantifies the RAMIFIED part precisely.
- **VS chain analogy** (already in `collatz_rapidity_s116n.py` §XII): both are
  iterated arithmetic maps with 2-adic-controlled attractors; the rapidity
  conservation law gives Collatz the kind of exact bookkeeping the VS chain has.

## Open threads / next steps

1. **Prove CONJ-defect-bound.** Reduce `sup D` to a feasibility statement over
   the small-odd predecessor tree of 1: which finite multisets of small odds can
   co-occur on one trajectory? A bound there gives `D* < ∞` unconditionally.
2. **Closed form for `D*`?** `0.225654…` — is it the value of an explicit sum
   over the optimal small-odd path, or just a sup with no closed form? (`ln(5/4)
   = 0.2231` is close but not it.)
3. **OEIS hooks.** `K+L` = total stopping time (A006577); `L` relates to
   A006667; the exact defect numerator `Π(3a_i+1)` and denominator `Π(3a_i)`
   are candidate new sequences (the "rapidity-defect" pair). Worth extending and
   submitting per the repo's OEIS practice.
