---
id: HYP-2059
status: PROVED (Pinch Lemma + r/s corollary; proof in reflection) + a clean provable target (N3)
source: opus-2026-06-02-S557
related:
  - HYP-2058
  - HYP-2055
  - THM-369
---

# HYP-2059: the LRC loneliness radius is r/s, pinned by a straddling pair

## Proved (the exact moments & conditions)

`f_S(t)=min_i||v_i t||`, `M(S)=max_t f_S(t)`. LRC(14) ⟺ `M(S)≥1/14` ∀S.

- **Pinch Lemma.** If `M(S)<½`, the max is attained at `t*` where two distinct
  runners `a,b` **straddle**: `frac(v_a t*)=M`, `frac(v_b t*)=1−M`. (Proof: the
  max of a min-of-tents is a breakpoint with a positive-slope binding runner on
  the left and a negative-slope one on the right; equal only at an apex ⇒ M=½.)
- **Radius = r/s.** `(v_a+v_b)t*≡0` ⇒ **`t*=m/(v_a+v_b)`**, `(v_a−v_b)t*≡2M`,
  and **`M(S)=r/s`** with `s=(v_a+v_b)/gcd(v_a,v_b)` (reduced pair-sum), `r≥1`.
  Verified exactly: 200/200 random 13-sets + every critical config.

## Exact reformulation

**LRC(14) ⟺ every 13-set has a pair `(a,b)` and integer `m` with
`‖v_j·m/(v_a+v_b)‖ ≥ 1/14` for all `j`.** The only times that matter are
`m/(v_a+v_b)`.

- **N1 (counterexample ⇒)** the optimal binding pair has reduced sum `s≥15`.
- **N2 (tight ⇒)** binding reduced sum `s≡0 (mod 14)`; at the floor, `s=14`
  (the two straddling runners satisfy `(v_a+v_b)/gcd=14`). This deductively
  yields the tight-witness lattice and subsumes: the sieve's "multiple of 14"
  (pair `(1,13)`), the spectral-gap value `2/27` (apex-doubled pair, `s=27,r=2`).
- **N3 (the provable target)** every 13-set has a pair with reduced sum `≤14`
  whose pinch time clears all other runners. Abundant candidate pairs
  (`(1,13),(2,12)[s=7],…,(6,8)`); content = one clears the rest.

## Next step

Attack N3 by pigeonhole over the `⌊n/2⌋` "small" pairs summing to `n`: show some
pinch time `m/n` avoids every other runner's (thin) danger arcs. This routes the
even-fold/odd-coupling residual (S554) through *which pair binds*.

**See:** `07-reflections/lrc-n14-the-exact-moments-pinch-pair-and-r-over-s-radius-s557.md`,
`04-computation/lrc_n14_pinch_pair_radius_s557.py` (+.out); HYP-2058, HYP-2055,
oracle-S552 (gap), THM-369.
