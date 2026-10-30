# Two aliquot residues are a reversible p-adic coordinate system

**Status: PROVED**, with independent finite divisor-sum and residue-permutation checks in `verify_aliquot_modular.py`. This refines the inherited exact reconstruction and retained-core dynamics with a reversible finite-residue coordinate system. It does not provide a fast method for evaluating divisor sums or factoring an integer.

## 1. Inherited exact inverse and the apparent missing coordinate

Write `s(n)=sigma(n)−n`. Fix a prime `p` and an integer `t≥1`. Put

    S_j = 1+p+...+p^(j−1),    j≥1.

The prime-dilation inverse and retained-core affine doubling update already appear in [Session 78](../../session78_fourth_power_aliquot_20261027/ALIQUOT.md); [Session 81](../../session81_counterexample_families_20261028/ALIQUOT_OBSERVERS.md) extends reconstruction to prime-power dilations. [Session 110](../../session110_overnight_20261030/CHECKPOINT_01.md), Section 1, gives the recent presentation used here: if `a=s(n)`, `b=s(p^t n)`, and `n=p^e u` with `p∤u`, then

    b−p^t a = S_t sigma(u),
    v_p(b−a) = e,
    n = S_(e+1) (b−p^t a)/S_t − a.                 (1)

Its proof uses `b−a=p^e S_t [u+p s(u)]`, with bracket and `S_t` units modulo `p`. At first sight, reducing the observations modulo `p^k` loses the full exponent `e`, so (1) might seem unusable when `e≥k`. The correction is that `S_(e+1)` has already stabilized modulo `p^k` once `e≥k−1`. The missing valuation is then irrelevant to the target residue.

The retained-core idea itself is inherited: Session 78 already uses `(s(n),sigma(oddpart(n))) -> (2s(n)+sigma(oddpart(n)),sigma(oddpart(n)))` on a doubling ray. Our refinement below is the universal residue-plane decoder, its Latin/isometry structure, and the finite conjugacy and precision boundaries. The earlier modular-trace collisions `27,35 mod 8` and `3,7 mod 4` are compatible: each pair has the same source residue. Recovering a source residue does not distinguish all unbounded integer sources.

## 2. Universal residue decoder

Fix `k≥1`, and let `A,B` be **any** residues modulo `p^k`. Set `d=B−A` and

    h = min(v_p(d), k−1),

where a zero residue has truncated valuation at least `k`. Define

    F_(p,t,k)(A,B)
      = S_(h+1) S_t^(-1) (B−p^t A) − A  (mod p^k).       (2)

The inverse of `S_t` exists because `S_t≡1 mod p`.

**Theorem A.** For every positive integer `n`, without a size or coprimality restriction,

    F_(p,t,k)(s(n),s(p^t n)) = n (mod p^k).

**Proof.** If `e<k−1`, the residue of `b−a` determines `e` exactly. If `e≥k−1`, then `S_(e+1)≡S_k mod p^k`, because all omitted summands are divisible by `p^k`. In both cases (2) is (1) reduced modulo `p^k`.

Consequently, if `1≤n≤B<p^k`, the two aliquot residues recover the **integer** `n` directly. One does not need enough CRT modulus to recover the potentially much larger values `s(n)` or `s(p^t n)` first. If the range is `1≤n≤p^k`, the residue zero uniquely denotes its upper endpoint; a range with two representatives of a residue needs additional information.

This is a native prime-power observer. The inverse denominator `S_t` can contain primes other than `p`, but none is bad for the `p`-primary observation. For example `p=2,t=2` gives `S_t=3`; it is invertible modulo every power of 2. Trying the same division modulo 3 fails.

## 3. Stronger result: the decoder is a permutation in its second observation

Rewrite (2), using `(p−1)S_(h+1)+1=p^(h+1)`, as

    F(A,B) = U_h d − p^(h+1) A,
    U_h = S_(h+1) / S_t,      d=B−A.                    (3)

Every `U_h` is a unit modulo `p^k`.

**Theorem B.** For each fixed residue `A`, the map `B↦F(A,B)` is a permutation of `Z/p^k Z` that preserves p-adic distances. Thus

    (A,B) ↔ (A,F(A,B))

is a bijection of the full `p^(2k)`-element residue plane, not merely an injection on pairs that happen to arise from integers.

**Proof.** On a nonzero difference shell `v_p(d)=h<k`, (3) is multiplication by a unit followed by translation by a multiple of `p^(h+1)`. It permutes that shell. For `h=k−1`, the translation vanishes modulo `p^k`; the conclusion remains true. The zero difference maps to zero. Distances between two differences in the same shell are preserved by the common unit multiplier. Differences in different shells have distance governed by the smaller valuation; their images retain those valuations, so that distance is preserved too.

The explicit inverse is particularly simple. Given `A` and desired source residue `r`, take `h=min(v_p(r),k−1)` and set

    B = A + (S_t/S_(h+1)) [r+p^(h+1) A]  (mod p^k).    (4)

The denominators in (4) are units. Formula (3) shows why the valuation of `r` equals that of `d`, so the same `h` is used in both directions. Each source residue therefore has exactly `p^k` formal observation pairs `(A,B)`, one for each `A`. This is not a claim that every formal pair is realized by an actual integer's aliquot sums.

### Both coordinates are reversible, not just the second

For fixed `B`, the map `A↦F(A,B)` is also an isometry and permutation. To see this, write `d=B−A`. Changing `A` by `epsilon` changes `d` by `−epsilon`. The vertical isometry gives an output change of norm `|epsilon|_p`; the remaining change in the explicit `A` coefficient has norm at most `p^(−1)|epsilon|_p`. The ultrametric inequality makes the larger term decisive, so the exact output distance remains `|epsilon|_p`.

Equivalently, given `B` and source residue `r`, with `h=min(v_p(r),k−1)`, the missing first observation is

    A = [S_(h+1) B − S_t r] / S_(h+t+1)  (mod p^k).     (5)

The denominator is again a unit. Formula (5) follows by rearranging (4) and using `S_(h+t+1)=S_(h+1)+p^(h+1)S_t`.

Thus `F` is a **Latin-square operation** on the residue ring: either observation is uniquely determined by the other and the decoded source residue. No associativity or group law is asserted.

The two-coordinate change `(A,B)↦(A,F(A,B))` is also an isometry for the sup p-adic metric. Pass first to `(A,d)=(A,B−A)`, a linear sup-metric isometry. The map `(A,d)↦(A,U_h d−p^(h+1)A)` is 1-Lipschitz because it is isometric in `d` and strictly contracting in `A` when `d` is fixed. Its explicit inverse has the same two properties. Mutual 1-Lipschitz bounds give equality of distances. In the inverse limit it preserves product Haar measure, since it permutes every finite residue plane.

**Sharp metric boundary.** For a fixed first observation `A`, a discrepancy of exact order `p^j` in the second observation changes the decoded source by exact order `p^j`. No extra p-adic precision can be inferred from a noisy second reading merely from this coordinate map. On genuine source pairs, a discrepancy in `A` is a separate input and must also be retained.

## 4. Compatible inverse limits and a finite arithmetic automaton

The decoders commute with reduction from modulus `p^(k+1)` to `p^k`. The only potentially changed valuation case is when the difference vanishes at the smaller precision; then its geometric sum is already stable there. Hence these finite bijections define a continuous coordinate change on `Z_p²`. For fixed `A`, it is an isometry in the second coordinate. This statement concerns the decoder; it does **not** assert that `n↦sigma(n)` is p-adically continuous or that the observation pair depends only on the source residue. For example, at `p=2,t=1,k=2`, sources `1` and `5` are congruent modulo 4, but their observation pairs are `(s(1),s(2))=(0,1)` and `(s(5),s(10))=(1,8)≡(1,0)`. Both decode to residue 1; the retained first coordinate distinguishes them.

There is also a genuine finite transition rule for aliquot observations along a prime-multiplication ray. Keep the state

    (r,A) = (n mod p^k, s(n) mod p^k).

Compute `B` by (4), then update to

    (p^t r mod p^k, B).

For actual integers this is exactly the state of `p^t n`; it extends as a defined transition on every formal pair. After at most `ceil(k/t)` steps, the first coordinate is zero. Formula (4) then gives `B=A`, so every orbit has reached a fixed state. In particular the quotient has no nontrivial cycles.

For an actual source `n=p^e u`, its eventual fixed second coordinate is

    sigma(u)/(1−p) mod p^k.

Indeed `s(p^(e+jt)u)=S_(e+jt+1)sigma(u)−p^(e+jt)u`, and once `e+jt≥k` the geometric sum is `1/(1−p)` modulo `p^k` and the second term vanishes. This fixed coordinate retains arithmetic information about the p-free core even after the visible source residue has become zero.

This is a useful distinction between **reversible observation coordinates** (Theorem B) and **irreversible dynamics** (multiplication by `p` modulo `p^k`). The two aliquot readings are not independent extra bits: they are exactly another coordinate system for the source residue together with its first aliquot residue.

### Exact conjugacy: the nonlinear update is scaling with a retained core

The finite dynamics has a stronger normal form. For any formal state `(r,A)`, define

    h=min(v_p(r),k−1),
    C=(A+r)/(1−p^(h+1))  (mod p^k).                    (6)

The denominator is a unit. The change `(r,A)↔(r,C)` is bijective, with inverse

    A=(1−p^(h+1))C−r.

**Theorem C.** Under this coordinate change, the aliquot-ray update is exactly

    (r,C) ↦ (p^t r,C).

**Proof.** Substitute `A=(1−p^(h+1))C−r` in (4), using `(1−p)S_t=1−p^t` and `(1−p)S_(h+1)=1−p^(h+1)`. This gives

    B=(1−p^(h+t+1))C−p^t r.

If `h+t≥k−1`, the coefficient is already 1 modulo `p^k`, matching the truncated valuation at the new state; otherwise its new valuation is `h+t`. The zero-residue case also follows, since all the `p^(h+1)` terms then vanish. This proves (6) is invariant and the conjugacy is exact.

Thus the apparently nonlinear finite arithmetic dynamics is a direct product of nilpotent multiplication on the source coordinate with an unchanged core coordinate. Each of its `p^k` fixed points has exactly `p^k` states in its basin. Per fixed core coordinate, the number of states reaching it in **exactly** `j≥1` steps is

    p^min(jt,k) − p^min((j−1)t,k).

Indeed reaching in at most `j` steps means `p^(jt)r=0 mod p^k`, which has `p^min(jt,k)` solutions. There is one depth-zero state per basin. For a nonzero residue of valuation `e<k`, the exact transient length is `ceil((k−e)/t)`.

For an actual source `n=p^e u`, equation (6) is precisely `C=sigma(u)/(1−p)` modulo `p^k`, because `A+r=sigma(n)=S_(e+1)sigma(u)`. The retained sidecar is therefore the scaled divisor sum of the p-free core. The theorem does not claim to recover `u` from that divisor sum.

### Shifted probes incur an exact precision debt

Suppose the two observations are at known positions `j` and `j+t` on the ray, both modulo `p^k`, where `0≤j<k`. Theorem A recovers `p^j n mod p^k`, hence `n mod p^(k−j)`. To guarantee the original residue modulo `p^k` by this method, observe the shifted pair modulo `p^(k+j)`.

The loss of `j` digits is sharp for actual integers, not just formal states. For `1≤j<k`, take

    n_1=p^(k−j),       n_2=p^(k−j+1).

At **every** observed position `ell≥j`, the values `s(p^ell n_1)` and `s(p^ell n_2)` are both `−1/(p−1) mod p^k`, because their occupied exponents are at least `k`. Nevertheless `v_p(n_2−n_1)=k−j`, so the sources differ already modulo `p^(k−j+1)`. Arbitrarily many later observations at the same precision cannot restore the missing digit.

The core-coordinate change is independent of the step `t`, so all positive step sizes are simultaneously conjugated to `(r,C)↦(p^t r,C)`. In particular these formal transitions obey the exact semigroup law `Phi_t Phi_u=Phi_(t+u)`. Retaining the invariant core permits continued reconstruction of the aliquot-observation coordinate `A`; forgetting it loses that reconstruction, although source-residue dynamics remains well-defined.

## 5. Connection to the paper and a strict limitation

The source paper's counting interface multiplies counts by a power of 2, losing all information in sufficiently small 2-primary observations. Here the paired arithmetic observer does the opposite: its geometric-series denominator is always a p-adic unit, and the apparently lost large valuation ceases to matter at exactly the needed precision. Both interfaces must be analyzed at the target quotient, not by the size of their integer intermediates.

Computing `s(n)` is not made efficient by these formulas. This is an exact information and dynamics theorem conditional on the two residues being available, not a new factoring algorithm or a resolution of the paper's headline claim. Composite multipliers remain outside the theorem; Session 110's multiplier-10 and multiplier-15 collisions persist even with exact unbounded observations.
