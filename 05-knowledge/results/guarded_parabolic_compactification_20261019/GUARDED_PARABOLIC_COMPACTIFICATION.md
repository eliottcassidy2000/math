# Guarded parabolic compactification and original-height Collatz repair

19 October 2026. Status: **PROVED, scoped**. The elementary general repetition identity is inherited from the repository's binary-fuel theorem; the new work here applies it to near-neutral bank compactifications, gives a positive repeated-word repair theorem without an episode cap, and characterizes fixed-offset uniform repair by an exact source-sensitive boundary condition. No general Collatz convergence or arbitrary pointwise bank-entry theorem is claimed.

## 1. What the near-cusp idea contributes

For the real quadratic map, shifting the fixed point by `u=x−1/2` gives exactly

    u ↦ u+u²+ε.

Full passage through the real bottleneck has leading clock `π/√ε`. Thus `ε=b^(−2m)` produces a clock of order `b^m`. This is a useful way to encode a large scale in a compact parameter. It does not, by itself, transport the parity guards of a Collatz integer. In particular, a bounded error in an escape clock does not determine an integer's residue class.

There is a genuinely exact parabolic object on the Collatz side. A guarded affine word `B(x)=s x+e`, compactified by `t=1/x`, becomes

    G(t)=t/(s+e t).                                      (1)

At `s=1` this is a parabolic Möbius map, with the exact translation coordinate `1/(e t)` when `e>0`. The word's arithmetic guard remains an additional coordinate. Theorems below explain both a useful consequence of retaining that coordinate and a sharp failure when it is omitted.

The earlier quadratic parameters must also remain distinct. At `c=−29/16`, the rational three-cycle `−7/4→5/4→−1/4→−7/4` has multiplier `35/8`; it is not parabolic. At `c=−7/4`, the algebraic three-cycle defined by `x³+x²/2−9x/4−1/8=0` has multiplier one. These are the exact inherited classifications in [the quadratic-cycle phase note](../quadratic_cycle_phase_20261015b.md); the neutral three-cycle there is not rational. No conjugacy between those scalar quadratic orbits and a Collatz integer orbit is asserted here.

## 2. A positive theorem: a repeated contracting word cannot overshoot

Use the Terras map

    T(x)=x/2                  for x even,
         (3x+1)/2             for x odd.

For a parity word `W` of length `n≥1`, with `h` odd bits, write its branch as

    B_W(x)=(a x+b)/M,    a=3^h, M=2^n, b≥0,
    s=a/M, e=b/M.

A word is *lawful* at x when it is the actual chronological parity word. Every word has one lawful residue modulo M. This follows inductively: the two lifts of a length-n cylinder differ by `2^n`, so their n-step endpoints differ by the odd integer `3^h` and have opposite next parities. Equivalently, its residue is `−b a^(−1) mod M`.

### Theorem A: square-guard descent

If `a<M`, x is a positive integer, and the concatenation `WW` is lawful at x, then

    T^n(x)<x   or   T^n(x)=x.                            (2)

Positive overshoot is impossible.

**Proof.** The two inputs `x` and `y=T^n(x)` both follow W, so `y−x` is a multiple of M. On the other hand, if `y≥x`, then

    0≤y−x=(s−1)x+e<e.

The carry always obeys `e<M` when `a<M`. Indeed, if the odd positions, numbered from zero, are `i_1<...<i_h`, then

    e=Σ_(k=1)^h 3^(h−k) 2^(i_k−n)
      ≤ (1/2) Σ_(k=1)^h 3^(h−k)
      = (a−1)/4 < M.

For h=0 the carry is zero. The only nonnegative multiple of M below e is zero. This proves (2). If equality holds, the actual orbit returns after n steps, so it is periodic with period dividing n. ∎

The periodic alternative cannot simply be removed: W=`10`, x=1 and W=`01`, x=2 are exact returns. The theorem makes no assertion that an arbitrary source repeats W.

### Corollary A1: no episode cap for square-guarded first crossings

Use exactly the inherited bank alphabet

    K=110110110110110101011011101010,
    F=10, H=110101000, V=K^4 F.

Let W begin with K or V and end at its first explicit F/H reset boundary with cumulative multiplier below one. Suppose `WW` is lawful at x. There is **no restriction on the number of expansion episodes in W**. Then either the W endpoint is already below x, or the existing first-crossing zero-return theorem applies: x is the minimum of a cycle of least period n, and a distinct nontrivial positive companion cycle contains `(x+1)/4<x`.

Consequently a least positive nonconverging integer cannot begin `WW` for any such W. This extends the restricted guard language without enlarging a numerical episode limit. It does not supply bank entry or repetition for every integer.

The inherited zero-return argument is used only in its actual scope: fully guarded K,V,F,H words, first crossing at the stated endpoint, the separately lawful adjacent companion, and the forced initial long K guard. In the equality case, the thirteen short descending branches in the guard-refinement note cannot occur, so the source begins `K^5` or `VK^3`, in particular K^4. The companion's first two steps are even. No smaller companion point is misidentified as an iterate of the lower orbit.

### Corollary A2: logarithmically many matching bits suffice

More generally, let ell be a nonnegative integer. If `a<M`, `2^ell≥e`, and the first ell parity bits at `x` and `T^n(x)` agree, then (2) holds. Equality of those bits is exactly equality of the two inputs modulo `2^ell`.

For an interleaved first-crossing bank word with m episodes, the inherited exact bound is `e<mE`, where

    E=83986881009841395726881741282576539890
      /157489106280555890461927692448373459.

Therefore it suffices to match

    ell=ceil(log_2(mE))

initial parity bits at the source and endpoint. This is a conditional self-matching guard, not an automatic consequence of the carry bound. The actual additional source information must be supplied. It retains the original source even when m is unbounded.

## 3. Exact arithmetic cost of the parabolic limit

### Lemma B: inherited repetition fuel in displacement form

Put `D_W=a−M`. For every k≥1,

    W^k lawful at x  ⇔  M^k divides D_W x+b.             (3)

In particular, if W^k is lawful, then

    B_W(x)−x ∈ M^(k−1) Z.                              (4)

When `D_W x+b≠0`, the maximum number of lawful repetitions is exactly

    floor(v_2(D_W x+b)/n).

**Proof.** The repeated affine branch is

    B_W^k(x)=[a^k x+b(a^k−M^k)/(a−M)]/M^k.

The nonempty word has `a≠M`, by unique prime factorization. Multiplication of its numerator by the odd number `a−M` gives

    a^k((a−M)x+b)−bM^k.

Since both `a` and `a−M` are odd, endpoint integrality is equivalent to (3). The unique lawful-cylinder property restores all intermediate parities. Dividing `D_W x+b` by M gives (4). The zero mismatch is the exact fixed-point/cycle case. ∎

This formula is not claimed as a new general pumping identity: it is precisely the forward binary mirror in [the binary-fuel theorem, Section 5](../collatz_guarded_pumping_memory_20261004.md#5-the-forward-mirror-and-a-sharply-located-refill-event). The following scale exclusion is its application here.

### Theorem B: lawful repeated sources avoid the nonlinear crossover

Suppose `s=1+δ>1`, `e>0`, and W^k is lawful at a positive integer x. Then

    δx+e ≥ M^(k−1).                                    (5)

If `M^(k−1)>e`, this implies

    e/(δx) ≤ e/(M^(k−1)−e).                            (6)

Consequently, for a family of words with `M→∞`, `δ→0+`, and bounded positive e, any sources lawfully repeating each word at least twice satisfy

    e/(δx)→0.

They are outside the parabolic crossover `x≈e/δ`, or equivalently `t≈δ/e`.

**Proof.** In the expansion case the one-block displacement is strictly positive. Apply (4). Rearrangement proves (6). ∎

This has a direct normal-form interpretation. Expanding (1) near t=0 gives

    G(t)=t/(1+δ+e t),

so the relative sizes of the parameter term `δt` and nonlinear term `e t²` are measured by `e/(δx)`. Equation (6) forces this ratio to zero for legal double repetitions. Formal iterability of a Möbius map is weaker than arithmetic iterability of its word.

For a contracting word `s=1−α`, the same calculation gives either `x=e/α` exactly, or, whenever `M^(k−1)>e`,

    x ≥ (M^(k−1)+e)/α.

The exact fixed-point case must be retained. Otherwise a lawful repeated contracting word lies above its fixed point and descends, recovering Theorem A.

### The actual near-neutral bank family

Retain

    P=3^19, Q=2^30, C0=4912087879, D=P−Q=88519643,
    c=C0/D, λ=P/Q, μ=3/4,
    θ=log(λ)/log(4/3).

The existing interleaved note proves θ irrational. Choose L→∞ with `{Lθ}→0+`, put `j=floor(Lθ)`, and use

    U_L=K^L F^j.

Its exact multiplier and carry are

    s_L=λ^L μ^j→1+,
    e_L=c s_L+1−(c+1)μ^j→c+1.                         (7)

Thus its reciprocal maps genuinely converge to

    G_0(t)=t/[1+(c+1)t],
    Ψ(t)=1/[(c+1)t],   Ψ(G_0(t))=Ψ(t)+1.              (8)

But `M_L=2^(30L+2j)→∞`, so Theorem B applies to any positive integer source that follows U_L twice. The arithmetic guard removes the fixed crossover scale from these repeated integer realizations.

There is an additional exact distinction even at the limiting parabolic parameter. With `u=−(c+1)t`, the map (8) becomes `u/(1−u)=u+u²+u³+...`, whereas the cusp at ε=0 is `u+u²`. These two germs are not conjugate by an invertible formal power series. Indeed, for `f(u)=u+A u²+B u³+...`, comparison of `h∘f=g∘h` with `h(u)=αu+βu²+γu³+...`, α≠0, gives `A_g=A/α` and `B_g=B/α²`. Thus `B/A²` is invariant: it is one for the Möbius germ and zero for the quadratic germ. This excludes formal/analytic local conjugacy, not every topological or asymptotic comparison.

Moreover, the deformation in (1) has a `δt` term and preserves t=0. It is different from `u↦u+u²+ε`, whose parameter is a constant term and whose fixed-point equation is `u²+ε=0`. A differentiable local coordinate change preserving a fixed parameter conjugates derivatives, so it cannot make a hyperbolic `s≠1` bank fixed point into a multiplier-one fixed point. No transfer of the `π/√ε` passage law follows.

For comparison, in the unguarded affine model set `X=δx/e`. Then

    X_(r+1)=(1+δ)X_r+δ,
    X_r+1=(1+δ)^r(X_0+1).

Moving between two fixed positive X thresholds has clock asymptotic to

    δ^(-1) log[(X_1+1)/(X_0+1)].                       (9)

This is a formal affine scaling law. Theorem B explains why legal repetitions of U_L do not remain in those fixed X windows.

## 4. The K cylinder has both a real scale and an arithmetic phase

For the repeated K word, (3) becomes

    K^L lawful at x ⇔ Q^L divides Dx+C0,
    K^L(x)+c=λ^L(x+c).                               (10)

Let `r_L` be the least nonnegative source residue modulo Q^L. Then

    Dr_L+C0=t_L Q^L,
    r_L/Q^L=t_L/D−c/Q^L.                              (11)

For L≥2, `1≤t_L≤D−1` and

    t_L ≡ C0 Q^(−L) (mod D),
    t_(L+1) ≡ Q^(−1)t_L (mod D).                      (12)

The exact checks give

    D=193·458651, gcd(C0,D)=1,
    ord_D(Q)=146768=2^4·9173.                         (13)

Thus the compact real mantissas `r_L/Q^L` have exactly 146768 subsequential limits `t/D`, one for each phase of the modular cycle. Meanwhile the unscaled source residues converge 2-adically to the **single** rational anchor `−c`.

Proof of the range in (12): `r_L<Q^L`, `Q^L>C0`, and `Dr_L+C0>0` imply `0<t_L<D+1`; the nonzero residue modulo D excludes D. The order certificate checks primality of the factors, `Q^146768=1 mod D`, and nonidentity at the two prime-divisor quotients. This is a finite exact arithmetic fact, not a floating period estimate.

In particular,

    r_L ≥ (Q^L−C0)/D.                                (14)

For an actual positive integer x the repeat budget obeys

    L_max=floor(v_2(Dx+C0)/30),
    Q^L_max≤Dx+C0.                                   (15)

A compact mantissa does not make the consumed binary precision disappear.

### Quantitative version of the rational-anchor repair criterion

Fix a positive integer d. On odd-denominator rationals use the same parity-dependent Terras formulas. Suppose `d−c` has a contracting parity prefix of length q, with branch `(a z+b)/2^q`, `a<2^q`. Set

    R_d=(a d+b)/(2^q−a).

Then **every** positive source in K^L has

    T^q(x+d)<x

whenever

    30L≥q  and  Q^L>D R_d+C0.                         (16)

Indeed, the first q parity bits of x+d agree with those of `d−c`, because `x≡−c mod Q^L`. Equations (14) and (16) put every source above the exact carry threshold R_d. Conversely, an all-lift offset repair using a prefix within a K^L guard forces a contracting prefix of `d−c`. This is the existing qualitative K∞ obstruction with an explicit source-height/precision bound.

The companion checker supplies a small declared control set, not a universal rational stopping theorem: d=4 has q=92 and L=4 suffices; d=12852 has q=135 and L=5 suffices. It also checks the fifteen exceptional offsets from the earlier 1000-episode forest. None of these checks asserts that all positive `d−c`, d≥56, have finite coefficient stopping time.

## 5. An exact compact boundary retaining the missing coordinate

There is a natural simultaneous real/2-adic compactification of positive integers:

    x ↦ (1/x,x) ∈ [0,1]×Z_2.

Its closure is exactly

    X = {(1/x,x):x is a positive integer} ∪ ({0}×Z_2).  (17)

Every 2-adic phase occurs at the same real-height boundary t=0. To prove this, for ξ∈Z_2 choose positive integers `x_k≥2^k` with `x_k≡ξ mod 2^k`. Then `(1/x_k,x_k)→(0,ξ)`. A limit with positive first coordinate has bounded integer sources and therefore is one of the original points.

The Terras map extends continuously to X by

    (t,ξ) ↦ (2t, Tξ)              if ξ is even,
             (2t/(3+t), Tξ)       if ξ is odd.          (18)

These formulas are asserted on X, not on an unrelated full rectangle. On the boundary they induce `(0,ξ)↦(0,Tξ)`. A finite parity word acts in the real direction by (1), and its derivative at t=0 is `1/s=M/a`.

For a positive integer d, the offset x↦x+d acts by

    (t,ξ) ↦ (t/(1+d t),ξ+d).                          (19)

It has the same real boundary value zero but changes the arithmetic phase. This is precisely why an offset can have a very different future parity word even when its relative real-height change is tiny.

The near-cusp scale is compatible with using a compact real coordinate. It does not identify the entire boundary `{0}×Z_2` with a single arithmetic state. The additional coordinate in (17) is a mathematically exact way to retain the information needed by the guard-refinement certificates.

## 6. A necessary-and-sufficient source-sensitive boundary criterion

Let S be any set of positive integer sources and fix a positive integer d. Define its unbounded-source 2-adic boundary by

    Λ∞(S)= intersection over H≥1 of closure_Z2({x∈S:x≥H}).

Equivalently, ξ belongs to Λ∞(S) exactly when, for every H and k, some `x∈S` satisfies `x≥H` and `x≡ξ mod 2^k`. This set is compact; it is empty when S is finite.

A 2-adic number z has *finite coefficient stopping time* if some finite parity prefix has multiplier `3^h/2^j<1`. This definition uses no real ordering of z.

### Theorem C: compact boundary criterion

The following are equivalent:

1. Every `ξ∈Λ∞(S)` has finite coefficient stopping time at `ξ+d`.
2. There exist finite integers J,H such that, for every `x∈S` with `x≥H`, some `1≤j≤J` satisfies `T^j(x+d)<x`.

Thus coefficient contraction on the actual arithmetic boundary is exactly what is needed for a fixed-offset, uniform finite repair at all sufficiently large sources. Small sources remain a finite separate task.

**Proof, 1⇒2.** If Λ∞(S) is empty, S is finite by compactness, and any J≥1 with H above its largest element works (take H=1 if S is empty). Otherwise each contracting prefix is constant on a finite 2-adic cylinder. These open cylinders cover Λ∞(S); compactness supplies finitely many. Let J be their largest prefix length. Their exact thresholds `(a d+b)/(2^j−a)` have a finite maximum R. All sufficiently large sources in S lie in their union: otherwise an unbounded sequence of uncovered sources has a 2-adic convergent subsequence with limit in Λ∞(S), outside that open union, a contradiction. Increase H so that `H>R`. The guarded affine inequality gives the claimed descent for every x≥H.

**Proof, 2⇒1.** Fix ξ in Λ∞(S). Choose a source x≥H matching ξ modulo `2^J`; arbitrarily large such sources exist. The first J parities from x+d and ξ+d agree. If none of those prefixes had contracting multiplier, then for every j≤J, positive carries and d>0 would give `T^j(x+d)≥x+d>x`, contradicting 2. Thus ξ+d has a contracting prefix of length at most J. ∎

The compactness argument supplies a criterion, not an effective algorithm for finding H from an arbitrary description of S. A practical certificate must still exhibit the finite cover and verify the excluded small/tail cases.

### Apply C to the actual bank frontier

For each positive integer d, let S_d be the positive sources admitting a fully guarded first-crossing bank word W with

    T^|W|(x)=x+d.

This set retains the endpoint equation as well as the source prefix. A target for the unbounded-episode problem is

    Λ∞(S_d)+d avoids the coefficient non-stopping set,
    for each d>0,                                     (20)

together with repair of the finitely many small exceptions for each d.

Every unbounded-source limit in Λ∞(S_d) is encoded by an infinite noncrossing bank path. To see this, choose an associated finite first-crossing parse for each source. Bounded-length parses supply only finitely many sources for a fixed d, by the endpoint equation below. Hence the lengths go to infinity on a subsequence. The finite alphabet allows successive common-prefix selection. No first crossing can occur at a fixed prefix on that subsequence, again because that fixed word permits at most one source for displacement d. The resulting infinite word is the claimed path.

The converse inclusion is **not** automatic. The current complete-prefix forests cover all symbolic bank paths, including the path K∞ with 2-adic source `−c`. Such a stronger certificate requires rational coefficient stopping for `d−c`. The source-sensitive condition (20) requires that test at `−c` only if `−c∈Λ∞(S_d)`. Establishing or excluding that membership is a precise additional arithmetic problem. An infinite formal source path is not itself a positive integer and need not be realized by the specified endpoint displacement.

### Exact endpoint filtering

For any fixed contracting word W with branch `(a x+b)/M`, a source with displacement d must be

    x_d=(b−Md)/(M−a).                                (21)

It exists exactly when this is a positive integer. If it is, the endpoint equation makes `a x_d+b` divisible by M and hence authenticates the complete W guard. Thus each fixed word contributes at most one source to S_d.

This gives a constructive refinement option at completed first-crossing leaves:

- discard the leaf for d if (21) is not a positive integer;
- otherwise retain that exact source and verify its subsequent original-height repair;
- do not promote an unrealizable endpoint offset into an actual counterexample.

It does not justify pruning an infinite pre-crossing path without a proof that no future extension realizes d. Equation (20) states exactly the missing limit-set obligation.

The criterion is stronger than individual eventual descent and is not a reformulation of global Collatz convergence. For example, for S equal to all positive integers, Λ∞(S)=Z_2 contains `−d−1`; after translation this is the all-odd fixed point −1. No bounded repair deadline can therefore work for all sufficiently large x. Explicitly, `x=2^J q−d−1>0` gives J initial odd bits at x+d and no descent below x through those J steps. This does not contradict possible later descent.

## 7. Relation to the current results and stopping boundary

The existing results remain unchanged:

- [original-height repair](../collatz_original_height_repair_20261019/ORIGINAL_HEIGHT_REPAIR.md) handles one expansion/reset episode with 92 extra steps;
- [interleaved repair](../collatz_interleaved_height_repair_20261019/INTERLEAVED_HEIGHT_REPAIR.md) handles up to 24 episodes with 92 steps;
- [source guard refinement](../collatz_guard_refined_repair_20261019/GUARD_REFINED_REPAIR.md) handles 25 episodes with 183 steps and its separate exact 1000-episode corollary with 273 steps;
- [height-controlled entry](../collatz_height_controlled_entry_20261019/HEIGHT_CONTROLLED_ENTRY.md) is a density-one, joint height/time theorem and supplies no exceptional-integer coverage.

The present positive addition is the unlimited-episode **self-matching guard** theorem. The structural addition is the arithmetic exclusion of the formal near-parabolic crossover for repeated words. The research criterion isolates the actual unbounded-source limit set rather than replacing source information by a real scaling limit or an almost-everywhere statement.

The exact distinction between coefficient stopping and actual descent is also standard in Rozier–Terracol, *Paradoxical behavior in Collatz sequences*, Definition 1.2 and its surrounding discussion: https://arxiv.org/html/2502.00948v5 . That paper is cited only for terminology and the distinction, not for any rational-source or pointwise bank-entry assertion here.

## 8. Reproduction

Run `python verify_guarded_compactification.py` and repeat with `python -O verify_guarded_compactification.py`. Both outputs must match [the verification receipt](verification-receipt.json). The checker uses only integer/Fraction arithmetic for acceptance tests, verifies the repetition and carry identities, performs declared finite all-lift controls for contracting squares, authenticates the K phase period, and checks the explicit rational-anchor examples. The all-size theorems and compactness criterion are proved above; finite controls are not their quantifiers.
