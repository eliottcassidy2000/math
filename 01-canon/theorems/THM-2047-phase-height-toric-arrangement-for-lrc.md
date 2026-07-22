---
id: THM-2047
title: The phase-height toric arrangement is an exact carrier for the lonely-runner max-min
status: PROVED (carrier, vertex law, deletion identity, Euler-characteristic detector, and Fejer-regularized relation-lattice formula). The proposed Wall-A localization theorem remains open.
source: codex-2026-07-21-LRC-arrangement-audit
depends_on:
  - THM-1002
  - THM-1142
  - THM-1017
  - THM-2043
related:
  - HYP-8830
  - HYP-7310
  - HYP-2986
  - HYP-3025
  - THM-752
  - THM-668-pair-sum-ruler-witness-structure
---

# THM-2047 -- the phase-height toric arrangement

> **Prior-art boundary.** HYP-2986 already constructs the threshold-`1/14`
> signed endpoint topes and boundary cocircuits. HYP-3025 already keeps the
> individual-arc Cech nerve and proves why the runner quotient needs a Betti
> defect sidecar. THM-1142 gives the exact essential deletion/replacement
> containment, THM-752 supplies the one-danger-tooth interval exit, and the
> pair-sum THM-668 plus its medial-axis reflection already identify the
> winding-circle contact geometry. This
> theorem does not rename those results as new. Its added content is their
> general-`delta` phase-height synthesis, the local top-vertex normal form and
> boundary-layer coefficient, and the exact distinction between the
> Fejer-regularized Fourier annihilator and an Orlik--Solomon layer formula.

Let `T=R/Z`, write `||x||` for distance to the nearest integer, and let `S` be
a nonempty finite set of positive integral speeds. Define

```text
f_S(t) = min_{v in S} ||v t||,
M(S)   = max_{t in T} f_S(t),
G_delta(S) = {t in T : f_S(t) >= delta},
E_S = {(t,delta) in T x [0,1/2] : delta <= f_S(t)}.
```

The cylinder `T x [0,1/2]` is cut by the `2|S|` signed character walls

```text
H_v^+ : v t - delta in Z,
H_v^- : v t + delta in Z.                         (1)
```

They are ordinary toric walls with character vectors `(v,-1)` and `(v,+1)`
after the height coordinate is viewed modulo one, restricted back to the
height half-cylinder. The important object is not the complement of these
walls: it is the **oriented cell subcomplex** `E_S`, including its boundary.

## 1. Exact carrier theorem

**Theorem.** For every `0<=delta<=1/2`,

```text
E_S intersect (T x {delta}) = G_delta(S) x {delta}.       (2)
```

Consequently

```text
M(S) = max{delta : E_S has a point of height delta},       (3)
```

and, for thirteen relative speeds,

```text
LRC(14) for S
  iff E_S intersects the horizontal slice delta=1/14.      (4)
```

*Proof.* Equations (2) and (3) are the definitions, and the maximum exists
because `E_S` is compact. Equation (4) is the weak lonely-runner inequality
`M(S)>=1/14`. The use of a closed oriented subcomplex is essential: a tight
witness can be an isolated boundary point. ∎

Thus `(owner,sign,t,delta)` is a lossless local address. Projecting away
`delta`, the integral lift of `t`, or the owner/sign labels is not justified
without a separate theorem. THM-2043 supplies an infinite exact warning: full
period-14 parity--Hasse data, all tests through `q=13`, `q_threshold`, and any
fixed finite 7-adic height truncation can agree while a resolved
`(q,a,margin)=(41,17,1)` exit differs.

## 2. Top vertices and the pair-sum law

**Theorem.** Every maximizer `t_*` is supported either

1. by two active walls of opposite sign, `H_{v_i}^+` and `H_{v_j}^-`, or
2. by the two sides of one cusp, the self-pair `v_i=v_j`.

In either case, if `t_*=a/q` in lowest terms, then

```text
q divides v_i+v_j <= 2 max(S).                             (5)
```

Moreover `M(S)=r/q`, where `r=min_v |a v|_q` is integral. Hence enumerating
every numerator on every pair-sum ruler computes `M(S)` exactly.

*Proof.* The function `f_S` is a continuous piecewise-linear lower envelope.
It is positive somewhere, so a maximizer is not a zero cusp. Away from the
cusps of the individual triangular waves, every active branch has nonzero
slope `+v` or `-v`. At a local maximum the active slopes cannot all have the
same sign: moving a sufficiently small distance in the improving direction
would increase every active branch, while the inactive branches remain slack.
Thus an active rising branch and an active falling branch meet. Write them as

```text
delta = v_i t - m,
delta = n - v_j t.
```

Adding gives `(v_i+v_j)t=m+n`. If the maximum occurs at a cusp of one active
wave, its two sides give `2v_i t in Z`, which is the same calculation with
`i=j`. Reduction of the resulting rational proves (5). At `t=a/q`, every
distance is `|av|_q/q`, proving the last assertion. ∎

This is the geometric form of THM-1002. It also explains why the correct
arrangement must include height: pair **sums** arise from intersecting a `+`
owner wall with a `-` owner wall.

### Local wedge and boundary-layer coefficient

At a top vertex `t_*`, let `A(t_*)` be the slopes of all active local affine
branches, including both slopes at an active self-cusp, and put

```text
s_-(t_*) = min A(t_*) < 0 < max A(t_*) = s_+(t_*).
```

For all sufficiently small `u`, the local normal form is

```text
f_S(t_*+u) = M(S) + min_{s in A(t_*)} s u
           = M(S) + s_- u,   u>0,
             M(S) + s_+ u,   u<0.                            (6)
```

Consequently, for all sufficiently small `epsilon>0`, normalized Haar length
satisfies the exact boundary-layer law

```text
|G_{M-epsilon}(S)|
 = epsilon sum_{t_*:f_S(t_*)=M}
     (1/s_+(t_*) + 1/(-s_-(t_*))).                           (7)
```

*Proof.* Inactive branches have positive slack at `t_*`, so on a small common
linear refinement only the active affine branches can realize the minimum.
For positive `u` the smallest slope wins; for negative `u` the largest slope
wins. The strict sign inequality is the local-maximality condition. The
piecewise-linear function has no constant interval, hence finitely many
isolated maximizers. Choose disjoint neighborhoods of them and use the compact
gap below `M` on their complement. In each neighborhood the superlevel set is
the interval `[-epsilon/s_+, epsilon/(-s_-)]`; summing its lengths gives (7).
∎

Only the extreme rising and falling owners enter (7); intermediate coincident
blocks are locally redundant. All constraints share one transverse variable.
Thus the proposed S209 Cartesian product of lower-rank layer complements has
the wrong local dimension and multiplicative order after restriction to the
LRC orbit. Formula (7), not a block product, is the exact first-order object.

## 3. Exact deletion/restriction identity

For one speed `w`, put

```text
B_w = {(t,delta): delta<=||w t||},
D_delta(w) = {t: ||w t||<delta}.
```

Then, exactly,

```text
E_{S union {w}} = E_S intersect B_w,                        (8)
G_delta(S union {w}) = G_delta(S) \ D_delta(w).             (9)
```

In particular, if a thirteen-speed set is written `C union {w}`, it fails the
weak `1/14` conclusion precisely when

```text
G_{1/14}(C) is contained in D_{1/14}(w).                    (10)
```

The left side is a closed finite union of arcs whose endpoints lie on the
signed walls of the twelve-speed core; the right side is a union of `w` open
danger arcs. This is a lossless deletion formulation, including isolated
boundary witnesses. THM-1017 solves (8) when `C=d{1,...,12}`. Proving that the
remaining compact-covering residual forces that core, or finding a cell-word
invariant that does, is still Wall A rather than a consequence of arrangement
topology alone.

## 4. Euler characteristic detects the boundary that volume misses

**Theorem.** If `0<delta<=1/2`, then `G_delta(S)` is a proper compact finite
union of closed arcs and points in `T`. Consequently

```text
chi(G_delta(S)) = number of connected components of G_delta(S),             (11)
G_delta(S) is nonempty iff chi(G_delta(S))>0.                               (12)
```

In particular, LRC(14) for `S` is equivalent to
`chi(G_{1/14}(S))>0`.

*Proof.* Each single-runner danger set `{t:||vt||<delta}` is a finite union of
open arcs with signed-wall endpoints. Its finite union over `v in S` has a
finite endpoint set, so the closed complement `G_delta(S)` is a finite union
of closed arcs and isolated points. It is proper because `t=0` is dangerous
when `delta>0`. Every component is therefore contractible and contributes one
to Euler characteristic. ∎

This is the rigorous topology transfer. Unlike Haar volume, Euler
characteristic sees an isolated tight phase. The signed cyclic endpoint word
computes it exactly and retains owner data, making it a plausible target for
deletion--restriction or a nerve calculation. No current argument forces its
positivity for every thirteen-speed set.

## 5. What the relation-lattice formula really says

Let `g_delta(x)=1[||x||>=delta]`. Its Fourier coefficients are

```text
g_hat(0)=1-2delta,
g_hat(k)=-sin(2 pi k delta)/(pi k),  k!=0.                  (13)
```

If `sigma_R g_delta` is the `R`-th Fejer mean, termwise integration gives

```text
integral_T product_j sigma_R g_delta(v_j t) dt
 = sum_{k dot v=0} product_j
     (1-|k_j|/(R+1))_+ g_hat(k_j).                         (14)
```

Fejer means lie in `[0,1]` and converge to `g_delta` away from finitely many
boundary phases on the orbit. Dominated convergence therefore yields the
exact regularized identity

```text
|G_delta(S)| = lim_{R->infinity} (right side of (14)).      (15)
```

This is the useful THM-1820 relation-lattice pairing. The lattice
`{k:k dot v=0}` is the annihilator of the one-parameter subtorus
`t -> (v_1t,...,v_nt)` under the dual character map. Equation (15) is a
Fourier/Poisson selection law with sinc weights; it is not, merely by being a
lattice sum, an arithmetic-Mobius sum over a De Concini--Procesi layer poset.

Nor does safe measure decide the weak conjecture: when `M(S)=delta`, the fiber
`G_delta(S)` may be nonempty but finite and hence have measure zero. The tight
AP at the lonely-runner threshold is the mandatory control. Positive measure
is a strict-exit certificate, not an iff criterion for weak loneliness.

## 6. Scope correction to the S209 arrangement proposal

Three distinctions are forced by the theorem.

1. A standard toric-arrangement complement removes codimension-one hypertori,
   normally a zero-measure set. `G_delta(S)` instead removes positive-width
   coordinate slabs and pulls the result back to one one-dimensional subtorus.
2. On a torus the affine offsets `0` and `1` in a Shi wall coincide. More
   importantly, LRC safety is bounded by the coordinate walls (1), not by the
   pair-collision walls `x_i-x_j in {0,1}`. The verified Shi point counts and
   parking-function region counts therefore do not count LRC safe phases.
3. The S209 statistic `N_R` is a bounded-box count of short vectors in one
   relation lattice, checked only for primitive triples in a small finite
   speed box. It is not an Orlik--Solomon Betti number or arithmetic-Mobius
   mass, and it proves no AP extremality for twelve-speed cores.

The surviving inspiration is precise: use the **signed phase-height toric
cell complex** and its deletion identity, not an unlabelled complement
cohomology. It retains exactly the resolved height and ownership that the
local Hasse and scalar-energy quotients lose.
