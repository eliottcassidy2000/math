---
id: THM-2047
title: The labelled phase-height cell complex is an exact carrier for the lonely-runner max-min
status: PROVED (carrier, top-vertex and boundary-layer laws, paired deletion, Euler detector, Fejer relation formula, and full-lattice reconstruction). Wall-A localization remains OPEN.
source: codex-2026-07-21-LRC-arrangement-audit
depends_on:
  - THM-1002-pair-sum-denominator-bound-and-the-bounded-gap-case
  - THM-1017
  - THM-1142
  - THM-2043
related:
  - HYP-8830
  - HYP-7310
  - HYP-2986
  - HYP-3025
  - THM-752
  - THM-668-pair-sum-ruler-witness-structure
---

# THM-2047 -- the labelled phase-height carrier

> **Prior-art boundary.** HYP-2986 constructs the threshold-`1/14` signed
> endpoint topes and boundary cocircuits. HYP-3025 retains the individual-arc
> Cech nerve and explains why a runner quotient needs a Betti-defect sidecar.
> THM-1142 gives the essential deletion/replacement containment, THM-752 the
> one-danger-tooth interval exit, and THM-668/1002 the pair-sum contact law.
> The new content here is the general-height synthesis, local top normal form,
> exact Euler/Fourier readings, and the sharp distinction between the full
> embedded relation lattice and lossy arrangement summaries.

Let `T=R/Z`, write `||x||` for distance to the nearest integer, and let `S` be
a nonempty finite set of positive integral speeds. Define

```text
f_S(t) = min_{v in S} ||v t||,
M(S) = max_{t in T} f_S(t),
G_delta(S) = {t in T : f_S(t)>=delta},
E_S = {(t,delta) in T x [0,1/2] : delta<=f_S(t)}.
```

The owner-constraint walls in the height half-cylinder are

```text
H_v^+ : v t-delta in Z,
H_v^- : v t+delta in Z.                              (1)
```

After periodizing `delta` modulo one, these are character hypertori with
vectors `(v,-1)` and `(v,+1)` in the full two-torus. Restricted back to the
strip, each connected wall segment retains its owner, sign, height, and which
adjacent cells satisfy the inequality. There is no globally ordered side of a
torus hypertorus. The literal boundary of `E_S` also contains the ambient
bottom `delta=0`. The carrier is this labelled, cellwise selected subcomplex,
not the ordinary complement of the wall set.

## 1. Exact carrier

For every `0<=delta<=1/2`,

```text
E_S intersect (T x {delta}) = G_delta(S) x {delta}.        (2)
```

Consequently

```text
M(S)=max{delta:E_S has a point of height delta}.           (3)
```

For thirteen relative speeds,

```text
LRC(14) for S
  iff E_S intersects the horizontal slice delta=1/14.     (4)
```

**Proof.** Equation (2) is the definition. The continuous function `f_S`
attains its maximum on the compact circle, and `E_S` is its subgraph, proving
(3). Equation (4) is exactly the weak lonely-runner inequality
`M(S)>=1/14`. The closed selected complex matters: at equality the safe fiber
may consist only of isolated boundary points. ∎

Thus `(owner,sign,t,delta,selected-cell)` is a lossless local address.
Projecting away height, owner/sign, or the resolved phase is legal only after
a separate preservation theorem. THM-2043 gives the hostile control: complete
period-14 parity--Hasse data and every fixed finite 7-adic height prefix can
agree while a resolved off-period strict exit differs.

## 2. Top vertices and pair-sum rulers

Every maximizer `t_*` is supported either by active walls of opposite slope or
by the two sides of one active half-integer cusp. If `t_*=a/q` is reduced, then

```text
q divides v_i+v_j <= 2 max(S)                              (5)
```

for active owners `v_i,v_j`, allowing `i=j`. Moreover

```text
M(S)=r/q,       r=min_{v in S}|a v|_q in Z.
```

Hence enumerating every numerator on every pair-sum ruler computes `M(S)`
exactly.

**Proof.** First, `M(S)>0`: choose a phase outside the finite union of the
zero sets of the waves. Refine the circle at all cusps and all crossings of
the finitely many affine branches of `||v t||`. On each open refined cell the
active branches have nonzero slopes `+v` or `-v`. Their lower envelope cannot
be constant on an interval, since one of finitely many affine branches would
then be constant on a subinterval.

At a local maximum the active smooth slopes cannot all have the same sign;
moving a sufficiently small distance in the common improving direction would
increase every active branch while inactive branches remained slack. Thus a
rising branch of owner `v_i` and a falling branch of owner `v_j` meet. For
integers `m,n`,

```text
v_i t_*-m = n-v_j t_*,
```

so `(v_i+v_j)t_*=m+n`. If an active branch is at a cusp, positivity excludes
the zero cusp; the half-integer cusp gives `2v_i t_* in Z`, the self-pair
case. Reduction proves (5). At `a/q`, each distance is `|av|_q/q`, proving the
last assertion. ∎

### Local wedge and boundary-layer coefficient

At a maximizer let `A(t_*)` be the slopes of all active affine branches,
including both slopes of an active self-cusp, and put

```text
s_-(t_*)=min A(t_*)<0<max A(t_*)=s_+(t_*).
```

For all sufficiently small `u`,

```text
f_S(t_*+u)=M(S)+s_-u,   u>0,
f_S(t_*+u)=M(S)+s_+u,   u<0.                              (6)
```

Therefore, for every sufficiently small `epsilon>0`, normalized Haar length
satisfies

```text
|G_{M-epsilon}(S)|
 = epsilon sum_{t_*:f_S(t_*)=M}
     (1/s_+(t_*)+1/(-s_-(t_*))).                          (7)
```

**Proof.** Inactive branches have positive slack, so on a common small linear
neighborhood only active branches can realize the minimum. For positive `u`
the smallest slope wins and for negative `u` the largest wins, proving (6).
The lower envelope has no constant interval, hence only finitely many isolated
maximizers. Choose disjoint neighborhoods and use the compact gap below `M`
on their complement. Each superlevel interval has the two displayed lengths;
summing gives (7). ∎

Only the extreme rising and falling owners enter the boundary coefficient;
intermediate coincident blocks are locally redundant. This additive
one-dimensional wedge, not a Cartesian product of layer complements, is the
correct first-order object on the LRC orbit.

### What the full-torus determinant does and does not count

For two linearly independent periodized characters `(v,sigma)` and
`(w,tau)`, their full-two-torus intersection index is

```text
|det((v,sigma),(w,tau))|=|v tau-w sigma|.                 (8)
```

Same signs expose `|v-w|`, opposite signs expose `v+w`, and the self-opposite
pair exposes `2v`. If the determinant vanishes, the intersection is
non-proper, not an intersection of multiplicity zero. The index also need not
equal the number of representatives retained in `0<=delta<=1/2`: the
opposite-sign characters `(1,+1)` and `(2,-1)` have full-torus index three but
only two points in the closed strip. Thus the determinant displays pair-sum
arithmetic; the top-slope proof, not a raw intersection count, selects the
actual LRC ruler.

## 3. Exact paired deletion

For a speed `w`, set

```text
B_w={(t,delta):delta<=||w t||},
D_delta(w)={t:||w t||<delta}.
```

Then exactly

```text
E_{S union {w}}=E_S intersect B_w,                        (9)
G_delta(S union {w})=G_delta(S)\D_delta(w).              (10)
```

Equivalently, deleting `w` removes its pair of signed walls and its single
inequality. For a thirteen-speed set `C union {w}`, failure of the weak
`1/14` conclusion is precisely

```text
G_{1/14}(C) subset D_{1/14}(w).                           (11)
```

The left side is a closed finite union of arcs and points; the right side is a
union of open danger arcs. This retains isolated witnesses. THM-1017 solves
the extension when `C=d{1,...,12}`. Nothing in (9)--(11) forces an arbitrary
maximum-deletion core to be an AP.

## 4. Euler characteristic detects what volume misses

If `0<delta<=1/2`, then `G_delta(S)` is a proper compact finite union of closed
arcs and isolated points in `T`. Hence

```text
chi(G_delta(S))=# connected components of G_delta(S),     (12)
G_delta(S) is nonempty iff chi(G_delta(S))>0.              (13)
```

In particular, LRC(14) for `S` is equivalent to
`chi(G_{1/14}(S))>0`.

**Proof.** Each danger set `{t:||vt||<delta}` is a finite union of open arcs
with endpoints on the signed walls. Their finite union has finitely many
endpoints, so its closed complement is a finite union of arcs and points. It
is proper because `t=0` is dangerous. Every component of a proper subset of
the circle of this form is contractible and contributes one to Euler
characteristic. ∎

The owner-labelled cyclic endpoint word computes this Euler characteristic
and sees an isolated tight phase, unlike Haar volume. No current argument
forces its positivity for every thirteen-speed set.

## 5. Fejer relation formula and the full embedded lattice

Let `g_delta(x)=1[||x||>=delta]`. Its Fourier coefficients are

```text
g_hat(0)=1-2delta,
g_hat(k)=-sin(2 pi k delta)/(pi k),   k!=0.                (14)
```

If `sigma_R g_delta` is the `R`-th Fejer mean, finite expansion and termwise
integration give

```text
integral_T product_j sigma_R g_delta(v_j t) dt
 = sum_{k dot v=0} product_j
     (1-|k_j|/(R+1))_+ g_hat(k_j).                        (15)
```

The Fejer means lie in `[0,1]` and converge to `g_delta` away from the finitely
many boundary phases on the orbit. Dominated convergence yields the exact
regularized identity

```text
|G_delta(S)|=lim_{R->infinity}(right side of (15)).       (16)
```

This is a Fourier/Poisson annihilator sum, not an arithmetic-Mobius sum over a
De Concini--Procesi layer poset. Positive measure is exactly a strict-exit
certificate: if `M(S)=delta`, the finite maximizer set can be nonempty with
measure zero; if `M(S)>delta`, continuity supplies a positive interval.

After ordering `S` as `mathbf v=(v_1,...,v_n)`, its full embedded annihilator

```text
L_mathbf_v={k in Z^n:k dot mathbf v=0}
```

is lossless up to common scale:

```text
(L_mathbf_v tensor Q)^perp=Q mathbf v.                    (17)
```

Thus the embedded labelled lattice determines the primitive positive speed
vector. Coordinate relabelling becomes an ambiguity only after coordinate
labels are forgotten. If `S=dS_0`, surjectivity of `t->dt` on `T` gives
`M(S)=M(S_0)`, so the omitted scale does not affect the maximum.

Compressed relation data can be genuinely lossy. The lattices for `(1,2)` and
`(1,3)` are abstractly both `Z`, while their maxima are `1/3` and `1/2`.
Given any coefficient cutoff `B`, choose even `N>B`. The vectors `(1,N)` and
`(1,N+1)` both have no nonzero annihilating relation of sup-norm at most `B`,
yet `M(1,N)<1/2` while `M(1,N+1)=1/2` at `t=1/2`. THM-2043 likewise shows
that a fixed residue/height packet need not determine a global strict exit.

## 6. Scope correction and open application

Three distinctions are now exact.

1. A standard toric complement deletes codimension-one hypertori. LRC safety
   deletes positive-width coordinate slabs and pulls them back along one
   one-dimensional orbit; its selected phase-height cells are extra state.
2. On a torus the Shi offsets `0` and `1` coincide, and Shi pair-collision
   walls are not the owner constraint walls (1). Braid/Shi region formulas do
   not count LRC safe phases.
3. A bounded-box count of short vectors in one relation lattice is neither an
   Orlik--Solomon Betti number nor an arithmetic-Mobius mass, and a finite
   triple census proves no AP extremality for twelve-speed cores.

The live Wall-A question is therefore precise: can paired deletion or a nerve
calculation that retains owner, sign, selected cell, height, and endpoint word
force the AP maximum-deletion core required by THM-1017? No ordinary
complement invariant is currently known to preserve that top-cell predicate.
