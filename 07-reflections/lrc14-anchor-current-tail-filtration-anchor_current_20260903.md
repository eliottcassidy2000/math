---
title: "LRC14 sharp-current layer cake and the critical septimal fibre transport bound"
status: >
  PROVED ANALYTIC + VERIFIED-EXACT, promoted as THM-4372; LRC(14) remains
  open.  The sharp absolute-current rebate is exactly a
  threshold-three atom plus a third binomial tail moment.  The 7-adic
  reverse martingale and the even-anchor wall combine to give a rigorous
  fibrewise convex-transport lower bound; composition with THM-4370 retains
  at least 4/7 of the deeper-shell rebate on the counterexample locus, and
  that coefficient is sharp for the retained typed-edge fibre data. A
  four-tail physical row is the minimal exact obstruction to the stronger
  naive 6/7-survival assertion.
source: anchor_current / LRC14 continuation session, 2026-09-03
artifacts:
  - 04-computation/lrc14_anchor_current_tail_filtration_anchor_current_20260903.py
  - 05-knowledge/results/lrc14_anchor_current_tail_filtration_anchor_current_20260903.out
canonical_result: 01-canon/theorems/THM-4372-lrc14-sharp-current-layer-cake-and-critical-fibre-transport.md
script_sha256: 5a79bcff44bfde4e9d655f1640fe3ee2078abd4e9a5ef65f3f5f201a44dc414f
output_sha256: 6276649efa7fcde418387cd397bae72a2e42e3d8818b29aa0f8076fb685f5cd2
hash_basis: raw LF bytes
audit: >
  PASS.  A separate integer dynamic program checks all 169 level/budget
  entries against the closed six-safe-vertex envelope.  Exact rational wall
  sweeps verify the theorem on 90 structured/pseudorandom profiles, the
  minimal four-tail hostile, the exact six-edge cancellation, and the first
  positive sharp local fibre, the exact global 4/7 five-tail equality, and a
  primitive twelve-tail 7/2/3 local embedding. Common-dilation tests check
  every current threshold 1,...,12 against the residual-wall formula at
  reduced remainders 0,1,6. Normal, optimized, and nondefault-hash-seed runs
  reproduce the frozen output byte-for-byte.
related:
  - THM-4345-lrc14-halfperiodic-anchor-strip-euclidean-remainder-and-current-envelope
  - THM-4346-lrc14-halfturn-current-brownian-kernel-and-cubic-moment-boundary
  - THM-4370-lrc14-septimal-wall-quadrature-and-valuation-reanchor
  - THM-2216-residual-capacity-hinge-gram-law
  - HYP-3223-lrc14-green-current-lorentzian-exchange-angles
---

# Sharp-current tails and the critical septimal fibre

## 1. Inheritance and concept board

The closest proved mechanisms are the exact half-periodic anchor-strip
transducer in THM-4345 and the exact divisor-square/reverse-martingale current
filtration in THM-4346.  The canonical hostile is the one-shell geometric
chain

```text
W_9={1,9,...,9^11},
```

whose quadratic current energy is below the variance threshold.  The corrected
near miss is that even exact anchor-strip quadratic energy is too coarse; on
the THM-4345 physical control, the sharp nonlinear current envelope is positive
where the variance and positive-part quadratic certificates are negative.  The
least-used sidecar is the Green-current interpretation of a current as vertex
divergence, together with the integer layer-cake precedent of THM-2216.

The live concept board was:

```text
sharp rebate g(|C|)       <-> nested current superlevel sets
7-adic reverse martingale <-> convex order, not tail stochastic order
anchor residual wall      <-> one deleted vertex on a septimal fibre
Green current             <-> unit-edge divergence and transport budget
common dilation           <-> exact residual-tail transducer control
```

The relevant META-PATTERNS card was “correct the object before sharpening the
technique”: replace scalar current energy by its nested tail profile.  The
controlled-forgetting card also applies: the deleted anchor vertex is the
sidecar lost by the full-circle shell norm.

All integrals below use ordinary `dt` on `[0,1/2)`.  Endpoints form a null set.

## 2. Exact binomial and tail-probability form of the sharp rebate

For an integer `d>=0`, THM-4345's sharp current-only envelope is

```text
g(0)=g(1)=g(2)=0,  g(3)=1/2,
g(d)=d(d^2-6d+11)/12                 (d>=4).          (1)
```

It has the simpler exact form

```text
g(d)=1/2 [1_{d>=3}+binom(d-1,3)],                     (2)
```

where the binomial term is zero below its natural range.  Equivalently, its
successive differences are

```text
Delta g(1)=Delta g(2)=0,
Delta g(3)=1/2,
Delta g(j)=1/2 binom(j-2,2)          (j>=4).          (3)
```

Therefore, on any measurable region `E`, if

```text
T_j(E;C)=integral_E 1_{|C|>=j} dt,
```

(After division by `|E|`, these are the corresponding tail probabilities.)

then the current rebate is exactly

```text
integral_E g(|C|)dt
 =1/2 T_3(E;C)
  +1/2 sum_(j=4)^12 binom(j-2,2) T_j(E;C).            (4)
```

Thus the nonlinear signal is not an unspecified higher moment.  It is one
exceptional threshold-three atom plus a positive third-binomial layer cake.
Equation (2) follows by expanding the cubic in (1), and (4) is the standard
integer layer cake applied to (3).

### Exact residual-wall transport of every tail

Let `B` be a finite odd family, `C_B=sum_(w in B) sigma_w`, and let the even
anchor be `2h`.  For an odd common multiplier `m`, use THM-4345's notation

```text
c=gcd(m,h), D=m/c=7q+r, H=2h/c,
R={x: ||Hx+q/2||<r/14}.
```

The indicator `1_{|C_B|>=j}` is half-periodic.  Applying the exact transducer
to it, separately for every `j`, gives

```text
T_j(strip; C_(mB))
 =[q T_j(half;C_B)+T_j(R;C_B)]/D.                    (5)
```

Summing (5) with the nonnegative weights (3) gives the same formula for
`integral g(|C|)`.  This is a specialization of THM-4345, not a new
common-dilation family theorem.  Its value here is that it identifies the
exact missing object as a nested residual-wall tail vector, rather than one
residual energy scalar.

The exact script checks (5) for all thresholds `1,...,12` on the THM-4345
base family at multipliers `1,13,49,127`.  At reduced degree `D=7`, every
tail, hence the whole nonlinear rebate, splits exactly in the ratio `1:6`
between strip and core.  As in THM-4345, the multiplier-`49` row is a
deliberately nonprimitive normalization control, not an unresolved primitive
LRC row.

## 3. Convex order from the 7-adic reverse martingale

Let `Phi` be the even piecewise-linear interpolation of `g(|n|)` from integer
points.  Its nonnegative-side slopes are the successive differences (3),
which are nondecreasing; reflection across zero shows that `Phi` is convex on
the real line.

For a finite odd family `W`, split its current by exact septimal height:

```text
C=sum_(e>=0) C_e,
M_e=sum_(j>=e) C_j=Q_e C,                            (6)
```

with `Q_e` the reverse-martingale conditional expectations from THM-4346.
Jensen gives the unconditional monotonicity

```text
integral Phi(C) >= integral Phi(M_1) >= integral Phi(M_2) >= ... . (7)
```

This is convex order.  It does **not** imply monotonicity of each individual
tail probability.

Now let

```text
G_h(t)=1-A(2ht),        a=nu_7(h).                    (8)
```

For `e+1<=a`, `G_h` is measurable for the sigma-field of `Q_(e+1)`, because it
is invariant under translation by `1/7^(e+1)`.  Weighted conditional Jensen
therefore gives

```text
integral G_h Phi(C) >= integral G_h Phi(M_a).         (9)
```

This safely deletes every shell below the anchor's septimal height.  The next
step, from `a` to `a+1`, is exactly where the anchor ceases to be measurable;
the residual-wall sidecar becomes load-bearing.

## 4. The critical fibre is a directed-edge transport problem

Put

```text
n_a = number of tails w with nu_7(w)=a,
Z=M_(a+1)=sum_(nu_7(w)>a) sigma_w.                   (10)
```

For `s,n>=0`, define

```text
u=min(n,6(s-2)_+),   u=6q+r, 0<=r<=5,
L_n(s)=(6-r)g(s-q)+r g(s-q-1),                       (11)
```

with the second term omitted when `r=0`, and with `g(k)=0` for every integer
`k<=2` when translated arguments occur.

> **Critical-fibre transport bound (THM-4372).** For any finite family of
> positive odd tails and any positive integer `h`, one has
>
> ```text
> integral_[0,1/2) G_h(t) g(|C(t)|)dt
>  >=(1/7) integral_[0,1/2) L_(n_a)(|Z(t)|)dt.       (12)
> ```

### Proof

Equation (9) reduces the left side to `M_a=C_a+Z`. Both `M_a` and `G_h` have
period `1/7^a`; on the corresponding quotient circle, a generic seven-point
fibre is

```text
t+j/7^(a+1),  j=0,...,6,                             (13)
```

On it, the high current `Z` is constant.  The anchor is dangerous at exactly one
vertex.  If `w=7^a u` with `7` not dividing `u`, the seven values
`wt+uj/7` run through a full `1/7` grid.  The lower-sheet danger interval has
length `1/7` and selects exactly one vertex; the upper-sheet danger interval
also selects exactly one, distinct vertex.  Thus every height-`a` runner is a
unit directed edge, and `C_a` is the vertex divergence of an `n_a`-edge
directed multigraph on the fibre.

More precisely, if `p` is its positive vertex and `w=7^a u`, then its negative
vertex `q` obeys

```text
u(q-p)=3 or 4                         (mod 7),        (13a)
```

with the choice determined by which side of the positive tooth centre the
base point occupies.  The relaxation below forgets this residue-labelled
edge type; equation (13a) is the next available faithful sidecar.

Let the common high current be `z`, change all signs if needed so `z>=0`, and
write the six safe divergences as `d_i`.  Each unit edge contributes one
negative endpoint, hence

```text
sum_(i safe) (-d_i)_+ <= n_a.                        (14)
```

If `r_i=(-d_i)_+`, monotonicity and evenness give

```text
g(|z+d_i|) >= g(max(z-r_i,0)).                       (15)
```

The marginal saving from lowering an integer level `k` to `k-1` is
`Delta g(k)`, and these savings increase with `k`.  Therefore a fixed budget
is optimally spent by lowering the six equal stacks as evenly as possible,
stopping at the zero plateau `0,1,2`.  Euclidean division of the spent budget
is exactly (11).  The sum over the six safe vertices is at least `L_(n_a)(z)`.
Averaging the seven-point fibres contributes the factor `1/7` in (12).

This proof explains why the Green-current language is faithful here.  The
useful electrical object is not effective resistance or `L2` energy; it is
the integer divergence together with its edge-transport budget and the
identity of the deleted anchor vertex.

Nothing in the fibre argument is peculiar to the displayed cubic beyond the
closed form (11).  For any even convex lattice cost `phi` that is
nondecreasing on the nonnegative integers, the same proof gives the minimum
of `sum_(i=1)^6 phi(z-r_i)` under a total reduction budget `n_a`; its value is
the six-bin water-fill of `phi`.  If `phi` has a zero plateau through `k`, the
budget stops at `6(s-k)_+`.  The sharp-current rebate is the case `k=2`.

Combining (12) with THM-4345's pointwise sharp-current envelope gives the
proof-facing q-cubic certificate

```text
integral G_h p(min(lower_depth,upper_depth))
 >=B_3(G_h)+(1/7) integral L_(n_a)(|Z|),              (12a)
```

where sheet symmetry makes the two one-sheet cubic integrals equal to
`B_3(G_h)`. On the half-base this equality is supplied by the reflection
`t -> 1/2-t`, which preserves `G_h` and swaps the sheets. Thus the transport
term can be inserted directly into the anchored cubic route; it is not merely
a lower bound for an auxiliary norm.

### Two simple consequences

First, every one of the `n_a` edges can erase at most one of the six initial
copies of `g(s)`, so (11) gives

```text
integral G_h g(|C|)
 >=(6-n_a)_+/7 integral g(|Z|).                      (16)
```

This is useful whenever the critical shell contains at most five tails.

The bound can be pushed farther up the reverse filtration without changing
its form.  Write `n_a=6q+r` and extend `g` by zero on the whole half-line
`x<=2`.  Then

```text
L_(n_a)(s)=(6-r)g(s-q)+r g(s-q-1).                  (16a)
```

This is a nonnegative combination of translates of a convex nondecreasing
function.  Hence the even extension `L_(n_a)(|x|)` is convex.  For every
`b>=a+1`, another conditional Jensen step gives

```text
integral G_h g(|C|)
 >=(1/7) integral L_(n_a)(|M_b|).                    (16b)
```

The choice `b=a+1` is strongest in convex order, while deeper `b` can be
arithmetically simpler.

Second, for twelve tails one has `n_a<=12-N_>` and `|Z|<=N_>`, hence
`n_a<=12-|Z|`.  Since `L_n` decreases with `n`, define

```text
Lambda(s)=L_(12-s)(s),  s=0,...,12.
```

Its exact values are

```text
s:       0 1 2 3  4    5  6  7  8     9   10  11  12
Lambda:  0 0 0 0  2  11/2 15 38 78 279/2 227 345 498. (17)
```

Consequently

```text
integral G_h g(|C|)
 >=(1/7) integral Lambda(|Z|)
 >=(2/7) integral g(|Z|)1_{|Z|>=4}.                 (18)
```

The first inequality in (18) has its own exact tail form.  The nonzero
successive weights of `Lambda` at thresholds `4,...,12` are

```text
2, 7/2, 19/2, 23, 40, 123/2, 175/2, 118, 153.       (19)
```

Thus (18) is already an explicit positive linear functional of the high-shell
tail probabilities.  Its sharp qualitative boundary is visible: level three
can disappear, but every level at least four forces some safe rebate under the
twelve-tail budget.

### Composition with the septimal wall restriction

THM-4370 meets this filtration at the same height `a=nu_7(h)`.  On a strict
counterexample it forces `n_a+N_> <=5`.  If `n_a>=3`, then `N_> <=2` and
`g(|Z|)=0`; otherwise (16) has coefficient at least `4/7`.  Hence every
surviving counterexample satisfies the sharper conditional estimate

```text
integral G_h g(|C|) >=(4/7) integral g(|Z|).        (19a)
```

A positive right side can occur only for `n_a<=2`, `N_> >=3`, and total
upper-cone count at most five.  In particular the six-edge level-three
cancellation below is a genuine sharp boundary for the transport theorem,
but THM-4370 excludes it from the critical shell of a counterexample.  The
remaining difficulty has moved to the at-least-seven lower-shell tails and
their one-sheet cubic debt.

The `4/7` cannot be improved from the typed critical edges alone. For

```text
h=1, W_a={1,3}, W_>={7,21,35},
```

every `0<t<1/490` gives, with the deleted vertex first,

```text
Z=(3,3,3,3,3,3,3),
C_a=(2,-1,0,-1,0,0,0),
C_a+Z=(5,2,3,2,3,3,3).
```

The edges of speeds `1` and `3` run respectively from `0` to `3` and from
`0` to `1`, so both obey `v(q-p)=3 mod 7` and have distinct safe negative
endpoints. The safe cost is `2=4g(3)`. The exact global values are

```text
integral g(|Z|)=1/70,
integral G_1g(|C|)=2/245=(4/7)(1/70).
```

There is also a primitive full-count local embedding at `h=7`:

```text
W_<={5,9,11,15,17,29,31},
W_a={7,21},
W_>={49,147,245}.
```

At `t=1/1000000` it has the required `7/2/3` valuation counts and the same
upper projection; its lower current makes the raw total
`(12,2,2,2,2,2,2)`, so all six raw safe rebates vanish. This is a local
control, not a strict counterexample. It shows that any stronger
counterexample-conditional coefficient must use global coherence or
one-sheet data, not a pointwise refinement of the typed critical fibre.

## 5. Hostiles, equality controls, and the strongest survivor

### Minimal exact failure of naive `6/7` survival

The tempting assertion

```text
integral G_h g(|C|) >=(6/7) integral g(|Z|)          (20)
```

is false even for a physical integer-speed row.  Take anchor `2` (`h=1`) and

```text
W={1,7,21,35}.                                       (21)
```

The exact wall sweep gives

```text
high-shell rebate                         =1/70,
core rebate of {7,21,35} alone            =3/245=(6/7)(1/70),
core rebate after adjoining speed 1       =1/98,
naive slack                               =-1/490.   (22)
```

The corrected bound (12) is exact on (21): `n_a=1` and its right side is
`1/98`.  The witness is minimal in tail count.  Since `g` vanishes through
current two, at least three high runners are needed; with only those runners,
the reduced degree is seven and exact independence holds.  One critical-shell
runner is the smallest possible perturbation.  Among distinct positive odd
speeds, `7,21,35` are also the three smallest high runners and `1` is the
smallest critical runner.

There is a twelve-tail version with the smallest possible maximum speed for a
three-runner high shell:

```text
W={1,5,7,9,11,13,17,19,21,29,31,35},
core rebate =23423/8369690,
slack in (20)=-79063/8369690.                        (23)
```

### Complete loss at level three and sharp level-four survival

At `t=1/10000` on the anchor's seven-fibre, the nine-tail family

```text
{1,3,5,7,9,11,13,21,35}
```

has exactly six critical-height runners, constant high current `3`, and
exact-shell/total fibre currents

```text
(6; -1,-1,-1,-1,-1,-1),
(9; 2,2,2,2,2,2).
```

All six safe rebates vanish. Hence six edges really suffice to erase the
level-three signal. The twelve-tail padding

```text
{1,3,5,7,9,11,13,15,17,19,21,35}
```

has constant high current `3`, while its total currents on the danger vertex
and six safe vertices are

```text
(12; 1,1,1,2,2,2).
```

Again all six safe rebates are zero, but this row is not needed to witness the sharp
edge-budget threshold.

With

```text
{1,3,5,7,9,11,13,17,19,21,35,49},
```

the corresponding high current is `4` and the total fibre is

```text
(12; 2,2,3,3,3,3).
```

Its safe cost is exactly `2=Lambda(4)`, so the first positive entry in (17)
is physically sharp at the fibre level.

### Canonical one-shell hostile

The THM-4346 family `W_9={1,9,...,9^11}` occupies one septimal shell.  For
`h=1`, the high projection `Z` in (10) is zero.  If `7|h`, the projection
`M_a` after the safe Jensen reductions is already zero.  Hence (12) is exactly
silent on this canonical hostile.  The result is a genuine mixed-height
certificate, not a closure of the arbitrary one-shell branch.

The strongest survivor of the failed implication (20) is therefore (12), not
a different universal survival fraction.  It retains exactly the coordinate
that (20) forgot: the number of critical-height unit edges that can transport
current into the deleted anchor vertex.

## 6. Scope and next sharp problems

The result supplies a positive anchored lower bound whenever either

```text
n_a<=5 and the high-shell rebate is nonzero, or
the high-shell current reaches absolute level at least four.
```

It does not prove LRC(14).  The one-sheet Bonferroni term can still be
negative, the high projection can vanish, and a lower-shell version of the
canonical one-shell hostile escapes.  The transport relaxation also forgets
the residue-labelled edge endpoints and their coherence as `t` moves.

The next plausible refinements are:

1. keep the one-sheet cubic debt jointly with the fibre divergence, rather
   than bounding only the rebate;
2. retain residue-labelled edge types to see whether physical edge coherence
   improves `L_n` away from the two sharp local witnesses;
3. iterate the nested tail vector through the residual-wall transducer before
   scalarizing it; and
4. combine a private high-shell Fourier mode with the threshold-four tail
   functional, since quadratic energy alone does not force that tail.

Connection contract:

```text
source:     7-adic current reverse martingale plus sharp g(|C|)
target:     anchor-core nonlinear rebate
map:        critical seven-fibre; shell runners become directed unit edges
preserved:  integer current, exact height, tail levels, deleted anchor address
destroyed:  residue-labelled endpoints and time-coherent edge motion
sidecar:    critical-shell edge count / divergence transport budget
test:       exact wall sweep, minimal four-tail hostile, D=7 dilation control
```

Reproduce with

```bash
python3 -B 04-computation/lrc14_anchor_current_tail_filtration_anchor_current_20260903.py
```
