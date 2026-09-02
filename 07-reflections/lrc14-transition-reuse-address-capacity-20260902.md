# LRC(14): addressed transition capacity and the exact reuse quotient

**Status:** PROVED elementary address-capacity and aggregation lemmas +
FINITE-EXACT controls.  Research reflection only; not canonized here and not a
proof of the minority branch or LRC(14).

## Inheritance pass and concept board

The closest proved mechanism is Section 7 of
[THM-4335](../01-canon/theorems/THM-4335-lrc14-owner-permutation-component-address-and-minority-renewal.md):
every covered nonspanning component of the `{2h}` anchor has a shortest greedy
chain of addressed odd-tail teeth, and each consecutive crossing pays an exact
`gcd`-quantized overlap.  The canonical hostile is its `h=420` pair of rows,
whose `730` nonspanning obligations contain only four actual renewal
components.  The corrected near miss is to charge one obligation to a speed
pair without retaining its two tooth addresses.  The least-used sidecar is the
whole handoff window between the incoming left wall and outgoing right wall.

The concept board was

```text
addressed edge | speed-pair quotient | determinant/q class
anchor address/owner bit | CRT phase orbit | aggregate slot capacity
```

The applicable method cards were “attack a proposed bound before extending
it,” “test structured adversaries,” and “turn certificate failure into an
address.”

## 1. Proper crossings and the handoff window

For an odd speed `w` and integer tooth address `n`, put

```text
a(w,n)=(14n-1)/(14w),       b(w,n)=(14n+1)/(14w).
```

An ordered addressed pair

```text
e:(u,n) -> (v,m)
```

is a proper crossing exactly when

```text
a(u,n)<a(v,m)<b(u,n)<b(v,m).                         (1)
```

Write

```text
D=um-vn,                 q=u+v-14D,                 (2)
g=gcd(u,v).
```

Direct subtraction gives

```text
a(v,m)-a(u,n)=[14D+v-u]/(14uv),
b(v,m)-b(u,n)=[14D+u-v]/(14uv),
b(u,n)-a(v,m)=q/(14uv).                              (3)
```

Consequently `(1)` is equivalent to

```text
|u-v|<14D<u+v.                                       (4)
```

In particular `D>0`, `D in g Z`, and, because `u/g,v/g` are odd,

```text
q in 2g Z_(>0).                                      (5)
```

Now suppose `e` is actually used consecutively by the farthest-reach greedy
chain on

```text
I_k=[(14k+1)/(28h),(14k+13)/(28h)].                  (6)
```

When `(u,n)` was selected, `(v,m)` already had the farther right wall.  It
could not yet have been eligible, so its left wall is at or to the right of
the preceding frontier, hence at or to the right of `L_k`.  Since `(u,n)` is
not the final chain tooth, its right wall cannot pass `R_k`.  Therefore every
actual greedy edge satisfies the sharper located condition

```text
H_e=[a(v,m),b(u,n)] subset I_k.                       (7)
```

The endpoint inclusions in `(7)` are nonstrict: an incoming open tooth may
start exactly at `L_k`, and an outgoing open tooth may end exactly at `R_k`.
Condition `(7)` is necessary, not sufficient in a larger tail bank; a third
tooth can dominate one of the greedy choices.

Because `H_e` has positive length and the interiors of the anchor components
are separated, `(7)` can hold for **at most one** `k`.  Thus a fixed addressed
transition is never reused across anchor components.  Reuse is created only
by the quotient

```text
(u,n;v,m) -> {u,v}.                                  (8)
```

This is the missing logical distinction in the renewal aggregation problem.

## 2. Exact component-address congruence

The unique candidate component in `(7)` can be read using integers only.  Set

```text
Y=2h(14m-1),               X=2h(14n+1),              (9)
Y=14v k_in+r_in,           0<=r_in<14v,
X=14u k_out+r_out,         0<=r_out<14u.
```

Then `(7)` holds if and only if

```text
k_in=k_out=k,
v<=r_in<13v,               u<r_out<=13u.             (10)
```

The strict/nonstrict pattern in `(10)` records exactly the open-tooth
semantics at the two handoff endpoints.

Writing `k=j+h epsilon`, `epsilon in {0,1}`, gives the MC2 owner bit of the
transition.  The quotient tooth integers are

```text
N_u=2n-epsilon u,          N_v=2m-epsilon v,          (11)
N_u=N_v=epsilon mod 2,
uN_v-vN_u=2D.
```

Thus component address, orientation, owner bit, determinant, and `q` all
survive in one exact state.

## 3. Exact speed-pair parametrization and count

Normalize

```text
u=gU,        v=gV,        gcd(U,V)=1.                 (12)
```

The possible normalized determinants are

```text
Dset(U,V)={d in Z_(>0): |U-V|<14d<U+V}.              (13)
```

Their number is

```text
N(U,V)=max(0,
 floor((U+V-1)/14)-floor(|U-V|/14)).                 (14)
```

For `d in Dset`, let `n_d` be the representative in `[0,U)` satisfying

```text
V n_d=-d mod U,
m_d=(d+V n_d)/U.                                     (15)
```

All period-one proper crossings of determinant `D=gd` are exactly

```text
n=n_d+Us,        m=m_d+Vs,        0<=s<g.            (16)
```

Hence an ordered speed pair has exactly `gN(U,V)` unfiltered addressed proper
crossings.  In particular

```text
gN(U,V) <= g ceil(min(U,V)/7).                        (17)
```

This already replaces an unspecified reuse factor by a finite pair-specific
one.

There is also an exact anchor filter.  Put

```text
q_d=g(U+V-14d),
lambda_d=h q_d/(7uv),
eta_(d,s)=h(14(m_d+Vs)-1)/(7v).                      (18)
```

Here `lambda_d` is the length of the handoff window after the coordinate
change `y=2ht`, and `eta_(d,s)` is its incoming left endpoint.  The exact set
of candidate component addresses for orientation `u -> v` and determinant
`gd` is

```text
K_d^(u->v)={ floor(eta_(d,s)) : 0<=s<g,
  {eta_(d,s)} in [1/14,13/14-lambda_d] }.             (19)
```

If `lambda_d>6/7`, this set is empty.  Formula `(19)` is equivalent to the
integer remainder test `(10)`.

For a closed count, define

```text
delta=gcd(g,2h),        G=g/delta,
rho_d={h(14m_d-1)/(7 delta V)}.                       (20)
```

The `g` translated windows form `delta` copies of a shifted `1/G` grid.
Therefore the number of **addressed events**, with multiplicity even when two
events land in one component, is

```text
A_d=delta max(0,
 floor(G(13/14-lambda_d)-rho_d)
 -ceil(G/14-rho_d)+1).                               (21)
```

The component congruences are equally explicit.  Put `H=2h/delta`.  For
`0<=r<G`, every live base event generates

```text
k_(r+G ell)=k_r+H ell,          0<=ell<delta.         (22)
```

Thus `(19)` is a union of arithmetic progressions modulo `H`.  If `g<=2h`,
the base addresses are distinct, so the number of component addresses is
also `A_d`.  For `g>2h`, the exact distinct count is obtained by deduplicating
the live base residues `k_r mod H`; this loss is real because several narrow
pair windows can occupy one anchor component.

Since `delta` is odd, one progression in `(22)` has `(delta+1)/2` events with
owner bit zero when its base residue is below `H/2`, and `(delta-1)/2`
otherwise; the remaining events have owner bit one.  Reflection sends

```text
(u,n)->(v,m), k, epsilon
  to (v,v-m)->(u,u-n), 2h-1-k, 1-epsilon,            (23)
```

preserving `D` and `q`.  Hence the two orientations have the same event
count, determinant histogram, and overlap histogram.

## 4. Global aggregation theorem

For each component whose shortest greedy chain has length `m_k>=2`, tag every
edge by its full addressed state `(u,n;v,m)`.  Section 1 makes these tags
injective across all components.  Therefore

```text
sum_(renewal k)(m_k-1)
 <= C_h(W):=sum_({u,v} subset W) 2 sum_(d in Dset) A_d.  (24)
```

In particular the number `R_h(W)` of renewal components obeys

```text
R_h(W)<=C_h(W).                                       (25)
```

Let `P_h(W)` be the addresses spanned by one tooth, as in THM-4335.  A strict
counterexample has neither a missing point nor an unspanned nonrenewing
component, so necessarily

```text
2h-|P_h(W)| <= R_h(W) <= C_h(W).                      (26)
```

Thus

```text
C_h(W)<2h-|P_h(W)|                                   (27)
```

is a rigorous, finite, exact safety certificate.  This is strong enough to
aggregate every renewal obligation without assuming independence or a
constant speed-pair reuse bound.

There is a `q`-weighted version:

```text
sum_(actual edges e) q_e/(14u_ev_e)
 <= sum_({u,v}) [2/(14uv)] sum_d A_d q_d.             (28)
```

The actual handoff windows are pairwise interior-disjoint: different
components are disjoint, and within one shortest chain a tooth two steps
later does not overlap a tooth two steps earlier.  Consequently the left
side of `(28)` is also at most `mu(G_(2h))=6/7`.  Formula `(28)` supplies the
exact global collision-toll sidecar for summing THM-4335's renewal budgets.

## 5. Exact controls

The companion artifacts are

```text
04-computation/lrc14_transition_reuse_address_capacity_probe_20260902.py
05-knowledge/results/lrc14_transition_reuse_address_capacity_probe_20260902.out
```

They compare `(13)--(21)` against direct wall inequalities on `7,220` small
`(h,u,v)` triples and independently verify reflection.  They also reconstruct
the literal greedy chains.

For `h=10`, `W={3,13}`, the entire candidate capacity is used:

```text
k=6:   (13,4)->(3,1),    epsilon=0, D=1, q=2,
k=13:  (3,2)->(13,9),    epsilon=1, D=1, q=2.         (29)
```

For both `h=420` controls (`P=1287` and `P=9009`), the actual renewal data are
identical.  There are six addressed edges on four components:

```text
k=147: 1365->525 (q=420), 525->11 (q=4), epsilon=0,
k=310: 11->1575 (q=4),                     epsilon=0,
k=529: 1575->11 (q=4),                     epsilon=1,
k=692: 11->525 (q=4), 525->1365 (q=420),  epsilon=1.  (30)
```

The pair `{525,1365}` is the important hostile.  It has `g=105`, one
determinant class, and **210** admissible addressed events on 210 components,
though the greedy chains use only two.  Its two orientations occupy the
progressions `k=3 mod 8` and `k=4 mod 8`.  Thus even inside the inherited
`h=420` data, a speed pair has no small uniform reuse constant.

The full candidate capacities are large:

```text
P=1287: 21,176 filtered events (28,344 unfiltered),
P=9009: 27,404 filtered events (35,628 unfiltered).   (31)
```

So `(27)` does not close either control: it removes the reuse ambiguity but
is intentionally an upper-capacity theorem, and competition among tails is
still doing most of the work.

## 6. Two unbounded-reuse counterfamilies

### 6.1 GCD/CRT all-or-none family

For odd `a`, take the pair `{3a,13a}`.  It has one determinant class,

```text
g=a,       D=a,       q=2a,
overlap=1/[7 lcm(3a,13a)]=1/(273a).                  (32)
```

At `h=10a`, every one of the `a` events in each orientation is admissible and
is literally used by the two-tail greedy cover.  Their addresses are

```text
13a->3a: k=6 mod 20,
3a->13a: k=13 mod 20.                                (33)
```

At `h=8a`, the same `2a` proper pair crossings all lie in anchor-danger gaps,
so the filtered capacity is zero.  This is an exact CRT phase switch, not a
sampling effect.

Common dilation alone would not be a satisfactory primitive hostile.  Add
the odd tail `1`.  On every component in `(33)` lying wholly in the speed-1
safe band `[1/14,13/14]`, the original pair chain is unchanged.  The number
of such components is `12a/7+O(1)` across the two orientations, hence is
unbounded, while

```text
gcd(20a,1,3a,13a)=1.                                 (34)
```

The exact probe records `46` surviving pair renewals already at `a=27`.

### 6.2 Primitive coprime family

Unbounded reuse is not only a large-gcd phenomenon.  For every `N>=1`, put

```text
u=14N+1,       v=14N+3,       h=42N.                 (35)
```

Then `gcd(u,v)=1`, and `(13)` consists of `d=1,...,2N`.  For every

```text
d=N+1,...,2N                                  (36)
```

the unique crossing in **each** orientation is literally a two-tooth greedy
cover of its anchor component.  Hence the primitive three-speed row
`{2h,u,v}` has at least `2N` actual renewals using the same unordered coprime
pair.

Here is a direct proof.  For orientation `u->v`, if `d=2r`, then

```text
n=u-r,       m=v-r,       k=84N-6r,
entry phase=6(3r-N)/v,    exit phase=6(r+N)/u.        (37)
```

If `d=2r-1`, then

```text
n=7N+1-r,    m=7N+2-r,    k=42N-6r+3,
entry phase=3(6r-2N-3)/v,
exit phase=3(2N+2r-1)/u.                             (38)
```

For the range `(36)`, these phases lie in `[1/14,13/14]`; subtracting the
scaled `u`-tooth width `12N/u` puts its left wall strictly before `1/14`, and
adding the scaled `v`-tooth width `12N/v` puts its right wall strictly after
`13/14`.  Thus the greedy chain is exactly `u->v`.  Reflection proves the
opposite orientation.  The probe reaches `62` actual pair renewals at `N=25`
and certifies the stated `50`-event subfamily.

Therefore no aggregation may replace a speed pair by `O(1)` renewal slots,
even after primitive normalization and even when the pair is coprime.

## 7. Scope and remaining frontier

```text
source:       shortest addressed cover chains on the {2h} anchor
target:       determinant/CRT-labelled component slot system
map:          edge -> handoff window -> y=2ht phase and component address
preserved:    actual cover implication, orientation, addresses, owner bit,
              D, q, overlap width, reflection
destroyed:    the speed-pair quotient forgets n,m and greedy competition;
              candidate capacity does not say an edge is selected
sidecar:      handoff window, delta/G/H orbit, and competing-tooth frontier
decisive test: compare C_h(W) with 2h-|P_h(W)|, then inspect only the
              surviving candidate slots under greedy domination
```

The fixed-address reuse problem is solved: its multiplicity is one.  The
fixed-speed-pair problem has an exact finite formula but unbounded families.
The next useful refinement is therefore not another pair-only bound.  It is a
competition theorem showing that most candidate windows are dominated by a
third tooth, as the `210`-candidate/`2`-actual `{525,1365}` control makes
plain.  Endpoint-tooth reuse in the summed width budgets also remains to be
bounded.  Neither `(24)` nor `(28)` alone closes the `420|h` minority wall.

## Reproduction

From the repository root:

```bash
python3 -B 04-computation/lrc14_transition_reuse_address_capacity_probe_20260902.py \
  | diff -u 05-knowledge/results/lrc14_transition_reuse_address_capacity_probe_20260902.out -
```

Normal, `python3 -O`, and `PYTHONHASHSEED=913` replays agree.  Raw-LF SHA-256:

```text
ccfce1ccdd0efcdcdfb254c1208adfdaeefd9834146b461a44a01e1907f25d16  script
77e412732397db357cacdfcfb73de02570c91c6fab510f46b20fe9035e41df25  output
```
