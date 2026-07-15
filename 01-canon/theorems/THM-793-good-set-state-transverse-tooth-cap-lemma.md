---
id: THM-793
title: Good-set-state transverse-tooth cap lemma and uniform lacunary terminal cones for every far-count flag
status: PROVED (elementary circle geometry plus THM-755; uniform corollaries use THM-780 and settled lower-runner LRC) + FINITE-EXACT (the 2,002-core factor-19 enhancement)
source: codex-2026-07-14-S2
renumber_note: claimed locally as THM-780 and then THM-784; renumbered to
  THM-793 after integrating the live THM-779--791 namespace.
depends_on:
  - THM-755   # capped-envelope criterion
  - THM-780   # strict-margin phase floor for the uniform corollaries
related: [THM-752, THM-761, THM-766, THM-792, THM-777, HYP-6815, HYP-6830]
verification: 04-computation/lrc14_lacunary_far_flag_certificate_codex_S2.py
  (+ 05-knowledge/results/lrc14_lacunary_far_flag_certificate_codex_S2.out)
---

# THM-793 - Good-set-state transverse-tooth cap lemma

## Statement

Let `B` be a nonempty finite set of distinct positive integer speeds. Write

```text
mu  = |G'_B|,
r_B = full topological component count of G'_B.
```

For a positive integer `N` not in `B`, put `P_N=B union {N}`, and write
`mu_N=|G'_{P_N}|` and `r_N=r_top(P_N)`. Then

```text
mu_N >= 6mu/7 - 2r_B/(7N),              (1)
r_N  <= N+r_B.                           (2)
```

This is the invariant-level form: the state `(|G'_B|,r_B)` controls a
transverse frequency insertion without retaining every endpoint.
Here `r_B` and `r_N` include isolated equality points. THM-755 needs only the
number `r_+` of positive-length interval components, so `r_+<=r_N` makes this
state deliberately conservative and every cap consequence below valid.
In particular, if `mu>0`, then

```text
liminf_(N->infinity) mu_N >= 6mu/7 > 0.                         (2a)
```

Thus one high-frequency insertion over a fixed positive-mass base can never be
a safe-mass-decay family. In the LRC14 use with `|B|=11`, settled lower-runner
cases give `M(B)>=1/12>1/14`, hence `mu>0` by continuity.

There is a useful marked-window form. If `G'_B` contains a connected closed
circle interval `J` of length `L>0`, and `S=sum(B)`, then

```text
mu_N >= 6L/7 - 2/(7N),                  (3)
r_N  <= N+S.                             (4)
```

If `mu>0`, then for any integer `a>=2` with all speeds in
`B union {N,aN}` distinct,

```text
r_N/(aN mu_N)
 <= (N+r_B)/(aN(6mu/7-2r_B/(7N)))       (5)
```

whenever the denominator is positive. In particular,

```text
limsup_(N->infinity) r_N/(aN mu_N) <= 7/(6a mu).   (6)
```

The marked-window version replaces `(mu,r_B)` in (5)-(6) by `(L,S)`.

There is a completely rational THM-755 threshold. Put `q=333/106<pi`. If

```text
q*a*(6mu/7) > 1                                                (7)
```

and

```text
N > r_B*(q*2a/7+1)/(q*6a*mu/7-1),                             (8)
```

then the peel `v=aN` satisfies

```text
pi*(aN)*mu_N > r_N.
```

Alternatively, one marked safe interval gives the sufficient conditions

```text
q*a*(6L/7) > 1,
N > (q*2a/7+S)/(q*6aL/7-1).                                  (9)
```

In either case THM-755 closes `P_N union {aN}` at the `1/14` threshold.
Covering and primitivity are not analytic hypotheses. To call the result a
13-speed LRC14 row, additionally require `|B|=11` and check those arithmetic
properties.

## Proof

The `N`-runner's strict danger set is the union of `N` disjoint open teeth,
each of length `1/(7N)`, centred at `k/N`.

First decompose the positive-length part of `G'_B` into closed circle intervals
`J_i` of lengths `L_i`; there are at most `r_B` of them and
`sum_i L_i=mu`. A tooth can meet `J_i` only if its centre lies in the circular
`1/(14N)`-neighbourhood of `J_i`. That expanded interval has length
`L_i+1/(7N)`. A mesh-`1/N` grid meets it in fewer than `NL_i+2` points. Bounding
each intersection by a full tooth gives

```text
mu_N >= sum_i (L_i-(NL_i+2)/(7N))
     >= 6mu/7-2r_B/(7N),
```

which proves (1). Applying the same argument only to a named interval `J`
proves (3).

Remove the `N` open teeth from `G'_B` one at a time. Removing one connected
open arc from a finite union of closed arcs and points can increase its number
of components by at most one: it either trims, deletes, or splits one retained
component. Hence `r_N<=r_B+N`, proving (2). Independently, the union of all
danger arcs for `P_N` has at most `sum(P_N)=N+S` components, so its complement
has at most that many components; this proves (4). These bounds count singleton
boundary components as well as nondegenerate interval components.

Division gives (5), and the limit gives (6). For the full-state exact cap,
multiply (1) by `q*aN` and use (2):

```text
q*aN*mu_N-r_N
 >= (q*6a*mu/7-1)N-r_B*(q*2a/7+1).
```

Conditions (7)-(8) make this strictly positive. The same calculation with
(3)-(4) gives (9). Since `q<pi`, either inequality implies
`pi*aN*mu_N>r_N=r_top(P_N)`. Since `r_+(P_N)<=r_top(P_N)`, this stronger
inequality implies THM-755's strict capped-envelope criterion. ∎

## Iterated transport and the order gauge

The state bound composes. Let `B_0=B`, and successively adjoin distinct positive
frequencies `N_1,...,N_k` not already present. Put

```text
B_j  = B_{j-1} union {N_j},
mu_j = |G'_{B_j}|,
r_j  = r_top(B_j).
```

Repeated application of (1)-(2) gives

```text
r_j <= r_0 + sum_(h<=j) N_h,                                  (10)

mu_j >= (6/7)^j mu_0
        - (2/7) sum_(i=1)^j (6/7)^(j-i)
            * (r_0+sum_(h<i) N_h)/N_i.                        (11)
```

Indeed, (10) is immediate by induction. Substituting its `j-1` case into the
one-step mass bound and iterating proves (11).

The final set `B_k` is independent of the insertion order, but the certified
lower bound (11) need not be. Therefore every permutation gives a valid lower
bound for the same `mu_k`, and their maximum is valid as well.

This maximum has a canonical order. On an enclosure state `(m,r)`, write

```text
T_N(m,r)=(6m/7-2r/(7N),r+N).
```

For `x<y`, direct subtraction of the mass coordinates gives

```text
[T_y(T_x(m,r))]_mass - [T_x(T_y(m,r))]_mass
  = (2/7)*(y-x)/(xy)*(x+y+r/7) > 0.                 (12)
```

The component coordinate is the same after either two-step order, so every
later raw transition preserves this strict improvement, multiplied by a positive
power of `6/7`. An adjacent-exchange argument therefore shows that the affine
enclosure (11) is uniquely maximized by inserting the distinct frequencies in
increasing order. This does not claim that exact safe mass, or a hybrid bound
maximized against THM-780's global floor, has the same unique order.
In a fixed-core four-far chart, a one-peel capped-envelope search starts with 24
conceivable gauges,

```text
4 choices of peel * 3! insertion orders = 24
```

but (12) reduces them to four canonical state certificates, one per peel. This
is a finite proof calculus, not a finite classification of the chart: all four
lower bounds may be nonpositive, and the compression deliberately forgets
endpoint owners and correlations.

## THM-792 specialization

For

```text
B={1,...,9,15,110},  L=1/1540,  S=170,  a=1092,
```

the marked-window threshold (9) has exact crossing

```text
11734415/9278 < 1265.
```

Thus THM-793 supplies the elementary infinite tail in THM-792. THM-792's
divisor-packet statement, fragmentation lower bound, and 176-prime finite base
remain separate exact inputs.

## Uniform consequences with the THM-780 phase floor

The state theorem combines with the strict-margin floor in two useful ways.

First, let `B` be any eleven-speed base and put `H=max(B)`. Settled LRC for
twelve runners gives a `1/12`-deep time. Applying THM-780 with `d=11`,
`beta=1/12`, and `alpha=1/14` gives

```text
|G'_B|>=84^(-11).
```

Also `r_top(B)<=sum(B)<=11H`. Therefore (1) gives the uniform varying-base
bound

```text
|G'_(B union {N})|
  >= 6/(7*84^11)-22H/(7N).                                  (13)
```

In particular, along any sequence with `N/H -> infinity`,

```text
liminf |G'_(B union {N})| >= 6/(7*84^11)>0.                  (14)
```

Thus even a varying base cannot create mass decay on a genuinely transverse,
scale-separated ray. Comparable-scale phase incidence and the peel rate remain
outside this corollary.

Second, there is a much cleaner explicit terminal region in every literal
four-far chart. Let `C subset {1,...,14}` have nine speeds, put `S=sum(C)`, and
take ordered distinct far speeds `14<n_1<n_2<n_3<n_4`. If

```text
n_1 >= 412,
n_2 >= 412*(S+n_1),
n_3 >= 412*(S+n_1+n_2),
n_4 >= 412*(S+n_1+n_2+n_3),                              (15)
```

then `C union {n_1,n_2,n_3,n_4}` satisfies LRC14.

Indeed, settled LRC for ten runners gives `M(C)>=1/10`. The minimum-clearance
function is `max(C)<=14` Lipschitz, so `G'_C(1/14)` contains a closed circle
interval of length

```text
L >= 2*(1/10-1/14)/14 = 1/245.
```

Put `a=6/7`, `R_1=S+n_1`, `R_2=S+n_1+n_2`, and
`R_3=S+n_1+n_2+n_3`. Apply the marked-window bound to `n_1` and then the state
bound to `n_2,n_3`. For `P=C union {n_1,n_2,n_3}` this gives

```text
|G'_P| >= a^3*L-(2/7)*(a^2/n_1+a*R_1/n_2+R_2/n_3)
         >= A-B/412,                                      (16)

A=216/(343*245),       B=254/343,
r_top(P)<=R_3.
```

With `q=333/106`, exact arithmetic gives

```text
q*(412*A-B)=4455873/4453855>1.                            (17)
```

For comparison `q*(411*A-B)=4419909/4453855<1`, so `412` is the least integer
constant certified by this uniform calculation.

The final condition in (15) therefore implies
`q*n_4*|G'_P|>R_3>=r_+(P)`. Since `q<pi`, THM-755 closes the row by peeling
`n_4`.

Equations (15)-(17) prove a uniform `412`-lacunary terminal cone in every one
of the 2,002 four-far charts. Covering and primitivity are not required. This
is an exact certificate face in the literal four coordinates, not a finite
sample or a compactness slogan.

There is a complementary finite-exact cone with a far smaller repeated factor.
Exact rational enumeration of all `binom(14,9)=2002` small cores gives

```text
min_(C subset {1,...,14}, |C|=9) |G'_C|
  =10601/114660,
```

uniquely at `C={1,2,3,5,7,8,9,11,13}`. Therefore every four-far row closes if

```text
n_1>=19*S,
n_2>=19*(S+n_1),
n_3>=19*(S+n_1+n_2),
n_4>=19*(S+n_1+n_2+n_3).                                  (17a)
```

Indeed, use the full-state recurrence three times. With
`mu_*=10601/114660`, the final cap comparison is

```text
q*(19*(6/7)^3*mu_*-254/343)
  =66520746/57900115>1.                                    (17b)
```

The same calculation at factor `18` gives
`55930347/57900115<1`, so `19` is least for this uniform finite-core floor.
This factor-19 cone and the analytic factor-412 cone are complementary: the
former asks more of the first coordinate relative to `S`, but vastly less of
the subsequent separations.

### Every fully lacunary far-count flag is terminal

The same calculation closes all higher far-count strata and shows that the
large constant is confined to small bounded cores.

Let `4<=f<=12`, let `C subset {1,...,14}` have `13-f` speeds, put `S=sum(C)`,
and let `14<n_1<...<n_f`. For `f=4,5,6`, define

| `f` | `R_f` |
|---:|---:|
| 4 | 412 |
| 5 | 405 |
| 6 | 394 |

If

```text
n_1>=R_f,
n_j>=R_f*(S+sum_(h<j)n_h),       2<=j<=f,                  (18)
```

then the row closes by peeling `n_f`. To see this, settled LRC on the
`13-f`-speed core and the `14`-Lipschitz bound give a marked safe interval of
length

```text
L_f = 2*(1/(14-f)-1/14)/14 = f/(98*(14-f)).               (19)
```

With `m=f-1` and `a=6/7`, insertion of `n_1,...,n_m` leaves mass at least

```text
a^m*L_f-(2/(7R_f))*sum_(j=0)^(m-1) a^j.                  (20)
```

The table entries are the least integers for which

```text
q*(R_f*a^m*L_f-(2/7)*sum_(j=0)^(m-1)a^j)>1,              (21)
```

so the final inequality in (18) makes THM-755 fire.

For `7<=f<=12`, the core has at most six speeds. The union bound is much
stronger than one marked interval:

```text
|G'_C|>=1-(13-f)/7=(f-6)/7,       r_top(C)<=S.
```

Use the same cumulative conditions, now including the initial component load,

```text
n_1>=R_f*S,
n_j>=R_f*(S+sum_(h<j)n_h),       2<=j<=f,                 (22)
```

with the least certified integers

| `f` | 7 | 8 | 9 | 10 | 11 | 12 |
|---:|---:|---:|---:|---:|---:|---:|
| `R_f` | 27 | 17 | 14 | 13 | 13 | 13 |

The proof is (11) with `mu_0=(f-6)/7`, followed by the same `q`-cap
comparison.

Finally, if all thirteen speeds are far, start with `{n_1}`. Its safe set has
mass `6/7` and `r_top=n_1`. If

```text
n_j>=13*sum_(h<j)n_h,       2<=j<=13,                     (23)
```

insert `n_2,...,n_12` and peel `n_13`. The exact comparison is

```text
q*(13*(6/7)^12-(2/7)*sum_(j=0)^10(6/7)^j)
  =948176665875/733588221653>1.                            (24)
```

Thus fully lacunary infinity is terminal in every far-count stratum. Any
unresolved unbounded family must contain a bounded-ratio cluster at some level
of its cumulative scale flag. The remaining geometry is a recursive cluster
tree, not an arbitrary escape to infinity.

## Structural meaning

A high frequency can create `Theta(N)` endpoint walls inside an old good set,
so raw fragmentation can diverge with no divisor packet. But its danger duty
cycle removes asymptotically only one seventh of each retained interval, while
a proportional peel supplies `Theta(N)` sampling resolution. The minimal state
used by this proof operation is

```text
(safe mass, safe component count, new wall rate, named peel rate).
```

No single coordinate suffices for the insertion/cap predicate. THM-780 controls
safe mass uniformly, but not peel-relative jump load or phase incidence. The exact
theorem-facing cap load is `r_+(P_N)/(aN mu_N)`; the state above certifies it by
the conservative upper load `r_N/(aN mu_N)`.

This loss is exact, even if one adds every component length. Put

```text
B={1,2,3,6,12},       C={1,2,4,6,12}.
```

Both bases have `mu=4/7`, `r_+=r_top=12`, no isolated points, and the identical
component-length multiset

```text
{1/168 x2, 1/24 x2, 3/56 x4, 11/168 x4}.
```

After adjoining the same frequency `59`, however, exact endpoint arithmetic
gives

```text
|G'_(B union {59})|=2425/4956,
|G'_(C union {59})|=41/84=2419/4956,
```

while both new sets have `r_+=r_top=44`. Thus mass, component count, and even
the whole length multiset do not factor exact insertion. A lossless recursive
state must retain phase placement relative to the pending tooth grid, for
example the owner-coloured endpoint word or its `N`-pushforward incidence
profile.

## Assumption challenge and Tournament Analysis

The useful vertices are base safe components, new danger teeth, and the named
peel obligation, not runners. The exact carrier is the bipartite
component-tooth incidence with lengths and peel rate as sidecars. For the
upper bounds above, it admits the sound enclosure `(|G'_B|,r_B,N,a)`, while
losing exact transition data as the canary above proves.
For several insertions, the remaining-frequency multiset is also retained;
insertion order is a proof gauge rather than object data, and (12) fixes its
canonical increasing representative.

A runner-pair tournament destroys how many teeth land in one component and
cannot reconstruct (1). A diagnostic tournament may orient candidate carriers
by predicate retention, with ties resolved

```text
component-tooth incidence
  -> (safe mass,component count,wall rate,peel rate)
  -> owner-coloured endpoint word
  -> peel-relative scalar
  -> divisor-support profile
  -> raw runner tournament.
```

This Hamiltonian path is a declared gauge, not the proof. The proof is the
metric incidence inequality.

## Scope

THM-793 proves that positive safe mass persists under one high-frequency
insertion over a fixed base, and proves eventual closure only when the
peel rate is large enough relative to it. The one-step inequalities (1)-(2)
alone do not give a uniform lower bound over arbitrary cores, classify
sublinear peels, or prove the global HYP-6830 splice. THM-780 supplies the
uniform floors used above and rules out absolute safe-mass decay for
LRC14-relevant bases of at most twelve speeds. Varying-base and multiscale
directions can still obstruct a practical sharp bound, a
sublinear/nonproportional peel, or owner/event transport erased by this
two-number state.
