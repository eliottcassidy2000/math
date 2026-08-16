# The visible U_full cyclic-component and point diagonals miss the frozen endpoint bridge

**Status: FINITE-EXACT bridge-only hostile and API/kernel result.**  This note
tests only the bridge between `q_H=(1,0,1)` and `q_q5=(1,0,0)`.  It does not
extend promoted THM-3479 to a physical current, construct a THM-2471 ancestry stalk for `U_full`,
supply a chronological arrival, evaluate either `K4`, remove a row, or prove
LRC(14).

## 1. Inheritance and the exact question

The inherited refined endpoint bank gives, in the certified split-prime chart,

```text
p=572252886246508880869,

A(q_H)  =320618948602619577408,
A(q_q5) =503604956476841920373,

A(q_H)-A(q_q5)=389266878372286537904 mod p.          (1)
```

The bridge in (1), and both later `K4` factors, are nonzero in all `72` role
charts.  Those are endpoint-aggregate graph facts.  The immediately preceding
reflection showed that an abstract five-atom positive realization is easy but
that the endpoint APIs return only separately integrated factors.

The inheritance pass is

```text
closest proved mechanism: THM-2471's linked-node Boolean ancestry stalk;
positive model:           THM-2594's one-base, before-integration contraction;
canonical hostile:        THM-2538's mixed-Haar transportation kernel;
corrected near miss:      MISTAKE-293, common ancestry is not one circle point;
least-used sidecar:       THM-2471's root/sheet/horizon coordinates (28)--(31).
                                                                    (2)
```

The cheapest exact question was whether a joint relation already visible in
the current interval geometry recovers (1).  There is no reason to evaluate
either `K4` unless this bridge gate passes.

## 2. Three geometric products, with the cyclic validity repair

Let `E_ell` be the shifted present Boolean set and `Q` the fixed delayed-word
set.  On the scaled circle the current engine computes

```text
AX_ell = sum_C AX_(ell,C),
BY_ell = sum_D BY_(ell,D),

gamma_ell = phase_ell (sum_C AX_(ell,C))(sum_D BY_(ell,D)),       (3)
```

where `C,D` range over maximal cyclic components of `E_ell`; the left factor
also intersects `Q`.  The endpoint frequencies are

```text
X=13,  Y=X+w_c,  w_c=742586.                                    (4)
```

The stored API exposes only the totals in (3), but its internal geometry
permits two cheap restrictions.

First, use equality of maximal cyclic `E` components:

```text
K_ell = phase_ell sum_C AX_(ell,C) BY_(ell,C).                    (5)
```

This is the unweighted/counting same-component restriction.  It uses the same
two endpoint factors and hence the same global two-factor unit as (3).  To make
(5) invariant under implementation partition and choice of circle cut, the
probe first merges adjacent linear segments and then glues the first and last
segments when they meet across `0~T`.  The full frozen universe happens to
contain neither adjacent stored segments nor origin-crossing components:

```text
raw segment count     =71,070,080,
maximal segment count =71,070,080,
cyclic origin glues   =0.                                         (6)
```

Thus the stored-segment and lawful cyclic-component banks coincide here.  The
algorithm nevertheless gates the two failure modes rather than assuming them
away.

Second, identify the two circle points.  Since the endpoint exponents are
`-X y` and `+Y y`, their product has combined frequency

```text
X-Y=-w_c=-742586.                                                  (7)
```

The resulting point bank is

```text
D_ell=phase_ell sum_[a,b] in components(R E_ell intersect Q)
        (g^(w_c a)-g^(w_c b)).                                    (8)
```

By MISTAKE-293, neither (5) nor (8) is a THM-2471 ancestry stalk.  Linked
nodes may share a base without occupying one component or one circle point.

## 3. Exact component verdict and cross-component mechanism

The full normalized inverse `F_13^3` transform of (5) gives

```text
K^(vee)(q_H)  = 99143203042879994518,
K^(vee)(q_q5) =130742602587392835137,

K^(vee)(q_H)-K^(vee)(q_q5)
 =540653486701996040250 mod p.                                    (9)
```

Its complete character-bank digest is

```text
1270078f3c019a3fa5ab507128e5010149ba93bbb6f74d4c6967cc2753df21ea. (10)
```

Therefore the unweighted same-cyclic-component restriction does **not**
recover (1).  More usefully, (3) has the exact decomposition

```text
gamma_ell = K_ell + R_ell,

R_ell = phase_ell sum_(C != D) AX_(ell,C) BY_(ell,D).              (11)
```

The inverse transform of the cross-component remainder is

```text
R^(vee)(q_H)  =221475745559739582890,
R^(vee)(q_q5) =372862353889449085236,

R^(vee)(q_H)-R^(vee)(q_q5)
 =420866277916799378523 mod p,                                    (12)
```

with bank digest

```text
529cf77d326ada80ab99bfbf4b7b6444d41794820f1fc4552555846a5f49c978. (13)
```

Equations (9) and (12) sum to (1) modulo `p`.  Thus the entire failure of the
named component diagonal is exactly the nonzero `C!=D` bank, in the inherited
two-factor normalization.  This is a localization statement, not a semantic
interpretation: it does not say that cross-component pairs are physical
ancestry.  Nor does it rule out component-dependent coupling weights such as
mass-normalized conditionals.

## 4. The finer point restriction also misses, with a normalization caveat

The inverse transform of (8) gives

```text
D^(vee)(q_H)  =   633668780131603861,
D^(vee)(q_q5) =405160484437854840264,

D^(vee)(q_H)-D^(vee)(q_q5)
 =167726070588785644466 mod p.                                    (14)
```

Its character-bank digest is

```text
771545a5cb1f0f03459b8d351de668ad950ece5fcb985fa61d599d643de3303f. (15)
```

The Cartesian bank has two omitted endpoint denominators while (8) has one,
so (14) is comparable to (1) only in the pinned raw endpoint-jump convention.
The exact discrete geometric-series restoration gives

```text
((1-g^(-X))(1-g^Y))/(1-g^w_c)
 =123588610788991450223,

scaled point bridge=24876622649422736677,                          (16)
```

which still misses (1).  For orientation only, the unsigned rational
frequency-magnitude probe `XY/w_c` gives

```text
scale=572122651179088191560,
scaled point bridge=284514977757516176864.                         (17)
```

Equation (17) is not a literal Fourier restoration: sign and phase depend on
the endpoint convention, and a finite split field has no canonical
Archimedean `i`.  An arbitrary fitted scalar could map any nonzero bridge to
another, so (14)--(17) do not exclude a rescaled or differently typed point
lift.  The component result (9)--(13) avoids this global denominator mismatch,
but still fixes only the unweighted equality relation on components.

## 5. The precise API/kernel obstruction

The actual signatures are

```text
fast_x_sweep(...)       -> (tuple(totals), overlap),
fast_endpoint_sum(...)  -> tuple(totals).                           (18)
```

Neither accepts or returns an atom, ancestry, stalk, joint-key, horizon, or
address field.  The left loop internally sees an `E` component, a delayed-word
component, and a wrap branch.  The right loop sees an `E` component.  All such
labels are discarded before the refined worker forms
`phase*AX_total*BY_total`.  The current returned scalar API therefore cannot
even recover (5), although the internal geometry can compute it.

After full atom spaces `U,V`, their augmentations, and a common corner have
actually been defined, the full row/column marginal map has mixed kernel

```text
ker(epsilon_U) tensor ker(epsilon_V),                              (19)
```

of dimension `(m-1)(n-1)`; this is the pinned proved THM-2538 mechanism.  Its
smallest direction is the checkerboard

```text
[[ 1,-1],
 [-1, 1]],                                                         (20)
```

whose row and column sums vanish while its diagonal functional is `2`.

The actual API does not literally land on a full marginal fibre product: no
such atom spaces or common corner have been defined.  Under any compatible
atomization its two character-weighted scalar returns are a further linear
image of the factor marginals, so they lose at least (19), and possibly more.
We do not claim that the Cartesian and diagonal banks are two joint lifts with
identical full marginals.  The exact computation rules out the named geometric
restrictions; (19)--(20) separately prove nonidentifiability of a lawful joint
lift from the returned scalars.

## 6. Cheapest lawful engine change

A future coupling must be one object before character evaluation, not a
different relation chosen separately for each `ell`.  The minimum typed record
is schematically

```text
Atom(omega, measure, left_factor, right_factor),

omega=(base,root,owner_sheet,word_sheet,source_sheet,
       left_horizon,right_horizon,address),                        (21)

r: Omega -> F_13^3 quotient.                                      (22)
```

It must carry projections from the common stalk to the linked left and right
nodes, the THM-2471 root/sheet/horizon data, and one `ell`-independent exact
address map `r`.  The whole character bank must then arise covariantly by
twisting the same records with `chi_ell(r(omega))`, multiplying the linked
factors before summation, and only afterward applying inverse DFT.  Per-`ell`
keyed records alone need not be the Fourier transform of one physical
coupling.

The cheapest next engine work is therefore:

1. retain factor-labelled cells before Boolean interval merging;
2. attach the common-stalk projections, ancestry sheets, horizons, and (22);
3. state an explicit `U_full` support relation;
4. multiply linked endpoint factors before summing; and
5. test only recovery of (1), before either `K4`.

THM-2594 is the positive implementation model because its factors live at
linked nodes on one base and are multiplied inside one integral.  It is not a
transplant to `U_full` and supplies none of the missing fields in (21)--(22).

## 7. Controls, connection contract, and strict boundary

The Cartesian positive control recomputes the entire inherited bank and
recovers its digest, both inverse values, and bridge exactly:

```text
1fabc5cfdbaa1455e10cd6bf9264488133616a7b0ff381623d729b4b4bfa9682. (23)
```

The optimized periodic sweep and independent boundary-refinement engine agree
at `ell=0` and the inherited hostile representative `v2`.  At those controls,
the point sums are respectively

```text
24175413190668994534,
486792831500322935265,                                             (24)
```

and independent public per-component aggregation agrees with (5).  The full
run covers all `2197` characters, not representatives or samples.

| field | exact content |
|---|---|
| source | frozen `U_full` `13^3` endpoint product bank |
| attempted target | atom-level `H-q5` common-ancestry bridge |
| first map | unweighted equality of maximal cyclic `E` components |
| component verdict | bridge `540653486701996040250`, not (1) |
| exact complement | `C!=D` bridge `420866277916799378523` |
| finer hostile | common-point bridge `167726070588785644466` in the raw convention |
| destroyed information | component values, factor labels, linked ancestry, horizons, common address |
| needed sidecar | one THM-2471-style stalk with projections and `ell`-independent map (22) |
| next decisive test | keyed bridge only; no `K4` until exact recovery of (1) |

Nothing here says that no lawful `U_full` ancestry coupling exists.  Nothing
here proves a grouped coefficient `C(a;X,m)`, all-unit projector `B(q)`,
physical current, chronological owner loop, scalar-row exclusion, `U_clock`
statement, THM-3479 promotion, or LRC(14).

## 8. Reproduction and replay scope

Run from the repository root:

```text
PYTHONHASHSEED=0 python3 -B \
  04-computation/lrc_endpoint_ufull_atom_bridge_kernel_probe_20260816.py
```

One full normal exact sweep generated the frozen banks.  Because that sweep
processed `71,070,080` components, a second full `-O` sweep was deliberately
not run.  After freezing, the script passed `py_compile`, an AST audit with
`assert_nodes=0` and no forbidden dynamic calls, and the two independent cheap
reference/API controls.  The separate diagonal audit checked the component
definition, cyclic-cut repair, normalization scope, kernel typing, and
consequence boundary.

The final sync promoted the inherited THM-3479 package and changed only
docstrings, status text, and dependency/semantic metadata in the three pinned
endpoint artifacts; the computational functions and frozen endpoint bank did
not change.  The dependency-sensitive digest was recomputed from the frozen
banks after first reconstructing the old digest exactly as a positive control.

LF hashes are

```text
script:   d1182de4d777bab20a8d423cf942151ac3149014b67d9c34883cbce37a7b0a9f
output:   a01826e37a674476d76089e2425cce45d02fb1c79134f1bf5a7430c2cef66a2a
semantic: b28f64e9beac1a9677d6562b606983b35824f54d9d2dc812d8678f8dbea47948
```
