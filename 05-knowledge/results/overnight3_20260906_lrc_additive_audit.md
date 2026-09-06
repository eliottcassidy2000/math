# Additive mixed-parity tails have a larger physical bulk

**Status: PROVED ANALYTICALLY + INDEPENDENT ROOT AUDIT PASS, with exact
finite controls.** Let `a<b`, `c=a+b`, and assume `a,b,c` are positive
integers, primitive as a triple, and nonzero modulo three. No oddness
assumption is imposed. For the literal scale-three distinct-sheet physical
comb defined below,

```text
|mu(F_(a,b,c))-9/98| <= 6/(7c).                       (1)
```

In particular every such additive triple with `c>=62` satisfies

```text
mu(F_(a,b,c))>6/77.                                  (2)
```

Thus the failure of extending the odd ternary-unit theorem to mixed parity
is a cofinite physical bulk phenomenon, not just a small isolated example.
This does **not** say that only additive triples can exceed the target:
the nonadditive `(2,11,20)` has mass `11/140>6/77`.

## 1. Parity-free physical address, derived from the sheets

For positive ternary-unit integer speeds `w=(w_1,w_2,w_3)`, define the
literal sheet sets by

```text
D_(w_i,s)={x mod1: ||w_i(x+s/3)||<1/14}, s in {0,1,2},
F_w=union_(pi in S_3) intersection_i D_(w_i,pi(i)).
```

The sheet assignment is unique wherever it exists, since the three sheet
families for one speed are disjoint. Put `y=3x` and `r=3/14`. Membership
means that there is an integer vector `n` with

```text
|w_i y-n_i|<r,
s_i=-w_i^(-1)n_i mod3,
{s_1,s_2,s_3}=F_3.
```

Shifting `y` by one shifts `n` by `w` and subtracts one from every sheet
label, preserving distinctness. The raw condition therefore has period one
in `y`. The three inverse x-branches each scale length by one third, so the
literal mass equals its raw y-mass over one period. No oddness is used in
this change of variables.

For primitive `w`, the map

```text
n mod Zw -> C=w cross n
```

is a bijection onto `L={C in Z^3:C dot w=0}`. Choose an integer Bezout
vector `z` with `z dot w=1`; then `n=C cross z` is an explicit inverse.
The kernel is `Zw`, since a rational multiple of primitive `w` that is
integral has integral multiplier.

Modulo three,

```text
C_i=w_j w_k(s_j-s_k).
```

All speeds are units, so the distinct-sheet condition is exactly
`C_i!=0 mod3` for every coordinate. For a fixed lift `n`, the three real
phase intervals have radii `r/w_i`. They have a common positive-length
intersection exactly when their three pairwise strict roof inequalities
hold. Their intersection length is

```text
L_w(C)=max(0,min( 2r/w_1,2r/w_2,2r/w_3,
                 r/w_1+r/w_2-|C_3|/(w_1w_2),
                 r/w_1+r/w_3-|C_2|/(w_1w_3),
                 r/w_2+r/w_3-|C_1|/(w_2w_3) )).
```

Different nearest-integer lifts have disjoint interiors because `r<1/2`.
Consequently the physical mass is exactly the sum of `L_w(C)` over the
full owner-live carrier dictionary. This paragraph rederives the address
and residue facts for the actual mixed-parity sheet definition; it does
not import an odd-only theorem outside its hypotheses.

## 2. Complete additive dictionary and its exact tent

For `w=(a,b,a+b)`, use the primitive relation `v=(1,1,-1)`. If a phase
is physically eligible and `e=n-yw`, then

```text
delta=v dot n=v dot e is an integer,
|delta|<r||v||_1=9/14<1.
```

Thus `delta=0`, even before invoking the owner residues. The identity
`v cross C=delta w` forces `C=k v`, with integer `k` by primitivity of
`v`. The three roof bounds reduce to the shortest pair sum, `a+b=c`.
The complete dictionary is therefore

```text
Lambda(w)={k(1,1,-1): k in Z, 0<|k|<3c/14, 3 does not divide k}.         (3)
```

The unit-speed condition also shows `a=b mod3`; their sum must not vanish
modulo three. In particular, the weighted relation residues are all equal.
There is no missing owner fiber or multiplier in (3).

Let `f(t)=L_w(t v)` for real `t`, extended evenly and by zero outside its
support. Put `q=3/(7c)`. For `t>=0`, direct comparison of the three pair
roofs gives

```text
f(t)=q,                                 0<=t<=3a/14;
f(t)=[3(a+2b)-14t]/(14bc),               3a/14<=t<=3b/14;
f(t)=[3c-14t]/(14ab),                    3b/14<=t<=3c/14;
f(t)=0,                                 t>=3c/14.                    (4)
```

For clarity, call the middle expression `h_1`, the pair roof omitting the
smallest speed. Throughout the live interval,

```text
h_1-q=(3a-14t)/(14bc),
h_2-h_1=(b-a)(3c-14t)/(14abc)>=0,
h_3-h_1=(3b-14t)/(14ac).
```

These identities give exactly the two switches in (4), and the values
match at both switches. The function is continuous, even, nonnegative,
and decreasing for positive `t`. At `t=3c/14` its value is zero, agreeing
with the strict endpoint exclusion.

## 3. Exact integral and the cofinite obstruction

The central rectangle, middle trapezoid, and final triangle in (4) give

```text
integral_0^infinity f(t)dt
 = q(3a/14)
   +[q+3/(14b)](3(b-a)/28)
   +[3/(14b)](3a/28)
 =27/392.
```

Hence the full integral is **`I=27/196`**, independent of `a,b`. By (3),

```text
mu(F_w)=sum_(k in Z)f(k)-sum_(k in Z)f(3k).            (5)
```

The term at zero cancels; it is not a live carrier. For an even decreasing
function, the elementary rectangle comparison proves

```text
|sum_(k in Z)f(hk)-I/h|<=f(0).
```

Apply it at `h=1` and `h=3` to (5). Since `f(0)=q`, this gives

```text
|mu(F_w)-(2/3)I|<=2q,
|mu(F_w)-9/98|<=6/(7c),
```

proving (1). The target gap is

```text
9/98-6/77=15/1078.
```

The lower estimate is strictly above the target whenever
`c>2156/35`, hence at every integer `c>=62`. The estimate is uniform over
the entire additive ratio range; it is not confined to `a=1` or to a fixed
ratio as heights grow. It proves convergence to `9/98` along every
unbounded eligible additive sequence.

Oddness blocks this mechanism at its source: three odd speeds cannot
satisfy `a+b=c`. In the universal odd coefficient proof, the same fact
appears as even relation norm, excluding the norm-three pattern `(1,1,1)`.
Its normalized slice-count slope is `2/7`, rather than the nonexceptional
`15/98` ceiling. This is a precise parity boundary, not evidence that the
audited odd theorem fails.

## 4. Exact controls and the next boundary

The parity-free integer-address verifier in the
[universal audit companion](../../04-computation/overnight3_20260906_lrc_universal_audit.py)
independently reconstructs these dictionaries and masses:

| Speeds | Physical mass | Raw carriers |
|---|---:|---:|
| `(2,5,7)` | `22/245` | 2 |
| `(1,7,8)` | `31/392` | 2 |
| `(1,10,11)` | `6/55` | 4 |
| `(1,61,62)` | `2467/26474` | 18 |
| `(1,121,122)` | `9599/103334` | 36 |

All exceed `6/77`; the last two directly check the cofinite conclusion.
The integral identity and the quadrature error are checked separately.
The theorem makes no assertion that every smaller additive triple exceeds
the target: for example `(1,4,5)` has only `k=+-1`, each of length `1/56`,
so its mass is `1/28<6/77`.

The tempting repair that **only additive triples exceed `6/77` is false**.
At `(2,11,20)`, the relation `v=(1,-2,1)` has norm four; the same direct
error estimate gives `|v dot n|<6/7<1`, so all carriers lie on its line.
Its complete live list is `+-v,+-2v`. The individual physical lengths are

```text
L(v)=L(-v)=3/140,
L(2v)=L(-2v)=1/56,
mu(F_(2,11,20))=11/140>6/77.
```

Thus the odd-speed norm-four exception also needs a new analysis after
dropping parity. No classification of every mixed-parity relation sector,
no sharp maximum over all additive triples, and no entry theorem is claimed
here. The parent's separate `(1,b,b+1)` investigation may sharpen that
family; it is not a dependency of the broader bulk law (1).

The connection contract is exact: the literal sheet union maps to its
nearest-integer error lift, then to the complete norm-three carrier line.
It preserves physical interval length and the mod-three deletion word.
The integral gives a uniform bulk value while the central height bounds
the sampling error. Forgetting parity admits this new relation, and
forgetting the mod-three word would give the wrong bulk constant.

Reproduction and frozen source/output hashes are recorded in the
[universal independent audit](overnight3_20260906_lrc_universal_audit.md).
This result remains a local physical-comb obstruction, not a counterexample
to LRC(14).

## Independent root review

Root independently rederived the nearest-integer address, Bezout inverse,
mod-three sheet residues, interval intersection formula, zero-defect line,
all three tent pieces, and the integral `27/196`. The two rectangle-rule
errors are both retained in the difference of sampled sums; the resulting
strict cutoff is exactly `c>2156/35`. The
[separate literal consumer companion](../../04-computation/overnight3_20260906_lrc_consumers.py)
reconstructs every eligible primitive additive pair with `1<=a<30` and
`a<b<60` from six physical sheet intersections and agrees with the tent.
Its independent adjacent-family formula also agrees for every `b=4,7,...,1000`.
That second derivation provides the sharp adjacent-family result without
being a dependency of the broad estimate. No odd-only carrier premise is
used in either physical derivation.
