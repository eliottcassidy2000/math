---
id: THM-2456
title: "Two-root replica uniform-offset boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For an
  invertible nonnegative cyclic convolution,
  the anchored delta-plus-replicas locus is exactly a uniform
  owner-minus-replica density offset. One Boolean root chart leaves
  only the trivial full-owner/empty-replica table, but two rational
  charts realize a nonconstant hostile whenever the convolution
  kernel has a missing support point. The actual seven source phases
  impose 7f+c<=1, yet an eight-chart model is sharp; singleton-only
  and paired-safe controls use 20 and 104 charts. Both ordinary
  danger and safe kernels lie on this hostile side. This is a finite
  coefficient boundary, not a canonical LRC packet or row exclusion.
source: codex-2026-07-26-two-root-replica-offset
depends_on:
  - THM-2449-coprime-owner-anova-and-delta-replica-boundary
related:
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2401-common-filter-endpoint-or-first-death-certificate
  - THM-2449-coprime-owner-anova-and-delta-replica-boundary
  - THM-2452-indicator-idempotent-aggregate-endpoint-restoration
  - THM-2457-complete-atom-root-cosupport-graph-and-semantic-word-hostile
script: 04-computation/lrc14_two_root_replica_uniform_offset_thm2456.py
output: 05-knowledge/results/lrc14_two_root_replica_uniform_offset_thm2456.out
script_sha256: 118d21f6d6218bd059aa7ee47bfacfb72f5138adc31662fa9240dce4e1f1059f
output_sha256: c7e5e9333a3e0fc900667f86359aa0361def44bbacb15ebf494dfb126bd55a35
hash_basis: working-tree bytes (LF)
---

# THM-2456 -- two-root replicas are uniform offsets before mixing

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2449 proves that the only mixed-zero owner table is one positive
source delta plus six identical target replicas. THM-2452, on the
other hand, supplies complete Boolean endpoint atoms only after
integration. This theorem identifies the exact gap between those two
facts:

```text
Boolean on one common root chart       -> replica branch is trivial;

Boolean before integration over charts -> a convex uniform-offset
                                           hostile remains.          (1)
```

The distinction is load-bearing. Complete truth bits do not turn an
integrated additive identity into a pointwise additive identity.

## 1. Invertible cyclic-kernel normal form

Let `G=C_p` with `p>1`, let

```text
h:G->Q_(>=0),                 kappa=sum_G h>0,
```

and use the cyclic convolution

```text
(Kf)(s)=sum_(r in G) h(r-s)f(r).                    (2)
```

Assume that `K` is invertible over `Q`; equivalently, every cyclic
Fourier eigenvalue of `h` is nonzero. Let

```text
0<=f_ell(r)<=1,               ell in C_7, r in G,

A_ell=K f_ell.                                      (3)
```

The `f_ell` are averaged root occupancies. They need not themselves be
Boolean.

Suppose that the untwisted target column is the positive owner anchor

```text
A_ell(0)=a 1_(ell=0),                  a>0,          (4)
```

and every THM-2449 rectangle vanishes. Then there is a nonnegative
target row `w`, with `w(0)=0`, such that

```text
A_0=a 1+w,

A_ell=w                         for ell!=0.          (5)
```

Put

```text
c=a/kappa.
```

Since `K1=kappa 1`, invertibility of `K` turns (5) into the unique
pre-gate density normal form

```text
f_1=...=f_6=f,

f_0=f+c 1,

Kf(0)=0,

0<=f<=1-c.                                           (6)
```

Conversely, every `f,c` satisfying (6) gives (4)--(5), with

```text
a=kappa c,                    w=Kf.                  (7)
```

Thus the mixed-zero condition forgets exactly one coordinate: whether
the uniform offset in (6) occurs on every common chart or only after
averaging different charts.

Because `h,f` are nonnegative, the anchor condition in (6) is
pointwise:

```text
f(r)=0 whenever h(r)>0.                              (8)
```

## 2. Sharp support boundary and the two-chart hostile

Within the hypotheses of Section 1, a nonconstant averaged hostile
exists if and only if

```text
supp(h) is a proper subset of G.                    (9)
```

If `h` is positive on all of `G`, equation (8) forces `f=0`, so
`w=0`. Conversely, suppose `h(r_0)=0`. On two charts of equal mass
take the six replica masks to be

```text
chart I:  {r_0},
chart II: empty,
```

and the owner masks to be

```text
chart I:  G,
chart II: {r_0}.                                    (10)
```

Their averaged occupancies are

```text
f=(1/2)1_{r_0},

f_0=f+(1/2)1.                                       (11)
```

Hence

```text
a=kappa/2,

w=(1/2)K1_{r_0},

w(0)=0.                                             (12)
```

The row `w` is nonconstant: otherwise invertibility would send a
constant vector to the nonconstant singleton `1_{r_0}`. Therefore
(10) is an exact nonconstant delta-plus-six-replicas table made only
from Boolean charts.

For the ordinary two-root danger kernel on `C_13`,

```text
(Kf)(s)=f(s)+f(s+1),

r_0=2,                                              (13)
```

the table is

```text
a=1,

w_1=w_2=1/2,                  w_s=0 otherwise,

A_0=1+w.                                           (14)
```

Each replica row has total mass `1`, the owner row has total mass
`14`, and the full seven-row table has total mass `20`.

This two-chart control is not yet typed by the pointwise-exclusive
seven source dangers. Section 4 installs that missing constraint.

## 3. One Boolean chart is rigid

Suppose in addition that the normal-form occupancies themselves are
Boolean:

```text
f_0,f in {0,1}^G.                                  (15)
```

The positive constant difference

```text
f_0-f=c1
```

can then have only the value `c=1`. Consequently

```text
f=0,                  f_0=1,

a=kappa,              w=0.                         (16)
```

Thus a single common Boolean root chart has no nonconstant replica
hostile. The failure in Section 2 is purely convex: Booleanity is lost
when the common chart is forgotten.

## 4. Seven exclusive source rows: the sharp eight-chart model

In the lawful THM-2449 table the seven translated source dangers are
pointwise exclusive. If the Boolean source masks on every chart obey

```text
sum_(ell=0)^6 b_ell(r)<=1,                          (17)
```

then averaging (17) and using (6) gives the exact extra constraint

```text
7f(r)+c<=1                         for every r.      (18)
```

Conversely, for rational `f,c`, conditions (6) and (18) are sufficient
for an abstract finite source-exclusive Boolean realization. At each
root choose independently among the labels

```text
none,0,1,...,6
```

with probabilities

```text
1-7f-c,        f+c,        f,...,f;                 (19)
```

the finite product distribution is rational and decomposes into
finitely many Boolean charts.

For the two-root kernel (13), (18) is sharp on eight equal charts.
Put

```text
c=1/8,                    f=(1/8)1_2.               (20)
```

An explicit realization is:

```text
chart 0: owner on all thirteen roots;
chart 1: owner only at root 2;
chart ell+1: replica ell only at root 2,  1<=ell<=6. (21)
```

At every chart/root pair at most one source row is present. The table
has

```text
a=1/4,

w_1=w_2=1/8,                 w_s=0 otherwise,

A_0(1)=A_0(2)=3/8,

A_0(s)=1/4                   otherwise,

sum_(ell,s)A_ell(s)=5.                              (22)
```

Eight is minimal among equal-weight source-exclusive chart models
with a positive anchor and a nonconstant replica. Indeed, on `N`
equal charts, nonzero `f` and positive `c` are both at least `1/N` at
some root. Equation (18) then gives

```text
8/N<=1.
```

## 5. Singleton and paired-safe controls

The all-root owner chart in (21) is not responsible for the failure.
A `20`-chart source-exclusive model uses only singleton root masks:

```text
thirteen owner charts, one at each root;
one extra owner chart at root 2;
six replica charts at root 2.                       (23)
```

It has

```text
c=1/20,                   f=(1/20)1_2,

a=1/10,                   w_1=w_2=1/20,

sum_(ell,s)A_ell(s)=2.                              (24)
```

The count `20` is minimal among **equal-weight** models under the
additional rule that every chart contains at most one source/root
pair: the uniform offset costs thirteen occupied owner pairs and a
nonzero common replica costs one more owner pair plus six replica
pairs.

The equal-weight qualifier is load-bearing.  With arbitrary rational
chart weights, the same table (24) has a `19`-support realization:
twelve owner-singleton charts away from root `2` have weight `1/20`,
the owner-singleton chart at root `2` has weight `2/20`, and the six
replica-singleton charts at root `2` each have weight `1/20`.  Thus
`20` is not minimal by support cardinality among weighted models.

Nor does the other member of a lawful target dipole automatically
remove the boundary. Let

```text
g_t(s)=1_(s-t notin {5,6}),              t in C_13. (25)
```

For `t=0`, use the eight charts in (21). For each `t!=0`, use eight
charts with only chart zero carrying the owner on all roots. Average
the resulting `13*8=104` source-exclusive charts, and multiply the
response on the `t`-block by `g_t(s)`. Since

```text
sum_t g_t(s)=11,
```

the exact paired table is again delta plus replicas:

```text
a=11/52,

w_1=w_2=1/104,                w_s=0 otherwise,

A_0(1)=A_0(2)=23/104,

A_0(s)=11/52                  otherwise,

sum_(ell,s)A_ell(s)=75/26.                         (26)
```

This is a finite coefficient model of an ordinary two-root gate paired
with a translated two-hole safe gate. It is not asserted to be a
canonical speed-row packet.

## 6. Both ordinary truth values lie on the hostile side

Let `S` be a nontrivial cyclic shift on `C_13` and let `J` be the
all-ones convolution. The ordinary danger and safe kernels are

```text
K_d=I+S,

K_g=J-I-S.                                          (27)
```

For a thirteenth root `zeta` their eigenvalues are

```text
K_d:  2                         at frequency 0,
      1+zeta^k                  at k!=0;

K_g:  11                        at frequency 0,
      -(1+zeta^k)               at k!=0.            (28)
```

Since `13` is odd, `-1` is not a thirteenth root. Both kernels are
invertible. Their supports have sizes `2` and `11`, respectively, so
both satisfy the proper-support hostile criterion (9). The danger
inverse is explicit:

```text
(I+S)^(-1)
 =(1/2)(I-S+S^2-...+S^12).                         (29)
```

The same statements hold for `S^d`, `d!=0 mod 13`, after an affine
root reindexing.

## 7. Quantitative fixed-kernel stability

Use the normalized `ell^2(C_13)` norm. For either kernel in (27), and
for one selected replica row, put

```text
e=f_0-f-c1,

Delta(s)=A_0(s)-A_ell(s)-a=(Ke)(s).                 (30)
```

The least singular value of both `K_d` and `K_g` is

```text
sigma_min=2 sin(pi/26),                             (31)
```

attained at the two frequencies nearest `-1`. Therefore

```text
||Delta||_2^2
 >=4 sin^2(pi/26)||e||_2^2.                         (32)
```

For `K_d`, the Gram operator is

```text
K_d^*K_d=2I+S+S^(-1).
```

For `K_g`, it differs by `9J`; this changes only the constant
eigenvalue and leaves the minimum in (31) unchanged.

Equation (32) is scoped to one fixed cyclic kernel and one common
density comparison. It gives no denominator-free amplitude floor
from Boolean truth bits alone, and it does not survive an unrecorded
mixture of different root-status kernels.

## 8. LRC application and exact nonclosure

On a frozen common `C_13` root chart, an ordinary danger bit with two
roots has kernel `K_d` after an affine reindexing; its safe complement
has kernel `K_g`. If all other target-moving factors have been proved
root-constant or have been retained in the fixed chart, Sections
1--7 classify the THM-2449 replica branch exactly.

That conditional statement does not close the current LRC endpoint
problem:

1. THM-2452 supplies a complete matched Boolean atom after integration,
   not a pointwise-in-base-chart THM-2449 identity.
2. A lawful target dipole also moves a blocker factor. The exact
   `104`-chart model (25)--(26) shows that this extra translated safe
   gate is not, by itself, a rectangle certificate.
3. The prescribed seven comb bits and terminal word may forbid the
   displayed convex mixtures, but no current theorem proves that
   finite realizability exclusion.
4. A fibrewise same-chart rectangle, an extremal capacity identity
   forcing Boolean densities almost everywhere, or a separately
   proved semantic owner/root service would remove this boundary.

Thus complete local truth bits alone do not force a rectangle defect.
This theorem supplies the exact finite-state test to run on any future
physical root atlas: solve for the averaged occupancies in (6), add
the source constraint (18), and determine whether the actual atom
types realize a nonzero `f`. No semantic THM-2305 word, canonical
THM-2401 root service, scalar-row exclusion, or LRC(14) conclusion is
asserted.

## 9. Exact companion

The companion

```text
04-computation/lrc14_two_root_replica_uniform_offset_thm2456.py
```

uses only integers and `fractions.Fraction`, with explicit `require`
checks. It:

1. verifies exact rank `13` for all twelve danger and safe shifts;
2. checks the full-support invertible control;
3. checks the normal form and one-chart rigidity;
4. constructs the two-, eight-, twenty-, and `104`-chart tables;
5. checks every anchor, rectangle, source-exclusivity condition, row
   sum, and displayed constant;
6. verifies the exact Gram identities underlying (31)--(32).

Reproduce with

```bash
python3 04-computation/lrc14_two_root_replica_uniform_offset_thm2456.py
python3 -O 04-computation/lrc14_two_root_replica_uniform_offset_thm2456.py
```

Both transcripts must equal

```text
05-knowledge/results/lrc14_two_root_replica_uniform_offset_thm2456.out
```

byte for byte. The script intentionally contains no Python `assert`,
floating-point, numerical trigonometry, or external algebra package.

An independent hostile audit rederived the normal form, converse,
proper-support boundary, Boolean-chart rigidity, source-exclusive
realization, all four finite controls, both kernel spectra, and the
fixed-kernel stability constant.  It also caught and repaired the
convolution-sign wording in (8) and the equal-weight qualifier on the
`20`-chart minimum before promotion.  Normal, optimized, and stored
transcripts and both LF hashes were independently reproduced.
