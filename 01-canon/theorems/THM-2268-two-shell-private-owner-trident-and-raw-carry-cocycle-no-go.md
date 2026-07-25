---
id: THM-2268
title: "Two-shell private-owner trident and raw carry-cocycle stopping boundary"
status: >
  PROVED + VERIFIED-EXACT. Every scalar five-unit/three-blocker cover with
  blocker valuations lambda_1<=lambda_2<lambda_3 has at least ten strict
  c_1-private and ten strict c_2-private guard-safe points on the primitive
  shell modulo 13^(lambda_3+1), and at least ten strict c_3-private points
  on the embedded shell of numerator valuation lambda_3-lambda_2. Thus its
  cyclic private-owner word contains all three labels, has at least thirty
  letters, and pays at least three owner switches. On the other hand, after
  quotienting the linear sheet slopes, the natural relative integer carry
  of two depth-two owners is a globally defined periodic potential, so its
  transition cocycle is exact and has zero loop holonomy. An exact profile
  (1,2,5) three-terminal diagram has a nonzero adjacent relative-carry jump
  of two, a common 68-obligation stalk, zero holonomy, and one persistent
  owner. Hence global owner switching is forced, but raw carry alone cannot
  locate or price it. An ancestry/gluing or support-gap sidecar is still
  required; no valuation profile is excluded and LRC(14) remains open.
source: codex-2026-07-25-depth-one-multi-terminal-cocycle
depends_on:
  - THM-2135-root-profile-invoices-and-first-deep-scalar-tail-closure
  - THM-2138-all-depth-unit-annulus-extremal-law
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
related:
  - THM-2211-carry-regime-root-transducer-and-infinite-autonomous-index
  - THM-2234-first-depth-one-private-two-owner-mass-and-two-step-expansion
  - THM-2246-depth-one-private-joint-two-step-fibre-cap
  - THM-2255-valuation-separated-pair-cap-and-exclusive-owner-mass
  - THM-2261-expiration-image-surjectivity-and-one-core-carrier-no-go
  - THM-2263-thirteen-adic-gap-pair-spectrum-and-profile-sharp-owner-floor
  - THM-2267-static-owner-coverage-is-flag-and-transition-holonomy-is-a-cut-kernel
script: 04-computation/lrc14_depth_one_multi_terminal_cocycle_stopping.py
output: 05-knowledge/results/lrc14_depth_one_multi_terminal_cocycle_stopping.out
script_sha256: 3806fee0aae6202cb738b8b44b150448be6c2468eb4c08d4f4c4d9683ff2a844
output_sha256: 6cc5671869ed92492731561d3d86f93b80afd58dc0d2bdc6b2a223b04c2b6f6c
hash_basis: working-tree bytes (LF)
---

# THM-2268 -- the private-owner trident and the exact raw-carry no-go

Use the scalar-cover notation

```text
D_a={x in R/Z:||ax||<1/14},
C_H={x in R/Z:||Hx||>1/7}.                           (1)
```

Suppose `H,q_1,...,q_5` are positive thirteen-units, `H` is odd, the
`q_i` are pairwise distinct, the three actual blockers `c_1,c_2,c_3` are
distinct positive multiples of thirteen, and

```text
C_H subset union_(i=1)^5 D_(q_i)
             union union_(j=1)^3 D_(c_j)             (2)
```

almost everywhere. Label the blockers as in THM-2198, so

```text
lambda_j=nu_13(c_j),
1<=lambda_1<=lambda_2<lambda_3.                      (3)
```

Then there are three pairwise disjoint finite sets `A_1,A_2,A_3` of strict
guard-safe torsion points such that

```text
|A_j|>=10,                                           (4)

x in A_j
  => x in C_H intersection D_(c_j)
     and x notin union_i D_(q_i)
     and x notin union_(k!=j)D_(c_k).                (5)
```

The sets `A_1,A_2` lie on the primitive numerator shell modulo
`13^(lambda_3+1)`. The set `A_3` lies on the same master grid, on the shell
of numerator valuation `lambda_3-lambda_2`.

Consequently, cyclically ordering `A_1 union A_2 union A_3` and writing the
unique owner label gives a cyclic word of length at least thirty, containing
all three labels. Its switch number is at least three.

This is the global positive statement. The second statement identifies the
corresponding local no-go. On a `169`-root terminal fibre, the raw relative
integer-carry transition of two depth-at-least-two blockers is always an
exact cocycle after the individual sheet slopes are removed. It can jump
across an adjacent chamber wall, but every closed overlap loop has zero
holonomy. The exact hostile diagram in Section 4 realizes a jump of two with
zero owner switches and a common full-simplex stalk.

## 1. The unit-annulus capacity debit is exactly ten

For `m>=3`, put `N=13^m`, `epsilon=(-1)^m`, and use THM-2138's primitive
guard annulus and one-mask maximum:

```text
U_N={z mod N:13 does not divide z and 7||z||_N>N},

|U_N|=(60N-130epsilon)/91,
M_m=(10N+130epsilon)/91.                             (6)
```

Every unit mask has size at most `M_m`. Every positive-valuation mask at
even `m` has size at most `M_m-20`; at odd `m` it has size at most `M_m`.
The parity arithmetic is

```text
m odd:   |U_N|-6M_m=10,

m even:  6M_m-|U_N|=10,
         |U_N|-[5M_m+(M_m-20)]=10.                  (7)
```

Thus five unit masks together with one positive-valuation mask miss at least
ten points of `U_N` at every `m>=3`.

The base `m=2` must not be inferred from the uniform formula without its
separate finite input. THM-2135 proves exactly that `|U_169|=110` and every
unit mask has size at most `20`; hence five unit masks miss at least ten
points. Its stronger exact union census in fact misses twenty-two.

No endpoint can hide this debit. Since `N` is a power of thirteen, neither

```text
7||z||_N=N
nor
14||rz||_N=N                                         (8)
```

can occur. Every guard and danger incidence on these grids is strict. If a
torsion point were missed by all eight danger combs, the failure would
persist on an open interval and contradict the almost-everywhere cover (2).

## 2. Two primitive-shell private sets

Set

```text
N=13^(lambda_3+1).                                   (9)
```

Reindex a physical numerator `r` as

```text
r=H^(-1)z mod N,              equivalently z=Hr mod N.        (9a)
```

This makes the guard coefficient one, keeps the five transformed
coefficients `q_i H^(-1)` as units, and preserves all three blocker
valuations.

On `U_N`, the blocker of valuation `lambda_3` is safe: its value is a
nonzero thirteenth root and therefore has norm at least `1/13>1/14`.
Omit `c_1`. The remaining possible danger masks are the five unit masks and
the positive-valuation mask from `c_2`. Equations (6)--(7) leave at least
ten points uncovered by those six masks. By (2) and the strictness audit
(8), `c_1` must own every such point. Choose ten of them for `A_1`.

Interchanging `c_1,c_2` gives ten points private to `c_2`; choose them for
`A_2`. Both sets are primitive modulo `N`. They are disjoint because a point
private to one blocker is, by definition, safe for the other.

Notice that this conclusion uses the global cover. THM-2246's sharp local
phase need not extend to such a cover and therefore need not contain either
private set.

## 3. The embedded deep-owner shell

Put

```text
N'=13^(lambda_2+1),
d=lambda_3-lambda_2>=1,                              (10)
```

and use the compatible reindexing (9a) modulo `N'`. On the primitive annulus `U_(N')`,
the blocker `c_2` is safe. If `lambda_1<lambda_2`, the other masks after
omitting `c_3` are five unit masks and one positive-valuation mask from
`c_1`; equations (6)--(7) leave at least ten points.

If `lambda_1=lambda_2`, both `c_1,c_2` are safe on `U_(N')`, leaving only
five unit masks. For `m=lambda_2+1>=3`, the formulas in (6) give

```text
|U_(N')|-5M_m
 =(10N'+780)/91    if m is odd,
 =(10N'-780)/91    if m is even,                    (11)
```

which is at least ten. When `lambda_2=1`, so `m=2`, use THM-2135's separate
`110-5*20=10` base certificate.

The omitted blocker `c_3` is divisible by `N'`, so it is dangerous at every
point of this grid. The global cover and (8) therefore turn ten of the
missed points into a set `A'_3` private to `c_3`.

Embed this shell in the master grid by

```text
iota: Z/N'Z -> Z/NZ,
iota(y)=13^d y.                                      (12)
```

The circle point is unchanged:

```text
iota(y)/N=y/N'.                                      (13)
```

For primitive `y`, the numerator `iota(y)` has exact valuation `d`. Define
`A_3=iota(A'_3)`. It has at least ten points and is disjoint from the
primitive sets `A_1,A_2`. This proves (4)--(5).

Order the at least thirty points in the union cyclically. After contracting
maximal constant-owner runs, the resulting cyclic word still contains all
three symbols, so it has at least three runs and therefore at least three
switches. This forced switch count is combinatorial; the argument does not
yet bound the arc length, mass, or carry cost between successive witnesses.

## 4. The raw relative carry is always exact

Let two blockers at a `169`-terminal be

```text
c=169a,                    d=169b.                   (14)
```

For a terminal phase `z` and sheet `0<=n<169`, their integer carry
potentials are

```text
p_c(n,z)=floor(c(z+n)/169)=an+floor(az),
p_d(n,z)=floor(d(z+n)/169)=bn+floor(bz).             (15)
```

After quotienting the forced sheet slopes `an,bn`, the terminal gauges are

```text
g_a(z)=floor(az),             g_b(z)=floor(bz).      (16)
```

The natural relative gauge which also quotients a change of circle lift is

```text
F_(a,b)(z)=a floor(bz)-b floor(az).                  (17)
```

Indeed,

```text
F_(a,b)(z+1)=F_(a,b)(z).                             (18)
```

Thus `F_(a,b)` is a globally defined integer potential on the terminal
circle. On any chamber-overlap edge `C->C'`, its transition

```text
omega(C,C')=F_(a,b)(C')-F_(a,b)(C)                  (19)
```

is an exact one-cocycle. Every closed overlap loop has

```text
sum omega=0.                                         (20)
```

This does not contradict THM-2211's infinite autonomous carry index. The
full carry stream is needed to propagate exact labelled masks. The present
statement says that one particularly natural scalarization of two terminal
carries has no cohomological obstruction until owner eligibility or another
sidecar restricts the available gauges.

## 5. Exact three-terminal hostile diagram

Take the profile `(1,2,5)` fixed speeds

```text
H=1,
c_1=26,                c_2=338,                c_3=742586,

q_i=1+169*700000i,                    1<=i<=5.       (21)
```

For a terminal phase `z`, put `x_n=(z+n)/169` and retain the obligations

```text
U(z)={n:
  ||x_n||>1/7,
  ||26x_n||>=1/14,
  ||q_i x_n||>=1/14 for every i}.                    (22)
```

The original THM-2246 hostile terminal is

```text
z_0=325007/700000.                                   (23)
```

It has `112` obligations, all owned by `c_2` and none by `c_3`.

Now take the simultaneous `c_3` carry wall

```text
z_*=157/338,
(c_3/169)z_*=4394z_*=2041.                           (24)
```

The nearest retained membership wall is exactly

```text
155/233248167061,                                    (25)
```

coming from `q_5` on sheet `72`; the nearest other raw carry wall is

```text
157/199927000338.                                    (26)
```

Choose

```text
epsilon=155/699744501183,
z_-=z_*-epsilon,
z_+=z_*+epsilon.                                     (27)
```

Hence `z_-` and `z_+` are adjacent in the fully retained local partition,
separated only by the simultaneous `c_3` carry wall. Exact enumeration gives

```text
U(z_-)=U(z_*)=U(z_+),
|U(z_-)|=68,
U(z_- ) subset U(z_0).                               (28)
```

Both `c_2,c_3` own all sixty-eight common obligations at the three phases in
(28), while `c_2` is also the unique owner of the `112` obligations at
`z_0`. Thus choosing `c_2` gives one persistent section over the entire
three-terminal common stalk.

The slope-quotiented gauges `(g_2,g_4394)` and relative potential are

```text
phase       gauge           F_(2,4394)

z_0         (0,2040)             4080
z_-         (0,2040)             4080
z_*         (0,2041)             4082
z_+         (0,2041)             4082.               (29)
```

Therefore the adjacent transition `z_- -> z_+` is nonzero:

```text
Delta gauge=(0,1),             Delta F=2.            (30)
```

But the overlap triangle

```text
z_0 -> z_- -> z_+ -> z_0
```

has edge values

```text
(0,2,-2)
```

and holonomy zero. Since the single legal `c_2` stalk services the whole
common universe, its coverability complex is the full `67`-simplex:

```text
2^68=295147905179352825856 faces,
no minimal nonfaces,
zero owner switches.                                 (31)
```

So even a nonzero adjacent relative-carry jump is not itself a switch tax.
The global trident in Sections 2--3 proves that a genuine cover cannot use
one persistent owner everywhere, but it does not force the three private
witness sets into one common ancestry fibre or prevent support gaps between
them. That ancestry/gluing problem is the next exact carrier.

## 6. Scope and reproduction

This theorem proves a global owner-switch **existence** statement and a local
raw-carry **no-go**. It does not assign positive Haar mass to the private
torsion witnesses, lower-bound the length of a switch arc, or show that a
switch intersects THM-2234's private measurable remainder. It excludes no
one of the `165` first-depth-one valuation profiles and does not prove
LRC(14).

The distinction from the measurable owner ledger is useful. THM-2255 and
THM-2263 give positive Haar mass to at least one exclusive labelled owner
and, in every strict profile, push one such stratum past the one-comb
threshold at expiration. The present theorem guarantees strict finite stalks
for **every** owner and forces a nonconstant global owner tour. Neither
result identifies which successor sheets receive the expanding measurable
stratum. THM-2267's transition cut kernel is therefore the exact next
consumer: it needs the common sheet/carry correspondence which (17)--(20)
and the hostile diagram show cannot be reconstructed from raw carry alone.

Reproduce the finite diagram with

```bash
python3 04-computation/lrc14_depth_one_multi_terminal_cocycle_stopping.py
python3 -O 04-computation/lrc14_depth_one_multi_terminal_cocycle_stopping.py
```

Both modes reproduce the stored transcript byte for byte. The companion
checks all `169` sheets, every membership boundary, the nearest other carry
wall, the sheet-potential quotient, owner labels, common overlap, relative
transition, loop holonomy, and full-simplex consequence. QED.
