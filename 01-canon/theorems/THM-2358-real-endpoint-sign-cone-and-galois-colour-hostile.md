---
id: THM-2358
title: "Real endpoint sign-cone boundary and Galois-colour hostile"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. For a
  nonzero real array, confinement to a phase cone of width less than pi
  is exactly sign-definiteness. Hence the THM-2355 cone route for
  THM-2344's real endpoint arrays is equivalent to proving both arrays
  sign-definite; if either support is nonsingleton this excludes the
  shifted-delta boundary. A rational N=91 nested-word example has
  nonzero fixed-colour correlations of both signs, proving that
  colourwise Galois straightening gives no common physical cone. More
  sharply, a literal centered (1,13,26), R=169 target carrier has a
  nonzero marked word/deep/bare triangle while its real inverse endpoint
  array has U_0>0 and U_1<0. Positivity, evenness, nonconstant word
  content, two target axes, and marked-triangle nonvanishing therefore do
  not imply an acute cone. The carrier is local, not a canonical
  nine-coordinate scalar cover; no target landing, row exclusion, or
  LRC(14) closure is proved.
source: codex-2026-07-25-real-endpoint-sign-cone
depends_on:
  - THM-2323-primitive-fixed-colour-cross-correlation-and-same-gauge-word-alignment
  - THM-2344-correlation-inverse-rigidity-and-aligned-tooth-twist-hostile
  - THM-2355-component-deletion-gram-and-twist-energy-phase-transport
related:
  - THM-2356-finite-field-chirp-gram-tomography-and-bockstein-pairing
script: 04-computation/real_endpoint_sign_cone_hostile_thm2358.py
output: 05-knowledge/results/real_endpoint_sign_cone_hostile_thm2358.out
script_sha256: 8aea0df87e4584ea79df9531f5faf41232b2823b0c186dc4432a9a1439f1b101
output_sha256: 5820d11eb06f38c12caf393f6d969990e36093d5ae7bb1081a00a22914f70315
hash_basis: working-tree bytes (LF)
---

# THM-2358 -- the real endpoint sign-cone boundary

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2355 gives a sharp sufficient condition for an endpoint
cross-correlation to have honest support: the two arrays must lie in phase
cones whose widths sum to less than `pi`. THM-2344 proves that the actual
endpoint residue arrays are real. These statements collapse the geometric
question to an exact sign question.

That collapse produces both a useful positive target and a sharp no-go.
Sign-definite endpoint arrays would finish the inverse-character
separation as soon as one support is nonsingleton. But positivity and
evenness of the physical interval factors do not imply sign-definiteness
after target inversion. This remains false even with a nonconstant
transported word and a literal nonzero marked triangle.

## 1. For real arrays, the cone condition is sign-definiteness

For a finite complex array `Z`, define its phase width to be the infimum
of the lengths of circular intervals containing the arguments of all
nonzero values of `Z`.

> **Real-cone lemma.** If `Z` is a nonzero real array, its phase width is
> less than `pi` if and only if `Z` is sign-definite:
>
> ```text
> Z>=0,                    or                    Z<=0.              (1)
> ```
>
> In that case its minimal width is zero. If both signs occur, its
> minimal width is exactly `pi`.

Indeed, the only nonzero phases of a real array are `0` and `pi`. One of
them occurs exactly in the sign-definite cases; both occur exactly when
the array has both signs, and they are antipodal.

THM-2344's endpoint arrays `U,V` are real, and `K(0)!=0` makes both
nonzero. Therefore THM-2355's hypothesis

```text
phase_width(U)+phase_width(V)<pi                    (2)
```

is equivalent in this application to

```text
U and V are both sign-definite.                     (3)
```

There is a sharp positive consequence. If nonzero `U,V` are
sign-definite, every nonzero summand in their cross-correlation has the
same sign. Hence

```text
support R_(U,V)=support(U)-support(V).               (4)
```

If the bad THM-2344 boundary held, this support would be the singleton
`{-p}`. For nonempty finite sets `A,B`,

```text
A-B is a singleton
 iff A and B are both singletons.                   (5)
```

To see the reverse direction in (5), fix `b in B`; then every `a in A`
must equal the unique difference plus `b`, so `A` is a singleton, and
symmetrically so is `B`.

Thus

```text
U,V sign-definite
 + max(|support U|,|support V|)>1

 => some nonzero full target survives.              (6)
```

This criterion is sharp: THM-2344's aligned point masses are
sign-definite singleton arrays and attain the inverse-character boundary.

## 2. Colourwise Galois acuteness is not one physical cone

The acute sector in THM-2323 is produced after a Galois automorphism that
depends on the chosen primitive colour. The following rational step
functions show why those sectors cannot be identified across colours.

Put `N=91`, choose rational

```text
0<epsilon<1/182,
```

and on the circle define

```text
A_0=(-epsilon,epsilon),

A_+=A_0+6/91,                 A_-=A_0-6/91,

f=1_(A_0),                    g=1_(A_0 union A_+ union A_-).       (7)
```

The three arcs are disjoint. Since

```text
6/91+epsilon<1/14,                                 (8)
```

after reducing `epsilon` if necessary, `f,g` are even rational-valued
step functions satisfying

```text
0<=f<=g<=1,

support(f),support(g) subset D_1.                  (9)
```

In THM-2323's fixed-colour cross-correlation, the only shifts are
`0,+6,-6`, all with the same positive weight

```text
c=91*(2epsilon)=182epsilon.                        (10)
```

For `zeta=exp(2*pi*i/91)`,

```text
C_k=c(1+zeta^(6k)+zeta^(-6k)).                     (11)
```

Both of the following colours are primitive:

```text
k=76:       6k=1 mod 91,

k=53:       6k=45 mod 91.                          (12)
```

Consequently

```text
C_76=c(1+2cos(2*pi/91))>0,

C_53=c(1-2cos(pi/91))<0.                           (13)
```

The second inequality uses `pi/91<pi/3`. None of the primitive values in
(11) vanishes. If `1+z+z^(-1)=0`, then `z` has order three; but for a
primitive `k`, `z=zeta^(6k)` has order `91`.

The same nested positive word therefore produces both physical signs and
has phase width exactly `pi`. THM-2323 remains valid: applying its
colour-dependent Galois automorphism to a hypothetical zero straightens
that one colour into an acute sector. What fails is an intertwiner turning
all those distinct algebraic automorphisms into one common archimedean
phase rotation.

## 3. A literal target-level centered-factor hostile

The preceding example acts at the fixed-colour cross-correlation level.
There is a stronger hostile directly in THM-2344's target-coordinate
representation.

Let

```text
d(x)=1_(||x||<1/14),                g=1-d,

w=(1,13,26),                        R=169,

X=26,                 c_3=26,       m=1,       Y=52.              (14)
```

For target shifts `s,t in F_13`, put

```text
L_(s,t)
 =integral_T
   d(13x+s/13)
   g(26x+t/13)
   g(169x)
   exp(-2*pi*i*26x) dx.                            (15)
```

At the untwisted point every factor is a literal centered danger or
complement. The last factor is a genuinely nonconstant transported
complement word.

Define the normalized inverse target array

```text
U(a,b)
 =1/169 sum_(s,t in F_13)
    zeta^(-as-bt)L_(s,t),

zeta=exp(2*pi*i/13).                                (16)
```

The atomic Fourier indices obey

```text
13u+26v+169 beta=26,

u+2v+13 beta=2.                                    (17)
```

Thus `U` is supported on the target line

```text
q_r=(2-2r,r),                    r in F_13,          (18)
```

and

```text
U_r:=U(q_r)
 =1/13 sum_(t in F_13)zeta^(-rt)L_(0,t).            (19)
```

All `U_r` are real: expansion into the centered even base Fourier
coefficients makes every atomic weight real, equivalently reflection
conjugates `(s,t)` with `(-s,-t)` in (16).

### The coefficient `U_1` is strictly negative

Split the transported complement as

```text
g(169x)=1-d(169x)
```

and accordingly write

```text
U=U^B-U^T.                                         (20)
```

On the bare term parameterize each tooth of `D_13` by

```text
x=(k+y)/13,       0<=k<=12,       |y|<1/14.         (21)
```

Then

```text
B_t
 =integral_(-1/14)^(1/14)
   g(2y+t/13) exp(-4*pi*i*y)dy.                    (22)
```

Let

```text
T(y)={t in F_13:||2y+t/13||<1/14}.                 (23)
```

Since the complete character sum is zero, (19) and `g=1-d` give

```text
U^B_1
 =-1/13 integral
    sum_(t in T(y))
    exp(-2*pi*i(2y+t/13))dy.                        (24)
```

Every summand in (24) has real part greater than `cos(pi/7)`, and
reflection makes the integral real. Exact rational interval intersection
gives the complete table, in denominator `364`,

```text
t= 0: (-13, 13),       length 26;
t= 1: (-26, -1),       length 25;
t=-1: (  1, 26),       length 25;
t= 2: (-26,-15),       length 11;
t=-2: ( 15, 26),       length 11;                  (25)

all other t: empty.
```

Therefore the total incidence is

```text
(26+25+25+11+11)/364=49/182,
```

and

```text
U^B_1<-49 cos(pi/7)/2366.                          (26)
```

For the tooth correction,

```text
|U^T_1|
 <=measure(D_13 intersection D_169).               (27)
```

Parameterize a `D_169` tooth by `x=(K+y)/169`. The `D_13` condition
holds on the tooth exactly when `K=0 mod 13`; the other twelve residue
classes miss even at the open boundary. Hence

```text
measure(D_13 intersection D_169)
 =13*(1/(7*169))=1/91=26/2366.                     (28)
```

Finally,

```text
cos(pi/7)>cos(pi/4)>2/3,

49*(2/3)>26.                                       (29)
```

Equations (20), (26)--(29), and the reality of `U_1` give

```text
U_1<0.                                             (30)
```

### The coefficient `U_0` and the marked triangle are positive

Put

```text
A(y)=1/13 sum_t g(2y+t/13).                         (31)
```

The table (25) shows that at most two danger roots are active, so

```text
A(y)>=11/13.
```

It is even. Formula (19) gives

```text
U_0
 =integral_(-1/14)^(1/14)
   g(13y)A(y)exp(-4*pi*i*y)dy.                     (32)
```

The imaginary part cancels by evenness, `cos(4*pi*y)>0` throughout the
interval, and `g(13y)` is nonzero on positive measure. Thus

```text
U_0>0.                                             (33)
```

So the real endpoint array already contains both signs.

The untwisted left coefficient is also positive:

```text
L_(0,0)
 =integral_(-1/14)^(1/14)
   g(2y)g(13y)exp(-4*pi*i*y)dy>0.                  (34)
```

The integrand has positive real part, and the interval

```text
(1/28,1/26)
```

is a literal positive-support witness. The bare right endpoint is

```text
F_(0,0)
 :=(d(13x)g(26x))_hat(52)

 =2 integral_(1/28)^(1/14) cos(8*pi*y)dy

 =[sin(4*pi/7)-sin(2*pi/7)]/(4*pi)>0,              (35)
```

because `sin(4*pi/7)=sin(3*pi/7)>sin(2*pi/7)`.
The deepest leg is

```text
(d(26x))_hat(26)=sin(pi/7)/pi>0.                   (36)
```

Equations (34)--(36) give a literal nonzero word/deep/bare Fourier
triangle at `X=26,Y=52,m=1`, while (30) and (33) show that its left
endpoint target array has phase width `pi`.

This carrier preserves:

```text
centered nonnegative danger/complement factors;
real/even base symmetry;
a nonconstant R=169 transported word;
two independent target axes;
and one nonzero marked triangle.                   (37)
```

It is a local three-coordinate factorized carrier, not a canonical
nine-coordinate scalar cover and not one of the `165` rows. It refutes
only an implication from the qualitative data in (37) to an acute cone.

## 4. Frontier effect

The cone route is now exact rather than qualitative:

```text
prove both grouped endpoint arrays sign-definite
  and prove one support nonsingleton;

or construct a genuinely coefficient-sensitive nonlinear probe,
  such as THM-2356's candidate planar chirp service.                 (38)
```

Physical factor positivity, evenness, a nonconstant word, an extra target
axis, the fixed-colour Galois argument, and marked-triangle nonvanishing
do not prove the first line. THM-2356 is only a candidate under audit and
is not used as a dependency here; its nonlinear-mask route is
mathematically complementary to this amplitude-level half-plane no-go.

No canonical sign theorem, planar target--jet coupling, target landing,
scalar-row exclusion, or LRC(14) closure is proved. The ledger remains
`165`.

## 5. Exact companion

The dependency-free companion uses `Fraction` arithmetic to:

- derive every interval in (25), the exact `49/182` incidence, and the
  maximum of two simultaneous danger roots;
- enumerate all `169` deep teeth, finding exactly `13` in the
  `D_13` intersection and the exact mass `1/91`;
- verify the rational strict sign gap and the thirteen-point support line;
- certify the positive support and angle orderings in (34)--(36); and
- instantiate (7) with `epsilon=1/1000`, checking both primitive sign
  colours and the order-three nonvanishing obstruction.

Reproduce with

```bash
python3 04-computation/real_endpoint_sign_cone_hostile_thm2358.py
python3 -O 04-computation/real_endpoint_sign_cone_hostile_thm2358.py
```

Both transcripts must match

```text
05-knowledge/results/real_endpoint_sign_cone_hostile_thm2358.out
```

byte-for-byte after LF normalization. Every executable check raises
explicitly under optimized Python.

Independent audit is pending. QED.
