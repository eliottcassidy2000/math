---
id: THM-2800
title: "Two-pole two-double-zero bitangent eliminant and complete Nielsen corridor"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  The complete balanced-response corridor e=2,h=2,p=(N-1,1) has
  exactly N-3 normalized maps in every degree N>=4.  Their two double zeros
  are the common-tangency points of a line to x^(N-1)(x-1), and an explicit
  degree-(N-3) recurrence polynomial Q_N is their exact Maxwell/bitangent
  eliminant.  The N=4 map has A4=V4 semidirect C3 monodromy.  At N=5 the
  two first Nielsen classes are the conjugate roots of 5a^2+2a+1 over Q(i)
  and have C5 semidirect C4 monodromy.  This is a response-layer theorem,
  not Keller-chart entry, JC(2), or DC(2).
source: root/jc-e2-h2-recurrence-2026-07-28
depends_on:
  - THM-2796-balanced-response-stieltjes-pade-normal-form-and-one-double-zero-classification
related:
  - THM-2799-one-pole-two-double-zero-chebyshev-response-classification
  - THM-2784-nonsplit-response-square-potential-divisor-and-infinity-classification
  - THM-2768-modular-c2-c3-a4-s4-bass-serre-quotient
script: 04-computation/jc2_two_pole_e2_bitangent_nielsen_thm2800.py
output: 05-knowledge/results/jc2_two_pole_e2_bitangent_nielsen_thm2800.out
script_sha256: 20b5b0e1c9383cfbb58dd53c317e3bbd55c23f7de89a07a7b6eef5f7f5fe4955
output_sha256: 055c62af3bf00c563a868719f95f6abf2529c16fdb9ea37ed9510f5cefd08389
hash_basis: LF-normalized bytes
---

# THM-2800 -- two-pole two-double-zero bitangent eliminant and complete Nielsen corridor

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

THM-2799 classifies the one-pole part of the first open `e=2` response
layer.  The present theorem classifies the adjacent two-pole corridor.  Its
missing Nielsen coordinate has a concrete algebraic meaning: it chooses a
bitangent of the binomial graph

```text
y=D(x)=x^(N-1)(x-1).                                (1)
```

The first passport multiplicity from THM-2796 is therefore not an
unstructured dessin accident.  It is the first self-intersection of the
dual/Legendre curve of `(1)`.

## 1. Exact bitangent eliminant

Work over an algebraically closed field of characteristic zero.  Fix

```text
e=2,                    h=2,                 p_poles=(N-1,1),
s_simple=N-4,           r=N-1,               N>=4.     (2)
```

Normalize the order-`N-1` pole to `0`, the simple pole to `1`, the
order-`N-1` point over the third branch value to infinity, and
`F(infinity)=1`.  Put

```text
D=x^(N-1)(x-1),                 T=x(x-1),
E=(x-gamma)(x-delta)=x^2-sx+p.                 (3)
```

The two roots of `E` are distinct and avoid `0,1`.  Since the third fibre
has partition `(N-1,1)`, a monic numerator `B=S E^2` satisfies

```text
F=B/D,                          B-D=A(x-z),       A!=0. (4)
```

At either double zero `u in {gamma,delta}`, equations `B(u)=B'(u)=0`
say

```text
A=-D'(u),
z=phi(u):=u-D(u)/D'(u)
 =u ((N-1)u-(N-2))/(Nu-(N-1)).                    (5)
```

Thus `gamma,delta` are two tangency points of the same line to `(1)`.
Equality of the two intercepts in `(5)` gives

```text
p=p_N(s):=((N-1)s-(N-2))/N.                        (6)
```

Define the complete-homogeneous recurrence

```text
h_0=1,                    h_1=s,
h_k=s h_(k-1)-p_N(s) h_(k-2).                       (7)
```

It is characterized by

```text
x^k = h_(k-1)x-p h_(k-2)             modulo E.      (8)
```

Equality of the two tangent slopes is

```text
[D'(gamma)-D'(delta)]/(gamma-delta)
 =H_N(s):=N h_(N-2)-(N-1)h_(N-3)=0.                (9)
```

The diagonal collision

```text
gamma=delta=(N-2)/N,
s=2(N-2)/N                                           (10)
```

is a root of `H_N`, because `D''((N-2)/N)=0`.  The recurrence shows
directly that the corresponding linear factor divides `H_N`.  Define

```text
Q_N(s)
 :=H_N(s)/(s-2(N-2)/N).                              (11)
```

Then

```text
deg Q_N=N-3.                                         (12)
```

Every normalized response in `(2)` gives a root of `Q_N`.

## 2. Algebraic converse

Before the global root count is used, let `s` be a root of `Q_N` satisfying
the evident squarefree/disjointness gates.  In the quotient by `Q_N`, set

```text
A=p [N h_(N-3)-(N-1)h_(N-4)],
z=p,
B=D+A(x-p),
S=B/E^2,
C=-(N-1)A.                                           (13)
```

Modulo `E`, equation `(8)` gives

```text
D'=-A,                    D=-A(x-p).                 (14)
```

Hence `E^2` divides `B`, and `S` is monic of degree `N-4`.  A direct
calculation using `(6)` gives

```text
D-(x-p)D'=-(N-1)x^(N-2)E.                           (15)
```

Consequently

```text
F=B/D,
G=C E/(2DT),
V=4 S D T^2/C^2                                     (16)
```

satisfy

```text
F'=2G,                    F=VG^2,
2VG'+V'G=2.                                         (17)
```

The provisional converse gate is exactly

```text
Disc(E) A Disc(S)
Res(S,E x(x-1)) Res(E,x(x-1)) !=0.                  (18)
```

The Nielsen count below proves that every root of `Q_N` passes `(18)`.

## 3. Exact Nielsen count and automatic admissibility

Put `L=N-1`.  Fix the third branch permutation

```text
sigma_1=(0 1 ... L-1)
```

and let it fix one additional sheet `*`.  The zero inertia is a product of
two disjoint transpositions.  Transitivity forces one transposition to
contain `*`.  Conjugation by the centralizer `<sigma_1>` normalizes that
transposition to

```text
(*,0).
```

Write the other one as `(b,c)`, with

```text
1<=b<c<=L-1.                                         (19)
```

For `tau=sigma_0 sigma_1`, a direct orbit trace gives the two cycle lengths

```text
c-b,                    L+1-(c-b).                  (20)
```

The pole partition is `(L,1)` exactly when

```text
c=b+1.                                               (21)
```

Therefore the normalized Nielsen classes are indexed by

```text
b=1,...,L-2,
```

and their exact number is

```text
N-3.                                                 (22)
```

The genus-zero permutation triples give normalized rational maps by the
Riemann existence correspondence.  The pole multiplicities and the high
third-fibre point make the source normalization unique.  A normalized map
determines `s`; conversely `(13)` makes the map unique once `s` is fixed.
Thus the `N-3` Nielsen classes inject into the roots of the degree-`N-3`
polynomial `Q_N`.  It follows that:

1. `Q_N` has exactly `N-3` simple roots;
2. every root passes all gates `(18)`;
3. `(13)--(17)` construct every map; and
4. the corridor `(2)` has exactly `N-3` normalized maps in every degree.

This is the complete all-degree classification.

## 4. Degree four: the exact `V4`/`C3` response

For `N=4`,

```text
s=1/2,                 p=-1/8,
E=x^2-x/2-1/8,
B=E^2,
C=-3/8,                                         (23)

F=E^2/[x^3(x-1)],
V=(256/9)x^5(x-1)^3.                            (24)
```

The branch partitions are

```text
(2,2),                 (3,1),                 (3,1).
```

The monodromy has order `12` and is `A4`.  Its double transpositions form
the normal Klein four group, and

```text
A4/V4 isomorphic to C3.                              (25)
```

This is an exact occurrence of the repo's binary/ternary motif on a
response map: an order-two inertia sits in the `V4` layer and the
order-three rail acts on it.  It is not the proposed `S4/V4 isomorphic S3`
quartic resolvent.  All branch permutations here are even, so the odd
inertia needed for `S4` is absent.  Nor does `(24)` assert entry into a
Keller chart.

## 5. Degree five: the first Nielsen pair over `Q(i)`

For `N=5`, put

```text
a=1-2s.                                              (26)
```

Then

```text
Q_5(s)=5s^2-6s+2,
5a^2+2a+1=0,                a=(-1 plus/minus 2i)/5. (27)
```

Define

```text
u=(a-1)/2,                    v=-(2a+1)/5,
E=x^2+ux+v,                   S=x-a,
C=4(a-2)/25.                                         (28)
```

The two maps are

```text
F=(x-a)E^2/[x^4(x-1)],
F-1=(2-a)(x-v)/[25x^4(x-1)],                        (29)

G=C E/[2x^5(x-1)^2],
V=4(x-a)x^6(x-1)^3/C^2.                             (30)
```

They are exchanged by the Galois involution of `Q(i)`.  Relative to the
fixed four-cycle `(0 1 2 3)`, exact representatives of the two zero
inertias are

```text
(1 2)(3 4),                 (1 4)(2 3).             (31)
```

Both monodromy groups have order `20` and are the Frobenius group

```text
C5 semidirect C4.                                    (32)
```

Thus the first lost Nielsen coordinate is literal: it records the placement
of a double-transposition involution against a `C4` rail.  The integer
passport forgets that placement; the quadratic bitangent coordinate
recovers it.

## 6. Exact controls

The companion uses two independent representations.

The algebraic side, through `4<=N<=10`, checks:

1. the recurrence, diagonal factor, and degree of `Q_N`;
2. nonzero discriminant and every converse gate `(18)`;
3. exact division `E^2 | B`;
4. the linear third-fibre defect and response identities `(15)--(17)`; and
5. the `N=4` and `N=5` specializations over exact number fields.

The permutation side checks all normalized pairs in `(19)` through `N=10`,
the cycle-length law `(20)`, the adjacent-chord criterion `(21)`, all eight
raw transitive `N=5` candidates, the two centralizer orbits, and monodromy
orders `12` and `20`.  The script contains no Python `assert` node.

Run:

```text
python 04-computation/jc2_two_pole_e2_bitangent_nielsen_thm2800.py
python -O 04-computation/jc2_two_pole_e2_bitangent_nielsen_thm2800.py
```

The bounded controls do not replace the all-degree recurrence and Nielsen
proofs.

## 7. Scope and next target

This is a complete abstract response-layer classification.  It does not
show that any member satisfies the inherited Faber flux equations or comes
from a polynomial Keller pair.  Hence it proves neither JC(2) nor DC(2).

Together THM-2799 and the present theorem close the `e=2` corridors with
one pole and with pole partition `(N-1,1)`.  The next exact chamber is
`e=2,h=2` with a genuinely variable pole partition `(d,N-d)`.  There the
bitangent curve is no longer the single binomial `(1)`; the cheapest test is
to retain the pole cross-ratio and ask whether its dual curve has an
equally finite Maxwell eliminant before imposing the Keller flux equations.
