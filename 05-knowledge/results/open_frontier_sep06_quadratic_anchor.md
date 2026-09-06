# The quadratic anchor excludes the escaping-phase obstruction

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.** This is a
structural restriction on the q=5,h=4 model with both C and D interlacing B.
It does not prove full-response negativity at every original first-row root,
and does not prove actual Laurent noncancellation.

## 1. Inheritance and the new constraint

The [linearly anchored hostile](open_frontier_sep06_laurent.md) is proved for
all real geometric ratios R>=384: e1=13 leaves a degeneration toward one
positive B root, and the original first phase escapes to infinity with a
positive full response. Its algebraic shape continuation also gives exact
model double cancellation. This is the canonical hostile for the present
repair, rather than a counterexample to the actual coefficient law.

Now impose the genuine second coefficient as well:

```text
e1(a)=13,  e2(a)=55,  equivalently sum a_i=13, sum a_i^2=59.
```

The incoming [coupled-window cone](creative_20260906_laurent_bridge.md) remains
a sufficient signed-response certificate only after the same-carrier equality
to its weights has been supplied. No such general equality is claimed here.
The five live objects are the two exact coefficient anchors; the compact
positive-root shape space; both contiguous interlacers; the original first
zero; and the sign of the complete carried response. The underused sidecar
is the number of roots surviving at the boundary of this shape space.

The corrected near miss is to infer negativity merely because a finite
anchored probe finds no hostile. The result below instead rules out the
specific escaping-phase mechanism analytically. Its scope includes weak
interlacing limits of the strict certificate class; it is not a classification
of every possible all-weight pencil without those interlacers.

## 2. Compactness forces a positive fourth coefficient

Let `0<=a_0<=...<=a_4`, with the two displayed anchors, and set

```text
B(v)=product(v-a_i)=v^5-13v^4+55v^3-e3 v^2+e4 v-e5,
C(v)=v^4-12v^3+45v^2-(2/3)e3 v+(3/7)e4,
D(v)=v^4-11v^3+36v^2-(5/12)e3 v+(1/7)e4.
```

Assume both monic degree-four polynomials C,D weakly interlace B. Let K be
this closed shape class, with repeated and zero B roots allowed. K contains
the closure of the strict positive-root class used by the earlier pencil
certificate. It is compact: the roots lie in [0,13], the coefficient map is
continuous, and the weak root-order inequalities are closed.

**Claim.** There is a constant kappa>0 such that `e4>=kappa` throughout K.
No numerical value for kappa is supplied here.

Suppose e4=0. Nonnegativity implies at most three positive roots. There
cannot be only two: Cauchy would give `13^2<=2*59`, which is false.
Hence

```text
B(v)=v^2 F(v),
F(v)=v^3-13v^2+55v-e3,
0<r1<=r2<=r3
```

are its remaining roots. The exact identity

```text
C(v) - v^2(v-5)^2/3 = (2/3)v F(v)                     (1)
```

gives `C(r2)=r2^2(r2-5)^2/3>=0`. Weak interlacing at the fourth ordered B
root gives `C(r2)<=0`. Therefore r2=5. The two anchors now force
`r1+r3=8` and `r1*r3=15`, so `(r1,r2,r3)=(3,5,5)` and e3=75.
Weak interlacing of D between the repeated B roots at 5 would require D(5)=0,
but direct evaluation gives `D(5)=-25/4`. This is a contradiction.

Thus e4 is strictly positive at every point of compact K, proving the
uniform positive minimum. In particular two B roots cannot simultaneously
approach zero in the strict class. The would-be boundary is a useful hostile:
`B=v^2(v-3)(v-5)^2` and `C=v(v-2)(v-5)^2` satisfy the first interlacing,
while D rejects it. The second contiguous row performs an essential step.

## 3. Every escaping original phase has a negative full response

Reverse B,C,D with alternating signs and retain the original shifts
`beta=t^-1 G_B`, `C_raw=t^-1 G_C`, `D_raw=t^-1 G_D`. Their first two
coefficients are automatically the genuine contiguous anchors:

```text
G_B=(1,13,55,e3,e4,e5),
G_C=(1,12,45,(2/3)e3,(3/7)e4),
G_D=(1,11,36,(5/12)e3,(1/7)e4).
```

Keep `O=sum binom(14,2j+1)t^j`, `E=sum binom(14,2j)t^j`, and define

```text
P=O star beta,
Q=(O^2+t^-1 E^2) star(beta^2+2t C_raw D_raw).
```

At the same original root t=-s<0, the complete first row is

```text
P(-s)=182-20020s+2002e3 s^2-3432e4 s^3+2002e5 s^4.
```

Consequently the following relation is exact at any such root:

```text
e5 s=(12/7)e4-e3/s+10/s^2-1/(11s^3).                 (2)
```

All Q coefficients are continuous functions on compact K. Its top two raw
coefficients, with the negative lower carry still retained, are

```text
Q_8=binom(28,18)e5^2,
Q_7=binom(28,16)(2e4e5+(6/49)e4^2).                  (3)
```

Now consider any sequence of shapes in K and their original positive phases
s_n with `P_n(-s_n)=0` and `s_n->infinity`. Pass to a convergent subsequence
of shapes. Write its limiting fourth coefficient as e4_*>0. Equation (2)
gives `e5_n->0` and `e5_n s_n->(12/7)e4_*`. Dividing the full response by
s_n^7 in (3), every lower term tends to zero, as does the top term
`e5_n^2 s_n=e5_n(e5_n s_n)`. Therefore

```text
Q_n(-s_n)/s_n^7 -> -(6/49)binom(28,16)e4_*^2 < 0.     (4)
```

The exponent-seven sign comes from the complete crossing coefficient;
discarding that coefficient or the Laurent shift would change the conclusion.

**Uniform tail consequence.** There is a finite S such that every shape in
K and every positive original root s>S satisfy `Q(-s)<0`. Otherwise a
sequence with s_n>=n and Q_n(-s_n)>=0 would contradict (4). This is a
qualitative uniform threshold, with no numerical S obtained. It does not
turn the remaining finite-phase problem into a verified finite census.

## 4. Exact nonvacuity and the remaining question

The simple rational shape

```text
a=(1,3,9,22,30)/5,
(e3,e4,e5)=(2127/25,27144/625,3564/625)
```

has both exact anchors and both strict interlacings. Its four original
positive phase roots and all four negative full responses are verified by
the accompanying exact control. The boundary tuple `(0,0,3,5,5)` preserves
the two coefficient anchors and the C interlacing, but is rejected by D.
These are controls for the proof and its necessary hypothesis, not a claim
that the full anchored class has been exhausted.

A bounded rational-sphere reconnaissance found no positive full response;
that observation is only a discovery diagnostic. The exact gain is the
positive e4 boundary restriction and the resulting uniform negative tail.
The next target is an explicit lower bound for e4, or a same-carrier sign
certificate on the remaining compact finite-phase class. Full Q-negativity
under two anchors remains **OPEN**.

## 5. Exact verification

The [standalone source](../../04-computation/open_frontier_sep06_quadratic_anchor.py)
and [frozen output](open_frontier_sep06_quadratic_anchor.out) pass **57 exact
gates**, with byte-identical normal and optimized output. The declared finite
universe is the displayed rational shape and all four of its original roots,
the rejected boundary `(0,0,3,5,5)`, and two exact evaluations of identity (1)
whose coefficient degree in e3 is proved to be at most one. The source uses
Fraction arithmetic, no repository producer, and no floating-point roots.

The four phase brackets are disjoint, each has an endpoint sign change,
and they exhaust degree four. The full response is strictly negative on
each entire bracket by rational interval evaluation. Independent ordinary
carriers reproduce the first-row and full carried response at rational phases.
The exact top two coefficients and the exponent -1 coefficient are retained.
These controls support the proof; compactness and the sequence argument
supply its uniform quantifiers.

```bash
python3 -B 04-computation/open_frontier_sep06_quadratic_anchor.py
python3 -B -O 04-computation/open_frontier_sep06_quadratic_anchor.py
```

SHA256:

```text
source   ebe6101b08f748d710b7bb3894ce064e476ce80eaf2758995c2b07bf1e5ae9d7
output   6e04d52afd8b8961cc9d95c498c0d36abf95b25c0d81ff7f50e9c1ac4c300fe4
semantic a2fa833d88b725985e8d071936f274498498b3707b295bdccbe6731373419eaa
```

Root independently audited the compact weak-interlacing closure, the
two-zero contradiction including repeated cases, the exact original-root
relation, both carried top coefficients and the uniform-tail consequence:
**PASS**. Root then independently read the complete source and reproduced
normal, optimized and frozen output agreement for all 57 gates: **PASS**.
The review includes the generic boundary identity, full top and lower carry,
four distinct root brackets exhausting degree four, and the separate ordinary
coefficient extraction. All three files are frozen. No earlier source, shared
navigation, theorem ID or Git state is changed by this bundle.
