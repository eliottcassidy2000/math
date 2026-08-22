---
id: THM-3497
title: "Berggren T4 calibration no-go, two word languages, and harmonic parity smoothing"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  The true mod-two linear
  sidecar blocks static and six-frame T4 calibration for Berggren generators
  A and C, while B has exactly 4 pinned, 8 with translation allowed, and 512
  six-frame calibrations per lawful translation.  Allowing a new translation
  for each composite word gives a minimal S3 x C2 language of density 1/3;
  retaining THM-3339's fixed drift gives a distinct minimal S4 x D4 language
  with alternating level densities 1/6 and 3/16, no shortlex natural density,
  and logarithmic density 17/96.  The root-to-child convention is fixed by a
  first-step hostile, and the two harmonic arguments are cylinder-uniform.
source: berggren-t4-three-language-independent-audit-2026-08-16
audit: >
  Independent standard-library reconstruction from THM-3339's integer
  matrices accepts the substantive claims of commits 62457a7db, 9459b7693,
  85f20bfde, and 4aae3772a.  It exhausts the finite gauges and both group
  images, checks minimality and recurrences, verifies all three Fibonacci rays
  through 48 periods, and proves the fixed harmonic cylinder gate from the
  96-state two-step kernel plus return lengths 2 and 3.  The sole repair is
  terminological: 17/96 and -1/96 are normalized polar coefficients, not
  analytic residues.  Normal and optimized replays match all stored outputs.
depends_on:
  - THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction
related:
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
  - THM-3359-modular-c-finite-supports-harmonic-density-and-periodic-scar
  - THM-3382-fibonacci-ray-dual-index-harmonic-bifurcation-and-ternary-heap-addresses
  - THM-3487-two-twenty-four-state-fibonacci-bundles-cycle-type-obstruction
script: 04-computation/berggren_t4_three_language_independent_audit_20260816.py
output: 05-knowledge/results/berggren_t4_three_language_independent_audit_20260816.out
script_sha256: 22e0cac6bc9294418df7ac667321dd8a091f0ec302d8be724c6593847478d891
output_sha256: 95fbcfa549921dc13cc65615a167dbb735625f97cbe604a4eda65c484779cfcf
semantic_sha256: 0c363166248508d0712542fd725df56b32769680ba9acbf8ee63b4a6115f88cf
full_branch_script: 04-computation/berggren_full_branch_t4_static_calibration_no_go_20260816.py
full_branch_output: 05-knowledge/results/berggren_full_branch_t4_static_calibration_no_go_20260816.out
full_branch_script_sha256: 557b18c515659ca6b7acdbf4c8d66dddeb07849f85e3514780d89a6a76c2c6a4
full_branch_output_sha256: 77b7a038cf4a6fdb7f1b609729f5e387af8ea8290ef2a976c31d2141e6e85dc3
variable_script: 04-computation/berggren_transplantable_word_language_harmonic_probe_20260816.py
variable_output: 05-knowledge/results/berggren_transplantable_word_language_harmonic_probe_20260816.out
variable_script_sha256: 86a7a891ada981b61b3ebc4b02cc1b02ba29fafbd05129e5a32a91fa4c134874
variable_output_sha256: 8690f9254744ac2fae5ad700e1932f45a7b996cb88be5153d864e6e9f1a02a9a
fixed_script: 04-computation/berggren_fixed_drift_word_language_parity_probe_20260816.py
fixed_output: 05-knowledge/results/berggren_fixed_drift_word_language_parity_probe_20260816.out
fixed_script_sha256: f503aa44b7287ea200589093c17c3e68f01e456b314f84b607c6695359249d4a
fixed_output_sha256: c8c87fa66bae84d1165dbdd00312e25d22fdae17841c5cc0d45aa205b7433d9b
hash_basis: LF-normalized bytes
---

# THM-3497 -- three quantified Berggren/T4 calibration problems

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

Equal four-point cardinalities do not supply a common connection.  They do,
however, produce two different regular languages once the calibration is
allowed to depend on the whole ancestry word.  This theorem keeps those three
quantifier levels separate and computes each one exactly.

## 1. Typed actions and the multiplication convention

Use THM-3339's parameter matrices

```text
A=[0 1;-1 2],       B=[0 1;1 2],       C=[1 0;2 1].       (1)
```

On

```text
P1(F_3)=([1:0],[0:1],[1:1],[1:2])
```

their projective permutations are

```text
A_3=(0 1 3), fixing 2;
B_3=(0 1 3 2);
C_3=(0 3 2), fixing 1.                                    (2)
```

On the three nonzero directions of `V4=F_2^2`, their true linear parts are

```text
Abar=Bbar=J=[0 1;1 0],              Cbar=I.               (3)
```

The word convention is chronological.  If `w=h_1...h_n` records successive
children from root to leaf, then

```text
rho_3(w)=h_(n,3)...h_(1,3),
L_2(w)=hbar_n...hbar_1.                                   (4)
```

Thus the displayed word `BA` means apply `B` and then `A`, so its matrix is
`AB`.  The first nontrivial Fibonacci ray fixes this convention sharply:

```text
AB(1,2)=(5,8),              BA(1,2)=(3,8).                (5)
```

Only the first value is THM-3339's `u_5`.  The reverse convention therefore
fails before any asymptotic or automaton calculation.

There are three distinct predicates below.

1. A static point gauge, or a six-frame family of gauges, must intertwine a
   specified generator while retaining (3).
2. In the **variable-translation language**, each accepted word may choose
   both a new point gauge and a new affine translation.
3. In the **fixed-drift language**, each accepted word may choose a new point
   gauge, but its affine target is the composite of THM-3339's already fixed
   lifts.

None says that one gauge intertwines the entire branch monoid.  In fact the
generator calculation below rules that out.

## 2. Static and six-frame calibration

For a generator `h`, a static calibration is a bijection
`f:P1(F_3)->V4` and a translation `t_h` satisfying

```text
f h_3 f^(-1)(u)=hbar u+t_h.                               (6)
```

The source cycle types in the order `A,B,C` are

```text
1*3,                         4,                         1*3. (7)
```

For linear part `J`, the four affine lifts have type `1^2*2` or `4`; for
linear part `I`, they have type `1^4` or `2^2`.  Hence `A` and `C` have no
solution of (6).  The two lawful four-cycle translations for `B` each have
four conjugating gauges, so the exact counts, after allowing `t_h`, are

```text
A:0,                         B:8,                         C:0. (8)
```

For THM-3339's pinned drift

```text
Atilde(u)=Ju+p,       Btilde(u)=Ju+p,       Ctilde(u)=u+p, (9)
```

the counts are `(0,4,0)`.  In the point order `(0,p,q,r)`, `B_3=Btilde`, so
the identity gauge is one of the four positive controls.

Let `Ord(V4*)` be the six frames of the three nonzero directions.  A
frame-dependent family asks for six point gauges `f_pi` satisfying

```text
f_(hbar pi) h_3=(hbar(-)+t_h) f_pi.                       (10)
```

On the resulting 24 states, the source/target cycle types are

```text
source A: 2^3 6^3;       target A: 2^12 or 4^6;
source B: 4^6;           target B: 2^12 or 4^6;
source C: 1^6 3^6;       target C: 1^24 or 2^12.          (11)
```

Thus `A` and `C` still have no family.  For `B`, the exact counts by
translation `(0,p,q,r)` are

```text
(0,512,512,0).                                             (12)
```

Indeed `J` has three two-cycles on the six frames.  On each return orbit the
source and target squares are double transpositions, with eight conjugating
seed gauges; the three choices are independent, giving `8^3=512` for each
lawful translation.  This proves the advertised positive counts `4/8/512`
without confusing their different quantifiers.

## 3. Variable translation: the twelve-state language

For a word `w`, call it variable-translation calibratable when there are a
word-dependent bijection `f_w` and translation `t_w` such that

```text
f_w rho_3(w) f_w^(-1)(u)=L_2(w)u+t_w.                    (13)
```

Let

```text
q:S4 -> S3                                                  (14)
```

be the action on the three perfect matchings of `K4`.  Its kernel is the
normal `V4`.  With `r` a three-cycle and `s` a reflection in `S3`, (2)--(3)
give

```text
q(A)=q(C)=r,                 q(B)=s,
epsilon(A)=epsilon(B)=1,     epsilon(C)=0.                (15)
```

Here `epsilon` records whether the true composite linear part is `J`.  The
word map

```text
Phi_var(w)=(q(rho_3(w)),epsilon(w)) in S3 x C2             (16)
```

is onto: `C` supplies `(r,0)`, `A` supplies `(r,1)`, their quotient supplies
the pure bit, and `r,s` generate `S3`.

The exact acceptance criterion is

```text
w satisfies (13)
iff Phi_var(w)=(1,0)
 or Phi_var(w)=(reflection,1).                             (17)
```

For bit zero, the affine targets are precisely the identity and three double
transpositions, namely `ker(q)`.  For bit one, the four translates of `J`
are two transpositions and two four-cycles, precisely the two `S4` conjugacy
types lying over the three reflections of `S3`.  This proves both directions
of (17), not just necessity.

The accepting set has four of the twelve states.  Its right stabilizer in
`S3 x C2` is trivial: a bit flip cannot preserve the bit-fibre sizes `1` and
`3`, and with bit zero the unique accepted identity forces the `S3` component
to be the identity.  Hence all twelve states are Myhill--Nerode inequivalent.
The minimal DFA and its syntactic transition group are both the regular
`S3 x C2` action.  Acceptance is not multiplicative (`B` and `BC` pass but
`BBC` does not), so this is a language, not a submonoid.

If `a_n` counts accepted words of length `n`, then

```text
a_0,a_1,...=1,1,3,11,25,81,251,715,2193,6593,...,         (18)

a_n=2a_(n-1)+2a_(n-2)+6a_(n-3)-9a_(n-4),                 (19)

sum_(n>=0) a_n x^n
  =(1-x-x^2-3x^3)/((1-x)(1-3x)(1+2x+3x^2)).              (20)
```

Writing `c_0=1`, `c_1=-1`, and
`c_n=-2c_(n-1)-3c_(n-2)`, one has the exact closed form

```text
a_n=(3^n+1+c_n)/3,
c_n=Re((-1+i sqrt(2))^n).                                 (21)
```

Therefore `a_n/3^n=1/3+O(3^(-n/2))`.

Enumerate all words in shortlex order, with the empty word indexed by `1`,
and let `S_var` be the accepted index set.  The full 12-state transfer matrix
is annihilated by

```text
(lambda-3)(lambda-1)(lambda+1)(lambda^2+2lambda+3).       (22)
```

The same `O(3^(d/2))` suffix error therefore holds from every prefix state.
A lexicographic initial segment is a union of at most two full ternary
cylinders at each remaining depth.  Summing the geometric errors, together
with the complete earlier levels, gives

```text
|S_var intersect [1,N]|=N/3+O(sqrt(N)).                  (23)
```

Partial summation yields

```text
sum_(m<=N, m in S_var) 1/m
   =(1/3)log N+C_var+O(N^(-1/2)).                         (24)
```

## 4. Fixed THM-3339 drift: the 192-state language

In the pinned order `(0,p,q,r)`, write

```text
G(u)=Ju+p=(0 1 3 2),             T(u)=u+p=(0 1)(2 3).    (25)
```

Then `G^4=T^2=1` and `TGT=G^(-1)`, so `<G,T>=D4`.  Let
`rho_fix(w)` be the chronological composite of the fixed targets

```text
A -> G,                  B -> G,                  C -> T. (26)
```

Call `w` fixed-drift calibratable when a word-dependent point gauge satisfies

```text
f_w rho_3(w) f_w^(-1)=rho_fix(w).                         (27)
```

Equivalently, the two permutations in (27) have the same cycle type.  The
paired generators generate all of `S4 x D4`.  This follows constructively:
chronological words `AAAA` and `AAAB` supply `(A_3,1)` and `(B_3,1)`, while
`AAA` and `CCC` supply `(1,G^3)` and `(1,T)`.  The first pair generates `S4`
and the second `D4`, so the image has order `192`.

The accepting subset has trivial right stabilizer in this 192-element group,
as an exact `34*192` census verifies.  Thus all 192 states are inequivalent,
and the minimal DFA and syntactic transition group are the regular
`S4 x D4` action.

There is a character `chi:D4->C2` with

```text
chi(G)=chi(T)=1.                                          (28)
```

Every letter flips `chi`, so length parity confines the walk to one of two
96-state chambers.  In normal forms `G^iT^j`, the even chamber has target
types `1^4`, `2^2`, and twice `1^2*2`; the corresponding `S4` class sizes
sum to

```text
1+3+2*6=16.                                               (29)
```

The odd chamber has twice type `4` and twice type `2^2`, giving

```text
2*6+2*3=18.                                               (30)
```

If `b_n` denotes the fixed-drift count, its first terms are

```text
1,1,1,4,13,46,113,421,1121,3667,9801,33166,... .         (31)
```

Its reduced generating function has numerator

```text
1+x-6x^2-3x^3-x^4+11x^5-56x^6-45x^8-99x^9
 -162x^10+243x^11-243x^12-729x^13                       (32)
```

and denominator

```text
(x-1)(x+1)(3x-1)(3x+1)(3x^2+1)
 (3x^2-2x+1)(3x^2+2x+1)(9x^4-2x^2+1).                  (33)
```

Equivalently, its minimal order-fourteen recurrence is

```text
b_n=7b_(n-2)+7b_(n-4)+71b_(n-6)+213b_(n-8)
    +189b_(n-10)+1701b_(n-12)-2187b_(n-14).              (34)
```

The numerator and denominator in (32)--(33) are coprime.  Exact residuals
vanish for more than 192 consecutive indices, and Cayley--Hamilton for the
192-state transfer operator propagates (34) to every `n>=14`; coprimality
proves minimality.

The normalized polar coefficients at the two dominant poles are

```text
lim_(x->1/3)  (1-3x) sum b_n x^n = 17/96,
lim_(x->-1/3) (1+3x) sum b_n x^n = -1/96.                (35)
```

They are not analytic residues.  Every remaining reciprocal pole has
modulus at most `sqrt(3)`, so

```text
b_n=(17/96)3^n-(1/96)(-3)^n+O(3^(n/2)).                 (36)
```

Consequently

```text
b_(2m)/3^(2m) -> 1/6,
b_(2m+1)/3^(2m+1) -> 3/16.                               (37)
```

At complete shortlex level boundaries, geometric summation gives two
different natural-density limits,

```text
11/64                    on even terminal levels,
35/192                   on odd terminal levels.          (38)
```

Thus the fixed-drift shortlex set `S_fix` has no natural density.

## 5. Why the fixed harmonic density nevertheless exists

The level totals (37) alone do not control harmonic weights inside a
shortlex level.  The needed cylinder-uniform statement follows from a second
finite-group gate.

Indeed, at an even level, placing the same asymptotic fraction `1/6` entirely
at the front or back of the level would give limiting masses `log(4/3)` or
`log(9/8)`, rather than `(1/6)log 3`.  Thus no argument from (37) alone can
prove the claimed harmonic coefficient.

Group the walk into ordered pairs of letters.  The nine two-letter increments
generate all `96` elements of `ker(chi)`, so the two-step walk is irreducible
on either parity chamber.  The paired generator orders are

```text
ord(A_3,G)=12,             ord(B_3,G)=4,
ord(C_3,T)=6.                                             (39)
```

Hence `B^4` and `C^6` are identity returns of two-step lengths `2` and `3`.
Their gcd is one, so the two-step walk is aperiodic.  Because every transition
is an average of permutations, the stationary measure on each chamber is
uniform.  Irreducibility and aperiodicity therefore give uniform convergence
to that measure from every prefix state; an optional final letter handles odd
suffix length.

At level `n`, the shortlex block starts at

```text
s_n=(3^n+1)/2.                                            (40)
```

A fixed word prefix cuts the lexicographic ranks into a triadic cylinder.
The mixing result says that within every such fixed cylinder, the accepted
proportion tends to `1/6` on even terminal levels and `3/16` on odd terminal
levels.  Approximating the continuous weight

```text
t |-> 1/(1/2+t),                 0<=t<=1,                 (41)
```

by functions constant on fixed-depth triadic cylinders gives the complete
level harmonic masses

```text
sum_(m in accepted level n) 1/m
 -> (1/6)log 3              for even n,
 -> (3/16)log 3             for odd n.                   (42)
```

Cesaro averaging of (42) over pairs of levels gives the mean

```text
(1/2)(1/6+3/16)=17/96.                                   (43)
```

There are `log N/log 3+O(1)` completed levels below `N`.  Any unfinished last
level lies in an integer interval of endpoint ratio at most three, so even
its total harmonic mass is `O(1)`.  Therefore

```text
sum_(m<=N, m in S_fix) 1/m=(17/96)log N+o(log N).         (44)
```

This proves logarithmic density `17/96` without assuming the nonexistent
natural density and identifies the cylinder-mixing input on which it depends.

## 6. The two Fibonacci clocks

THM-3339's three root-to-child rays are

```text
(BA)^r,                  A(BA)^r,                  C(BC)^r. (45)
```

Applying (17) gives

```text
(BA)^r passes       iff r is even;
A(BA)^r passes      iff r is odd;
C(BC)^r passes      iff r is odd.                         (46)
```

These are exact power laws in the finite quotient, not extrapolations from a
bounded ray prefix.

Thus the variable-translation Fibonacci indices satisfy exactly

```text
k mod 6 in {0,1,2}.                                      (47)
```

Applying the fixed-drift criterion (27) instead gives

```text
(BA)^r passes       iff r=0 mod 4;
A(BA)^r passes      iff r=3 mod 4;
C(BC)^r passes      iff r=3 mod 4,                        (48)
```

Again the period-four calculation is an exact power computation in the paired
finite group.

or exactly

```text
k mod 12 in {0,1,2}.                                     (49)
```

The corresponding reciprocal-index sums have logarithmic coefficients
`1/2` and `1/4`.  Reciprocal Fibonacci values and reciprocal triple
coordinates instead converge exponentially; no geometric mass statement is
being inferred from the index clock.

## 7. Boundary

The theorem classifies finite calibration predicates on Berggren ancestry
words.  It does not construct one branch-monoid-equivariant point gauge, a
physical tournament current, an LRC word-current, a D5 cohomology class, or a
Jacobian invariant.  The perfect-matching quotient is intrinsic and useful,
but it deliberately forgets the actual word-dependent conjugating gauge.

QED.
