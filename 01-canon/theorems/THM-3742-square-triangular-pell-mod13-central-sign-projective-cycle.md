---
id: THM-3742
title: "Square-triangular Pell mod-thirteen central-sign and stereographic projective cycles"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The integral
  half-step 1+sqrt(2) interlaces the square-triangular norm-one family and
  the T_n=2T_m norm-minus-one family.  Modulo thirteen its even step is a
  sharp C14 clock on every nonzero norm fibre.  Ordinary vector
  projectivization kills the central sign and gives two quadratic-character
  C7 orbits; a different, stereographic parametrization of the norm-one
  conic is bijective with P1(F13) and conjugates the same update to one
  Singer-type C14 cycle.  Every nonzero scalar linear observer has image
  exactly seven or eight, with an exact dual-norm criterion.  The companion
  exhausts all fibres, observers, and projective states.  The closest LRC
  offset shell lies on the wrong norm torsor, and the Pell matrix is excluded
  from every constant-right exposed Cohn closure.  No LRC(14), planar-JC, or
  constant-Cohn closure follows.
source: root + jc_arithmetic_archaeology + lrc_khinchin_archaeology / 2026-08-23
audit: >
  PASS.  An independent derivation reconciled the apparent C7/C14 conflict,
  checked the stereographic inverse including its removable t=6 boundary,
  computed N^7 and N^14, recovered the Singer conjugacy, and corrected the
  raw-observer coefficient norm to a^2-5b^2.  A second audit separated the
  square and nonsquare direct-projective torsors and the THM-3713 five-offset
  near-match.  Normal and optimized exact companions are byte-identical.
depends_on:
  - THM-3335-square-triangular-pell-markov-pythagorean-selector
  - THM-2619-psl2f13-seven-edge-norm-minimal-projective-kernel-and-retained-target-obstruction
  - THM-3726-automorphic-cohn-constant-sl2-orbit-broughton-classification
  - THM-3736-automorphic-cohn-complete-constant-sl2-polynomial-exposure-classification
related:
  - THM-3234-singer-owner-compactification-and-pointed-heisenberg-carrier-gate
  - THM-3713-lrc-orbitwise-deep-offset-detector-and-successor-quotient-hostile
  - THM-778-centered-christoffel-endpoint-skew-product
  - THM-3570-universal-pell-conic-target-graph-factor-compiler
  - THM-3573-polynomial-target-graph-pell-parameter-descent-classification
  - THM-3575-coprime-pell-target-coordinate-euler-quotient-no-go
  - HYP-2453-triangular-tower-moment-bridge
  - HYP-2529-triangular-tower-repunit-carrier-tournament
script: 04-computation/square_triangular_pell_mod13_projective_cycle_thm3742.py
output: 05-knowledge/results/square_triangular_pell_mod13_projective_cycle_thm3742.out
script_sha256: 3232e2185a8b05af3caddb0a9ecc43d61030b76acd2969d599ccaa4642570299
output_sha256: c4623ddab191b2a1e1106adf2197332721b165729a9e61e6792158f455a6b19f
hash_basis: raw LF bytes
---

# THM-3742 -- a signed Pell clock, and exactly what its quotients forget

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem
does not introduce the square-triangular recurrence, the projective
seven-clock, or the Singer cycle.  Its repository-new payload is their exact
typed connection: an integral signed Pell lift, the norm-character and
central-sign invoices lost by ordinary projectivization, a stereographic
fourteen-state repair, and sharp scalar-observer and target-specific
boundaries.

## 1. Inheritance pass and notation

The closest proved arithmetic mechanism is
[THM-3335](THM-3335-square-triangular-pell-markov-pythagorean-selector.md),
whose square-triangular update is

```text
M=[3,8;1,3]                                                   (1)
```

in coordinates `(x,q)` on `x^2-8q^2=1`.  The closest finite clock is
[THM-2619](THM-2619-psl2f13-seven-edge-norm-minimal-projective-kernel-and-retained-target-obstruction.md),
which already proves the nonsplit `PSL_2(F_13)` seven-cycle and its retained
target-transition obstruction.  The canonical hostile is therefore not a
failed orbit calculation; it is the loss of the norm torsor, central sign,
physical root, owner, and target when the orbit is projectivized.

The corrected near miss is MISTAKE-369: two Pell families in the same
quadratic field, or with the same recurrence, need not lie in the same typed
norm/parity torsor.  The least-used sidecar is the signed norm value before
projectivization.  It will distinguish every lawful statement below.

Put

```text
Q(x,s)=x^2-2s^2,
S_0=[1,2;1,1],
A=S_0^2=[3,4;2,3],
D=[1,0;0,2],               M=D^(-1) A D=[3,8;1,3].    (2)
```

The symbol `S_0` is reserved for the integral Pell half-step.  It is not the
Singer matrix in Section 5.

## 2. The two classical triangular families are alternate half-steps

### Theorem 2.1 (signed Pell interlacing)

For every `k>=0`, write

```text
(x_k,s_k)=A^k(1,0),          (u_k,v_k)=A^k(1,1).       (3)
```

Then there are unique nonnegative integers `n_k,q_k,N_k,m_k` with

```text
x_k=2n_k+1,       s_k=2q_k,
u_k=2N_k+1,       v_k=2m_k+1,                          (4)
```

and

```text
T_(n_k)=q_k^2,              T_(N_k)=2T_(m_k).          (5)
```

These are the complete nonnegative solution branches of the two displayed
equations, with their shared seed included.  Moreover

```text
A^k(1,1)=S_0 A^k(1,0),                                  (6)
```

so the norm `+1` and norm `-1` families are precisely the even and odd powers
of the same half-unit `1+sqrt(2)`:

```text
S_0^(2k)(1,0)=A^k(1,0),
S_0^(2k+1)(1,0)=A^k(1,1).                              (7)
```

The first two scalar projections are

```text
q_k: 0,1,6,35,204,1189,...,
v_k: 1,5,29,169,985,5741,... .                         (8)
```

The second sequence is the odd companion appearing in HYP-2529.  Retaining
only it preserves the recurrence but loses the first coordinate, depth, and
norm-coset address.

*Proof.*  Direct multiplication gives

```text
S_0^T diag(1,-2) S_0=-diag(1,-2),
A=S_0^2,
A^T diag(1,-2) A=diag(1,-2).                          (9)
```

Thus `(3)` has norms `+1` and `-1`, respectively, and its parity is preserved.
Substituting `(4)` gives exactly

```text
(2n+1)^2-8q^2=1       iff T_n=q^2,
(2N+1)^2-2(2m+1)^2=-1 iff T_N=2T_m.                  (10)
```

The matrices in `(2)` commute, and `S_0(1,0)=(1,1)`, proving `(6)--(7)`.

For completeness, the standard inverse-unit descent proves exhaustion.  On
the norm-one branch, THM-3335 descends every positive nonseed solution by
`A^(-1)`.  On the norm-minus-one branch, a solution other than `(1,1)` has
`v>=5`; then

```text
(u',v')=(3u-4v,-2u+3v)                                (11)
```

is again positive of norm `-1` and has `0<v'<v`.  Indeed
`u/v=sqrt(2-v^(-2))` lies strictly between `4/3` and `3/2`, while `u>v`.
Descent terminates at `(1,1)`.  This proves the two complete branches.  QED

The word *interlaces* is typed: `S_0` reverses the norm.  It does not identify
the two norm fibres as groups.  The norm-one branch is a unit group; the
norm-minus-one branch is only its torsor until a basepoint is chosen.

## 3. Modulo thirteen: one vector C14 on every nonzero norm fibre

Work in `F=F_13` and put

```text
C_e={(x,q) in F^2:x^2-8q^2=e},          e in F*.       (12)
```

Both `2` and `8` are nonsquares in `F`.  In
`E=F[sqrt(8)]`, the matrix `M` is multiplication by

```text
z=3+sqrt(8),                  Norm(z)=1.               (13)
```

Exact powers are

```text
M^7=-I,                       M^14=I.                  (14)
```

Consequently `z` generates the fourteen-element norm-one subgroup of
`E*`.  The norm map `E* -> F*` is onto and every nonzero fibre has fourteen
points.  Multiplication by `z` is therefore sharply transitive on every
`C_e`:

```text
|C_e|=14,                       C_e is one M-orbit.     (15)
```

After seven steps every vector is its negative.  This central sign is a real
coordinate of the vector orbit, not cosmetic duplication.

### Ordinary vector projectivization

The map

```text
pi_lin:C_e -> P^1(F),             (x,q) |-> [x:q]      (16)
```

identifies `(x,q)` with `(-x,-q)`.  Each fourteen-cycle therefore becomes a
seven-cycle.  The quadratic character of the norm is well-defined on a line,
because scalar multiplication changes the norm by a square.  It gives the
complete two-orbit split.  In the affine slope `q/x`, the two cycles are

```text
O_square   ={0,4,5,6,7,8,9}       ={0,+/-4,+/-5,+/-6},
O_nonsquare={infinity,1,2,3,10,11,12}
            ={infinity,+/-1,+/-2,+/-3}.               (17)
```

Since `-1=5^2 mod 13`, both integral families from Section 2—norm `+1` and
norm `-1`—land in `O_square`.  Moving to `O_nonsquare` changes the torsor and
is not licensed by their common Pell recurrence.

The companion enumerates all `169` vectors, all twelve nonzero norm fibres,
and the central-sign pairing in each orbit.

## 4. Stereographic completion retains all fourteen conic states

There is a second appearance of `P^1(F)` which must not be confused with
`(16)`.  On the norm-one conic define

```text
theta(x,q)=q/(x+1),              theta(-1,0)=infinity. (18)
```

Because `8` is nonsquare, the finite inverse has no pole:

```text
theta^(-1)(t)=((1+8t^2)/(1-8t^2), 2t/(1-8t^2)),
theta^(-1)(infinity)=(-1,0).                            (19)
```

Thus `theta:C_1 -> P^1(F)` is a bijective stereographic parametrization of
the genus-zero conic.  It is not induced by a linear quotient of the ambient
two-dimensional vector space.

Substitution in `(18)` gives

```text
theta M theta^(-1)(t)=(4t+1)/(8t+4),
N=[4,1;8,4],                         det N=8.           (20)
```

At `t=6` the cancelled affine expression has a removable `0/0`; the projective
matrix sends it lawfully to `infinity`.  Starting at zero, `N` has the single
cycle

```text
0,10,9,1,8,11,6,infinity,7,2,5,12,4,3.               (21)
```

Its exact half- and full-period powers are

```text
N^7=[0,12;5,0],                 N^14=8I.              (22)
```

Projectively, the half-period is

```text
t |-> 5/t=1/(8t).                                      (23)
```

Equation `(19)` shows that `(23)` is exactly
`(x,q)|->(-x,-q)`.  Hence there is no C7/C14 contradiction:

- `pi_lin` is two-to-one and kills the antipode;
- `theta` is bijective and retains it as the nontrivial involution `(23)`.

This distinction is the central validity gate of the theorem.

## 5. Both finite clocks were already present, but without the Pell typing

Let

```text
g_6=[0,1;-1,6],              P=[11,5;5,0].             (24)
```

Then, modulo thirteen,

```text
det P=1,                     P g_6=M P.                (25)
```

Thus the direct projective action of `M` is precisely an `SL_2` lift of the
THM-2619 nonsplit seven-clock.  THM-2619 already proves that this projective
clock does not by itself retain the target transition needed by LRC(14).

For the stereographic clock, use the distinct Singer matrix and conjugator

```text
G=[1,4;2,1],                 J_0=[0,1;2,0].            (26)
```

Exact multiplication gives

```text
J_0 G^5=6N J_0,
[N]=[J_0][G]^5[J_0]^(-1) in PGL_2(F).                 (27)
```

So `(21)` is a Singer-cycle generator of the kind already present in
THM-3234.  The new information is not “there exists a projective fourteen-
cycle”; it is that the classical integral square-triangular unit maps to it
through `(18)`, together with the precise sign and norm invoices.

## 6. Every scalar linear observer sees exactly seven or eight residues

Fix `e!=0` and a nonzero raw linear form

```text
L_(a,b)(x,q)=a x+b q.                                  (28)
```

Its dual norm is

```text
c(a,b)=a^2-8^(-1)b^2=a^2-5b^2.                       (29)
```

This is never zero for `(a,b)!=(0,0)`, because `5` is nonsquare.  With the
quadratic character extended by `chi(0)=0`, every output fibre has the exact
size

```text
#{(x,q) in C_e:L_(a,b)(x,q)=y}
  =1-chi(y^2-e c(a,b)).                                (30)
```

Consequently

```text
|L_(a,b)(C_e)|
 =7+(1+chi(e c(a,b)))/2
 ={8 if chi(e c(a,b))=+1,
   7 if chi(e c(a,b))=-1}.                            (31)
```

In the size-eight case there are two tangent singleton fibres and six double
fibres.  In the size-seven case all seven nonempty fibres are double.  On
each `C_e`, exactly `84` of the `168` nonzero coefficient vectors give each
size.

*Proof.*  If `a!=0`, eliminate `x=(y-bq)/a`.  The discriminant of the resulting
quadratic in `q` is

```text
4a^2[8y^2+e(b^2-8a^2)]
 =32a^2[y^2-e(a^2-5b^2)].                              (32)
```

Because `8` is nonsquare, the number of roots is exactly `(30)`.  When `a=0`,
substitution of `q=y/b` gives the same formula directly.  Counting the
values for which `(30)` is positive proves `(31)`; the tangent alternatives
give the stated multiplicities.  The `84/84` split is an exhaustive finite
classification over the fourteen projective coefficient directions, replayed
for every nonzero scalar and every `e`.  QED

The coefficient convention matters.  For raw `ax+bq`, the correct expression
is `a^2-5b^2`, not `a^2-8b^2`.  For example `L=x+q` has eight images because
`1-5=9` is square, although `1-8=6` is nonsquare.  If instead the observer is
normalized as `ax+8bq`, its dual norm becomes `a^2-8b^2`.

### Square and triangular scalar palettes are even smaller

Modulo thirteen,

```text
squares     ={0,1,3,4,9,10,12},
triangulars ={0,1,2,3,6,8,10}.                        (33)
```

Their union misses `{5,7,11}`.  Thus neither an untyped linear compression nor
the raw square/triangular scalar palette can saturate all twelve nonzero LRC
colours.  The full rational denominator in `theta` is doing real work.

## 7. Exact LRC resonance: the right shell on the wrong torsor

The five active offsets in THM-3713's non-cover control satisfy the literal
set identity

```text
{1,2,-3,-2,-1}=O_nonsquare \ {infinity,3}.             (34)
```

This is a sharper connection than a numerical resemblance: the current five
deep offsets occupy five of the seven rays of one Pell projective orbit.
But the square-triangular and `T_n=2T_m` families live in `O_square`, and a
change of affine chart can move the two deleted labels.  Equation `(34)` is
therefore a **PROVED exact near-match and a proved typing obstruction**, not
an LRC carrier theorem.

There is also an exact polarity match.  Put

```text
J_1=diag(1,-1),
R=M^(-1)J_1=[3,8;12,10].                               (35)
```

Then

```text
R^2=I,                       RMR=M^(-1).               (36)
```

If `v_r=M^r(1,0)` indexes the norm-one vector orbit, then

```text
R v_r=v_(-1-r)               (indices mod 14),         (37)
```

and the same affine reflection holds modulo seven after the sign quotient.
This is exactly the formal reflection `r|->-1-r` in THM-778's centred
Christoffel token calculus.  What is not preserved is decisive: THM-778
updates the actual event owner by an inverse-speed increment, whereas the Pell
clock advances uniformly.  No canonical map identifies those physical steps.

A lawful future test must retain at least

```text
(norm-character torsor, antipodal sign, affine chart and deleted rays,
 physical root/orientation, owner, inverse step, midpoint phase and ties,
 temporal word).                                      (38)
```

The cheapest decisive experiment is to place only the fourteen dihedral
gauges of `O_nonsquare` over the THM-3713/3718 exact ledger, demand that the
two deleted rays be the physical root and excluded target, and then test the
actual owner increment and THM-2305 temporal word inside every surviving
fibre.

## 8. The same constant matrix cannot enter the exposed Cohn route

Now interpret the integer matrix `M` over a characteristic-zero field as the
constant right factor `R=[a,b;c,d]` in THM-3726/3736.  Here

```text
(a,b,c,d)=(3,8,1,3),             epsilon=2a-d=3,
epsilon^2-1=8!=0.                                    (39)
```

The raw THM-3726 constant lower equations force `h=0` and then `0=8`; the
upper equations already force `0=-1`.  Hence there is no constant exposed
closure.  THM-3736 proves that a nonconstant lower closure requires `c=0`
and a nonconstant upper closure requires `b=0`.  Both fail.  Therefore:

```text
No polynomial h closes either exposed row of M_0 M.    (40)
```

This is a complete nonentry theorem for the constant Pell matrix, not a
bounded search.  Any useful Pell route into the planar Cohn program would
need a genuinely polynomial right factor or a longer nonconstant word.

The other Pell syntax in THM-3570/3573/3575 is also distinct.  Setting
`b=0,a=-1/6` in THM-3570 reproduces the bare conic
`W^2-2x^2=1`, but collapses the target base; its involution `q|->-2/q` is
branch exchange, not Pell time.  The mod-thirteen orbit preserves none of the
required target packet

```text
(a,b,phi,q/H,P,Q, source coordinate, collision, Euler quotient).          (41)
```

Thus “quadratic norm equation plus involution” is the full preserved
predicate.  No target-graph, polynomial-descent, Keller, collision, or
Euler-quotient conclusion transfers.

## 9. Exact audit and scope

Reproduce with

```bash
python3 -B 04-computation/square_triangular_pell_mod13_projective_cycle_thm3742.py
python3 -B -O 04-computation/square_triangular_pell_mod13_projective_cycle_thm3742.py
```

The companion checks the integral quadratic-form identities and depths
`0..12`; all `169` finite-field vectors and twelve nonzero norm fibres; both
ordinary projective cycles; all fourteen stereographic states; the antipodal,
Singer, `g_6`, and dihedral identities; every one of the
`12*168*13` scalar-observer output fibres; the residue palettes; and the raw
Cohn gates.  Its `30,503` exact checks are unchanged by `-O`.

The theorem supplies an exact address and a precise list of lost coordinates.
It supplies no physical LRC owner/root/word transport, no semantic arrival,
no LRC(14) proof, no planar Keller map, and no planar-Jacobian conclusion.
