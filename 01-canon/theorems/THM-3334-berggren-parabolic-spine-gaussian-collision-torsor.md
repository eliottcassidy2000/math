---
id: THM-3334
title: "Berggren fixed-cusp Farey fan, parabolic spine, complete descendant-angle law, and Gaussian collision torsors"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED, with cited
  antecedents separated from repo deductions.  The consecutive-parameter
  Berggren U-spine is the Gaussian-square image of a chain of Farey faces
  around the fixed cusp (1,1).  Its Q-values are one scalar evaluated on the
  cusp and chain, not one orbit containing the exceptional 2.  Every
  Berggren descendant triangle has one fixed angle with cosine 8/9 and is
  acute/U-obtuse/D-obtuse according to the exact 8/9 and 9/8 leg-ratio
  walls.  Equal-hypotenuse parents form a Gaussian factor-choice torsor
  F_2^r/<1>; along the U-spine these Boolean fibres are unbounded by a
  discriminant -4 CRT construction.  The first affine V4 fibre is c=1105,
  whose four ancestry words and three prime-XOR perfect matchings are
  explicit.  XOR supplies edge colours, not a canonical tournament.
source: codex-2026-08-12-berggren-spine-collision-torsor
audit: >
  Two independent read-only audits rederived the matrix, Farey, plane,
  multiplicity, fixed-cusp, angle, and boundary formulas; one separately
  traced the primary literature and corrected the fixed-hypotenuse scope,
  Price-branch wording, figure attribution, and exceptional-term status.
depends_on:
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2632-farey-v4-theta-channel-and-hurwitz-crt-parity-sidecar
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
  - THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone
related:
  - THM-1310-jacobian-counterexample-fiber-geometry-S3-resolvent-jelonek
  - THM-3536-berggren-angle-languages-signed-c4-and-harmonic-density
script: 04-computation/berggren_parabolic_gaussian_torsor_thm3334.py
output: 05-knowledge/results/berggren_parabolic_gaussian_torsor_thm3334.out
script_sha256: bc201467b685fc0ba3654b17726ddbaae403a4389cbcf7259613e47d3366063a
output_sha256: 3c2416ecb78c1b25a744a1cc70b1ed2eafa9fcb0552a628f298590ce2b23403c
hash_basis: working-tree bytes
---

# THM-3334 -- one fixed cusp, an unbounded Boolean family of planes

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem
combines cited published antecedents with new elementary deductions.  No
literature-priority claim is made for the synthesis or the angle law.

## 1. Inheritance, conventions, and the fixed-cusp theorem

[THM-2596](THM-2596-modular-free-factor-farey-gram-owner-cocycle.md) proves
that the three Berggren branches are infinite-order `PGL_2(Z)` reduction
maps, not one `C_3` action.  [THM-3333](THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md)
uses the ordered raw Gaussian-square lift

```text
Phi(m,r)=(m^2-r^2, 2mr, m^2+r^2)                         (1)
```

and proves

```text
<Phi(u),Phi(v)>_L=2 det(u,v)^2.                          (2)
```

Thus unimodular spinor pairs, not unordered primitive triangles, are the
faithful Farey carrier.

Use the Berggren matrices

```text
U=[1 -2 2;  2 -1 2;  2 -2 3],
A=[1  2 2;  2  1 2;  2  2 3],
D=[-1 2 2; -2  1 2; -2  2 3],                           (3)
```

and the parameter matrix

```text
g=[2 -1; 1 0].                                           (4)
```

For every integer spinor `u`, direct expansion gives the intertwining law

```text
U Phi(u)=Phi(gu).                                         (5)
```

Put

```text
s=(1,1),                 u_n=(n+2,n+1),        n>=0.      (6)
```

Then

```text
gs=s,                     (g-I)^2=0,
gu_n=u_(n+1),             u_(n+1)=u_n+s,                 (7)
det(s,u_n)=-1,            det(u_n,u_(n+1))=1.            (8)
```

The third determinant in the face also has absolute value one.  Hence

```text
{[s],[u_n],[u_(n+1)]}
```

is a Farey triangle for every `n`.  The U-spine is therefore more precisely
the Gaussian-square image of a chain of Farey faces in the fan about one
fixed cusp.  This does not identify the full Berggren tree with the Farey
graph: a general Berggren edge can be a longer Farey prefix, as THM-3333
shows.

Writing `P_n=Phi(u_n)`, (5)--(7) give

```text
P_n=U^n(3,4,5)^T
   =(2n+3, 2n^2+6n+4, 2n^2+6n+5).                       (9)
```

If `P_n=(a_n,b_n,c_n)`, define

```text
Q_n=a_n^2+2.
```

Then

```text
Q_n=a_n^2+2=2b_n+3=2c_n+1=4n^2+12n+11,                 (10)
Q_(n+1)-Q_n=8n+16,                    Delta^2 Q_n=8.     (11)
```

The indexing in the motivating list is `t=n+1`:

```text
q_t=4t(t+1)+3,  t>=1,       giving 11,27,51,83,... .     (12)
```

## 2. Why the growth is quadratic: parabolic horocycles

Let `N=U-I`.  Exact multiplication gives

```text
N^2!=0,                    N^3=0,                        (13)
U^n=I+nN+binom(n,2)N^2.                                  (14)
```

Moreover

```text
U^T diag(1,1,-1) U=diag(1,1,-1).                        (15)
```

Thus `U` is a unipotent integral Lorentz transformation.  Its fixed space is
the line spanned by the primitive null vector

```text
v=(0,1,1).                                                (16)
```

There is a full horocycle statement behind the special spine.  The linear
functional `c-b` is U-invariant.  For every nonzero scalar `delta`, the
intersection of the null cone with `c-b=delta` is

```text
P_delta(a)=
 (a, (a^2-delta^2)/(2delta), (a^2+delta^2)/(2delta)),     (17)
```

and substitution gives

```text
U P_delta(a)=P_delta(a+2delta).                           (18)
```

The positive primitive spine is the integral section `delta=1`, with odd
`a>=3`.  Formula (18), not merely sequence fitting, is the mechanism for the
quadratic coordinate growth and constant second difference.

Projectively,

```text
P_n/c_n -> v.                                             (19)
```

So the fixed null ray is the cusp at infinity of the invariant parabola.

## 3. The exceptional 2 and the degenerate 3

The raw cusp has

```text
Phi(s)=Phi(1,1)=(0,2,2)=2v.                              (20)
```

Its content two is exactly THM-3333's odd/odd spinor parity case.  Define on
spinors

```text
q(m,r)=(m^2-r^2)^2+2.                                    (21)
```

Then

```text
q(s)=2,                    q(u_n)=Q_n.                    (22)
```

Consequently

```text
2,11,27,51,...
```

is one scalar evaluated first at the fixed cusp and then on its adjacent
fan chain.  It is **not** one orbit: `s` is fixed while the `u_n` form a
separate orbit.  The Euclidean identity `||v||^2=2` is also true after the
primitive lattice normalization, but Euclidean norm is not U-invariant and
does not supply an orbit law.

The uniform quadratic progression would instead begin with `3`.  Indeed

```text
U(1,0,1)^T=(3,4,5)^T,                                    (23)
```

and `q_0=3` reconstructs the degenerate null point `(1,0,1)`.  Its three
putative descendants are not distinct, so it is not a positive parent and
has no descendant triangle.

A scalar `Q` reconstructs the displayed parent only under the typing

```text
Q-2 is an odd square,  Q>=11.                             (24)
```

For example `Q=13` is the minimal odd hostile to the unqualified statement
that every odd `Q>=11` reconstructs a spine parent.

## 4. The descendant plane and the complete angle classification

Kőszegyová--Csókási--Hirjak prove that the three descendants of a primitive
parent `P=(a,b,c)` are noncollinear and lie in

```text
2x+2y-3z+c=0;                                             (25)
```

they also give the area formula below.  See their
[published 2024 paper](https://doi.org/10.29020/nybg.ejpam.v17i3.5323).
All formulas in this section also follow directly from (3).

Write the child vertices as `V_U,V_A,V_D`.  Exact subtraction gives

```text
u=V_A-V_U=(4b,2b,4b),
v'=V_D-V_U=(-2a+4b,-4a+2b,-4a+4b),
w=V_D-V_A=(-2a,-4a,-4a).                                 (26)
```

Hence

```text
|u|=6b,
|w|=6a,
|v'|=2 sqrt(9a^2-16ab+9b^2),                             (27)
u x v'=(8ab,8ab,-12ab),
area=2 sqrt(17) ab.                                      (28)
```

Put

```text
L=sqrt(9a^2-16ab+9b^2).
```

The three angle dot products reduce to

```text
(-u).w =32ab,
u.v'   =4b(9b-8a),
v'.w   =4a(9a-8b).                                      (29)
```

Therefore every Berggren descendant triangle satisfies

```text
cos(angle A)=8/9,
cos(angle U)=(9b-8a)/(3L),
cos(angle D)=(9a-8b)/(3L).                               (30)
```

The equality walls would force a positive integer Pythagorean hypotenuse
proportional to `sqrt(145)`, so they never occur.  The complete classification
is

```text
b/a < 8/9          iff the descendant triangle is U-obtuse,
8/9 < b/a < 9/8    iff it is acute,
b/a > 9/8          iff it is D-obtuse.                   (31)
```

Minimal positive controls are respectively

```text
(15,8,17),          (21,20,29),          (3,4,5).         (32)
```

The 2024 paper proves only that the triangle is nonright/nonisosceles and
asks whether it is acute or obtuse.  Equations (30)--(31) answer that stated
question for the Berggren tree; they are repo deductions from the published
vectors, with no global-priority claim.

On the U-spine,

```text
b=(a^2-1)/2,            b/a>9/8 for every odd a>=3,       (33)
```

so every spine descendant triangle is D-obtuse.  Moreover

```text
angle U -> 0,
angle A = arccos(8/9),
angle D -> arccos(-8/9).                                 (34)
```

In terms of `Q`, (25) and (28) become

```text
4x+4y-6z+Q-1=0,
distance from origin=(Q-1)/(2 sqrt(17)),
area=sqrt(17)(Q-3)sqrt(Q-2).                             (35)
```

These `Q` formulas are scoped to the spine representation, not every parent
in the plane.

## 5. Equal-plane fibres are Gaussian factor-choice torsors

For odd `c`, let `X_c` be the set of ordered primitive Pythagorean parents
with hypotenuse `c`.  Tripathi's fixed-integer theorem, as imported by the
2024 paper, states

```text
|X_c| = 2^(omega(c)-1)   if c>=3 and no p=3 mod 4 divides c,
|X_c| = 0                if c=1 or such a prime divides c. (36)
```

See [Tripathi, 2008/09](https://doi.org/10.1080/00150517.2008.12428142).
The prime condition is load-bearing; `c=3` is the minimal hostile to the
unqualified power-of-two formula.

There is a structural refinement of (36).  Suppose

```text
c=product_(j=1)^r p_j^(e_j),               p_j=1 mod 4.  (37)
```

Choose the conjugate Gaussian prime pair

```text
p_j=pi_j conjugate(pi_j).
```

If `z=m+in` has norm `c` and `gcd(m,n)=1`, then unique factorization in
`Z[i]` forces the whole exponent `e_j` to occur on exactly one side:

```text
pi_j^(e_j)  or  conjugate(pi_j)^(e_j).                   (38)
```

Splitting the exponent between both sides would make `p_j` divide both
`m` and `n`.  Thus the raw choices form a torsor for `F_2^r`.  Global complex
conjugation complements all bits and gives the same positive ordered parent.
Consequently

```text
X_c is an affine torsor for
V_c=F_2^r/<(1,...,1)>,                    dim V_c=r-1.    (39)
```

This proves (36) and records the missing coordinate: a plane label retains
`c` but loses the Gaussian factor allocation.

The descendant-plane theorem identifies the set of triangles in (25) with
`X_c`.  Distinct parents in one fibre are pairwise incomparable in the
Berggren tree.  Indeed every child has strictly larger hypotenuse:

```text
c_U-c=2(a-b+c)>0,
c_A-c=2(a+b+c)>0,
c_D-c=2(-a+b+c)>0.                                      (40)
```

Thus an equal-plane fibre is an ancestry antichain, not a tree level or a
set of branch edges.

The smallest global collision is

```text
c=65=8^2+1^2=7^2+4^2,                                  (41)
```

with parents `(63,16,65)` and `(33,56,65)`.  Its scalar plane label is
`Q=131`.  This is THM-3333's norm-65 representation hostile in the present
plane language.

## 6. Discriminant -4 makes the spine fibres unbounded

The U-spine hypotenuse is

```text
c_t=t^2+(t+1)^2=2t^2+2t+1,                    t>=1.       (42)
```

The quadratic polynomial in (42) has discriminant

```text
disc=-4.                                                   (43)
```

For an odd prime `p`, it has a root modulo `p` exactly when `-1` is a square:

```text
c_t=0 mod p has two roots  iff p=1 mod 4,
c_t=0 mod p has no root    iff p=3 mod 4.                 (44)
```

This also proves directly that every prime divisor of every `c_t` is
`1 mod 4`; equivalently, the represented sum of two consecutive coprime
squares is automatically admissible in (36).

Now choose any `r` distinct primes `p_j=1 mod 4` and one of the two roots in
(44) for each.  The Chinese remainder theorem gives one residue class
modulo `product p_j`, hence infinitely many positive `t`, for which

```text
product p_j divides c_t.                                  (45)
```

Therefore

```text
omega(c_t)>=r,             |X_(c_t)|>=2^(r-1).            (46)
```

The Boolean dimensions and equal-plane multiplicities on this one parabolic
branch are unbounded.

The exact companion finds the first record rungs

| `omega(c_t)` | first `t` | `c_t` | distinct prime support |
|---:|---:|---:|---|
| 1 | 1 | 5 | `5` |
| 2 | 6 | 85 | `5*17` |
| 3 | 23 | 1105 | `5*13*17` |
| 4 | 223 | 99905 | `5*13*29*53` |
| 5 | 1081 | 2339285 | `5*13*17*29*73` |
| 6 | 12131 | 294346585 | `5*13*17*41*73*89` |

Only the displayed first-record claims are finite-exact.  Unboundedness is
the CRT proof (44)--(46).

## 7. The first external K4: c=1105

At `t=6`, `c=85=5*17` gives the first collision along the spine:

```text
(13,84,85),  word U^5;
(77,36,85),  word AD.                                   (47)
```

Thus `Q=171` is the first **spine-selected** shared plane, not the first
shared plane globally.  The 2024 paper's finite Figure 1 displays the second
row in (47), not both rows; completeness of the tree, rather than that figure,
supplies both.

At `t=23`,

```text
c=1105=5*13*17                                           (48)
```

is the first spine fibre of size four.  Choose

```text
5 =(2+i)(2-i),
13=(3+2i)(3-2i),
17=(4+i)(4-i).                                           (49)
```

Modulo global conjugation, the four factor choices give

| bit representative | Euclid parameters | parent | Berggren word |
|---|---|---|---|
| `000` | `(32,9)` | `(943,576,1105)` | `DAUD` |
| `100` | `(31,12)` | `(817,744,1105)` | `UDUA` |
| `010` | `(33,4)` | `(1073,264,1105)` | `DADDD` |
| `001` | `(24,23)` | `(47,1104,1105)` | `U^22` |

The plane quotient has collapsed the ancestry depths `4,4,5,22`.  This is
the explicit branch transplant: it preserves the hypotenuse and descendant
plane while destroying the legs, Gaussian allocation, depth, and word.

Because this fibre is an affine `V_4` torsor, its six unordered pairs are the
edges of `K_4`.  XOR difference partitions them into three perfect matchings,
canonically labelled here by the three split primes:

```text
p=5:   (47,1073) and (817,943),
p=13:  (47,817)  and (943,1073),
p=17:  (47,943)  and (817,1073),                          (50)
```

where each number denotes the odd leg of its row in the table.  Equations
(48)--(50) are a literal four-vertex/six-edge/three-matching carrier of the
kind classified abstractly in [THM-2606](THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin.md)
and [THM-2753](THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration.md).

It is not intrinsically a tournament.  For any two torsor points `x,y`,
translation by `y-x` swaps them.  A translation-invariant orientation would
therefore have to orient the same pair both ways.  The complete enumeration
finds

```text
0 of the 64 tournaments on the four points are V4-translation-invariant. (51)
```

XOR supplies a symmetric edge colour, not an arrow.  Ordering the odd legs,
slopes, depths, or words would give an added transitive-tournament gauge and
would break the affine symmetry.

## 8. The other K4: generalized Fibonacci boxes

For one U-spine parent at parameter `t`, put

```text
E_t=(1,t,t+1,2t+1).                                      (52)
```

This is a four-term generalized Fibonacci window:

```text
1+t=t+1,                 t+(t+1)=2t+1.                   (53)
```

Its three `K_4` perfect matchings encode the one parent:

```text
1(2t+1)=a,                  t(t+1)=b/2,
1(t+1)+t(2t+1)=c,
(t+1)(2t+1)-1t=c.                                       (54)
```

This **internal K4 on four box entries** is different from the **external
K4 on four equal-plane parents** in Section 7.  Cardinalities four and six
do not identify their vertices, operations, or predicates; no canonical map
between the two K4s is supplied here.

The U-spine is also not the true golden Fibonacci path.  On parameters its
step is the parabolic matrix `g` from (4), which preserves `m-r=1`.  The
golden step is

```text
h=[1 1;1 0],                                               (55)
```

with spectral radius `phi`.  Consecutive Fibonacci pairs have unit gap only
at `(2,1)` and `(3,2)`, so the two paths share the first two triples and then
separate.  Modulo two,

```text
g: (0,1)<->(1,0),
h: (0,1)->(1,0)->(1,1)->(0,1).                           (56)
```

Thus the parabolic spine remains in the raw primitive parity stratum, whereas
the golden path visits the odd/odd content-two channel every third step.
Fibonacci labels alone carry no LRC or tournament predicate.

## 9. The nonlinear index identity is not a tree embedding

If depth is indexed from zero, then

```text
a_j=2j+3,
b_n=2(n+1)(n+2),
Q_n=a_(b_n).                                              (57)
```

This is an exact self-index identity.  It is not a graph, monoid, or
semigroup embedding.  The adjacent depths `0,1` map to depths `4,12`; their
spinors `(6,5)` and `(14,13)` have determinant

```text
6*13-5*14=8,                                              (58)
```

not one.  Hence (57) destroys Farey adjacency already on the first edge.

## 10. Transfer contract and cross-frontier boundary

| source | target/map | preserved | destroyed or required sidecar |
|---|---|---|---|
| fixed-cusp spinor fan | `Phi` | slope, Farey faces, Lorentz pairing | raw parity scale if normalized |
| U horocycle | `Q=a^2+2` | one quadratic height | cusp versus orbit, spinor/leg data |
| parent `P` | plane label `c` | descendant plane and its multiplicity | Gaussian representation, word, depth |
| Gaussian factor choices | quotient by conjugation | affine Boolean torsor | an origin and factor orientations |
| four-parent fibre | six XOR-coloured pairs | three perfect matchings | tournament orientation |
| one generalized-Fibonacci box | matching formulas (54) | one Pythagorean parent | external equal-plane fibre |

The discriminant `-4` in (43) is also sharply different from the cubic
discriminant in
[THM-1310](THM-1310-jacobian-counterexample-fiber-geometry-S3-resolvent-jelonek.md):

```text
disc_x(N)=-4 Q^2 L.                                      (59)
```

Here `-4` is the discriminant of a fixed quadratic polynomial and controls
Gaussian prime splitting.  In (59), the variable square class `-L` controls
the cubic monodromy; `-4Q^2` is only its square factor.  Matching the visible
constant does not classify the dimension-three Keller counterexamples or
transport their fibres to (39).

Likewise, the fixed-cusp Farey fan is a genuine instance of THM-3333's
carrier, but it supplies none of the LRC owner, phase, clock, endpoint-word,
or global-exit data.  No LRC row is removed.  The Boolean fibre grade is a
factor-choice rank, not the Keller composition-degree monoid.  No Jacobian,
Dixmier, or LRC statement follows.

## 11. Literature boundary

The published inputs are:

- Berggren/Barning/Hall: unique ternary generation of primitive Pythagorean
  triples;
- Janičková--Csókási (2023) and Kőszegyová--Csókási--Hirjak (2024): the
  U-power formula, descendant vectors, plane (25), area (28), and
  nonright/nonisosceles statements;
- Tripathi (2008/09): the fixed-hypotenuse count (36).

The fixed-cusp Farey fan, horocycle packaging, complete angle law, Gaussian
torsor formulation, unbounded CRT theorem, first `V_4` ancestry atlas,
prime-XOR matchings, two-K4 separation, and Fibonacci-path hostile are proved
here as elementary deductions.  No claim of external novelty or priority is
made.

## 12. Exact reproduction

Run

```bash
python 04-computation/berggren_parabolic_gaussian_torsor_thm3334.py
python -O 04-computation/berggren_parabolic_gaussian_torsor_thm3334.py
```

Both executions byte-match the stored transcript.  The companion checks:

- the Lorentz matrices, nilpotence, 501 U-spine powers, fixed cusp, every
  Farey face, all `Q` formulas, and the self-index hostile;
- all `1,314` primitive Euclid rows with `m<=80`, including planes, areas,
  side vectors, angle signs, and all three angle chambers;
- all `5,000` odd `c<=10,000`, with `1,073` nonempty Gaussian torsors built
  independently from factor choices and direct sum-of-two-squares search;
- every spine factorization through `t=12,500`, the first six Boolean record
  grades, five CRT split-prime sets with four lifts each, and inert-prime
  hostiles;
- the complete Berggren tree through `c=1105`, both ancestry fibres, all
  prime-XOR matchings, and all `64` tournament orientations;
- `2,000` generalized-Fibonacci boxes, the parabolic/golden path split, and
  both parity cycles.

The proof is the displayed integer/Gaussian algebra and CRT argument; the
computation audits the consequence universe and finite first-record claims.
