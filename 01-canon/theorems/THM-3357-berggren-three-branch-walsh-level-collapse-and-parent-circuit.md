---
id: THM-3357
title: "Berggren three-branch Walsh level collapse and parent circuit"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The three Berggren
  parameter children form an exact positive circuit whose minors and Lorentz
  shells recover the parent Pythagorean triple.  Four signed branch sums obey
  a Walsh quartet: the whole ternary level and its three branch-parity twists
  collapse to even powers on the middle, left, right, and constant rays.
  Branch-parity sectors, Pell/Markov and square-triangular level masses, a
  triple-valued second-moment recurrence, a prime-five sibling clock, and an
  exact factorial-moment law follow.  The sibling hypotenuses carry a
  transitive T3 and an owner-free convex-gate Horn implication.  Word order,
  nonlinear Gaussian squaring, determinant maxima, owners, and LRC phase
  height do not survive the relevant aggregate quotients; no LRC(14),
  EW-design, JC, FC(3), or arbitrary ternary-tree equivalence follows.
source: codex-2026-08-14-berggren-three-branch
audit: >
  three independent structural/hostile audits rederived the local circuit,
  Walsh/parity formulas, Pell indexing, norm tournament, convex-gate Horn
  rule, prime-five quotient, quadratic/energy recurrences, even/odd square
  split, factorial integral, and JC/FC/LRC scope boundaries; normal,
  optimized, and stored exact outputs match byte-for-byte
depends_on:
  - THM-2056-kelvin-polar-farey-defect-certificate
  - THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone
  - THM-3335-square-triangular-pell-markov-pythagorean-selector
  - THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction
  - THM-3341-u-spine-square-hypotenuse-transplant-and-triangular-plane-torsors
related:
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
  - THM-3300-factorial-gaussian-torus-bridge-and-the-archimedes-no-go
  - THM-3334-berggren-parabolic-spine-gaussian-collision-torsor
  - THM-3345-prime-xor-ancestry-path-groupoid-and-source-dependent-berggren-cost
  - THM-3346-u-spine-prime-toggle-root-atlas-and-conjugation-monodromy
  - THM-3353-split-prime-parabolic-branch-transplant-and-unary-transducer-compiler
  - THM-3356-primitive-affine-determinant-shells-parabolic-orbits-and-prime-clock-resultants
  - THM-3358-admissible-composite-parabolic-compiler-and-hensel-normal-offset-atlas
script: 04-computation/berggren_three_branch_walsh_level_thm3357.py
output: 05-knowledge/results/berggren_three_branch_walsh_level_thm3357.out
script_sha256: 84d65e03afae17f5c1c24b84b4f646805b5931cabf2dd8686f2111813af4ca68
output_sha256: 4868cb22c216c2a35327d79977de0bc93d0a5085cc368a75af85959bb9d188db
hash_basis: working-tree bytes (LF)
---

# THM-3357 -- all three Berggren branches at once

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This is a repository synthesis and proof interface; no literature-priority
claim is made.  The point is to replace three unrelated unary spines by one
typed ternary calculation without pretending that the branch labels form a
cyclic group.

## 1. Inheritance, hostile controls, and convention

The closest proved mechanisms are [THM-3334](THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md),
which identifies one parabolic Berggren spine, and
[THM-3339](THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction.md),
which writes all three parameter branches but selects only three special
Fibonacci rays.  The canonical hostile is [THM-2596](THM-2596-modular-free-factor-farey-gram-owner-cocycle.md):
the three branch maps have infinite order and are not one `C3` action.  The
corrected near miss is that three selected rays do not classify the full
ternary tree.  The least-used sidecars taken up here are the signed sibling
minors, branch-count parity, and the Lorentz second moment.  THM-3345 remains
the warning that branch counts do not retain source-dependent ancestry words.

Use positive Euclid parameters

```text
u=(m,n)^T,                 0<m<n,                        (1)
Psi(u)=(a,b,c)=(n^2-m^2,2mn,n^2+m^2).                  (2)
```

When `gcd(m,n)=1` and `m+n` is odd, (2) is a primitive Pythagorean triple.
In the parameter convention of THM-3339, the left, middle, and right children
are

```text
L=[ 0 1;-1 2],       M=[0 1;1 2],       R=[1 0;2 1],   (3)
Lu=(n,2n-m),          Mu=(n,2n+m),       Ru=(m,2m+n).   (4)
```

Their slopes lie respectively in `(1,2)`, `(2,3)`, and `(3,infinity)`, so the
labels in (3) are intrinsic in this positive chamber.  The root is

```text
u_0=(1,2),                 Psi(u_0)=(3,4,5).             (5)
```

Matrix products act on columns.  A root word records successive branch
choices; its matrix product therefore reads rightmost action first.

## 2. The sibling circuit reconstructs its parent

Put

```text
x=Lu,                      y=Mu,                      z=Ru. (6)
```

Direct expansion gives the complete oriented minor table

```text
det(x,y)=b,                det(x,z)=c,                det(y,z)=a, (7)
```

and hence the unique positive circuit

```text
a x+b z=c y.                                             (8)
```

Thus the three labelled sibling spinors reconstruct the parent's entire
Pythagorean triple from their pairwise determinants.  This is a weighted
oriented two-simplex, not merely three vertices with a cosmetic tournament.

Let `X=Psi(x)`, `Y=Psi(y)`, and `Z=Psi(z)`.  THM-3333's Lorentz polarization
turns (7) into

```text
<X,Y>_L=2b^2,             <X,Z>_L=2c^2,             <Y,Z>_L=2a^2. (9)
```

Consequently

```text
det[X Y Z]=-4abc,                                      (10)
det Gram_L(X,Y,Z)=16a^2b^2c^2,                        (11)
det(x,y)^2+det(x,z)^2+det(y,z)^2=2c^2.                (12)
```

Equation (12) is the local determinant energy.  The signs in (7), or the
slope labels, are load-bearing: the absolute Lorentz shells alone retain
only `a^2,b^2,c^2`.

There is an exact outer linear intertwiner.  With

```text
H=[-1 1;1 1],                  H^2=2I,                  (13)
```

one has

```text
HL=RH,                         HM=MH,                  HR=LH. (14)
```

For a primitive opposite-parity `u`, the coordinates of `Hu=(n-m,n+m)` are
coprime and odd, and

```text
Psi(Hu)=2(b,a,c).                                      (15)
```

Thus `H` swaps the outer branches and the two legs only after a raw-content
two normalization.  Since `det(H)=-2`, it is not an integral tree
automorphism, and (14) does not manufacture a three-cycle.

## 3. The sibling norm tournament and a strict determinant-gate Horn rule

Write

```text
N_L=||Lu||^2,              N_M=||Mu||^2,              N_R=||Ru||^2. (16)
```

The exact differences are

```text
N_M-N_L=4b,                N_M-N_R=4a,                N_R-N_L=4(b-a), (17)
N_L+N_R=6c.                                               (18)
```

Hence `N_M` is always largest.  The three values are distinct: `a=b>0`
would make `c^2=2a^2`.  Orienting an arc toward the larger child therefore
gives an intrinsic transitive tournament on the three siblings.  Its only
variable arc is `L--R`, oriented by the sign of `b-a`.

This sign has a two-state branch automaton.  If

```text
sigma(u)=sign(b-a),                                      (19)
```

then

```text
sigma(Lu)=+1,              sigma(Mu)=-sigma(u),         sigma(Ru)=-1. (20)
```

Thus `L` and `R` reset the state and `M` flips it.  Moreover

```text
|N_R-N_L|=4|b-a|>=4.                                    (21)
```

By THM-3341, equality in (21) occurs exactly on the positive middle-only ray
of consecutive-leg triples

```text
(3,4,5),(21,20,29),(119,120,169),... .                  (22)
```

The outer arc alternates on that ray.  The tournament records this bit but
forgets the magnitudes and the circuit (8).

There is also a useful LRC-facing implication.  Fix one ordered rank-two
deck `C={c_i}` in the sense of THM-2056 and define

```text
Delta(v)=max_i |det(v,c_i)|,
G_kappa(v)=||v||^2-kappa Delta(v),              kappa>=0. (23)
```

The support function `Delta` is positively homogeneous and subadditive.
Applying it to (8) gives

```text
c Delta(y)<=a Delta(x)+b Delta(z).                       (24)
```

The norm identities above give the strictly positive correction

```text
K(u)=cN_M-aN_L-bN_R
    =2m(n-m)(3n^2+4mn-m^2)>0.                           (25)
```

Combining (24)--(25) proves, for every fixed deck and every `kappa>=0`,

```text
c G_kappa(Mu)
 >=a G_kappa(Lu)+b G_kappa(Ru)+K(u).                    (26)
```

In particular, if both outer siblings have nonnegative score, the middle
sibling has strictly positive score.  At `kappa=91`, (26) is an owner-free
Horn rule for THM-2056's LRC(14) determinant gate:

```text
outer L and R gate-pass  ==>  middle M gate-passes strictly. (27)
```

This is a real search-pruning implication, but only for the sufficient
determinant gate in one fixed basis/deck.  It does not recover the active
owner, phase height, or runner labels, and it supplies no outer seed or full
LRC-row certificate.  It does certify the middle direction for THM-2056's
sufficient gate once the two outer gate seeds are known.  The standard
Euclidean norm in (23) is basis-dependent exactly as in THM-2055/2056; only
the algebraic circuit and Walsh identities below are freely conjugation
covariant.

### 3A. The only local shared-prime edge is a prime-five clock

For every primitive opposite-parity parent in (1), the three child
hypotenuses have a rigid gcd graph:

```text
gcd(N_L,N_M)=5  iff  5|m,              and is 1 otherwise,
gcd(N_M,N_R)=5  iff  n=m mod 5,        and is 1 otherwise,
gcd(N_L,N_R)=1.                                         (27a)
```

In the positive cases the minimum of the two `5`-adic valuations is exactly
one, so the gcd is `5`, even though one child may contain a higher power.  For the first row, use
`N_M-N_L=8mn` and primitivity.  A common odd prime must divide `m` or `n`;
reducing `N_L` modulo those coordinates leaves only `p=5` with `5|m`.
Writing `m=5k` and reducing both norms modulo `25` excludes a second common
factor of five.  For the second row, use
`N_M-N_R=4(n-m)(n+m)`: the `n+m` factor would force an odd prime to divide
`2m^2`, while the `n-m` factor leaves only `5`; the mod-25 argument is the
same.  Finally a common odd divisor of the outer norms divides
`N_L+N_R=6c`.  Reduction modulo `c=m^2+n^2`, with `p=3` checked separately,
gives a contradiction.

This is a symmetric shared-support edge, not a tournament arc.  Each parent
has either no collision edge or exactly one of `L--M` and `M--R`.

There is a three-state **unlabelled level-count** clock.  On the projective slope
`t=n/m in P^1(F_5)`, the branches act by

```text
L:t->2-1/t,                  M:t->2+1/t,               R:t->t+2. (27b)
```

Partition the six residues into

```text
Gamma={infinity,1},          Xi={0,4},             Upsilon={2,3} (27c)
```

The two residues in `Gamma` are exactly the collision conditions in (27a).
After forgetting the outgoing branch label, these classes have
transition-count matrix and root state

```text
[Gamma_(d+1)]   [1 2 0][Gamma_d]                    [Gamma_0]   [0]
[Xi_(d+1)]    = [0 0 3][Xi_d],                       [Xi_0]    = [0]
[Upsilon_(d+1)] [2 1 0][Upsilon_d]                    [Upsilon_0] [1]. (27d)
```

The number `c_d=Gamma_d` of depth-`d` parents with a shared-five sibling edge is

```text
0,0,6,6,24,96,222,726,2256,...,                         (27e)
c_(d+3)=c_(d+2)+3c_(d+1)+9c_d,
sum_(d>=0)c_d t^d=6t^2/((1-3t)(1+2t+3t^2)).             (27f)
```

Therefore `c_d/3^d -> 1/3`; among all three sibling-pair slots, the collision
edge density tends to `1/9`.  This prime-five clock loses height, exact norm,
and ancestry word.  The classes in (27c) are equitable only after summing over
the three outgoing branch labels; the full six-state labelled automaton
suffices for wordwise prediction.  Also the sibling-collision class is
`Gamma={infinity,1}`, whereas `5` divides the **parent** hypotenuse exactly in
`Upsilon={2,3}`.  This is not the fixed-hypotenuse prime-XOR torsor of
THM-3345/3346, THM-3356's affine-shell resultant fingerprint, or THM-3358's
exact-grade Hensel-normal compiler: those retain equal-norm allocation, a
fixed affine shell, or a selected normal lift, whereas (27a) is a local
sibling shared-support event.

## 4. One weighted transfer compiles every ternary level

For commuting weights `x,y,z`, put

```text
F(x,y,z)=xL+yM+zR
        =[ z, x+y; -x+y+2z, 2x+2y+z].                  (28)
```

If `S_d(x,y,z;u)` is the sum over all depth-`d` descendants, weighting a word
by `x^(#L)y^(#M)z^(#R)`, distributivity gives

```text
S_d(x,y,z;u)=F(x,y,z)^d u.                              (29)
```

This is an identity of polynomials, so coefficient extraction gives the
aggregate spinor for any prescribed triple of branch multiplicities.  A
fixed numeric weight specialization is computable in `O(log d)` matrix
multiplications rather than by enumerating `3^d` nodes.  Expanding all formal
coefficients has `Theta(d^2)` support and is a different complexity claim.

Writing

```text
Sigma=x+y+z,                  delta=x^2-y^2+z^2,         (30)
```

one has

```text
tr F=2Sigma,                 det F=delta,                (31)
S_(d+2)=2Sigma S_(d+1)-delta S_d,                       (32)
sum_(d>=0) S_d t^d=(I-tF)^(-1)u.                        (33)
```

The quotient in (29) destroys ancestry order.  Already at the root, the two
words with one `L` and one `M` give

```text
MLu_0=(3,8),                    LMu_0=(5,8).             (34)
```

Their coefficient in (29) retains only the sum.  THM-3345's source-dependent
path word therefore cannot be reconstructed from branch multiplicities.
Here (34) uses matrix-product notation; a chronological branch string is
read in reverse matrix order under the convention stated after (5).

## 5. Four Walsh characters collapse the full ternary tree to four rays

The three matrices satisfy the exact quartet

```text
L+M+R =M^2,              L+M-R =L^2,
-L+M+R=R^2,              L-M+R =I.                      (35)
```

Equivalently,

```text
M^2=2M+I,       L^2=2L-I,       R^2=2R-I,       L+R=M+I. (36)
```

Specializing (29) at the four sign characters gives, for every root `u` and
every `d>=0`,

```text
sum_(|w|=d) w u                         =M^(2d)u,
sum_(|w|=d) (-1)^(#R(w)) w u            =L^(2d)u,
sum_(|w|=d) (-1)^(#L(w)) w u            =R^(2d)u,
sum_(|w|=d) (-1)^(#M(w)) w u            =u.             (37)
```

Thus all `3^d` nodes collapse unsigned to one point twice as deep on the
middle ray.  The other three character sums land twice as deep on the left
ray, twice as deep on the right ray, and back at the root.  This is the
simultaneous three-spine law.

At the standard root these three rays are already closed forms:

```text
M^(2d)u_0=(P_(2d+1),P_(2d+2)),
L^(2d)u_0=(2d+1,2d+2),             R^(2d)u_0=(1,4d+2).  (37a)
```

Here the Pell numbers `P_r` are defined in (44).  Thus every parity-sector
sum below is an explicit quarter-combination of one Pell state, two linear
states, and the root.

There is an exact branch-parity refinement.  Let

```text
p=(p_L,p_M,p_R) in F_2^3,       p_L+p_M+p_R=d mod 2,    (38)
```

and let `V_(d,p)(u)` be the sum of those descendants whose three branch
counts have parity `p`.  Fourier inversion on the four-point affine plane
(38) gives

```text
4V_(d,p)(u)
 =M^(2d)u+(-1)^p_R L^(2d)u
          +(-1)^p_L R^(2d)u+(-1)^p_M u.                (39)
```

The number of nodes in the sector is

```text
#V_(d,p)={3^d+(-1)^p_L+(-1)^p_M+(-1)^p_R}/4.            (40)
```

Parity triples violating (38) are empty.  Over `Z[1/2]` (or any coefficient
ring in which `2` is invertible), the four characters in (37) are the full
Walsh basis on (38); they retain branch-count parity, not letter order.  Over
`Z`, (39) remains a cleared identity with divisibility by four.  In
characteristic two the characters coalesce, so no Fourier-recovery claim is
made there.

Equations (35)--(37) and the cleared identity (39) are not peculiar to
positive integers.  They base-change to endomorphisms over any commutative
coefficient ring satisfying the quartet (35), with the inversion qualification
just stated, and they survive every simultaneous conjugation.  Formula (40)
is instead an integer word-count identity.

```text
(L,M,R)->(CLC^(-1),CMC^(-1),CRC^(-1)).                  (41)
```

This is the exact port to another ternary matrix presentation: test its four
signed branch sums.  It is not a theorem about every ternary tree.  For the
unimodular control

```text
R_bad=[1 1;2 3],                                        (42)
```

one has

```text
L+M+R_bad=[1 3;2 7]!=M^2,          L-M+R_bad!=I.        (43)
```

Thus arity three, unimodularity, and positivity do not force Walsh collapse.
The control in (42) is a branch-matrix triple, not a claim about a genuine
disjoint-chamber enumeration tree.  Conversely, quartet closure compiles a
multiset of branch words; it does not by itself prove injectivity, positivity,
disjoint chambers, or treehood.

## 6. The full-level mass is the square-triangular Pell/Markov clock

Let the Pell numbers be

```text
P_0=0,          P_1=1,          P_(r+2)=2P_(r+1)+P_r.   (44)
```

At the root, the unsigned case of (37) is

```text
sum_(|w|=d) w u_0=M^(2d)u_0=(P_(2d+1),P_(2d+2)).       (45)
```

Define

```text
alpha_d=P_(2d+1),       s_d=P_(2d+2),       beta_d=P_(2d+3). (46)
```

Pell Cassini and THM-3335 give

```text
alpha_d beta_d-s_d^2=1,
alpha_d^2+beta_d^2+4=6alpha_d beta_d,                   (47)
n_d=(alpha_d+beta_d-2)/4,        q_d=s_d/2,
T_(n_d)=q_d^2.                                          (48)
```

Thus the first-coordinate masses of consecutive full levels are adjacent
fixed-two Markov coordinates, while the intervening second-coordinate mass
is their Cassini carry.  The first rows are

| `d` | full-level spinor mass | next first mass `beta_d` | square-triangular index `n_d` |
|---:|---:|---:|---:|
| 0 | `(1,2)` | 5 | 1 |
| 1 | `(5,12)` | 29 | 8 |
| 2 | `(29,70)` | 169 | 49 |
| 3 | `(169,408)` | 985 | 288 |
| 4 | `(985,2378)` | 5741 | 1681 |
| 5 | `(5741,13860)` | 33461 | 9800 |

Every coordinate in (45)--(46) obeys

```text
X_(d+2)=6X_(d+1)-X_d,                                  (49)
```

and `n_(d+2)=6n_(d+1)-n_d+2`.  Equations (45)--(49) are an
`O(log d)` ring-arithmetic/matrix-multiplication compiler for a global
statistic of a `3^d`-node level.  This is not a bit-complexity claim—the
outputs have `Theta(d)` bits—and it is not a nodewise bijection.

At depth two, the nine-node second-coordinate mass is the exceptional
cannonball value:

```text
70^2=sum_(r=1)^24 r^2
    =sum_(r=1)^12 (2r)^2+sum_(r=1)^12 (2r-1)^2
    =2600+2300.                                         (50)
```

The last line splits the **classical range** `1,...,24`; it is not the split
of the even and odd Euclid parameters at ternary depth two.

THM-3335's global selector/cannonball classification implies that `d=2` is
the unique positive full-level depth whose `s_d^2` is a positive square
pyramidal number.  The same level compiler supplies two different scalar
addresses already separated in THM-3335:

```text
square-arc tournament order:     n_d+1=2,9,50,289,...,
skew-EW arithmetic candidate:    s_d^2+2=6,146,4902,... . (51)
```

The second list satisfies THM-476's necessary square gate; it supplies no
EW design.  Apart from the trivial order two, THM-3335 proves that the two
order predicates in (51) cannot hold for the same order.

## 7. Nodewise Gaussian squaring needs a second transfer

The quadratic map `Psi` does not commute with the linear collapse.  Its three
induced Berggren matrices are

```text
T_L=[ 1 -2 2; 2 -1 2; 2 -2 3],
T_M=[ 1  2 2; 2  1 2; 2  2 3],
T_R=[-1  2 2;-2  1 2;-2  2 3].                         (52)
```

For commuting weights,

```text
det(xT_L+yT_M+zT_R)
 =(x-y-z)(x+y-z)(x+y+z).                                (53)
```

The unsigned triple transfer is

```text
K=T_L+T_M+T_R=[1 2 6;2 1 6;2 2 9].                    (54)
```

Put

```text
(A_d,B_d,C_d)=sum_(|w|=d) Psi(wu_0)=K^d(3,4,5).        (55)
```

Then

```text
A_d-B_d=(-1)^(d+1),                                    (56)
sum_(d>=0)(A_d+B_d)t^d=(7-3t)/(1-12t+3t^2),
sum_(d>=0)C_d t^d=(5-t)/(1-12t+3t^2).                  (57)
```

Each of `A_d,B_d,C_d` obeys

```text
X_(d+3)=11X_(d+2)+9X_(d+1)-3X_d,                       (58)
```

while `A_d+B_d` and `C_d` obey

```text
Y_(d+2)=12Y_(d+1)-3Y_d.                                (59)
```

The first rows are

```text
d=0: (3,4,5),             d=1: (41,40,59),
d=2: (475,476,693),       d=3: (5585,5584,8139),
d=4: (65587,65588,95589).                               (60)
```

The smallest hostile is already the first level:

```text
sum Psi(child)=(41,40,59),
Psi(sum child)=Psi(5,12)=(119,120,169).                 (61)
```

Thus the Pell collapse (45) cannot be pushed through Gaussian squaring.

The nodewise triple transfer also gives an internal even/odd square law.
Every primitive node has one even and one odd parameter.  Let `Eeven_d` be
the sum of the squares of all even parameters at level `d`, and `Eodd_d` the
corresponding odd-parameter sum.  Give a node sign `rho=+1` when its second
coordinate is even and `rho=-1` when its first coordinate is even.  Modulo
two, `L` and `M` swap the coordinates while `R` fixes them.  Hence the signed
quadratic transfer is

```text
-T_L-T_M+T_R=[-3 2 -2;-6 1 -2;-6 2 -3],                (61a)
```

whose characteristic polynomial is `(lambda+1)^2(lambda+3)`.  Its first
coordinate from `(3,4,5)` is `(-1)^d(4*3^d-1)`.  Since that coordinate is
exactly the even-square total minus the odd-square total,

```text
Eeven_d+Eodd_d=C_d,
Eeven_d-Eodd_d=(-1)^d(4*3^d-1),                         (61b)
Eeven_d={C_d+(-1)^d(4*3^d-1)}/2,
Eodd_d ={C_d-(-1)^d(4*3^d-1)}/2.                       (61c)
```

For example,

```text
(Eeven_0,Eodd_0)=(4,1),       (Eeven_1,Eodd_1)=(24,35),
(Eeven_2,Eodd_2)=(364,329).                              (61d)
```

The `364+329=693` split is the actual depth-two Euclid-parameter square
partition, typed separately from `2600+2300=70^2` in (50).

There is nevertheless an exact quadratic sidecar.  Let

```text
E_d=sum_({v,w} subset level d) det(v,w)^2.              (62)
```

Cauchy--Binet applied to `sum_v vv^T` gives

```text
E_d=(C_d^2-A_d^2-B_d^2)/4.                             (63)
```

Consequently

```text
E_d: 0,50,7012,967510,133454152,18407967994,...,        (64)
E_(d+4)=142E_(d+3)-564E_(d+2)+450E_(d+1)-27E_d,        (65)
sum_(d>=0)E_d t^d
 =(50t-88t^2+6t^3)/(1-142t+564t^2-450t^3+27t^4).       (66)
```

For mechanism rather than prefix fitting, let `lambda_+,lambda_-` be the
roots of `lambda^2-12lambda+3`.  Equations (56)--(57) express (63) using
only products of the modes `lambda_+,lambda_-,-1`.  Its possible modes are

```text
lambda_+^2,       lambda_+lambda_-=3,       lambda_-^2,       1,
```

whose characteristic polynomial is
`(T-1)(T-3)(T^2-138T+9)`.  This is exactly recurrence (65) and denominator
(66).

This is the total pairwise Lorentz-shell energy of one full ternary level.
The scalar does not itself exhibit an extremal shell or owner; recovering
either requires the determinant distribution.  Section 8 gives the exact
non-recovery boundary on same-level subpackets, not on complete levels.

### 7A. The triple-transfer determinant has a closed factorial moment law

Let

```text
P(x,y,z)=det(xT_L+yT_M+zT_R)
        =(x-y-z)(x+y-z)(x+y+z),                         (66a)
```

and let the three-variable factorial functional be

```text
calL(x^i y^j z^k)=i!j!k!.                              (66b)
```

Then, for every `r>=0`,

```text
calL(P^r)=0                              if r is odd,
calL(P^r)=(3r+2)!/[2(r+1)^2]            if r is even.   (66c)
```

The exponential-integral representation of `calL`, followed by
`t=x+y+z`, separates the radial integral

```text
integral_0^infinity exp(-t)t^(3r+2)dt=(3r+2)!           (66d)
```

from the standard two-simplex.  On that simplex put

```text
u=2x-1,                       v=1-2z.                   (66e)
```

The simplex becomes `-1<=u<=v<=1`, the Jacobian is `1/4`, and `P=uv`.
Since `(uv)^r` is symmetric, its integral over the ordered triangle is one
half its integral over the square.  The simplex factor is therefore

```text
(1/8)(integral_(-1)^1 u^r du)^2,
```

which is zero for odd `r` and `1/[2(r+1)^2]` for even `r`.  This proves
(66c).  The first values are

```text
1,0,2240,0,1743565824,... .                             (66f)
```

This is an FC(3)-facing detection benchmark: every odd factorial moment
vanishes, but `calL(P^2)=2240>0`.  It is not an FC(3) counterexample and
proves no general factorial-conjecture case.  Its value is an exact infinite
moment formula produced by all three branch matrices, not a finite
coincidence.

## 8. Exact loss ledger and hostile maximum

Inside the actual depth-three Berggren level, the two three-node subpackets

```text
P={LLM=(4,11),MLR=(5,18),MRL=(9,16)},
Q={LMM=(8,19),LMR=(3,14),LRL=(7,12)}                    (67)
```

have the same full first and quadratic moments:

```text
sum P=sum Q=(18,45),
sum_(v in P) vv^T=sum_(v in Q) vv^T=[122 278;278 701],
sum Psi(P)=sum Psi(Q)=(579,556,823),       E(P)=E(Q)=8238. (68)
```

But their sorted absolute determinant patterns are

```text
(17,35,82)                    and                  (37,55,62). (69)
```

Thus even the full first and second aggregate moments do not recover the
maximum pairwise shell on same-level subpackets.  This still does not assert
that two **complete** Berggren levels have the same aggregate state.

There is also a fixed-candidate owner hostile.  Put `d=(2,3)`.  Both packets

```text
P'={(1,4),d,(2,9)},              Q'={(1,6),d,(2,7)}      (69a)
```

have first moment `(5,16)` and energy `170`.  Against the common candidate
`d`, their deck determinant lists are `(5,12)` and `(9,8)`: the maxima are
`12` and `9`, with different owners.  Finally, THM-3335's labelled clock
benchmark keeps the same spinor, Gaussian point, determinant maximum, signed
owner, and Kelvin defect across tail placements `j=7,...,12`, while the exact
phase maxima are `1/7,...,1/12`.  Every spinor-only Walsh datum is identical
there.  The labelled column/clock placement is therefore indispensable for
LRC phase height.

| source | target/map | exactly preserved | destroyed information | needed sidecar / decisive hostile |
|---|---|---|---|---|
| full ternary level | `F^d u` | weighted first moment, branch multiplicities | word order and source path | noncommutative word or THM-3345 groupoid; (34) |
| full level | four Walsh rays | branch-count parity sums | ordering within a parity class | ordered language; (34) |
| sibling triple | norm `T3` | hypotenuse order and sign state | circuit weights and LRC phase height | three gaps or minors; (17) |
| same-level spinor subpacket | first and quadratic moments | full quadratic moment and total shell energy | maximum shell and owner | determinant distribution; (67)--(69a) |
| spinor sum | `Psi(sum)` | Gaussian square of the aggregate | sum of nodewise triples | separate quadratic transfer; (61) |
| root level mass | Pell/Markov selector | recurrences, square-triangular index | nodewise square-leg claim | depth/type label; (45)--(51) |
| simultaneous conjugate branches | (41) | quartet and Walsh identities | standard Euclidean/LRC metric unless transported | conjugated metric/deck |
| arbitrary unimodular ternary system | same proposed map | nothing automatic | quartet itself | signed-sum test (43) |

## 9. Frontier consequences and stopping boundaries

### LRC(14)

Equation (27) is the positive gain: for a fixed rank-two determinant deck,
one need not treat the middle sibling as an independent gate failure once
both outer siblings pass.  Equations (67)--(69a) show that the displayed
moment tuple cannot be a universal recovery rule for arbitrary primitive
packets, even for subpackets selected within one Berggren level; they do not
assert a collision between complete levels.  THM-3335's labelled-plane
hostile separately shows that identical spinor/Kelvin scalar state can have
different phase height.  No open LRC row is closed here.

### Tournament structure

The sibling hypotenuse comparison is an intrinsic binary relation with no
ties, so the transitive `T3` of Section 3 is lawful.  Its middle vertex is
always the champion and its outer arc is the two-state automaton (20).  A
whole Berggren level has no corresponding intrinsic pairwise orientation
from (29) or (37); branch parity is a grading, not a tournament.

### Efficient sequences, square sums, and skew-EW arithmetic

Equations (32), (37), (49), (58), and (65) replace exponential level
enumeration by fixed-order recurrences or logarithmic matrix powering.
Equations (48)--(51) explain how square-triangular, even/odd square-sum, and
skew-EW arithmetic sequences arise as **level masses**.  They do not turn
every node at that level into a selector node or a design.

### Planar JC and FC(3)

The unsigned collapse uses the even decimation `M^(2d)` of the norm-`-1`
silver recurrence, so it retains the positive Pell magnitude while erasing
the alternating norm sign emphasized in the planar-JC metallic route.  The
Walsh sidecar records branch-count parity, not that lost chronological
parity/owner circuit.  Concretely, the linear map from parent triple
coordinates `(a,b,c)` to `(N_L,N_M,N_R)` has determinant `48`, so it is
three-dimensional rather than planar.  The natural planar projection

```text
(m,n)->(N_L,N_R)
```

has Jacobian `48(m^2-2mn-n^2)`, which is nonconstant and has a complex zero
curve.  It is not a Keller map.  Section 7A is the positive FC(3) contact: the
triple-transfer determinant has an exact factorial moment formula and is
detected at exponent two.  It does not settle FC(3).  For the direct
Gaussian/Pythagorean port, THM-3300's factorial invariant quotient keeps
norm-type moments but not the represented angle and determinant shells in
(7)--(9); that port still needs a signed/weight-two sidecar.  No planar-JC
circuit theorem follows.

### A procedure for another ternary tree

For another triple of branch endomorphisms `(B_0,B_1,B_2)`, the cheapest
useful audit is:

1. fix one common lattice, action convention, and positive chamber;
2. compute the four signed sums corresponding to (35);
3. if the quartet closes, use (29) and Fourier inversion before enumerating;
4. compute the three sibling minors and ask whether they form a positive
   circuit such as (8);
5. transport the quadratic form separately before making norm or LRC claims;
6. test one pair of anagram words and one same-moment/different-maximum packet.

Passing steps 2--3 gives a lawful branch-word-multiset Walsh compiler even
outside the present positive Berggren realization; treehood and injectivity
remain separate tests.  Failing the quartet, as in (43), is an exact stopping
reason rather than evidence against every other statistic on that system.

The main open pulls are to classify integral positive-chamber ternary
quartets satisfying (35), to find the least order-aware sidecar refining
(39), and to replay the Horn rule (27) on actual THM-2056 residual decks.

## 10. Reproduction

Run

```bash
python3 04-computation/berggren_three_branch_walsh_level_thm3357.py
python3 -O 04-computation/berggren_three_branch_walsh_level_thm3357.py
```

The exact companion checks:

- `1714` primitive parents in an explicit parameter box for (7)--(21);
- the convex-gate inequality for `2170` parent/`kappa` cases including
  `kappa=91`;
- all `125` integer weight triples in `[-2,2]^3` for the structural formulas
  (28), (30)--(32), with direct enumeration versus (29) at six signed and
  nonsigned specializations;
- direct tree enumeration against all four Walsh characters and parity
  sectors through depth eight;
- the Pell/Markov compiler through depth ten and the cannonball split;
- the triple transfer, determinant factorization, recurrences, and
  Cauchy--Binet energy by an independent path, including the internal
  even/odd parameter-square split;
- the prime-five sibling gcd law and collision recurrence through depth ten;
- the factorial determinant formula through exponent eight by exact
  polynomial expansion, independently of its simplex proof;
- the nonlinear, word-order, same-moment/different-maximum, conjugation, and
  arbitrary-ternary hostile controls.

The bounded checks verify the implementation and constants.  The global
quantifiers in Sections 2--7 follow from the displayed matrix, determinant,
Fourier, Cayley--Hamilton, and Cauchy--Binet identities.
