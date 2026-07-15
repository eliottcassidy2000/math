---
id: THM-830
title: Gap-tournament suspension and the weighted deletion deck classify the B3 faces, half cube, and radius-one metagraph
status: PROVED (all-size coordinate, face, deletion, seam, half-cube, line-defect, and radius-one laws) + FINITE-EXACT (full tilings through n=7; root invariants through n=9; formula replay through n=14)
source: codex-2026-07-14-S14
depends_on: [THM-002, THM-442, THM-513, THM-549, THM-550, THM-553, THM-559, THM-796, THM-801]
related: [THM-781, THM-809, THM-811, THM-814, THM-818, THM-825, THM-828, THM-832, THM-838, THM-839, THM-841, HYP-2685, HYP-3234, HYP-6880]
verification:
  - 04-computation/b3_gap_tournament_deletion_deck_codex_S14.py
  - 05-knowledge/results/b3_gap_tournament_deletion_deck_codex_S14.out
---

# THM-830 - the gap tournament and weighted deletion deck

The fixed-Hamiltonian-path tiling at size `n` is exactly a labelled tournament
on the `n-1` gaps of that path.  This is the missing tournament meaning of the
middle `B` face.  It also exposes the one bit that ordinary internal deletion
must retain: the shortcut seam made adjacent when a path vertex is removed.

The same coordinates give three further exact descriptions.

1. The full `A+B+C-D-E-F+G` recursion is the complete-homogeneous
   polynomial recurrence of a three-axis address simplex.
2. A full tiling is an arrow between two half-tilings with the same fixed
   boundary.  Blue is the identity arrow, black is a nonzero additive defect,
   and the defect alphabet at size `n` is literally the blue cube at size
   `n-1`.
3. On the entire one-flip shell, the tile, ordinary node, converse-merged node,
   complement line, reflection orbit, colour, `C3`, and `H` are classified by
   one deletion profile.

These are tournament/metagraph theorems.  They do not preserve the LRC metric
predicate, so they reorganize the `n=14` information frontier without proving
the fourteen-runner case.

## 1. A fixed-path tiling is a tournament on path gaps

Use THM-796's cube

```text
X_n=F_2^{S_n},
S_n={(a,b):1<=b<a<=n, a-b>=2},
m_n=|S_n|=binom(n-1,2).
```

Put `N=n-1`.  For `t in X_n`, define a labelled tournament `G_t` on the
ordered gap vertices `[N]` by

```text
(G_t)_(i,j)=t_(i+1,j),             1<=j<i<=N.              (1.1)
```

Every pair `(i,j)` occurs once, so

```text
Gamma_n:X_n -> LabTour_(n-1),      t |-> G_t               (1.2)
```

is a literal cube bijection.  Its inverse `Sigma_n(G)` fixes the marked path
`n->...->1` and uses the gap arc `(a-1,b)` to orient every nonpath arc `(a,b)`.
This is a path suspension, not an induced subtournament construction.

Let `kappa` be all-tile antipode, `sigma` staircase reflection, and
`rho(i)=N+1-i`.  Directly from (1.1),

```text
Gamma_n(kappa t)=G_t^op,
Gamma_n(sigma t)=rho(G_t^op).                              (1.3)
```

Consequently the complement-line sort has the exact tournament model

```text
L_n=X_n/<kappa>
   ~= LabTour_(n-1)/{G~G^op}.                              (1.4)
```

An upper line is blue exactly when

```text
G=rho(G^op).                                               (1.5)
```

Reflection sends the line `e_G={G,G^op}` to `e_(rho G)`.  A fixed line would
require `rho G=G` or `G^op`.  The first would make `rho` a nontrivial
tournament automorphism of order two, which is impossible; hence the fixed
lines are exactly (1.5), recovering blue versus black without using endpoint
nodes.

This separates three operations that were historically blurred:

```text
all-tile antipode     path-relative reversal T -> T^(dagger_p)
grid reflection      true converse plus reversal of the marked path
merged-node map      forget p, quotient by isomorphism, then merge converse.
```

If `Phi_n(G)=pi_n(Sigma_n(G))`, then the complete three-sort pushforward is

```text
K_n^c(u,v)=#{G:Phi_n(G)=u, Phi_n(G^op)=v, colour(G)=c}.     (1.6)
```

It is symmetric.  If `Fib_n(u)={G:Phi_n(G)=u}`, summing over endpoint `v` and
both colours gives `|Fib_n(u)|`; at fixed colour the row sum is the corresponding coloured
half-edge count.  For `u!=v`, `K_n^c(u,v)` is the literal line multiplicity,
while `K_n^c(u,u)/2` is the loop multiplicity.
Thus nodes are fibres of `Phi`, tilings are labelled gap tournaments, and
edges are their converse orbits.  No relation among the three sorts is left
implicit.

There is also a second node map

```text
gamma_n:L_n -> M_(n-1)^gap                             (1.7)
```

which sends `e_G` to the ordinary converse-merged class of `G`.  If `G` is a
self-converse class, the fibre contains

```text
(n-1)!/(2|Aut G|)
```

upper lines; a non-self-converse merged pair contributes

```text
(n-1)!/|Aut G|.
```

This is the all-orderings counterpart of THM-781's `HP(T)/Aut(T)` fibre.

## 2. Every A/B/C face is a tournament operation

Apply `Gamma_(n-1)` after each THM-801 face.  For `p>q` in `[N-1]`,

```text
A(G)_(p,q)=G_(p,q),
B(G)_(p,q)=G_(p+1,q),
C(G)_(p,q)=G_(p+1,q+1).                                  (2.1)
```

Thus `A` and `C` delete the maximum and minimum gap vertices.  The middle
face is a tournament on cuts: cut `p` has lower incarnation `p`, upper
incarnation `p+1`, and the omitted adjacent arc `G_(p+1,p)` is its self-loop
or seam bit.  The statement in THM-801 that `B` is not deletion of one fixed
original vertex remains true; (2.1) says what it is instead.

The whole commuting face monoid has a single formula.  For nonnegative
`alpha,beta,gamma`, put `r=alpha+beta+gamma`.  Then

```text
D_(alpha,beta,gamma)(G)_(p,q)
  =G_(p+beta+gamma,q+gamma),       p>q in [N-r].            (2.2)
```

It retains precisely the arcs with

```text
j>gamma,       i<=N-alpha,       i-j>beta.                 (2.3)
```

Reversal exchanges `alpha,gamma` and fixes `beta`.  Equations (2.2)-(2.3)
give exact tournament semantics to the three pair overlaps and the triple
overlap, not only to their cardinalities.

Iterating the gap face peels successive span words:

```text
B^r(G)_(p,q)=G_(p+r,q),
LabTour_N ~= LabTour_(N-r) x product_(h=1)^r F_2^(N-h).    (2.4)
```

The discarded factor at step `h` is the arc-span word
`(G_(j+h,j))_j`.  At the marked-path level, `B^r` uses length-`r` path
windows as vertices and compares nonadjacent windows by their outer endpoint
arc.  This supplies an exact recursive ledger for what every gap contraction
forgets.

## 3. The address simplex is the full seven-term object

For a tile `(a,b)`, use positive barycentric coordinates

```text
(u,g,w)=(b,a-b-1,n-a+1),          u+g+w=n,                (3.1)
```

or shifted coordinates

```text
(x,y,z)=(w-1,u-1,g-1),            x+y+z=n-3.              (3.2)
```

Here `x,y,z>=0`, reflection swaps `x,y`, and the full three-face membership is

```text
A iff x>0,             C iff y>0,             B iff z>0.  (3.3)
```

Thus `S_n` is the lattice-point simplex of weak compositions of `n-3` into
three parts.  For `n>=4`, at least one coordinate is positive, the three faces
cover, and the seven exact support strata have census

```text
A,B,C                         1 each,
AB,AC,BC                      n-4 each,
ABC                           binom(n-4,2).                (3.4)
```

The seam layer is the whole edge `z=0`, of size `n-2`.

The slot-preserving tile enumerator is the complete homogeneous polynomial

```text
T_n(X,Y,Z)=sum_(x+y+z=n-3) X^x Y^y Z^z=h_(n-3)(X,Y,Z).
```

Its generating function and recurrence are

```text
sum_(n>=3) T_n t^(n-3)=1/((1-Xt)(1-Yt)(1-Zt)),             (3.5)

T_n=(X+Y+Z)T_(n-1)
    -(XY+XZ+YZ)T_(n-2)
    +XYZ T_(n-3).                                         (3.6)
```

Here `T_3=1`; (3.6) starts at `n=6`, with `T_4,T_5` supplied by the
generating function.

This is the exact polynomial `A+B+C-D-E-F+G`, with the three letters and
all overlaps still typed.  Setting `X=Y=Z=1` produces THM-442's scalar third
difference and discards which geometric direction carried each bit.

There is a useful tournament-root chart

```text
d=w-u=x-y=n-(a+b-1),
1<=g<=n-2,
|d|<=n-g-2,
d == n-g (mod 2).                                         (3.7)
```

Reflection is `d->-d`.  On their domains the faces act by

```text
A:(g,d)->(g,d-1),
B:(g,d)->(g-1,d),
C:(g,d)->(g,d+1).                                         (3.8)
```

Therefore the full Laurent root spectrum

```text
P_n(q,y)=sum_(a,b) q^g y^d
```

satisfies the stronger seven-slot law

```text
P_n=(q+y+y^-1)P_(n-1)
   -(qy+1+qy^-1)P_(n-2)
   +qP_(n-3),                                             (3.9)
```

with

```text
sum_(n>=3) P_n t^n
 =q t^3/((1-qt)(1-yt)(1-y^-1 t)).                         (3.10)
```

Equation (3.9) starts at `n=6` with the coefficients of (3.10) as initial
data.

## 4. Internal deletion creates a repaired card and one seam bit

Delete path vertex `j`.  Relabel the remaining ordered vertices.  All new
adjacent path arcs are still forced except the shortcut between old vertices
`j+1` and `j-1`.  For `2<=j<=n-1`, define

```text
s_j(t)=t_(j+1,j-1)                                        (4.1)
```

and let `c_j(t) in X_(n-1)` retain every other nonadjacent arc.  Then
`(c_j(t),s_j(t))` determines the actual induced ordered tournament `T-j`.
The seam value says whether the inherited order is still a Hamiltonian path.

For endpoint deletions there is no seam.  For a lower tile `(p,q)`,

```text
c_j(t)_(p,q)=
  t_(p,q)         if j>p,
  t_(p+1,q)       if q<j<=p,
  t_(p+1,q+1)     if j<=q.                                (4.2)
```

Define the crossing window

```text
U_j={(p,q) in S_(n-1):q<j<=p}.                            (4.3)
```

On `U_j`, every repaired card agrees with the middle face:

```text
c_j(t)|_(U_j)=d_B(t)|_(U_j).                              (4.4)
```

The `U_j` cover `S_(n-1)`, their restrictions agree, and hence `d_B(t)` is
their unique Cech gluing.  Exact incidence counts are

```text
|U_j|=(j-1)(n-j)-1,
|U_j cap U_k|=(j-1)(n-k)                for j<k,
#{j:(p,q) in U_j}=p-q,
sum_j |U_j|=binom(n,3)-(n-2).                            (4.5)
```

Equivalently, the full gap tournament splits as

```text
X_n ~= X_(n-1)^B x F_2^(n-2 seam bits).                  (4.6)
```

The lower `B` face is its nonadjacent-arc chart; the seam word is its adjacent
arc row.  Complement flips both factors, and reflection reflects the lower
face and reverses the seam word.

Let `Blue_n=Fix(sigma_n)` and `Def_n=im(1+sigma_n)`.  The seam split gives

```text
0 -> F_2^(n-2) -> L_n -> L_(n-1) -> 0,
0 -> F_2^f_n   -> Blue_n -> Blue_(n-1) -> 0,
0 -> F_2^floor((n-2)/2) -> Def_n -> Def_(n-1) -> 0.       (4.7)
```

For lines, the exact lower-to-upper colour channel is

```text
BB=2^f_n b_(n-1)=b_n,             BK=0,
KB=(2^(n-2)-2^f_n)b_(n-1),
KK=2^(n-2)k_(n-1).                                      (4.8)
```

Here the first letter is upper colour.  The fresh seam defect has polynomial

```text
2^ceil((n-2)/2) (1+z)^floor((n-2)/2).                    (4.9)
```

## 5. Weighted deletion deck and mirror current

For one chord `(a,b)`, the nonendpoint deletion positions split into

```text
above a                 J_A=n-a=w-1=x,
strictly inside         g=a-b-1,
below b                 J_C=b-1=u-1=y.                   (5.1)
```

If `g>=2`, every internal deletion leaves a mutable lower tile and belongs to
the `B` role.  If `g=1`, the unique internal deletion turns this chord into
the seam, so define

```text
J_B=g 1[g>=2],          J_S=1[g=1].                       (5.2)
```

Then, cell by cell,

```text
J_A+J_B+J_C+J_S=n-2.                                    (5.3)
```

For a tiling, sum these weights over its one-bits.  Equations (4.2) give the
coordinatewise integer deck sum

```text
sum_(j=1)^n c_j(t)_(p,q)
 =(n-p)t_(p,q)+(p-q)t_(p+1,q)+q t_(p+1,q+1),             (5.4)
```

and the seam-corrected Hamming checksum

```text
sum_j |c_j(t)| + sum_(j=2)^(n-1) s_j(t)
 =(n-2)|t|.                                               (5.5)
```

At the carrier level this is

```text
(n-2)m_n=n m_(n-1)+(n-2).                                (5.6)
```

The exact basis totals are

```text
sum J_A=sum J_C=binom(n-1,3),
sum J_S=n-2,
sum J_B=binom(n,3)-(n-2),
sum (J_B+J_S)=binom(n,3).                                 (5.7)
```

The gap and signed boundary current are

```text
g=J_B+J_S,
d=J_A-J_C=n-(a+b-1).                                     (5.8)
```

Thus the full three-role deletion profile is equivalent to the tile address;
folding keeps `(g,|d|)`, and the sign of `d` is precisely the lost converse
sheet.  For general tilings put

```text
D(t)=sum_(a,b) d_(a,b)t_(a,b).                            (5.9)
```

Then

```text
D(sigma t)=-D(t),             D(kappa t)=-D(t),
sum d^2=2binom(n,4),          Var_uniform(D)=binom(n,4)/2,
sum g d=0,
sum g^2=(n-2)(n-1)^2 n/12.                                (5.10)
```

Blue implies `D=0`; the converse first fails at `n=6`, so signed current is a
useful edge coordinate but not a colour codec.

There is also an exact deletion-incidence matrix.  Put

```text
A_(j,p)=1[j is not an endpoint of chord p].                (5.11)
```

Every column sum is `n-2`.  If `e_j` is one at a path endpoint and zero
otherwise, then

```text
|row_j(A)|=m_(n-1)+1-e_j,
(AA^T)_(j,k)=m_(n-2)+2-e_j-e_k-1[|j-k|=1]   (j!=k).       (5.12)
```

Over `F_2`, for `n>=4`,

```text
rank(A)=n        for n odd,
rank(A)=n-1      for n even.                               (5.13)
```

For even `n`, the sole row relation is the all-card parity sum because every
chord survives an even `n-2` cards.  For odd `n`, complement-of-path
connectivity kills every left-null word.  For
`K_j(t)=sum_p A_(j,p)t_p`, the number of selected chords surviving in card
`j`, the selected-chord deck inventory is

```text
sum_t product_j z_j^(K_j(t))
 =product_(p=(a,b)) (1+product_(j notin {a,b}) z_j).       (5.14)
```

## 6. The half cube is a pair groupoid and defect recursion

Let `R_n` be the root half `d>=0` and `Bd_n` its fixed boundary `d=0`.  The
strict half maps literally to the previous half by

```text
(g,d>0) |-> (g,d-1),
R_n\Bd_n ~= R_(n-1).                                     (6.1)
```

Write

```text
h_n=|R_n|=floor((n-1)^2/4),
f_n=|Bd_n|=floor((n-1)/2).
```

Then

```text
m_n=h_n+h_(n-1),              f_n=h_n-h_(n-1).            (6.2)
```

Let `C_n=F_2^{R_n}` and let `tr:C_n->F_2^{Bd_n}` be fixed-boundary trace.
Restriction to the two mirror halves gives the canonical fibre product

```text
X_n ~= C_n x_(F_2^{Bd_n}) C_n.                            (6.3)
```

For each trace word `z`, a tiling is an arrow `(u,v,z)` between two interior
half-words.  Reflection is inversion `(u,v,z)->(v,u,z)`, blue tilings are the
identity arrows `u=v`, and composition adds defects

```text
delta(u,v)=u+v in F_2^{R_(n-1)}.                          (6.4)
```

Complement translates both endpoints and the trace by the all-one word.
After choosing the first half, (6.3) splits as

```text
X_n ~= Blue_n x Blue_(n-1),
L_n ~= (Blue_n/<1>) x Blue_(n-1).                         (6.5)
```

Thus the entire previous blue cube is the current defect alphabet.  Defect
zero is the new blue slice; every nonzero old-blue word labels an equal-mass
black slice.  Exact line populations are

```text
l_n=2^(h_n+h_(n-1)-1),
b_n=2^(h_n-1),
k_n=(2^h_(n-1)-1)2^(h_n-1),
P(blue)=2^(-h_(n-1)).                                     (6.6)
```

If defect weight counts mismatched reflection pairs, then

```text
#{lines of defect weight q}=2^(h_n-1) binom(h_(n-1),q).   (6.7)
```

The universal reflection-orbit line carrier has

```text
|Q_n|=2^(h_n-2)(2^h_(n-1)+1).                             (6.8)
```

At `n=14`,

```text
m_14=78,       h_14=42,       h_13=36,       f_14=6,
l_14=2^77,
b_14=2^41,
k_14=(2^36-1)2^41,
|Q_14|=2^76+2^40.                                         (6.9)
```

Equivalently, the `n=14` carrier is 64 trace layers of a
`2^36 x 2^36` arrow table.  The diagonal cells are blue; every nonzero
additive matrix diagonal is black.  Binary colour alone discards 36 defect
bits.

## 7. The master node-edge-tiling tensor

In the coordinates of Section 6, define

```text
c_n(z,u,delta)=pi_n(the tiling with halves u and u+delta and trace z). (7.1)
```

Converse merging gives

```text
c_n(z,u,delta)=c_n(z,u+delta,delta).                       (7.2)
```

Node mass and exact-defect line multiplicity are

```text
M_n(a)=#{(z,u,delta):c_n(z,u,delta)=a},                    (7.3)

E_(n,delta)({a,b})
 =1/2 #{(z,u):
        {c_n(z,u,delta),c_n(z+1,u+1,delta)}={a,b}}.         (7.4)
```

The factor `1/2` removes the two complement endpoints; loops retain their
usual two half-edges.  Summing (7.4) over `delta` reconstructs the whole merged
metagraph.

For a node `a`, let

```text
Def(a)={delta:c_n(z,u,delta)=a for some z,u}.               (7.5)
```

Then

```text
Def(a)={0}                  pure blue,
0 in Def(a), Def(a)!={0}    mixed,
0 notin Def(a)              pure black / non-self-converse. (7.6)
```

Every nonzero-defect multiplicity in a node or edge fibre is even, by the
free translation `u->u+delta`.  The zero-defect node multiplicity is odd exactly
for a self-converse node.  This recovers THM-796's blue/black parity theorem
as the fixed-point formula for groupoid inversion.

For black `delta!=0`, complement and reflection act on `(z,u)` by

```text
kappa:(z,u)->(z+1_Bd,u+1_I),
sigma:(z,u)->(z,u+delta).
```

They generate a Klein four action.  If `(r,s)` is a Walsh frequency on the
boundary and interior words, its two character parities are

```text
(<r,1_Bd>+<s,1_I>, <s,delta>) in F_2^2.                   (7.7)
```

Tilings retain all sectors.  The reflection quotient retains only
reflection-even data, the complement quotient retains only complement-even
data, and their joint quotient `Q_n` retains only the `(0,0)` sector.  The
merged-node map `pi_n` factors through the reflection quotient and then
forgets still more marked-path and tournament-label data; it need not be
complement-even because complement exchanges the two line endpoints.  The
`(1,1)` sector is simultaneous endpoint/converse holonomy.  This is a literal
four-component carrier, not a four-dimensional visual metaphor.

## 8. The half recursions are invariant-ring covers

Reflection swaps the simplex variables `X,Y` and fixes `Z`.  Therefore

```text
Q[X,Y,Z]^(X<->Y)=Q[s=X+Y,p=XY,Z],
deg(s,p,Z)=(1,2,1).                                       (8.1)
```

Its Hilbert series is

```text
1/((1-t)^2(1-t^2))=1/((1-t)^3(1+t)),                      (8.2)
```

whose degree-`n-3` dimension is `h_n`.

If `n` is even, `n-3` is odd, so every invariant monomial is divisible by
`s` or `Z`.  Their two-set monomial cover is exactly

```text
h_n=2h_(n-1)-h_(n-2),                    A+B-C.            (8.3)
```

If `n` is odd, the missing pure monomial `p^((n-3)/2)` forces a third corner.
The monomial cover by `s,Z,p` is

```text
h_n=2h_(n-1)-h_(n-2)+h_(n-2)
    -2h_(n-3)+h_(n-4),                  A+B-C+D-E-F+G.     (8.4)
```

The overlap `C=sZ` and corner `D=p` have the same scalar degree but are
different ideals.  Cancelling their sizes erases the exact geometry.

The weighted folded root spectrum

```text
H_n(q,y)=sum_(u<=w) q^g y^(w-u)                            (8.5)
```

has generating function

```text
sum_(n>=3) H_n(q,y)t^n
 =q t^3/((1-t^2)(1-qt)(1-yt)).                            (8.6)
```

Hence

```text
H_n=(q+y)H_(n-1)+(1-qy)H_(n-2)
    -(q+y)H_(n-3)+qyH_(n-4).                              (8.7)
```

This homogeneous recurrence starts at `n=7`; smaller coefficients come
directly from (8.6).

At even output sizes it additionally factors to

```text
H_n=(q+y)H_(n-1)-qyH_(n-2).                              (8.8)
```

Equation (8.8) holds for every even `n>=4`.

Equation (8.7) is the weighted odd seven-slot chart: the `+H_(n-2)` corner
and `-qyH_(n-2)` overlap become equal only after setting `q=y=1`.

The scalar curvature ledger is

```text
Delta^2 h_n=(1-(-1)^n)/2,
Delta^2 f_n=(-1)^(n+1),
Delta^2 m_n=1,                                            (8.9)

sum m_n t^n=t^3/(1-t)^3,
sum f_n t^n=t^3/((1-t)^2(1+t)),
sum h_n t^n=t^3/((1-t)^3(1+t)).                           (8.10)
```

Thus the full staircase is two adjacent half phases,
`m_n=h_n+h_(n-1)`; adjacent summation kills the half law's mirror eigenvalue
`-1`.

## 9. Exact radius-one node, tiling, and edge classification

Let `e_(a,b)` be the tiling with exactly one flipped tile.  THM-513 proves
that its ordinary tournament class determines the tile.  In (3.7), converse
sends `(g,d)` to `(g,-d)`.  Therefore

```text
ordinary one-flip node       (g,d),
merged one-flip node         (g,|d|),
reflection-orbit line packet (g,|d|).                     (9.1)
```

Moreover,

```text
C3(e_(a,b))=g,
H(e_(a,b))=1+2^g,
b=(n-g-d)/2,
a=(n+g-d+2)/2.                                            (9.2)
```

At fixed gap `g`,

```text
ordinary tilings/nodes       n-g-1,
merged nodes/Q packets       floor((n-g)/2),
blue/SC packet               1[n-g even],
black packets                floor((n-g)/2)-1[n-g even].  (9.3)
```

On the mirror spine `d=0`, the canonical complement line is one blue line.
Off the spine, the two roots `(g,+d),(g,-d)` give two black parallel lines
with the same merged-node boundary.  Hence the complete atomic sub-metagraph
has

```text
m_n literal root lines,
h_n merged root nodes and reflection packets,
f_n blue singleton packets,
h_n-f_n black double-line packets.                         (9.4)
```

This is an exact recursive relation among the user's three tracked sorts.
The only radius-one projected loop is the `n=4` apex blue line; every
radius-one line at `n>=5` has unequal endpoint `C3`.

At `n=14`, the merged packet counts by `g=1,...,12` are

```text
6,6,5,5,4,4,3,3,2,2,1,1,                                (9.5)
```

and the six even-`g` central packets are blue.  There are 78 literal root
lines, 42 merged nodes/packets, 6 blue lines, and 36 black packets containing
72 literal lines.  Their `H` values are `1+2^g`.

The universal `n=8` fixed-layer collision in THM-814 swaps tiles `(7,2)` and
`(6,3)`.  Both have `d=0` and lie in `ABC`, but their deletion profiles are

```text
(J_A,J_B,J_C,J_S)=(1,4,1,0),
(J_A,J_B,J_C,J_S)=(2,2,2,0).                              (9.6)
```

Thus the weighted deck separates all sixteen listed xor-`0x02080` collisions
immediately.  It identifies the first positional moment in THM-814/825 as the
internal-deletion coordinate, not an arbitrary extra statistic.

The aggregation boundary is sharp.  At `n=5`,

```text
{(4,2),(5,2)} and {(4,1),(5,3)}
```

have the same total deck vector `(1,2,2,1)` but `C3=2,3`.  At `n=6`,

```text
{(6,2),(5,3)} and {(6,4),(5,1)}
```

have the same deck vector `(1,3,3,1)` and `C3=4`, but score multisets

```text
(0,2,3,3,3,4) and (1,1,2,3,4,4).                         (9.7)
```

So labelled root positions or full cards must be retained; additive role
totals do not determine a general node.

## 10. Tournament packet formulas inside metagraph lines

For distinct roots `p,q`, let

```text
eta_(p,q)=-1  if they share an upper or lower endpoint,
eta_(p,q)=+1  if the upper endpoint of one is the lower endpoint of the other,
eta_(p,q)=0   otherwise.                                  (10.1)
```

The fixed-path triangle count has the exact quadratic form

```text
C3(t)=sum_p g_p t_p + sum_(p<q) eta_(p,q)t_p t_q.         (10.2)
```

The signed packet graph has

```text
negative edges     2binom(n-1,3),
positive edges     binom(n-2,3).                           (10.3)
```

If `w_p=1[b=1]+1[a=n]`, then

```text
sum_(q!=p) eta_(p,q)=w_p-2g_p,
C3(kappa t)-C3(t)=n-2-sum_p w_p t_p.                      (10.4)
```

Thus the nonlinear packet interaction collapses along a complement line to
the Hamming mass on the two boundary legs: `2n-5` leg coordinates, with the
apex coordinate counted twice.  Under the uniform measure on `t in X_n`, define
`Delta C3=C3(kappa t)-C3(t)`.  Expanding (10.2) in independent bits, or
equivalently applying Parseval to its linear and quadratic parts, gives the
line ANOVA

```text
Var(C3)=(n^3-7n^2+20n-16)/32,
E[(Delta C3)^2]=(n-1)/2,
Var((C3(t)+C3(kappa t))/2)=(n-3)(n-2)^2/32.               (10.5)
```

These equations are another direct tournament-to-edge dictionary: the odd
Walsh part is endpoint motion on a line, while the even quadratic part is its
midpoint energy.

## 11. Historical reconciliation and corrected boundaries

At least five different operators have shared the seven-slot word:

```text
1. size finite difference of m_n;
2. Cech inclusion-exclusion for the A/B/C face cover;
3. anchored Boolean vertex-deletion sieve;
4. Walsh derivatives on the arc cube;
5. divisor-lattice Mobius inversion in exact-period LRC packets.
```

The sign word alone never identifies the carrier, restriction maps, or
preserved predicate.

For a deletion-compatible support-additive invariant

```text
phi(T)=sum_R w_T(R)
```

whose weights satisfy `w_(T-A)(R)=w_T(R)` whenever `R` is disjoint from `A`,
and deletable set `U=V\{anchor}`, the anchored sieve is exactly

```text
sum_(empty!=A subset U)(-1)^(|A|+1) phi(T-A)
 =phi(T)-sum_(R:U subset R)w_T(R).                         (11.1)
```

Hence it reconstructs a `k`-local invariant when `|U|>k`.  For `c3`, it is
exact from `n=5` onward.  At `n=4` it **undershoots**, rather than overshoots:

```text
sieve=c3-[the all-deletable triple is cyclic].             (11.2)
```

The exhaustive audit gives residual histogram `0:48, 1:16` over the 64
tournaments.  For `H`, the residual is not vague nonlinearity: in the
Grinberg-Stanley odd-cycle-packing expansion it is precisely the weighted
sector whose support contains every deletable vertex.

Likewise, `F_2^{h_n}` parametrizes reflection-fixed tilings, not the orbit set
`X_n/<sigma>`.  The correct orbit count is

```text
|X_n/<sigma>|=(2^m_n+2^h_n)/2.                            (11.3)
```

Weighted orbit enumeration uses `(2^m_n+2^h_n)/2` representatives rather than
`2^m_n`, for the exact speedup

```text
2/(1+2^(h_n-m_n)),
```

which approaches two but is smaller because reflection-fixed tilings remain.

The corrected odd geometric role order is `(A,B,D,C,E,F,G)`, with ordinary
`+++---+` Boolean Mobius signs.  The alphabetical `++-+--+` match with a
quadratic-residue word mod seven is a slot permutation mnemonic; it is not a
new cyclotomic theorem, and the positive `G` slot has no counterpart in
`chi_7(0)=0`.

## 12. Tournament Analysis and the open recursive classifier

Tournament Analysis uses the exact radius-one merged nodes `(g,|d|)` as
vertices, challenging the assumption that tournament vertices must be
runners or original arcs.  The pair observable is lexicographic comparison.
The two gauges put `g` or `|d|` first; the other coordinate breaks ties and
is the tie Hamiltonian path.  Each gauge is transitive: score histogram
`{0:1,...,h_n-1:1}`, no directed cycles, singleton SCCs, and one Hamiltonian
path.  The gauges nevertheless flip

```text
n       4  5  6  7  8  9  10  11  12  13  14
flips   1  3  9 19 38 66 110 170 255 365 511.             (12.1)
```

Thus cyclic depth and mirror/deletion position define non-equivalent,
non-monotone orders even on the shell where nodes, tilings, and line packets
coincide.

The exact next object is the tensor `c_n(z,u,delta)` in Section 7, or equivalently
the coloured pair-groupoid table.  THM-818's kernel-pair relation is the first
equality-witness layer needed when the nonlinear node labels are not closed
under arrow composition.  The sharp research question is whether a feasible
refinement makes these labels a coherent configuration with constant
composition numbers.  Binary blue/black colour cannot: two black arrows
compose to blue exactly when their full nonzero defects agree.

The incoming THM-828/832/838 chain supplies the first nontrivial action test
on this proposal.  Its 58 `n=9` reflection doubletons form a punctured
rank-four defect cube; every face still sees that cube, and the centered
`n=9 -> n=10` coordinate copy preserves literal rank four while exposing only
two dimensions to raw `S2`.  Literal `Q` and the endpoint-coupled node pair
remain injective on that bank.  This is consistent with the master tensor:
the full defect is transported, while a histogram is only a projection.  It
does not prove closure for arbitrary continued-fraction words.

THM-839 identifies the missing nonlinear stalk more precisely.  Reflection
has 388 endpoint-swap fixed post-B bases, but only 58 also satisfy every local
ordered-state balance character.  Each of the three nonlinear punctures is
killed by one two-toggle parity obligation.  Thus groupoid inversion symmetry
is necessary but not sufficient for the raw ordered sidecar: a recursively
closed node label must retain these transverse local characters or an
equivalent ordered-state witness.  THM-841 reserves the corresponding
rolewise continuation test, especially whether the gap face is distinguished
from the endpoint faces; its reserved title is not used here as a theorem.

The barycentric carrier has a formal `S_3` triality permuting endpoint,
internal-gap, and opposite-end roles.  Only the endpoint swap is tournament
converse.  Transport kernels for the other two swaps measure exactly how the
node quotient breaks endpoint/cut duality and are a new finite classifier to
compute.

The preservation boundary is explicit.  The gap tournament, seam tower,
defect, node tensor, and root packets preserve all fixed-path tournament and
metagraph relations proved here.  They destroy runner speeds, metric scale,
owner, threshold, carry, and continued-fraction transport.  Those LRC fields
must be attached as a separate stalk before any conclusion about loneliness
at `n=14` is valid.
