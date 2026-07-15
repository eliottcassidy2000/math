---
id: THM-812
title: Centered-Christoffel coordinate replication embeds complement lines and transports projected coloured edge cells
status: PROVED (general coordinate-copy and Mobius-pullback lemma) + FINITE-EXACT (the (q,p,s)=(3,2,1) action from n=5 to n=6, including coloured edges and first truncation obstruction)
source: codex-2026-07-15-S13
depends_on: [THM-778, THM-796, THM-801, THM-808]
related: [THM-781, HYP-6880]
verification:
  - 04-computation/continued_fraction_metagraph_edge_transport_codex_S13.py
  - 05-knowledge/results/continued_fraction_metagraph_edge_transport_codex_S13.out
  - 05-knowledge/results/continued_fraction_metagraph_edge_transport_codex_S13.json
---

# THM-812 — centered-Christoffel transport on the metagraph stalk

There is an explicit centered-continued-fraction coordinate replication

```text
Phi:X_5 -> X_6
```

that embeds all 64 fixed-path tilings, commutes with complement and staircase
reflection, and consequently embeds all 32 complement lines while preserving
blue/black colour.  It does not descend from the 10 bare source nodes to the
target nodes.  It does, however, descend injectively on all 20 projected
coloured edge cells.  Thus the edge fibre, not the bare isomorphism-class
node, is the first functorial metagraph carrier for this continued-fraction
action.

## 1. The first nontrivial centered schedule

THM-778's centered height and increment word are

```text
F^s_(q,p)(i)=ceil(qi/p+q/(2p)-s/2),
d_i=F^s_(q,p)(i+1)-F^s_(q,p)(i).
```

For `(q,p,s)=(3,2,1)` they are

```text
(F(0),F(1),F(2))=(1,2,4),       d=(1,2).                (1)
```

The Euclidean partial quotient is one, the transported phase is zero, and
there is no midpoint tie.  Write an `X_5` tiling as

```text
a=t(5,1),       c=t(4,2),
(h0,h1)=(t(3,1),t(4,1)),
(l0,l1)=(t(5,2),t(5,3)).                                (2)
```

Use `d` on the low leg and its reverse on the reflected high leg, and extend
the core by the palindromic word `(c,a,c)`.  Explicitly,

```text
target apex:       t(6,1)=a
target high leg:  (t(3,1),t(4,1),t(5,1))=(h0,h0,h1)
target low leg:   (t(6,2),t(6,3),t(6,4))=(l0,l1,l1)
target core:      (t(4,2),t(5,2),t(5,3))=(c,a,c).        (3)
```

In explorer tile order the resulting target-to-source coordinate map is

```text
rho=(0,1,2,2,3,0,4,5,4,5).                              (4)
```

It is surjective and satisfies

```text
rho sigma_6 = sigma_5 rho.                               (5)
```

## 2. Coordinate-copy lemma

Let `S` and `T` be finite binary-coordinate sets, let `rho:T->S` be
surjective, and define

```text
Phi_rho(x)_j=x_(rho(j)).                                 (6)
```

Then `Phi_rho` is injective and commutes with all-bit complement.  If `S,T`
carry involutions and `rho sigma_T=sigma_S rho`, it also commutes with their
coordinate reflections.  Therefore, for staircase tilings, it embeds
complement lines and preserves reflection-fixed (blue) versus nonfixed
(black) colour.

The proof is immediate but structurally useful.  Surjectivity recovers every
source bit from one of its copies.  Complement commutes with copying.  The
intertwining identity gives

```text
(sigma_T Phi(x))_j=x_(rho(sigma_T j))
                    =x_(sigma_S rho(j))
                    =(Phi sigma_S(x))_j.
```

Equations (4)--(5) prove every general claim above for (3).  The verifier
checks all tilings independently: image size `64`, with zero complement,
reflection, colour, or intertwining failures.

## 3. Exact coloured-edge transport

For a complement line `e={t,kappa t}` define its projected coloured edge cell

```text
P_n(e)=(colour(e), unordered node pair of its endpoints). (7)
```

Although (3) spreads the 10 source nodes across 23 target nodes, with support
size histogram

```text
{1:3, 2:1, 3:2, 4:1, 5:1, 6:2},                        (8)
```

the composite `P_6 Phi` is constant on every fibre of `P_5`.  Moreover its
20 values are distinct.  In atlas-rank coordinates the exact map is

```text
B(0,5)->B(0,15)     B(1,6)->B(7,31)
B(2,8)->B(2,24)     B(5,6)->B(24,25)
B(5,8)->B(23,25)    B(6,9)->B(23,32)
B(7,9)->B(7,23)     B(8,9)->B(16,32)

K(1,3)->K(9,13)     K(1,4)->K(4,12)
K(1,5)->K(26,26)    K(1,8)->K(26,28)
K(3,5)->K(3,15)     K(3,6)->K(10,11)
K(3,7)->K(13,22)    K(3,8)->K(3,24)
K(5,7)->K(9,29)     K(6,6)->K(22,33)
K(6,8)->K(10,24)    K(8,8)->K(21,22).                  (9)
```

The eight blue source cells are singleton literal lines.  Each of the twelve
black cells contains exactly two lines, and those two lines form one
staircase-reflection orbit.  Hence (9) transports the complete source edge
quotient, including the loop cells, without pretending that a bare node has
a target.

This part is finite-exact, not an all-size functor theorem: constancy in (9)
uses the `n=5,6` tournament atlases.  The all-size statement is the tiling and
line embedding supplied by the coordinate-copy lemma.

## 4. Pullback of the Boolean Mobius stalk

Write a multilinear target function as

```text
f(y)=sum_(B subset T) mu_f(B) product_(j in B)y_j.
```

After substituting (6), repeated copies square-free collapse.  Therefore the
anchored Boolean Mobius coefficients satisfy the exact general law

```text
mu_(f o Phi)(A)=sum_(B subset T: rho(B)=A) mu_f(B),       (10)
```

where `rho(B)` is the set image.  Thus CF replication acts on coefficients
by a finite pushforward on subset support; coefficient degree can fall when
several target variables copy the same source variable.

For (3), fix the literal core bit `c` and order the five source variables as

```text
(l0,l1,h0,h1,a).                                        (11)
```

Pull back every target-node indicator to this five-cube.  The numbers of
distinct node indicators retained through anchored degrees `0,1,...,5` are

```text
c=0:  2,5,10,13,13,13       (13 target nodes),
c=1:  2,5,10,14,15,15       (15 target nodes).           (12)
```

Hence degree at most three is exact over `c=0`, but over `c=1` it has one
collision, between target nodes 12 and 31.  Their truth sets in the binary
indexing (11) are respectively

```text
node 12: {27,29},               node 31: {15}.           (13)
```

Exactly three degree-four coefficients first separate them:

```text
subset                 bidegree (low,high,apex)   mu_12  mu_31
l0 l1 h0 h1                    (2,2,0)                0      1
l0 l1 h1 a                     (2,1,1)                1      0
l0 h0 h1 a                     (1,2,1)                1      0. (14)
```

This is a sharp obstruction to the degree-three anchored truncation for this
transport.  It is not a collision of the full stalk and does not refute the
`Omega+B2` codec: THM-809 proves that codec exact at `n=8`.

## 5. Odd/odd tie sanity

For `(q,p,s)=(3,1,1)`, the centered word is `d=(3)` and has the predicted
midpoint tie.  Copying each one-bit `X_4` leg three times and setting all three
new core bits equal to the apex gives a second map `X_4->X_6`.  It is injective
on all 8 tilings and reflection/complement equivariant; its four complement
lines form three projected edge cells, and their target cell is again a
function.  This confirms that the coordinate-copy algebra survives the
smallest tie case.  It does not eliminate THM-778's tie sidecar from a metric
wall decoder: simultaneous owner blocks still have to be retained there.

## 6. Preservation boundary and LRC meaning

The natural vertices in this calculation are not runners or tournament arcs.
We explicitly compared literal tilings, complement lines, projected coloured
edge fibres, fixed-core node obligations, and Mobius coefficient subsets.
The quotient in (9) preserves precisely the source coloured edge cell for this
one replication.  Bare node projection destroys the chosen path/core lift;
coefficient truncation destroys a quartic cross-leg/apex interaction.

Nothing here preserves the metric LRC loneliness predicate, owner assignment,
wall position, common gcd, or prime-sheet carry.  The contribution is instead
a rigorous bridge between two previously separate descriptions: a centered-
CF digit supplies a coordinate substitution, and (10) supplies its exact
action on the metagraph stalk.  A future LRC transport must take the fibre
product of this path/core action with THM-808's owner-labelled redundancy-root
action, rather than attach either scalar to a bare isomorphism node.
