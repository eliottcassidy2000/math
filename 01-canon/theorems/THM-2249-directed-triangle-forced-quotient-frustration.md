---
id: THM-2249
title: "Directed-triangle forced quotient frustration"
status: >
  PROVED + VERIFIED-EXACT. In THM-2221's equal-block pinned transport
  model with common quotient the directed triangle, every transport X has
  a forced core-reversal energy
  F(X)=sum_(i,j)X_ij X_(i+1,j-1). It is a literal subset of the core cost:
  G(X)>=F(X). Hall decomposition into permutation layers gives an exact
  rotation/reflection formula. The zero layer consists precisely of the
  three whole-block cyclic rotations. Every other transport pays the sharp
  floor phi(1)=3, phi(N)=3(N-1) for 2<=N<=4, and phi(N)=2N+1 for N>=5.
  This yields an exact rotation-reduction certificate for the full pinned
  response and a sharper singleton-pinning criterion. It does not compute
  the residual internal-block kernel above F or handle unpinned exchange.
source: codex-2026-07-25-directed-triangle-transport
depends_on:
  - THM-2221-tournament-context-cut-metric-and-pinned-transport-response
related:
  - THM-2183-order-join-is-an-exact-tournament-metric-product
  - THM-2195-transitive-quotients-exactly-control-universal-substitution-products
  - THM-2242-tournament-complement-transport-and-knot-kernel-green-rigidity
script: 04-computation/tournament_directed_triangle_frustration_thm2249.py
output: 05-knowledge/results/tournament_directed_triangle_frustration_thm2249.out
script_sha256: 64718d49b4588e0a46d60a506229f8fb2d6c84b7928fd6a6b8f94bb50b6b85b3
output_sha256: f48e0616f6cd7df45d09715feb087737b4e6d71e7e8b10167fe3871d476417cf
hash_basis: working-tree bytes (LF)
---

# THM-2249 -- directed-triangle transport has a sharp frustration floor

Use THM-2221's pinned equal-order model with common quotient the directed
triangle. Label its vertices by `Z/3Z`, with

```text
i -> i+1.                                                   (1)
```

Thus

```text
P=C_3[T_0,T_1,T_2],       Q=C_3[S_0,S_1,S_2],
|T_i|=|S_i|=N>=1.                                          (2)
```

For a core bijection `pi`, let

```text
X_ij=#{x in T_i:pi(x) in S_j}.                             (3)
```

Every row and column sum of `X` is `N`. Let `g(pi)` be its full core
reversal cost and

```text
G(X)=min_(pi:X_pi=X) g(pi)                                 (4)
```

as in THM-2221.

The transitive prefix tax in THM-2242 cannot survive literally: the directed
triangle has three nonidentity-labelled cyclic transports with zero quotient
cost. The correct replacement is a quadratic frustration energy whose zero
set retains exactly those automorphisms.

## 1. The forced quotient energy

All indices below are modulo three. Define

```text
F(X)=sum_(i,j in Z/3Z) X_ij X_(i+1,j-1).                   (5)
```

Then

```text
g(pi)>=F(X_pi),
G(X)>=F(X).                                                (6)
```

Indeed, fix the source quotient arc `i->i+1`. A source vertex from `T_i`
mapped into `S_j` and a source vertex from `T_(i+1)` mapped into `S_(j-1)`
form a reversed pair, because

```text
j-1 -> j                                                   (7)
```

in the target quotient. There are exactly

```text
X_ij X_(i+1,j-1)                                          (8)
```

such pairs. The nine pair families in (8) are disjoint. They include only
pairs whose source labels and target labels are both distinct, so their
orientations are forced by the two quotient triangles and are independent
of every internal block orientation. All omitted pair classes contribute
nonnegative cost. This proves (6).

Thus `F` is not an analogy or a relaxation imported from another problem.
It counts a literal, canonically specified subset of the reversals in every
bijection with transport `X`.

## 2. Hall layers and the six-type table

View `X` as the multiplicity matrix of an `N`-regular bipartite multigraph
between source and target labels. Hall's theorem gives a perfect matching.
Delete it and repeat. Hence

```text
X=sum_(t=1)^N P_(sigma_t),       sigma_t in S_3.            (9)
```

This is the integral Birkhoff decomposition, proved here by repeated
matching rather than by a convexity claim.

For two permutation layers define

```text
f(sigma,tau)
 =#{i:tau(i+1)=sigma(i)-1}.                               (10)
```

Expanding (5) through (9) gives

```text
F(X)=sum_(s,t) f(sigma_s,sigma_t).                        (11)
```

Write the three rotations and reflections as

```text
R_a(i)=i+a,             H_a(i)=a-i,        a in Z/3Z.     (12)
```

The complete interaction table is

```text
f(R_a,R_b)=3 if b=a+1, and 0 otherwise;
f(R_a,H_b)=f(H_b,R_a)=1;
f(H_a,H_b)=3 if a=b, and 0 otherwise.                     (13)
```

Each line follows directly from the single congruence in (10). Suppose a
decomposition (9) contains `r_a` copies of `R_a` and `s_a` copies of `H_a`;
put `r=sum_a r_a`, `s=sum_a s_a`, so `r+s=N`. Equations (11)--(13) give the
exact layer formula

```text
F(X)
 =3(r_0 r_1+r_1 r_2+r_2 r_0)
  +2rs
  +3(s_0^2+s_1^2+s_2^2).                                (14)
```

Although a transport can have more than one matching decomposition, the
right side always equals the intrinsic expression (5). No uniqueness of
the layer counts is asserted or needed.

## 3. Exact zero layer and sharp positive floor

Equation (14) first gives

```text
F(X)=0
 iff X=N P_(R_a) for one a in Z/3Z.                       (15)
```

If `s>0`, the last two terms in (14) are positive. If `s=0`, any two
positive rotation multiplicities contribute to one of the three products
in the first term. Hence exactly one `r_a` equals `N`, proving (15).

There is also a sharp scale-bearing gap above zero. Define

```text
phi(1)=3,

phi(N)=min(3(N-1),2N+1)             for N>=2
      =3(N-1)                       for 2<=N<=4,
      =2N+1                         for N>=5.              (16)
```

Then every transport outside the three rotations satisfies

```text
F(X)>=phi(N),                                               (17)
```

and equality is attained for every `N`.

To prove the lower bound, first suppose `s=0` but at least two `r_a` are
positive. With their sum `N`, the first term of (14) is at least

```text
3(N-1).                                                     (18)
```

Equality is attained by `N-1` copies of one rotation and one copy of a
different rotation.

Now suppose `s>=1`. Dropping the first term in (14) gives

```text
F(X)>=2s(N-s)+3 sum_a s_a^2.                               (19)
```

For `s=1`, this is `2N+1`. For `s>=2`, use
`3 sum_a s_a^2>=s^2`. The right side is at least

```text
2Ns-s^2,                                                    (20)
```

which is at least `2N+1` at both endpoints `s=2,N` for `N>=3`, and hence
throughout the interval because (20) is concave. The remaining boundary
`N=s=2` is direct from (19), which gives at least six. Equality
`F=2N+1` is attained by `N-1` copies of one rotation and one reflection.
For `N=1`, the only nonrotations are reflections and have `F=3`. Comparing
(18) and (19) proves (16)--(17) with sharpness.

Combining (6) and (17),

```text
X notin {N P_(R_0),N P_(R_1),N P_(R_2)}
  implies G(X)>=phi(N).                                    (21)
```

This is the first positive universal core tax for the nontransitive
directed-triangle quotient. It grows linearly with block size even though
the three automorphism transports remain exactly free at quotient level.

## 4. Exact rotation-reduction certificate

Let `D=D_mu` be any nonnegative pinned exterior cut semimetric from
THM-2221. For each rotation put

```text
A_a=sum_i d_iso(T_i,S_(R_a(i))),

U_a(D)=A_a+N sum_i D(i,R_a(i)),
U_min(D)=min_a U_a(D).                                     (22)
```

Because `R_a` is an automorphism of the quotient, THM-2221's whole-block
formula is exact:

```text
G(NP_(R_a))=A_a.                                           (23)
```

Define the computable nonrotation certificate

```text
L_D=min_(X not a rotation) [F(X)+<X,D>].                   (24)
```

Then the full pinned response obeys

```text
min(U_min(D),L_D)<=R_(P,Q)(mu)<=U_min(D).                  (25)
```

The upper bound uses the best rotation. For every nonrotation, (6) gives
the lower affine-quadratic cost in (24), while the rotation costs are
exactly (22). Consequently

```text
L_D>=U_min(D)  implies  R_(P,Q)(mu)=U_min(D).              (26)
```

Since `D` is nonnegative, (17) gives the cheaper universal test

```text
U_min(D)<=phi(N)  implies  R_(P,Q)(mu)=U_min(D).           (27)
```

Thus the first nontransitive core does not require classifying all of
`G(X)` whenever either (26) or (27) fires. Only three labelled cyclic
reassignments survive.

There is also a sharpened uniform-pinning corollary. Let `mu_K` contain `K`
copies of every singleton cut word, so

```text
D_(mu_K)(i,j)=2K          for i!=j.                         (28)
```

Put `A_0=sum_i d_iso(T_i,S_i)`. Every nonidentity integer transport has at
least two off-diagonal units. Hence identity transport is optimal if

```text
A_0<=phi(N)+4K,
A_0<=A_a+6KN              for a=1,2.                       (29)
```

If all inequalities are strict, identity is the unique minimizing
transport matrix. The first line improves THM-2221's core-blind nonidentity
bound by the exact additive frustration `phi(N)`; the second line is the
unavoidable separate check on the two zero-frustration rotations.

## 5. Boundary and connection contract

The theorem deliberately stops short of claiming `G=F`.

```text
source:
  arbitrary split transports over a common directed-triangle quotient;

target:
  the quadratic forced-pair energy F and its Hall-layer formula;

map:
  retain only pairs with distinct source and target quotient labels;

preserved:
  a literal subset of every core reversal count, the complete zero layer,
  the sharp first positive layer, and certificates (26)--(29);

destroyed:
  source-internal pairs, target-internal pairs, and all internal block
  tournament structure;

needed sidecars:
  the three rotation assignment costs A_a, the cut metric D, and the
  residual kernel G-F when the certificate does not fire;

cheapest hostile controls:
  mixed rotations attain 3(N-1), while one reflection layer attains 2N+1.
                                                                  (30)
```

The faithful replacement for the transitive prefix order is therefore not
a cosmetic tournament on transports. It is a nonnegative quadratic
cycle-frustration energy plus three automorphism-sector sidecars. Unpinned
core/context exchange and larger strongly connected quotients remain open.

## 6. Independent exact audit

The companion enumerates every `3 x 3` nonnegative integer transport with
all margins `N` for `1<=N<=12` (`13,558` tables total). It verifies (14),
the three-element zero layer, and the sharp floor (16). It also checks all
`720` bijections of two six-vertex, three-block cores and all `64` choices
of source/target internal orientations (`46,080` full-cost checks), proving
computationally that the pairs counted by `F` are always present in the
actual reversal cost. Normal and optimized runs reproduce the stored
transcript exactly.

QED.
