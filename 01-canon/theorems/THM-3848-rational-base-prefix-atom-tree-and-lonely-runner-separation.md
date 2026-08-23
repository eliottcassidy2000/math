---
id: THM-3848
title: "Rational-base prefix atom tree and lonely-runner separation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY PROOF-AUDITED, with an explicitly
  CONDITIONAL ABC addendum.  Every coprime rational base p/q has an exact
  p-ary finite-prefix wall-atom tree and a saturated adjacent relation
  lattice.  Its mixed-power speed row has maximum loneliness
  floor((p+q)/2)/(p+q).  For p/q=3/2 this is exactly 2/5 at every nontrivial
  level, with exactly two maximizing times; the level-12 row is a primitive
  thirteen-speed positive control, not an LRC(14) reduction.  The distinct
  closed formal safe-tail shift has an exact renewal law, entropy log(3/2),
  binary-ultrametric dimension log_2(3/2), and is nonsofic; this is not the
  unknown Z-language.  Assuming ABC gives radical growth at odd carry steps
  and radical saturation of an explicit denominator-19 safe-prefix family,
  but no horizon or terminal bound.  Mahler's Z-number problem and LRC(14)
  remain open.
source: root + lrc-rank-transfer + thm3833-referee + encoding-wildcard + safe-language-referee / 2026-08-23
audit: >
  PASS.  The companion performs 3,083 active exact gates: it exhausts five
  coprime bases through four atom levels, checks every parent/remainder child
  fibre, audits adjacent relations and maximal-minor gcds through level six,
  independently maximizes fifteen small mixed-power rows by the complete
  triangular-wave breakpoint set, classifies the 3/2 maximizers through level
  five, and checks the level-12, denominator-19, even-run, and radical/valuation
  hostiles, generates the greedy equality boundary, and verifies the safe-word
  renewal through depth 22.  Normal and optimized runs byte-match the frozen
  output.  Separate read-only proof audits checked the all-level maximizer
  classification, THM-2228 boundary interpretation, equality-tail uniqueness,
  renewal decomposition, nonperiodicity, nonsoficity, and dimension scope.
depends_on:
  - THM-2228-mahler-three-halves-carry-tail-and-integral-stabilization
related:
  - THM-2047-phase-height-toric-arrangement-for-lrc
  - THM-2050-period14-top-germs-do-not-determine-global-loneliness
  - THM-2352-q-adic-prefix-residue-collision-spectrum
  - THM-3743-lonely-runner-polyhedron-khinchin-flatness-relation-reduction
  - THM-3833-abc-conditional-cube-radical-and-hyperbolic-power-finiteness
script: 04-computation/mahler_rational_base_prefix_lrc_atom_tree_thm3848.py
output: 05-knowledge/results/mahler_rational_base_prefix_lrc_atom_tree_thm3848.out
script_sha256: 4510dc9a122e5d22f1e1ed7134472b753cfd94740fcd93aa27e7d6880a961ab5
output_sha256: 3dfba51c1680e1d0ec94adb18cd604774f4ca441f3baa2551e5af7a9afe7d748
hash_basis: raw working-tree bytes (LF)
---

# THM-3848 -- the rational-base prefix tower is exactly solvable, but its two target predicates separate

**PROVED + VERIFIED-EXACT + INDEPENDENTLY PROOF-AUDITED, with a
CONDITIONAL ABC addendum.**  The word *prefix* is load-bearing throughout.
Nothing below constructs a Mahler Z-number or a counterexample to LRC(14).

## 1. The mixed-power prefix carrier

Fix coprime integers `p>q>=2`.  For `N>=0` put

```text
v_(N,n)=q^(N-n)p^n,              0<=n<=N,
V_N=(v_(N,0),...,v_(N,N)),
H_N={t in R/Z : {v_(N,n)t}<1/q for every n<=N}.       (1)
```

If `xi` is read modulo its prefix period `q^N` and `t=xi/q^N`, then

```text
{xi(p/q)^n}={v_(N,n)t},           0<=n<=N.             (2)
```

Thus `H_N` is exactly the normalized finite one-sided rational-base prefix
set.  Let `D_q(t)=qt mod 1` and `L_m={t:{mt}<1/q}`.  Since

```text
v_(N+1,n)=q v_(N,n)       (n<=N),
v_(N+1,N+1)=p^(N+1),
```

one has the exact pullback law

```text
H_(N+1)=D_q^(-1)(H_N) intersect L_(p^(N+1)).           (3)
```

For a fixed `t in H_N`, its `q` inverse sheets are `(t+epsilon)/q`.  The old
phases agree on all sheets, while the appended phases differ by
`p^(N+1)/q`.  Because `p` is a unit modulo `q`, those phases hit the `q`
half-open arcs `[j/q,(j+1)/q)` once each.  Exactly one sheet lies in
`[0,1/q)`.  Consequently

```text
D_q:H_(N+1) -> H_N is a set-theoretic bijection,
mu(H_N)=q^(-(N+1)).                                    (4)
```

The bijection is piecewise, not one affine inverse branch.  Its switching is
what creates the finer atom tree.

## 2. The full `p`-ary wall-atom tree and natural ordinal

Every boundary in (1) lies on the common grid

```text
Delta_N=p^N q^(N+1).                                  (5)
```

Membership is constant on each open grid atom.  Equations (4)--(5) show that
exactly `p^N` of the `Delta_N` open atoms are safe.  More strongly, every safe
level-`N` atom has exactly `p` safe children.

Indeed, write a parent index as `k` and subdivide its `q` pullback sheets at
the next wall grid.  For each `r in {0,...,p-1}`, the `q` candidates have
indices

```text
j=p k+p Delta_N epsilon+r,           0<=epsilon<q.      (6)
```

They differ by translation through `1/q`; the terminal mask in (3) selects
exactly one `epsilon`.  The residue `r=j mod p` labels that child, while
`k=((j-r)/p) mod Delta_N` recovers its parent.  Hence the rooted safe-atom
tree is canonically the full `p`-ary tree.

For a path word `r_1...r_N`, `0<=r_i<p`, its reversible level rank and global
shortlex ordinal are

```text
rank_N=1+sum_(i=1)^N r_i p^(N-i),
rho=1+sum_(d=0)^(N-1)p^d+sum_(i=1)^N r_i p^(N-i).      (7)
```

This realizes the user's proposed principle: the ordinal is the discrete
address, not the numerical evaluator.  The selected sheet `epsilon`, the
level, and the map (6) remain reconstruction data.  In the original `xi`
coordinate, one period has `p^N` safe atoms, each of width `1/(qp^N)`, and
total length `1/q`.

The unique point lift in (4) and the `p` atom children in (6) are not a
contradiction.  The parent is refined into `p` remainder subatoms; adjacent
subatoms may use the same sheet, and the selected sheet changes only where
the terminal mask requires it.

## 3. The complete adjacent relation lattice

For `0<=n<N`, let

```text
r_n=p e_n-q e_(n+1).                                  (8)
```

Then `r_n dot V_N=0`, the `N` rows are independent, and they generate the
full integer annihilator of `V_N`.  To see the last assertion, the maximal
minors of the bidiagonal row matrix, up to sign, are

```text
q^N, p q^(N-1), ..., p^(N-1)q, p^N.                  (9)
```

Their gcd is one.  The row lattice is therefore saturated.  It has the same
rank and rational span as `ker_Z(V_N)`, so the two lattices agree.  Each row
has `l1` norm `p+q`; for `p/q=3/2`, every adjacent relation has norm five.

This is a lossless relation coordinate for the mixed-power ray.  It is not a
badness certificate: THM-3743 says a hypothetical LRC(14) counterexample has
one short relation, but safe rows can have a whole short basis.

## 4. Exact loneliness of every mixed-power row

For a positive integer row `W`, write

```text
M(W)=max_(t in R/Z) min_(w in W)||wt||.                (10)
```

For every `N>=1`,

```text
M(V_N)=floor((p+q)/2)/(p+q).                           (11)
```

For the upper bound, the first two coordinates are a common scale times
`{q,p}`.  The lower envelope of the two triangular waves can have a maximum
only at a cusp or where their distances agree.  Equality of the distances
means `(p-q)t` or `(p+q)t` is integral.  The difference branch is no larger
than the sum branch; the latter has maximum
`floor((p+q)/2)/(p+q)`.  A simultaneous cusp gives the same value only when
`p,q` are both odd.  This proves the pair upper bound and hence the upper
bound in (11).

Put `s=p+q` and `h=floor(s/2)`.  Since `q` is a unit modulo `s`, choose `a`
with

```text
a q^N = h mod s.                                       (12)
```

As `p=-q mod s`, one has

```text
a v_(N,n)=(-1)^n h mod s.
```

Every coordinate at `t=a/s` therefore has distance `h/s`, proving the lower
bound.

### The `3/2` tower and all maximizing times

For

```text
W_N={2^(N-k)3^k:0<=k<=N},             N>=1,           (13)
```

equation (11) gives `M(W_N)=2/5`.  Let `b_N in {1,2,3,4}` be defined by

```text
2^(N-1)b_N=1 mod 5.                                    (14)
```

Then the complete maximizing set is

```text
G_(W_N)(2/5)={b_N/5,-b_N/5};
b_N=1,3,4,2,1,3,4,2,... .                              (15)
```

Every adjacent pair is `g_k{2,3}` with
`g_k=2^(N-k-1)3^k`.  The pair `{2,3}` reaches `2/5` only at `+/-1/5`.
Thus a full maximizer satisfies `5g_k t in Z` for every `k`.  Since the
endpoint scales `2^(N-1)` and `3^(N-1)` are coprime, `5t` is integral; the
first pair and (14) leave exactly the two times in (15).  At either time the
ordered phases alternate `2/5,3/5`.

At `N=12` this gives the primitive thirteen-speed row

```text
(4096,6144,9216,13824,20736,31104,46656,
 69984,104976,157464,236196,354294,531441),            (16)
```

with twelve independent norm-five relations and exact maximum `2/5`, attained
at `t=2/5,3/5`.  It has `3^12=531441` one-sided safe wall atoms and normalized
safe measure `2^-13`.  None of its speeds is divisible by five, so the
mod-five witness makes it a sieve-trivial positive control, far from the
current covering residual for LRC(14).

## 5. Why the exact LRC solution is not a Mahler solution

The centered LRC evaluator retains only `||phase||`; Mahler's finite prefix
retains the oriented half-circle predicate `phase in [0,1/2)`.  The maximizing
word `2/5,3/5,2/5,3/5,...` loses precisely that orientation coordinate: half
of its entries are on the wrong side.

The same boundary is visible through THM-2228.  The periodic carry address

```text
c=(01)^infinity
```

has

```text
Phi(c)=-2/(3^2-2^2)=-2/5,
Y_even=(2/3)^2/(1-(2/3)^2)=4/5,
Y_odd =(2/3)  /(1-(2/3)^2)=6/5.                       (17)
```

The canonical residues of `-2/5 mod 2^N` are

```text
A_N=(2^N b_N-2)/5.                                    (18)
```

They never stabilize: `5A_N+2=2^N b_N` is unbounded.  Moreover
`Y_odd>1`, so the strict one-sided tail test fails.  The algebraic cancellation
`Phi(c)+Y_even/2=0` must not be read as THM-2228's real reconstruction formula;
its safety and stabilization hypotheses are absent.  The carry word `01` is
also not the ordinary least-significant-bit word of `-2/5`, which is
`0110`-periodic.

Thus the same finite carrier has two exact evaluators:

```text
source:       normalized mixed-power prefix atom
LRC target:   minimum centered distance
Mahler target:all oriented phases lie in the lower half-circle
preserved:    speeds, phase word, adjacent relation lattice, prefix level
destroyed by centered quotient:
              orientation of each phase
additional Mahler sidecar:
              strict real tails plus eventual ordinary residue stabilization. (19)
```

This is a proved separation, not merely a failed analogy.

## 6. The other prefix tree: a nonsofic formal safe-tail shift

The torus wall atoms above must not be confused with the binary carry-word
language.  Put `r=2/3`, and for a binary sequence `c=c_0c_1...` define

```text
Y_i(c)=sum_(j>=0)c_(i+j)r^(j+1).                       (S1)
```

Let `P_m` contain the length-`m` binary words satisfying every strict
truncated suffix inequality

```text
2 sum_(j=0)^(m-i-1)c_(i+j)2^j 3^(m-i-1-j)<3^(m-i),
0<=i<m.                                                (S2)
```

Write

```text
K={c:Y_i(c)<=1 for every i},
S={c:Y_i(c)< 1 for every i}.                           (S3)
```

Bond `P_(m+1)` to `P_m` by deleting the last letter.  Appending zero
preserves every strict suffix inequality, so these bonding maps are onto.

Then

```text
inverse_limit_(last-letter truncation) P_m=K=closure(S),
K\S={c in K:sigma^n c=d for some n>=0},               (S4)
```

where `d=d_1d_2...` is the greedy equality word

```text
d=101000001001001010000000001000000100...,
1=sum_(k>=1)d_k(2/3)^k.                                (S5)
```

In particular, `K\S` is countable, but finite prefixes cannot detect its
removal.

To prove (S4), a finite sum in (S2) can never equal one: after clearing its
power of three, its numerator is even and its right side is odd.  Every prefix
of `K` is therefore in `P_m`, while every `w in P_m` extends as
`w0^infinity in S`.  Monotone convergence gives the inverse-limit equality
and density.

For equality-tail uniqueness, start at `z_0=1` and set

```text
d_k=floor((3/2)z_(k-1)),
z_k=(3/2)z_(k-1)-d_k.                                  (S6)
```

The digit is forced except at state `2/3`.  Inductively `z_k` is an odd
integer divided by `2^k`, so it never equals `2/3`.  Thus any admissible suffix
with value one is exactly `d`, proving the second part of (S4).

### Renewal, entropy, dimension, and nonsoficity

Let `a_m=|P_m|`, `a_0=1`, and `D(z)=sum_(k>=1)d_k z^k`.  First disagreement
with `d` gives the unique decomposition

```text
P_m={d_1...d_m}
    disjoint_union_(k<=m,d_k=1)
      {d_1...d_(k-1)0 v:v in P_(m-k)}.                 (S7)
```

An upward first disagreement at a zero of `d` violates the current suffix
bound; after a downward disagreement, any safe suffix remains safe.  Hence

```text
a_m=1+sum_(k=1)^m d_k a_(m-k),
A(z)=sum_(m>=0)a_m z^m=1/((1-z)(1-D(z))).              (S8)
```

The exact counts begin

```text
1,2,3,5,8,12,18,27,40,60,90,134,201,302,452,... .    (S9)
```

The word `d` is not eventually periodic.  If its tail had period `L`, some
state `z_N` in (S6) would have reduced denominator dividing the odd integer
`3^L-2^L`.  But `z_N` is dyadic and lies strictly between zero and one, an
impossibility.

Every admissible prefix `d_1...d_n` has a follower set whose unique
lexicographic maximum is `sigma^n d`.  Equal follower sets would make two
suffixes of `d` equal and hence make `d` eventually periodic.  The follower
sets are all distinct, so the closed shift `K` is **nonsofic**.

The lexicographic maximum assertion uses the same first-disagreement gate:
every `x in K` satisfies `x<=d` lexicographically.  An upward first
disagreement at `d_k=0` has current tail at least `2/3`, while the greedy tail
is strictly below `2/3`, forcing `Y_0(x)>1`.

Finally, `D` is analytic on `|z|<1` and `D(2/3)=1`.  The bound
`|D(z)|<=D(|z|)` excludes zeros of `1-D(z)` inside that circle; equality in
the triangle inequality, together with `d_1=1`, makes `z=2/3` its only zero
on the circle.  Since `D'(2/3)>0`, the zero is simple.  Singularity extraction
from (S8) gives

```text
a_m ~ C(3/2)^m,
C=9/(2D'(2/3))=1.5510451884...,
h_top(K)=log(3/2),
dim_H(K)=dim_H(S)=log_2(3/2),                           (S10)
```

where dimension uses the standard prefix-cylinder formula for a closed binary
subshift in the binary ultrametric.  The union of safe level-`m` cylinders has
fair Bernoulli/Haar mass `a_m/2^m~C(3/4)^m`; continuity from above makes `K`
Haar-null.  The strict set `S` has the same finite language and dimension
because only the countable boundary (S4) was removed; `S` is not itself a
closed shift.

This theorem does **not** call the Mahler Z-language nonsofic.  That language
is `S` intersected with the ordinary positive-integer stabilization condition
of THM-2228 and may still be empty.  At depth four, the carry-residue map sends
the eight safe words to

```text
{0,4,5,6,8,9,13,14} mod 16,                            (S11)
```

four ordinary intervals rather than carry-lex adjacency.  Cardinality,
prefix bonding, entropy, and ultrametric dimension survive; Euclidean
adjacency and terminality do not.  This safe-word tree has `a_m` nodes at
level `m`, whereas the mixed-power torus partition has `3^m` atoms.  No
conjugacy between them is asserted.

## 7. What ABC lawfully adds to the carry dynamics

This section assumes the standard ABC schema with constant `K_epsilon`, as
typed in THM-3833.  It is **CONDITIONAL**.

For every odd positive integer `a`, put

```text
b=ceil(3a/2)=(3a+1)/2.                                (20)
```

The packet `3a+1=2b` is pairwise coprime and `gcd(a,b)=1`.  ABC gives

```text
2b <= K_epsilon rad(6ab)^(1+epsilon)
   <= K_epsilon (6 rad(ab))^(1+epsilon),

rad(ab) >= (2b/K_epsilon)^(1/(1+epsilon))/6.           (21)
```

Hence, for each `delta<1`, only finitely many odd carry steps satisfy
`rad(ab)<=b^delta`.  Every positive ceiling orbit, defined by
`a_(n+1)=ceil(3a_n/2)`, is strictly increasing and has infinitely many odd
states: an eventually even orbit would make one positive integer divisible
by every power of two.  Therefore, conditionally,

```text
liminf_(n -> infinity, a_n odd)
 log rad(a_n a_(n+1))/log a_(n+1) >=1.                (22)
```

This is a radical floor on carry-one transitions, not a Z-number
obstruction.  The exact hostile `a_0=2^H` has `H` consecutive zero-carry
steps

```text
a_n=2^(H-n)3^n,             0<=n<=H,
rad(a_n a_(n+1))=6,         0<=n<H.                   (23)
```

The normalized equation on these steps has no nonzero additive term, so ABC
cannot bound `H`.

### An unbounded safe-prefix family forced to be radical-saturated

For every positive multiple `m` of 18 define

```text
D_m=(2^m-1)/19,
Q_m=(3^m-2^m)/19,
alpha_m=9*2^m/19,
A_m=floor(alpha_m)=9D_m.                               (24)
```

Since `2^m=1 mod 19` and `3/2=11 mod 19`,

```text
{alpha_m(3/2)^n}=[9*11^n]_19/19
                 in {9/19,4/19,6/19}<1/2,
0<=n<=m.                                               (25)
```

Its first `m` carries are `(100)^(m/3)`, and their cleared carry polynomial
is `9Q_m`.  Thus (24) supplies arbitrarily long unconditional finite safe
prefixes.  The two primitive ABC packets

```text
19D_m+1=2^m,
2^m+19Q_m=3^m                                          (26)
```

give

```text
rad(D_m) >= 2^(m/(1+epsilon))/(38 K_epsilon^(1/(1+epsilon))),
rad(Q_m) >= 3^(m/(1+epsilon))/(114 K_epsilon^(1/(1+epsilon))). (27)
```

Since `D_m~2^m/19` and `Q_m~3^m/19`, ABC conditionally forces, along
positive multiples of 18,

```text
log rad(D_m)/log D_m ->1,
log rad(Q_m)/log Q_m ->1.                              (28)
```

The conclusion is radical saturation, not finiteness: the prefixes in (24)
exist for every such `m`, and their integer starts vary with `m`.  Assuming
the disputed IUT-to-ABC chain would yield only the same conditional consumer;
this theorem assumes ABC directly and has no IUT dependency.

## 8. Radical support is also too coarse for LRC(14)

Let

```text
P=product_(prime ell<=121) ell,
L=lcm(1,...,121),
U=84P,                  V=84L,
S_U={1,...,11,13,U},    S_V={1,...,11,13,V}.           (29)
```

Both rows are primitive and have identical componentwise radical vectors,
because `rad(U)=rad(V)=P`.  Yet `U=55 mod 121`, while `V=0 mod 121`.  At
`t=10/121`, the distance numerators for `S_U` are

```text
10,20,30,40,50,60,51,41,31,21,11,9,55.               (30)
```

Thus `S_U` is strictly lonely there, with margin `9/121>1/14`; the last
runner of `S_V` is at zero.  In fact every denominator at most 121 divides
`V`, so no member of that denominator bank can witness `S_V`.

The radical address preserves prime support but destroys valuation depth and
residue phase.  In this hostile, `v_11` is the first decisive missing
coordinate.  This does not compare the two rows' global maxima; it proves
only the stated failure of radical support to preserve a named phase or a
bounded-denominator witness bank.

## 9. Boundaries, reproduction, and open consequences

- `N=0` is excluded from (11): its singleton row has maximum `1/2`.
- The atom theorem counts open wall-grid atoms, not connected components;
  adjacent safe atoms may coalesce at included boundaries.
- Half-open `[0,1/q)` is load-bearing for the unique point lift.
- An infinite path in the normalized prefix tree can encode a nonterminal
  `q`-adic state.  It need not come from one fixed positive real `xi`.
- The `p`-ary atom tree is not the safe carry-word tree of THM-2228; they use
  different partitions and have different level counts.
- ABC constants are unknown, and Section 7 is conditional even if every
  displayed finite identity is exact.

Reproduce the finite audits with

```bash
python3 -B 04-computation/mahler_rational_base_prefix_lrc_atom_tree_thm3848.py
python3 -B -O 04-computation/mahler_rational_base_prefix_lrc_atom_tree_thm3848.py
```

The proofs of Sections 1--6 and 8 are exact and all-level; Section 7 is
explicitly conditional on ABC.  The companion is a hostile and independent
finite audit.  No section proves existence or nonexistence of
a Mahler Z-number, bounds a Z-number candidate beyond the cited literature,
excludes any live LRC(14) covering row, proves ABC, validates IUT, or improves
the `1/14` frontier.  **QED.**
