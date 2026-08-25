---
id: THM-4057
title: "Stern--Brocot depth pullback and rational-edge tournament gauge"
status: >
  PROVED + VERIFIED-EXACT.  Reduced positive rationals are exactly the
  ordered arcs of the coprimality graph, with reciprocal as arc reversal;
  this raw object is not a tournament.  The reciprocal is the global
  Stern--Brocot mirror and preserves path depth.  Pulling that depth through
  the exact Berggren/Farey prefix code gives D(wx)=D(x)+2|w|-#A and turns the
  A-branch Walsh character into the Stern--Brocot checkerboard character.
  An explicit scale-invariant depth gauge gives a lawful nontransitive
  tournament, with infinite Pell and Fibonacci directed-cycle families.
  Trace/norm/signed-sheet and Gaussian-square formulas identify exactly what
  reciprocal reversal preserves and loses.  None of this proves an
  irrationality statement, Duffin--Schaeffer, LRC(14), or a tournament
  characterization of Berggren ancestry.
source: codex-khinchin-ds-rational-tournament-20260824
audit: >
  The exact companion checks the rational-edge typing, reciprocal path mirror,
  canonical continued fractions, trace/norm/sheet identities, Gaussian
  content, all Berggren words through depth eight, five additional rational
  roots through depth five, three distinct reflections, the tournament on
  initial natural-number segments through 60, and finite Pell/Fibonacci
  controls.  An independent audit rederived the all-depth cocycle and Walsh
  ray, caught the gcd/height loss and the distinction from HYP-2925's
  Euclidean quotient-count parity, and verified the stated scope.
depends_on:
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
  - THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone
  - THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit
related:
  - HYP-2212-quadratic-discriminant-carrier-pi-e
  - HYP-2628-lrc14-euler-copy-squarefree-profile
  - HYP-2925-lrc14-magnitude-aware-tournaments
  - THM-3509-reduced-fraction-harmonic-k4-face-and-fibonacci-unit-cassini-ray
  - THM-3744-pell-prefix-loneliness-constant-carry-exact-formula
  - THM-4056-divisor-phase-compiler-and-duffin-schaeffer-lcm-clock
script: 04-computation/stern_brocot_rational_edge_tournament_thm4057.py
output: 05-knowledge/results/stern_brocot_rational_edge_tournament_thm4057.out
script_sha256: 98406973d1a2c60dfe30b498d36c8cd30e35b8088e851d612b1db28dde9e1d27
output_sha256: 13102a7220e5e02adfc4eae5dc920ea7acb63ba284216c004e0b3fe9e4473692
---

# THM-4057 -- reciprocal rational edges, Berggren depth parity, and a lawful tournament gauge

**PROVED + VERIFIED-EXACT.**

The user's rational-edge picture has an exact core, but its faithful raw
object is a coprimality graph, not a tournament.  The nontrivial theorem is
what happens after the rational edge is retained long enough to compare its
three natural structures: Stern--Brocot path depth, Gaussian/Berggren
transport, and an antisymmetric tournament gauge.

## 1. Rational numbers are ordered coprime arcs, not tournament edges

Let `G_copr` be the graph on `N_(>0)` with

```text
{a,b} in E(G_copr)  iff  a!=b and gcd(a,b)=1.           (1)
```

### Theorem 1.1 (exact rational-edge dictionary)

Reduction gives a bijection

```text
Q_(>0) \ {1}
  -> {ordered arcs of G_copr},
a/b |-> a -> b.                                        (2)
```

Reciprocal is exactly arc reversal:

```text
a/b |-> b/a,             (a->b) |-> (b->a).            (3)
```

Quotienting `(2)` by reciprocal gives the undirected edges of `G_copr`.
On the first `N` natural vertices,

```text
|E(G_copr[N])|=sum_(q=2)^N phi(q),                      (4)
```

because the smaller coprime neighbors of `q` are the `phi(q)` positive
residues below it.

The dictionary proves the tournament obstruction immediately:

- keeping all reduced rationals gives both directions on each coprime pair;
- a noncoprime pair has neither direction;
- `1/1` is a loop rather than an edge;
- retaining only `a<b` gives the increasing, acyclic orientation of
  `G_copr`;
- completing every missing pair by the same order gives only the transitive
  tournament.

On a finite pairwise-coprime vertex set the `a<b` restriction is a complete
transitive tournament, but it contains no directed-cycle information.  Any
nontransitive tournament therefore requires an additional antisymmetric
observable.

There are also two distinct coprimality levels:

```text
gcd(a,b)=1                                               (5)
```

makes one rational node `a/b` primitive, whereas

```text
|ad-bc|=1                                               (6)
```

is Farey adjacency between two rational nodes `a/b` and `c/d`.  Conflating
`(5)` and `(6)` deletes the Stern--Brocot geometry one is trying to study.

## 2. The global Stern--Brocot mirror and its depth

Give each positive rational its canonical simple continued fraction

```text
x=[a_0;a_1,...,a_r],             a_r>=2 if r>0,         (7)
```

and define

```text
D(x)=a_0+a_1+...+a_r-1.                                 (8)
```

This is the edge depth in the full Stern--Brocot tree rooted at `1`; it is
also the number of unit subtractions in the subtractive Euclidean algorithm.

The coordinate swap

```text
J=[0 1;1 0],                 J(a,b)=(b,a)               (9)
```

acts projectively by `x|->1/x`.  It swaps every left/right path symbol, so

```text
D(1/x)=D(x).                                            (10)
```

Geometrically `(9)` is reflection across `a=b` in endpoint coordinates,
`t|->-t` for `t=log x`, and reflection in the unit semicircle in the upper
half-plane model.  It is not the affine line reflection `x|->2-x`.

Two repo-local reflections must remain separate:

```text
j(x)=1-x                                                (11)
```

is the depth-preserving mirror of the `(0,1)` subtree rooted at `1/2`, while

```text
H=[-1 1;1 1],       H^2=2I,       det(H)=-2,
tau(x)=(1-x)/(1+x)                                    (12)
```

is THM-3357's Pythagorean leg/outer-branch reflection.  It needs a
content-two normalization and is not an integral tree automorphism.

## 3. Trace, norm, and the signed reciprocal sheet

For an ordered endpoint pair `(p,q)`, put

```text
T=p+q,                P=pq,                Delta=q-p.   (13)
```

Then

```text
T^2-4P=Delta^2.                                          (14)
```

Reciprocal reversal fixes `(T,P)` and sends `Delta|->-Delta`.  The square in
`(14)` is not enough to retain orientation; the signed sheet `Delta` is the
missing coordinate.  Any two of `T,P,Delta` reconstruct the endpoint pair up
to the corresponding finite branch:

```text
p=(T-Delta)/2,              q=(T+Delta)/2,
{p,q}=roots of X^2-TX+P,
{p,-q}=roots of X^2+Delta X-P.                          (15)
```

Also

```text
gcd(p,q)=1  iff  gcd(T,P)=1.                            (16)
```

Indeed, a prime dividing `T` and `P` divides both endpoints, and the converse
is immediate.

Gaussian squaring becomes

```text
Psi(p,q)=(q^2-p^2,2pq,p^2+q^2)
        =(T Delta,2P,T^2-2P).                           (17)
```

Reversal negates the first signed Gaussian coordinate and fixes the other
two.  If `gcd(p,q)=1` and `p,q` have opposite parity, `(17)` is a signed
primitive Pythagorean triple and the first coordinate is the signed odd leg;
taking its absolute value gives the positive primitive triple.  If both are
odd, the raw triple has content two; coprimality alone is insufficient.

For two spinors `u,v`, THM-3333 gives the Lorentz pairing

```text
<Psi(u),Psi(v)>_L=2det(u,v)^2.                          (18)
```

Simultaneous reciprocal reflection negates `det(u,v)` and preserves `(18)`.
Thus the positive/unordered Pythagorean image loses exactly the orientation
needed for a tournament chirality; a signed spinor sheet is mandatory.

Equations `(13)--(15)` are the same trace/norm/discriminant functor used in
[HYP-2212](../../05-knowledge/hypotheses/HYP-2212-quadratic-discriminant-carrier-pi-e.md)
for `e` and `pi`.  That is an exact shared quotient mechanism, not evidence
that `e+pi` is irrational: the individual status of `e+pi` remains open.

## 4. Stern--Brocot depth is the Berggren A-Walsh character

Use the THM-2596/THM-3357 positive parameter convention

```text
u=(m,n)^T,                    0<m<n,
x=m/n,                        Psi(u)=(n^2-m^2,2mn,n^2+m^2),
A=[0 1;-1 2],     B=[0 1;1 2],     C=[1 0;2 1].         (19)
```

On `(0,1)`, THM-2596 proves the exact Farey-prefix identities

```text
A=rho,
B=lambda o rho o j,
C=lambda^2,                                            (20)
```

where `lambda,rho` are the two positive Farey children and `j(x)=1-x`.
Each child adds one Stern--Brocot edge and `j` preserves depth.  Therefore

```text
D(Ax)=D(x)+1,
D(Bx)=D(x)+2,
D(Cx)=D(x)+2.                                          (21)
```

### Theorem 4.1 (all-root depth cocycle)

For every reduced `x in (0,1)` and every chronological Berggren word `w` of
length `h`,

```text
D(wx)=D(x)+2h-#A(w).                                   (22)
```

This follows by adding the three increments in `(21)`; word order is
irrelevant only for this additive depth statistic, not for the node itself.

THM-3357 equation (37), with `A,B,C` corresponding to its `L,M,R`, proves

```text
sum_(|w|=h)(-1)^(#A(w)) w u = C^(2h)u.                 (23)
```

Taking parity in `(22)` gives the exact pullback

```text
sum_(|w|=h)(-1)^(D(wx)) w u
  =(-1)^(D(x)) C^(2h)u.                                (24)
```

Thus the previously branch-defined `A`-Walsh character is precisely the
Stern--Brocot checkerboard depth character pulled back to the Berggren tree.

At the standard root `u_0=(1,2)`, `D(1/2)=1` and

```text
D(w(1/2))=1+2h-#A(w),
sum_(|w|=h)(-1)^(D(w(1/2)))w u_0
  =-C^(2h)u_0=(-1,-4h-2).                              (25)
```

Since the difference between the numbers of words with even and odd `#A` is
`(2-1)^h=1`, the standard level has

```text
# odd Stern--Brocot depth =(3^h+1)/2,
# even Stern--Brocot depth=(3^h-1)/2.                   (26)
```

This grading retains only depth parity.  It does not recover branch order,
Berggren ancestry, the sibling minors, or THM-3357's intrinsic norm
tournament.

## 5. A lawful scale-invariant tournament completion

For distinct positive natural numbers `a,b`, set `g=gcd(a,b)` and define

```text
sigma_D(a,b)
 =sgn(b-a)(-1)^(D((a/g)/(b/g))).                        (27)
```

Orient `a->b` iff `sigma_D(a,b)=+1`.

### Theorem 5.1

Equation `(27)` defines a tournament on every finite subset of positive
naturals and an infinite tournament on `N_(>0)`.  It is invariant under common
positive scaling.  Its first directed triangle on the initial natural-number
segments is

```text
2 -> 1 -> 5 -> 2.                                      (28)
```

Indeed, depth is reciprocal-even by `(10)`, while `sgn(b-a)` is
antisymmetric, so

```text
sigma_D(b,a)=-sigma_D(a,b) in {+1,-1}.                  (29)
```

Reduction by `g` proves scale invariance.  Finally

```text
D(2/1)=1,       D(1/5)=4,       D(5/2)=3,              (30)
```

which gives `(28)`.  There is no directed triangle on `{1,2,3,4}`, as the
finite exact audit verifies.

This is a mathematically lawful completion, but it is an **added gauge**.  It
changes existing natural-order coprime arcs—already `1/2` becomes `2->1`—and
reducing by `g` destroys the original gcd and height.  Those quantities must
remain sidecars for any arithmetic target.

The gauge is also distinct from the historical `CF-parity` experiment in
HYP-2925.  That code counts divisions in the ordinary Euclidean algorithm;
`(27)` uses subtractive/Stern--Brocot path length `sum a_i-1`.  Shared
continued-fraction vocabulary does not identify the two tournaments.

The companion records, as **FINITE-EXACT** data only, directed-triangle counts
on orders `3..20`:

```text
0,0,1,2,2,9,17,22,32,38,52,61,74,98,122,157,176,213.   (31)
```

No closed form or reconstruction claim is made for `(31)`.

## 6. Infinite Pell and Fibonacci cycle families

Let

```text
P_0=0, P_1=1, P_(n+1)=2P_n+P_(n-1),
F_0=0, F_1=1, F_(n+1)=F_n+F_(n-1).                     (32)
```

Their adjacent ratios have canonical continued fractions, for `j>=2` and
`n>=2`,

```text
P_j/P_(j-1)=[2;2,...,2]       (j-1 digits),
F_(n+1)/F_n=[1;...;1,2]       (n-2 copies of 1, then 2),
D(P_j/P_(j-1))=2j-3,
D(F_(n+1)/F_n)=n-1.                                   (33)
```

The recurrences prove the continued fractions by induction.  Modulo two,

```text
P_j is even iff j is even,
F_n is even iff 3|n.                                  (34)
```

For every integer `M>1`, `(8)` gives `D(M)=M-1`; hence `(27)` orients

```text
1->M when M is odd,          M->1 when M is even.       (35)
```

Combining `(33)--(35)` proves two infinite families.

### Corollary 6.1 (Pell cycles)

For every odd `j>=3`,

```text
1 -> P_j -> P_(j-1) -> 1.                              (36)
```

The first rows are

```text
(1,5,2), (1,29,12), (1,169,70), (1,985,408),...        (37)
```

The first row is exactly `(28)` and uses the first three nonzero Pell
numbers.  It links the new tournament gauge to the Pell-prefix family solved
in [THM-3744](THM-3744-pell-prefix-loneliness-constant-carry-exact-formula.md),
without transferring that theorem's loneliness maximum.

### Corollary 6.2 (Fibonacci cycles)

For `n>=2`,

```text
n=5 mod 6:   1 -> F_n     -> F_(n+1) -> 1,
n=0 mod 6:   1 -> F_(n+1) -> F_n     -> 1.             (38)
```

Thus every node `F_(6k)` with `k>=1` (which is even) lies in the two adjacent
cycles

```text
1 -> F_(6k-1) -> F_(6k) -> 1,
1 -> F_(6k+1) -> F_(6k) -> 1.                          (39)
```

The first four are `(1,5,8)`, `(1,13,8)`, `(1,89,144)`, and
`(1,233,144)`.  This is a new depth-parity cycle law, not the unit-Cassini
matching quotient of THM-3509.

## 7. Khinchin, Duffin--Schaeffer, Apéry, and `e+pi`

The exact connections now separate cleanly.

### Khinchin

Reciprocal inserts or removes only the leading continued-fraction digit/zero.
The positive digit word is unchanged, and any asymptotic geometric mean is
unchanged by one finite insertion.  Khinchin's almost-everywhere geometric
mean is therefore reciprocal-even.  So are depth `(8)` and finite products
of all positive digits.  None can orient an edge without an odd sidecar such
as `sgn(b-a)` in `(27)`.  MISTAKE-231 separately forbids treating a finite
digit mean as a universal entropy.  The irrationality of Khinchin's constant
itself remains open; this theorem supplies no approximation sequence for it.

### Duffin--Schaeffer

At denominator `q`, the primitive centres `a/q`, `a<q`, are exactly the
incoming natural-order coprime edges `a->q`, numbering `phi(q)`.  Thus a
Duffin--Schaeffer layer is a weighted denominator star in `G_copr`, and
THM-4056 compiles finite collections of those stars into LCM clocks.
Reciprocal reversal changes the denominator/radius; it is not the native
circle reflection for the approximation event.  The endpoint graph here is
also not the bipartite GCD graph used in the Koukoulopoulos--Maynard
second-moment proof.

### Apéry-style irrationality

An approximant `p_n/q_n` is an oriented projective edge.  The consecutive
Casoratian

```text
p_n q_(n-1)-p_(n-1)q_n                                (40)
```

is its signed Farey area; reciprocal negates it.  But irrationality needs
more than a nonzero path or fast convergence: after clearing denominators by
`L_n^m`, the integer linear form must tend to zero while remaining nonzero.
The LCM clock of THM-4056 supplies the valuation address, not that decay
estimate.  MISTAKE-209 remains the firewall.

### `e+pi`

The map `(p,q)->(T,P,Delta)` is the same finite quotient as
`(e,pi)->(e+pi,e*pi,pi-e)`.  It explains why a signed discriminant sheet is a
faithful orientation coordinate.  It does not specialize the integer
tournament to `(e,pi)` or prove that any individual shadow is irrational.
The proved survivor remains HYP-2212's statement that at least two of
`e+pi,e*pi,e-pi` are transcendental.

## 8. Transfer contract and stopping boundary

| source | target/map | preserved | destroyed / required sidecar |
|---|---|---|---|
| reduced rational `a/b` | ordered coprime arc `a->b` | endpoints, scale-free ratio, reciprocal | common scale/gcd if starting from arbitrary pair |
| rational node | full Stern--Brocot path | exact Euclidean run word, depth, and reconstructible reduced endpoints | none |
| full Stern--Brocot path | scalar depth `D` | path length and checkerboard parity | run word and endpoints |
| reciprocal pair | `(T,P)` | unordered endpoints | signed sheet `Delta` |
| spinor `(p,q)` | Gaussian triple `(17)` | norm-cone point | reciprocal orientation after positive normalization |
| Berggren word | depth cocycle `(22)` | `#A` parity and total depth | word order, sibling minors, ancestry |
| depth coloring | tournament `(27)` | antisymmetry, scaling, depth parity | gcd, height, raw numerator-to-denominator direction |
| Khinchin statistic | reciprocal quotient | asymptotic digit mean | finite phase, order, owner, odd orientation |
| DS denominator star | LCM clock in THM-4056 | primitive centres and denominator label | infinite overlap independence |

Theorem 4.1 is the positive bridge: it identifies one existing Berggren Walsh
character with one exact Stern--Brocot grading.  Theorem 5.1 is a new lawful
tournament, not a claim that tournaments are faithful for every arithmetic
target.  The raw rational-edge graph, the sibling norm tournament of
THM-3357, and the depth-gauged tournament are three different objects.

## 9. Reproduction and hostiles

Run

```bash
python3 -B 04-computation/stern_brocot_rational_edge_tournament_thm4057.py
python3 -B -O 04-computation/stern_brocot_rational_edge_tournament_thm4057.py
```

Normal and optimized transcripts are byte-identical.  Hostile controls check:

- the missing noncoprime pair `{2,4}` and reciprocal two-arcs;
- full reciprocal `J`, subtree mirror `j`, and content-two `H` separately;
- odd/odd coprime parameters, whose Gaussian image has content two;
- all `3^0+...+3^8=9841` standard-root Berggren words and five nonroot
  chambers through depth five;
- the first directed cycle and the fact that `(27)` reverses the raw `1->2`
  edge;
- normal versus optimized execution without truth-bearing `assert` paths.

The computation establishes its declared finite consequences.  The all-word
and infinite Pell/Fibonacci claims are proved symbolically by `(20)--(24)` and
`(33)--(35)`, not inferred from the finite prefixes.

Nothing here proves irrationality or transcendence of Khinchin's constant or
`e+pi`, the Duffin--Schaeffer theorem, a Berggren-tree tournament
classification, LRC(14), JC(2), or any transfer among those open problems.
