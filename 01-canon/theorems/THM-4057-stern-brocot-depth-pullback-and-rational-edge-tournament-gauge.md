---
id: THM-4057
title: "Stern--Brocot depth pullback and rational-edge tournament gauge"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Reduced positive
  rationals are ordered coprime arcs, not tournament edges; a tournament on
  positive naturals is exactly a scale-decorated reciprocal selector.
  Reciprocal is the global Stern--Brocot mirror and preserves its continued-
  fraction depth. Pullback through the Berggren prefix code gives
  D(wx)=D(x)+2|w|-#A and the A-branch Walsh character. Separately, the
  Calkin--Wilf heap gives a natural ordinal: primitive Pythagorean parameters
  are exactly k=3,5 mod 6, with explicit Berggren transducers and radix
  inverses. A lawful depth-gauged tournament has infinite Pell/Fibonacci
  cycle families; branch-triangle gcds, parity counts, and the unique root
  shadow collision are exact. No LRC, Jacobian, irrationality, or arbitrary
  tournament-invariant transfer is asserted.
source: codex-khinchin-ds-rational-tournament-20260824
audit: >
  The primary companion checks rational-edge typing, reciprocal Stern paths,
  canonical continued fractions, trace/norm/sheet identities, Gaussian
  content, Berggren depth words and roots, the lawful tournament through
  finite segments, and Pell/Fibonacci controls. The secondary and independent
  companions separately check Stern versus Calkin--Wilf reflection, the heap
  involution, mod-six Pythagorean section, Berggren ordinal transducers and
  inverses, reduction fibres, scale selectors, branch gcds, parity counts,
  the extreme-edge cycle rule, and the unique root collision.
depends_on:
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
  - THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone
  - THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit
  - THM-3509-reduced-fraction-harmonic-k4-face-and-fibonacci-unit-cassini-ray
  - THM-3756-odd-square-ordinal-berggren-affine-descent
  - THM-4056-divisor-phase-compiler-and-duffin-schaeffer-lcm-clock
related:
  - HYP-2212-quadratic-discriminant-carrier-pi-e
  - HYP-2628-lrc14-euler-copy-squarefree-profile
  - HYP-2925-lrc14-magnitude-aware-tournaments
  - THM-2294-anchored-plucker-tournament-and-kakeya-address-bank
  - THM-3382-fibonacci-ray-dual-index-harmonic-bifurcation-and-ternary-heap-addresses
  - THM-3744-pell-prefix-loneliness-constant-carry-exact-formula
script: 04-computation/stern_brocot_rational_edge_tournament_thm4057.py
output: 05-knowledge/results/stern_brocot_rational_edge_tournament_thm4057.out
script_sha256: 98406973d1a2c60dfe30b498d36c8cd30e35b8088e851d612b1db28dde9e1d27
output_sha256: 13102a7220e5e02adfc4eae5dc920ea7acb63ba284216c004e0b3fe9e4473692
secondary_script: 04-computation/calkin_wilf_berggren_ordinal_thm4057.py
secondary_output: 05-knowledge/results/calkin_wilf_berggren_ordinal_thm4057.out
secondary_script_sha256: 0ad8497d50ab4ec59568489b6dc8c8f4ed9243b64c36e6608ab369137e99cfa2
secondary_output_sha256: 00f13ea27cfb7ec1bcc44dd1375623e4f8803de642b19679f64419fa3d354e2b
independent_audit_script: 04-computation/calkin_wilf_berggren_ordinal_thm4057_independent_audit.py
independent_audit_output: 05-knowledge/results/calkin_wilf_berggren_ordinal_thm4057_independent_audit.out
independent_audit_script_sha256: a65dc37d2fd4a2a592ee3fcd0095ec5ba8ff0a42a46e88fd32bd169baee7fa82
independent_audit_output_sha256: b6eef90db4c150e3570454c1236145a6de9f5ef533aef2e282beb59d64c3a335
hash_basis: raw bytes
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

After deleting only a leading zero, the full projective Euclidean
coefficient word is reciprocal-invariant. The standard finite Khinchin word
instead deletes the integer part on both sides and can lose or gain one
positive digit; `3/5=[0;1,1,2]` versus `5/3=[1;1,2]` is the exact hostile.
Consequently an infinite asymptotic geometric mean is unchanged by the finite
prefix modification, but a finite positive-digit word or product need not be.
Stern--Brocot depth `(8)` is reciprocal-even. Neither it nor an asymptotic
digit mean can orient an edge without an odd sidecar such as `sgn(b-a)` in
`(27)`. MISTAKE-231 separately forbids treating a finite digit mean as a
universal entropy. The irrationality of Khinchin's constant itself remains
open; this theorem supplies no approximation sequence for it.

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

## 10. Complementary Calkin--Wilf ordinal and scale-selector layer

Three coordinates must remain distinct. Stern--Brocot depth is the sum-of-
partial-quotients path statistic used by the Berggren cocycle above.
Calkin--Wilf heap position is a breadth-first natural ordinal with a different
reciprocal involution. Finally, the dilation or gcd of a natural-number edge
is erased by reduction. The results below show how the three coordinates
interact without identifying them.

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** The theorem
implements the proposed “rational number as an edge” only after separating
four objects which otherwise collapse: completeness, direction, primitive
ratio, and dilation. No literature-priority or global-novelty claim is made.

### 1. Inheritance and type repair

[THM-2294](THM-2294-anchored-plucker-tournament-and-kakeya-address-bank.md)
is the orientation guardrail: antisymmetric signs orient, while symmetric
characters are colours. [THM-3509](THM-3509-reduced-fraction-harmonic-k4-face-and-fibonacci-unit-cassini-ray.md)
supplies the exact reduced-fraction/Pythagorean leg-swap fibre.
[THM-3756](THM-3756-odd-square-ordinal-berggren-affine-descent.md) supplies
the closest natural-number chart and the canonical warning that outer odd
rank is not a node address. [THM-3382](THM-3382-fibonacci-ray-dual-index-harmonic-bifurcation-and-ternary-heap-addresses.md)
is the index-map hostile: density and harmonic class change under a different
injection of the same tree into the naturals.

There are three distinct rational-edge objects.

1. All positive reduced rationals `p/q`, oriented `p->q`, give both arcs on
   every coprime pair, plus the loop `1->1`. This is a bidirected coprimality
   graph, not a tournament.
2. Restricting to `p<q` gives the upward acyclic orientation of the
   coprimality graph. Noncoprime pairs such as `{2,4}` are missing, so it is
   still not a tournament.
3. Orient every natural pair by size, `a->b` for `a<b`, and attach its
   reduced ratio and gcd. This is a complete transitive tournament whose
   rational label has a scale sidecar.

The first two constructions are useful graphs. Calling either a tournament
would erase a load-bearing failure of completeness or antisymmetry.

### 2. Every tournament is a scale-decorated reciprocal selector

Write every unordered natural pair uniquely as

```text
{a,b}=g{p,q},       g=gcd(a,b),       gcd(p,q)=1,       p<q. (CW1)
```

For a tournament `T` on `N_(>0)`, define

```text
epsilon_T(p,q;g)=0  if gp->gq,
epsilon_T(p,q;g)=1  if gq->gp.                            (CW2)
```

Then

```text
T  <->  arbitrary bit functions epsilon(p,q;g)             (CW3)
```

is a bijection. Tournament converse toggles every bit, and the selected
reduced rational changes from `p/q` to `q/p`.

A direction on the reduced rational independent of representative scale is
well-defined if and only if

```text
epsilon(p,q;g) is independent of g.                         (CW4)
```

Without (CW4), two multiples of one rational ray can select opposite reciprocal
nodes. Thus an arbitrary tournament is not a bit function on reduced
rationals alone; it is a **scale-decorated** selector.

If a directed edge `u->v` is labelled `u/v`, then every directed path has

```text
product_(i=0)^(r-1) v_i/v_(i+1)=v_0/v_r.                  (CW5)
```

Every directed-cycle product is therefore one. Rational magnitudes form an
exact multiplicative coboundary and cannot detect tournament curl; curl lives
in the selector bits.

#### Coprimality masks impose no finite graph-pattern restriction

Every finite simple graph occurs, up to relabelling, as an induced
coprimality graph. Given a desired graph, assign a private prime `p_ij` to
both vertices `i,j` for each desired **nonedge**, and give every vertex an
additional private prime. The products at two vertices have gcd greater than
one exactly at the prescribed nonedges and are all distinct. Consequently a
coprimality mask alone cannot force a finite tournament pattern; direction,
height, prime labels, or Euclidean ancestry is additional data.

#### Minimal completion hostile

On `[4]` the upward coprime graph misses only `{2,4}`. Completing it as
`2->4` gives the transitive tournament with

```text
(H,c_3)=(1,0),
```

while `4->2` gives a source over a directed triangle with

```text
(H,c_3)=(3,1).                                            (CW6)
```

The same coprime-edge data therefore admits different Hamiltonian-path and
cycle counts.

### 3. Exact Stern--Brocot reciprocal reflection

Represent a reduced fraction by its primitive column `(numerator,denominator)`.
At a Stern--Brocot interval address `w`, let `P_w` have the left and right
boundary columns. Put

```text
J  =[[0,1],[1,0]],
U_L=[[1,1],[0,1]],
U_R=[[1,0],[1,1]].                                        (CW7)
```

Then

```text
P_empty=J,       P_(wL)=P_w U_L,       P_(wR)=P_w U_R.     (CW8)
```

Since `J U_L J=U_R`, induction gives

```text
P_(bar w)=J P_w J,                                        (CW9)
```

where `bar w` interchanges every `L` and `R`. The mediant vector is
`P_w(1,1)^T`, so (CW9) swaps its coordinates:

```text
p/q  ->  q/p.                                             (CW10)
```

Thus arc reversal is literally reflection across the center of the
Stern--Brocot tree. It is not, by itself, a vertex permutation of a finite
tournament.

The positive pair monoid obeys the same mirror identities

```text
L(p,q)=(p,p+q),       R(p,q)=(p+q,q),
S L=R S,              S R=L S,                            (CW11)
```

with `S(p,q)=(q,p)`. Equation (CW11) does not identify the Stern--Brocot and
Calkin--Wilf parent maps.

### 4. The distinct Calkin--Wilf natural ordinal

The two trees enumerate the same positive rationals, but the same mixed word
need not name the same node:

```text
word LR:       Stern--Brocot 2/3,       Calkin--Wilf 3/2.  (CW12)
```

Let Stern's diatomic sequence be

```text
s(0)=0, s(1)=1,
s(2n)=s(n),       s(2n+1)=s(n)+s(n+1).                    (CW13)
```

The Calkin--Wilf heap node `k>=1` is

```text
C(k)=s(k)/s(k+1).                                         (CW14)
```

Write the binary expansion of `k` as `1w`, with `L=0,R=1`, and put
`ell=floor(log_2 k)`. Letterwise reciprocal reflection has the affine ordinal

```text
k*=3*2^ell-1-k.                                           (CW15)
```

Indeed, (CW15) complements the `ell` suffix bits and the pair recursion in
(CW11) swaps `s(k),s(k+1)`.

For the natural-order tournament on `[N]`, every edge has the unique lossless
label

```text
(a,b)=g(s(k),s(k+1)),
k>=2 even,       g s(k+1)<=N.                             (CW16)
```

Thus

```text
E(T_N) <-> {(k,g): k>=2 even, g s(k+1)<=N}.                (CW17)
```

Forgetting `g` gives the fibre `floor(N/s(k+1))` proved in
[THM-4056](THM-4056-divisor-phase-compiler-and-duffin-schaeffer-lcm-clock.md).
Reflection preserves maximum height and this fibre while reversing the arc.

Reciprocal rationals share the **projective Euclidean coefficient word**
obtained by deleting only a leading continued-fraction zero. They do not
generally share the standard finite Khinchin word, which always discards
`a_0`; THM-4056 gives the exact `3/5` versus `5/3` hostile.

### 5. Primitive Pythagorean parameters are `k=3,5 mod 6`

Two elementary inductions from (CW13) give

```text
s(k) is even  iff  3|k,
s(k)>s(k+1)   iff  k is odd,       k>1.                   (CW18)
```

Adjacent Stern values are coprime. Therefore the standard primitive
Pythagorean parameters

```text
m>n>0,       gcd(m,n)=1,       m-n odd                    (CW19)
```

are exactly

```text
(m,n)=(s(k),s(k+1)),       k=3 or 5 mod 6.                 (CW20)
```

The opposite-parity condition is essential: `(3,1)` is coprime but generates
the nonprimitive triple `(8,6,10)`.

For the standard Berggren parameter branches

```text
A(m,n)=(2m-n,m),
B(m,n)=(2m+n,m),
C(m,n)=(m+2n,n),                                           (CW21)
```

the induced natural-number maps are

```text
A(k)=2k-1,
B(k)=4k*+3=12*2^floor(log_2 k)-4k-1,
C(k)=4k+3.                                                 (CW22)
```

If the Calkin--Wilf word of `k` is `w=uR`, the word laws are

```text
A: uR -> uLR,
B: w  -> (bar w)RR,
C: w  -> wRR.                                              (CW23)
```

The child branch is visible in the terminal binary digits:

```text
A-image: k=1 mod 4,
B-image: k=3 mod 8, except the root k=3,
C-image: k=7 mod 8.                                       (CW24)
```

Consequently the inverse Berggren descent on the set (CW20) is

```text
k -> (k+1)/2                         if k=1 mod 4,
k -> ((k-3)/4)*                      if k=3 mod 8 and k>3,
k -> (k-3)/4                         if k=7 mod 8.          (CW25)
```

This is a complete radix descent to the root `k=3`. Under THM-3756's branch
convention, `A,B,C` are respectively its `L,M,R`.

#### Exact depth pullback

The root has Calkin--Wilf depth one. An `A` step adds one bit; a `B` or `C`
step adds two. At Berggren depth `d`, the nodes at Calkin--Wilf depth

```text
1+d+j
```

are therefore counted by

```text
C(d,j)2^j,       0<=j<=d.                                 (CW26)
```

The ordinal image (CW20) has natural density `1/3`, Dirichlet series for
`Re(s)>1`

```text
6^(-s)[zeta(s,1/2)+zeta(s,5/6)],                           (CW27)
```

and, as `X->infinity`, harmonic sum

```text
sum_(k<=X,k=3,5 mod 6)1/k
 = (1/3)log X + gamma/3 - pi*sqrt(3)/12
   + log(2)/3 - log(3)/12 + O(1/X).                         (CW28)
```

These are properties of the Calkin--Wilf ordinal. The canonical ternary heap
of THM-3382 gives the same abstract tree a different density and harmonic
class.

### 6. Odd-square coordinates and the Vieta light-cone carrier

The exact bridge from (CW19) to THM-3756's odd roots is

```text
q=m+n,       d=m-n,
m=(q+d)/2,  n=(q-d)/2.                                    (CW29)
```

It bijects primitive opposite-parity `(m,n)` with coprime odd `q>d>0` and
preserves the full primitive triple:

```text
(m^2-n^2,2mn,m^2+n^2)
 =(qd,(q^2-d^2)/2,(q^2+d^2)/2).                           (CW30)
```

For the oriented edge `d->q`, put

```text
Sigma=q+d,       Pi=qd,       Delta=q-d.                   (CW31)
```

Reciprocal reflection fixes `Sigma,Pi`, negates `Delta`, and

```text
Sigma^2-4Pi=Delta^2,
A=Pi,       2B=Sigma*Delta,       4C=Sigma^2+Delta^2.      (CW32)
```

Thus the generic sum/product/discriminant carrier from the `e,pi` work has an
exact Pythagorean realization: the symmetric coordinates retain the unordered
edge, while the discriminant sign is the orientation sidecar. Reversal sends

```text
(A,B,C)->(A,-B,C).                                        (CW33)
```

If `g=gcd(q,d)>1` is odd, (CW30) has exact content `g^2`; reducing `d/q`
collapses that scaled Berggren component. For an odd denominator `q`, exactly
`phi(q)/2` incoming reduced numerators are odd, so THM-3756's odd/odd chart is
precisely half of that denominator shell.

### 7. Berggren branch triangles and the minimal cycle obstruction

For a parent `(m,n)` satisfying (CW19), join the parent edge, one child edge in
(CW21), and the remaining closure edge on their three natural endpoints. The
closure gcd is

```text
gcd_A=gcd_B=gcd(n,2),       gcd_C=gcd(m,2).                (CW34)
```

Indeed,

```text
gcd(n,2m-n)=gcd(n,2m+n)=gcd(n,2),
gcd(m,m+2n)=gcd(m,2n)=gcd(m,2).                            (CW35)
```

At Berggren depth `d`, the counts of nodes with `n` odd/even are

```text
(3^d+(-1)^d)/2,       (3^d-(-1)^d)/2.                     (CW36)
```

This follows from the two-state parity transfer: `A,B` toggle the state and
`C` preserves it. Hence the marked branch incidences split as

```text
primitive K3       = (3^(d+1)+(-1)^d)/2,
scale-two closure  = (3^(d+1)-(-1)^d)/2.                  (CW37)
```

In the natural-order tournament each branch triangle is transitive. Flipping
the unique edge joining its least and greatest natural vertices creates a
directed 3-cycle; flipping either short edge alone leaves a transitive
triangle. The extreme edge is the closure in branches `A,B` and the child
edge in branch `C`. Global reflection flips all three and is again transitive.

#### Unique unordered-Pythagorean collision in a branch triangle

Put `t=n/m in (0,1)`. The three reduced edge labels in the `A,B,C` triangles
are respectively

```text
A: t, 1/(2-t), t/(2-t),
B: t, 1/(2+t), t/(2+t),
C: t, t/(1+2t), 1/(1+2t).                                (CW38)
```

By THM-3509, two reduced fractions in `(0,1)` have the same primitive
normalized unordered Pythagorean triple exactly when they are related by

```text
tau(x)=(1-x)/(1+x).                                         (CW39)
```

Solving the three pair equations in each row of (CW38) gives

```text
A: t^2-4t+1=0,  2t-1=0,  t^2-3t+1=0;
B: t^2+2t-1=0, t^2+t-1=0, no solution;
C: 3t^2-1=0,   t=0,       t=0.                            (CW40)
```

The only rational root in `(0,1)` is `t=1/2`. Primitivity forces the root
`(m,n)=(2,1)`, branch `A`, with labels

```text
(1/2,2/3,1/3),
U(1/2)=U(1/3)={3,4,5},       U(2/3)={5,12,13}.             (CW41)
```

Reflecting only the closure `1->3` changes the natural transitive triangle
into

```text
1->2->3->1                                                   (CW42)
```

while preserving every gcd and the entire unordered-Pythagorean multiset.
Equation (CW41) is the unique additional shadow collision among all marked
branch incidences. Direction was already lost by reciprocal folding; the
collision shows that even distinct edge identities can merge at the root.

### 8. Validity boundary and generated frontiers

- The natural-order completion is transitive. Nontrivial tournament structure
  enters only through explicit selector bits; coprimality is not an intrinsic
  total orientation.
- Stern--Brocot reflection (CW9) and Calkin--Wilf ordinal (CW15) are compatible
  reciprocal mirrors but different address systems. Formula (CW22) belongs to
  Calkin--Wilf.
- The quotient to an unordered primitive triple loses reciprocal direction,
  leg order, and raw scale. Retaining ordered legs repairs the `tau` collision;
  retaining scale repairs the infinite dilation fibre.
- No LRC(14), planar Jacobian, `e+pi`, or Duffin--Schaeffer consequence follows
  from the shared Vieta grammar alone.
- Live generated tasks are: classify finite-state dilation-invariant
  selectors on Stern--Brocot words; combine branch closure signs with
  THM-3357's parity/Walsh characters; determine distinct rather than marked
  branch-triangle overlaps; and find the smallest extension of THM-3756's
  ordinal carrying both reciprocal direction and the `tau`-parity bit.

### 9. Replay

From the repository root:

```text
python -B 04-computation/calkin_wilf_berggren_ordinal_thm4057.py
python -B -O 04-computation/calkin_wilf_berggren_ordinal_thm4057.py
python -B 04-computation/calkin_wilf_berggren_ordinal_thm4057_independent_audit.py
python -B -O 04-computation/calkin_wilf_berggren_ordinal_thm4057_independent_audit.py
```

The primary path checks `1,048,575` Calkin--Wilf nodes and every Berggren node
through depth eleven. The independent path reconstructs word reflection,
all `1,024` tournaments on `[5]` from scale selectors, and the branch/gcd/
collision atlas through depth eight. Both normal/optimized pairs are
byte-identical. **QED.**
