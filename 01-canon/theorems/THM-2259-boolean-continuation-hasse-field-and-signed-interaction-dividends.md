---
id: THM-2259
title: "Boolean continuation Hasse field and signed interaction dividends"
status: >
  PROVED (abstract commutative nonexpansive metric-monoid theorem) + CITED
  APPLICATION (Brittenham--Hermiller). On every labelled packet and at every
  fixed diagonal continuation, the multiway defect is a normalized
  superadditive Boolean game. Its increments form a nonnegative curl-free
  field on the Boolean Hasse diagram and reconstruct every subset defect
  along any ordering. Boolean Möbius inversion gives unique signed
  interaction dividends; all disjoint-cut sums are nonnegative, and a
  minimal nonface is exactly a pure positive dividend, but higher dividends
  can have either sign. Two explicit translation-invariant integer word
  metrics have identical singleton and pair response data at every context,
  and the same weighted pair defect graph and zero-defect complex, while
  their triple dividends are respectively -1 and +1. Thus even
  context-indexed weighted pair relations cannot classify compositional
  interaction. For T(2,7) and its mirror, the Brittenham--Hermiller edge is
  persistent at every knot context with defect at least one. No new
  unknotting number, higher-order knot symbiont, or positive Gordian catalyst
  is proved.
source: codex-2026-07-25-knot-interaction-cube
depends_on:
  - THM-2176-gordian-continuation-profile-and-interaction-cocycle
  - THM-2220-fixed-context-stable-response-and-catalyst-complexity
  - THM-2248-higher-interaction-defect-complex-and-tropical-trace-spectrum
related:
  - THM-2191-catalytic-localization-of-the-gordian-metric
  - THM-2221-tournament-context-cut-metric-and-pinned-transport-response
  - THM-2242-tournament-complement-transport-and-knot-kernel-green-rigidity
external:
  - "Mark Brittenham and Susan Hermiller, Unknotting number is not additive under connected sum, arXiv:2506.24088v2."
---

# THM-2259 -- the finite interaction carrier is a flat Hasse field

THM-2176 identifies the full min-plus continuation kernel and its positive
binary defect. THM-2220 proves that every fixed diagonal context gives a
subadditive response. THM-2248 turns the root defects of a labelled packet
into a zero-defect simplicial complex and shows that its minimal nonfaces can
have arbitrary arity.

The complex records where interaction first appears, but not how already
active shortcuts reinforce or screen one another. The missing finite
coordinate is the signed Boolean Möbius spectrum. Equivalently, one may retain
the nonnegative integrable edge field on the Hasse diagram of partial sums.

## 1. Diagonal continuation games

Let `(M,+,0)` be a commutative monoid with a metric `d` satisfying joint
nonexpansivity

```text
d(a+c,b+d) <= d(a,b)+d(c,d).                         (1)
```

Fix a labelled packet

```text
x_1,...,x_r in M.
```

Labels remain distinct even if two packet entries are equal. For
`S subset [r]`, put

```text
x_S=sum_(i in S)x_i,             x_emptyset=0.       (2)
```

For each fixed diagonal context `z in M`, define

```text
rho_z(x)=d(x+z,z),                                    (3)

Delta_z(S)=sum_(i in S)rho_z(x_i)-rho_z(x_S).         (4)
```

THM-2220's fixed-context argument gives

```text
rho_z(x+y)<=rho_z(x)+rho_z(y).                        (5)
```

Indeed, insert `y+z` and translate the first leg:

```text
d(x+y+z,z)
 <=d(x+y+z,y+z)+d(y+z,z)
 <=d(x+z,z)+d(y+z,z).                                (6)
```

Consequently

```text
Delta_z(S)>=0,
Delta_z(emptyset)=Delta_z({i})=0.                    (7)
```

More precisely, for disjoint `A,B`,

```text
Delta_z(A union B)-Delta_z(A)-Delta_z(B)
 =rho_z(x_A)+rho_z(x_B)-rho_z(x_A+x_B)
 =:sigma_z(A,B)>=0.                                  (8)
```

Thus `Delta_z` is a normalized superadditive set function on disjoint
coalitions. Taking `B` one vertex at a time also makes it monotone under
inclusion. Therefore

```text
K_z={S subset [r]:Delta_z(S)=0}                      (9)
```

is an abstract simplicial complex for every context `z`. At `z=0`, this is
exactly THM-2248's root interaction complex.

There are two canonical context envelopes:

```text
K_all = intersection_(z in M) K_z,
K_some= union_(z in M) K_z.                         (10)
```

Both are simplicial complexes, because arbitrary intersections and unions of
downward-closed families are downward closed. They separate three behaviors:

```text
S in K_all:                 additive at every diagonal context;
S in K_some\K_all:          context-sensitive;
S not in K_some:            positive defect at every context.   (11)
```

No finiteness or attainment claim is hidden in (10).

## 2. The nonnegative Boolean Hasse field

For `i notin A`, define the context-edge weight

```text
c_z(i|A)
 =Delta_z(A union {i})-Delta_z(A)
 =rho_z(x_A)+rho_z(x_i)-rho_z(x_A+x_i)
 =sigma_z(A,{i}).                                    (12)
```

Equation (5) proves

```text
c_z(i|A)>=0.                                         (13)
```

These weights live on the intrinsically oriented Hasse edges

```text
A -> A union {i}
```

of the Boolean lattice. They are not orientations of pairs of atomic
vertices. Every square is flat: if `i,j notin A` and `i!=j`, then

```text
c_z(i|A)+c_z(j|A union {i})
 =Delta_z(A union {i,j})-Delta_z(A)
 =c_z(j|A)+c_z(i|A union {j}).                       (14)
```

Hence the field is curl-free. Conversely, its potential is recovered by
path integration. For every ordering `i_1,...,i_m` of `S`, with
`S_k={i_1,...,i_k}`,

```text
Delta_z(S)=sum_(k=1)^m c_z(i_k|S_(k-1)).             (15)
```

The first term vanishes because `c_z(i|emptyset)=0`. Flatness makes the sum
independent of the ordering. Thus either of the following is an exact finite
carrier for all pinned subset defects at the chosen context:

```text
the potential table Delta_z(S);
the nonnegative flat Hasse field c_z(i|A).            (16)
```

The second representation explains the recursion: vertices are partial
composites, not only the original objects.

### Continuation-kernel formula

Let

```text
P_x(a,b)=d(x+a,b)
```

be THM-2176's faithful min-plus kernel. With tropical convolution `tensor`,

```text
P_(x+y)=P_x tensor P_y.                              (17)
```

Write

```text
P_S=tensor_(i in S)P_(x_i)=P_(x_S).
```

Then

```text
rho_z(x_S)=P_S(z,z),                                 (18)

c_z(i|A)
 =P_A(z,z)+P_(x_i)(z,z)
  -(P_A tensor P_(x_i))(z,z).                       (19)
```

So the Hasse field is the defect in composing diagonal kernel responses. It
retains every diagonal subset response in the chosen packet, but still
forgets all off-diagonal entries `P_S(a,b)` with `a!=b`. The full kernel
remains the faithful operation carrier.

## 3. Boolean Möbius/Harsanyi interaction dividends

Define the signed interaction dividend of a nonempty face `A` by Boolean
Möbius inversion:

```text
h_z(A)
 =sum_(B subset A)(-1)^(|A|-|B|)Delta_z(B).          (20)
```

Here and below `subset` in a summation permits equality. The inverse formula
is

```text
Delta_z(S)=sum_(A subset S)h_z(A).                   (21)
```

Since `Delta_z(emptyset)=Delta_z({i})=0`, the empty and singleton dividends
vanish. For `|A|>=2`, the modular singleton term in (4) has zero Möbius
transform, and hence

```text
h_z(A)
 =-sum_(B subset A)(-1)^(|A|-|B|)rho_z(x_B)

 =-sum_(B subset A)(-1)^(|A|-|B|)P_B(z,z).          (22)
```

Thus `h_z` is the unique signed inclusion-exclusion spectrum of the pinned
diagonal continuation kernels.

Subtracting (21) on the two endpoints of a Hasse edge gives

```text
c_z(i|A)=sum_(B subset A)h_z(B union {i}).           (23)
```

In particular, individual higher dividends need not be nonnegative. What
metric subadditivity forces is the nonnegativity of every cumulative edge
sum in (23). More generally, for disjoint nonempty `A,B`, equations
(8) and (21) give the exact Boolean cut inequality

```text
sigma_z(A,B)
 =sum_(
      T subset A union B,
      T intersection A != emptyset,
      T intersection B != emptyset
    ) h_z(T)
 >=0.                                               (24)
```

The signed spectrum is therefore constrained by nonnegative sums across
every disjoint coalition cut, not by coefficientwise positivity.

### Minimal nonfaces are pure positive dividends

Let `A` be a minimal nonface of `K_z`. Then

```text
Delta_z(A)>0,
Delta_z(B)=0                  for every proper B subset A.   (25)
```

Equation (20) immediately yields

```text
h_z(B)=0                     for every proper B subset A,
h_z(A)=Delta_z(A)>0.                                  (26)
```

Conversely, if every proper dividend on `A` vanishes and `h_z(A)>0`, then
(21) gives (25). Thus:

> Minimal nonfaces are exactly the pure positive interaction dividends.

Negative dividends can occur only above interaction that has already
appeared on proper faces. They measure screening or redundancy among lower
shortcuts; a positive nonminimal dividend measures additional higher-order
reinforcement.

For three labels the formula is especially transparent:

```text
h_z(123)
 =Delta_z(123)-Delta_z(12)-Delta_z(13)-Delta_z(23)

 =c_z(1|23)-c_z(1|2)-c_z(1|3).                     (27)
```

It is the failure of the response to label `1` to split additively across
the two-label context.

## 4. Exact pair-shadow collision: screening versus reinforcement

Higher dividends genuinely require signs. They are not the weights of a
nonnegative hypergraph, and even the complete weighted pair graph does not
determine them.

Let `M=Z^3` with standard basis `e_1,e_2,e_3`, and put

```text
p=e_1+e_2,
q=e_2+e_3,
g=e_1+e_2+e_3.                                      (28)
```

Define two symmetric weighted word metrics.

For `d_-`, give the generators

```text
+/-e_1,+/-e_2,+/-e_3                 cost 2,
+/-p,+/-q                             cost 3.        (29)
```

For `d_+`, use the same generators and add

```text
+/-g                                  cost 3.        (30)
```

Let `ell_-(v)=d_-(v,0)` and `ell_+(v)=d_+(v,0)`. Both are
translation-invariant integer metrics, and the triangle inequality for word
length gives joint nonexpansivity (1).

For the labelled packet `(e_1,e_2,e_3)`, the complete relevant length table
is

| `S` | singleton | `12` | `13` | `23` | `123` |
|---|---:|---:|---:|---:|---:|
| `ell_-(e_S)` | `2` | `3` | `4` | `3` | `5` |
| `ell_+(e_S)` | `2` | `3` | `4` | `3` | `3` |

Here is a direct proof of the table. The displayed generators give every
upper bound. A singleton has length two because the least generator cost is
two. The vectors `p,q` have length three: neither is a cost-two generator,
and two letters cost at least four. The vector `e_1+e_3` has length four:
it is not any generator, while `e_1+e_3` is a two-letter word of cost four.
In `d_-`, the vector `g` has the cost-five word `p+e_3`; no one-letter word
is `g`, and the only words of cost below five after the one-letter cases are
two cost-two coordinate letters, which cannot sum to a vector with three
nonzero coordinates. In `d_+`, `g` is a cost-three generator. Adding it does
not lower any earlier entry: none of those vectors is the new generator, so
the same lower bounds below four remain valid.

The two defect tables therefore agree on every singleton and pair:

```text
Delta_-(12)=Delta_+(12)=1,
Delta_-(13)=Delta_+(13)=0,
Delta_-(23)=Delta_+(23)=1.                           (31)
```

But

```text
Delta_-(123)=1,              h_-(123)=-1,
Delta_+(123)=3,              h_+(123)=+1.            (32)
```

Both zero-defect complexes are literally

```text
{emptyset,{1},{2},{3},{1,3}},                        (33)
```

and both weighted minimal-nonface ledgers consist of the two edges
`12,23`, each of weight one. Nevertheless their triple defects and the signs
of their triple dividends differ.

The Hasse field locates the lost coordinate:

```text
c_-(1|{2})=c_+(1|{2})=1,

c_-(1|{2,3})=0,             c_+(1|{2,3})=2.          (34)
```

In `d_-`, the shortcut through `p` becomes redundant once the shortcut
through `q` is present. In `d_+`, the new generator `g` creates additional
reinforcement.

Because both metrics are translation invariant,

```text
rho_z(e_S)=ell(e_S)             for every z in Z^3.  (35)
```

Thus the two examples have identical singleton and pair diagonal response
data at **every** context, not merely at the root. Any context-indexed family
of pair graphs, weighted binary symbiosis relations, or orientations
constructed from those pair data gives the same answer on the two models,
while the triple response differs. This is an exact stopping boundary for
pairwise continuation shadows.

## 5. Knot application: a persistent Brittenham--Hermiller edge

Now take `M` to be oriented knot types under connected sum and `d=d_G`.
For a context knot `J`,

```text
rho_J(K)=d_G(K#J,J).                                  (36)
```

On a labelled packet of prime-summand occurrences, the tables `Delta_J`,
`c_J`, and `h_J` are canonical up to relabelling by Schubert prime
decomposition. They refine THM-2248's root interaction complex to a stalk
over all diagonal continuations.

Let

```text
K=T(2,7),                  Kbar=mirror(K).            (37)
```

Brittenham--Hermiller prove

```text
u(K)=u(Kbar)=3,
u(K#Kbar)<=5.                                       (38)
```

For every knot `J`, half the signature difference gives the lower bound

```text
d_G(K#J,J)>=3,
d_G(Kbar#J,J)>=3.                                   (39)
```

Translation of a three-change unknotting path gives the reverse inequalities.
Hence

```text
rho_J(K)=rho_J(Kbar)=3                  for every J. (40)
```

Likewise, translating the published at-most-five crossing-change path gives

```text
rho_J(K#Kbar)
 =d_G(K#Kbar#J,J)
 <=u(K#Kbar)
 <=5.                                               (41)
```

Therefore the pair dividend is uniformly positive:

```text
h_J({K,Kbar})
 =Delta_J({K,Kbar})
 =6-d_G(K#Kbar#J,J)
 >=1                         for every knot J.       (42)
```

Thus `{K,Kbar}` lies outside `K_some`: it is a persistent interaction edge,
not merely a defect at the unknot context. This conclusion uses only the
published upper-bound certificate and the classical signature calibration.
It does **not** determine the exact value in (42), even at `J=U`.

## 6. What binary relations and tournaments lose

For a packet of atomic objects, the pair dividend

```text
h_z({i,j})=sigma_z({i},{j})                         (43)
```

is symmetric, nonnegative, weighted, and legitimately zero. There is no
intrinsic tournament orientation in this observable. Forcing one discards
ties and still leaves all higher coefficients absent.

The exact hierarchy is:

```text
pair support:              which two-body shortcuts occur;
weighted pair graph:       their magnitudes;
zero-defect complex:       first activation faces of every arity;
signed dividends h_z:      reinforcement and screening by inclusion-exclusion;
flat Hasse field c_z:      nonnegative recursive response on partial sums;
full kernel P_x(a,b):      exact continuation through arbitrary intermediates.
                                                               (44)
```

The pair of word metrics in Section 4 agrees on the first three rows of
(44), and on every context-indexed pair response, while disagreeing on the
fourth and fifth. The signed dividend table and flat Hasse field are
equivalent finite encodings of all subset defects through (21) and (23).
Neither recovers the off-diagonal kernel.

The natural finite directed relation is therefore the Boolean Hasse DAG

```text
A -> A union {i}
```

with weight `c_z(i|A)`, not a tournament on atomic summands. Its direction
comes from inclusion, its weights are nonnegative, and its square relations
(14) are the associativity audit. The Möbius spectrum is the unique signed
coordinate change which exposes pure higher interactions and lower-order
screening.

## 7. Scope and failure boundary

1. The context-indexed construction is abstract and applies to every
   commutative nonexpansive metric monoid. The two sign examples are weighted
   word metrics, not Gordian metrics.
2. Equation (42) proves a persistent pair edge, not its exact weight. It
   computes no previously unknown unknotting number.
3. No higher-order knot minimal nonface, negative knot dividend, or positive
   Gordian translation catalyst is exhibited.
4. Infimizing over contexts does not commute with the signed Möbius
   transform. The diagonal stalk must not be replaced by separately
   minimized catalytic scalars.
5. Even the complete diagonal stalk forgets off-diagonal common
   intermediates. THM-2176's full min-plus kernel remains the faithful
   continuation object.

The structural reframe is:

```text
atomic binary relation
  -> partial-composite Hasse field
  -> signed Möbius interaction spectrum
  -> full min-plus continuation kernel.             (45)
```

QED.
