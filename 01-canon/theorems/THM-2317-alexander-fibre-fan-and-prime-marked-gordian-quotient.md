---
id: THM-2317
title: "Alexander fibre fan and the prime-marked Gordian quotient"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE UNDER INDEPENDENT AUDIT.
  Maximal-ideal Alexander-fibre dimensions of a finite knot packet form a
  finite integral slope set whose support function is a polyhedral lower
  envelope for ordinary and homogenized unknotting number on the physical
  cone. A direct Fox calculation gives
  A_(T(2,7))=Lambda/(Delta), and the mirror has the same module, so the
  complete two-knot fibre fan has only slopes (0,0) and (1,1): its exact
  lower envelope is P+Q and cannot sharpen THM-2308. The factor-marked
  Gordian graph is a Cartesian product with l1 distance; connected-sum
  realization is a symmetric path-functor and a nonexpansive metric
  quotient. Its marked move corollas compose by fibre product, whereas
  support labels, binary relations, and shortest-path scalarization do not.
  Brittenham--Hermiller nonadditivity is exactly strict quotient contraction,
  and THM-2308's stable mixture defect is its asymptotic two-colour
  distortion. No exact unknotting number or positive catalyst is obtained.
source: codex-2026-07-25-stable-knot-relation
depends_on:
  - THM-2191-catalytic-localization-of-the-gordian-metric
  - THM-2308-mirror-double-nakanishi-floor-and-sharp-stable-mixture-profile
related:
  - THM-2176-gordian-continuation-profile-and-interaction-cocycle
  - THM-2248-higher-interaction-defect-complex-and-tropical-trace-spectrum
  - THM-2294-anchored-plucker-tournament-and-kakeya-address-bank
  - THM-2315-marked-target-gain-corolla-and-pairwise-composition-boundary
external:
  - "Mark Brittenham and Susan Hermiller, Unknotting number is not additive under connected sum, arXiv:2506.24088v2."
---

# THM-2317 -- the scalar is a quotient of a marked path object

**RESERVED / PROVISIONAL PROOF CANDIDATE UNDER INDEPENDENT AUDIT.**

Unknotting number is subadditive and is now known not to be additive under
connected sum. There are two different information losses behind that scalar
statement:

```text
lower certificates:
  one additive Alexander fibre -> maximum over fibres;

upper witnesses:
  factor-marked crossing paths -> forget the factor lift -> shortest path.
```

The first loss produces a finite tropical fan. The second is a metric
quotient. For the first Brittenham--Hermiller pair, the entire fibre fan
collapses to one ray, while the path quotient is already strict. This
separates the lower-bound obstruction from the actual geometric shortcut.

## 1. Every finite knot packet has an Alexander fibre fan

Put

```text
Lambda=Z[t,t^(-1)].
```

For a knot `J`, let `A_J` be its Alexander module. For a maximal ideal
`q` of `Lambda`, write

```text
k_q=Lambda/q,
nu_q(J)=dim_(k_q)(A_J tensor_Lambda k_q).             (1)
```

Fix a finite labelled packet of knots

```text
K_1,...,K_s.
```

Its **Alexander fibre-slope set** is

```text
V_A={v(q):q maximal in Lambda},
v(q)=(nu_q(K_1),...,nu_q(K_s)) in N^s.               (2)
```

This set is finite. Indeed, THM-2308 gives

```text
0<=nu_q(K_i)<=m_N(K_i)<=u(K_i),                      (3)
```

so `V_A` lies in a finite integer box. Let

```text
P_A=conv(V_A),
alpha_A(x)=max_(v in V_A) <v,x>
          =max_(v in P_A) <v,x>,    x in R_(>=0)^s.  (4)
```

There is also a signed group-completion envelope:

```text
P_A^+-=conv(V_A union -V_A),
beta_A(z)=max_(v in V_A)|<v,z>|
         =max_(v in P_A^+-)<v,z>,       z in R^s.    (4a)
```

Thus `alpha_A` is the support function of an integral polytope and `beta_A`
is the support function of its central symmetrization. They are continuous,
convex, positively homogeneous, subadditive, and piecewise integral-linear.
Their domains of linearity form finite normal fans. Since every slope has
nonnegative coordinates,

```text
beta_A(x)=alpha_A(x)               for x in R_(>=0)^s. (4b)
```

The reason `beta_A` controls formal differences is a local crossing-change
lemma.

> **Fibre Lipschitz lemma.** For every maximal ideal `q`,
>
> ```text
> |nu_q(J)-nu_q(L)|<=d_G(J,L).                        (4c)
> ```
>
> Consequently `nu_q` extends from the connected-sum monoid to an additive
> integer-valued Gordian-`1`-Lipschitz functional on its Grothendieck group.

To prove the one-crossing case, use the standard band-twist form of a
crossing change. Simultaneous Seifert surfaces for `J` and `L` have the same
genus and, in suitable homology bases, Seifert matrices `V_J,V_L` with

```text
V_L-V_J=+E_(rr)             or             -E_(rr). (4d)
```

One way to see (4d) is to apply the Seifert disk-band construction
simultaneously, take the core of the changed band as the last basis element,
and observe that reversing the crossing changes only that band's
self-linking by one; all mutual linkings are unchanged.

The square Alexander presentation matrices `tV-V^T` therefore differ by

```text
+/-(t-1)E_(rr),                                      (4e)
```

a matrix of rank at most one after reduction to `k_q`. Over a field, the
nullities, equivalently the cokernel dimensions, of two same-size matrices
whose difference has rank at most one differ by at most one. This proves the
one-crossing case. A Gordian path proves (4c). Direct-sum additivity proves
the last assertion.

For `n=(n_1,...,n_s) in N^s`, put

```text
K(n)=#_(i=1)^s (#^(n_i) K_i).                        (5)
```

> **Fibre-fan lower envelope.** For every `n in N^s`,
>
> ```text
> u(K(n))>=alpha_A(n),
> u_hash(K(n))>=alpha_A(n).                          (6)
> ```
>
> More generally, on the rationalized group span,
>
> ```text
> u_hash(sum_i z_i[K_i])>=beta_A(z).                 (6a)
> ```
>
> The stable inequalities extend by homogeneity and continuity to the
> corresponding real span.

### Proof

Connected sum gives a direct sum of Alexander modules. Hence for every
fixed maximal ideal,

```text
nu_q(K(n))=sum_i n_i nu_q(K_i)=<v(q),n>.             (7)
```

Equation (3), applied to `K(n)`, gives

```text
u(K(n))>=nu_q(K(n)).
```

Maximizing over `q` proves the first inequality in (6).

For any positive integer `r`,

```text
#^r K(n)=K(rn),
alpha_A(rn)=r alpha_A(n).
```

Apply the first inequality to `rn`, divide by `r`, and let `r` tend to
infinity. This proves the stable inequality. THM-2191 extends `u_hash` to a
seminorm `p` on the rationalized knot group. The fibre Lipschitz lemma makes
each additive extension of `nu_q` one of THM-2191's admissible stable dual
functionals. Therefore

```text
p(sum_i z_i[K_i])
 >=max_q |sum_i z_i nu_q(K_i)|
 =beta_A(z).                                         (6b)
```

first for rational `z`. Positive homogeneity and finite-dimensional
continuity give the real extension. Equation (4b) recovers the physical-cone
bound. QED.

Each individual slope `v(q)` is additive under connected sum. The maximum
over `q` is not additive in general because the maximizing fibre may switch.
The fan records exactly that selector change. It is a concrete
operation-compatible lower envelope, not a complete Alexander invariant and
not an equality claim for unknotting number.

## 2. The `T(2,7)` fibre fan is exactly one nonzero ray

Let

```text
K=T(2,7),
Delta=t^6-t^5+t^4-t^3+t^2-t+1.                     (8)
```

The cyclicity of `A_K` can be proved directly, rather than inferred from its
Alexander polynomial. Use the standard torus-knot presentation

```text
pi_1(S^3\K)=<x,y | x^2 y^(-7)>                      (9)
```

and the abelianization

```text
x -> t^7,                 y -> t^2.
```

Fox differentiation gives the Alexander chain maps

```text
Lambda --d_2--> Lambda^2 --d_1--> Lambda,

d_1=(t^7-1, t^2-1),
d_2=(1+t^7, -sum_(j=0)^6 t^(2j))^T.                (10)
```

Set

```text
a=1+t+...+t^6,
b=1+t,
Delta=t^6-t^5+t^4-t^3+t^2-t+1.                     (11)
```

Then

```text
d_1=(t-1)(a,b),
d_2=Delta (b,-a)^T,                                 (12)
```

because

```text
b Delta=1+t^7,
a Delta=1+t^2+...+t^12.
```

The two polynomials `a,b` generate the unit ideal: reducing `a` modulo
`b=t+1` gives `a(-1)=1`. Since `Lambda` is a domain,

```text
ker d_1=Lambda (b,-a)^T.                             (13)
```

Equations (12)--(13) prove

```text
A_K=H_1(cover_tilde(K);Z)
   isomorphic to Lambda/(Delta).                     (14)
```

The mirror twists the `Lambda` action by `t -> t^(-1)`. Since

```text
Delta(t^(-1))=t^(-6)Delta(t),
```

the ideal `(Delta)` is invariant up to a unit, and

```text
A_(Kbar) isomorphic to Lambda/(Delta).               (15)
```

Consequently, for every maximal ideal `q`,

```text
nu_q(K)=nu_q(Kbar)
 =1 if Delta in q,
 =0 otherwise.                                      (16)
```

Both cases occur. THM-2308 gives the explicit maximal ideal

```text
q_1=(2,t^3+t+1)
```

containing `Delta`, while

```text
q_0=(2,t+1)
```

does not contain `Delta`, because `Delta(1)=1` in `F_2`. Therefore

```text
V_A(K,Kbar)={(0,0),(1,1)},
P_A=conv{(0,0),(1,1)},                               (17)

alpha_A(P,Q)=P+Q                  for P,Q>=0,
beta_A(P,Q)=|P+Q|                 for P,Q in R.      (18)
```

Equation (18) proves a sharp stopping statement for this method:

> **Fibre-bank saturation.** The complete bank of maximal-ideal Alexander
> module dimensions for `T(2,7)` and its mirror gives exactly THM-2308's
> common-fibre floor `P+Q`. Switching maximal ideals, taking their convex
> hull, passing to signed formal differences, or using the whole normal fan
> cannot improve that physical-cone floor or distinguish the two chiralities.
> On the group plane it sees only the total coordinate `P+Q` and is blind to
> the anti-mirror coordinate `P-Q`.

This does not rule out stronger information in the full Blanchfield form,
higher-order Alexander data, concordance invariants, or actual crossing
witnesses. It says that the dimension-fibre quotient itself is exhausted.

There is a useful exact two-coordinate synthesis. THM-2308's half-signature
functional has values

```text
S(K)=3,                     S(Kbar)=-3.
```

Together with any nonzero fibre functional from (16), it gives the explicit
recognizable stable lower seminorm

```text
gamma_A,S(P,Q)=max(|P+Q|,3|P-Q|)
 <=u_hash(P[K]+Q[Kbar]).                             (18a)
```

This is the support function of

```text
conv{+/-(1,1), +/-(3,-3)}.
```

THM-2308 proves equality with the second term throughout the
opposite-sign chambers `PQ<=0`. Thus the two classical slopes already
reconstruct those chambers exactly. Their common kernel geometry also
pinpoints the residual: the Alexander slope is blind to `P-Q`, the signature
slope is blind to `P+Q`, and neither determines the same-sign interior above
the diagonal floor.

## 3. Factor marking restores an exact Cartesian product

Let `Gamma_G` be the undirected Gordian graph: its vertices are oriented
knot types and two vertices are adjacent when one crossing change relates
them. Its path metric is `d_G`.

For a finite label set `I`, form the Cartesian product graph

```text
Gamma_I=box_(i in I) Gamma_G.                        (19)
```

A vertex is a labelled tuple `x=(x_i)_(i in I)`, and one edge changes
exactly one coordinate by one Gordian edge. The label is a sidecar: when
the `x_i` are prime factors it refines Schubert's unordered prime
decomposition by retaining the factor names.

> **Marked product theorem.** For all `x,y in Gamma_I`,
>
> ```text
> d_I(x,y)=sum_(i in I)d_G(x_i,y_i).                 (20)
> ```
>
> The connected-sum realization
>
> ```text
> Sigma_I:Gamma_I -> Gamma_G,
> Sigma_I(x)=#_(i in I)x_i                           (21)
> ```
>
> is a graph map, is symmetric under relabelling, and is monoidal under
> disjoint union of label sets. In particular,
>
> ```text
> d_G(Sigma_I(x),Sigma_I(y))<=d_I(x,y).              (22)
> ```

### Proof

Every path in the Cartesian product changes one coordinate per step. Its
number of `i`-steps is at least `d_G(x_i,y_i)`, so every path has length at
least the right side of (20). Concatenating one coordinate geodesic for
each `i` attains this sum.

If an edge changes the `i`th knot by one crossing change, perform the same
crossing change inside the `i`th connected-summand ball. This gives one
Gordian edge between the connected sums, so `Sigma_I` is a graph map and
(22) follows. Relabelling summands changes the connected sum only by
isotopy, and

```text
Sigma_(I disjoint_union J)(x,y)
 =Sigma_I(x)#Sigma_J(y),
```

which proves the symmetry and monoidal assertions. QED.

At the root `0_I=(U,...,U)`, define the **quotient contraction defect**

```text
delta_I(x)
 =d_I(x,0_I)-d_G(Sigma_I(x),U)
 =sum_i u(x_i)-u(#_i x_i)>=0.                        (23)
```

Thus connected-sum nonadditivity is not a failure of the marked local
geometry. The marked geometry is exactly `l_1`; nonadditivity is the
contraction created by forgetting the factor lift and then minimizing in
the larger unmarked Gordian graph.

Equivalently, strict defect means that no unmarked geodesic from
`Sigma_I(x)` to the unknot lifts to a coordinatewise path from `x` to
`0_I` of the same length. This is the precise path-lifting obstruction.

## 4. The marked crossing-move corolla composes functorially

For `x in Gamma_I`, let `E_x` be its set of oriented outgoing marked edges.
An element retains

```text
(changed coordinate i, actual Gordian edge e, target tuple y).           (24)
```

There are source, target, and coordinate maps

```text
V(Gamma_I) <-s- E_I -t-> V(Gamma_I),
lambda:E_I->I.                                       (25)
```

This is the **marked crossing-move corolla**. Length-`r` marked path
witnesses are the iterated fibre product

```text
E_I x_(t,s) E_I x_(t,s) ... x_(t,s) E_I.            (26)
```

The canonical bijection between the two parenthesizations of (26) is the
associativity of path composition. The graph map `Sigma_I` sends every
marked edge to the corresponding local crossing edge and therefore induces
a functor on path categories. It preserves concatenation and the length of
each chosen marked path.

Forgetting `(e,y)` and retaining only the coordinate word

```text
lambda(e_1)...lambda(e_r)                            (27)
```

does not preserve the terminal knot. The collision is present already in
one step: in one coordinate, the unknot has Gordian edges to both a trefoil
and the figure-eight knot, while the word in (27) is the same one-letter
word. A binary relation on factor labels forgets still more.

Taking shortest-path length is a second, different quotient. Although every
chosen marked path maps functorially, an unmarked geodesic need not lie in
the image of any marked geodesic. Hence

```text
path witnesses compose exactly;
their coordinate support words do not determine targets;
terminal minimization after forgetting the lift can strictly contract.  (28)
```

This is the knot analogue of THM-2315's span boundary, but the proof here
uses the actual Cartesian Gordian product and does not depend on the
provisional LRC theorem.

## 5. Brittenham--Hermiller is strict quotient contraction

Return to

```text
K=T(2,7),                  Kbar=mirror(K).
```

The one-body costs are

```text
u(K)=u(Kbar)=3.
```

For the marked two-factor point `x=(K,Kbar)`,

```text
d_{ {1,2} }(x,(U,U))=6.                              (29)
```

Brittenham--Hermiller give

```text
u(K#Kbar)<=5,                                        (30)
```

so `Sigma_{ {1,2} }` is strictly distance-contracting on this pair. The
shortcut is exactly a nonliftable unmarked path, not a failure of the
coordinate product law.

For integers `P,Q>=0`, let `x_(P,Q)` be the labelled tuple of `P` copies of
`K` and `Q` copies of `Kbar`. Then

```text
delta_(P,Q)
 =3(P+Q)-u((#^P K)#(#^Q Kbar)).                     (31)
```

Replicate the entire labelled packet `n` times. THM-2308's stable mixture
defect satisfies the exact quotient-distortion identity

```text
D(P,Q)
 =lim_(n->infinity) delta_(nP,nQ)/n
 =3(P+Q)-u_hash((#^P K)#(#^Q Kbar)).                (32)
```

Thus its one-variable concave profile is the asymptotic distortion profile
of the two-colour connected-sum realization. Section 2 supplies the exact
Alexander-fan floor

```text
u_hash((#^P K)#(#^Q Kbar))>=P+Q,                    (33)
```

but the fan has no further slope with which to constrain the interior of
that profile. The remaining information must enter before the fibre
dimension quotient or after it through a different marked witness.

## 6. Why this is not a tournament

The scalar two-factor contraction

```text
w(J,L)=u(J)+u(L)-u(J#L)                              (34)
```

is symmetric. Its positive support is an undirected weighted symbiont
graph, not an orientation. For `K,Kbar`, mirror and factor swap preserve
the entire datum in (34). Any head selector using only this symmetric datum
and equivariant under the swap must return a tie; choosing a direction adds
an external gauge.

For larger packets, pair weights are only the two-dimensional faces of the
metric quotient. THM-2248 gives exact abstract metric-monoid controls in
which every pair is additive while a higher packet has strict defect, so no
pair graph, and therefore no tournament built from it, classifies quotient
contraction in the allowed axiomatic category.

The lawful relation is instead the span (25):

```text
move witness -> changed factor,
move witness -> actual target tuple.                 (35)
```

Its path fibre products retain simultaneity and composition. The
connected-sum map then deliberately forgets the lift. Directional catalytic
observables from THM-2176 can define different pairwise data, but they do
not turn the symmetric connected-sum defect into a tournament.

## 7. Relation and loss ledger

| source | map | preserved | first information lost | repair |
|---|---|---|---|---|
| Alexander modules | `A_J -> (nu_q(J))_q` | direct-sum fibre dimensions | extension, pairing, and chirality data | full module or pairing sidecar |
| fibre vectors | maximum over `q` | best common-fibre lower floor | which one fibre should persist under mixture | normal-fan selector |
| labelled factor tuple | `Sigma_I` | every chosen local crossing path | factor lift of an unmarked path | tuple plus marked edge |
| marked move corolla | coordinate word | move-support sequence | actual edge and target | the span `(s,t,lambda)` |
| quotient metric | pair defect graph | two-factor contractions | higher packet distortion | full subset/cost profile |
| symmetric pair weight | attempted orientation | magnitude only | swap symmetry | honest tie or a new directional observable |

The theorem produces three concrete stopping boundaries:

1. maximal-ideal fibre dimensions cannot sharpen the `T(2,7)` mirror floor
   beyond `P+Q`;
2. marking factors restores exact `l_1` composition but cannot force the
   unmarked shortest path to lift;
3. pairwise contraction is intrinsically symmetric and higher interaction
   can be invisible on every pair.

No exact value of `u(K#Kbar)` or `u_hash(K#Kbar)` follows, no positive
Gordian catalyst is produced, and no tournament classification of knots is
claimed.
