# Two sentences, one schema, and the missing witness transport

**PROVED elementary reformulations and coding lemmas / CITED foundational
equivalence and reflection / FINITE-EXACT controls / OPEN Collatz.**
Date: 2026-09-26. This is a research sidecar, not a new independence result
or a claim of novelty for Ackermann coding or finite set theory.

The proposed small axiom package has a precise realization for arithmetic.
There is also a reflection presentation of ZFC that can be packaged in the
same syntactic shape. The substantive distinction is which induction or
reflection schema is available, not the number of primitive symbols or
displayed sentences. Transitive containment, even with closure under taking
subsets of members, is already compatible with the finite-set world.

## 1. Inheritance and the operative connection

Closest proved mechanism: [the guarded affine ports](creative_descent_20260925.md),
which keep both the source residue and the actual endpoint. Canonical hostile:
an encoded Collatz step need not decrease the well-founded rank of its code.
The example `3 -> 5` below gives a two-state witness. Corrected near miss:
the [creation decoder](creation_decoder_20260925.md) transports a legitimate
elliptic reconstruction, but its height need not decrease under the distinct
Collatz operation. Least-used sidecar: preservation of an **existential
witness**, not just preservation of the objects mentioned in a formula.

The concept board is finite syntax, unbounded depth, exact arithmetic coding,
structural induction, witness reflection, and selected-return certificates.
The anchor is an exact certificate interface; the niche is finite set theory;
the wildcard is reflection as a disciplined comparison of a large structure
with a smaller one. Targeted searches found no prior PA/HF theorem in the
current results/canon routes. The recovered
[THM-4090 / two-sort matching-logic obstruction](../../01-canon/theorems/THM-4090-two-sort-matching-logic-global-completeness-obstruction.md)
is a useful methodological warning about proof transport, but its particular
calculus obstruction is not an incompleteness theorem about PA or Collatz.

## 2. An exact arithmetic package in the language of one binary relation

Use ordinary first-order logic with equality and the sole nonlogical symbol
`in`. The following are abbreviations, not new primitives:

```text
Empty(x)     := forall u, not (u in x).
Adj(x,y,z)   := forall u, (u in z iff (u in x or u=y)).
```

Take two fixed sentences:

```text
(E)  forall x,y, [(forall u, u in x iff u in y) -> x=y].
(A)  forall x,y, exists z, Adj(x,y,z).
```

Take one schema, for every formula `phi(v,p)` with arbitrary parameters `p`
(universally close the displayed free parameters):

```text
(I_phi)
 [forall x, (Empty(x) -> phi(x,p))
  and forall x,y,z,
      (Adj(x,y,z) and phi(x,p) and phi(y,p) -> phi(z,p))]
 -> forall x, phi(x,p).
```

**PROVED reformulation.** This theory is definitionally equivalent to the
usual structural-induction theory HF. Here is the entire reduction.

First, instantiate the schema with a contradiction, such as `v != v`.
If there were no empty object, both premises would hold and the conclusion
would contradict the nonempty domain of first-order logic. Thus an empty
object exists. Sentence E makes it unique; A and E make every adjunct unique.
Introduce their names `0` and `x triangle y` by definitions. Then the scheme
becomes induction from `0` under adjoining one object to another, with both
arguments satisfying the induction predicate.

Conversely, in the expanded HF presentation use the two sentences

```text
z=0              iff Empty(z),
z=x triangle y   iff Adj(x,y,z),
```

and its structural induction. Induction proves that every object is either
`0` or some adjunct: the base is immediate and every induction step is an
adjunct by construction. Two objects with the same elements therefore
coincide: use the first displayed biconditional for an empty object, or the
second for its adjunct representation. This proves E; A follows from the
function's totality. Definitional elimination transports the induction
instances in both directions.

**CITED connection to PA.** S. Swierczkowski gives the two-axiom/one-schema
HF presentation in section 1 and its definitional equivalence with PA in
section 10 of [*Finite sets and Godel's incompleteness theorems* (2003)](https://www.impan.pl/shop/en/publication/transaction/download/product/88124).
The published signature includes `0` and adjunction as well as membership;
the elementary elimination above supplies the literal one-relation version.
Thus this is an affirmative, precise reading of the user's proposed package.

The schema is **structural adjunction induction**, not merely induction
along membership. ZF with Infinity also supports membership induction.
Moreover HF is a first-order theory: its nonstandard models should not be
identified with the externally standard universe of hereditarily finite sets.

## 3. What transitive containment does, and does not, add

Write

```text
Trans(t) := forall y in t, forall z in y, z in t.
TC       := forall x, exists t, [x subseteq t and Trans(t)].
```

In the structural-induction HF presentation, finite union and transitive
containers are provable by induction. For the inductive construction, if
`t_x` and `t_y` are transitive and contain all members of `x` and `y`, then
`t_x union t_y union {y}` is a transitive container for `x union {y}`.
The required finite-union operation is itself constructed by adjoining the
members one at a time under structural induction. TC is consequently not
an additional unmotivated operation to impose on the HF presentation above.

For other presentations the qualification matters. Kaye and Wong analyze
inverse interpretations between PA and finite ZF and emphasize the
transitive-containment condition, or a suitable membership-induction
equivalent in their base theory. See [*On interpretations of arithmetic and
set theory*, NDJFL 48 (2007), 497-510](https://research.birmingham.ac.uk/en/publications/on-interpretations-of-arithmetic-and-set-theory/),
DOI `10.1305/ndjfl/1193667707`. Do not replace this with an unqualified
equivalence for every theory called `ZF minus Infinity`.

The following explicit standard-model calculation isolates a more surprising
point. Define the Ackermann bijection recursively by

```text
A(0) = empty,
A(n) = {A(i) : binary digit i of n is 1}.
```

Membership is one binary digit test. Adjunction is `n | (1 << m)`.
The subsets of `A(n)` correspond exactly to binary submasks of `n`.

**PROVED container lemma.** For every ordinary natural number `n`,

```text
B_n = {A(0), A(1), ..., A(n)} = A(2^(n+1)-1)
```

contains `A(n)`, is transitive, and contains every subset of every member.
Indeed a member code of `A(i)` is less than `i`, and a submask of `i` is
at most `i`. Both therefore lie among `0,...,n`.

Thus two superficially similar closure conditions must be distinguished:

```text
SC(B): forall x in B, forall y, (y subseteq x -> y in B).

PC(B): forall x in B, exists p in B,
       forall y, (y in p iff y subseteq x).
```

SC includes each subset separately. PC includes the **set of all subsets**
as one of the container's members. Every `B_n` satisfies SC. No nonempty
finite HF container satisfies PC: choose a member with maximal code `h`.
The code of its power-set object contains binary position `h`, because
`A(h)` is one of its own subsets. That power-set code is at least `2^h>h`,
contradicting maximality. This is a quantifier/aggregation distinction,
not a numerical accident in the first few examples.

The Ackermann map preserves equality, membership, adjunction and subsets
exactly. It does not preserve Collatz edges as membership edges. Its needed
sidecar for that target is the actual arithmetic operation and its guards.

## 4. A genuine ZF/ZFC reflection package

The relevant classical phenomenon is **reflection**, which adds a statement
about formulas to a statement about containers. Bagaria's
[*Large cardinals as principles of Structural Reflection*, pp. 2-3](https://arxiv.org/pdf/2107.01580)
recalls the characterization of ZF by Extensionality, Separation, Foundation,
and reflection in transitive sets closed under subsets of members.
Levy's [1960 paper, Theorem 6](https://msp.org/pjm/1960/10-1/pjm-v10-n1-p14-p.pdf)
gives a related reflection equivalence over its explicitly stated base S;
in particular its proof reflects both a formula and its existential closure
to obtain Replacement. These are primary research sources, not an assertion
that bare transitive closure is equivalent to ZF.

Here is a convenient explicit package, with a self-contained reverse
derivation. For each finite list `Phi` of formulas and each finite list of
parameters `a`, let `Ref(Phi;a)` assert the existence of `M` such that

```text
every a_i is in M;  Trans(M);  SC(M);
for each phi(x_1,...,x_r) in Phi,
  forall x_1,...,x_r in M, [phi(x_1,...,x_r) iff phi^M(x_1,...,x_r)].
```

`phi^M` means syntactically restrict every quantifier to `M`. The list Phi
is a metalinguistic finite list, not a variable ranging over formulas with
a universal truth predicate. Unlisted free parameters of phi are included
among its displayed arguments; the finite parameter list can be padded.

Take Extensionality and Foundation, all Separation instances, and all
these Reflection instances. **PROVED: this is equivalent to ZF.**

For the reverse direction, use Separation with the false predicate on any
object to obtain the empty set. For arbitrary `a,b`, choose a container
with `a,b in M`, and separate `{a,b}` from M. If `a in M`, transitivity
makes `union a` a subset of M, so Separation gives Union. SC makes every
subset of `a` an element of M, so Separation gives the power-set object.

For Replacement suppose `phi(x,y,p)` defines a unique y for every x in a.
Reflect **both** `phi(x,y,p)` and `exists y phi(x,y,p)` in one M containing
`a,p`. Each x in a belongs to M by transitivity. Reflection of the
existential formula supplies a witness y in M satisfying the relativized
relation; reflection of the relation makes it an actual witness. Separate
the image from M. This is the load-bearing witness transport.

For Infinity let `Succ(x,y)` say `y=x union {x}`, whose totality was just
proved from Pairing and Union. Reflect both Succ and its existential closure
in M containing the empty set. The same two-formula argument shows that M
is closed under actual successor. Thus M is an inductive set, as required.
Foundation and Extensionality were assumed, completing ZF.

Conversely, finite-formula reflection in ZF supplies a parameter-containing
`V_alpha` with the required equivalences. Such a `V_alpha` is transitive
and satisfies SC. This is the cited classical direction; no single set
reflecting **every** formula at once is required or asserted.

For an exact **two sentences plus one schema** ZFC package, choose:

```text
sentence 1: Extensionality;
sentence 2: Foundation and Choice;
one schema, indexed by (psi,Phi,a): Sep(psi) and Ref(Phi;a).
```

Choice can be any conventional pure-membership first-order formulation,
with ordered pairs and functions expanded by their definitions. Conjunction
and indexing combine the two schema families honestly; this is an explicit
repackaging, not a claim that the conjunction has become weaker or that
Choice follows from reflection. Omitting Choice gives ZF. Both packages
use only the binary relation of membership besides logical equality.

## 5. Exact integers, small membership graphs, and the orbit obstruction

There is an exact version of the small-carrier/unbounded-depth intuition.
Let `u_0=0` and `u_(k+1)=2^u_k`. Then

```text
A(u_(k+1)) = {A(u_k)}.
```

The recursively reachable membership graph is a chain with k+1 vertices
at stage k; its rank is k. Codes begin
`0,1,2,4,16,65536,2^65536`. A 65,537-bit integer therefore has this
particular seven-vertex representation. The constructor is fixed while
the possible depth is unbounded; this is a family of growing graphs, not
one fixed finite graph that stores unlimited information without a sidecar.

Every membership edge strictly lowers the Ackermann integer code and HF
rank. But under the shortcut Collatz map

```text
T(n)=n/2 if n even;  T(n)=(3n+1)/2 if n odd,
```

the edge `3 -> 5` raises HF rank from 2 to 3. Neither endpoint is an
Ackermann member of the other. Even taking a transitive closure need not
lower the numerical code: the members of `A(27)` have codes `{0,1,3,4}`;
their transitive closure adds code 2 and is `A(31)`.

**PROVED computability distinction.** Every prescribed finite orbit prefix
is computable and can be stored with its exact values and affine carries
as an HF object. The relation saying that a finite code is an actual
first-descent certificate is decidable by replay. Neither fact establishes
totality of the unbounded search for such a certificate.

For all positive n>1, existence of some k>=1 with `T^k(n)<n` is equivalent
to Collatz termination: convergence supplies such a k; conversely strong
induction on n applies to the smaller actual iterate. This equivalence
identifies the remaining universal coverage obligation without making an
independence claim about PA or ZFC.

## 6. A concrete replacement question: reflect the certificate witness

Let `D(n,c)` mean that c codes a finite, correctly replayed orbit segment
starting at n and ending below n. Let `q` be a proposed finite-state
summary and let `C(q(n))` be a decoded family of candidate certificates,
potentially carrying unbounded counters or recursive structure.
An exact analogue of the useful reflection step is

```text
exists c D(n,c)
    iff exists c in C(q(n)) D(n,c).                         (W)
```

The right-to-left direction is automatic once D includes exact replay.
The left-to-right direction says that the representation loses no genuine
witness. Even W does not prove a witness exists. One still needs a
construction producing the right-hand witness for each n in scope, with
a proved terminating certificate-construction rank.

This splits a confused obstacle into three testable obligations:

1. **Soundness:** accepting a code really supplies an actual descent.
2. **Witness retention:** representation/quotient preserves an existing
   certificate or an adequate smaller replacement.
3. **Coverage:** every source in the stated domain has a constructed code.

The [bounded carry-defect result](nextforest_20260926_boundary.md) already
does something substantive of this kind: its finite bank plus exact
source-dependent lifts covers an all-height class of inputs. It does not
need full set-theoretic reflection; its arithmetic lifting theorem supplies
the relevant witness transport directly. The next useful question is
whether a proposed larger forest has an analogous lifting theorem and a
coverage argument, including its residual branches.

One should not assume that every finite summary is inadequate. A finite
controller with an unbounded exact counter may retain what is needed.
Nor can one infer coverage from the existence of the carrier or from its
well-founded membership relation: the certificate-construction operation
must actually decrease a specified rank. Defining that rank by the
unknown remaining Collatz stopping time would only move the open obligation.

## 7. Reproduction and scope

Run from the repository root:

```text
python 04-computation/experiments/reset_20260926_logic.py
python -O 04-computation/experiments/reset_20260926_logic.py
```

[The script](../../04-computation/experiments/reset_20260926_logic.py) uses
explicit exceptions, so optimized Python does not remove its checks.
The independent paths are recursive literal finite sets versus bitmask
formulas. It checks 4,096 round trips and transitive closures, 4,096
adjunction pairs, 32 powersets, 128 containers, the rank/closure hostiles,
and seven singleton-chain levels: **57,809 exact checks**. The frozen
[output](reset_20260926_logic.out) records the universes and counts.
These finite controls audit coding identities and boundary examples.
They do not prove a schema, a metatheoretic equivalence, or Collatz.
