# Five indecomposable modules arise from the actual interface algebra

**PROVED elementary classification and transfer interpretation; CITED classical
algebra identification. Chan–Pak remains OPEN.** This is a different, intrinsic
five-element collection. The five tournament cores are operators, not five
modules. No bijection between these two collections is asserted.

## 1. Inheritance and the new connection

The [concurrent Session252 algebra](../../session252_fivefold_bridges_20261111/algebra/PROOF.md)
identifies the polynomial residual-path interface with

    M(d,l,o,i,j) = [[d,o,j],[0,l,i],[0,0,d]],
    M(C -> D) = M(C) M(D).

At any fixed path-weight value these matrices belong to the five-dimensional
algebra A over the chosen field. Session252 already identifies A with
End_R(R plus k), R=k[epsilon]/(epsilon^2). That algebra identification and
the five-dimensional statement are inherited. The additional connection here
is to its complete five-object module classification, with an elementary
rank reconstruction and a concrete interpretation of the lost central channel.

This is a classical algebra: the Auslander algebra of the dual numbers, a
Nakayama algebra with projective lengths 3 and2. Its five indecomposable
modules are stated in René Marczinzik, *Dominant dimensions of algebras*
(2017), Example4.3.5, printed page43,
[primary thesis](https://elib.uni-stuttgart.de/server/api/core/bitstreams/4a8f7dec-f65d-4712-8bc1-9066e85529db/content).
The classification below is proved directly; it is not claimed as new
representation theory. The useful new bridge is its application to this
specific tournament interface and the earlier actual unit switches.

## 2. The five modules and a direct completeness proof

Let e1=E11+E33, e2=E22, a=E12, b=E23. Matrix multiplication gives

    ba=0, ab=E13=:c, c^2=0.

Together e1,e2,a,b,c are a basis of A. A finite-dimensional left A-module
is therefore precisely two vector spaces V,W with maps

    a: W -> V,       b: V -> W,       ba=0.                (1)

There are exactly five indecomposable types. The labels describe the
successive nonzero maps, with every displayed arrow an identity between
one-dimensional spaces:

|Type|Dimensions (V,W)|Nonzero chain|
|---|---:|---|
|S1|(1,0)|1|
|S2|(0,1)|2|
|L12|(1,1)|1 -> 2|
|L21|(1,1)|2 -> 1|
|L121|(2,1)|1 -> 2 -> 1|

In L121 use V=(v_top,v_bottom), W=(w), with b(v_top)=w,
a(w)=v_bottom, and b(v_bottom)=0. This fixes the direction convention.

For a general pair (1), put ra=rank(a), rb=rank(b), r=rank(ab).
Its multiplicities in this list are exactly

    m121 = r,
    m21  = ra-r,
    m12  = rb-r,
    m1   = dim(V)-ra-rb,
    m2   = dim(W)-ra-rb+r.                               (2)

All quantities are nonnegative: im(a) is contained in ker(b), and
rank(ab)>=ra+rb-dim(W).

Here is a basis construction proving existence, not just dimension matching.
Choose W3 complementary to im(b) intersect ker(a) inside im(b), and let
Wb=im(b) intersect ker(a). Choose Wa so that a(Wa) complements a(W3)
inside im(a), and W0 complementary to Wb inside ker(a). Then

    W = W3 plus Wb plus Wa plus W0.

Indeed W3+Wb=im(b); Wa can be chosen outside im(b)+ker(a), and W0
fills the remaining kernel. Lift bases of W3 and Wb through b to V.
Together with the images a(W3),a(Wa), they are independent; complete
them by V0 inside ker(b). The paired basis vectors give respectively
L121, L12, L21, S2 and S1 summands. Formula(2) follows.
The five blocks are indecomposable: the simples are immediate; either
length-two nonzero arrow forces its two one-dimensional spaces into the
same summand; in L121 both arrows are nonzero and link all three basis
positions (equivalently its endomorphisms form k[epsilon]/epsilon^2,
which has no nontrivial idempotent). Thus the list is complete and the
five rank quantities determine the isomorphism class over every field.

## 3. What this five does, and does not, say about the tournament five

The natural three-coordinate matrix carrier is L121: take V=span(e1,e3)
and W=span(e2). EVERY coloured tournament core acts on this same module.
Choosing K1,C3,S4,S5a,S5b chooses five actual operators in A, with integer
polynomial entries and severe counting constraints; it does not choose
one operator from each module type. For example S1 observes only d=H(C),
so it cannot distinguish the two H9 cores. The module classification
supplies no choice of a core, cut or address for a prescribed target N.

The radical is span(a,b,c), its square is span(c), and its cube is zero.
The four smaller indecomposables all kill c. The length-three carrier is
the only indecomposable retaining the two-step response ab, precisely the
central j coordinate. Thus discarding that channel has an exact module
meaning, rather than being a vague loss of information.

Under transpose together with exchange of the first and third matrix
coordinates, a and b are exchanged. The induced duality fixes S1,S2,L121
and swaps L12,L21. Its orbit pattern is 1+1+1+2. This differs both from
Platonic duality's 1+2+2 and from converse fixing all five small cores.
It must not be used as an equivariant identification of those collections.

## 4. A physical transfer: central unit switches

The actual pair in
[Session229 rich modules](../../session229_combined_constructions_20261109/rich/PROOF.md)
uses the strong core B95=(64,121,35,37,13,81,28) and apex outmasks12,20.
That packet proves equality of the six old endpoint coordinates and a
two-unit difference in the remaining joint-cover coordinate. The literal
children have counts351 and353; its additional star context gives681,683.
These witnesses and their propagation theorem are inherited. The full
polynomial endpoint channels I(z),O(z) agree for the named pair, and

    J20(z)-J12(z) = -2z^2-6z^3+z^4+3z^5+2z^6.

Thus centrality holds before choosing z, not merely after evaluating at
plus/minus one. These coefficient equalities are explicitly checked.

In the current seven-coordinate compiler, the two cuts have identical
(d,L,i,o,s,e), and h differs by2. Consequently every pair of arbitrary
binary flank addresses gives consecutive odd outputs. `check_modules.py`
independently checks these coefficient equalities and small actual graphs.
At the matrix level, their difference in each evaluated interface is a
scalar multiple of c. Since M c=c M=d c, fixed left/right contexts multiply
that difference by H(left)H(right). Transitive contexts have H=1 and
preserve it. This recovers the old central-defect mechanism inside A.

This is the constructive content of the five-module bridge: the largest
indecomposable identifies exactly the response that can store a unit switch,
and the algebra proves its legal propagation. The four smaller modules
would all erase it. The inherited hostile remains decisive: k copies with
independent central switches produce only k+1 scalar levels, spaced by
2*d^(k-1), because c^2=0. They do not produce a complete radix digit set.
The absolute lower endpoint is not freely prescribed; the arbitrary-N
lifting step remains open.

## 5. Reproduction and exact scope

`python -I -B check_modules.py` independently enumerates every pair of
binary-field maps with dim(V),dim(W)<=3 and dim(V)+dim(W)<=5, filters ba=0,
and compares their full change-of-basis orbits with the five-block normal
forms. It also checks the matrix relations, reconstructs the physical
B95 cut pair, and recounts its small flank contexts. The field census is
a control for the all-field basis proof, not its justification. No claim
about arbitrary multi-apex observers or a five-state witness machine follows.
