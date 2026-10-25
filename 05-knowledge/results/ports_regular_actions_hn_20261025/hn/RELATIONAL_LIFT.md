# Keeping the same exterior colouring through every local test

2026-10-25. **PROVED:** an exact relational lift of the path/C6 interface,
its separator rule, and an unpartitioned finite-menu certificate for black
detours. **FINITE-EXACT:** the small controls below. These are certificate
constructions, not a new chromatic lower bound or all-order connectivity.

The inherited mechanisms are the [q+2 path automaton](README.md), the
[circle-port decomposition](../../circle_ports_polygon_palettes_20261021/README.md),
and the [affine ledger compiler and menu criterion](../../trijection_hourglass_ports_20261025/trijection/README.md).
The canonical hostile is an exterior colouring that fails to extend even
though another exterior colouring succeeds. The corrected near miss is
replacing correlated lists by their independent marginal possibilities.
The least-used sidecar is a *relation* between the unfinished exterior
variables and the local automaton state. The live concepts are shared
variables, existential projection, universal failure, signed constraints,
and separator size.

## 1. An exact finite relation for the universal quantifier

Let G be a finite graph and C a named palette of size q>=2. Partition its
vertices into an exterior B and an interior U, where G[U] is a disjoint
union of paths and six-cycles. This includes the radius-one neighbourhood
of a planar point: put the centre itself in B and its actual neighbours
in U. Let H=G[B]. Define

    R_H = {beta in C^B : beta is a proper colouring of H},
    L_v(beta) = C minus {beta(b) : b in B adjacent to v}.

The centre colour is one coordinate of beta. Fixing it is an additional
relation filter; allowing every centre colour means retaining that
coordinate throughout. Exterior constraints on non-neighbouring vertices
are included in R_H and are not replaced by per-port colour sets.

For each component J of G[U], define A_J to be the relation on B consisting
of those beta for which the lists L_v(beta), v in J, are colourable. Compute
this relation by the q+2 path automaton or its anchored C6 version.

**Relational extension lemma.**

    G is q-colourable  iff  R_H intersect intersection_J A_J is nonempty.

Equivalently, universal failure of extension means

    for every beta in R_H, there exists J whose local automaton rejects.

**Proof.** Restriction of a proper colouring of G supplies one beta and a
list colouring of every component. Conversely, for one beta in the stated
intersection, choose a list colouring of each component. Components have
no edges between them, and each list excludes precisely the colours of
its exterior neighbours. Their union with beta is therefore proper. The
same beta is used for every component. QED.

For C6, the local relation can be obtained by introducing an anchor a in
its first list, starting the path state at forced(a), and requiring that
the final attainable colours contain one different from a. Existentially
project a only after that closing-edge test. An empty relation means that
every possible anchor failed, with no exchange of quantifiers.

## 2. A unit-distance hostile to independent list products

Take a regular unit pentagon, with cyclic vertices o,r,u,v,s. Put
U={r,s} and B={o,u,v}. The two rim components are isolated vertices; the
only exterior edge is u--v. Use C={0,1} and fix beta(o)=0. Its two exterior
colourings give exactly the two list pairs

    (empty, {1}),   ({1}, empty).

The exact relation contains no pair on which both components extend.
Each component separately extends under some exterior colouring, but the
successful exterior colourings are different. The independent product of
the two list projections has four pairs and introduces the false pair

    ({1}, {1}).

In automaton-state terms, the exact relation is

    {(dead, forced(1)), (forced(1), dead)},

whereas its product completion also contains
(forced(1),forced(1)). The obstruction is a lost correlation even though
each local state is exact. This is a genuine planar unit-distance example
at q=2; it makes no six-colour claim. The q=3 positive control has 30
proper colourings.

Because a product of projections is an overapproximation, proving failure
on that entire product would still be sufficient. Its defect is loss of
completeness: an accepted product row can be fictional. Concluding that a
common extension exists merely because each component has some extension
is the invalid quantifier exchange.

## 3. A small local interface, with explicitly measured correlation memory

The relations need not be stored as full tables on every exterior vertex.
Write the problem as a finite constraint network containing exterior edge
inequalities, list-determination factors, automaton transitions, and cycle
anchor/closure factors. Process factors in any order compatible with the
chosen path scans. At a cut retain every variable occurring in both the
processed and unprocessed portions, together with the active automaton
state and any unclosed cycle anchor. Store the *joint relation* on these
variables, not their separate allowed values.

**Separator projection lemma.** If R(x,z) is the processed relation and
all future factors S(z,y) are independent of x, eliminating x now is
exact for every future completion:

    exists x [R(x,z) and S(z,y)]
        iff [exists x R(x,z)] and S(z,y).

Apply this identity at each join-and-project step. Two histories may be
merged exactly when they have the same retained assignment and local
automaton state; omitted histories are joined existentially. The invariant
is that the table contains precisely the assignments extendible through
the processed factors. Induction proves the final empty/nonempty answer.

If at most w exterior colour variables remain active and components are
scanned one at a time, an intermediate table has at most

    q^w (q+2) rows for a path,
    q^(w+1) (q+2) rows for an anchored six-cycle.

These are bounds for a specified factor order, not constant-size bounds
for arbitrary exteriors. Deriving R_H itself has not been made free.
After filtering a completed component for acceptance, its state and its
anchor can be projected away; the remaining boundary correlation stays.

This gives an exact finite certificate language: initial finite factor
tables, successive relational joins and projections, deterministic DFA
transitions, closing-edge filters, and a final empty table. Each step is
locally checkable by finite table operations. Conversely, for every finite
instance the full-variable tables provide such a derivation whenever the
final relation is empty. The certificate can be exponential; no bounded
width or polynomial-size unsatisfiability certificate is asserted.

There is an elementary lower bound on the required sidecar at a prescribed
cut. Let the exterior be w disjoint paths x_i--z_i--y_i and use two named
colours. Projecting out z_i gives the exact relation

    (x_1,...,x_w)=(y_1,...,y_w).

After reading all x coordinates, its 2^w different assignments have
distinct continuation languages: the suffix y=x distinguishes each one.
These disjoint length-two paths admit unit-distance realizations. Thus a
constant local DFA does not bound the memory of arbitrary boundary
correlations. A different factor order can avoid this particular large
cut; the claim is about the stated cut, not an invariant lower bound over
all processing orders.

## 4. The exact transfer to black-detour menus

Fix a literal source family F of labelled tournaments described by unary
or binary parity equations on arc bits. Each candidate ledger j consists
of Hamiltonian words and has already passed the word-permutation and
literal-target XOR checks. Let

* L_j be the conjunction of its actual intermediate join equations;
* B_(j,i) be the full affine system saying its i-th word is blue.

The ledger succeeds at a source x precisely when

    L_j(x) and, for every i, not B_(j,i)(x).

The source x is shared across every ledger, every join and every blue
test. Intermediate flips are compiled into constants as in the inherited
compiler. Define a failure-alternative family W_j consisting of

1. one single negated join equation for each equation in L_j;
2. one whole system B_(j,i) for each move i.

Negating a unary parity equation just toggles its constant. Each member
of W_j is therefore still a conjunction of unary/binary parity equations.

**Unpartitioned menu lemma.** The candidate ledgers cover every source in
F by legal all-black detours iff, for every choice W_j in the respective
failure-alternative families, the system

    F together with all equations in all the chosen W_j

is inconsistent.

**Proof.** A ledger fails exactly when at least one of its joins is false
or at least one of its moves is blue. An uncovered source therefore
chooses one failure alternative for each ledger and satisfies their common
intersection. Conversely, a source in one such intersection defeats every
ledger. Distributing the finite unions gives precisely the criterion.
The source must be common; separate witnesses against separate ledgers
do not demonstrate an uncovered source. QED.

All the selected systems are signed pair-constraint graphs, including a
constant vertex for unary equations. An intersection is inconsistent iff
it contains a cycle whose equation constants XOR to one. Thus one such
signed cycle for each selection is a complete finite coverage certificate.
Self-loops and parallel edges must remain available in this constraint
graph. They are proof constraints, not tournament loops.

This extends the inherited all-legal menu criterion to variable legality
without first partitioning by every boundary-join bit. If there are t_j
distinct join equations and m_j moves in ledger j, the direct expansion
has at most product_j(t_j+m_j) selections. Factoring, shared tables, or
redundant-branch elimination may reduce it; no uniformly small expansion
is promised. When all joins are already implied by F, the negated-join
alternatives are immediately inconsistent and the criterion reduces to
the inherited product_j m_j blue-choice test.

For the actual [endpoint-four Theorem A](../ports/README.md), put
x=T(2,4), y=T(2,6), z=T(3,5). Its four guarded ledgers occupy

    x=0;  x=1,y=0;  x=1,y=1,z=1;  x=1,y=1,z=0.

These guards partition every source family, whether or not x,y,z are
independent. Each guarded leaf has a proof of all joins, all three
blackness conditions, and the same literal target. This is a concrete
relational decision certificate, retaining every port and every shared
endpoint. It establishes this family, not the existence of an all-order
menu for arbitrary port sets.

## 5. What the connection preserves

| Source | Target representation | Preserved predicate | Information that must remain |
|---|---|---|---|
| Exterior colouring plus path/C6 lists | Boundary relation joined to local DFA | Existence of one common full colouring | Shared exterior colours, active state, cycle anchor |
| Literal tournament source plus detour menu | Joint failure relations on arc bits | Existence of one source defeating every ledger | Exact ports, intermediate flip constants, complete blue systems |
| A processed relation and future constraints | Projection onto their separator | Every possible future satisfiability test | Every variable used on both sides |

The map between the two research lanes is compilation into the same
relation/join/projection calculus. It is not a graph map from a unit graph
to a tournament, and it does not identify palette colours with ports.
The transferable rule is exact: preserve shared variables until their
last use; only then eliminate them. Counting multiplicities and clone
transport require further coordinates and do not follow from Boolean
satisfiability preservation.

The next concrete target is a detour menu whose failure-relation factors
have uniformly controlled separators or a repeating signed-cycle
obstruction. The corresponding HN target is a realizable exterior whose
joined acceptance relation is empty. Neither target is solved merely by
having a small local automaton.

## 6. Small independent controls

The checker below uses literal proper-colouring enumeration for C5, the
exact correlated list relation, and 300 synthetic parity-source/menu
families. For the latter it compares direct source-by-source coverage
against the unpartitioned criterion and independently checks all 19,200
selected signed systems by parity propagation and literal assignments.
The menu examples audit the algebraic lemma, not new tournament families.

Run from this directory:

    python -c "from pathlib import Path; t=Path('RELATIONAL_LIFT.md').read_text(); exec(t.split(chr(96)*3+'python\n')[1].split(chr(96)*3)[0])"

```python
from itertools import product
from random import Random
import json

def need(p,m):
    if not p: raise RuntimeError(m)
def col(q):
    return [x for x in product(range(q),repeat=5)
            if all(x[i]!=x[(i+1)%5] for i in range(5))]
need(len(col(2))==0 and len(col(3))==30,'pentagon controls')
R=set()
for u,v in product(range(2),repeat=2):
    if u!=v:R.add((frozenset({0,1}-{0,u}),frozenset({0,1}-{0,v})))
need(len(R)==2 and not any(all(L for L in row) for row in R),'joint relation')
P=set(product({row[0] for row in R},{row[1] for row in R}))
need(sum(all(L for L in row) for row in P)==1,'marginal hostile')
X=list(product((0,1),repeat=3))
def satisfies(x,eqs):
    z=x+(0,)
    return all(z[a]^z[b]==v for a,b,v in eqs)
def consistent(eqs):
    adj=[[] for _ in range(4)]
    for a,b,v in eqs:adj[a].append((b,v));adj[b].append((a,v))
    seen={}
    for start in [3,0,1,2]:
        if start in seen:continue
        seen[start]=0;stack=[start]
        while stack:
            a=stack.pop()
            for b,v in adj[a]:
                if b in seen:
                    if seen[b]!=seen[a]^v:return False
                else:seen[b]=seen[a]^v;stack.append(b)
    return True
rng=Random(30030);intersections=0
families=[[],[(0,1,0)],[(0,1,1),(1,2,0)]]
for run in range(100):
    menus=[]
    for j in range(3):
        legality=[(rng.randrange(3),3,rng.randrange(2)) for _ in range(2)]
        blues=[[(rng.randrange(4),rng.randrange(4),rng.randrange(2))
                for _ in range(2)] for _ in range(2)]
        menus.append((legality,blues))
    alternatives=[ [[(a,b,1-v)] for a,b,v in L]+B for L,B in menus ]
    for F in families:
        covered=all(any(satisfies(x,L) and not any(satisfies(x,B) for B in Bs)
                        for L,Bs in menus) for x in X if satisfies(x,F))
        expanded=True
        for choice in product(*alternatives):
            eqs=F+[e for option in choice for e in option]
            c=consistent(eqs)
            need(c==any(satisfies(x,eqs) for x in X),'parity consistency')
            intersections+=1
            if c:expanded=False
        need(covered==expanded,'unpartitioned menu equivalence')
print(json.dumps({'status':'PASS','C5_colourings_q2_q3':[0,30],
    'joint_list_rows':len(R),'independent_product_rows':len(P),
    'phantom_accepting_rows':1,'menu_families':300,
    'signed_systems':intersections},indent=2))
```
