# Ternary ports, golden defects and coefficient identifications

**Status: PROVED elementary identities and scoped family classification;
FINITE-EXACT geometric tests; CITED classical star-triangle context;
OPEN six-colour forcing and theorem transfer.**
Session `clock_duality-20261018j`. The results are self-audited; no priority
claim or new canon ID is assigned.

The useful common object is a labelled boundary relation. A graph move,
unit walk or quotient that preserves connectivity can destroy that
relation. Three-port elimination makes the loss explicit: at six colours
the missing term is `-19/9` on the all-equal boundary state. Its numerator
is a golden-field norm. At seven colours the corresponding parameter is
exactly the previously studied `-29/16`, and seven is the unique integer
palette in this family whose associated quadratic has a rational
three-cycle. These are exact algebraic connections, not implications for
the chromatic number of the plane or for Collatz convergence.

## 1. Inheritance before adding another graph analogy

The incoming [graph/Pythagorean synthesis](graph_pythagorean_20261018_synthesis.md)
already proves the Seidel `3,4,5` triple, its infinite conference-family
continuation, and the ordered boundary-vector construction. This session
uses those stronger results rather than rederiving their scalar shadows.

- [THM-261, Petersen root orthogonality](../../01-canon/theorems/THM-261-petersen-root-orthogonality.md):
  `D(K5)=Petersen`, with `D(G)=complement(L(G))`.
- [THM-4116, boundary-state gluing](../../01-canon/theorems/THM-4116-boundary-state-gluing-and-ap-odd-shell-tree-synchronizers.md):
  ordered partial Tait-colouring tensors glue by a dot product. Petersen's
  two sides can each be colourable while their supports are disjoint.
- [THM-4529, Petersen and Heawood in Paley coordinates](../../01-canon/theorems/THM-4529-petersen-and-heawood-families-in-paley-coordinates-anti-circulant-all-odd-tournaments.md):
  the Petersen Delta-Y family is the class of K6, not a class containing
  bare K5 and K3,3. Its Fano/Pasch construction and parity limitations stand.
- [THM-4690, the 19-adic copy obstruction](../../01-canon/theorems/THM-4690-first-arrival-ladder-port-and-integer-clock-decompositions-and-the-nineteen-adic-copy-obstruction.md)
  excludes particular step alphabets; it does not separate all copies in
  the enlarged field. Its corrected sheet scope is essential.
- [THM-4700, unit-path ring and reflection links](../../01-canon/theorems/THM-4700-the-unit-path-ring-of-the-six-colour-field-short-reciprocal-walks-and-the-heron-census-of-coefficient-triangles.md)
  supplies the reflection construction. Its global classification is
  still labelled PRE-AUDIT; the links used here are checked directly.

The old [Kuratowski/Tait note](procgen_kuratowski_20260925_tait_kempe_triples.md)
already shows why a minor obstruction cannot prove Collatz entry: a
functional graph has at most one cycle per component. The graph is planar
whether or not every positive component reaches ROOT.

## 2. An explicit ternary coordinate chart for the graph triple

Represent Petersen by the two-subsets of `{0,1,2,3,4}`, with disjointness
as adjacency. Set `v={3,4}` and Seidel-switch across its three neighbours.
Then v is isolated. Label the other nine vertices by `F3^2`:

```text
{0,1,2} minus {i}  -> (i,i),
{i,3}             -> (i+1,i-1),
{i,4}             -> (i-1,i+1).                         (3)
```

The switched graph is the rook graph: two points are adjacent iff one
coordinate agrees. It is `L(K3,3)`. The invertible ternary map
`(a,b)->(a+b,a-b)` interchanges this relation and its complement.
This is a labelled verification of the incoming synthesis's isomorphism.

There is also a norm interpretation: in `F9=F3[i]`, `i^2=-1`, adjacency is

```text
Norm((a+bi)-(c+di))=(a-c)^2+(b-d)^2=1.
```

Thus the two-coordinate ternary grid, the line graph of K3,3, and this
finite-field unit-circle graph are literally isomorphic. The Petersen
step requires switching; it is not an isomorphism of the original graphs.
The positive `A4`-root presentation of K5 edges remains another exact
representation, with its original pair labels.

**Preserved:** the declared adjacency or switched adjacency, including
nonedges. **Lost:** original vertex degrees, planarity, geometric unit
realization and colour-extension tensors. **Sidecar:** the switch shore,
the edge-to-vertex map and the ordered ports. The independent script checks
all pairs in (3). These losses forbid applying the finite-field colouring
or Seidel spectrum directly to the planar unit-distance problem.

The role of three is exceptional in this chart. Over any odd finite
field F_q, the norm-one graph on F_(q^2) has degree q+1; the q by q rook
graph has degree 2(q-1). These degrees agree only at q=3. The Paley graph
on F_(q^2), used in the inherited infinite conference family, instead
joins differences whose norm is a nonzero square in F_q, and has degree
(q^2-1)/2. At q=3 the only nonzero square is 1, so all three descriptions
coincide. Beyond q=3, one must change the adjacency predicate: simply
inflating the ternary chart does not preserve unit-norm adjacency. The
norm's fibres have q+1 elements, proving the degree statements. Exact
controls at q=3,5,7,11,19 check both norm predicates and this failure
boundary. This distinguishes a genuine recursive family from an
accidental equality of its smallest members.

## 3. The exact three-port term lost by a star-triangle replacement

Let a,b,c be boundary colours from q labelled colours. Write
`delta_ab=[a=b]` and `delta_abc=[a=b=c]`. Eliminating the centre of a
properly coloured three-leaf star gives

```text
Y_q(a,b,c)=q-|{a,b,c}|
         =q-3+delta_ab+delta_bc+delta_ca-delta_abc.       (4)
```

For q>3 put `C=q-3`. The unique symmetric pair-only product matching the
all-distinct and exactly-two-equal entries is
`C product_(pairs)(1+delta_pair/C)`. Its full correction is

```text
Y_q = C product_(pairs)(1+delta_pair/C)
      - [(q^2-3q+1)/(q-3)^2] delta_abc.                 (5)
```

**Proof.** Products of two or three equality deltas equal `delta_abc`.
The pair-only product has triple coefficient `3/C+1/C^2`, whereas (4)
has coefficient -1. Their difference is
`-(C^2+3C+1)/C^2`. QED.

At q=6, the three star entries are 5,4,3, all positive. Equation (5) has
the correction `-19/9`. An ordinary triangle instead permits only three
distinct boundary colours. A bare Delta-Y move therefore does not
preserve six-colour boundary relations. If the centre itself is retained
as a port, elimination is invalid without its conditional tensor.

This is the classical two-site/three-site interaction issue in Potts
star-triangle transformations. The primary source
[Wu and Guo, Physical Review E 86 (2012), Eqs. 9--12](https://physicsfaculty.bnu.edu.cn/Public/web/application/research/statistics_group/research/papers/pre2012_wu_guo.pdf)
gives the three boundary weights, the three-site interaction, and the
same polynomial `q^2-3q+1` at the honeycomb zero-temperature endpoint.
No novelty is claimed for that polynomial or the local transform.

For comparison, a unit path of length two has boundary matrix

```text
(J-I)^2=(q-2)J+I.
```

At six colours every endpoint assignment has 4 or 5 extensions. More
generally every path of length at least two extends every endpoint
assignment for q>=3. Consequently subdividing every edge of **any** graph
does not force a relation on its old vertices, if the only added
constraints are the path edges and their internal vertices are disjoint.
Additional authenticated contacts can change this conclusion.

This pinpoints the obstacle with long coefficient connectors: a unit
walk pays connectivity, while a forcing gadget must pay an actual zero
of the boundary-extension tensor.

## 4. Golden norms connect 11, 19, 29 and the old 76

Let `phi^2=phi+1`. The exact identity is

```text
D(q)=q^2-3q+1=Norm_(Q(sqrt5)/Q)(q-phi^2).
```

For q=4,5,6,7,8,9 its values are 5,11,19,29,41,55. This is a norm
identity, not a claim that all values are prime. The vanishing parameters
are `phi^2` and `phi^-2`. Their chromatic/Potts role is classical; it does
not transfer a golden identity at a nonintegral parameter to six-colour
extendability. In particular (5) requires a signed three-port correction
at q=6, not a positive graph substitution.

The older [golden phase recursion note](difference_diagonals_golden_phase_recursion_20261004.md)
already records `Norm(phi^7-1)/16=-29/16` and the ideal
`(29,phi-24)`. There is now an exact unit-and-conjugation comparison:

```text
phi^5-1 = phi^3 (5-phi^2),
phi^7-1 = phi^3 conjugate(7-phi^2),
phi^9-1 = 2 phi^5 (6-phi^2).                            (6)
```

Since `Norm(phi)=-1`, these give norms -11, -29 and -76. In particular
the six-colour numerator 19 and the inherited 76 differ by the norm of
the rational factor 2 in the last identity. This compares actual
algebraic integers, not merely their decimal values. For q=7 the
conjugation in (6) is necessary: `7-phi^2` vanishes at `phi=6 mod 29`,
while `phi^7-1` vanishes at the conjugate residue 24.

The 19-adic coefficient obstruction and the first conjugation-fixed
prime 29 in THM-4700 have their own local-field proofs. Equations (5)--(6)
do not identify those valuation predicates with a colour-extension
predicate. They give a concrete common arithmetic place to investigate,
with the missing map stated explicitly.

The polynomial also occurs in the owner-tail algebra of
[THM-4033](../../01-canon/theorems/THM-4033-prime-sector-ap-cover-eventual-owner-tail.md),
where `q=h+2` gives `n_-=q^2-3q+1`. That prior formula is
retained as an algebraic overlap; no source/owner map to (4) has been
constructed.

## 5. A scoped dynamical classification of the correction family

Define, for integer q>=4,

```text
c_q=-(q^2-3q+1)/(q-3)^2,       f_q(x)=x^2+c_q.
```

**Proposition T (PROVED).** A rational point of exact period three for
f_q exists if and only if q=7.

*Proof.* For an exact three-cycle of `x^2+c`, its trace s satisfies
`s^2+s+c+2=0`. This follows by reducing the trace polynomial
`s=x+f(x)+f^2(x)` modulo `(f^3(x)-x)/(f(x)-x)`; the remainder of
`s^2+s+c+2` is identically zero. Thus `-4c-7` must be a rational square.
Put `d=q-3>=1`. Here

```text
-4c_q-7 = (-3d^2+12d+4)/d^2.
```

For d>=5 this is negative. At d=1,2,3,4 its numerator is 13,16,13,4,
leaving only d=2 and d=4. For d=2, `c=-11/4`: valuation growth forces
every rational periodic point to have exact denominator 2. The escape
bound `|x|<5/2` leaves only `+-1/2,+-3/2`, all of which reach `-5/2`
and escape. For d=4 the cycle is

```text
-1/4 -> -7/4 -> 5/4 -> -1/4,        c=-29/16.
```

This proves the stated family classification. QED.

The cycle, its multiplier `35/8`, and the entire rational preperiodic
graph were already proved more strongly in
[THM-4139, rational three-cycle and order-six lift](../../01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md).
The extension here is the exhaustive integer-q *star-correction family*,
not a rediscovery of that cycle or a general rational-cycle theorem.
The parabolic parameter `-7/4` corresponds to the nonintegral values
`q=5+-4sqrt(3)/3`; it is distinct from the rational cycle containing
the point `-7/4` at parameter `-29/16`.

No dynamics on graph colourings has been conjugated to f_q. The map
proved here is the scalar parameter extracted from the exact boundary
tensor, followed by a separate arithmetic classification.

## 6. Authenticated coefficient contexts

Use the inherited five independent K-coordinates
`(1,u,v,uv,Y)`, where `K=Q(zeta15,sqrt(-11))`, `Y=e^2`, and

```text
u^2=a+b u+c Y,    a=conjugate(1+zeta5).
```

The starting factors are `A+aM`, `C+bM`, `C`, `C`, `C+cM`, with the
original eight ports at `(p,0,0,0,0)`. The 570 patterns are all six-colour equality patterns on the eight
original ports, up to colour permutation, satisfying the inherited port
edges: 01,12,23,34,40,05,56,67,70. All results below are for the axes
and the explicitly named unit contact families, **not a complete unit
difference census** in this five-coordinate space.

### Reflection links

Writing `h=s/k`, the inherited reflection construction becomes the
Möbius map

```text
R_s(q)=(q-h)/(1-conjugate(h)q).
```

For unit q and s, with `Norm(k)!=1` and denominator nonzero, it has unit output
and `|s+k(R_s(q)-q)|=1`. Its inverse is `R_(-s)`. The script checks both
the output norm and the actual cross-copy edge in characteristic zero.
The independent audit also proves the two norm identities and the inverse
identity formally, substituting only q*conjugate(q)=s*conjugate(s)=1.
Tangencies may return q itself; the selected new points are checked to
be distinct. This normalized formula exposes a reversible generator
action without claiming that its unmarked orbit retains port relations.

Add one new reflection label for k=b and one for k=c. **FINITE-EXACT:**

```text
factor sizes:          472,90,10,10,90
implicit vertices:    382,320,000
axis edge counts:     1430,137,15,15,137
boundary patterns:    all 570 extend
direct witness checks:145,482,240
```

As a hostile control, requiring the b/c colour factors to ignore their
copy labels blocks 123 of those patterns in the restricted additive
model. The unrestricted additive model extends every one. Thus copy
labels are active information; that restriction is not a geometric
forcing certificate.

### Literal overlaps

For a stronger connection, add the four common coefficient labels

```text
1/b, 1/b+1, 1/c, 1/c+1
```

to M. For example `0+b(1/b)=1+b(0)`, so the b-axis copies literally share
a vertex. The pair `1/b,1/b+1` also has unit difference, creating another
u-square contact. These identities bypass the need to construct a long
unit walk to the copy. All old ports are retained in their old positions.

**FINITE-EXACT:**

```text
factor sizes:          566,106,10,10,106
implicit vertices:    635,957,600
axis edge counts:     1720,175,15,15,175
boundary patterns:    all 570 extend
direct witness checks:332,516,340
```

This kills the completely copy-blind colouring shortcut: assigning the
same base colour to all representations of the shared vertex would
identify the colours of base points 0 and 1, which are adjacent. It still
does not kill the general additive witnesses. Connectivity and literal
identification are now paid; palette compatibility remains unpaid.

### Changing the cubic generator

The cubic relation is `u^3=m2 u^2+m1 u+1`. Hence

```text
u^-1=(a-m1)+(b-m2)u+cY.
```

This supplies another authenticated three-coordinate unit circuit,
sharing its Y coefficient with u^2. It is a concrete changing-generator
test, not an extra independent Cartesian axis. **FINITE-EXACT:** factor
sizes `688,130,10,10,70`, hence 626,080,000 implicit vertices. All 570
patterns extend, with 268,113,180 direct witness checks. The extra circuit
has 189,000 constraints. Receipt: `golden_inverse.json`.

### Extension-complete pieces cannot create a new port restriction

**Proposition G (PROVED).** Let each graph G_i have marked boundary B_i
and a declared boundary graph H_i, with every proper six-colouring of H_i
extending to G_i. Glue any finite collection by identifying boundary
vertices, retaining disjoint interiors and adding no other edges. Then
its six-colour boundary relation is exactly that of the glued union of
the H_i. The same holds with additional edges solely between boundary
vertices, provided those edges are included in the declared union.

*Proof.* Every colouring of the glued graph restricts to its boundary
union. Conversely restrict a proper colouring of that union to each H_i
and extend it inside G_i. Interior disjointness makes the choices
compatible; the already fixed boundary colours agree at every identified
vertex. QED.

The saved 570 witnesses, followed by global colour relabelling, establish
this premise for each selected graph above. Thus finite repetition and
port permutation of these pieces alone cannot strengthen their boundary
relation. This does not bound the chromatic number of the plane, nor does
it cover extra unit contacts between different interiors. Those contacts
are an explicit way around the obstruction. Petersen's complementary
five-pole supports show exactly what is absent here: locally nonempty
but globally incompatible boundary relations.

## 7. Reusable rule and remaining frontier

The common research move is **eliminate only with the exact interface**:

| Source | Exact retained interface | What a scalar shortcut loses |
|---|---|---|
| Collatz first-arrival packets | source-labelled Eisenstein polynomials, entry port and integer clock | phase orientation, arrival timing and moving mass |
| Three-port graph centre | equality-partition tensor (4) | the triple term in (5) |
| Cubic coefficient copies | coordinate identity and colour-extension relation | copy-dependent palette choices |
| Petersen decomposition | ordered Tait boundary tensor | compatibility despite positive local counts |

The [ROOT companion](root_eisenstein_first_arrival_20261018j.md) proves a
uniform linear mass bound but leaves quadratic concentration open.
Here the connector and overlap tests leave all 570 patterns extendable.
The next useful geometric target is a compatible multi-generator
boundary relation that defeats the saved colourings, with every claimed
edge independently norm-checked. Another bare unit path cannot do that.

Reproduce the finite controls with `ternary_interfaces_20261018j.py`.
The geometry runner is `golden_reflection_context_20261018j.py`, with
`--cubic --links 2`, `--cubic --links 2 --copy-blind`,
`--cubic --links 0 --overlap 2`, or `--cubic --links 0 --inverse`.
The JSON receipts save every extending additive colouring for direct
verification. The [manifest](../../04-computation/experiments/frontier_interfaces_20261018j/manifest.json)
records dependencies and reproduction commands. No additive UNSAT result is promoted to full colouring
UNSAT, and no six-colour forcing claim is made.
