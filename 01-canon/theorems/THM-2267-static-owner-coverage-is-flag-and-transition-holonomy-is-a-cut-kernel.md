---
id: THM-2267
title: "Static owner coverage is flag and transition holonomy is a cut kernel"
status: >
  PROVED + VERIFIED-EXACT. For any finite family of measurable owner-service
  sets, the Boolean redundancy game has Hasse increments equal to overlap
  with the accumulated service union, signed dividends equal to alternating
  multiple intersections, and nonnegative coalition-cut sums. Its
  zero-redundancy complex is always flag: every minimal nonface is already a
  positive pair overlap. Thus replacing its owner pair graph by this full
  zero complex adds no obstruction. The first genuinely new finite coordinate
  is the exact gluing of obligations across chambers. On any
  finite transition graph, minimum owner-switch energy dominates every
  binary owner-cut min-cut; for two owners this is equality, and on a unit
  cycle every positive cut holonomy is at least two. A four-chamber
  identity-versus-swap pair has identical complete static service/Hasse data
  but switch energies zero and eight. This is a transition-carrier theorem
  and a static-data no-go, not an exclusion of an LRC(14) profile.
source: codex-2026-07-25-owner-transition-holonomy
depends_on: []
related:
  - THM-2211-carry-regime-root-transducer-and-infinite-autonomous-index
  - THM-2221-tournament-context-cut-metric-and-pinned-transport-response
  - THM-2246-depth-one-private-joint-two-step-fibre-cap
  - THM-2255-valuation-separated-pair-cap-and-exclusive-owner-mass
  - THM-2256-automorphism-contact-dichotomy-for-quotient-frustration
  - THM-2259-boolean-continuation-hasse-field-and-signed-interaction-dividends
script: 04-computation/owner_coverage_hasse_transition_cut_thm2267.py
output: 05-knowledge/results/owner_coverage_hasse_transition_cut_thm2267.out
script_sha256: f659e812379367b6013c930e7db821da86542c9600012328e17dc1b50723ce11
output_sha256: 8a3edef1c0c5c58e1e93a3efca8f2996e8e56a6726b9865491518b908550c2ef
hash_basis: working-tree bytes (LF)
---

# THM-2267 -- static owner data is flag; gluing carries the holonomy

THM-2259 identifies a nonnegative flat Hasse field and signed Boolean
dividends as the exact finite diagonal interaction carrier in a general
nonexpansive metric monoid.  The literal union-redundancy game of owner
coverage is more rigid than that general setting: its zero-defect complex can
never have a genuinely higher minimal nonface.  THM-2256 gives a
complementary warning about forcing the remaining pair data into a
tournament.

The missing coordinate is not a more elaborate static complex.  It is the
matching of partial obligations across phase, root, or carry transitions.
Once that matching is retained, its first intrinsic quotient is a symmetric
cut kernel and its first obstruction is owner-switch holonomy.

## 1. The static owner game

Let `(Omega,mu)` be a finite measure space, let `I` be a finite set of owner
labels, and let

```text
B_i subset Omega,                         i in I,     (1)
```

be measurable service sets.  For `S subset I`, put

```text
U_S=union_(i in S)B_i,

Delta(S)=sum_(i in S)mu(B_i)-mu(U_S).                 (2)
```

The quantity `Delta(S)` is the total service redundancy.  Pointwise, if

```text
m_S(x)=#{i in S:x in B_i},
```

then

```text
Delta(S)=integral (m_S-1)_(+) dmu.                   (3)
```

In particular `Delta` is nonnegative, vanishes on the empty set and on
singletons, and is monotone.

For `i notin S`, its Boolean Hasse increment is exactly

```text
c(i|S)
 =Delta(S union {i})-Delta(S)
 =mu(B_i intersection U_S)
 >=0.                                               (4)
```

Every Boolean square is flat because `c` is the gradient of `Delta`.
For disjoint `A,C subset I`,

```text
Delta(A union C)-Delta(A)-Delta(C)
 =mu(U_A intersection U_C)
 >=0.                                               (5)
```

Thus the general THM-2259 coalition-cut inequality becomes literal overlap
of the two accumulated service unions.  Equivalently,

```text
sum_(
  T subset A union C,
  T intersection A nonempty,
  T intersection C nonempty
) h(T)
=mu(U_A intersection U_C)>=0.                       (5a)
```

Let `h(T)` be the Boolean Möbius dividend of `Delta`.  Inclusion--exclusion
in (2) gives the complete signed spectrum:

```text
h(T)=0,                                      |T|<=1,

h(T)=(-1)^|T| mu(intersection_(i in T)B_i),   |T|>=2. (6)
```

Consequently the higher dividends may be negative, but their signs and
magnitudes are not free interaction data: they are the alternating
multiple-intersection ledger of one set system.

### The flag theorem

Define the zero-redundancy complex

```text
K={S subset I:Delta(S)=0}.                           (7)
```

Equation (3) proves

```text
S in K
 iff mu(B_i intersection B_j)=0
     for every distinct i,j in S.                   (8)
```

Indeed, pairwise null intersections make the finite union essentially
disjoint, while any positive pair intersection makes the integrand in (3)
positive there.  Hence `K` is the clique complex of the graph whose edges
are the null-overlap owner pairs.  Equivalently:

> Every minimal nonface of a static owner-coverage redundancy complex has
> cardinality two.

This is sharper than the general interaction hierarchy in THM-2259.  Higher
minimal nonfaces are possible for metric composition, but impossible for
literal union coverage.  Passing from the owner-overlap graph to its full
zero complex adds no obstruction.

## 2. Exact transition energy and its binary cut quotients

Now let `G=(V,E)` be a finite undirected graph with strictly positive edge
weights `w_e`.  Its vertices are partial obligations and its edges are
**actual continuation correspondences**: the same sheet, carry stalk, or
labelled piece before and after one transition.  Let `L` be a finite owner
set, and give every vertex a nonempty eligibility set

```text
A_v subset L.                                       (9)
```

Define the minimum owner-switch energy

```text
omega_G(A)
 =min_(ell(v) in A_v)
    sum_(uv in E) w_(uv) 1_[ell(u)!=ell(v)].        (10)
```

This is intrinsic and symmetric.  No orientation of an owner pair has been
introduced.

Because all edge weights are positive,

```text
omega_G(A)=0
 iff intersection_(v in C)A_v is nonempty
     for every connected component C of G.          (11)
```

The forward implication follows because a zero-energy labelling is constant
on each component; the reverse implication chooses one common owner on each
component.

For every nonempty proper owner cut `Q subset L`, project (9) to the two
sides of the cut:

```text
A_v^Q contains 0 iff A_v intersection Q is nonempty,
A_v^Q contains 1 iff A_v intersection (L minus Q) is nonempty. (12)
```

Let `kappa_Q` be the binary version of (10) for these projected eligibility
sets.  Every legal owner labelling projects to a legal binary labelling, and
an owner-side change can occur only on an owner-change edge.  Therefore

```text
omega_G(A)>=kappa_Q

and hence

omega_G(A)>=max_(empty!=Q proper subset L) kappa_Q. (13)
```

The right side is an exact, efficiently computable cut obstruction.  Put

```text
V_0={v:A_v^Q={0}},             V_1={v:A_v^Q={1}}.    (14)
```

Then `kappa_Q` is the minimum weight of an edge cut separating `V_0` from
`V_1`; vertices with eligibility `{0,1}` may lie on either side.  This is a
standard source--sink min-cut after attaching the forced sets by
infinite-capacity edges.  If `|L|=2`, projection loses no label information,
so (13) is equality.

For a unit-weight cycle, every cyclic label word has either zero changes or
at least two.  Thus

```text
omega_G(A)=0 or omega_G(A)>=2.                       (15)
```

For a binary owner cut, delete the flexible `{0,1}` vertices and read the
remaining forced bits cyclically.  If no forced bit occurs, or all forced
bits agree, then `kappa_Q=0`.  Otherwise

```text
kappa_Q
 =#{cyclic alternations in the remaining forced word},    (16)
```

an even integer at least two.  Each unequal consecutive pair of forced bits
needs a switch in the intervening arc, and placing exactly one switch in
each such arc realizes (16).

### Min-plus transition kernel

For a unit-weight path with vertices `v_0,...,v_(m-1)`, let

```text
D(i,j)=0 if i=j, and D(i,j)=1 otherwise,

Q_A(i,j)=0 if i=j in A, and infinity otherwise.      (17)
```

Use tropical matrix multiplication.  The exact minimum switch cost from
initial owner `i` to terminal owner `j` is the `(i,j)` entry of

```text
Q_(A_v0) tensor D tensor Q_(A_v1) tensor ...
             tensor D tensor Q_(A_v(m-1)).           (18)
```

On a cycle, append one final `D` and take the tropical trace.  This recovers
`omega_G` exactly.

The diagonal projectors in (18) are the static eligibility data.  The
off-diagonal switch kernel and, crucially, the order in which the projectors
are glued are additional coordinates.  This is the finite transition
analogue of THM-2259's distinction between diagonal Hasse stalks and a full
min-plus continuation kernel.

## 3. Exact static-data collision

The need for gluing appears in the smallest useful hostile model.  Take four
cyclic chambers `r in Z/4Z`, two chamber-local obligations

```text
Omega_r={a_r,b_r},
```

and two owners.  In every chamber set

```text
B_(r,0)={a_r},                 B_(r,1)={b_r}.        (19)
```

Thus the complete service incidence, every value of `Delta`, every Hasse
increment, every Möbius dividend, and the zero complex are identical in all
four chambers.  In fact `Delta` is identically zero.

There are two gluings:

```text
identity: a_r -> a_(r+1),       b_r -> b_(r+1),

swap:     a_r -> b_(r+1),       b_r -> a_(r+1).      (20)
```

Include the closing transition from chamber three to chamber zero.  Under
identity gluing, each continuation component has one constant eligible
owner, so

```text
omega_identity=0.                                   (21)
```

Under swap gluing, all eight transition edges join obligations with opposite
singleton owner sets.  Every edge is forced to switch, so

```text
omega_swap=8.                                       (22)
```

The two systems have the same complete per-chamber static owner data but
different transition energy.  Therefore no invariant factoring through the
unordered collection of static service/Hasse stalks can determine switch
holonomy.  The missing information is exactly the chamber-to-chamber
obligation correspondence.

## 4. Why a static tournament is not the repair

The pair observable in (4)--(5) is symmetric, weighted, and may vanish.  It
does not orient owner pairs.  There is also an equivariance obstruction to a
universal forced orientation.  If two owner labels have identical service
sets and identical sidecars, their transposition is an automorphism of the
input data.  An intrinsic construction must preserve that automorphism, but
no tournament admits the transposition of the endpoints of one arc as an
automorphism.  Hence a tournament cannot be assigned equivariantly on all
static owner systems without adding a genuine asymmetry.

There is a second LRC-specific boundary.  In the strict valuation profiles,
orienting blockers by increasing `13`-adic depth produces a transitive
tournament.  THM-2256 proves that its forced quotient-frustration floor is
exactly one at every blow-up size, so that bare quotient envelope has no
macroscopic scale.  In the repeated-first-depth profiles, depth itself leaves
a tie and does not define a tournament at all.

This does not prohibit a tournament on richer transition states.  It says
that the vertices would have to be owner-carry obligations with an intrinsic
asymmetric transition observable.  The static owner labels are insufficient.

## 5. LRC(14) instantiation and honest remaining implication

At a fixed `13`-adic scale and a phase chamber, let `Omega_v` be the finite
root-sheet obligations left by the guard and the five unit masks.  For a
blocker label `j`, let

```text
B_(v,j)={obligations serviced by blocker j}.         (23)
```

The scalar cover hypothesis says that these service sets cover
`Omega_v`.  Section 1 applies verbatim, with counting weights on a torsion
fibre or cell/Haar weights on a finite chamber partition.

THM-2255 proves that every one of the remaining `165` profiles has a
positive exclusive-owner stratum.  In the `150` strict profiles, at least
one labelled stratum expands through its owner's expiration to mass

```text
88159/415800>1/7,                                   (24)
```

so it cannot remain inside one successor danger comb.  Before expiration,
its owner and nonowner bits transport exactly.

Equations (10)--(16) identify a lawful consumer for this information, but do
not manufacture its missing hypothesis.  To obtain a positive cut tax one
must still build:

```text
vertices: exact transported pieces of the exclusive flow;
edges:    their common sheet/carry correspondences through expiration;
A_v:      the named successor owners or absorbers available to each piece.
                                                               (25)
```

Only then may one exhibit a cut `Q` with `kappa_Q>0` and convert its edge
weights into a measure loss.  The mass expansion (24) forces a split or
handoff, but does not by itself identify the target service sets in (25).
Therefore THM-2267 excludes none of the `165` profiles.

THM-2246's hostile `112`-root fibre gives the sharp zero control.  One deep
blocker services every obligation.  Any transition graph confined to that
persistent-owner stalk has a common eligible owner, so (11) gives
`omega=0`; on the retained one-owner label set its static zero complex is a
simplex.  Neither that redundancy complex nor a raw winding number can turn
this specific local carrier into a positive holonomy.

## 6. Typed connection and loss ledger

The useful map is

```text
source:
  multi-chamber labelled root sheets with exact carry continuation;

target:
  a weighted eligibility graph (G,(A_v)) and its cut energies kappa_Q;

map:
  make every transported obligation piece a vertex, join exactly
  corresponding pieces, and retain its complete set of available owners;

preserved:
  legal owner sections, every forced owner switch, componentwise persistent
  owners, binary owner-cut obstructions, and tropical composition;

destroyed:
  phase distance inside a chamber, mask endpoint geometry, service
  multiplicity unless put in edge weights, and all data not encoded in the
  chosen owner labels;

mandatory sidecars:
  the exact sheet/carry gluing, births and deaths at endpoints, and explicit
  absorber labels when an obligation leaves the blocker residual;

cheapest decisive test:
  on an actual THM-2255 exclusive flow, compute one binary min-cut across
  owner expiration and ask whether it has positive Haar-weight scale.      (26)
```

For knots, THM-2259's Boolean Hasse stalk similarly retains diagonal subset
responses while the full continuation kernel retains arbitrary
intermediates.  Formula (18) is not a Gordian theorem, but it makes the
structural correspondence exact: static potentials are diagonal projectors,
and transition gluing lives in the off-diagonal min-plus product.  The
identity-versus-swap collision is the finite proof that diagonal stalks do
not recover that product.

## 7. Independent exact audit

The companion exhausts all `512` three-owner service systems on three
obligations.  It checks (4)--(8), all Boolean squares, all Möbius dividends,
and all disjoint coalition cuts.  It also exhausts all `2,401` nonempty
three-label eligibility systems on a four-cycle, checking (13), (15), and
the binary formula (16), and verifies (21)--(22).

Reproduce with

```bash
python3 04-computation/owner_coverage_hasse_transition_cut_thm2267.py
python3 -O 04-computation/owner_coverage_hasse_transition_cut_thm2267.py
```

Normal and optimized transcripts are byte-identical to the stored output.
QED.
