# What576 stores, and how easily the prescribed input leaks

**PROVED elementary** fibre decoder and fresh-colour obstruction.
**FINITE-EXACT** one-port perturbations. No new six-chromatic graph.

## The correct arithmetic and the actual fibre

625+125+25+1=776, not576. Useful exact identities are
```
576=24^2=(25-1)^2=625-50+1=625-49,     7^2+24^2=25^2.
```
The24 in our [composed gate](../first_exit_complement_games_20261020/PALETTE-GLUE.md)
is an actual extension count, not initially a Euclid parameter.
In the inherited16-point gate, the three monochromatic input triples
use a palette A of three colours. Every noncentral interior vertex sees
all three colours A, and must use its complementary three-colour palette B.
For an output colour in A, the remaining six-vertex graph consists of
two triangles joined by one edge. Its colourings have the explicit decoder
```
bridge endpoint:3 choices; opposite endpoint:2 choices;
the two remaining triangle vertices:2 orders on each side.
```
Thus there are3*2*2*2=24 extensions per allowed output. With an output in B,
all seven core vertices would use B, contradicting the Moser spindle's
three-colour obstruction. This recovers the old gate relation structurally.

Two copies sharing only the output have independently chosen internal
extensions, so the count per shared allowed colour is24^2=576. This is the
same product mechanism as the incoming q-by-q Hamilton array, but its
consumer differs: here matching named output colours, there matching full
subset choices under an anti-involution. No graph-minor claim transfers.

There is also a parity obstruction to identifying the two symmetry models.
The incoming pointed-set involution has one fixed point on an odd set of
size q. An involution on the24-element fibre cannot have exactly one fixed
point, since its nonfixed points occur in pairs. The numerical triple7,24,25
does not supply that missing action. In the incoming orbit construction,
24 instead occurs as (7^2-1)/2, with Euclid parameters4,3. The binary family
q=2^(k-1)+1 itself does not contain q7; that is a general odd-set comparison,
not a member of the claimed tournament family.

## Every proper one-port change escapes this bad input

Take the31-point composition, with left input palette012 and right345.
The prescribed full boundary is proper and has no extension. Perturb
exactly one of the nine left ports to one of the other five colours.
There are45 changes. Four violate an existing port edge; the other41 are
proper, and **every one of those41 admits extensions of the whole graph**.

The checker computes all six output counts on the actual16-point graph.
The right unmodified copy then contributes24 per output in345, so the
complete extension number is24 times the sum of those three output counts.
For example recolouring port7 from0 to1 changes the left row to
(60,84,84,8,8,8), giving576 extensions of the full composition.
Changing it from0 to3 gives (32,48,48,8,4,4), giving384.

This is a finite exact sensitivity certificate at the displayed boundary.
It does not classify all perturbations or prove a universal robustness
law. It does show that this proposed synchronizer must rule out41 adjacent
proper boundary states, not just the one scalar palette intersection.

## Why literal input equality is an expensive target

If a graph is k-colourable and q>=k+s, any selected s vertices can be
recoloured with s distinct unused colours, preserving all edges. Thus in
a four-colourable carrier at palette6, any two ports can be assigned two
fresh distinct colours. More generally, forcing equality of two distinct
ports in every proper q-colouring of a colourable graph requires chromatic
number at least q: a (q-1)-colouring could recolour one port freshly.

Consequently a nonvacuous context forcing even one of the input equalities
at palette6 already needs a six-chromatic graph. Forcing that entire bad
input and obtaining no six-colouring would seek an even stronger bound.
The known four-colour field carrier cannot supply such a context, however
many copies are glued within it. A full relation can still be useful, but
externally prescribed equality must not be treated as a free constructor.

The constructive reroute is the incoming five-colour circle-saturation
target: retain every possible boundary partition, and seek a geometric
fragment whose every five-colouring uses all five colours on neighbours
of a new centre. Adding the centre would require a sixth colour. This
is an OPEN target, not a construction supplied here. Complete alternative
input patterns, a carrier escape and actual unit contacts are essential.

Our receipt checks all45 changes and36 fresh-two-port recolourings.
It imports only the pinned exact edge/colour receipt; no geometric
identification or extra unit edge is silently introduced.

