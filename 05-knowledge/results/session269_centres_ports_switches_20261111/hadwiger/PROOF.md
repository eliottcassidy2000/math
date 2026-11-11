# Geometric two-ear suspension and the connector-layer support boundary

**Status: PROVED elementary geometric and counting theorems; FINITE-EXACT examples.** No improved bound for the chromatic number of the plane is claimed. No literature-priority claim is made. Here a unit-distance graph is the induced graph on a finite set of distinct points, and every geometric construction below respects that stronger convention.

## Inheritance and concept board

The direct source is Session254, `05-knowledge/results/session254_four_frontiers_20261111/hadwiger_nelson/PROOF.md`, inherited at origin/main `1ca9c5f83b55bcac7dacc0a6f5def38e02f245e7`. It realizes independent four-edge paths with six-colour conditional response `(105^r,104^r)` and asks for cross-path interactions that kill a port state. The hostile is positivity of every finite product. The recovered niche is the actual equilateral-diamond/spindle geometry in `04-computation/hadwiger_dehn_niven_spindle_s685.py`; only its elementary diamond mechanism is used, not its broader interpretive claims. The useful operation lesson from Session268's `automata/PROOF.md` is that a seam must reach the information-bearing layer to break a quotient. Its tournament formula is not transplanted into graph colouring.

Live concepts: exact induced Euclidean gluing; middle-layer interactions; external-field colouring counts; connector-layer state deletion; finite bad-angle avoidance. The source-to-target map here sends a realized core to a realized two-port graph and preserves its full core edges. The missing information in a scalar path count is the core's labelled colour pattern. Retaining that pattern yields both a positive construction and a stronger no-go.

## 1. An arbitrary interacting core can be suspended geometrically

Let W be a finite set of n distinct plane points, and H its induced unit-distance graph. Assume

`U = intersection over w in W of {p : |p-w|<2}`

is nonempty. Equivalently, W lies in some open disk of radius 2. We do NOT rescale unit distances.

**Theorem 1 (exact two-ear suspension).** There are distinct ports A,B and distinct new points x_w,y_w for each w in W such that the entire induced unit-distance graph consists exactly of:

* all edges of H;
* `A-x_w-w-y_w-B` for each w;
* no other edges.

In particular the core's arbitrary existing edges become lawful interactions between the middle vertices of these parallel four-edge paths. The graph has 3n+2 vertices and |E(H)|+4n edges. Indeed almost every pair (A,B) in U squared works with either prescribed choice of the two circle-intersection branches for each ear.

**Proof.** Remove the finite set W from U. For p in the remainder and w in W put v=w-p and choose either sign in

`e_w(p)=(p+w)/2 ± sqrt(1/|v|²-1/4) Jv`,

where J is rotation by 90 degrees. These are distinct real points with distance 1 from p and w. Locally each branch is analytic. We show the forbidden events occupy measure-zero semialgebraic subsets, so their finite union cannot exhaust any nonempty open port domain.

First consider one port p and its ears. A port-core edge puts p on one of finitely many circles. An ear coinciding with another core point z requires |z-w|=1 and |p-z|=1, again a circle condition. If an ear intended for w is also adjacent to z≠w, that ear belongs to the intersection of the two fixed unit circles about w,z. The intersection has at most two points, so p lies on at most two fixed unit circles. Two ears coinciding have the same restriction.

The only less immediate event is an unwanted unit edge between ears for distinct w,z. Because both ears have distance 1 from p, such an edge forms an equilateral triangle. Write the ears as p+u and p+Ru, where |u|=1 and R is rotation by plus or minus 60 degrees. The remaining intended adjacencies say

`|p-(w-u)|=1,   |p-(z-Ru)|=1`.

For each u, these two circles in p have at most two intersections, except when their centres coincide. Centre coincidence is `(R-I)u=z-w`; since R-I is invertible, there is at most one such u, and its possible p form one circle. Thus the whole exceptional set is semialgebraic of dimension at most one: a one-dimensional circle parameter with finite fibres, together with at most one exceptional circle. Its plane measure is zero. This argument applies irrespective of ear branch.

Choose A outside these exceptional sets. Its ears and core are now distinct with exactly the intended one-port incidences. Choose B outside its analogous exceptional sets. Additional unwanted incidences with the finite already chosen A-side set also exclude only measure-zero sets: for an old point t≠w, adjacency of y_w(B) to t places y_w in the at-most-two intersections of the circles about w,t and hence B in finitely many circles. Equality with t puts B on the circle about t (or is impossible). A unit edge from B to any old point is directly a circle condition; equality excludes one point. The case t=w is the intended edge and is not excluded. Avoid A=B as well. The remaining choices form a full-measure subset of U. Setting x_w=e_w(A), y_w=e_w(B) proves the assertion. The same sequential argument and Fubini give the almost-every-pair statement. ∎

This is a genuine induced Euclidean construction, not an abstract graph gadget presumed to have a planar embedding. Edge crossings are allowed; graph planarity is not required by Hadwiger–Nelson.

## 2. Exact response: interacting middle vertices still cannot kill a port state

Write Col_q(H) for the proper labelled q-colourings of H. Fix port colours a,b in {1,...,q}. The exact number of extensions is

`Z_H(a,b)= sum over c in Col_q(H) product over w in W [(q-2+1[c(w)=a])(q-2+1[c(w)=b])].`       (1)

Indeed x_w must avoid exactly the colours of A and w, and similarly y_w avoids those of B and w. These choices are mutually independent after the core colouring is fixed, precisely because Theorem 1 excludes unintended edges.

**Corollary 2 (support saturation).** For every q≥3, every pair (a,b), and every eligible core,

`Z_H(a,b)>0 iff H is q-colourable.`

More quantitatively, writing P_H(q)=|Col_q(H)|,

`(q-2)^(2n) P_H(q) ≤ Z_H(a,b) ≤ (q-1)^(2n) P_H(q)`.

All summands in (1) have positive integer factors for q≥3. This proves both claims. Colour permutation symmetry reduces (1) to equal and unequal port channels, but does not make their counts equal.

**Consequence for Session254.** Adding ANY unit-realizable middle-layer interaction core, not merely a bounded number of cross edges, cannot eliminate unequal-port support while retaining equal-port support at six colours (or any q≥3). The interactions need to touch the connector layer, shorten an attachment, or change the construction in some other substantive way. A non-q-colourable core simply imports an existing obstruction and kills both channels. This is a sharp location-of-interaction obstruction, not an assertion that all geometric interactions fail.

The q≥3 boundary is essential: the explicit square-core example has P_H(2)=2 but both channels zero. At q=2 a connector can run out of colours.

For a one-vertex core, (1) recovers 105 and 104 at q=6, exactly the Session254 primitive. An edgeless n-vertex core recovers their nth powers. Hence the new operation is an actual extension of the old one.

### 2.1 A retained Gram coordinate: exactly when the interaction is invisible

Let q≥3 and put `f_a(c)=product_w(q-2+1[c(w)=a])`. Then equation (1) is the Gram matrix `Z_H(a,b)=sum_c f_a(c)f_b(c)`. By colour symmetry every diagonal entry is D and every off-diagonal entry is O. Therefore

`D-O = (1/2) sum_c (f_a(c)-f_b(c))² ≥ 0` for a≠b.

**Corollary 2a (exact balanced-colour kernel).** For a q-colourable core, D=O if and only if every proper labelled q-colouring of the core has all q colour classes of exactly the same size. In particular q not dividing n forces D>O.

Indeed `f_a(c)=(q-2)^n ((q-1)/(q-2))^N_a(c)`, where N_a counts vertices of colour a. The ratio is strictly greater than 1, so a summand vanishes exactly when N_a=N_b. If D=O this holds for the selected pair in every colouring; colour permutation then gives it for every pair. The converse is immediate. This proof is independent of geometry once the exact induced suspension has been realized.

Thus the core's *colour-class size distribution*, not its unmarked chromatic count alone, controls the bias. Middle interactions can erase or increase a count bias but cannot delete unequal support. The Gram identity gives an exact nonnegative measure of the missing coordinate. At q=3 a unit triangle is a genuine geometric balanced-core control and has D=O=24; at q=6 the square has D>O. An abstract K_q is a useful algebraic control but is not a plane realization for q≥4.

## 3. Positive connector-layer calibration: a continuous family of equality gadgets

A unit diamond has points O=(0,0), P=(sqrt(3)/2,1/2), Q=(sqrt(3)/2,-1/2), T=(sqrt(3),0), with edges OP,OQ,PQ,PT,QT. Its ports O,T have distance sqrt(3). The cross-connector edge PQ forces equal port colours in every proper 3-colouring. Fixed equal ports admit 2 extensions; fixed unequal ports admit zero. At general q its two responses are `(q-1)(q-2)` and `(q-2)(q-3)`. This is inherited elementary spindle forcing, used as a positive control rather than presented as a new lower bound.

**Lemma 3 (generic rooted rotation).** For any finite induced unit-distance point graph containing O, the union of it and a copy rotated about O has, for all but finitely many angles, no coincident vertices other than O and no new cross-copy unit edges, apart from the existing edges incident to O.

For any nonzero original vectors u,v, both `|u-R_theta v|=0` and `|u-R_theta v|=1` have at most two angle solutions: the squared distance is `|u|²+|v|²-2|u||v|cos(theta+phi)`, a nonconstant trigonometric polynomial. There are finitely many pairs. Incidences with O are exactly the intended rooted incidences. ∎

Apply this lemma to the diamond and take as new ports T and R_theta T. The resulting seven-point graph is exactly two diamonds joined at O. Both ports must have the colour of O, so the fixed-port 3-colour response is `(4,0)`. Their distance is

`d(theta)=2 sqrt(3) |sin(theta/2)|`.

Consequently, for every d in (0,2sqrt(3)) except at most finitely many values, a seven-vertex induced unit-distance equality gadget exists at port distance d. In particular such gadgets exist at arbitrarily small positive distances. This supplies a positive, variable-length geometric constraint compiler, not merely count amplification. It is only a three-colour calibration; no five- or six-colour equality gadget is obtained.

An exact small-distance instance is

`cos(theta)=107/108, sin(theta)=sqrt(215)/108`.

It has squared port distance 1/18, seven distinct vertices, exactly ten unit edges, and fixed-port responses `(4,0)` at q=3 and `(48,32)` at q=4. Thus the same geometric gadget also provides a hostile control against silently promoting three-colour equality to four or more colours.

## 4. Exact experimental scope

Run `python3 replay.py` and `python3 equality.py` from this directory (Python + SymPy). The first constructs a square core with ports A=(1/3,2/5), B=(7/5,4/3) and the positive-sign ear formula of Theorem 1. All 91 unordered pairs of the 14 points are classified by exact radical simplification and algebraic sign determination, with no numerical-distance tolerance. Exactly 20 edges occur. The equal/unequal response pairs for q=3,4,5,6 are respectively:

* (162,128)
* (55584,50720)
* (2894292,2772306)
* (57807360,56451584).

The response enumeration covers every q^4 labelled square-core assignment and sums the proved connector choices. Controls include the q=2 failure, the Session254 one-core-vertex response, and abstract K4 non-3-colourability. The K4 control is expressly NOT asserted to have a unit-distance plane realization.

The second script checks all 21 point pairs of the small-distance equality gadget exactly and enumerates all 3^7 and 4^7 labelled colourings subject to the two fixed port choices. Its exact equality-zero and higher-colour positivity controls are saved separately.

## 5. Frontier and cross-lane connection

The constructive gain is a general induced-geometric gluing theorem with arbitrary middle-core interactions. The obstruction gain is its all-core, all-q≥3 support theorem. The positive gain is an exact connector-sensitive variable-distance equality family at q=3, with a four-colour hostile exposing its boundary.

The operation-sensitive connection to Session268 is precise but limited: changes internal to a layer can move counts without exposing the support coordinate; coupling the actual connector layer can expose that coordinate. Source: labelled proper core colourings. Target: port extension counts. Map: equation (1). Preserved predicate: existence of a core q-colouring when q≥3. Lost information under the two-number quotient: individual core colour classes and spatial locations. Necessary sidecar for further construction: the connector adjacency pattern and exact Euclidean incidence. No tournament-to-unit-distance embedding or general chromatic obstruction follows.

Next substantive target: a five-colour gadget where a connector is constrained by at least five already forced distinct neighbour colours, with all those relations realized at exact unit distances. A single unit circle's unit-distance graph is bipartite (unit chords advance by 60 degrees), so a naive clique around that connector cannot supply this palette. The missing relations must therefore be enforced by a larger external geometric network. Theorem 1 does not supply that network. Hadwiger–Nelson remains open.
