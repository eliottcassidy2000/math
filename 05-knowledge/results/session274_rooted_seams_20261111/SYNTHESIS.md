# Rooted seams are a recursive coordinate system, not just a count operation

**Session274, 11 November 2026. Chan–Pak remains OPEN.**

The fixed Hamiltonian path is useful boundary data, not a cosmetic drawing.
Once the seam root is kept in that path, physical vertex births, tiling
coordinates, metagraph lines and path-count responses fit into one exact
framework. The session supplies all-order identities and structural tests,
not a larger numerical coverage range.

## 1. A reusable physical seam with a complete three-experiment decoder

Add a new vertex b which beats the root a and loses to every other old
vertex. Insert b immediately before a in the fixed Hamilton path. This always
works, including either endpoint. Repetition produces a transitive fan before
a and never leaves the fixed-path chart.

The whole path-cover polynomial of the k-fold fan is

    Q_k(t)=A(t)(t+1)^k+B(t)t^k+C(t)(t−1)^k.

The coefficients are explicit rooted-colour boundary responses, not fitted
parameters. Therefore a cubic recurrence generates every fan size; observing
the full profiles at k=0,1,2 determines the entire infinite continuation
sequence. This is a precise finite description of an infinite construction.

At t=1, the zero mode disappears after the first birth. This explains why
scalar Hamiltonian counts lose information even when every future scalar fan
count is known. The seam renews indefinitely, but its scalar excess doubles
rigidly. **Renewal is not independent target control.**

[Physical construction, proof and initial cases](constructions/PROOF.md)

## 2. Every framed tournament has a unique seam-and-gauge growth history

Keep the path arcs fixed and view all other arcs as tile bits. Their graph is
G_n=K_n minus the path. The inherited cut-space theory classifies gauge orbits
by binary cycle charges. A cut-seam flips the free rectangle across a rooted
path cut; these moves generate the cut space.

There is a stronger recursive bridge: a rooted birth, followed by a child
cut gauge whose labels at the old root and newborn are fixed to zero,
is a **bijection** from

    parent tiling × {0,1}^{n−1}  ->  child tiling.

It works at every size and at every root position. An explicit reverse map
reads the new chord bits, undoes the gauge, and deletes the newborn. Repeating
first-root births from a singleton gives a lossless recursive description
of every framed tournament. This recovers the known tiling count rather than
proving a new A000568 formula: different frames and ordinary isomorphism
classes are still not identified.

For n≥4, one parent seam sector births into every compatible child cycle-charge
sector, with exactly two birth images in each. Thus new geometric choices
are supplied at every step, rather than drawn from a finite stock.

[Exact geometry and reverse maps](geometry/PROOF.md)

## 3. The metagraph connection depends on where the root sits

For parent order n≥4, a metagraph complement line flips ALL free tiles, not one arc. Birth at the
first path vertex intertwines complement lines modulo cut gauges. Birth at
an interior or final path vertex has a nonzero transport defect: the old path
edge from predecessor to root becomes a free chord, creating a detectable
cycle charge. This is why a seam that is perfectly legal as a tournament
operation need not behave identically as a metagraph operation.

The obstruction has an explicit cycle witness, not an appeal to a visual
analogy. Root position is the missing coordinate.

## 4. Forest observables cannot see cycle charge; correlated paths can

Every seam sector has exactly the same mean **complete path-cover profile**.
A path cover uses a forest of free edges; any required orientations on a
forest impose independent cut-switch equations. The cycle charge therefore
cannot affect the mean. This is an exact all-order conservation mechanism.

For higher moments, unions of several paths can contain cycles. The q-th
moment is an exact sum over compatible q-tuples of paths, with compatibility
determined by the sector's cycle parities. This gives a concrete hierarchy
of observables that progressively retains information discarded by averaging.

At order5, the four sectors all have mean H=79/8 but four distinct variances.
At order6, variance gives only19 values for32 sectors. The latter is an
explicit warning against promoting the first successful decoder into a
universal one.

The common first moment also has an explicit three-state run automaton and
factorial readout. This connects the old finite-state spine machinery and
factorial-Hankel perspective without pretending either settles Chan–Pak.

[Moment formula, recurrence and controlled examples](moments/PROOF.md)

## 5. What this changes for a constructive Chan–Pak strategy

Session273 transported hidden resources while preserving full unmarked
profiles, then exposed them through one seam. This session separates three
questions that were previously easy to conflate:

1. **Can the interface be renewed?** Yes, the named seam and frame survive
   arbitrary fan length.
2. **Can recursive growth reach all geometric configurations?** Yes, birth
   plus the specified cut gauge has an explicit inverse and spans all tilings.
3. **Can these choices be scheduled to realize each prescribed odd count?**
   Still OPEN. Neither uniform means, all charge sectors, nor a universal
   geometric grammar supplies the required integer selector.

The next target is consequently more precise: carry the rooted response
through the gauge stage of the bijection, and find a target-bearing invariant
whose admissible inverse step decreases a prescribed odd count. A fixed
root's three-mode fan packet is sufficient only before arbitrary gauges;
it must not be assumed closed under them. The moment collision at order6
also rules out replacing that obligation by variance alone.

Names for the reusable principles:

- **Frame-anchored interface:** retain a witnessing path so operations have
  canonical coordinates and explicit inverse maps.
- **Forest blindness / loop sensitivity:** one-path linear observations
  forget cycle charge; correlated histories can detect it.
- **Geometric completeness versus arithmetic selection:** a universal
  construction grammar still needs a count-controlling scheduler.

These are descriptions of the proved mechanisms above, not additional grand
theorems. No proof of Chan–Pak, eventual threshold, or complete ordinary
isomorphism classification is claimed.

## Reproduction and scope

See [REPRODUCE.md](REPRODUCE.md), [review](review/REVIEW.md), and
[open obligations](OPEN_QUESTIONS.md). Historical claims were used as scoped
inputs; the independent review concerns the new bridge, not a re-audit of
the repository's old theorem inventory.
