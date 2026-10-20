# Independent audit of the fifteen-probe boundary certificate

Status: **AUDITED / PASS** for the mathematical statements and finite
controls listed here. This audit read [boundary.md](boundary.md) and
[boundary_check.py](boundary_check.py), independently reconstructed the
support cover without importing the producer, and ran the producer's
complete checker successfully. It does not claim an independent
enumeration of every physical pole or audit of the cited external paper.

## 1. Independent cover reconstruction and normalization

I generated the 105 labels by choosing a triple and independently pairing
the four remaining ports, then tested all 2187 ternary words. After the
explicit odd-multiplicity parity filter, the fourteen displayed labels
have this exact support table:

| In target J | Hit by at least one zero probe | Labelled words |
|---|---|---:|
| yes | no | 54 |
| no | yes | 492 |

No other row occurs. The numbers of words hit by 0,1,3,5,9 probes are
respectively 54,384,96,6,6. These sum to 546. All admissible words use
all three colours, so the colour-permutation action is free: division
by six gives 9 target orbits and 82 off-target orbits. In particular the
normalization is 54m for labelled self-closure and 9m for the orbit
inner product. Neither 54 nor the factor six can be dropped.

Nonnegativity makes simultaneous vanishing of the fourteen closure
counts equivalent to support contained in J. This implication uses the
explicit cover, not a rank argument and not cancellation of signed data.

## 2. Physical support forces uniform multiplicity

The Kempe proof is valid. At each cubic internal vertex the subgraph of
two selected colours has degree two; its boundary vertices have degree
one. The component from one terminal of a monochromatic target wire
is consequently a path with exactly one other boundary endpoint.
If that endpoint were not its prescribed partner, switching the two
colours on the component would violate the equality on this wire while
remaining a proper colouring. The global support hypothesis rules this
out.

Switching the component therefore changes precisely the chosen wire's
boundary colour. It is an involutive bijection, fixes the ordered
junction palette and the other wire, and permits every one of the three
wire colours. Applying the same argument to the other wire equates all
nine assignments. Global colour permutations equate the six junction
palettes. Thus the physical extension tensor is mJ with m a nonnegative
integer; the zero tensor is included. Vertex-free open wires do not
invalidate this component argument. Closed wire circles only multiply
all extension counts by three per circle.

The fifteenth measurement fixes m, and hence the FULL tensor. This does
not hold for a general nonnegative tensor merely supported on J: its
nine invariant coordinates may differ. The note correctly retains the
physical hypothesis and separately proves nonnegative uniqueness using
all 105 readouts. That proof's two necessary inequalities, c_i<=0 and
c_i+c_j>=0 on unsupported singleton pairs, force every coefficient to
vanish. The added positive hostile is exact: weight nine on one supported
colour orbit has labelled self-closure 6*9=54 and passes the fourteen
zero probes, but differs from J. Its producer check also passes.

## 3. Independent check of the monomer factor and parity

For a pole with no vertex-free components, fix a boundary singleton
colour at port i. The edges with that colour form a perfect matching
of the internal graph after deleting the incident vertex v_i. Parallel
edges must remain distinct choices. Removing this matching and the
port-i semiedge leaves three paths and some cycles, with degree two at
every internal vertex.

An alternating two-colour path has differently coloured boundary edges
exactly when its number of internal vertices is odd: it has one more
edge than internal vertices. Every admissible path has two colourings;
each even cycle has two, and an odd cycle has none. The factor is

    3 singleton colours * 2^3 path orientations * 2^c cycles
    = 24 * 2^c.

This independently verifies equation (3), its sign and its divisibility
by 24. It also explains the inherited signed identity: if only one of
the three paths has odd internal order, the two even-order paths give
equal numbers of (3,3,1) and (5,1,1) boundary assignments, which cancel.
If all three have odd internal order, all eight assignments contribute
to the nonnegative defect count. The number of odd-order paths is odd,
since the total internal vertex count of a seven-pole is odd.

Scope precision communicated to the author and incorporated: equation (3) excludes
vertex-free OPEN wires and CLOSED wire circles. A separately adjoined
closed circle multiplies N_i by three, not by the residual two-colour
cycle weight two. The producer's matching engine checks simple internal
graphs; its prose correctly specifies edge-copy matchings for the
multigraph extension. No untested parallel-edge census is claimed.

## 4. Remaining exact checks and scope boundaries

The producer's full checker passed in this audit. In particular:

* Independent edge-colouring and matching/path-parity engines agree on
  all 19 specified connected simple five/seven-vertex seven-poles and
  the connected 17-vertex control. The latter has tensor 4J and labelled
  self-closure 216=54*4.
* The aggregate identity has the correct sign and factor six:
  11 C_path630 - 112 C_junction105 = 180 sum_i N_i. Its equality face
  is larger than the individual target family mJ.
* The positive equal-readout hostile is explicitly an ensemble, not a
  silently assumed connected selector gadget. The golden collision uses
  two actual connected poles with distinct defect vectors.
* The golden Gram calculation keeps its field, loop parameter, and
  planar carrier explicit. The block multiplication matrix for a+b*phi
  is correct; specialization phi->4 modulo11 satisfies the defining
  relation and supplies valid nonzero-minor lower bounds. The crossing
  and ternary detectors genuinely refute transport of that null vector.

No minor branch sets, disjoint rooted linkage, bridgelessness, general
Petersen reduction, all-order Gram-rank formula, tournament identity, or
Collatz implication follows from the fifteen measurements. The note
states these limits explicitly. No mathematical correction to the
fifteen-probe certificate or its physical uniformity proof was found.
