# The same two norm-seven ideals give density without chromatic forcing

**FINITE-EXACT construction, exact chromatic number3; recovered earlier work.**
This is an independent replay of the mechanism behind MISTAKE-575 and
`04-computation/experiments/oai2_20261006_lattice_u21.py`, not a new extremal
unit-distance theorem.

Use omega²+omega+1=0, alpha=3+omega, beta=2−omega, both norm7. Let B be
the center plus the six Eisenstein units and P={0,1,1+omega}. The sumset

    S = alpha P + beta B

has exactly21 distinct lattice points and57 pairs at squared distance7.
The entire graph has a proper coloring a+b modulo3 and contains triangles,
so its chromatic number is exactly3. Its degree multiset is 5^18,8^3.
Dividing all coordinates by alpha turns every displayed edge into a unit
edge while retaining all points in Q(√−3).

The two factors are precisely the conjugate norm-seven ideals used by the
seven-color torus construction. Therefore the SAME arithmetic carrier can
serve a dense distance graph, a perfect code and a seven-color quotient,
without giving the source graph chromatic number6 or7.
The earlier false 'no single lattice' inference is not revived.

The repository cites Alexeev–Mixon–Parshall for global extremality at21;
this control does not independently prove or need that external extremal
claim. Its verified consequence is the exact57-edge/three-color witness.

Run `python3 check_dense_flower.py` and `python3 -O check_dense_flower.py`.
The finite universe is every unordered pair and triple of the explicit21
points. No bounding-box search or numerical distance comparison is used.

## Irrational rotation does not mean lattice escape

For O=Z[omega] and K=Q(omega), a complex unit z has zO commensurate with O
if and only if z lies in K. If commensurate, some positive multiple of z
belongs to O. Conversely clear rational coefficients of z and z^−1 in
the basis(1,omega) to obtain finite-index containments. Thus the rotation
commensurator is the infinite group K intersect S1; the exact rotational
stabilizer of O is only its six units.

Here alpha/beta=(11+5sqrt(−3))/14 has infinite order: if its angle were a
rational multiple of pi, twice its real part11/7 would be a rational
algebraic integer, which it is not. Nevertheless alpha O and beta O are
commensurate index-seven lattices; as distinct prime ideals their sum is O.
The associated finite pattern remains three-colorable. The step from
finite D12 symmetry to an irrational-angle rotation has NOT escaped the
quadratic field or its coloring bound. This repairs an inference, not an
assertion that every irrational rotation is harmless.
