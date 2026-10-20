# Nodal-cylinder period obstruction

**Status:** proved; independently reviewed.

The exact generating-pair identity
\[
\int_{-1}^{1}P Q_y\,dy=\frac8{15}\det L
\]
is constant for
\(\mathbb C[x,P,Q]=\mathbb C[x]+(y^2-1)\mathbb C[x,y]\).
Consequently this nodal cylinder admits no generating Keller pair.

Combined with the earlier actual-source conductor classification, this
excludes generic conductor degree two. Under that packet's exact hypotheses,
a noninvertible Keller map must instead have generic conductor degree at
least four. This is a conductor-degree bound, not a sheet-number bound, and
does not solve the complex Jacobian conjecture.

- [Proof, classical input, and hostile controls](NODAL_CYLINDER_PERIOD_OBSTRUCTION.md)
- [Source verification and dependency record](SOURCE_VERIFICATION.md)
- [Exact control script](verify_nodal_period.py)
- [Finite results](results.json)
- [Completed independent audit](review/AUDIT.md)

Reproduce with Python 3 and SymPy:

```sh
python3 verify_nodal_period.py --output normal.json > normal.out
python3 -O verify_nodal_period.py --output results.json > optimized.out
cmp normal.json results.json
cmp normal.out optimized.out
```

Producer result: **285 exact checks across 17 explicitly declared generating
families**. Independent result: **491 checks across 36 families**. Normal and
optimized execution are byte-identical for both. These are finite controls,
not a computational proof of the imported all-embedding theorem.

No relative ambient automorphism over \(\mathbb C[x]\) is assumed. The proof
uses fiberwise ambient equivalence plus a global constant conormal determinant.
