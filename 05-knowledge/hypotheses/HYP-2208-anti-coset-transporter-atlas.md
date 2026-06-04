# HYP-2208: Anti-Cosets Are Transporters Across Complement, Converse, and Conjugation Carriers

**Status:** OPEN synthesis, with an exact finite atlas in S632.

## Claim

The anti-coset pattern is not special to self-converse tournaments.  It is the
generic transporter pattern for any carrier with an opposite/converse/
conjugation operation:

```text
Aut(X)  = Iso(X, X)
Anti(X) = Iso(X, JX).
```

If `Anti(X)` is nonempty, then it is a torsor/coset over `Aut(X)` by
composition.  The mathematical work is therefore not to pick a favorite
anti-map, but to name:

1. the carrier `X`;
2. the dual/opposite operation `J`;
3. the automorphism group that acts after one anti-map is chosen;
4. the side channel destroyed when a quotient turns `Anti(X)` into an ordinary
   stabilizer or forgets it entirely.

## Evidence

S632 adds `04-computation/anti_coset_everywhere_s632.py` and stores
`05-knowledge/results/anti_coset_everywhere_s632.out`.

The script checks three finite faces.

### Self-Converse Tournaments

Through `n=7`, every generated self-converse tournament satisfies the coset
law:

```text
|Anti(T)| = |Aut(T)|,
Anti(T) = Aut(T) * tau = tau * Aut(T) for any tau in Anti(T).
```

No coset failures occur.  `H=7` and `H=21` remain absent from the checked
self-converse classes.  The global anti-cycle-type histogram includes long
anti-cycles `(6,)` and `(6,1)`, reinforcing THM-409's warning that the stable
object is the coset/perspective carrier, not a vertex-level involution.

### LRC Pair-Sum Shells

For even `n=6,8,...,20`, the action of `<2,-1>` on folded shells modulo
`C=2n-1` has exact transporters from a shell to its doubled shell.  In every
checked row, the transporter size equals the stabilizer size.

At `n=14`, `C=27`, the thirteen folded shells collapse to three orbits:

```text
gcd 1: [1,2,4,5,7,8,10,11,13]
gcd 3: [3,6,12]
gcd 9: [9]
```

This recovers the THM-407/S599g prime-3 shell split as an anti-coset-style
transporter statement.  The important caveat is that folding `j ~ -j` turns
reflection anti-data into a stabilizer; the quotient has erased the orientation
side channel.

### Eisenstein Unit-Distance Prefixes

For Eisenstein/triangular prefixes under the finite `D6` lattice action, full
centered hex shells keep the full orientation-reversing anti-coset:

```text
n=1:  Aut+ = 6, Anti = 6
n=7:  Aut+ = 6, Anti = 6
n=19: Aut+ = 6, Anti = 6
```

Partial shell prefixes beyond `n=7` usually keep only one reflection or none.
The script's prefix `n=21` has `47` unit edges and a single anti-map, not the
full centered-hex coset.  This is a finite symmetry version of the warning that
the first complete Eisenstein shell does not control the later exact/Moser
`n=21` carrier, where HYP-2206 instead records `57 = 20 + 37`.

## Tournament Analysis

S632 uses proof/carrier lenses as vertices:

- `SC Anti(T)=Iso(T,T^op)`;
- LRC shell transporter;
- Eisenstein `D6` conjugation;
- coimage sign-reversing involution;
- raw complement merge;
- raw scalar `7/21` match.

The pairwise observable is which lens retains transporter, side-channel, and
quotient-loss data.  The gauges are SC torsor theorem, LRC shell action, and
Eisenstein `D6` conjugation.  The majority tournament is transitive with
Hamiltonian-path count `1`, ranking the explicit SC anti-coset first and raw
scalar matching last.

## Assumption Challenge

Alternate vertices considered: tournament vertices, anti-maps, rooted
perspectives, residue shells, folded shell orbits, unit-distance points, `D6`
symmetries, coimage cancellations, LRC proof obligations, and raw scalar
coincidences.

Chosen vertices: carrier lenses.  They preserve the predicate "`Anti(X)` is a
transporter to `JX`, and quotienting can hide it."  They destroy exact
embeddings, complete LRC certificates, full higher-`n` tournament spectra, and
raw equality among visible integers.

The challenged assumption is that complement/converse/conjugation should be
studied by selecting one canonical anti-map.  The coset is the object; a chosen
map is a coordinate.

## Next Tests

1. Replace raw complement merges in older tournament scripts by explicit
   `(X,J,Aut,Anti)` records.
2. Recompute the LRC `64` fixed-boundary classes as transporters before and
   after converse merging, recording which labels are lost.
3. Audit exact `n=21` unit-distance/Moser cores under orientation-preserving
   and orientation-reversing carrier symmetries, not only edge counts.
4. Search coimage/sign-reversing-involution scripts for hidden transporter
   sets: what is `X`, what is `JX`, and what stabilizer acts on choices?
