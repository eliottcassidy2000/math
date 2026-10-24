# Validation boundaries and checks

This is an additive structural packet. It does not modify the independently
frozen nonregular-order-eleven audit or its producer.

Inputs are the source-pinned elementary tournament helper implementation and
the source-pinned nonregular/SCC audit constants9166 and19480. No almost-regular
ten-vertex atlas, parent count, regular-eleven atlas, regular-eleven class count,
or A000568(11) value is read or used for calibration.

The exhaustive mathematical coverage comes from full automorphism orbit
partitions, with explicit elimination of all but11,5+5+1,3+3+1^5,3+1^8,1^11.
Candidate enumeration is deliberately small:

- Forty noncoordinate-line-invariant degree-three5-by5 matrices, quotiented
  into eight C5^2 orbits. Full group checks reject one hidden vertex-transitive
  Aut55 case as an error22 class.
- Five q7 presentations: one uniform crossing, and two C5 attachment orbits
  times two matching signs for a nonuniform crossing.
- Two q9 presentations: the perfect matchings between the named degree-one
  vertices of the low strong4 and degree-two vertices of the high strong4.
- Thirty-two circulant connection sets, quotiented into four multiplier orbits.
  Full pairwise isomorphism tests independently certify that none merge under
  a more general eleven-vertex relabeling.
- The elementary6-fixed-vertex exclusion is additionally checked on all32768
  six-vertex tournaments; all64 four-vertex tournaments verify the low/high
  internal block classification. Neither is an eleven-vertex census.

Nineteen full automorphism computations and ordinary pairwise isomorphism
checks leave eighteen nonrigid regular-eleven classes. One candidate, Aut55,
appears in the q3 presentation family and in the circulant family and is
explicitly identified and counted once. All153 pairs of the final eighteen
records are tested for full ordinary isomorphism. All regularity, strongness,
full orbit-size, and error claims are checked directly on literal rows.

The normal and optimized outputs are byte-identical. Six hostile changes,
each run normally and with-O, are rejected: a missing stabilizer slope,
missing negative matching sign, missing second q9 matching, truncated cyclic
universe, double-counted hidden Aut55 control, and confusion between pointed
orbits and labelled vertices. All twelve mutated scripts also fail the manifest
byte check. The final replay copies only manifest-listed files and repeats the
checks from an unrelated working directory.

The source-pinned helper contains additional unused main-audit routines; this
packet never invokes them or reads their data inputs. Its input contract is
therefore exactly the one stated above, not the larger main-audit universe.
