# Independent audit of G105 and normalized phase

Audit status: PASS, scoped to the seven-letter single-bridge language in
`NOTE.md`. No promotion to the full two-port observer or the Chan–Pak
conjecture is justified. Audited by the separate H5-compiler workstream.

## Proof audit

G105's necessity is valid: each letter has nonzero determinant supported
on primes 3,5,7; if a product matrix has every entry divisible by a prime p,
its determinant is divisible by p^3. Thus no outside prime can divide every
entry. This is not an inference from a finite search. Sufficiency follows
from the three displayed integer word identities M(Wp)=p Up: concatenating
the words with exponents a,b,c gives 3^a5^b7^c times an integral matrix.
The Up need not commute; their order is retained in the product. The empty
word correctly handles modulus1.

The note's stronger uniform-H obstruction on actual nonempty inputs is
also valid. The three displayed actual states form a determinant2 basis,
which spans over every odd prime. An invertible word matrix cannot send
all three H outputs to zero, since its H column is the nonzero image of
(1,1,0). At prime2, oddness of tournament H supplies the separate boundary.
This argument concerns uniform operations, not individual count55 or105
realizers.

Gphase's valuation equality and inverse-phase assertion are valid at every
exponent and precision: det Up is a p-adic unit, so Up and its inverse
preserve the common-coordinate p-adic ideal. The exact integer gcd claim
for W5 follows separately from det U5=1. Merely knowing that M(Wp) is
entrywise divisible by p would not imply these conclusions; the explicit
normalized determinants are essential. The q4 q1 hostile correctly has
normalized determinant3 and fails the primitive-unit hypothesis.

The observer scope is explicit and adequate. The word, original names,
and new active port determine a physical tournament under the inherited
all-order bridge theorem. The three coordinates do not encode arbitrary
joint attachments at two old ports. Vanishing raw residues modulo49 or343
therefore supply no proof of an H-neutral reset, and the note does not
claim one. The 81/27 control is correctly identified as automorphism
orders, not modular precision.

## Independent finite checks

The separate companion is `../h5/audit_golden.py`. It reads the seven
integer letter matrices and word lists, then recomputes their products
and Smith data using direct determinant and gcd-of-minor arithmetic.
It does not call the author's multiplication, Smith, RREF, graph search,
or subset-DP state-counting functions. All three normalized determinants
and the unimodular W5 claim agree with `CHECK.json`.

For modular reachability it uses an independent representation: each
subspace is stored as its entire set of vectors over the field, and each
letter sends every vector to its literal image. This avoids RREF entirely.
The complete breadth-first search gives:

| prime | reachable subspaces | dimensions0,1,2,3 | shortest zero word length |
|---:|---:|---|---:|
|3|28|(1,13,13,1)|2|
|5|64|(1,31,31,1)|7|
|7|116|(1,57,57,1)|6|
|11|1|(0,0,0,1)|none|

The first three counts equal the total number 2p^2+2p+4 of subspaces of
F_p^3. This independently verifies the stated exhaustive universe and
shortest letter counts, without claiming a vertex-cost optimum.

The two H55 graphs are recounted by literal permutation paths and
literal subset path-cover sums, independently of the author's subset DP.
Their states are (34,21,61) and (32,23,63); right multiplication by U7
gives normalized residues (3,0,2) and (2,3,5). Both raw W7 outputs vanish
mod7, while their normalized phases remain different, as claimed.

The author's normal and optimized receipts are byte-identical. The
independent audit's normal and optimized receipts are also byte-identical;
source hashes pin the audited NOTE.md, check.py and CHECK.json.

    python3 05-knowledge/results/chanpak_55_105_golden_20261110/h5/audit_golden.py
    python3 -O 05-knowledge/results/chanpak_55_105_golden_20261110/h5/audit_golden.py

Results: `../h5/audit_golden.json` and `../h5/audit_golden.optimized.json`.
No mathematical repair was required. The large 19/22/24-vertex words
continue to depend on the inherited all-order bridge theorem; this audit
does not misrepresent them as independently enumerated Hamilton paths.
