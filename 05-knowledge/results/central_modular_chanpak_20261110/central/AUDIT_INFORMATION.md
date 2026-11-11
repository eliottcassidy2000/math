# Independent audit of the golden content classification

**PROVED scoped audit; FINITE-EXACT independent products.** This audit reads
[information/PROOF.md](../information/PROOF.md), the inherited
[golden RECURRENCE](../../golden_hamilton_transfer_20261103/RECURRENCE.md),
and the inherited
[G105/Gphase note](../../chanpak_55_105_golden_20261110/golden/NOTE.md).
It imports no executable code from either the current information lane or
the older golden experiments. It does not independently re-prove the entire
inherited tournament transfer theorem.

## Verdict and completed quantifier correction

The mathematically correct classification is:

> For an integer r>=1, there exists a golden word w satisfying
> cont(vM(w))=r cont(v) for every integer row v if and only if
> r=5^b 9^k 21^c for nonnegative integers b,k,c.

The proof establishes this statement, and the live theorem now states it
with the required existential quantifier. The originally audited sentence
“A golden word w has M(w)=rU ... if and only if r belongs to the monoid”
has been replaced by the existential wording above. For a fixed arbitrary w
the converse is false: w=s cannot produce scalar5, even though5 belongs to
the monoid. This is a quantifier correction, not a failure of the construction
or determinant obstruction. The finding was sent to the author and parent;
this audit retains it as the correction lineage after the live text's repair.

## The content-to-unimodularity implication is sound

For nonsingular integral M and integer r>=1, testing standard basis rows
shows every row has content r, hence every entry is divisible by r. Therefore
A=M/r is integral before any argument about its determinant is used. There
is no hidden rational-normalization assumption.

If a prime p divides det A, choose a nonzero row in its left kernel over Fp.
Scale one nonzero coordinate to1, and lift it with that coordinate exactly1
over the integers. This gives a primitive integer row v whose image vA is
divisible coordinatewise by p. Nonsingularity ensures vA is not zero, so its
content is at least p, contrary to universal content preservation. Thus
det A=+/-1. Conversely integral A and A^(-1) preserve the coordinate ideal,
so unimodularity suffices.

Even if the hypothesis were initially stated for a positive real r, basis
rows would force r to be a positive integer. The text already explicitly
assumes integer r>=1. Its universal quantifier is over all integer rows,
not merely physically realized nonnegative tournament states; that stronger
domain is essential to this proof and is correctly retained in the theorem.

## Independent word products and determinant obstruction

The independent script [audit_information.py](audit_information.py) transcribes
only the seven matrices from the inherited recurrence, multiplies on the
right in displayed processing order, and evaluates determinants directly
by the three-term cofactor formula. It obtains:

| Word | Added vertices | Entry gcd | Determinant | Determinant after division |
|---|---:|---:|---:|---:|
| W5 | 22 | 5 | 125 | 1 |
| W9 | 16 | 9 | 729 | 1 |
| W21 | 34 | 21 | 9261 | 1 |

The displayed U9 and the full displayed M(W21) agree entry for entry with
the independent products. The seven letter determinants are
(1,1,1,9,5,9,21), in order (s,c,q0,q1,q2,q4,q5).
As additional hostile/positive controls, all343 integer rows in[-3,3]^3,
including zero and negative rows, satisfy the predicted content scaling
for each displayed word. These finite controls support the implementation;
the unimodular theorem supplies the universal quantifier.

For necessity, let a,b,c be the valuations of r at3,5,7. No other prime
can divide r. If n_5 counts q5, n_2 counts q2, and n_9 counts q1/q4,
det M=r^3 forces n_5=3c, n_2=3b and2n_9+3c=3a. Thus a-c is nonnegative
and even. The proof does not assume singular directions align merely because
the determinant permits them. Sufficiency is separately certified by the
three actual products and their concatenation. Noncommutativity is harmless:
the scalar factors commute and the ordered product of integer unimodular
normalizations remains integer unimodular. The empty word realizes r=1.

## Physical word and observer sidecars

Each displayed letter belongs to the inherited seven-letter language.
The processing order is explicitly the reversed source-component order;
the inherited physical charts also reverse their vertex labels as specified
there. The old marked source, literal word, chart identities and fresh vertex
names determine the actual graph: attach each block with all old-to-new arcs
except its first vertex pointing to the previous marked last vertex, then
mark the new last vertex. No selection of a fractional matrix or illegal
inverse operation is used in the construction.

Consequently the all-order recurrence supplies actual tournament outputs
for these words, adding16 and34 vertices, respectively, to the supplied
nonempty source. Counts of a source plus34 new vertices are not independently
enumerated in this audit. This use of a proved recurrence is stated openly
rather than replaced by a finite small-order check.

The normalized state vU is generally not itself a realized tournament
state, and the inverse U^(-1) is not asserted to be a legal bridge word.
The readout is r(vU)_1+r(vU)_2, not generally rH. Thus scalar21 in this
classification does not realize Hamiltonian-path count21. The word and
original graph retain attachment provenance; the three numerical coordinates
do not reconstruct arbitrary-cut responses or the two-port observer.
The information proof correctly states these boundaries.

## Reproduction

From the repository root run
`python3 -B 05-knowledge/results/central_modular_chanpak_20261110/central/audit_information.py`.
The JSON records each full product and its normalized matrix. Repeat with
`python3 -B -O`; the output is identical. No assertion-dependent checks,
floating arithmetic, imported information-lane implementation, or broad
word search is involved.
