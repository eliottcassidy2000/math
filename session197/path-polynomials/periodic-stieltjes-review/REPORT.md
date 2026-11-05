# Independent review: periodic physical Stieltjes words

**PASS, no correction.** The all-order Hankel total-nonnegativity theorem is valid for a fixed parity-returning physical word repeated from K1 or Tr2. The output-array argument, infinite-matrix locality, physical mirror decomposition and stated scope boundaries pass.

## Frozen binding

The reviewed source is `session197-odd-cycle-path-polynomials/`:

- `PERIODIC_STIELTJES_WORDS.md`: `82b941c9ce53b41ca73ee0224b240fe59c61df701df5fde6f9d29283f0fad166`.
- `verify_periodic_stieltjes_words.py`: `6464caf1620bb77e88924cfe7e43485384c7090828783cfa86108bc058a1f391`.
- Both frozen source receipts: `25b76f0827e580aff425cdf69ad17189d0f59c037f68f0554c96cc06e2368460`.

The named `PROOF.md`, `STIELTJES_VACUA.md` and `verify.py` dependencies were copied unchanged with the packet. Their hashes match the source's Section7. The source verifier was run normally and under `-O` in this copied snapshot; both outputs exactly reproduce the frozen receipt. No original source was changed.

## All-order proof review

The column matrices are correctly oriented. Ordered padding by two vertices is lower bidiagonal in either central-factorial basis. A triangle wrapper has subdiagonal1 and superdiagonal b_k. Consequently its exact L-diagonal-U factorization has lower entries1/ell_(k-1) and upper entries b_k/ell_(k-1). Their multiplication gives a_k on the diagonal, so no transpose or reciprocal is missing.

The pivot proof from the named vacuum packet is sufficient at every index: the odd induction gives ell_n>=(n+1)^2; the even finite boundary through8 is followed by the exact residual (3n+5)(n-9). Thus the factors are nonnegative bidiagonal/diagonal matrices. Every specified minor of their finite product has only finitely many intermediate indices, because the factors are banded. Finite Cauchy-Binet proves total nonnegativity without an analytic limit or an unproved positivity assertion for an infinite operator.

The output-array induction is valid. If selected row0 lacks selected column0, its minor vanishes. Otherwise expansion has positive sign because both are first in increasing order. Removing them decreases the pair (sum of row indices, minor size) lexicographically. If row0 is absent, subtracting1 from every selected row index strictly decreases the first coordinate. Row-finiteness makes the intermediate Cauchy-Binet sum finite: each earlier output row has finite support. This supplies a genuine well-founded induction for arbitrary minors, not only leading minors.

Both T and its transpose are row-finite because T is banded. The identity H=O(T) O(T^T)^T has the stated order: its summand is (T^i)_(0,k)(T^j)_(k,0). Any specified minor uses only the finite union of the corresponding output-row supports. Thus another finite Cauchy-Binet application proves every Hankel minor nonnegative. The shifted principal matrices are submatrices of this Hankel matrix; applying the classical Stieltjes moment criterion is therefore justified.

The cross-parity X matrices are identity from q to r and the positive padding bidiagonal matrix from r to q. An even number of X letters returns to the original coefficient space, so one fixed product T is iterated. For an odd number of X letters, the stated conclusion is only for W squared and the even subsequence; the proof does not silently iterate a matrix with the wrong parity.

The physical mirror decomposition is correct. Larger side labels in I form the initial transitive padding; larger labels in O form final padding. Centrally I_j points to the module (I_<j, input, O_<j), that module points to O_j, and the reversed edge is O_j->I_j. This gives exactly P^(t-1-j) D P^j, in column-action order.

The empty/padding-only periods give constant moments, so nonnegativity rather than universal strict positivity is the appropriate conclusion. The one-time-prefix hostile remains a real boundary: Tr4 followed by repeated D gives1,17,241 and determinant-48, whereas repeating PD from Tr2 repeats its padding and gives1,17,913. No arbitrary initial-profile statement follows.

## Independent implementation controls

`independent.py` does not import source code. Its monomial triangle operator is computed directly by the binomial identity for (x+1)^k-(x-1)^k. A separate central-factorial vector implementation reconstructs whole polynomials to check matrix orientation and parity, rather than checking only scalar readouts.

Normal and optimized runs agree on:

- 184 distinct P/D or parity-returning D/X words through length6, including overlap removal.
- 1,840 full virtual-basis column identities in both parities.
- 368 independently computed moment series through seven repetitions.
- 2,622 production minors and 2,622 Hankel minors, including all minor sizes through4 in the tested corners.
- 200 exact rational pivot controls.
- 60 actual physical graph counts by literal directed-word recursion, independent of the source's subset counter.
- 84 literal mirror block equalities, including nontransitive inputs and all marked positions for side sizes1–6, checking5,460 arcs.
- The exact one-time-prefix negative determinant.

These are finite implementation controls; the all-order assertion rests on the preceding proof review.

Independent checker SHA256: `1c8801a77c49749c10caf3b745429dc4b3a53549be8c9b913c1d214e8c8dfcaf`.
Independent normal/optimized receipt SHA256: `b967c4742685b7109b4a35224b39032f02b797876932c6b615aa4eff40a344f2`.

Run `python3 independent.py` or `python3 -O independent.py`. The copied source replay is `python3 source/verify_periodic_stieltjes_words.py`; all code uses only the Python standard library.
