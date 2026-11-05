# Primitive normalization does not preserve physical sources

Session 199. **PROVED:** dividing all four response coordinates by their
greatest common divisor can destroy physical source realizability. There is
an infinite explicit family with this failure. **FINITE-EXACT, complete by
cover height:** the first failure has U=39 and is unique at that height.

This is a hostile test of a proposed Euclidean normalization, not a claim
that any new Hamilton-path count is absent from the unrestricted spectrum.

## 1. Why this normalization was a plausible bridge

The historical Gaussian/Berggren work has a genuine primitive reduction:
THM-3336, `01-canon/theorems/THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation.md`,
divides the two coordinates of a Gaussian product by their common content.
The quotient remains a primitive integral spinor, and its content cocycle
records the change of norm and determinant. THM-3353's split-prime compiler
then performs a Euclidean inverse on positive coprime odd parameter pairs,
where every inverse step is proved to remain in the parameter domain.

The attempted transfer was to normalize a large modular source witness by
its four-coordinate content, thereby lowering ordinary height while retaining
its projective direction. The following small example identifies the first
failed implication: an integral primitive response vector need not be a
physical source. There is no general source-space analogue of the saturated
primitive-spinor domain.

## 2. A literal source with a nonphysical primitive quotient

A source is an actual rooted tournament with one Hamilton path starting at
the root, with response (D,E,U,d) and the usual conventions. Start at the
two-vertex transitive source t=(1,0,2,1). The actual chronological letters are

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D-E+U,d,2U-d,d),
    J(D,E,U,d)=(D-E+U,D-E,D+U,D).

The complete source grammar and all physical interpretations are in the
frozen dependency. The word OOZOJ has the successive responses

    (1,0,2,1), (1,0,4,1), (1,0,8,1),
    (9,1,15,1), (9,0,30,1), (39,9,39,9).

Thus s=(39,9,39,9) is an actual seven-vertex source, with content3. Its
primitive integral quotient is

    p=(13,3,13,3).                              (1)

**Theorem 1.** The profile p is not physically realizable.

**Proof.** The complete balanced-source normal form says that every
nontrivial source with D=U is J(O(W)) for an actual source W. If its balanced
profile is (h,e,h,e), then W has

    D_W=e,       U_W=(h-e)/2.

For (1), W would have D_W=3 and U_W=5. But the entire D=3 source slice has
U=3*2^k, k>=0. To see this, the initial strong component of a weak start-one
source is followed only by singleton components. Its path count is therefore
the path count of that initial component. A strong tournament with H=3 has
at most three vertices, because cutting a Hamilton cycle at each vertex gives
at least its order many distinct paths. It must be the directed triangle,
whose source cover count is3. Every later singleton doubles that count.
Consequently U_W=5 is impossible. QED.

The quotient p satisfies all nonnegative source-coordinate and parity
conditions, including genuine odd deletion count. Its failure is exact
physical membership, not a negative coordinate or a parity violation.

## 3. An infinite family with the same failure

For k>=1 let s_k=O^k(s). Then

    s_k=(39,0,39*2^k,9),
    s_k/content(s_k)=(13,0,13*2^k,3).             (2)

Both formulas are exact physical constructions before division, and their
content is exactly3.

**Theorem 2.** Every quotient in (2) is nonphysical.

**Proof.** Every weak source is a strong source core, or a singleton,
followed by a nonempty sink tail. A strong source has

    1 <= U_core/D_core < 2.

For example, the positive source grammar gives U<=2(D-E), and E>=1 on a
strong source. A singleton has ratio1. Each sink doubles U without changing
D. Thus if U/D=2^k, the unique weak-tail length is exactly k and its core
has U=D. Applied to (2), a putative quotient would have a balanced core
with (D,E,U,d)=(13,3,13,3): sink deletion preserves d, and balance forces
E=d. Theorem 1 excludes it. QED.

Thus common-divisor removal fails at arbitrarily large physical cover
counts; the failure is not an isolated small exception. It does not rule
out a different, separately constructed primitive representation of a
scalar target.

## 4. Exact first-failure certificate

The complete O/Z/J grammar has 97 nontrivial profiles with U<=39. Exactly
17 have content greater than one. Dividing each by its content and testing
against the same complete bank finds exactly one nonphysical quotient:

    (39,9,39,9) / 3 = (13,3,13,3).

This certifies minimality in **source cover count U**, not in graph order,
Hamilton count alone, or word length. The cap is complete because U strictly
increases under every physical letter. Each quotient has no larger U than
its input, so checking against the cap39 bank is an exact all-order
membership test, not a graph-size truncation.

The verifier regenerates the bank, checks all17 quotients, and independently
counts the seven-vertex graph of OOZOJ. It also checks the displayed source
families and the elementary cone for the false quotient. The all-order
nonphysicality proofs are Sections2-3, not extrapolation of these tests.

The surviving positive question is to identify a **physical** normalization
or replacement that supplies a bounded-height witness on the target fiber.
Arithmetic gcd division alone cannot be that operation. Universal odd
coverage remains OPEN.
