# Independent mathematical audit of the matrix boundary packet

2026-11-10. **PASS after a minor citation correction.** No mathematical
counterexample or required change to a displayed formula was found. This
review distinguishes all-order arguments from the finite checks actually
performed. It does not certify external novelty or tournament completeness.

## Construction

Reviewed `construction/PROOF.md` and `construction/check.py`. The exact law

    H(G(k,S)) = 1 + sum_(i in S) 2^(k-1-i) + [0 in S]

holds by the stated compression bijection and the unique path avoiding
a->b. The strongness criterion0 in S is both necessary and sufficient.
The interval, direct inversion, logarithmic vertex bound and exact missing-
pair budget follow. Their fixed-order scope is the displayed family, not
the entire oriented spectrum. The tournament slice S=all has H=2^k+1;
the construction of other totals genuinely uses missing pairs.

The independent `compact_audit.py` imports no producer code and has already
checked525 targets (every3 through512, plus selected binary boundaries
through4099), including63 independent permutation checks. It additionally
checks the generic exact2n+1 compiler on all4165 loopless labelled digraphs
through order4 and verifies the unsourced correction n*C_Ham(D). Those
normal and optimized receipts agree byte for byte.

One minor provenance issue was reported and repaired: the four-vertex H2
rows(2,4,9,2) belong to Session235's `audit/AUDIT.md`, whereas PB5 gives a
different five-vertex witness. The four-vertex graph is valid: its two
Hamilton paths are0123 and3120. No mathematical claim was changed.

## Weighted enumeration and cancellation fibres

Reviewed `enumeration/NOTE.md` and `enumeration/check.py`.

- The pair-orbit lists are correct: each vertex cycle of length L supplies
  floor((L-1)/2) preserving orbits of length L and, if L is even, one
  reversing orbit of length L/2. Two distinct vertex cycles of lengths
  a,b supply gcd(a,b) preserving orbits of length lcm(a,b).
- Preserving orbits contribute x^ell+y^ell+2z^ell; reversing orbits admit
  only the symmetric pair states and contribute x^ell+y^ell. Thus weighted
  Burnside gives the stated class polynomial, with centralizer denominator
  z_lambda. Killing both symmetric states correctly removes every even
  vertex cycle in the tournament specialization.
- Loop bits are constant on vertex cycles for fixed unmarked graphs,
  giving prod(1+t^L). They become independent only after a whole ordered
  Hamiltonian path has fixed all vertex positions.
- Cancellation fibres over an actual oriented O are exactly decorations
  of its missing pairs and vertices, modulo Aut(O). This gives the stated
  two-variable orbit polynomial. A full relabeling between two preimages
  must preserve their common cancellation O, so no extra quotient is lost.
- H(D)/|Aut(D)| and H(C(D))/|Aut(D)| count path-marked types because the
  action on ordered spanning paths is free. The fixed-path factors are
  respectively(y+z)^(n-1) and z^(n-1); all other pairs contribute
  x+y+2z. These weighted path masses are explicitly distinguished from
  the unweighted number of isomorphism classes.
- The expected values refer to uniform labelled matrices in the declared
  pair-state palettes. No independence between different paths is used.

The producer's declared universes match its code: all loopless labelled
matrices through4;139 pair-orbit profiles through10;104 fibre polynomial
checks meaning52 oriented types times the two loop conventions.
The saved normal and optimized receipts are byte-identical.

Fresh independent controls in `final_audit_check.py` reconstruct all raw
class coefficients at orders1–3, with1,4,64 labelled inputs and1,3,16
types. They explicitly construct the C4 fibre's4 loopless and64 looped
labelled preimages, giving3 and22 raw types. The three ordinary H values
are4,4,8. The script uses literal permutation paths and full relabeling,
without importing the producer. No fresh complete order-four census or
external sequence-source audit is claimed by this second review.

## Endpoint contractions and resolution overlap

Reviewed `completion/PROOFS.md` and `completion/check.py`.

CT1 is an exact bijection: selected unordered pairs must form a linear
forest; each nontrivial component has two possible traversal directions,
with illegal directed choices removed. Contracting a legal choice retains
the final and initial vertex of every block, so the quotient may have either
defect state even when the original graph is a tournament. Isolated blocks
are counted once. The empty case is consistent with H(empty)=1.

CT2's2m-vertex lift is valid for every loopless digraph. The pair controlling
D[i,j] is physically distinct from the pair controlling D[j,i], so no
orientation prescriptions conflict. Requiring every named matching arc
gives a path-compression bijection. Strong input of order at least2 gives
strong output by lifting directed walks. The universal response bound is
correct, and applies to a selected matching response, not the full H.

For CT3, the indexed resolution counts(5,5) and the entire unmarked
resolution type-list agree in the two order-four examples, but their
intersections are4 and2 and their unions are6 and8. Exact labels across
the resolutions are the lost information. The stated minimality has the
correct restriction to semicomplete graphs with exactly one reciprocal
pair; the complete order-two/order-three check pays that claim. Ordered
padding preserves the collision, and the note correctly avoids claiming
that padded graphs remain strong.

Fresh independent checks in `final_audit_check.py` list the actual paths
and canonical types of the two collisions. They also enumerate every
permutation of the eight-vertex H7 lift: the full tournament has79 paths,
exactly7 of which use its entire selected matching. The ten-vertex H21
lift is proof-reviewed and its saved producer evidence inspected, but this
second review does not freshly enumerate all its permutations.

The producer code matches the stated complete universes:4166 digraphs at
orders0–4 and262666 selected-pair responses;70 literal lift controls at
orders0–3; and separate named target and padding controls. The normal and
optimized receipts agree byte for byte. The order-three examples realize
the two defect types separately; publication wording now explicitly says
that each phenomenon occurs at order three. There is no claim that a
two-vertex quotient contains two different unordered pairs. No equation or
finite count needed repair.

## What remains outside the audited conclusions

The identities generate exact class counts and evaluate exact marked
responses. They do not infer the support of H from those aggregate counts,
preserve H under raw-code reinterpretation, or turn a constrained response
into the full count of a tournament. The independent receipt records its
source hashes and the smaller universe actually rerun here. The tournament
assertion that7 and21 are the only missing odd totals remains open.
