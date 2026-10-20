# Independent audit: divisor/tetrahedral and F5/icosahedral bridges

Date: 20 October 2026 (UTC).

## Verdict

**APPROVED for the exact eight-file manifest preserved in `exact_manifest.json`.**
The main claims are correct with their stated markings, quotients, and scope.
There is no mathematical publication blocker. Two optional wording improvements
are listed below; neither changes a map, a count, or a proof.

This is an independent review of the staged `session17-prime-shapes` packet.
No repository or remote writes were made. Review files live only in this audit
directory. The working repository reported `main...origin/main` without dirty
paths when inspected. Remote publication is not claimed by this audit.

The original eight manifest members and manifest itself are frozen under
`frozen/`. All eight byte counts and SHA-256 values match the supplied manifest.
The seven entries in `SHA256SUMS` are the seven other members, intentionally
excluding the checksum file itself. That checksum file is, in turn, hashed by
the publication manifest. The optimized duplicate receipts are not publication
members and are not needed for completeness.

## Reproduction and independent coverage

Producer checkers, executed without importing or changing them:

```
python3 -B session17-prime-shapes/verify.py
python3 -B -O session17-prime-shapes/verify.py
python3 -B session17-prime-shapes/verify_field.py
python3 -B -O session17-prime-shapes/verify_field.py
```

Results: 1,968 and 1,036 active checks respectively. Each normal result agrees
byte-for-byte with its supplied receipt and its optimized result.

Independent checker:

```
python3 -B session17-prime-field-audit/independent_check.py
python3 -B -O session17-prime-field-audit/independent_check.py
```

Both runs pass 3,124 active checks and produce identical JSON. The checker
imports no producer modules. It uses six-bit graph sets, row-major field
matrices, literal directed paths, recovered source coordinates, and actual
two-complement moves, rather than relying on the supplied checkers' outputs
as proof. Complete distance rows, divisor-to-complement-pair rows, and weighted
legal-return rows are retained in `independent-normal.json`.

The finite universes are small and fully declared: 64 four-vertex graph/tournament
bit patterns, 24 four-label permutations, 625 field matrices, 12 signed-vector
classes, 120 five-label permutations, 120 determinant-one matrices, 5,040
seven-label vertex words, and the three named source-specific Collatz words.
No large census, all-order computational extrapolation, or new conjectural
dependency was used.

## 1. Arithmetic definitions and complete classification

The note correctly excludes both 1 and N from the divisor carrier. Therefore
F=product(a_i+1)-2, S=2^r-1-[N squarefree], and U=r-[N prime]. In particular
U is not Omega and is not always omega. The empty identity at a prime is
handled correctly, and the optional N=1 boundary is explicitly excluded.

For squarefree composites F=S while U>0. In the nonsquarefree case the
necessary and sufficient equation is product(a_i+1)=2^r+r+1. At support
one it forces exponent three; at support two the required product seven
cannot be split into two integers at least two; at support three the only
three-factor product twelve is 3*2*2. At support at least four the minimum
3*2^(r-1) exceeds 2^r+r+1. Both directions therefore give exactly
p, p^3, and p^2qr, with distinct prime labels in the last shape.

This correctly inherits DB1--DB3 rather than claiming rediscovery. It also
excludes every square greater than one. The table at 25,125,625,776,576,
the arithmetic identity with difference p^2(p+3), and 776=sigma(5^4)-5 all
check exactly. The producer independently compares direct divisors with the
formulas at its ten named inputs and the classification at every N=2..1500.

## 2. Literal arithmetic fibre map

For p^2qr, the nonsquarefree proper divisors are exactly p^2,p^2q,p^2r.
The inherited atom map ell -> p*lcm(p,ell) sends p,q,r to these three
objects. Pairing these images with their atoms and fixing pq,pr,qr,pqr
gives an involution with cycle type 1^4 2^3. Its quotient has seven
fibres, with precisely the three prime-labelled fibres doubled.

Thus 10=7+3 is a true fibre-size identity despite U being a subset of S.
At p^3 the same construction is the two-point swap, and at a prime the
carrier is empty. Divisor complement at p^2qr instead has cycle type
2^5. The distinctions are essential and are stated correctly.

The failure of divisibility preservation is concrete: at N=60, 2 divides
6, but tau(2)=4 does not divide tau(6)=6. This also proves that the new
geometric involution is not an unmentioned order symmetry.

## 3. Tetrahedral carrier, map, and S4 obstruction

The twenty three-edge K4 graphs consist of twelve paths, four stars, and
four triangles plus isolates. Complement is fixed-point-free and gives
ten unordered pairs. Opposition commutes with complement because it is
a permutation of the six edges.

The map Phi is well-defined on all four path-order descriptions of a
complement pair. Reversal has even parity. Complementary order
(c,a,d,b) has odd relative parity and swaps endpoint and middle edges,
so the formula returns the same letter. For a target edge there are
exactly two positive endpoint orders, related by reversal; this supplies
the inverse. Corners have the unique star/triangle inverse.

The literal opposite-edge path (b,a,d,c) has even relative parity and
swaps the target edge with its opposite. Hence Phi intertwines opposition.
Four star/triangle pairs are fixed; the remaining six form three pairs.

For every g in S4 the correct identity is
Phi(gx)=J^parity(g) g Phi(x). This is checked for all 240 pairs (g,x),
and its proof is just the parity change in the defining path order.
Thus fixed-orientation equivariance is A4-equivariance; allowing the
orientation marking itself to transform restores ordinary marked
S4-equivariance. Reversing the chosen orientation gives J Phi.

There is no equivariant bijection for the two specified unmarked S4
actions. A transposition fixes two source pairs but four target letters.
The alternate stabilizer explanation also checks: the path-pair
stabilizer has element orders 1,2,4,4, whereas the edge stabilizer has
orders 1,2,2,2. Equal orbit size does not repair this distinction.

## 4. Changed divisor table and what it preserves

Every literal row of B is distinct, covers the ten letters, and satisfies
B tau=J B. The entire composite Phi^(-1) B is independently checked,
not merely the two conjugacies separately. The full divisor-to-graph
table is retained in the independent receipt.

The squarefree set becomes the four corners plus the three edges incident
to zero; its prime subset is exactly those three incident edges. The
quotient and its singleton/doubleton families are preserved. Arithmetic
multiplication, divisibility, size, and divisor complement are not.

The proper divisor poset has exactly the q/r swap as its nonidentity order
automorphism: adjoining bottom and top restores the lattice, whose
join-irreducibles are a distinguished two-element chain and two isolated
elements. No abstract three-cycle of letters becomes an order automorphism.

## 5. Actual balanced source and the 625/576 operation

The quoted source file defines T_(3,3) as the transitive seven-vertex
tournament with only 0--3 and 3--6 reversed. Its source decoder is
P(I,K)=I, sigma(A\K), 3, A\I, B\sigma(A\K), each displayed run increasing,
with infinity interpreted as A={0,1,2}. All 25 legal paths were independently
enumerated from the literal tournament and then decoded in the reverse
direction. They agree with both this packet's receipt and all 25 rows of
`session17-balanced-fibres/marked_f5_paths.json`.

The standard increasing path is the field zero. Actual positional
anti-reversal is coordinate swap. Choice complementation is simultaneous
negation under the printed five-state marking, with exactly the standard
path fixed. The tournament is independently verified rigid, so these
path-set permutations do not claim additional tournament automorphisms.

Two source paths as columns are all M2(F5), of size 625. Deleting the
standard-path row and column of the pair-choice grid gives precisely
matrices with both columns nonzero: 625-25-25+1=576. Deleting the diagonal
instead gives 600. The distinction is maintained in the main proof.

There are 24*4=96 dependent nonzero-column matrices. For each nonzero
determinant there are 24*5=120, yielding 480 invertible matrices and
120 determinant-one matrices. Independent column sign changes act freely,
so the quotient has 144 ordered pairs, not a projective-line quotient.

## 6. Determinant graph and literal icosahedral certificate

The signed-vector quotient has twelve elements. Squared determinant is
well-defined on independent sign classes. For a fixed class, determinant
zero gives itself and [2v], and the other ten split five/five between
squares one and four. Multiplication by two is a fixed-point-free
involution because its square is multiplication by -1.

The printed N,S,U_i,L_i are twelve distinct classes. Their determinants
give exactly the two-pentagon/two-pole adjacency. In particular
det(U_i,L_j)=2(i-j)+1 gives j=i or i+1. The L-L condition follows from
det(L_i,L_j)=4(i-j); the pole incidences are equally direct.

Independent BFS confirms every vertex has distance profile 1,5,5,1.
The distance-three vertex is precisely [2v], and determinant square four
means exactly distance two. Thus the ordered-pair decomposition
144=12+12+60+60 has all its geometric descriptions justified, not inferred
from the count alone.

## 7. Real tournament return graph, not just a Hamming coincidence

All twelve printed codes are regular five-tournaments, distinct, and on
the same orientation-parity sheet (the odd bit-count sheet). The source
theorem in the repository has SHA-256
db2bd384ca2ad4251870723f2e8cea7a872c65ec487ac10fb0a35e84bbef001b,
identical to the corrected session-16 v3-2 publication payload.

The inherited theorem's endpoint-dipole criterion, minimum-six-arc
argument, regular-five structural classification, and explicit legal
paths establish the legality behind the six-bit criterion. The new note
correctly cites that proof rather than treating XOR count as sufficient
without a source theorem.

For additional independence, this audit constructs every actual directed
P, every actual directed Q in K_P(T), and then K_Q(K_P(T)) for each of the
twelve codes. Regularity is equivalent to matching endpoints. Each code
has 15 identity presentations and five distinct nonidentity neighbours,
each of weight two, for total mass 25. Those five legal neighbours are
exactly its determinant-graph neighbours. No inherited producer code is
imported. This verifies the entire claimed isomorphism at the legal-move
level, including the distinction from the original tournament metagraph.

## 8. SL2 action, printed generators, and actual reflection

The kernel proof is sound: fixing the two coordinate-vector classes makes
a matrix diagonal with diagonal entries +/-1, and fixing their sum's
class makes the signs equal. Thus the kernel is exactly {I,-I}, and the
induced group has order sixty.

Every one of the 120 determinant-one matrices was evaluated independently.
The resulting sixty permutations equal exactly the sixty even relabelling
actions on the twelve tournament codes. Both printed generator-to-vertex-
permutation correspondences agree at every code. Independently generated
closure of the two printed matrices has 120 elements, supplementing the
elementary-row-reduction proof that they generate SL2(5).

Multiplication by two corresponds to converse at every code. For actual
source reflection R(x,y)=(y,x), the matrix 2R has determinant one, and its
action agrees with (0 1)(3 4). Consequently R transports to that even
relabeling followed by converse, exactly as printed. Although 2 is not
its own inverse in F5, it is its own inverse on the sign quotient, which
is why this factorization is correct on the stated carrier.

R-fixed classes satisfy Rv=+v or Rv=-v. The two nonzero sign classes on
each of y=x and y=-x give four fixed classes. Only the first two arise
from individually R-fixed source paths; the second two pair black paths
that R and J both exchange. Adding zero recovers five R-fixed J-orbits.
The text correctly avoids lifting a fixed quotient object to a fixed path.

## 9. Marking boundary and contextual examples

There are 2!*2^2=8 sign-equivariant five-state markings fixing infinity/zero.
The hostile h(x)=x^3 is allowed: it commutes with negation and permutes F5.
The determinant of (1,1),(2,1) is -1, whereas after cubing it is -2;
their squares are 1 and 4. Therefore J, even together with the pointed
choice set and coordinate swap, does not canonically determine this edge
relation. The packet correctly labels it a marked transport and does not
infer a native black-move relation or a tournament automorphism group.

The four-vertex contextual counts are independently confirmed: eleven
unlabelled simple graphs and four unlabelled tournaments. The Royle,
Praeger, Glasby, Freedman, Devillers primary abstract was rechecked at
https://arxiv.org/abs/2204.01947 on 20 October 2026. Its even-graph condition
concerns the edge-reversal sign of automorphisms, not Eulerian degrees;
the stated equinumerosity is up to isomorphism. The cited journal metadata
are J. Algebraic Combin. 57 (2023), 515-524. The new note's scope is correct.

The group-order example loses the action, as claimed: C3^3 and C27 have
order 27 but exponents 3 and 27. The inherited first-descent examples were
also independently recomputed by rational block-envelope DP: 27 first
descends to 23 after 37 odd steps and needs four blocks; 63 first descends
to 61 after 34 and needs five; 75 reaches 1 after three and needs one.
All smaller claimed budgets fail at those same prefixes. These examples
do not support a universal Collatz bound, and none is asserted.

## 10. Nonblocking editorial notes

1. PRIME_BALANCE_TETRAHEDRAL_BRIDGE.md section 5 says the old table placed
   “all six squarefree pair-products” at midpoints. The six divisors are
   p,q,r,pq,pr,qr. Prefer “the six non-top squarefree divisors” or list them.
   The actual old and new tables are correct; this is only terminology.
2. README's “zero-row/zero-column deletion” can be read as deletion of
   rows/columns of individual 2x2 matrices. Prefer “deleting the standard-
   path row and column of the 25x25 pair grid,” giving matrices with both
   columns nonzero. The full note states the operation correctly.

Neither suggestion is required to repair a theorem. If the producer edits
either file, refresh the manifest and checksums, and verify the revised
bytes before using this exact-byte approval for publication.

## Scope retained at close

The packet proves two different connection contracts: a marked involution/
fibre bijection and a marked graph/action transport. It does not prove a
universal trijection, unmarked S4-equivariance, divisor-order preservation,
a canonical path adjacency, extra source-tournament automorphisms, or any
Collatz/LRC conclusion. These limits are part of the approved result.
