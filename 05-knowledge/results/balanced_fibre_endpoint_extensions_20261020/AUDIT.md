# Independent audit: balanced deletion fibres and full endpoint transport

20 October 2026. **Verdict: PASS for the stated mathematical results, with the explicit coordinate and finite-set conventions incorporated in the publication overlay.** No theorem, numerical, printed-path, or marked-table counterexample was found. Infinite claims are supported by the arguments below; finite checks are evidence at their declared bounds, not substitutes for those arguments.

## Audited version and minimal edits

The frozen producer proof is `BALANCED_DELETION_EXTENSIONS.md`, SHA-256
`baeafd394f839681a775fe0c1786d3e87e6fdecb1979ac4781298b3d460f11a0`.
The initially supplied hash `e73e5577e530c88b4950bd2d66f1b776987805dd34660b98b5bb71101ccd1d1c`
is the hash of **check_fibres.py**, not the Markdown proof. This was an identification mismatch, not evidence of an unapproved proof replacement. The complete six-file original packet was frozen before testing and left unchanged.

The overlay makes two narrow wording repairs, preserved in `SUGGESTED_NOTATION_PATCH.diff`:

1. Section 1 initially defines the finite right coordinate as the raw subset `B minus J`. The balanced diagonal and field marking use the **reflected** coordinate `K=sigma(B minus J)`, with infinity fixed. The inherited microcosm has this convention, but the standalone note must state the change explicitly. Without it, applying ell to the raw right-shore subset is a type mismatch. No count or kernel changes.
2. The arbitrary-surjection identity is explicitly scoped to **finite sets**, as required by its integer subtraction and counting interpretation.

The original proof remains in `frozen/`. The overlay also supplies a README cross-link to the separately audited prime-field packet. No claim of unmarked canonical determinant adjacency is introduced.

## 1. Rectangular deletion and all-depth fibres

Independently derive the decoder by removing the centre c=a. All remaining arcs increase. Hence both the prefix before c and suffix after c increase, and each splits into its A and B portions. If an A vertex follows c, its first is 0, the only A out-neighbour of c; if a B vertex precedes c, its last is a+b, the only B in-neighbour of c. Thus I is a subset avoiding 0 or the full A, while J is empty or contains a+b. Conversely these conditions make every join legal for a,b>=2.

Deleting only interior vertices keeps c,0,a+b. The four runs restrict and their admissibility conditions persist. Each finite coordinate can independently contain or omit each removed interior vertex; each exceptional full/infinity coordinate has one extension. Thus the four fibres are exactly `2^(r+s),2^r,2^s,1`, for every permitted deletion depth. This proves actual path deletion and surjectivity, not merely equality of counts. The upper and lower parameter bounds and the requirement to retain the three distinguished vertices are essential.

For reflected balanced additions, write `m=q-1` and let each finite state have D lifts. Then `q'=Dm+1`. The total is `(Dm+1)^2`, giving `D^2 H-2D(D-1)q+(D-1)^2`. The black count is `Dm(Dm+1)`, which equals `D^2 Delta-D(D-1)(q-1)`. Its direct disjoint decomposition is

- finite unequal parents: `m(m-1)D^2` black lifts;
- mixed parents: `2mD` black lifts;
- finite equal parents: `mD(D-1)` black lifts;
- exceptional equal parent: no black lifts.

This also checks the direction carefully: a black parent has only black lifts, but a black child may delete to a blue parent. The same diagonal/off-diagonal classification proves the finite nonuniform-fibre formula. A noninjective map cannot be a homomorphism between loopless complete graphs because an edge inside one fibre collapses to a vertex.

## 2. Endpoint closure and its boundary

For finite nonempty I, the path start is min(I). Infinity starts at 0. Empty I starts in J, or at the centre if J is empty. The right coordinate has the reflected statements at the end. Consequently endpoints determine exactly the three types infinity, empty, and finite nonempty in each coordinate, but not the entire subsets.

The new left label L is larger than every old left finite label, and the new right label R is smaller than every old right label. A finite nonempty subset therefore retains its endpoint. For empty I, choosing its new bit makes the start L; with that bit absent, R becomes the start exactly if the right bit is absent. The symmetric argument gives the stated end rule. This proves the kernel on each individual path, with multiplicities, not just after summation.

Blue parent paths have equal full normalized coordinates, and the child is blue exactly when the two new bits agree. Restricting the kernel is therefore valid only on separately retained blue counts. At k=3, these two literal paths share endpoints (1,5):

- blue: `(1,4,6,3,0,2,5)`;
- black: `(1,6,3,0,2,4,5)`.

Both were checked against actual arcs and positional anti-reversal. Thus the endpoint-only colour guard is substantive. The three-type transfer matrix, its tensor square, the stratum polynomial, and the three rational formal series follow directly from one infinity lift versus two finite lifts. None closes arbitrary tournament insertion.

## 3. Full endpoint unimodularity

For the upper-right bulk block, the two minimum choices are independent, giving the stated outer product. In a side block, min(S)=alpha and max(A minus S)=beta force:

- the terminal interval when beta=alpha-1;
- alpha in S and beta out of S, with exactly beta-alpha-1 free interior bits when beta>alpha;
- no solution otherwise.

This proves c(alpha,beta), including the zero diagonal. The three remaining isolated supports are precisely the empty/infinity boundary cases. The supports are disjoint, so their number is `k^2+2*k(k-1)/2+3=2k^2-k+3`.

For a nonzero permanent term, column 0 forces row 1; successively columns 1 through k-2 force rows 2 through k-1. Previously possible alternatives have already been removed. Column k-1 forces row k+1, and column k forces row 0. The reflected chain then forces rows k+2 through 2k against columns k+1 through 2k-1, leaving row k for column 2k. Every entry is one. Hence there is exactly one perfect matching with weight one.

The row-to-column permutation is the single `(2k+1)`-cycle displayed in the proof. Its sign is positive, so both determinant and permanent are exactly 1 for every k>=2. This proves integer invertibility and full rank in every characteristic, despite the rank-one bulk. It supplies no replacement for the distinct-component tensor required by general insertion.

## 4. Homogeneous twins and printed-witness regression

The Markdown itself was parsed to obtain its label assignment, five outgoing bit rows, source path, and printed remainder. The parsed rows are `(10,8,19,20,3)`; the source is `(2,1,3,4,0)`. They form a tournament with homogeneous vertices 0,1 and a valid Hamiltonian path. Deleting vertex 1 gives the printed `(2,3,4,0)` with the unique bad join `2->3`, whose reverse is present.

For minimality, an orphan from internal deletion requires consecutive a,v,b and the cyclic join b->a. The other homogeneous twin u then satisfies a->u->b. On four vertices u can occur only before a or after b in that path, forcing respectively u->a or b->u, a contradiction. End deletion cannot create a join. Thus five is minimal.

The audit deliberately mutates the **printed text** to `(v,a,b,c,u)` and runs the same parsing/validation routine. It is rejected specifically because the parsed source is not a Hamiltonian path. This regression protects against the prior failure mode where a correct generic program did not certify a wrong printed example. The supplied 25-entry F5 JSON is likewise read and checked as data, not silently regenerated as a replacement.

## 5. Literal join and marked F5 carrier

All cross arcs in A->B point forward, so a path cannot return to A after entering B. Restriction and concatenation give the actual product bijection. For S=T_(3,3), the 25^2 paths in S->S are therefore literal concatenations, and excluding the standard component path independently removes a row and column, not the diagonal: `625-25-25+1=576`.

The literal directed cycle `(0,1,2,4,5,6,3,0)` proves S is strong. Its score classes prove rigidity, and reflection is its unique anti-automorphism; all 5040 vertex permutations independently confirm both symmetry counts. In the ordered join, SCC order forces an anti-automorphism to exchange the two components. Its action on component paths is `(P,Q)->(r(Q),r(P))`. Fixed points are exactly `(P,r(P))`: 25 total and 24 with neither component standard. Thus the full counts are blue/black 25/600 and the 576 bulk counts are 24/552, leaving boundary counts 1/48.

The five-state label table intertwines complement with negation. Applying it to the **normalized** left/right coordinate pair yields all 25 field vectors, with the standard increasing path at zero, actual anti-reversal as swap, and the independently defined complement operation J as negation. There are exactly `2!*2^2=8` sign-equivariant five-state markings. The nonzero sign quotient has 12 classes; it is not the six-point full projective-line quotient by all nonzero scalars. Its ordered square has 144 classes, four representatives each in the 576 bulk.

`cross_packet_receipt.json` independently compares the balanced table with the prime-field packet's receipt: all 25 literal paths, state masks, and field vectors agree exactly. The shared marking is therefore consistent. The determinant relation remains an additional marked graph relation, not a canonical adjacency of source paths and not an automorphism action on the rigid source tournament.

## 6. Exact checks and reproducibility

The producer script was run in separate replay directories, normally and with `-O`; both outputs and both generated marked tables are byte-identical to the frozen results. The independent audit imports no producer code. Its distinct methods are adjacency-matrix/set recursion, exhaustive vertex permutations at small orders, inverse coordinate extraction from literal paths, rational Gaussian elimination, and forced-column matching elimination.

Independent bounds:

- All a,b=2,...,5 and every allowed interior deletion: 676 maps and 112,896 individual source-path deletions, all target fibres.
- Every parent and literal child for balanced k=2,...,6: both full and blue endpoint kernels, black-parent lift direction, and actual colour/coordinate agreement.
- Full endpoint coefficients against literal paths for k=2,...,7; subset coefficients, support, unique unit perfect matching, and exact determinant for k=2,...,12.
- Every labelled tournament on two through four vertices, all homogeneous pairs, both twin deletions and all paths: 412 valid controls, plus the parsed five-vertex hostile.
- All 625 literal join paths and all eight total/bulk/boundary colour counts; the source Hamiltonian cycle and all 5040 source-vertex permutations for rigidity and anti-automorphism uniqueness.
- All 25 supplied marked paths; legality, inverse states, vectors, blue status, R and J; all 120 possible five-state bijections for the eight-marking count.
- 1,600 all-depth arithmetic controls, including q=1 and arbitrary D=1,...,16; 1,364 positive nonuniform fibre vectors.
- Both a malformed printed source-path regression and a same-endpoints/different-colours hostile.

Run from this directory:

```
python3 check_fibres.py > check.out
python3 -O check_fibres.py > check.optimized.out
cmp check.out check.optimized.out
python3 audit_balanced.py . > publication-audit-results.json
python3 -O audit_balanced.py . > publication-audit-results.optimized.json
cmp publication-audit-results.json publication-audit-results.optimized.json
python3 audit_balanced.py > audit-results.json
python3 -O audit_balanced.py > audit-results.optimized.json
cmp audit-results.json audit-results.optimized.json
```

The final two default runs audit the preserved original proof/table in `frozen/`; the explicit `.` runs audit the actual publication overlay. Every validation uses explicit exceptions, not removable assertions.

## Scope, provenance and publication boundary

The repository baseline inspected read-only was `f3e658d7cc105d5e367e4dd565812a81e127a6fc`. The agent contract, bounded startup packet, relevant current frontier, active guardrails, research protocol, repaired deletion mechanism, and incoming two-return/microcosm proofs were read. No literature-priority claim is made. The optional message-status reader was unavailable because that checkout has no `.machine-id`; no registration or repository mutation was attempted. This does not affect mathematical reproduction.

All audit/staging writes are outside the repository. Producer files are preserved byte-for-byte. The proposed publication is add-only under `05-knowledge/results/balanced_fibre_endpoint_extensions_20261020/`, with exact source/destination hashes in the external publication manifest. Repository and remote writes remain with the designated publisher after root review.
