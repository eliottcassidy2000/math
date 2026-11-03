# Independent audit: boundary pairing and parity beat endpoint vertex blocks

## Finding and status

**PROVED, all finite orders:** the inherited endpoint-block carrier is correct but can be compressed further. A selected old path component needs its two endpoints and the parity of its number of interior vertices. Its particular interior vertex set is unnecessary, provided the set of all interior vertices, the isolated vertices, and the pairing of endpoints remain known. The corresponding signed sums form an independent, attachment-closed rational polynomial basis containing every fixed source-oblivious Hamiltonian-count response.

**FINITE-EXACT, independently enumerated:** this basis has dimension 571 at order 6, versus the inherited 661 endpoint-block coordinates. The first 6 dimensions, orders 0 through 5, agree: 1,1,2,7,28,126. Thus the inherited sharpness claims at orders 3,4,5 survive unchanged. The inherited paper explicitly calls 661 only an upper bound, so its statement is not erroneous. Extrapolating its small-order minimum to all orders would be false.

**FINITE-EXACT, independently certified:** the new dimension571 is the exact minimum at order6 for the entire fixed-cut attachment language. The discovery search supplied571 exterior coefficient contexts; this audit independently reconstructed their coefficient matrix and proved full rank modulo101 and103. **OPEN:** all-order minimality, and whether the narrower S/B command language attains571 at order6.

## Source pins and scope

Read-only source commit: `7b76044ca13e6684d1746923285b83dbf83828d2`.

- `05-knowledge/results/primegame_port_spectrum_20261103/FOREST_RESPONSE_RECURRENCE.md`, SHA-256 `0bb6a6b0624495b9ccc4123a2db544f2764f92c01138c4b0b241d6ef2ad73b0b`.
- `05-knowledge/results/primegame_port_spectrum_20261103/GROWTH_FAITHFULNESS.md`, SHA-256 `22c8093fbb11523e112f22183f1cd8d4cc8845ca5f98e67683cd4c3e3174d71c`.

The finite universe is the complete cube of labelled tournaments on `[n]={0,...,n-1}`, with one Boolean bit `a_ij` for every lexicographically ordered pair `i<j`. No Hamilton-path chart, strongness condition, isomorphism quotient, or legality filter is imposed. All means refer to arbitrary rational probability laws on that cube. Linear dimensions are over Q; a known constant removes one coordinate for normalized laws. Arbitrary nonlinear encodings and individual deterministic-state cardinalities are different problems.

Operations preserve every old edge and its label. New arcs are fixed independently of the old source. All finite sequences of arbitrary cuts are allowed. Fixed relabellings and ordinary converse are also allowed. Adaptive source-dependent cuts, conditioning, and shared-source products are not automatically included as first-moment updates.

## 1. Definitions

Put `y_ij=a_ij-1/2` for `i<j`, and `y_ji=-y_ij`. For an undirected labelled linear forest F, orient each nontrivial component from its smaller endpoint to its larger endpoint. Let

- `s(F)` be the product of `+1` for increasing and `-1` for decreasing directed edges in these canonical traversals;
- `Z(F)` be the isolated vertex set;
- `M(F)` be the matching formed by its unordered endpoint pairs;
- `epsilon_F(u,v)` be the parity of the number of interior vertices in that component.

The new signature is `sigma(F)=(Z(F),M(F),epsilon_F)`. Its interior vertex set is implicitly `[n]` minus its isolates and endpoints. Define

`R_sigma(y) = sum_{F: sigma(F)=sigma} s(F) product_{e in F} y_e`.

All forests in one signature have the same number of edges: `r+k`, where r is the number of interior vertices and k is the number of endpoint pairs. Their monomial supports are disjoint between signatures, so the R polynomials are linearly independent. Multilinear interpolation makes them independent as functions on the entire Boolean orientation cube.

The inherited Q basis fixes the vertex set of every component as well as its endpoints. Each R is exactly the sum of the Q coordinates with the corresponding new signature. There are no extra numerical weights in this sum.

## 2. All-context coefficient theorem

**Theorem.** Every Hamiltonian-count response after a fixed source-oblivious attachment belongs to the span of the R polynomials.

**Proof.** Expand a Hamiltonian word in the old variables y, retaining fixed Boolean new arcs. A nonzero old monomial has linear-forest support F. Every component of F must appear consecutively in the word. Contract these components to tokens, retaining their entry endpoint, exit endpoint, and the remaining isolated old and new vertices.

Take F and F' with the same signature. Replace the interiors of each token by the corresponding path in F', in the same endpoint-to-endpoint direction. This is a bijection between the Hamiltonian vertex words contributing to the coefficient of F and those contributing to the coefficient of F'. Interior vertices may move between tokens: their total set is unchanged, each is used once, and no internal vertex has an adjacency outside its selected token. All old-old adjacencies not selected in the monomial occur between token endpoints or isolated vertices, so their factors 1/2 remain unchanged. Every adjacency involving new vertices is likewise unchanged.

For canonical traversal, the selected-edge sign is s(F) or s(F'). Reversing a token changes its sign by `(-1)^(number of edges)`. This factor is the same before and after replacement because each corresponding component has the same interior parity. Therefore the ratio of selected signs is independent of all traversal directions, and is `s(F')/s(F)`.

Hence the coefficient of F multiplied by s(F) is constant on each signature. All nonforest coefficients vanish. These are exactly the coefficient constraints defining the R span. QED.

This proof does not replace actual directed paths in a tournament. It compares formal coefficients at the fair old-edge point. The distinction is essential.

## 3. Closure and gauges

**Theorem.** Pullback under any fixed source-oblivious attachment maps the enlarged R span to the old R span.

**Proof.** Expand one enlarged R polynomial as a signed sum of enlarged forests, then substitute the fixed values of all new variables. Its old supports are linear forests. Fix two old forests F,F' with the same old signature. In every contributing enlarged forest, contract the selected old paths. Each selected old interior vertex already has degree two, so no exterior edge can touch it. Replacing those selected paths by F' leaves the contracted graph unchanged. Thus it preserves acyclicity, degrees, the enlarged endpoint matching and isolated set. The parity of each enlarged path is preserved because every replaced old path has the same parity. The resulting enlarged forest belongs to the same enlarged R signature. Reversing this replacement is the inverse. The canonical-sign ratio is s(F')/s(F), exactly as above. New-edge factors are unchanged. Therefore the pullback satisfies every old signature relation. QED.

A fixed permutation transports the signature and multiplies its basis polynomial by a known sign: for each endpoint pair whose canonical order is reversed, the multiplier is `(-1)^(epsilon+1)`. It does not require the internal vertex allocation. Ordinary converse sends `R_sigma` to `(-1)^(r+k)R_sigma`. These are signed rational coordinates, not nonnegative probabilities.

No closure under conditioning, pointwise products, source-dependent re-marking, or a native fixed-path complement follows from these proofs.

## 4. Exact coordinate count

Choose 2k endpoint labels, pair them, choose r interior labels, and leave the other labels isolated. For a prescribed parity vector of Hamming weight h, an allocation of r distinct interior labels exists exactly when `h<=r` and `r-h` is even. The remaining `r-h` vertices can be allocated in pairs to any block. Thus for k>=1 the number of possible parity vectors is

`A(k,r)=sum_{0<=h<=min(k,r), h congruent r mod 2} binom(k,h)`.

Consequently

`c_n = 1 + sum_{k=1}^{floor(n/2)} binom(n,2k) (2k)!/(2^k k!) sum_{r=0}^{n-2k} binom(n-2k,r) A(k,r)`.

The initial sequence is

`1,1,2,7,28,126,571,2668,13028,65908,347626,1902231,10777900`.

One exact EGF, without an unwarranted claim of a simpler form, is

`C(x)=exp(x+x^2/2) sum_{h,j>=0} (x^2/2)^h x^(h+2j)/(h! (h+2j)!) - exp(x)(cosh(x)-1)`.

The subtraction removes the impossible signatures with no endpoint pair and a positive even number of pooled interior vertices.

The extra term 1 counts the all-isolated forest. A parity vector is a signature, not a count of how many interior allocations realize it; counting those allocations would restore redundant coordinates.

At n=6 the old 661 coordinates become 571: 481 new classes contain one old coordinate, while 90 contain two. The reduction can arise either from two 3-vertex paths exchanging their single interiors, or from 2+4 and 4+2 path sizes with the same endpoint matching and even interior parity. No such ambiguity exists at n<=5.

## 5. Cheapest strict hostile and positive source laws

On six vertices, consider the forests

- F: `0-4-1` and `2-5-3`;
- G: `0-5-1` and `2-4-3`.

Their old component vertex sets differ. Their new signature is identical: no isolates, endpoint pairs `(0,1),(2,3)`, and interior parities `(1,1)`. Both canonical signs are +1. In lexicographic edge-bit order their masks are `10376` and `5392`.

Write `chi_F(a)=(-1)^(sum_{e in F} a_e)`, and similarly for G. The strictly positive laws

`mu_plus/minus(a) = [1 +/- (chi_F(a)-chi_G(a))/4]/32768`

have weights `1/65536`, `1/32768`, or `3/65536`. They are normalized because both nonconstant characters have uniform mean zero. Every fixed source-oblivious attachment has identical means under the two laws, by the coefficient theorem and the equal canonical signs.

The inherited Q coordinates do distinguish them: the Q for F has means `+1/64` and `-1/64`, while the Q for G has the opposite means. Their sum R has mean zero under both. Hence the old Q vector truly contains unobservable first-moment information.

A shared-source second moment sees the difference. Append vertex 6 with cut mask17, so it beats old vertices `{0,4}`. Append vertex7 with cut mask69, so it beats `{0,2,6}`. All other arcs involving each new vertex point toward it. For this fixed two-birth context:

- both mean Hamiltonian counts are `627/2`;
- second moments are `15043377/128` and `15038127/128`;
- the second-moment difference is `2625/64`.

The exact graph counter verifies these values over all 32768 sources. Independent sampling of separate source copies would see only products of the identical means; shared-source products are required for this particular distinction.

## 6. Minimality and operation-class hostiles

The quotient alone proves an upper bound, not observability. The independently verified rank certificate below now supplies the matching lower bound at n=6. Thus the rational source-law kernel has dimension exactly `32768-571=32197`. The inherited bound gave only at least32107. Exactly90 independent directions in the old Q coordinate space are newly discarded, with disjoint Fourier supports. No further linear information within the R span is lost at order6 for the full attachment language.

The original exact minima 7,28,126 at orders3,4,5 are unaffected because the two bases coincide there. They refer to unnormalized linear laws; for normalized laws the counts are6,27,125 measured coordinates.

**Smallest nonlinear operation-class warning, n=3.** Use uniform even/odd old-edge-parity laws on all8 sources. Their H means are both3/2, and after fixed B0 their H means are both7/2. Nevertheless:

- conditioning on `a01=1`, an event of probability1/2 under both laws, gives H means2 and1;
- the adaptive rule “use B0 on even source parity, otherwise append a sink” gives means7/2 and3/2.

Thus neither conditioning nor adaptive source-dependent attachments can be added to the operation language while retaining the first-moment quotient. The restriction is mathematical, not bookkeeping.

Deterministic source faithfulness is not contradicted: a separating family of functions can fail to linearly separate probability laws. The inherited full labelled-source cardinalities remain a different claim.

## 7. Independent minimum-dimension certificate at order six

For an exterior edge set G on N old-plus-new vertices, with every edge incident to at least one new vertex, extract its coefficient in the exterior variables of H_N. This gives a polynomial in the old variables. Because H_N is multilinear in every exterior variable, this extraction is a rational linear combination of its evaluations at all exterior orientations y_e=+/-1/2. Every such orientation is a legitimate fixed-cut attachment. Therefore each extracted row is in the actual attachment-response span, even though it is not itself a single attachment count.

For any edge set W, its coefficient in the full H_N is zero unless W is a linear forest and every nontrivial component has an even number of edges. Otherwise its coefficient is

`s(W) 2^k c!/2^(c-1)`,

where k is its number of nontrivial components and c its number of all components, including isolates. To prove this, contract the selected paths. Their c blocks have c! orders. A path of odd edge length cancels under reversal; an even path has two equal-sign orientations. Each of the c-1 joins contributes1/2. This yields the formula.

The discovery search supplied `rank_contexts.txt`, a frozen bank of571 exterior coefficient contexts. SHA-256: `19f3e917d64e15d6cc79d15bf11d5f1fa3d3cac78215096226715cb21636987f`. These contexts use0..12 new vertices and0..16 selected exterior edges. The bank is finite; no random search is required to replay it. Its rows are **coefficient-extraction contexts**, not a claim that only571 individual attachment evaluations were executed.

`verify_exterior_rank.py` generates its own canonical representative of every new signature directly from the full old forest census, independently reconstructs every row from the coefficient formula, and uses rectangular modular Gaussian elimination. It finds:

- modulo101: rank571, determinant40;
- modulo103: rank571, determinant84;
- deleting the final row: rank570 in each field;
- replacing the final row by a duplicate of the first: rank570 in each field.

The checker also rejects old-old edges, duplicate edges and loops in an alleged exterior context. It independently checks the closed coefficient formula by enumerating all Hamiltonian vertex words and all edge monomials for N=1,...,5. Normal and `python -O` outputs are byte-identical, with1680 explicit checks.

A nonzero determinant modulo either odd prime implies a nonzero rational determinant, since all denominators are powers of2. Thus the actual attachment span has dimension at least571. The all-order compression theorem gives at most571 at n=6. Equality follows. For normalized laws the minimum number of nonconstant coordinates is570. This certificate does not establish the same minimum using only S and named B commands, nor all-order minimality.

## 8. Exact verification and reproduction

The independent audit imports no discovery-checker code.

1. `audit_carrier.py` independently classifies every undirected edge subset at orders0..6, including all32768 at n6. It checks both counts, every partition into monomial signatures, and direct substitution of every possible one-birth cut at orders2..5. It checks20540 enlarged basis pullbacks and142618 explicit requirements. Checks use raising functions, not Python assert.
2. Normal and `python -O` outputs are byte-identical: `audit_normal.json`, `audit_optimized.json`.
3. `actual_responses.cpp` independently counts actual Hamiltonian words by subset dynamic programming on589824 actual Boolean graphs:18 fixed contexts on all32768 six-vertex sources. It transforms each complete response row by exact integer Walsh transform. All forbidden monomials and every signed new-signature relation pass579565 requirements. The output is `actual_responses.json`.
4. `verify_exterior_rank.py` independently reconstructs and verifies the571-dimensional certificate as detailed above. It requires NumPy and its sibling `audit_carrier.py`, but no discovery-checker representatives or generated JSON as input. `verify_egf.py` independently checks the exact EGF coefficients through order12 using rational arithmetic and rejects the formula obtained by omitting the no-endpoint-pair boundary subtraction.
5. `operator_hostiles.py` independently enumerates directed words on all8 three-vertex sources, checks every one-birth cut, and records the adaptive/conditioning witnesses. Normal and `-O` outputs are byte-identical.

Reproduce inside this directory:

```sh
python -B audit_carrier.py > audit_normal.json
python -B -O audit_carrier.py > audit_optimized.json
cmp audit_normal.json audit_optimized.json
python -B operator_hostiles.py > operator_hostiles_normal.json
python -B -O operator_hostiles.py > operator_hostiles_optimized.json
cmp operator_hostiles_normal.json operator_hostiles_optimized.json
python -B verify_exterior_rank.py > rank_normal.json
python -B -O verify_exterior_rank.py > rank_optimized.json
cmp rank_normal.json rank_optimized.json
python -B verify_egf.py > egf_normal.json
python -B -O verify_egf.py > egf_optimized.json
cmp egf_normal.json egf_optimized.json
g++ -O2 -std=c++17 actual_responses.cpp -o actual_responses
./actual_responses coefficient_map.txt > actual_responses.json
```

`coefficient_map.txt` is the explicit full list of1486 forest masks, their new group indices and canonical signs, generated from the independently enumerated basis. The C++ evaluator never substitutes the quotient to calculate H; it uses actual adjacency matrices. Its Fourier tests depend on the map only for the relations being tested, not for the Hamiltonian counts.

## Connection boundary

There is a genuine boundary-state/gluing mechanism: contract selected paths, retain their endpoints and reversal parity, and compare contexts. This is a signed linear-observation quotient. It does not preserve actual internal directed paths, embeddings, knot crossings, or arbitrary source-conditioned behavior. A knot/tangle transfer would need a specified tangle composition, boundary permutation, and crossing/framing data before this theorem could be exported. Neither the count571 nor a shared formal recurrence supplies a knot theorem, Ramsey bound, or arithmetic realization claim.
