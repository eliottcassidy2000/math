# Independent audit: Paley palettes and rounded carries

Audit target: `../session23b-palettes/PALEY_MULTICOLORING_AND_CARRY.md`.

Verdict: **PASS for the mathematical claims as stated.** No substantive mathematical defect was found. The proof is self-contained and correctly distinguishes the rational palette-budget series from the nonrational depth series. The finite controls below supplement, rather than replace, the general arguments.

## Exact palette law and its generating function

For the Paley tournament with difference set `{1,2,4}`, an independent source-removal test finds 7 legal singleton supports, 21 pairs, 21 transitive triples, and no larger legal support. In particular the maximum legal support size is 3.

The three displayed base palettes are correct, including the deliberately repeated support `{0,1,2}` in the 2-fold palette. Each repeated occurrence is a different color. Exhaustive breadth-first search of exact vertex-incidence states, allowing repeated legal supports, gives minimum budgets `0,3,5,7,10` for `k=0,1,2,3,4`. Direct verification of the displayed periodic construction succeeds for every `k=0,...,100`. Combined with the support-size lower bound, the proof establishes the law for every nonnegative integer k, not merely those checked.

Exact polynomial arithmetic confirms

`sum chi_k z^k = (3z+2z^2+2z^3)/((1-z)(1-z^3))`.

Numerator and denominator are coprime. At a primitive cube root w, the numerator is w, which is nonzero; at 1 it is 7. The reduced denominator degree is 4. Thus the claimed minimal eventual constant-coefficient recurrence order is exactly 4. The written recurrence holds already for k>=4.

## Substitution

The occupied-block support criterion is correct. Necessity follows by choosing one representative per nonempty block. Sufficiency follows by ordering the occupied base vertices transitively, then concatenating each block's transitive order. This also proves the general identity `chi_d(T[U]) = chi_{chi_d(U)}(T)` for nonempty finite tournaments T and U, with k-fold coloring defined by distinct colors at each vertex and repeated supports permitted.

The lower-bound construction correctly chooses any k distinct colors used in each U-block; discarding other incidences preserves legality. The upper-bound construction correctly relabels a fixed k-coloring separately in each block. Therefore the proposed specialization to P7 and the value `chi_d(P7[P7])=7` are sound, without enumerating the 49-vertex product.

Independent computational control used all 121 ordered pairs of labeled tournaments of orders 1, 2 and 3. For each pair it compared exact dichromatic coloring against exact multicoloring and checked every induced subset's support criterion. It also checked every induced subset of `P7[K2]`. Total: **51,410 induced-subset checks**. No large lexicographic power was enumerated.

The fractional value is correct under the intended fractional dichromatic interpretation. For maximal clarity, define this explicitly or write `chi_{d,f}(P7)=7/3` rather than leaving `chi_f` potentially confused with fractional coloring of the underlying graph.

## Nonrationality of the depth sequence

The general argument for coprime A>B>1 and positive integer x_0 is valid.

1. Rationality over Q implies an eventual homogeneous recurrence over Q for x. Its linear combination `e_n=B*x_(n+1)-A*x_n` then obeys an eventual homogeneous recurrence too.
2. The bounded integer sequence e has a finite set of possible recurrence-length windows. After the recurrence begins, the next window is determined by the current window. A repeated window forces eventual periodicity. This proof does not need the recurrence coefficients to be integers or a theorem about algebraic-integer roots.
3. If the carry has eventual period L, iteration over one period yields `x_(n+L)=q^L*x_n+c` with q=A/B and fixed rational c on a given residue class. Its fixed point is `alpha=c/(1-q^L)`. For this ceiling recurrence, c>=0 and hence alpha<=0.
4. Choose an integer M clearing alpha's denominator and put `N=M*(x_n-alpha)`. Then `M*(x_(n+tL)-alpha)=A^(Lt)*N/B^(Lt)` is an integer for every t. Since gcd(A,B)=1, every B^(Lt) divides N. Thus N=0. In fact this directly contradicts x_n>0 and alpha<=0; alternatively it gives the boundedness contradiction used by the candidate.

The positive initial value and reduced denominator greater than 1 matter: x_0=0 is the zero sequence, while an integer ratio has zero carry and a rational geometric generating function. The stated theorem correctly excludes both exceptions.

Independent rounded-sequence controls tested **432 parameter/initial-value cases**, using all coprime A>B>1 with A<=20 and x_0 in {1,2,3,4}, for 80 steps each. They check carry bounds, strict growth, the normalized telescoping identity, and every 1-through-8-step affine identity on the tested blocks. These are finite checks of identities, not a purported numerical proof of aperiodicity.

For the Paley ratio, the initial values agree with the candidate. The limiting normalized value exists because the nonnegative increments sum absolutely, and its bound `1 <= limit <= 3/2` is correct. The 200-step computation gives approximately **1.3658806634197658**, with an exact rational certified interval in `audit_results.json`; this decimal is not an exact-value, algebraicity, or transcendence claim.

Only nonrationality is audited here. No non-D-finiteness claim is being added.

## Fano, distinguishing colors, and transfer boundaries

The seven translated out-neighborhoods `{i+1,i+2,i+4}` form a Fano triple system: each unordered pair occurs once. Every one is a directed triangle. Across all 128 binary vertex colorings, the counts of monochromatic Fano lines have histogram `{1:70, 2:42, 4:14, 7:2}`. In particular every binary coloring has at least one.

Enumeration of all 7! relabelings finds 21 automorphisms of P7. A coloring with black set `{0,1}` and all other vertices white has trivial color-preserving stabilizer, verifying distinguishing number 2. A short direct proof is also available: the unique black edge 0->1 fixes both endpoints; their adjacency signatures fix vertices 2,4,6 and leave only the pair {3,5}; adjacency to vertex 2 separates that remaining pair. Nontrivial translations rule out distinguishing number 1.

The distinction among acyclic, proper and distinguishing colorings is maintained. No proper plane unit-distance lower bound follows. No native-source, guard, clock, arrival-budget or other Collatz/LRC structure is supplied, so no such theorem transfer follows. For maximum precision, any informal reference to a finite-state mechanism should mean an autonomous finite-state mechanism; a finite-state transducer supplied a nonperiodic external input is a different notion.

## Literature provenance check

The referenced primary records have the correct subjects. The publisher record for V. Neumann-Lara's 1994 paper is *The 3 and 4-dichromatic tournaments of minimum order*, Discrete Mathematics 135, 233-243, DOI `10.1016/0012-365X(93)E0113-I`. Its abstract identifies the minimum orders and characterizes the minimal examples. This audit did not inspect the complete proof.

https://www.sciencedirect.com/science/article/pii/0012365X93E0113I

A more directly accessible primary source for the product theorem is Costa and Silva's author preprint, *Acyclic Coloring of Products of Digraphs and of Digraphs with Bounded Treewidth* (2022), arXiv:2204.14212. Its abstract explicitly states the replacement of H in G[H] by the complete bidirected digraph on chi_a(H) vertices. The 2024 journal record is *Acyclic coloring of products of digraphs*, Discrete Applied Mathematics 349, 59-69, DOI `10.1016/j.dam.2024.01.020`.

https://arxiv.org/abs/2204.14212

https://flore.unifi.it/handle/2158/1358211

https://www.sciencedirect.com/science/article/abs/pii/S0166218X24000283

These abstract/metadata checks substantiate the provenance disclaimer, not priority or a complete literature audit. Direct initial publisher-page opens returned HTTP 403; primary indexed publisher material and the accessible author/institutional records supplied the bibliographic checks. No source-derived proof is needed for the candidate's exact certificates or carry argument.

## Reproduction and scope

Run `python independent_audit.py` in this directory. It writes `audit_results.json` and prints a concise summary. Python standard-library exact arithmetic is used throughout; SymPy is used only for formal polynomial/series verification.

Files written only in `session23b-palettes-review`; no repository files were changed and no Git operation was performed.

## Final revision binding and optimized-mode rerun

The final proof was reread after its status update and the three optional precision edits: explicit fractional dichromatic definition, the word “autonomous” in the finite-state statement, and the author-preprint link. The displayed formulas and mathematical arguments remain unchanged. The optional precision recommendations above are now resolved. **Final verdict: PASS.**

Final proof SHA-256:
`fe009bd8cdb5a7ab03fb8758a7c969594295ff403e7c6da910378847499afa2d`

The independent checker now uses explicit `need` calls instead of Python `assert` statements, so all checks remain active under optimization. Both `python independent_audit.py` and `python -O independent_audit.py` pass **329,752 exact checks** and produce byte-identical `normal.json` and `optimized.json` receipts. The checker pins the final proof hash above and records its own source hash. Independent checker SHA-256:
`429c6aa540088169cc2ef4f8467a67e5812f647e5a993ccc4363d2b43b9ad81e`

The root's `check_palettes.py` was also independently rerun normally and under `-O`. Both runs pass **2,266 checks** and produce byte-identical receipts, also matching the root's existing receipt. Copies are `root_checker_normal.json` and `root_checker_optimized.json` in this review directory. Root checker SHA-256:
`ac510c6f287d22d84678e1df78d1ecab677a4cd7386e23ff3aa4b817d1efbaf0`

`final_binding.json` binds this report, the proof, checker sources, and normal/optimized receipts by SHA-256. This revision made no changes outside the review directory and performed no Git operation.
