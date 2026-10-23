# Session41 synthesis: compatible odd-cell attachments and both Fibonacci branches

23 October2026. PROVED statements below have declared carriers. COMPLETE FINITE tests are not promoted to unrestricted tournament classification or an infinite Goldbach/Collatz result. The session follows the user's instruction to investigate both ends of discrepancies and identify larger universes when possible.

## 1. A shared chart is stronger than separate compatible attachments

The [actual attachment criterion](attachments/SHARED_GAUGE_CRITERION.md) compares two named cellular tournaments using one root internal isomorphism, transported along a spanning tree. A candidate survives exactly when:

1. every transported map is an allowed internal isomorphism preserving all retained marks;
2. all chord equations intertwine the two actual holonomy tuples;
3. all remaining block templates, polarities, singleton incidences and actual path flags agree.

The word “one” matters: independently selecting a valid chart on each edge can give inconsistent choices at the same middle cell. On actual cyclic triples, self-comparison may simplify internal checks because C3 is normal in S3, but arbitrary two-system comparison must retain tree chirality. Extra decorations can shrink the allowed domain further.

Two concrete hostiles are certified:
- A three-R5-cell TREE has individually equivalent attachment double-cosets but no simultaneous compatible middle chart. Thus checking only cycle holonomy is insufficient once the internal chart domain matters.
- In a native order31 example, twisting an unused triple-chord block preserves the supplied actual Hamiltonian path, regular scores, endpoint refinement and each local block orbit. The common-outneighbor spectrum changes from465 pairs of value7 to393 of7,36 of6 and36 of8. The tournaments are nonisomorphic; the common-root test rejects the false assembly.

The complete finite carrier check covers930 two-mark cases of the particular D(D(P7)) source,933 actual construction controls and1517 source-automorphic alignments. This is not a census of every order31 tournament or every Hamiltonian path.

## 2. Removing the path clock changes the reconstruction problem

The [indexed flag theorem](attachments/ONE_OWNED_ARC_AT_ORDER31.md) resolves the remaining alignment in this structured source with one appropriately chosen INDEXED actual path arc: its position and both vertex owners are retained.

Dropping the position is not harmless. Native control697 has21 source-automorphic candidate paths and no occurrence-unique single edge. Exactly13 pairs of noncontiguous actual edges isolate the path. A contiguous unindexed window needs five arcs, whereas one indexed arc suffices. The two-owner observation is not claimed to cost only one bit.

The [infinite doubled-Paley family](native-arc-family/INFINITE_ARC_FLAG_GAP.md) proves the same distinction for every prime q≡3 mod4, q≥7, inside an explicit orbit of q(q−1)/2 Hamiltonian paths of D(D(P_q)). One indexed flag identifies the path; two unindexed arcs suffice and one never does. The [window theorem](window-extension/WINDOW_THEOREM.md) proves five consecutive unindexed arcs are necessary and sufficient in that orbit.

These positive assertions are orbit-relative, not claims about every Hamiltonian path outside the specified orbit. Their negative assertions are actual native examples, not arbitrary feasible-margin counterfeits.

## 3. Residual-cycle potentials need the correct coupled integer fibre

The [integer-trade theorem](residual-cycles/README.md) gives an all-size exact criterion for a supplied bounded integer model Ax=b, l≤x≤u. A linear projected wiring observation is constant on that fibre if and only if every conformally primitive kernel move feasible at one base point has zero projection.

Proof mechanism: every integer difference decomposes conformally into primitive moves. Each summand is individually feasible at the base point because its coordinates lie between the base and the target. Therefore any changed observation has a feasible primitive witness. Ordinary network cycles are the special case for a network matrix, not a substitute for the primitive moves after shared-owner equalities are added.

In the smallest coupled-table example, neither individual2×2 switch is allowed, but their simultaneous switch is. A test restricted to separately feasible residual cycles would falsely declare uniqueness. The same coupling can also enforce cancellation for another projection. Fixed labelled mod8 linear guards can be retained using bounded integer auxiliaries; nonlinear electrical invariants are not declared linear by analogy.

An executable finite-domain compiler retains shared cell/frame variables rather than multiplying independent marginal choices. Controls include33750 basepoint tests,32768 encoding assignments and324 actual six-vertex frame comparisons.

The [owner-indexed generating functions](residual-cycles/OWNER_INDEXED_ENUMERATORS.md) give a cheaper exact method on genuine forests or a supplied tree decomposition: contract child messages over the SAME separator state. This counts compatible choices and certifies projected min/max equality. It does not assume that an arbitrary native support has small treewidth.

## 4. The Fibonacci discrepancy now has two distinct resolutions

The [Fibonacci grammar](fibonacci/RESULTS.md) uses A=F_-1, B=F_1, C=F_2. The zero slot is forced unselected.

- D deletes the zero position and reconnects A,B,C into consecutive slots.
- R retains its position, so A and B are nonconsecutive and A+B is legal.

Deleting zero preserves lengths but removes an exact representation sheet. The complete D fibre consists of the applicable supports Z(N), A+Z(N−1), B+Z(N−1); R adds A+B+Z(N−2) when the ordinary unit is absent from that tail. Proper self-decomposition of F_n is not unique: R has two proper decompositions for every n≥3, whereas D has one for odd n and two for even n.

Selecting the highest-index proper support and recursively expanding gives the canonical substitution

C→CA, A→B, B→CA.

Its row R_N is obtained from R_(N−1) by adding exactly one unit on the left. It matches the current16-letter prefix, but gives C at35 where the older supplied stripe records B. That historical value is retained, not silently corrected.

The [branch classification](fibonacci/BRANCH_ADDENDUM.md) resolves what each model permits:
- With fixed canonical blocks, old35B fails at36 in D. R repairs36 but fails37.
- With independent proper expansion choices per occurrence, synchronized in their visible absolute-position colors, the complete D family permits canonical B→A/B.
- The complete R family additionally permits canonical C→B/C. Canonical A remains A.

Thus old35B has an infinite coherent continuation in R;36 is forced A and37 can be A or B. The canonical branch also remains valid. This is a genuine encompassing universe, not a choice to discard one end.

A SINGLE globally reused noncanonical expansion dictionary, or a shared parse tree across all rows, is a stronger synchronization requirement and remains separate. Occurrence-dependent compatibility must not be relabelled as solving that problem.

## 5. Exact grid counts, branch generating functions and aperiodicity

The [diagonal-grid note](COLORED_GRID_AND_PHASE_QUOTIENT.md) groups CA into a long tile and B into a short tile, recovering the binary Fibonacci substitution. With alpha=phi^-2, the exact canonical prefix counts are

c(H)=floor((H+2)alpha), a(H)=floor((H+1)alpha), b(H)=H−a(H)−c(H).

The kth positions are C_k=floor(k phi²)−1, A_k=floor(k phi²), and B_k=2floor(k phi)+k. These formulas are proved, not fitted.

Put S(H)=sum_(m=1)^H floor(m alpha). In the entire triangular grid through row H,

T_C=S(H+2), T_A=S(H+1), T_B=H(H+1)/2−T_C−T_A.

The exact recurrence with K=floor(H alpha),

S(H)=K(H+1)−3K(K+1)/2+S(K),

gives logarithmic-depth counting with integer-square-root arithmetic, tested on501-digit inputs. The boundary excess T_C−T_A equals c(H).

The [branch partition functions](BRANCH_PARTITION_FUNCTIONS.md) are

Z_D=x_A^a x_C^c (x_A+x_B)^b,
Z_R=x_A^a (x_A+x_B)^b (x_B+x_C)^c.

Hence globally extendable prefix counts are2^b and2^(b+c), with exponential rates phi^-3 log2 and phi^-1 log2. Pinning old35B halves the retained prefix count after35 but changes neither rate.

All attainable limiting frequency triples are described explicitly. Nevertheless every branch's A-indicator is aperiodic and has a nonrational ordinary generating function: forced-A and forced-not-A Beatty positions both visit every residue class infinitely often. This remains true even for branches with rational frequency limits or no frequency limits.

Retain actual horizons: the first5000 canonical A occurrences hit only95 of144 residues modulo144; first complete coverage occurs at occurrence11001. An asymptotic density argument is not a finite clock certificate.

## 6. A real parity/phase dictionary, with its lost information exposed

The [binary phase dictionary](BINARY_PHASE_DICTIONARY.md) encodes the three F3 owners as the three nonzero F2² vectors. Translation becomes T=[[1,1],[1,0]], reflection becomes [[1,1],[0,1]], giving the exact Aff(F3)≅GL2(F2) action.

Connected parallel-owner constraints become homogeneous binary linear equations with the ZERO solution excluded. Nullity d≤2 gives2^d−1 owner sections. A system can have no nonzero fixed vector but still have a three-element commuting symmetry group: fixed sections and centralizers are different carriers.

The same T occurs as a quotient of the Fibonacci incidence matrix, but original letters C and B have already been identified there. It is not a cyclic permutation of the original three colors. The dictionary is multiplicative in a matrix group and therefore does not contradict the absence of additive C3 maps into a binary group.

T has exact order3·2^(s−1) modulo2^s. This clock, the Pascal path clock, and the native electrical mod8 guard are distinct objects.

## 7. Additive number theory: both a positive bridge and precise obstructions

The [additive packet](additive/RESULTS.md) retains the inherited charge

b(n)=floor((n+1)alpha), Q(n)=(n−2b(n),b(n)).

It proves the exact all-arity carry range and addition-tree independence. Shared-support generating functions use2^k digit-subset states for k summands, or k+1 states when supports must be disjoint. Carry zero is not disjointness:3+3 is the smallest positive overlap-with-zero-carry hostile.

A common quadratic-algebra family contains both golden charge and Gaussian norms. Their differences matter:
- Modulo2 the golden algebra is the field F4, while the Gaussian algebra has a nonzero nilpotent. They cannot be identified as rings.
- Modulo3 the explicit map t→2+i IS an isomorphism preserving multiplication, conjugation and norm. The obstruction is prime-local, not an absolute rejection of the analogy.
- Numerical evaluation Q(n)→n is not multiplicative in the golden ring; its exact product defect is retained.

There is also an actual infinite norm bridge: F_(2k+1)=F_k²+F_(k+1)². The associated Euclid triple has content2 exactly for k≡1 mod3 and content1 otherwise, with its degenerate initial case excluded from triangle claims.

Finite tests partition every even Goldbach target through2000 by actual canonical owner colors and carries. They prove no uniform prime-pair lower bound. Same owner color plus the same residue modulo4 still fails to determine primality or sum-of-two-squares membership; explicit hostiles are retained.

The larger stripe family allows arbitrary subsets of canonical C positions to remain C. Retaining C precisely at the prime canonical-C positions would therefore be a definition, not a discovery. An arithmetic claim must fix its selector independently or quantify over the full family. Both ends of the discrepancy have meaningful uses: the canonical branch supplies a rigid arithmetic register, while the larger family exposes what grammar alone cannot force.

## 8. Recent repo work is inherited, not duplicated

Incoming0429c7d8fdf and508cd019aa0 are preserved. Their [odd-cell gluing](../odd_cell_gluing_20261023/README.md) supplies functional partial-chart propagation and chart-conditioned residual potentials. Their [zero-slot note](../odd_cell_gluing_20261023/ZERO-SLOT.md) already proves the four-sheet fibre and arbitrary C→B repair of35; this session credits those convergent results and sharpens them to the complete language, necessity, counts and frequency region.

Their affine group-pair clock retains unbounded prime-power translation memory after an eight-step slope return. It is not a period of the color word or a Collatz trajectory. The same companion-matrix grammar can act on different carriers without identifying their meanings.

### Late incoming completion at 3f1a9a2c

The [late-source integration review](audit-late-incoming/INTEGRATION.md) incorporates [odd-cell completion](../odd_cell_completion_20261023/README.md). On a supplied named equitable partition into cells of sizes 1, 3, 5 and 7 with actual blocks, source-only canonicalization is now available, with at most 21 root charts per seven-cell informative component. A sufficient actual-owner marking budget is c3+c5+2c7, counting isolated nonsingleton cells as components; it need not be minimal. This genuinely advances beyond comparing two supplied arrays, but does not discover such a bounded partition or recover unknown native path data.

For fixed local block types with full internal groups C_p, the exact joint alignment coordinate is W/image(D), not a collection of independent edge orbits. Its dimension on a connected network is beta+r−1+a, where beta is cycle rank, r counts rigid edges, and a records the all-symmetric balanced case. Rigid tree attachments can therefore retain alignment information even without cycles. This finite-field quotient is distinct from our integer residual-move lattice.

The inherited network Burnside count retains named cells and full C_p internal groups, without converse or cell-permutation quotient. Nine-cell and prime-eleven hostiles prevent extending these conclusions from connectedness or prime size alone. The integration review checks proof scope; its numerical censuses were not independently rerun here.

The constructive combined route is to authenticate shared charts, retain joint alignment coordinates, and only then compile compatible states for coupled integer-move certificates or owner-indexed contraction. Unknown larger-cell kernels and lossless native-port observers remain substantive obligations.

## 9. Remaining obligations

- Prove that a compact attachment model exhausts the actual native port choices at arbitrary order, rather than merely classifying a supplied finite-domain model.
- Retain common cell charts, actual path positions, converse and electrical guards when moving from local certificates to full wiring.
- Classify the stronger globally synchronized noncanonical Fibonacci dictionary model separately from occurrence-dependent choices.
- For additive questions, retain shared digits and arithmetic filters, then prove a genuinely uniform admissible representation bound. The current formulas and finite prime tests do not prove Goldbach or Collatz.

Independent audits accompany the first, second, third, fourth and fifth packets. The native arc family and window theorem have standalone witnesses and explicit orbit scope. All executable receipts distinguish proof, exhaustive finite testing, and inherited results.
