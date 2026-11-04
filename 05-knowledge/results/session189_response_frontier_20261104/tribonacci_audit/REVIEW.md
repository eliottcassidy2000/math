# Independent review of the sole-exit tribonacci tower and finite root-response coverage

## Verdict and exact scope

PASS for the all-order tower law in the targeted source PROOF.md at commit 7002daf0f2b97257df431a3805715791b9283807, for the independently reproduced finite coverage below, and for the ten explicit ready witnesses. INCOMING.json pins every targeted upstream file by URL, Git blob and SHA-256.

The finite root-response computations are independently recalculated, not merely replayed from the incoming triple files. They reproduce 1,782,468 rooted bases and 1,728,779 distinct (h,e,B,strong) rows over the retained ordinary tournament representatives of orders 4–9. The source-defined union attains all 351,061 odd values in [69,702189]. In fact every one has a retained STRONG-base construction, yielding a strong tournament of order at most 19.

This is finite coverage within these particular operators. It is not an infinite selector, a new global record for known tournament-spectrum coverage, an asymptotic density theorem, a closure of the two-digit joint-port grammar, or a solution of the all-strong name-retaining +4 branches.

## 1. The all-order three-state law

For a nonempty rooted tournament(T,t), B_t counts pairs of a possibly empty path and a second path starting at t, jointly covering all names. Flattening a Hamilton path and cutting immediately before t injects all h(T) Hamilton paths into these pairs, so B_t≥h(T).

Append x with its unique exit x→t. A Hamilton path either ends at x, giving h−e_t possibilities, or uses x→t and is a unique counted pair with x inserted between its paths, giving B_t possibilities. Thus

    h'=h−e+B,   e'=h−e.

For the new rooted response B_x, a second path starting at x either continues to t, recovering B_t, or is the singleton x, leaving an arbitrary Hamilton path of the whole old tournament. Hence B'=B+h. Every case includes the correct empty-prefix convention and preserves all physical names.

The state vector(h,e,B) therefore evolves by the integer matrix

    M=((1,-1,1),(1,-1,0),(1,0,1)).

Direct integer multiplication gives M³=M²+M+I. This proves the tribonacci recurrence at every height, with h_0=h, h_1=h−e+B and h_2=h+2B. The incoming A,C,D initial sequences and D_k=C_(k+2) then follow uniquely. No finite interpolation, positivity assumption on a formal inverse, or marked-cycle bijection is needed.

The law is valid for every root, not only a sole exit. At a sole exit a→r, subtracting the two one-exit formulas gives B_a(T)−h(T)=e_r(T−a), exactly the readiness digit. This identity does not turn an arbitrary root into a joint-port-ready root.

## 2. Near-join and strongness, checked independently

Let U=T→TT_k with k≥1 and keep root t in T. Ordered concatenation gives h(U)=h(T), while e_t(U)=0. In the definition of B_t(U), the k added transitive names may independently be put in either path; each induced transitive suborder has one Hamilton path. Therefore B_t(U)=2^k B_t(T). Appending a sole-exit z→t gives

    H(T→TT_k→z, with z→t)=h+2^k B.

This proves the near-join law directly from the tower lemma, without relying on an uninspected external proof.

If the base T is strong and has at least three names, adding a sole-exit vertex preserves strongness: the new vertex reaches T through its root, and every old name reaches an old name different from the root and hence the new vertex. Thus every positive-height tower over the strong bases used here is strong. The near-join is also strong: z reaches all of T, T reaches every chain name, and every chain name reaches z. The single-vertex base at tower height one is not strong and is not used for the finite strong-base assertion.

For unbounded arbitrary rooted bases, the near-join value h+2^kB can also be written as a height-two value over T→TT_(k−1), whose response is2^(k−1)B. This identifies those two numerical families when base order is unrestricted. It does not preserve a strong base or a sole-exit/joint-port guard, and does not prove equivalence with the older strong sole-exit ready-pair selector.

## 3. Independent finite construction bank

The physical inputs are retained as column-major tournament codes, together with the source census's self-converse flags. They come from the already reproduced pinned deletion census, whose full upstream provenance is recorded in CODES.json. Both orientations of each non-self-converse representative are used. The finite ATTAINMENT result does not require any assertion of census completeness: every code is an actual tournament, and one independently checked construction witness is supplied for every claimed output.

coverage.cpp uses a predecessor endpoint recurrence, independently of the incoming successor implementation. It computes h and e, all induced-subset Hamilton counts, the exact root response B, and strongness for every physical base. The observed ordinary class counts 4,12,56,456,6880,191536 and all root-response row counts agree with the incoming source. Distinctness counts and first missing values refer to this finite input bank and the specified height ranges, not to all tournaments at all orders.

The exact reproduced intervals and value-set sizes are:

- Height two, all bases: [69,15449],8,637 values.
- Height two, strong bases only: [69,15449],8,594 values.
- Near-join parameters1–8: [69,382687],674,533 values.
- Tower heights3–10: [69,77],801,304 values over its complete range.
- Tower heights1–10: [69,432485],801,327 values.
- Source-defined union of height two, near-join1–8, and tower3–10: [69,702189],832,898 values.
- The same source-defined union using only strong bases: [69,702189],832,195 values.

All lower endpoints69 are explicitly tested. The next missing value for the union is702191 within this finite bank and operator range.

The incoming coverage.py starts its interval at 69 without checking membership of 69. This audit repairs that validation gap; it does not find a missing endpoint. The script also omits height one, although a README/table row describes heights1–10. Adding height one contributes exactly five extra low values 1,11,13,47,67 to the height-two-plus-higher-tower set. It leaves the stated upper interval unchanged. Accordingly the source-defined union has the five listed low exceptions in addition to7 and 21; adding height one gives 832,903 values and removes those five exceptions. These bookkeeping differences do not invalidate the reported positive interval.

## 4. A strong witness for every odd value, not only a set cardinality

WITNESSES.json specifies the packed format. Across its ordered parts, line i supplies the physical base code, root, converse choice, operation, parameter, and strong flag for H=69+2i. Every one of the351,061 lines uses a strong base.

witness_check.cpp independently decodes and regroups these records by physical base, checks strongness, and recomputes h,e,B with a SUCCESSOR recurrence. It verifies the chosen formula at every line. This is a separate direction of endpoint DP from the coverage generator. There are 44,306 distinct physical base graphs among the chosen witnesses. Every constructed output has at most 19 names. No isolated missing output is inferred away by density or a count comparison.

The proof in Sections1–2 certifies the actual graph construction for every witness. In addition, six concrete resulting tournaments are built and counted directly, including H=69,12875,12877,38629,38635 and 702189; their physical rows, counts and strongness are retained. These numerical realisers do not imply a legal H-neutral rearming from a prescribed lineage state or preserve its name budget.

## 5. Tower controls and explicit order-fourteen ready witnesses

The separate physical counter checks 2,859 actual graphs: 321 rooted base/tower controls, 963 near-joins, the explicit re-rooted tower sequence, and ten supplied ready memberships with deletions and selected actual C2 children. The base controls include every labeled rooted tournament through order four and deterministic additional bases through order seven. All small-base(h,e,B) values are also checked by literal flattened permutations with a cut immediately before the root.

The re-rooted tower family gives the stated digit-one sequence5,13,31,77,193,489,1255,3253,8493,22293 at orders 4–13. This direct construction requires no classification of all ready classes.

Both claimed strong order-fourteen digit-one realisers are independently verified:

- H=12875, sole exit0→10, digit e_10(T−0)=1.
- H=12877, sole exit0→10, digit e_10(T−0)=1.

Their actual two-birth children have order 16 and Hamilton counts 38627 and 38633 respectively. The other six order-fourteen memberships and the order-twelve digit-two 12877 and order-thirteen digit-one38615 witnesses also pass. Complete rows and exact deletion checks are retained in membership_controls.json.

These positive membership checks do not independently certify the order-fourteen ready generator's completeness, floor, or negative membership claims. Nor do digit-one witnesses give the +4 children38629 and 38635 from those same leaf states. The all-strong name-retaining order-fifteen digit-two rearmings remain open. The pinned7002daf0 source had only the smaller digit-two 12877 witness. A later targeted2cfc7bba update supplies a strong order-thirteen digit-two 12875 witness; this has been independently checked in the separate two-cycle comparison packet, so both +4 branches now have the stated weak-allowed route. This later positive update does not turn transitive-prefix padding into a strong name-retaining preparation.

## 6. The earlier marked-bijection correction still matters

The latest ready README continues to state the overly broad marked correspondence from THM-4862. Its four-name counterexample and corrected image condition have not been removed by the new tower theorem. The exact marked image must require that deleting the marked sole exit leaves a strong tournament. The canonical admissible-mark argument repairs the unmarked correspondence for digits1 and 2 at orders at least four, as proved in the separately retained repair packet.

The present tower and near-join laws, the finite root-response bank, and the explicit membership tests do not use that disputed marked correspondence. Therefore the earlier correction does not invalidate these accepted results. Conversely, this review does not silently approve the unqualified marked statement, external enumeration hypotheses, or an all-order ready-extension heredity theorem.

The new all-rooted selector is a broader candidate class than the older strong sole-exit ready-pair selector. The incoming note's description of the two as the same conjecture is not supplied with a reduction preserving the ready guard or its strong/name constraints. No such equivalence is used here. In particular, the702189 window concerns a union of three operations; the finite height-two selector alone is verified only through 15449.

## 7. Release and limitations

Normal and optimized Python run the same explicit checks. All retained upstream bytes, independent sources, input code parts, complete witness parts, actual graph controls and receipts are covered by PUBLICATION_FILES.json and SHA256SUMS. Executables and concatenated duplicates are excluded. Standard-library Python and installed g++ suffice to reconstruct the entire finite bank.

No source file was changed, no repository-wide scan was performed, and no package installation was needed. The audit's mathematical scope is exactly the tower/near-join identities, finite construction coverage, and positive named witnesses described above. There is no infinite-density or global-closure inference.
