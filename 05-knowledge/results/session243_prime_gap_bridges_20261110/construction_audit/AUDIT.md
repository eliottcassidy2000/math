# Independent hostile audit: physical prime-gap bridge constructions

**Verdict: PASS for the explicitly scoped mathematical statements below, with inheritance and boundary corrections communicated before publication.** The final producer proofs and both compact canon statements were read and accepted; their inspected bytes are pinned in SOURCES.json. No producer verifier was imported or used to generate the independent receipts. Chan–Pak completeness, Polignac, infinitude of prime values, and arbitrary target placement remain OPEN.

## 1. Scope and decisive findings

1. The four literal rooted sources do give the synchronized strong outputs
   `6Q-H, 6Q+H, 6Q+3H, 6Q+5H` for every nonempty tournament gadget with pair `(H,Q)`. The common difference is exactly `2H`.
2. For a transitive gadget of order `m>=1`, these are the four consecutive odd numbers `6*2^m-1, +1, +3, +5`. The supplied strong witnesses have orders `m+4,m+4,m+5,m+5`, respectively. They are all of order **at most** `m+5`. Adding a universal sink equalizes order but destroys strongness; no same-order strong conclusion is justified by that padding.
3. The retained-root map `(a,b)->(Ha,Qb)` is correct. Its scalar product form is already the singleton-root specialization of Session227's path-cover diagonalization, and repeated retained-root grafting is literally one graft with an ordinal-sum gadget. This is inherited structure, not a new diagonalization theorem.
4. The fixed finite source bank obstruction is valid only for the genuine final intercept `a=1`. It does **not** forbid coincidental Euclid-number outputs. The strong source `[20,1,2,7,14]` rooted at `3` has `(a,b)=(-5,18)`; a singleton graft gives `31=-5+2*18=2*3*5+1`.
5. The compact binary compiler realizes every positive even difference `2N` as a **single-vertex extension** between two strong tournaments, with smaller count at most `3N`, larger count at most `5N`, and larger order at most `ceil(log2 N)+4`. The deletion identity is literal, not merely equality of counts.
6. The words `JOJOZ` and `JOOOZ` have total-count rows `(9,-3,9,0)` and `(9,-1,9,0)`. Applied to a strong source with endpoint count `E=N`, they give same-order strong outputs separated by `2N`. Increasing the near-transitive base in the inherited endpoint selector preserves `E=N` and gives unbounded pair counts. **Existence and infinitude of all even differences were already proved in Session199**, including a single-arc realization, and Session221 supplies a later logarithmic-order construction. This is an alternative endpoint-controlled mechanism, not a newly solved difference-spectrum question.

The compact extension theorem and endpoint infinitude theorem are separate: the bound `5N` is not claimed for infinitely many distinct pairs at a fixed `N`.

## 2. Independent reconstruction of the graft law

All paths are directed vertex words. Empty induced tournaments have Hamilton count one. Write the rooted source state as `(D,E,U,d)`, where `U` counts ordered spanning two-path covers with one designated path starting at the root and the other possibly empty. Swapping the two path bins does not change this convention's count.

Set

    f=U-D-d+E,   a=d-f=D+2d-E-U,   b=U-d=D-E+f.

The literal graft retains source root `r`, puts every gadget vertex before `r`, and puts all old nonroot vertices before every gadget vertex. Its scalar count is the inherited

    N=(D-E)H+dH+(U-d)(Q-H)=aH+bQ.

A path ending at `r` must be an old-deletion Hamilton path, then a gadget Hamilton path, then `r`; hence `E'=dH`. Deleting `r` leaves the ordinal sum `(A-r)=>B`, giving `d'=dH`.

For `U'`, covers whose root-starting old path is just `r` contribute `dH`. Each of the other `U-d=b` old covers allows the gadget to be distributed as the two terminal path pieces, in exactly `Q` ways. Thus

    (D',E',U',d')=(aH+bQ,dH,bQ+dH,dH).

Consequently `a'=Ha` and `b'=Qb`. The law includes sink-root sources, but assumes a nonempty gadget. For an empty gadget the original source is unchanged, whereas `a+b=D+d-E` need not equal `D`. The actual seed `(D,E,U,d)=(5,1,9,3)` has `D=5` but `a+b=7`.

The strongness criterion is inherited from Session226: the graft is strong iff the old source has at least two vertices and its root reaches all source vertices, independently of gadget connectivity. All four sources are strong, so the criterion applies to every nonempty gadget, including transitive ones.

An independent source response computation enumerates permutations of the old vertices and eligible gadget-run slots. Its coefficients `(c1,c2)` give `c1*H+c2*(Q-2H)`. For the four sources these are `(11,6),(13,6),(15,6),(17,6)`, independently recovering the four `(a,b)` pairs without using `U` or the proposed transfer row.

## 3. Literal sources and independent finite universe

Outgoing bit rows and zero-based roots:

| Rows | Root | `(D,E,U,d)` | `(a,b)` |
|---|---:|---|---|
| `[8,1,3,6]` | 3 | `(5,1,7,1)` | `(-1,6)` |
| `[8,1,3,6]` | 2 | `(5,1,9,3)` | `(1,6)` |
| `[20,1,2,7,14]` | 1 | `(9,3,9,3)` | `(3,6)` |
| `[20,1,2,7,14]` | 2 | `(9,3,11,5)` | `(5,6)` |

`verify_independent.py` uses a first-vertex recurrence, independent literal permutation/split enumeration, and the slot-count method above. Its exact universe is:

- All 285 rooted labelled sources of orders one through four, **every root**, against every labelled gadget of orders one through three: 3,135 literal full-state and strongness tests. No properness or strongness filter hides boundary cases.
- The four designated sources against all 1,099 labelled gadgets of orders one through five: 4,396 literal full-state tests, all exact common-gap identities, and strongness.
- Transitive gadgets of orders one through eight, plus all four seeds and eight small graft outputs checked by literal permutations of every subset.
- Three-step retained-root constructions using a singleton, directed triangle, and directed edge; their fused ordinal-sum graphs agree arc for arc, and their state products agree.
- Explicit empty-gadget, changed-root, sink-padding, and accidental-Euclid hostiles.
- `H=1` iff transitive, and its response `Q=2^n`, checked on all 1,099 labelled gadgets through order five. The all-order argument is independent of this finite check.

The output is `CHECK.json`; normal and optimized runs are identical. Counts, states, source rows, and consequence objects are recorded, rather than only injected formula values.

## 4. Genuine unit intercept versus growing primorial support

For a nonempty sequence of retained-root grafts with gadget pairs `(H_i,Q_i)`,

    D_final=a_0*product(H_i)+b_0*product(Q_i),
    a_final=a_0*product(H_i), b_final=b_0*product(Q_i).

If `a_final=1`, positivity and integrality of each `H_i` force `a_0=1` and every `H_i=1`. A nontransitive tournament contains a directed triangle; the elementary vertex-insertion injection shows its full Hamilton count is at least the triangle's three. Thus `H_i=1` forces a transitive gadget. Conversely a transitive tournament has a unique Hamilton path. Hence `Q_i=2^{m_i}` and

    D_final=1+b_0*2^M,   M=sum(m_i)>=1.

For a fixed finite source bank, the odd prime factors of `D_final-1` therefore lie in the union of the fixed prime supports of the nonzero `b_0`. Growing primorials eventually introduce a new prime, so they cannot all occur in this genuine-intercept-one language. A primitive-looking numerical formula is not a substitute for the final intercept condition.

There is an even sharper primorial gate. For `P_n=2*3*5*...*p_n`, equality `D_final=P_n+1` forces

    v2(b_0)+M=1.

Since `M>=1`, it follows that `b_0` is odd, `M=1`, and `P_n=2b_0`. Thus only one singleton gadget can have been used. A seed with `b_0=0` gives count one and no Euclid value. This refinement is valid only for the same genuine unit-intercept, nonempty-word scope. Root changes and refreshed source banks are excluded from the argument.

## 5. Compact one-vertex extension differences

Start with the rooted edge `(D,E,U,d)=(1,0,2,1)`. The fixed-root moves act as

    O:(D,E,U,1)->(D,0,2U,1),
    Z:(D,E,U,1)->(D+U-E,1,2U-1,1).

Every integer `U>=2` is attained: invert an even `U` to `U/2`, and an odd `U>=3` to `(U+1)/2`. The inverse strictly decreases until two. The resulting word has exactly `ceil(log2 U)-1` letters. This is the inherited Session199 binary response compiler.

For `N>=2`, build its word at `U=N` and then append `O`; for `N=1`, use the original edge. The common source then has `U=2N,E=0,d=1`, and `1<=D<=N`. The upper bound follows from `f>=0`, which gives `D<=U-1+E` before the extra `O`, and from `E` being zero or one.

Appending respectively `Z` and `OZ` gives

    L=D+2N,  R=D+4N,  R-L=2N,
    2N<L<=3N,  R<=5N.

Both graphs are one-apex tournaments whose spine starts with an `O` vertex and ends with a `Z` vertex. The root reaches the first spine vertex, the chain reaches the last, and the last reaches the root, proving strongness. In the upper graph, delete the penultimate newborn vertex, the extra `O` immediately before final `Z`. Every surviving arc is exactly that of the lower graph after increasing-order relabeling. This supplies the promised actual induced extension.

The upper order is `ceil(log2 N)+4`, including `N=1`, and the lower is one smaller. For `N=1`, the counts are `3,5`, consistent with THM-4876's all-order restriction on a deletion difference of two. The result does not give such a gap from an arbitrary supplied parent count.

## 6. Endpoint suffixes and infinitude, with credit and primality boundary

Direct exact multiplication of the inherited O/Z/J maps gives

    H(source.JOJOZ)=9D-3E+9U,
    H(source.JOOOZ)=9D-E+9U.

Every intermediate `J` or final `Z` closes a strong graft whenever the current root reaches the old vertices, and `O` preserves that reachability. Hence both outputs are strong for a strong initial source, with the same order `|source|+5`.

THM-5001's endpoint selector gives `E=N`: use the transitive n-vertex base with only its extreme arc reversed, and add a new root with incoming neighbors selected by the binary weights `(1,2,...,2^{n-3})`. Taking any larger `n` still supports the same binary digits, preserves endpoint count `N`, and remains strong. Its total count is at least the induced base count `2^{n-2}+1`, which tends to infinity. The lower suffix count is at least `9D-3N`, so the resulting pairs are unbounded.

This proves the alternative mechanism, but not new universality or infinitude: `session199/physical-interpolation/PROOF.md` Section 4 already gives infinitely many strong same-order single-arc pairs for every even difference. Its frozen pre-audit label is superseded by its README and passing `audit/INDEPENDENT.json`. Session221 and the earlier apex-wedge packet also contain all-even difference constructions.

The lower endpoint-suffix output is always divisible by three and is greater than three. Therefore this specific infinite family does **not** produce prime pairs. The synchronized four-value window likewise contains `6Q+3H`, always divisible by three. Difference coverage for Hamilton counts is not Polignac, primality, consecutive-prime recurrence, or Chan–Pak absolute coverage.

## 7. Gap replay and reproduction

`verify_gaps_independent.py` imports only this audit's independent graph routines, never a producer script. It checks:

- 10,000 compact target differences by exact word arithmetic, including the precise order and `5N` bounds;
- both literal full states and strongness for `N=1,...,128`, and the induced penultimate-vertex deletion arc for arc;
- both suffixes on all 285 rooted labelled sources through order four;
- the endpoint selector and both literal suffix outputs for `N=1,...,12`, at three increasing base orders each.

Normal and optimized outputs are required to agree. Reproduce from repository root:

    python3 -B 05-knowledge/results/session243_prime_gap_bridges_20261110/construction_audit/verify_independent.py
    python3 -B -O 05-knowledge/results/session243_prime_gap_bridges_20261110/construction_audit/verify_independent.py
    python3 -B 05-knowledge/results/session243_prime_gap_bridges_20261110/construction_audit/verify_gaps_independent.py
    python3 -B -O 05-knowledge/results/session243_prime_gap_bridges_20261110/construction_audit/verify_gaps_independent.py

Both scripts have `--fail`; it must return nonzero in normal and optimized mode. No floating point, random sampling, external solver, or primality oracle supplies any proof. The finite tests are implementation controls; the all-order conclusions follow from the arguments above.

## 8. Inheritance search and publication limits

The audit read the exact current Git blobs for the sparse inherited packets. `SOURCES.json` records the inspected source paths and hashes. The important prior mechanisms are:

- Session225 graft recovery, full scalar row and O/Z/J meanings.
- Session226 strong-graft iff criterion.
- Session227 multiroot proof, Section 3, scalar diagonalization and ordinal-cover multiplicativity.
- Session199 physical interpolation, Section 4, all-response binary inverse and infinite all-even differences.
- Session221 physical exclusion-tape construction, explicit all-even single-arc differences.
- THM-5001, complete strong endpoint spectrum.
- THM-4876, actual deletion-increment restrictions, compatible with the compact boundary examples.

No absence-of-prior-art theorem is asserted. At the pinned audit baseline `9332291ed45937bfd949b10111a59246cdc201f4`, THM-5109 was RESERVED and was not used as a proved dependency. Concurrent newer checkpoint `be5ce9a33c369e99327b29cd4b3417ddb087ca95` renames that theorem and labels statements C–F PROVED in a pre-audit checkpoint; that newer file was not re-audited here and is not a dependency of these results. The finite-support obstruction is not an unrestricted fixed-menu or root-changing obstruction. No output divisor is silently promoted to a Hamilton count, and no conditional primality statement is promoted to surjectivity.

## 9. Final-prose and arithmetic supplement

The final `constructions/PROOF.md`, `constructions/GAP_UNIVERSALITY.md`, and compact canon statements THM-5111 / THM-5112 pass. The producer incorporated the inherited-diagonalization and prior gap-infinitude credits, nonunit Euclid hostile, exact strong-order distinction, and explicit prime limitations before this acceptance. THM-5111 now explicitly qualifies the product identity by a nonempty gadget word, as requested by this audit, preserving the checked empty-gadget boundary.

`verify_arithmetic_independent.py` separately checks every displayed local sieve row, the full combined period 27,720 against the reduced 330-period test, its 40 surviving classes, all listed prime/gap/consecutiveness instances for m=1,...,30, both m=23 factorizations, the E6 graph pair `(17,224)`, and the physical Brown root counts 5,11,71 plus their actual ordinal-square graphs with counts 25,121,5041. Root counts also use direct full permutation enumeration. It independently reproduces all 5,044 nonsink rooted instances through order five, their 43 full states and 39 coefficient pairs. It checks the `5N` equality classification through N=10,000 and all 121 O/Z/J words of length at most four for the claimed absence of an endpoint-only difference row `(0,2,0,0)`.

The exact factorization-free compact equality argument is sound: before the final common O, the deficit `U-D` is zero exactly after a nonempty all-Z word. Such words have U=2^k+1, k>=1. Any O creates deficit at least two; the first following Z can reduce it by at most one, and further consecutive Zs leave it unchanged. This is sharpness only in the declared constructor, not a minimum over all possible graphs.

The supplementary reproduction commands are `python3 -B construction_audit/verify_arithmetic_independent.py` and the same with `-O`, from the packet directory. All three audit scripts use explicit exception gates; `FAILURE_GATES.json` records six expected nonzero exits. Normal/optimized receipt pairs are byte-identical. This audit introduced no mathematical correction to the final claims beyond the scope and inheritance clarifications recorded above.
