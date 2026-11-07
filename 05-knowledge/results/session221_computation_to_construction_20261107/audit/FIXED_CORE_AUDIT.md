# Session221 fixed-apex compression and inverse audit

**Verdict: PASS** for the final `session221-automata/PROOF.md` and `automata.py` reviewed here. The result is exact fixed-core membership with optional apex endpoint constraints, not all-target acceptance or polynomial time in the binary target length.

## Mathematical review

The proof uses the **discrete** Moon–Busch floor, not only its coarser exponential consequence. THM-4744, `THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md`, Setting, cites A. H. Busch, EJC 13 (2006), N3, Corollary 1 and states the three residue formulas. They are equivalent to `min{3^a 5^b:2a+3b=s-1}` for nontrivial strong order `s`. This inherited cited theorem was checked as a dependency, not independently reproved from the external paper.

Delete precisely singleton chain SCCs, retaining every labeled apex. Every remaining chain vertex belongs to a nontrivial SCC containing at least one apex. Thus its chain population is at most `s-1`. Multiply the Busch lower bounds across SCCs and add their exponents to obtain the exact integer horizon

    L(m)=max{2a+3b:3^a 5^b<=m}.

The mask test is correct: prior masks' union must be contained in the current mask, which must be contained in the later masks' intersection; the current mask must also be an outward apex cut. These are exactly the three possible triangle obstructions. The test is linear in tape length at fixed apex number.

For endpoint preservation, retain a singleton chain global source/sink when present. These are exactly the first mask 0 and last full mask. The SCC concatenation bijection then preserves each apex's start/end status jointly as well as separately, with at most two extra vertices. No higher path-cover or continuation equivalence follows. In particular the source's `[3]` versus empty prefix hostile correctly prevents prefix normalization in the inverse search.

The inherited positive gap transfer is exact by assigning each chain vertex to its gap between an apex permutation. The dimension is `k! * 2 * 3^k`; multiplicities are retained. The start, end, and joint functionals have the stated empty-boundary-gap interpretation.

For fixed `k>=1`, enumerating words to `L(m)` or `L(m)+2` has `O_k(m^(3k log_5 2))` nodes. Carrying the transfer state gives fixed-k many integer additions per edge, on `O_k(log(m+1))` bits. DFS uses `O_k(log^2(m+1))` working bits. Thus the advertised numerical-target time bound is justified without a general exponential subset-DP evaluator. Output storage for a full census is additional.

Requested clarifications were applied: empty graph is distinct from empty word (`T_A(empty)=A`); arbitrary endpoint input has its input-reading cost included; count bounds and array dimensions are checked; JSON lists and tuple constraints are normalized consistently.

## Independent executed verification

`audit_fixed_core.py` uses the candidate implementation only as the object being tested. The comparison path is an independently written backwards endpoint recursion, returning ordinary H, apex starts, apex ends, and the full joint endpoint matrix.

Completed:

- 778 literal instances: all labeled apex cores for k=0,1,2 through tape length 4, plus 60 seeded k=3/4 instances;
- transfer outputs and all endpoint functionals on every instance;
- both scalar and endpoint-preserving normal forms on every instance, with exact output preservation and the L/L+2 bounds;
- 54 distinct physical endpoint profiles, independently reconstructed from short words, passed as JSON lists to the constrained compiler; every returned graph independently recounted;
- two impossible endpoint-range rejections;
- complete scalar target exclusions 7 and 21 as inherited hostile controls.

No external claim or larger coverage range is inferred from these checks. The global compression idea and inherited saturation/transfer are credited; the assessed contribution is the fixed-core physical deletion test, endpoint refinement, and resulting smaller exact search bound.

## Release and replay

Release `FIXED_CORE_AUDIT.md`, `audit_fixed_core.py`, `fixed_core_receipt.json`, and `FIXED_CORE_SHA256SUMS` alongside the candidate `automata` directory. Run:

    python3 -B audit_fixed_core.py /path/to/automata

The script defaults to a sibling `automata` directory and supports the original session staging path as a fallback. No unresolved correctness blocker remains.
