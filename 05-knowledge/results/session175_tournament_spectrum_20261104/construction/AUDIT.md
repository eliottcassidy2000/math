# Independent audit: one-feedback synthesis

**Verdict: PASS after the corrected modulo-eight example.** Audited the stable README and checker on 2026-11-04. No global spectrum-completeness claim is justified or made. The all-order exclusions of 7 and 21 are inherited; these results give explicit constructors and exact restrictions on one construction class.

## Proof review

1. **Binary formula and streaming state.** Deleting the distinguished vertex leaves two increasing sequences, uniquely determined by a subset and its complement. The end cases, adjacent cut, and reversed endpoint pair exhaust all paths. Appending 0 adds the binary-prefix weight and previous final bit. Appending 1 balances the terminal contribution against a new 01 adjacency. The empty prefix state is correct. PASS.
2. **Complete finite inverse decision.** Initial 0 vertices are sources and final 1 vertices sinks, so removing them preserves the Hamilton count. Any surviving core starts 1 and ends 0 and is strong. Its extreme inversion contributes 2^(n−2), separately from the two endpoint paths. This bounds the core length independently of padding. Thus the finite gap lists are complete for the entire deletion-one-transitive class, not merely graphs below a selected order. PASS.
3. **Four-run products.** The registers after 10^a and then 1^b give the stated product exactly. Because the word begins 1 and ends 0, its witness is strong. PASS.
4. **Asymmetric cyclic kernel.** With one singleton feedback vertex, the six skeletons have run multiplicities (1,1) with multiplicity three, (2,1), (1,2), and (2,2). They yield 3ac+q_A(2)c+a q_B(2)+q_A(2)q_B(2)=(F_A(2)−a)(F_B(2)−c)+2ac. PASS.
5. **Explicit positive family.** The four-vertex seed has H=5 and F(2)=32. A transitive r-vertex module yields 27·2^r−17. Every nonempty cyclic three-module substitution is strong. For r=3k+2, division by 7 modulo 7 gives 6+3k, exactly the stated exceptional class k≡5. This establishes infinitely many odd 7-adic valuations, not infinitely many primes. PASS.
6. **Orphan correction.** Legal insertion positions on an old path number 1+#10. Every remaining new Hamilton path has a unique split, last-left endpoint, and first-right endpoint with backward old arc. The proposed endpoint-sensitive subset sum counts each exactly once. PASS.
7. **Scope.** The asymmetric cyclic wrapper is broader than the mirror family but still lies within uniform C3 substitution. The separate obstruction packet shows that this broader architecture still excludes 7 and 21, and indeed many other small values; this does not conflict with its H=91 family. The binary feedback constructor is a different class and its gaps are not unrestricted gaps. PASS.

An initial sentence used H=11 as an example escaping H≡3 (mod8), but 11≡3 and its equal modules are mirrored. Author corrected this before the stable audit: H=23 is the appropriate residue escape; H=91 is the odd-7-valuation escape. No remaining mathematical correction is requested.

## Independent computations

Run `python3 audit_independent.py > audit_independent.json` and the same with `python3 -O`, then compare. Both runs passed and outputs agree exactly.

- 1,023 binary words of lengths 0–9: an independently written memoized suffix recursion agrees with the author's streaming and literal formulas. The graph for each word is built independently.
- Every odd gap through 199 agrees with the complete inverse solver, whose exact core-length limit is nine. This validates the complete list independently of the author's bound-10000 census.
- 1,098 extensions of all labelled tournaments on one through four vertices: explicit permutation enumeration validates the insertion-plus-orphan partition.
- The family at r=1,…,8 gives 37,91,199,415,847,1711,3439,6895 under the independent Hamilton counter.
- One hundred k values validate the residue formula for the 7-adic boundary.

The audit imports only the author module for comparisons and elementary encoders in the orphan/family tests; its Hamilton recursion and permutation oracle are independently implemented. The symbolic proofs, not these finite ranges, establish the general assertions. The author's endpoint-subset formula is also checked against explicit orphans in its 1,098 extension cases.

No prohibited path was read, no software was installed, and no source file outside this packet was changed.
