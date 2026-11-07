# Complete square-root inverse for fixed-core templates

**Status: PROVED exact all-length membership algorithm; independently audited PASS.** The fixed three-core coverage conjecture remains OPEN.

Fix a legal core/prefix, with A0=qP+A>0 and tau0=P G+h f(P)+T+S(P mod2). The physical cover inclusion gives 0<=tau0<=A0, and A0>=h for this core bank. For a suffix of length b write z=rev_b(R). Reversal preserves the pair-distance function f, hence the exact count equation is

    X=D-N=A0 z+h f(z)-2h rev_b(z)+E(z mod2),
    D=2^(b+1) A0+h-tau0.

Maintain the residual equation

    X=K u+h f(u)+S(u mod2)-B rev_m(u)+C,

where B=2h, initially (K,S,m,C)=(A0,E,b,0). Strip one low bit:

    bit0: (K,S,m,C) -> (2K,0,m-1,C);
    bit1: (K,S,m,C) -> (2K+h,h,m-1,C+K+S-B 2^(m-1)).

After t bits, K>=A0 2^t and 0<=S<=h. If l is the stripped low-bit integer, C=A0 l+h f(l)+E(l mod2)-B rev_b(l), so C>=-B 2^b. Therefore x=X-C<=4 A0 2^b+h, using A0>=h and tau0>=0.

For 0<=u<2^m, use 0<=f(u)<=u to obtain the necessary interval

    max(0, ceil((x-h)/(K+h))) <= u
    <= min(2^m-1, floor((x+B(2^m-1))/K)).

The unclipped real interval has length at most

    6h 2^b/(A0 4^t)+1/4^t+1/2^t.

Choose t=min(b,max(2,ceil(log_4(24 h 2^b/A0)))). If t<b, this width is strictly below1; every branch leaves at most one integer to test in the exact equation. If t=b, simply test its sole u=0. All 2^t branches suffice. Cost per fixed P,b is O(1+sqrt(h 2^b/A0)) exact arithmetic evaluations, each O(log N) by direct f recurrence.

For suffix lengths 2^b>=sqrt(N), P<=N/(q 2^b), and summing over all P gives O(sqrt N) branches per b: the additive1 costs at most O(sqrt N); the square-root term sums using sum_{P<=M}(qP+A)^(-1/2)=O(sqrt(M)) plus the legal P=0 term. Since 2^b<=N/A0, the P=0 term is also O(sqrt N). The fixed 80-state bank contributes a constant.

For 2^b<sqrt N, use the prior prefix parity recursion with fixed R, whose coefficient K0-h>=q 2^b. After t low prefix bits, Kt-h>=2^t q 2^b. Its constant C_t is at least its original C>0, because K-S>0 throughout. The necessary prefix interval width is at most hN/(q^2 4^(b+t)); t=max(0,ceil(log_4(2hN/(q^2 4^b)))) leaves at most one candidate per branch. The 2^b suffixes times 2^t branches cost O(2^b+sqrt N)=O(sqrt N).

Thus all legal suffix lengths, bounded by log2(N/h), admit a complete exact inverse in O(sqrt N (log N)^2) elementary arithmetic steps. This is an algorithm for template membership and a way to construct recipes when present. It does not prove that every odd N>=65 is present.

## General fixed-core scope

The mathematical algorithm applies to any fixed finite bank of nonempty strong tournament cores, including every cut on each core. Its constants may depend on that bank. The executable and retained80-state checks implement only the three specific Session216 cores.

Indeed q>=2h because the ordered cover sum includes the two placements with an empty path. Since G excludes covers with first path empty, G<=q-h. Also S,E<=h. For any nonempty cut, fix a selected vertex and split each Hamilton path immediately before that vertex: the first path may be empty and the second starts in the cut. Concatenation recovers the original path, proving A>=h. For P>0, A0=qP+A>=q>=2h; for P=0 the legal nonempty cut yields A0>=h. The cover subset relation tau0<=A0 survives padding. These inequalities supply every hypothesis used above.

A strong core and the stated cut condition ensure the closed graph is strong. For arbitrary nonstrong cores the arithmetic remains valid, but the source-component reachability condition must replace the simple nonempty-core-cut rule. This extension does not assert that any fixed bank represents every large odd integer.

## Implementation, integer domain, and reproduction

The complete mixed algorithm is inverse.cpp. It uses integer-only threshold selection, signed128-bit products in threshold and final count tests, and an explicit input domain 1<=N<=10^15. This is a conservative implementation domain, not a restriction on the abstract arbitrary-precision theorem. The claimed operation count is arithmetic-operation complexity; it is not a unit-cost bit-complexity claim. All arithmetic has O(log N)-bit operands for a fixed bank, so a corresponding bit bound can be obtained by multiplying by the cost of O(log N)-bit arithmetic.

Build and reproduce:

    g++ -std=c++17 -O3 inverse.cpp -o inverse
    ./inverse CORE_STATES.txt 1 7 21 23 47 63 65 265 287 1237 4295971 9948407 1000000007 1000000000039 > INVERSE_CASES.jsonl
    python3 check_suffix.py > SUFFIX_CHECK.json
    python3 verify_combined.py

The exact membership statuses are template-relative. verify_combined.py independently enumerates every parameter tuple allowed by the physical alpha bound for targets<=3000, checks the complete support against the mixed inverse, reconstructs literal graphs, and counts60 small outputs by vertex-subset DP. Large outputs are reconstructed and checked for strong connectivity, but their Hamilton counts are certified by the inherited formula rather than whole-graph DP. CASE_GRAPHS.json includes those literal outgoing masks.

For N=1000000000039 the algorithm finds core-cut index0, P=16231501406, b=0, R=0, a legal strong recipe because its prefix is nonempty. This improves that target's prior balanced-band witness, without any global coverage consequence.

### Conservative fixed-bank integer safety bounds

For the implemented bank h<=9 and q<=84, M=2^b<=N/5. On the long-suffix branch every enumerated prefix has A0 M<=N, and h<=A0. Thus K_t<=(A0+h)2^t<=2N, |C_t|<=4N+h by its explicit low-word formula, and all interval numerators have absolute value below10N+10. On the short-suffix branch M<sqrt N. Its initial K<=168sqrt N and C<=186sqrt N+9. If t>0, minimal threshold selection gives 2^t<=2sqrt(2hN)/(qM), hence every K_t<=17sqrt N; if t=0 the initial bound applies. The sum of constant increments is bounded by twice the largest K_t, so |C_t|<=522sqrt N+9. Reconstructed prefix candidates satisfy qP<=N/M+O(sqrt N), with an absolute fixed-bank bound below20N+100, before the exact alpha test. After that test the formula terms are individually bounded by a small constant times N; its leading product is additionally promoted to signed128-bit arithmetic. All remaining signed64 operations, including shifts, are bounded conservatively by1000N+1000. This is below2^63 for N<=10^15. Threshold products are evaluated in signed128-bit arithmetic. The executable's decimal-input domain check is part of its contract.

The final executable traverses both low-bit recursion trees depth first. It stores only O(log N) node records instead of materializing an entire level, and stops at the first verified witness. Thus extra working storage is O(log N) machine words in the capped implementation, or O((log N)^2) bits for an arbitrary-precision implementation with O(log N)-bit node coordinates. A worst-case complete negative search still visits the stated square-root-times-polylogarithmic number of arithmetic steps; early success may visit far fewer.
