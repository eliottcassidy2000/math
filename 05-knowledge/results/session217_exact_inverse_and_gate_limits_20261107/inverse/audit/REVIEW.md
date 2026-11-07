# Independent audit: exact three-core membership algorithm

Verdict: PASS for the mathematical algorithm and stated arithmetic-operation
complexity. All-target occupancy is not proved. Review of the 80-state C++
implementation found no correctness defect; its numeric domain is 1<=N<=10^15,
not arbitrary-precision inputs. Runtime or memory feasibility at every input
in that domain is not asserted by this audit.

Reviewed source: session217-audit/SUFFIX_INVERSE.md, check_suffix.py,
inverse.cpp, and verify_combined.py as supplied on this session. Source hashes
at audit completion are in INPUT_SHA256.txt. Later explanatory edits need
comparison against that freeze, not an assumption that a new file was read.

## Proof checks

1. Reversal preserves f because its defining pair weight depends only on the
pair distance. Inverting R via z=rev_b(R) therefore gives precisely
X=A0 z+h f(z)-2h rev_b(z)+E(z mod2).
2. Stripping low z bits gives the stated K,S,C recursion. After t bits,
C is obtained by setting the remaining high word to zero, so the claimed
closed form in the low word l follows. In particular C>=-2h 2^b.
3. The necessary suffix interval uses x<=(K+h)u+h and
x>=Ku-2h(2^m-1); both directions are correct even when x is negative.
Its real width is
h*x/[K(K+h)]+2h(2^m-1)/K+h/(K+h). The bounds A0>=h,
K>=A0 2^t, x<=4A0 2^b+h give the displayed three-term bound.
The depth t>=2 and 4^t>=24h 2^b/A0 make it at most 9/16<1;
when t=b the residual word is simply zero.
4. Prefix recursion retains K-h>=2^t q2^b and C_t>=C_0>0.
Its positive-prefix interval is thus less than one unit at the chosen depth.
The residual u=0 branch must be tested separately. Both implementations do so.
Short original prefixes are valid: their high residual is zero, not a missing
branch requiring a padding-length promise.
5. For each large b, summing 1+sqrt(h2^b/A0) over legal P gives O(sqrt N).
The legal P=0 term is bounded separately with A>=h. For each small b,
2^b times 1+sqrt(hN)/(q2^b) is O(sqrt N). There are O(log N)
possible suffix lengths. Exact function evaluation costs O(log N) arithmetic
steps, establishing O(sqrt N (log N)^2) arithmetic operations. This is not
a bit-operation bound with unit-cost unbounded multiplication.
6. Legality is preserved: the empty-core-cut P=0 case is excluded; P>0 uses
its minimal binary length. A0>=h implies all b are bounded by log2(N/h).
The C++ loop b<=log2(N/5) is a harmless superset for this bank.
7. The proof extends to any fixed finite bank of nonempty strong cores and
specified cuts. For a nonempty cut, splitting every Hamilton path immediately
before one fixed selected vertex injects h paths into A covers, so A>=h.
Also q>=2h, G<=q-h and S,E<=h. These are the only bank inequalities needed;
constants may depend on the fixed bank. The C++ input contract remains the
supplied 80 states, not arbitrary banks.

## Independent implementation and controls

`audit.py` was independently written and imports none of the producer code.
It returns every candidate recipe rather than stopping at the first one.
A separate direct enumeration uses the all-length alpha<=N-h bounds. The
compared objects are entire (core,P,R,b) recipe sets, not just hit labels.

Every target 1..1500 agrees, covering 8,366 recipes. A further 106 targets
through 10,000 agree against an expanded all-length brute universe. These
include even targets, 47 and 63, 265, 287, 1237, and dyadic-neighbor controls.
The receipt records the exact scope. The original inherited graph formula is
used by both evaluators; the producer's separate literal graph DP checks
cover that different question. No finite test is used to justify the
asymptotic bound or to infer universal target membership.

Reproduce from any working directory:

    python3 session217-occupancy/inverse-audit/audit.py

The independent checker uses Python integers, so its arithmetic is not subject
to the C++ input cap. It deliberately runs only its declared bounded controls.

## C++ integer audit

All depth thresholds use __int128, and all loops and ceilings use integer
arithmetic. No floating logarithm can omit a branch. Correct signed floor
and ceiling division are supplied. The principal output product is promoted
to __int128. In the long-suffix branch, M<=N/5 and A0*M<=N, hence the exact
node constant C=A0*l+h*f(l)+E(l mod2)-2h rev_b(l) is O(N), and the node
coefficient obeys K<=(A0+h)2^t<=N+hM. Candidate prefix sizes are O(N),
short-suffix coefficients/constants are O(sqrt N), and all shift counts remain
well below 63 at N<=10^15. These bounds leave ample signed-64 headroom for
the smaller unpromoted products. The release proof should retain explicit
bounds and the numerical cap rather than describe this executable as
arbitrary precision.

## Final depth-first revision and audit lineage

The initial breadth-first source freeze is preserved in INPUT_SHA256_BFS.txt
(code hash 88e9b41fe217aba8f55a46872412a3a1b30e5714e512f9d643e410abb28eb467).
The final code hash is
ad7349e751ac2acdbadd01368ce46faa57144253ba1d43b60f05e4e409b16b1c.
Both recursion trees now use depth-first traversal. Direct line review confirms
that their transitions, depth thresholds, low-word coordinates, residual-zero
tests, interval pruning, and final witness tests are unchanged. Each call has
only one active child, so O(log N) live node records suffice. Arithmetic-node
sizes are O(log N) bits, giving O((log N)^2) extra bits in the abstract
arbitrary-precision implementation. This memory claim concerns the algorithm,
not the deliberately breadth-first independent Python checker.

The final code was independently compiled with undefined-behavior sanitizer:

    g++ -std=c++17 -O1 -g -fsanitize=undefined -fno-sanitize-recover=all session217-audit/inverse.cpp -o /tmp/session217_inverse_audit_ubsan
    python3 session217-occupancy/inverse-audit/check_final_cpp.py

Its support agrees with the independent all-length brute universe for every
N=1..3000. Seven extra targets include 10^12+39 and exact integer-depth
threshold cases 392,1568,6272. All returned witnesses pass the independent
arbitrary-precision formula check, with zero sanitizer findings; see
DFS_RECEIPT.json. This checks final implementation behavior separately from
the earlier all-recipe-set proof checker.

One minor exposition correction was requested in the new signed64 paragraph:
minimal threshold selection gives 2^t<=2sqrt(2hN)/(qM), not strict <,
because the C++ loop continues when equality holds. At h=9,q=84,b=0 and
N=392*4^s the upper equality occurs. The following 17sqrt N bound and all
conclusions are unchanged. The final input freeze records the corrected prose.
