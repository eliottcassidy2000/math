# Session 57 integer biclique and Lagrange-family audit

Verdict: PASS for the stated mathematical claims and their precise bounds. No claim of a globally least K4,4 witness or a general MSOS impossibility is warranted or made.

## Verified results

1. The common-neighbor map is a bijection between positive same-parity factor pairs p<q of b-a and admitted x=((q-p)/2)^2-a. Positivity, the universe bound, and x not in {a,b} are essential and implemented. Strict decrease of (delta/p-p) makes the factor list injective. Delta congruent to 2 modulo 4 is excluded.
2. Complete enumeration proves the first ordinary K2,4 occurs at N=274 and the first ordinary K3,3 at N=673. The checker independently compares factor neighbors with direct integer-square tests. Every biclique has a seed pair; all right tuples of its common neighbors are enumerated. Taking the smallest permitted left labels is valid for minimizing the maximum label. Simple-graph adjacency prevents side overlap.
3. The factor-route K4,4 exclusion through N=10000 passes normal and optimized Python runs. A separately replayed adjacency/square-sum closure enumeration confirms the same exclusion, with 195219 eligible pairs and 731523 right quadruples. The adjacency route also reproduces K3,3 minimum 673, 660224 eligible triples' seed pairs, and 1836298 distinct right triples.
4. The primary-family box scan passes normal and optimized runs: 97856 primitive pairs in [-200,200]^2, 97780 valid positive disjoint normalized pairs, best maximum 1112833 at (-11,5). A separately implemented gauge minimizer searches the first nonforbidden integer on each side of the convex objective's center and reproduces the same result. The provided radius-17 method is also mathematically sound because at most 16 integer gauges are forbidden. The positive witness and all 16 exact roots match. Only bounded-family optimality is proved; the ordinary K4,4 global minimum lies in [10001,1112833].
5. The three symbolic A balances exactly factor as -v F12, F13, -v F14. The verifier checks all 48 projective points for each residual form over F47, including infinity. Independent plain-integer evaluation on every nonzero vector in F47² gives 6624 nonzero evaluations. Homogeneity and reduction of primitive integer coordinates exclude every rational projective zero, including rational parameters with denominators divisible by 47.
6. Therefore A is balanced exactly when v=0. At v=0, A=(2916,36,2916,36)u^4 and B=(0,220,2268,2268)u^4. The B balances are (-4316,-220,-220)u^4, so no nonzero rational parameter pair is doubly balanced. The origin is the collapsed degenerate exception. Permutations, part exchange, gauge translation, and nonzero square scaling cannot repair failure of double balance. This excludes this family from a decorated locus requiring both balances, not other families or all MSOS configurations.

## Source and attribution

The formulas and numerical example were checked against the supplied source at https://arxiv.org/html/2508.07806v1, introduction and equation (1.1). The source credits the parametric family to Jean Lagrange through Guy; the numerical example is presented through Guy without needing a stronger separate inventor attribution. The modular obstruction and computational optimizations are independently derived, not attributed to Choudhry. No full article is included in this audit. The source receipt is sufficient provenance for the final obstruction package.

## Reproduction and failure behavior

The copied scripts are in integer/ and family/. Run checker.py, verify_k44.py, primary_construction.py, and family/verify.py both normally and with python -O. Run integer/closure_search.py 10000 for the independent adjacency route. Run independent.py normally or with -O for the alternate gauge scan and full-vector modular check. Outputs are recorded beside this audit, with exact agreement between normal and optimized results apart from timing fields.

All theorem guard checks use explicit exceptions, not assertions. Mutation tests changed the expected K2,4 minimum and the expected boundary A tuple; both errors caused nonzero exits under normal and optimized Python. The family verifier does not catch its exception for pretty output, but correctly terminates nonzero. The closure enumerator is a generator rather than a claim validator and prints the actual search outcome; its null K4,4 result was explicitly inspected.

No original repository files or remotes were changed. Only audit-owned files were written. Input hashes are in input_receipts.json. Full article copies and Python cache files are not part of the deliverable.
