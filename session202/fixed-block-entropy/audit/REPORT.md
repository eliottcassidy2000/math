# Independent audit: fixed strong blocks with arbitrary root incidences

**PASS, without mathematical correction to the frozen proof.** The all-order support bound and its multiplicative-closure consequence are valid for the stated fixed-block family. They do not obstruct arbitrary internal singleton spacers, unrestricted block banks, or general tournament constructions. No eventual odd-count threshold is proved.

## Frozen inputs

The audited `PROOF.md` hash is
`96aac520797bfc6f8aef07ca9959de4cc8d51c87378b711cf92530ca5f542620`.
The main verifier hash is
`223e8ee44fb4610c93fed0c9462349cfc35fabdb3dfc424b03ede7f4f7b8e67b`.
The seed verifier hash is
`9f54e3f7164547fcf37fe291f26292a88b0fc2f54f2fea3eda5ca5b5eb028e30`.
These files remained unchanged during review. `INPUTS.json` records the frozen files and the inherited root-amalgamation formula dependency.

## All-order proof review

1. **Two-bin constant and endpoint paths.** The identity defining beta counts exactly ordered, possibly empty path pairs. Every induced tournament contributes at least one path, and the two extreme subset terms each equal h. The proof that every vertex of a strong tournament can end a Hamilton path is valid: an outside vertex with an arc into a longest path ending there can be inserted before the first vertex it beats, without changing the endpoint. Otherwise the path dominates its nonempty complement, contradicting strongness. Thus beta>=2^s+2(h-1)>2^s for the fixed strong block of order s>=3.

2. **Physical injection.** A strong output has the required root incidence at each extreme block. The displayed Hamilton word uses every vertex exactly once. All nonroot joins have increasing block index, and the two root joins use the prescribed arcs. The two bins of each interior block are recovered from its vertices before and after the unique root position in the full path, so distinct ordered covers cannot collide. The lower bound beta^(m-2) is therefore an actual injection, valid for every root-incidence pattern that gives a strong output. For m=2 the empty product correctly gives one path.

3. **Weak SCC factorization.** Each strong block lies wholly in one SCC. The root SCC contains a consecutive interval of blocks, since an intermediate block receives a path from the first and has a path to the last. An SCC not containing the root cannot contain two blocks: every path crossing backward between their indices would pass through the root. Hence all other SCCs are exactly whole copies of the fixed block. Hamilton paths traverse tournament SCCs in their forced order, so H=h^(m-ell)H(C) is exact. The ell=0 case has singleton root component and H(C)=1.

4. **Support estimate and constants.** For root components with ell>=2, the injection bounds ell by2+log_beta X. The ell=1 case also obeys that bound for X>=1. Summing the at-most2^(s ell) incidence presentations gives the stated constant2^(3s)/(2^s-1). This deliberately overcounts scalar values. Summing over the weak padding exponent converges because h>=3 and rho>0. The pure powers contribute only O(log X), which is absorbed by O(X^rho). No order cutoff is incorrectly imposed on neutral external source/sink padding.

5. **Arbitrary ordered-join closure.** The generator support bound implies convergence of sum a^(-t) for every t>rho by dyadic grouping. Because a>=2, that also implies convergence of the displayed Euler product. Its expansion counts exponent vectors and overcounts distinct numerical products, so collisions never invalidate the inequality. It does not count ordered multiplication words as independent values. Choosing rho<t<1 proves density zero. Taking t=1 proves convergence of the reciprocal sum; the conclusion that infinitely many odd primes are absent follows as stated. Alternatively it follows directly from density zero and multiplicative closure, because inclusion of all but finitely many odd primes would imply inclusion of a positive-density set of odd integers.

6. **Singleton scope.** The excluded singleton has h=1,beta=2 and exactly two incidence choices, so the strict entropy gap is absent. Neither the fixed-block presentation bound nor its weak-padding summation can simply be applied after arbitrary internal singleton decisions are added. External universal sources/sinks preserve H and are harmless. The proof does not identify these two very different padding operations.

7. **Inherited arithmetic benchmark.** Section7's coefficients, index order, scalar recurrence and graph-order formula agree exactly with `ROOTED_ORDERED_AMALGAMATION.md`, Sections4-5, whose hash is0166889f66a6363506d336dcf7ae3d34a768f4ae79be5d387cfb426d78fe2883. Prefix q is read most-significant first with O=1,Z=0, while suffix t is read least-significant first with Z=1,O=0. The prime-atom sufficiency proof handles all cases of7^a3^b r correctly:49 and343 handle every a>=2, a represented7p handles a=1,r>1, and63 handles a=1,r=1,b>=2. Counts7 and21 remain the explicitly excluded cases. This is a sufficient OPEN occupancy obligation, not a consequence proved by the entropy theorem.

## Independent exact checks

The independent `check.py` imports no author code or data. It uses a separately written graph compiler, Tarjan SCCs instead of transitive closure, memoized suffix path counts, and permutation/cut enumeration of the two-bin constant. Its normal/optimized outputs agree.

It verifies all896 declared incidence presentations, including the exactly specified deterministic38-presentation sample, and independently matches every frozen author summary. In addition, it literally constructs and validates6,034 distinct injected Hamilton paths across the strong controls, rather than only comparing counts with the lower bound.

It confirms h,beta=(3,12) for C3 and(5,32) for the four-vertex block, every root-component interval, all nonroot SCC blocks, both SCC-product formulas, the singleton beta=2 boundary, and literal Hamilton counts49,63,343 for the three supplied graph words. The exponents are therefore rho=3log2/log12 and rho=4/5 respectively.

The author main verifier was replayed normally and under optimization; the receipts agree with each other and the frozen normal receipt. The seed verifier was also replayed and matches its frozen receipt. Finite checks do not supply the all-order result; the proofs above do.

Reproduction from this audit directory:

    python3 check.py > check-normal.json
    python3 -O check.py > check-optimized.json
    cmp check-normal.json check-optimized.json
    cmp check-normal.json NORMAL.json

The separately proposed CONTENT companion is outside this frozen audit.
