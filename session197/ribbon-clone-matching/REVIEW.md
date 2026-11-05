# Independent review: exact clone ranks and finite-seed obstruction

**PASS, no correction.** The all-order proof, nontrivial strong-seed hypotheses, exact cycle skeleton normalization, leading bases and counting obstruction are sound. No claim about the unrestricted Hamilton-path spectrum follows from the clone-only language.

## Frozen source binding

In `session197/ribbon-clone-matching/`:

- `PROOF.md`: `f1f147491e1a4b8f24fe8965e683235b607288147993773ca97f4f6657c190ff`.
- `graphs.py`: `889c0bfd92e7e788443929b65688b22b40c447a1843b57ef87d44af8d812ab8a`.
- `verify.py`: `162194bc4d5172fe05bc1c0736265e215c7a4b027d80680ba67a8f1835ec93f7`.
- Both validation receipts: `a009f57868315d19b32f2caf788b299848c1700c2715a9c25c88fcaff6ff2604`.

## All-order audit

1. The induced-path insertion rule is always available and preserves the complete old word. For cycles, both endpoint tournaments are strong. A mixed outside vertex inserts in an old cyclic edge. If none exists, the two uniform classes must both be nonempty and strongness forces an arc from the class dominated by the current tournament to the dominating class. The resulting two-vertex ear inserts into any cyclic edge. Choose new vertices deterministically from the current vertex set; deleting all added vertices recovers the original cyclic order. Each enlarged induced set is strong, so the process continues. No monotonicity through a weak intermediate tournament is assumed.

2. Fixing an old vertex as the cycle anchor removes the rotation factor before skeleton counting. Each skeleton has an ordered list of its clone gaps. An ordered nonempty partition of the m distinct clone vertices into these gaps gives exactly one graph cycle because every gap must be internally increasing. Its number is j! S(m,j), with no extra division by j or m. The same argument holds for the linear path gaps, including the two exterior gaps. All skeleton coefficients are nonnegative integers.

3. Every internal old fragment from one clone run to the next starts in O and ends in I. Its first change from O to I supplies one actual crossing edge. Different fragments are disjoint, giving a matching. Conversely, disjoint crossing edges provide a path with nu+1 clone runs and a cycle with nu runs on the selected rooted-triangle subgraph. Arbitrary old-vertex insertion preserves the former run separations. Strong cycle extension preserves the latter cyclic separations. The selected triangle union is strong, so the cycle extension hypothesis holds. This establishes exact, rather than merely upper, ranks. Since j! S(m,j) has leading base j with positive coefficient, no cancellation can remove the maximum exponential base.

4. The explicit cycle lower bound fixes one return pair and orders the others in (nu-1)! ways, then assigns the clones to the nu gaps. Distinct data have distinct recoverable cycles. At m<nu the bound is zero. The matching-one graph is a nonempty star under the strongness assumption. Its cycle coefficient is constant; each old cycle supplies a distinct two-run path skeleton, proving A>=C(T). With rank at least two, using just two disjoint crossing edges gives C(T_m)>=2^m-2, including the claimed six-cycle value at m=3.

5. Replacing any vertex within one transitive module by another transitive module preserves that original module's total order and its uniform incidences to all other original modules. Induction therefore keeps every descendant a positive-multiplicity blow-up of the original seed, even when the chosen vertex changes. Every strong seed vertex lies in a rooted triangle. Retaining its whole clone module and singleton representatives of the other two triangle vertices gives H=2^m+1. Induced-path monotonicity bounds every multiplicity by floor(log2(X-1)). Counting the positive multiplicity vectors proves the stated finite sum bound and zero odd density. Weak singleton padding is explicitly excluded from this coordinate bound; a singleton seed contributes only 1.

## Independent controls

The independent `check.py` does not extract skeletons from interpolated subset-DP counts. It enumerates old-vertex permutations and multiplies the local allowed gap factors directly. Linear gaps include both endpoints; cyclic gaps use a fixed old anchor. It independently obtains crossing rank by augmenting paths and counts actual clone Hamilton paths/cycles by literal depth-first word enumeration without memoization.

- All 570 labelled strong tournaments through order five, with all 2,822 roots, pass exact path/cycle skeleton degree, nonnegativity, cycle lower bounds and binary amplitude checks.
- 9,306 literal clone counts agree with the directly enumerated skeletons.
- The regular seven-vertex rank-three control has path coefficients (175,227,98,6) and cycle coefficients (17,16,2); five literal clone sizes agree, including both exact leading degrees.
- Independent normal and optimized output is identical.

Artifacts in this directory:

- `check.py`: `7410dc8bfadf4e7c90208d6da70edd7e9874f651c8570e076eb856bcfeb990bc`.
- `independent.json` and `independent-optimized.json`: `22841a4084b7254a246dd49b5852cc9272e013044786e35d24b96a0177d7e3e7`.

The source verifier and graph helper were copied to `source/` before replay because the verifier writes its receipt. Both source modes reproduce the exact frozen validation bytes. The original source hashes remain unchanged.
