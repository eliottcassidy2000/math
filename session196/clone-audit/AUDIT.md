# Independent audit: transitive cloning and complete star synthesis

**Status: PASS.** The unrestricted `{7,21}` Hamilton-path spectrum remains open. This audit establishes no literature-priority claim and makes no inference from a finite census to an all-order theorem.

## Frozen target

The reviewed files are copied under `frozen/`. `SOURCE_SHA256SUMS` pins the author's proof, verifier, receipt, and normal/optimized outputs. The proof SHA-256 is `31982897f0092fafa040c9d60452fb2b69c128d4b7173b5d642d3ff56a622499`; the verifier SHA-256 is `7c87d7ad228a88b966384915ef0263ebef9290b3d6d3a9d0bf776d6de39bc0e6`.

The inherited named mechanisms were checked against the Session 195 rooted-interface and zero-slack packets. The latter defines the same quantity `R+e-h-d`, so its nonnegativity applies to the synthesis's `sigma`. The new coefficient identity does not depend on that inequality.

## All-order proof review

1. **Literal expansion.** Maximal clone runs have unique compression, and an ordered nonempty set partition of the clone set has exactly one transitive ordering within each run. Thus the coefficient is `k! S(m,k)`, not `S(m,k)` or a composition count. The formula is for positive m. The added deletion warning is necessary: a triangle gives `2^m+1` at positive m, while deletion has count one.
2. **Matching rank.** The cross graph uses exactly `r -> O -> I -> r`. Every gap between root symbols contains an O-to-I transition; disjoint gaps produce a matching. Conversely, disjoint matched arcs give a legal partial path with separated clones. Ordinary vertex insertion preserves their relative order and hence their separation, irrespective of the clones' internal order. This proves existence and positivity, including a cross graph consisting only of disjoint edges.
3. **Binary iff and return.** Matching number one is equivalent to a nonempty star, allowing isolated vertices. The top Stirling coefficient is positive, so distinct exponential bases cannot cancel the largest base. Matching number zero is the separate constant case. Transitive expansion adds only isolated vertices to each clone's cross graph. Further cloning therefore enlarges the same transitive block and preserves the offset, including its sign when the offset is a unit.
4. **Two-piece synthesis.** The converse operation correctly interchanges the two possible locations of the star center. With center c in the root's out-neighborhood, the remaining cross edges force `A -> B`, yielding precisely the claimed construction. `U=R-d` counts covers with a nonempty path after c and is positive exactly when c beats some A vertex. The five path-run cases in the proof are disjoint and exhaustive, including empty A or B. Their sum is exactly `UL 2^m + dg - sigma t`. Coprimality is only a necessary condition for a unit offset; actual rooted pieces are still required.
5. **Fixed-template bound.** Every vertex of a nontrivial strong component lies on a triangle. Deterministic insertion injects induced Hamilton words into full Hamilton words. The induced cyclic blowup therefore forces `H >= 2^(m_v)+1` for every active multiplicity separately. Singleton strong-component blowups do not affect the count. The stated polylogarithmic bound follows by counting active multiplicity vectors, not by assuming distinct vectors give distinct counts.
6. **Ordered products and primes.** Positive cloning preserves each original nontrivial strong component's template, and ordered joins only concatenate strong-component lists. Prime total count requires one prime-count component. The finite set of component templates still supplies only polylogarithmically many counts up to X. The dyadic reciprocal-prime/Euler-product argument rules out eventual coverage of all primes without using the prime number theorem. This does not apply to operations that create new strong templates.
7. **Scope.** The opposite-sign commutation obstruction fixes the two clone classes and their signs. It is not a prohibition on adaptive root changes or other operations. None of the small-order amplifier tables is promoted to an all-order exclusion.

No coefficient, root-incidence, positivity, or quantifier correction was required. During review the author added the explicit positive-m/deletion boundary clarification; this audit targets that final version.

## Independent computation

`audit.cpp` does not import or call the author's verifier. It constructs each expanded graph afresh from an ordered owner list, counts Hamilton words by removing their first vertex, finds matchings by recursive subset search, and computes quotient coefficients by a marked-word state DP with a reusable, nonconsecutive root symbol. The coefficient DP also includes `a_0=H(T-r)` so deletion is verified separately rather than by extrapolating a positive-m exponential formula.

The complete unfiltered labeled universe comprises all tournaments of orders one through six and every root:

- 33,867 graphs and 202,013 rooted graphs;
- matching-degree identity and both base cases `a_0=H(T-r)`, `a_1=H(T)`;
- 808,052 actual expanded-graph checks at m=1,2,3,4.

For synthesis, rooted profiles are counted by literal permutations and root-position cuts. The independent assembler uses the common-root model directly. All 285 rooted labeled inputs of orders one through four are paired:

- 81,225 ordered pairs;
- 243,675 actual synthesis checks at m=1,2,3;
- 27,722 pairs for which the formal positive-m exponential evaluated at zero differs from the actual empty-Z graph. These are intentional boundary controls, not theorem failures.

Two additional disjoint-edge cross-graph controls at orders seven and nine have matching ranks three and four. Each is checked on actual blowups at m=1 through five, extending beyond the maximal ranks possible in the exhaustive order-six universe.

The author's frozen Python verifier also passes independent normal and optimized replays. Both outputs and both generated receipts are byte-identical to each other and to the frozen author files. Its explicit exceptions remain active under `python3 -O`.

## Reproduction

From this audit directory:

```sh
g++ -std=c++17 -O2 -Wall -Wextra -pedantic audit.cpp -o audit
./audit > independent.stdout
cd frozen
python3 verify.py > ../replay.normal.stdout
cp receipt.json ../receipt.normal.json
python3 -O verify.py > ../replay.optimized.stdout
cp receipt.json ../receipt.optimized.json
cmp ../replay.normal.stdout ../replay.optimized.stdout
cmp ../receipt.normal.json ../receipt.optimized.json
cmp normal.stdout ../replay.normal.stdout
```

All work was performed in this dedicated scratch directory. The shared dirty repository checkout was neither edited nor staged. The parent session owns publication.
