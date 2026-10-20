# Infinitely many mixed-prime symmetry obstructions

**Status:** PROVED AND INDEPENDENTLY AUDITED: exact deficit one for the infinite family, using the inherited relative-even Burnside theorem and one CITED classical analytic-number-theory input for infinitude. Structural checks are included, together with an independently completed quotient census at order 22. No full census is claimed at order 38.

**Main result.** For the odd affine-Paley product group below, `tau_E(H)=tau_T(H)-1`. The exact-gap proof is §9; §§1–6 first establish the strictly weaker inequality without needing graph fusion-freeness.

This is a separate continuation. The frozen order-ten `C15` counterexample and its artifacts are unchanged.

## 1. A general no-fusion and sign-escape criterion

Let `H<=S_n` be odd-order and `N=N_(S_n)(H)`. Write `Fix_T(H)` and `Fix_G(H)` for the labelled invariant tournaments and graphs. A graph is `N`-even if its edge-reversal sign is trivial on `N∩Aut(X)`.

Use the relative-even theorem, THM-4721 Theorem R, from

`01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md`:

`#(N\Fix_T(H)) = #(N\{N-even graphs in Fix_G(H)}).                (1)`

Call the tournament side **fusion-free** if every full tournament isomorphism type meeting `Fix_T(H)` is exactly one `N`-orbit. This is a statement about the orbit partition, not that every individual isomorphism normalizes `H`.

**Criterion.** If the tournament side is fusion-free, then

`tau_E(H) <= tau_T(H)`.

The inequality is strict if any `H`-fixed graph is `N`-even but globally odd.

*Proof.* Full-even graphs form an invariant subfamily of the `N`-even graphs. Each full-even type contributes at least one `N`-orbit. A globally odd `N`-even graph contributes an additional orbit that cannot belong to any globally even type. Apply (1), identifying its left side with `tau_T(H)` by the no-fusion hypothesis. QED.

There is an exact nonnegative defect formula. Let `m_X` count the `N`-orbits of `H`-fixed labelled representatives in an even full graph type `[X]`, and let `q_X` count the `N`-orbits of `N`-even `H`-fixed representatives in an odd full graph type. Then, in the fusion-free tournament case,

`tau_T(H)-tau_E(H) = sum_([X] even, m_X>0)(m_X-1)
                     + sum_([X] odd) q_X`.                    (2)

Thus graph-type fusion and sign escape are separate nonnegative defects. This specializes THM-4721's already-proved class identity; it is not a new Burnside principle.

## 2. The affine-Paley groups

Let `a<b<2a` be primes with `a=b=3 mod 4`. Put

`k=b-a`, `n=a+b+k=2b`,

and split the vertices into sets `A,B,R` of sizes `a,b,k`. For a prime `p=3 mod4`, define

`H_p={x -> ux+v on F_p : u is a nonzero square, v in F_p}`.

It has order `p(p-1)/2`, which is odd. Let

`H=H_a × H_b`,

acting on `A` and `B` separately and fixing `R` pointwise. The two distinct prime divisors `a,b` make this a mixed-prime odd-order group.

**Normalizer lemma.**

`N_(S_n)(H)=AGL(1,a) × AGL(1,b) × S_k`.                        (3)

*Proof.* A normalizer permutes `H`-orbits preserving their cardinalities. Thus it preserves `A`, `B`, and the fixed-point set `R`. Its restrictions normalize `H_a` and `H_b`. In `H_p`, the translation subgroup of order `p` is the unique Sylow `p`-subgroup: the number of Sylow subgroups divides `(p-1)/2<p` and is one modulo `p`. It is characteristic. Its normalizer in `S_p` is the affine group, since normalizing the translation `x -> x+1` implies `f(x+1)=f(x)+c`, hence `f(x)=cx+d`. Conversely the affine group normalizes the square-multiplier subgroup. Independent root permutations are unrestricted. QED.

## 3. Internal invariant graphs and tournaments

The group `H_p` has one orbit on unordered pairs. Indeed, translations reduce a pair to `{0,d}`, square multipliers act transitively on square differences and on nonsquare differences, and forgetting order identifies `d` with `-d`. Since `-1` is a nonsquare, these two classes merge.

Therefore its invariant graphs are exactly the empty graph and `K_p`. Its invariant tournaments are the Paley tournament

`x -> y iff y-x is a nonzero square`

and its converse. The normalizer element `x -> -x` interchanges these two orientations.

The same negation is an odd automorphism of `K_p`: it has `(p-1)/2` transpositions on vertices, an odd number. On a complete graph the edge-reversal character is the ordinary vertex permutation sign. As all external adjacency of an invariant block is uniform, this odd internal automorphism extends by the identity outside that block.

Consequently every `N`-even `H`-fixed graph has both `A` and `B` independent. Its remaining data is an arbitrary graph on a weighted quotient with two distinguished vertices of weights `a,b` and `k` ordinary roots.

## 4. Prime-order vertex-transitive tournament blocks are prime

The following elementary proof avoids importing a modular-decomposition theorem.

**Lemma.** Every vertex-transitive tournament of prime order is prime: it has no proper nontrivial module.

*Proof.* Vertex transitivity makes the tournament regular. If `M` is a module, all its vertices have exactly the same external outneighbors, so its induced tournament is regular as well. A nontrivial regular tournament is strongly connected. For example, if a source strongly connected component has size `s<m`, its vertices have average outdegree

`(s-1)/2+(m-s) > (m-1)/2`,

contradicting regularity on `m` vertices.

Every strongly connected module is overlap-free. Suppose modules `M,W` overlap with `I=M∩W`, `U=M\W`, and `W\M` all nonempty. Choose `w in W\M`. For every `x in U` and `u in I`,

`x -> u iff x -> w iff u -> w`.

Because `w` lies outside the module `M`, the last orientation is independent of `u`. This gives a uniform cut between the nonempty sets `U,I`, contradicting strong connectivity of `T[M]`.

Now use a transitive subgroup acting on the prime number `p` of vertices. The translates of a proper nontrivial module have the same size; by overlap-freeness they are equal or disjoint. They cover the vertex set by transitivity, so the module size divides `p`, impossible for a proper nontrivial module. QED.

In particular each Paley block on `a` or `b` vertices is prime and strongly connected. Translation invariance supplies vertex transitivity.

## 5. Tournament fusion is absent

Every `H`-fixed tournament is a substitution of the two Paley blocks, each in one of its two converse orientations, and `k` singletons into a `(k+2)`-vertex quotient tournament. Use negation within either block to fix both internal orientations. This normalization lies in `N`.

Both large blocks are prime strong modules. Any other prime strong module must be disjoint from them unless equal to one: strong modules do not overlap, and a prime module cannot properly contain or be properly contained in another nontrivial prime module. Since only `k<a<b` singleton vertices remain, no other prime module of order `a` or `b` is possible. Thus the two large blocks are intrinsically unique by their orders.

Every full isomorphism between two such normalized tournaments preserves `A` and `B` and induces a permutation of `R`. The original isomorphism's restrictions to the Paley blocks need not lie in the affine normalizer. This is harmless: replace those restrictions by identity maps between the two fixed internal representatives. Cross-block and root-block arcs are uniform, so the replacement is still an isomorphism. Together with the induced root permutation it lies in `N`.

Undoing the initial internal normalizations proves that every full type meeting `Fix_T(H)` is a single `N`-orbit. This establishes the precise no-fusion hypothesis of §1 without assuming that all Paley automorphisms are affine.

## 6. A sign-escape graph and strict inequality

Take the complete bipartite graph with parts

`A∪R` and `B`.

Both parts have size `b`, so this is `K_(b,b)`. It is `H`-fixed. Every element of `N` preserves the two parts, since `N` preserves `A`, `B`, and `R` separately. Orient all edges from `A∪R` to `B`; these normalizer elements reverse no edge. Thus the graph is `N`-even.

A full graph automorphism can interchange the two parts. In that reference orientation it reverses all `b²` edges. Since `b` is odd, its sign is minus one. Hence `K_(b,b)` is globally odd.

The criterion now gives:

**Infinite-family theorem, strict form.** For every pair of primes

`a<b<2a`, `a=b=3 mod4`,

and the group in §2,

`tau_E(H) <= tau_T(H)-1`.                                      (4)

This is a family of genuine full-isomorphism-type inequalities. No graph fusion-free hypothesis is needed for the strict theorem.

## 7. An exact tournament-side formula and infinitely many parameters

Only the `k` singleton quotient vertices may be permuted. Thus `tau_T(H)` counts tournaments on `k+2` vertices with two individually distinguished vertices. For an all-odd partition `lambda` of `k`, append two parts of size one to obtain `mu`. Write

`z_lambda=product_l l^{m_l}m_l!`,

`e(mu)=sum_i (mu_i-1)/2 + sum_(i<j) gcd(mu_i,mu_j)`.

The elementary tournament Burnside formula yields

`tau_T(H)=sum_(lambda partition of k, all parts odd) 2^{e(mu)}/z_lambda`.  (5)

Permutations with an even cycle fix no tournament; for an all-odd cycle type, each unordered-pair orbit is oriented independently, giving the exponent above.

For `(a,b)=(7,11)`, `k=4`, formula (5) and the exact-gap theorem (6) give

`tau_T(H)=1408`, `tau_E(H)=1407`, with `n=22` and `|H|=1155`.

For `(a,b)=(11,19)`, `k=8`, they give

`tau_T(H)=874138624`, `tau_E(H)=874138623`, with `n=38` and `|H|=9405`.

These tournament counts are exact formula evaluations; the graph equalities use the audited exact-gap theorem in §9. The order-22 case is additionally checked by the independent quotient census recorded below. The order-38 value is a theorem consequence, not a full census.

**CITED classical input for infinitude.** NIST DLMF, equation 27.11.9, states the reciprocal-prime asymptotic in any reduced arithmetic progression. With residue 3 and modulus 4 it implies that `sum_(p=3 mod4)1/p` diverges. Source checked 2026-10-20: [DLMF 27.11.E9](https://dlmf.nist.gov/27.11.E9), with hypotheses `gcd(h,k)=1`, `k>0` stated immediately below the equation.

List these primes increasingly as `p_i`. If only finitely many consecutive pairs satisfied `p_(i+1)<2p_i`, then eventually `p_(i+1)>=2p_i`; the tail reciprocals would be bounded by a convergent geometric series. This contradicts divergence. There are therefore infinitely many pairs with `a<b<2a`, proving infinitely many examples in (4) and (6). This argument needs only the cited divergence, not a uniform prime-gap estimate.

## 8. Bounded structural checks and stopping condition

The standard-library script `infinite_family_structural_checks.py`:

- verifies all proper nontrivial vertex subsets are nonmodules for the Paley tournaments of orders 3,7,11,19;
- verifies regularity, pair transitivity, converse normalization, and the odd complete-graph negation witness;
- preserves a hostile: the regular cyclic tournament on nine vertices with differences `{1,3,4,7}` has three proper modules `{0,3,6}`, `{1,4,7}`, `{2,5,8}`, so the prime-order assumption matters;
- builds the pair orbits for `(a,b)=(7,11),(11,19)` from explicit generators;
- checks every normalizer generator preserves `K_(b,b)` with positive sign and checks the full part-swap witness has negative sign;
- evaluates (5) by exact rational arithmetic.

All 75 explicit control statements pass. The module scans examine respectively 3, 119, 2,035, and 524,267 candidate subsets. Normal and optimized Python outputs are byte-identical. This structural script performs no full census: its finite controls test the proof's building blocks, while the theorem handles every permitted prime pair.

**Independent exact-gap audit and order-22 control.** The separate reviewer accepted the normalizer replacement, role recovery, even false-twin corrections, and absence of graph fusion. Its proof record is `EXACT_GAP_FAMILY_ADDENDUM.md`; its separate standard-library checker is `exact_gap_family_check.py`, with matching normal/optimized results in `exact_gap_family_check.json` and `.out`. For `(a,b,k)=(7,11,4)` it exhausts all 32,768 normalized weighted graph presentations and obtains 1,992 full graph types, exactly 1,992 normalizer orbits, 1,408 normalizer-even types, 1,407 globally-even types, and 1,408 tournament types. The sole extra odd graph is `K_(11,11)`. These independent census counts are separate from the 75 structural controls above. There is no order-38 census claim.

Commands:

`python3 infinite_family_structural_checks.py > infinite_family_structural_checks.json`

`python3 -O infinite_family_structural_checks.py > infinite_family_structural_checks_optimized.json`

`cmp infinite_family_structural_checks.json infinite_family_structural_checks_optimized.json`

This closes the bounded extension: it establishes an infinite obstruction mechanism, the exact tournament-side count, and an independently audited graph-side deficit of exactly one. No unbounded search is required.

## 9. Exact-gap theorem: PROVED AND INDEPENDENTLY AUDITED

**Theorem.** For every prime pair and product group in §2,

`tau_E(H)=tau_T(H)-1`.                                        (6)

The proof below was independently accepted, with the additional order-22 exhaustive control in §8.

**Proof.** Every `N`-even graph has independent large blocks. A true-twin pair cannot involve a large-block vertex, since another vertex in that independent block distinguishes the pair. Any true-twin pair among roots gives an odd root transposition lying in `S_k<=N`. Hence all nontrivial twin classes of an `N`-even graph are false twins, and all internal class permutations are even.

If the two large blocks are in distinct twin classes, their sizes are

`a+r_A<=a+k=b`, `b+r_B>=b`, `r_A+r_B<=k`.

Root-only classes have size at most `k<a`. Therefore the two large roles are distinguished by class size, except when both classes have size `b`. Equality forces `r_A=k,r_B=0`; the two classes exhaust all vertices, and if they are distinct they must be joined completely, giving exactly `K_(b,b)`. If the large blocks share a twin class, that class has size at least `a+b` and is uniquely large.

Outside `K_(b,b)`, every full isomorphism preserves the twin class containing each original block. Compose with permutations inside the target false-twin classes to send the original `A` and `B` sets back onto themselves. These corrections are even. Then replace its internal restrictions to `A` and `B` by identities, which remain automorphisms because the two blocks are independent and uniform to their complements. The remaining root permutation lies in `S_k<=N`.

These operations prove both that every full-even graph type is one `N`-orbit and that every `N`-even graph other than `K_(b,b)` is globally even. For the latter assertion, an alleged odd full automorphism is corrected only by even false-twin permutations, so it would yield an odd element of `N`, contradicting `N`-evenness. The single exception `K_(b,b)` itself is one `N`-orbit: all `k` roots must join the original `a`-block, so its weighted quotient presentation is forced. Equation (2) therefore gives exact deficit one. QED.
