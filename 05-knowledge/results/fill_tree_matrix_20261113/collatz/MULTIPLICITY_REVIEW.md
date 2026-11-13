# Independent referee review of the neutral-fill proof

**ACCEPTED: PROVED elementary consequences of the cited repaired neutrality
theorem and inherited complete atom inventories; FINITE-EXACT controls have
their declared scope.** This review independently rederived the arguments in
[multiplicity/PROOF.md](../multiplicity/PROOF.md), read its checker, and read
[THM-1900, the insertion-response calculus](../../../../01-canon/theorems/THM-1900-the-insertion-response-calculus.md)
and
[Session24's actual-count multiplicity proof](../../actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md).
It does not claim another full tournament census.

1. **Neutrality.** The repaired THM-1900 argument counts all insertion
   positions into old Hamilton paths and leaves orphan paths as nonnegative
   additional contributions. For a split strong component, a Hamilton cycle
   has a P-to-complement transition; deleting another cycle edge leaves that
   transition in a spanning path. For a union of components that is not a
   terminal segment, the component-bit word also has such a transition.
   Thus every nonterminal attachment strictly increases H, while a terminal
   SCC cut preserves it. This proves the all-order implication required by
   the packet; it does not reuse the false deletion-of-a-path premise.

2. **Full neutral components.** A neutral birth inserts a singleton SCC at
   a condensation cut and cannot alter a nontrivial strong component.
   Conversely, an H-neutral deletion is the inverse of such a birth, so the
   deleted vertex was a singleton SCC. Therefore fixed ordered words W are
   exactly the components for these birth/deletion edges. SCC order is
   intrinsic, including when some component types repeat. The gap vector
   gives an injective and surjective ordinary-isomorphism coordinate system.
   There really are g_i+1 different named cuts over successor g+e_i.

3. **Perfect fill with retained multiplicity.** With a fixed named initial
   core, the number of coordinate orders ending at g is d!/product(g_i!).
   The product of successive cut choices in those orders is product(g_i!).
   Their product is d!, and the common probability of each literal history
   has denominator (r+1)...(r+d). Distinct histories cannot merge as named
   final graphs: deleting vertices in reverse birth-name order recovers the
   history. Core automorphisms do not identify these named objects. Hence
   the quotient cells have exactly the stated uniform probability. Choosing
   quotient directions uniformly loses the g_i+1 weights and yields the
   correctly stated binomial hostile at r=1,d=2.

4. **Primitive differences.** The zero-truncated binomial kernel has series
   z^m/(1-z)^(r+1), so r+1 backward differences isolate its coefficient at m
   at every order, including the onset. The rank r must remain in the state:
   the H9 example gives (1-z)^3 F_9(z)=2z^5-z^6, exactly as claimed. The only
   clarification requested in review is now applied: the auxiliary
   empty-tournament convention is N_(0,0)(1)=1 and all other order-zero
   slices vanish. This is the
   convention already used by the inherited series and the r=0 urn control;
   it does not alter any positive-order formula.

5. **Four-target completeness.** Busch gives H(S)^3>=5^(|S|-1); the integer
   inequality 5^6>21^3 therefore excludes every strong order>=7 for all four
   targets5,7,17,21. Exhausting orders1..6 is consequently sufficient.
   Prime5 and17 allow only one nontrivial SCC. For21 the only potential
   nontrivial factors are3,7,21, and the missing7/21 atoms make a product21
   impossible. The rigid one/two positive atoms give N_n(5)=max(n-3,0)
   and N_n(17)=2max(n-5,0). The checker uses actual Hamilton counts and
   ordinary classes; it does not divide labelled counts by n! without first
   certifying rigidity of the relevant atoms.

6. **Thirteen-to-twenty-seven obstruction.** Any induced strong H13 core
   lies inside one SCC of an H27 host. Components of count3 or9 are excluded
   by induced-count monotonicity. The five strong H27 atoms are complete by
   5^7>27^3, and the three converse representatives require precisely
   6+21+21 five-subset tests. Converse preserves Hamilton counts. Their
   failure to contain H13 therefore excludes every host order and every
   singleton padding. This is stronger than a single-source birth-menu
   failure, while leaving old-arc surgery and deletion outside its scope.
   Independently, the SCC series is

       (2z^6+3z^7)/(1-z)^2 + 4z^8/(1-z)^3 + z^9/(1-z)^4.

   Its coefficients simplify to n(n-5)(n-4)/6 for n>=6, agreeing with the
   formula in the proof. There are ten core words, not ten strong atoms.

7. **Global coverage boundary.** Source padding makes the supports S_n
   nested. Thus pointwise eventual attainment of every allowed odd target,
   complete attainment of each finite target prefix at some common order,
   and divergence of the least unfilled allowed target are equivalent.
   The finite first-gap table does not prove this divergence. Neutral
   uniform fill only propagates existing words W and cannot supply a core
   over an empty target fibre. The proof keeps this quantifier boundary.

No new computational concern required another full run. The substantive
review is of proofs, types, equality boundaries, named-versus-quotient
histories, and the actual consequences of the inherited finite inventories.
