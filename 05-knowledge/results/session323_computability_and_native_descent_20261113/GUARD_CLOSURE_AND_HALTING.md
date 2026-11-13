# Guard closure, cyclic priority, and the precise Chan–Pak machine boundary

**Proved:** the divisibility closure representation, exact priority inverse, and optimal shielding in the stated guard model. **Inherited:** the bounded tournament-core theorem and the compiled-machine translations. No new Collatz or Chan–Pak resolution is claimed.

## 1. Historical connections and the actual reduction

The [machine/automata/FRACTRAN/tournament dictionary](../machines_automata_fractran_tournaments_dictionary_20261101/README.md) distinguishes exact machine simulations, actual finite path certificates, and a deliberately limited analogy between path guards and applicable fractions. [Session 135](../session135_fractran_tm_structure_20261101/README.md) gives explicit register/FRACTRAN/Turing translations and complexity costs. A standard FRACTRAN program with reduced denominators \(d_i\) is an exact partial residue-controlled multiplicative map: its chosen fraction is determined by the least matching denominator, hence by the state modulo \(\operatorname{lcm}_i d_i\). This describes one step; the residue need not determine the next residue, so a finite selector is not a finite-state abstraction of the whole orbit. The repository's [guard precision note](../stability_priority_moduli_20261102/FRACTRAN_PRECISION.md) gives a concrete difference between knowing prime support and knowing valuations.

For Chan–Pak there is an actual counterexample machine, stronger than only a witness search. [Session 133's core theorem](../session133_fractran_spectrum_refinement_20261101/core_reduction/README.md) bounds the non-singleton core of every tournament with count at most \(M\) by
\[
K(M)=\max_{\substack{b=0,1,2\\ 5^b\le M}}
 \left(4b+3\left\lfloor\log_3(M/5^b)\right\rfloor\right).
\]
Deleting singleton strongly connected components preserves the count by ordinal factorization. Retain the convention \(H=1\) for the empty residual. Thus to decide exact membership at \(h>1\), enumerate the finitely many tournaments through order \(K(h)\), count their Hamiltonian paths exactly, and return a yes or a no. This is decidability, not a fast algorithm in the bit length of \(h\). A fixed TM can enumerate positive odd \(h\ne7,21\), run this complete decision for each, and halt with the first absent \(h\). Its nonhalting is equivalent to Chan–Pak (using the inherited exclusions). The general TM-to-FRACTRAN compilation transports this **particular search**, not the proof of its nonhalting.

For ordinary Collatz, the direct fixed-input machine halts exactly if its supplied input reaches 1; universal convergence has the natural form \(\forall n\,\exists t\), a totality statement. That is distinct from Chan–Pak's \(\Pi^0_1\) normal form after the core bound. An arithmetical-hierarchy comparison is not an impossibility theorem for a theorem-specific equivalence between two fixed statements. It does rule out treating *general* machine totality as equivalent, under an effective uniform reduction, to single-machine nonhalting. For reference, the earlier [Camion–Busch comparison](../camion_busch_gaps_polyhedra_collatz_20261002.md) already identifies several missing bridges. Neither a shared syntax nor a literal reduction of finite checks gives a geometric reduction between the conjectures.

**No autonomous guard-visible residue quotient for a genuine denominator.** Suppose a positive modulus \(M\) is required both to determine the chosen rule and to determine the residue of the next input, whenever a rule is selected. If the reduced fraction \(a/d\) is selected at \(n\), it is also selected at \(n+M\). Integrality of both outputs forces \(d\mid M\). The output difference is \(aM/d\); it is zero modulo \(M\) iff \(d\mid a\), hence iff \(d=1\). Conversely if all selected rules have denominator one and the current residue determines the rule, multiplication does descend. This excludes every fixed modulus for a guard-visible autonomous dynamics whenever some reachable rule has a genuine denominator. It does not exclude noncongruence abstractions, abstractions restricted to a special invariant input set, or state enlargements that retain higher valuations.

## 2. Exact applicability closure

Fix reduced positive denominators \(d_1,\ldots,d_m\), and put \(D(n)=\{i:d_i\mid n\}\) for positive inputs. For \(S\subseteq[m]\) write \(L_S=\operatorname{lcm}(d_i:i\in S)\), with \(L_\varnothing=1\). Then
\[
\operatorname{cl}(S)=D(L_S)
\]
is a closure operator, and the attainable applicability sets are exactly its fixed sets. To see idempotence, the denominators of every added index divide \(L_S\), so the lcm does not change. If \(S=D(n)\), then \(L_S\mid n\), so every index added by closure is already in \(D(n)\). Conversely every closed set is witnessed by its own lcm.

**Representation theorem.** Every finite intersection-closed family \(\mathcal C\) on \([m]\) containing \([m]\) arises in this fashion with squarefree denominators. For each \(C\in\mathcal C\) take a distinct prime \(p_C\) and let \(d_i=\prod_{C\not\ni i}p_C\). A prime \(p_C\) divides \(L_S\) iff \(S\not\subseteq C\). Hence \(d_j\mid L_S\) iff there is no closed set \(C\) with \(S\subseteq C\) and \(j\notin C\). This is exactly membership in the intersection closure of \(S\). Full ground set makes that intersection well-defined. The representation is not promised polynomial-size if the closed sets are not already listed.

This is a genuine native selection compiler: the same closure family is realized on all positive inputs, with an explicit smallest divisibility witness for each attainable set. It is not a compiler of an arbitrary state-transition diagram, because the fraction numerator and future residue constraints were not specified.

## 3. Total cyclic priority and when it linearizes

A tournament on the rule indices chooses the unique source of the induced tournament on the nonempty applicable set, when one exists. No fraction applicable is normal halting; a nonempty set with no source is a conflict failure. Conflict freedom is equivalent to having a source on every nonempty closed set, a finite condition by §2.

For a conflict-free choice rule, define the **revealed priority digraph** by drawing an arrow from the selected rule of each attainable nonempty set to every other applicable rule in that set. There exists one fixed linear priority on the same rule list making the same choices at all inputs iff this digraph is acyclic. Every successful linear order must obey all its arrows; conversely a topological ordering puts the specified choice before all its competitors on each set.

A concrete witness uses the reduced denominators \(2,3,5,30\), with numerator one. Orient the first three as a directed cycle, and let the fourth beat all three. The triple is applicable only when 30 is applicable, so every nonempty attainable set has a source. Each pair of the first three is attainable alone (witnesses \(6,10,15\)), exposing every arrow of the priority cycle. The revealed priority digraph therefore has a cycle, and no fixed linear order on the same four guards realizes the same rule choices. All four rules actually win on some input (for example \(2,3,5,30\)). Each chosen update strictly decreases the positive integer, so the entire resulting system still halts. Thus neither cyclic priority nor an irreducible difference from one linear ordering gives computational universality. Conversely this says nothing against simulations using additional rules or state.

## 4. Minimal triangle shielding with a protected sublanguage

Start with one original rule per distinct prime \(p_i\) and an arbitrary tournament on these rules. The original applicable sets are then all subsets. Add **shield rules** with reduced denominators that are squarefree products of original primes. Require shields never to become applicable when the original applicable set induces a transitive tournament. Whenever shields occur, let them beat all original rules, with a linear order among shields.

**Theorem.** The minimum number of shields making the entire priority rule conflict-free is exactly the number of directed triangles in the original tournament. One shield per triangle, with denominator the product of its three guard primes, achieves the optimum.

For sufficiency, a nontransitive applicable subtournament contains a directed triangle, and therefore enables at least one shield. The first enabled shield wins. A transitive nonempty applicable set triggers no shield and has its original source. For necessity, on the input \(\prod_{i\in C}p_i\) of an original directed triangle \(C\), an enabled shield must have its support contained in \(C\). The protection condition rules out every proper subset of \(C\), since all those induce transitive tournaments. A distinct shield for the full support \(C\) is therefore necessary for each directed triangle. All original pair choices remain exposed under this completion; if the original priority is nontransitive, no fixed order on the extended rule list can realize all its winners.

This completion is the exact analogy to **repairing a failed local selector while preserving every legal protected case**. It is stronger than merely noticing that odd many tournament guards may fire. It also isolates the limitation: the shield repairs a finite guard conflict and controls no downstream endpoint, exponent budget, or termination unless those are separately proved. A Chan–Pak repair chart needs the extra native state that proves its output is both an exact count and a legal next source. The guard closure can model finite availability; it cannot manufacture that invariant.

Replay: python3 replay_guard_closure.py checks a finite universe of closure systems and priority examples. The theorems are not inferred from the replay.
