# A native 4+6 prime cut, its exact selector, and the missing path interface

**Proved in the stated finite height window.** This is a structural distinction between a component selector and a Hamiltonian selector, not a moonshine derivation of the square-sum graph or a Chan–Pak proof.

## 1. Three genuinely different partitions

The [October 22 Monster-level note](../monster_levels_square_sum_7_21_20261022.md) gives the inherited exact prime-support filtration
\[
\pi(\mathrm{Fi}_{22})=\{2,3,5,7,11,13\},\quad
\pi(B)\setminus\pi(\mathrm{Fi}_{22})=\{17,19,23,31,47\},\quad
\pi(M)\setminus\pi(B)=\{29,41,59,71\}.
\]
It also separates the order quotient from the involution class index \([M:2.B]\); the latter, not \([M:B]\), is the group-theoretic index. These are inherited group/level facts, not rederived here. [Session 32](../session32_square_monster_20261023/SYNTHESIS.md) finds the full 15-prime square-sum tree, identifies the two components on Baby support, and exhibits failure of a same-cardinality block map. [THM-4860/4861's prime square-sum work](../fermat_monster_hamiltonian_counts_prime_classes_20261102/README.md) supplies the inherited residue restrictions. None identifies a prime-support classification with tournament Hamiltonian counts.

The first ten primes are
\[
S=\{2,3,5,7,11,13,17,19,23,29\}.
\]
They are also the first ten entries of the *Monster* prime list. They are not the first ten *Baby* primes: the latter replace 29 with 31, and the full Baby support has eleven primes. On \(S\), the group filtration has sizes \(6+3+1\); Baby membership alone gives \(9+1\).

Give a prime set the native square-sum edges. The first ten's graph is
\[
(3-13-23-2-7-29)\ \sqcup\ (5-11)\ \sqcup\ (17-19).
\]
The **6** is the connected component of the marked prime 2; the **4** is its complement, itself two disjoint pairs. This is a real \(6+4\) cut, but not a division into two connected blocks. On replacing 29 with 31 (the first-ten Baby list), the root component loses its pendant 29 and the pair \(5-11\) gains the pendant \(31-5\); hence the components have sizes \(5,3,2\). The native prime 4+6 is therefore not stable under that moonshine-inspired selection. For comparison, the parity bipartition on the original first ten, with 2 joined to the \(1\bmod4\) side, is \(5+5\), another distinct invariant.

## 2. Exact selection at height 71

Allow any subset \(U\) of additional primes \(29<p\le71\). Call a choice **connecting** if the ten marked roots \(S\) all lie in a common component of the graph on \(S\cup U\); it need not use every optional vertex in that component.

**Theorem.** The exact criterion is
\[
47\in U\quad\land\quad
\bigl(53\in U\;\lor\;\{41,59\}\subseteq U\bigr).
\]
Thus the inclusion-minimal choices are exactly \(\{47,53\}\) and \(\{41,47,59\}\). If choices are restricted to Ogg/Monster primes, the second is the only possibility. In the ordinary prime stream the marked roots first coalesce at 53; in the Ogg stream they first coalesce at 59. If only the added primes are required to lie in Baby support, the ten marked roots do not coalesce. Requiring the whole chosen graph, not just the ten roots, to be connected additionally excludes \(37,43,67\).

**Proof.** All edges incident with an added prime can be listed by subtracting each square up to \(11^2\) from its label and keeping only smaller primes. The entire lower-neighbor ledger is:

| added prime | smaller neighbors |
|---|---|
| 31 | 5 |
| 37 | none |
| 41 | 23 |
| 43 | none |
| 47 | 2, 17 |
| 53 | 11, 47 |
| 59 | 5, 41 |
| 61 | 3 |
| 67 | none |
| 71 | 29 |

No other added prime can attach the \(17-19\) component, so 47 is forced. With 47 present, the only ways to connect the \(5-11\) component to the root component are the path \(11-53-47\), or the path \(5-59-41-23\). All remaining optional vertices are leaves or isolated within this height window. These observations give both necessity and sufficiency, including the exact minimal choices. Removing Monster-only \(29,41,59,71\) from full Monster support leaves Baby component sizes 8 and 3; in the full tree, \(41,59\) are the connecting path and \(29,71\) extend a pendant path. The old Session 32 full-tree statement is consistent with, but does not replace, the arbitrary-subset criterion. ∎

## 3. Why connected selection is still not a path selector

**Obstruction.** No subset \(U\) of the allowed added primes yields a Hamiltonian path on \(S\cup U\), even if it satisfies the exact connecting criterion.

There are at least three different forced degree-one vertices. The marked vertex 19 always has only neighbor 17. If 61 is absent, 3 has only neighbor 13; if 61 is present, 61 is a degree-one vertex attached to 3. Independently, if 71 is absent, 29 has only neighbor 7; if 71 is present, 71 is a degree-one vertex attached to 29. A spanning path has only two endpoints, so it cannot cover all three leaves. The proof applies even if a different connector is picked and even if isolated optional vertices are excluded.

The structural lesson is explicit: collapsing the three initial components to contacts gives a complete Boolean selector for connectedness, but forgetting endpoint capacity loses the decisive obstruction to a spanning path. The additional resource is not a denser collection of connecting primes; it is a different boundary choice that removes or joins an unavoidable leaf.

## 4. Contrast with the other 4+6

The [size-25 overlap packet](../session321_selective_returns_and_fano_wiring_20261113/SQUARE_OVERLAP_PORTS.md) splits ten *Hamiltonian paths* into counts 4 and 6 using two marked ears \(11-25-24-1\) and \(11-25-24-12\). Here 4 and 6 count *vertices* in different components of a prime support. Their common structural pattern is a marked binary selection whose answer changes when the native interface changes. There is no map identifying the two finite sets by numeric labels, or taking Monster membership to an ear. In fact one carrier has ten full paths and the other can have none even after every allowed connectivity repair. This provides a direct hostile to using matching split cardinalities as a constructive equivalence.

For a tournament reservoir, the transferable step is the separation of a **selection certificate** (which marked choices meet the desired scalar condition) from a **continuation certificate** (which simultaneously retain endpoint readiness). The former is exact and simple in this prime model; the latter fails for every member of its fiber. Thus one must prove the analogue of the continuation gate, rather than assuming an exact connectivity or coverage selector proves it.

Replay: python3 replay_prime_selection.py (checks every optional subset).
