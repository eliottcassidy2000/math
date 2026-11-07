# Positive ribbon compiler and a marked local-edit observer

**Status: PROVED elementary application of inherited constructions; FINITE-EXACT executed fixed points. Independent audit PASS WITH TWO SCOPE CORRECTIONS, both applied.** This is a reusable exact compiler, not a new Hamiltonian-path spectrum theorem, a new finite-state method, or a resolution of Chan–Pak. The all-bound theorem below is the standard truncated-nonnegative-semiring argument, already implemented for other tournament grammars in Session203 and Session216, explicitly applied here to the repository's ribbon matrices. The machine and Rule 30 provenance is in `PROVENANCE.md`.

## 1. Inputs, source identity, and exact observable

Fix a ribbon span `w>=1` and an integer budget `B>=1`. The alphabet consists of subsets `a` of `{1,...,w}`. A word `a_1...a_n` gives the tournament on the **ordered labelled** vertices `1,...,n` in which `q -> q-s` for `s in a_q`, `s<q`, and every other arc points forward. Bits referring to nonexistent vertices do nothing. An actual graph can be recovered from the word without a counting oracle.

We use THM-4880(A)'s finite nonnegative profile transfer

    x(ua) = M_a x(u),       H(u) = f x(u),

where `x(empty)=e_empty`, and `f` sums the one-fragment states. The shipped `sources/ribbon.py` is an unchanged snapshot from commit `38fdf9f234159aaedb2d0ba1e0f7bcdefb49e817`, path `04-computation/experiments/ribbon_spectrum_20261104/ribbon.py`. It has `d=7,27,119` profile states for `w=1,2,3`. Its empty-word readout is zero, so this packet explicitly starts at the one-vertex word and discusses **nonempty tournaments**; the shortest nonempty witness for `H=1` has order one. No empty-word convention is silently changed.

The implementation searches all tournaments in this ordered bounded-span family, not only strong tournaments. It does not quotient by graph isomorphism. The signed reduced representations of dimensions 14/62 are deliberately not used for saturation.

## 2. The inherited saturation lemma

Put `C=B+1` and `tau(z)=min(z,C)` for nonnegative integers, entrywise for vectors and matrices. For `a,b>=0`,

    tau(a+b)=tau(tau(a)+tau(b)),
    tau(ab)=tau(tau(a) tau(b)).

For multiplication, a zero factor annihilates on both sides; otherwise an overflowing factor forces an overflowing product. Hence `tau` is a homomorphism into the finite semiring `{0,...,C}` with saturated addition and multiplication. In particular,

    tau(M_a tau(x)) = tau(M_a x),
    tau(f tau(x)) = tau(f x).

This is exactly Session135 `acceleration/THEORY.md` §4, specialized to the nonnegative linear ribbon maps. It was already used as a sound saturation principle in THM-4878. The nearest executable predecessors are `session203/saturated-profile-membership/PROOF.md` and `05-knowledge/results/session216_overnight_constructions_20261107/inverse/PROOF.md`: they already prove capped BFS termination, shortest witnesses, resource-safe incompleteness, and output-monotone pruning for their positive tournament grammars. We do not claim the theorem or generic algorithm is new; the release specializes it to the ribbon profile and adds the marked two-sided edit observer.

A coordinate cap cannot instead be put on a signed representation: with `C=6`, `10-9=1` but `tau(10)-tau(9)=0`. Nor does it preserve arbitrary TM/FRACTRAN microsteps: `12` and `18` have the same cap at five, but after eight decrements their caps are four and five.

## 3. All-length termination and exact shortest selection

**Theorem.** For every fixed `w,B`, breadth-first exploration of the capped profile vectors, starting with the one-vertex vector and using all ribbon letters, terminates. It decides whether every `h<=B` is attained by a nonempty span-`w` ribbon tournament. Keeping the first incoming letter and predecessor of each new vector supplies a shortest-order witness within this family for every attained `h<=B`.

A valid worst-case bound is `(B+2)^d` states and `2^w(B+2)^d` transitions. Every attained target has a witness of order at most `(B+2)^d`. These bounds are deliberately crude, exponential in the profile dimension, and polynomial only in the *numeric* budget at fixed width. There is no logarithmic-in-`h` selector, efficient uniform-width algorithm, or all-odd coverage conclusion.

**Proof.** Semiring compatibility proves by induction that the state of every word is exactly `tau(x(word))`, and its readout is `min(H(word),B+1)`. There are at most `(B+2)^d` states. Equal capped states have equal capped futures under every continuation. Therefore each new state can retain only its first predecessor, with no loss of any target up to `B`. Breadth-first order makes that predecessor word shortest. If a shortest witness revisited a capped state, removing the intervening block would preserve its target readout, a contradiction; this gives the length bound. The first letter is always equivalent to the empty subset because there are no earlier vertices. Thus starting with that one-vertex word loses no nonempty tournament. QED.

**Safe pruning used by the code.** Reject a state once its readout is `B+1`. Rédei insertion gives `H(T)<=H(T+v)`: every directed Hamiltonian path admits insertion of `v`, and deleting `v` recovers the old path, so choosing one insertion per path is injective. Consequently no appended continuation can return below the budget. This is pruning along vertex addition, not a claim that arbitrary local edits are monotone.

If a further regular language of words is wanted, take the product with its DFA. Initialize with every permitted first-letter DFA successor (or start at the empty word with an explicit nonempty marker): arbitrary DFAs may distinguish first letters even though their singleton tournaments coincide. The same argument gives at most `r(B+2)^d` nonempty states for an `r`-state DFA. THM-4880(B)'s strong-connectivity DFA is one option, but this optional product is **not implemented or counted here**.

**Difference from THM-4881(C).** That theorem already proves the all-length span-2 `H<=99` result by the exact fixed point of 2,865 uncapped reduced states at length twelve. The current cap compiler supplies a guaranteed finite state universe for every budget and span; its `B=99` run uses 5,102 positive capped profile states. Neither the `7,21,47` exclusions nor this range of spectral coverage is new. It is not asserted that the cap reduces the practical search.

## 4. Controlled local writes and the needed two-sided sidecar

A fixed ordered ribbon word is also a finite tape over an alphabet of `2^w` symbols. Replacing letter `a_j` by `b` toggles exactly the arcs `{j,j-s}` whose allowed bits differ, at most `w` arcs. One existing-offset bit flip changes exactly one arc; nonexistent offsets change none. Thus a TM write on this alphabet has a literal **controlled marked local arc-rewrite** image. Retain the ordered vertices, head position, control, program, and finite-window boundary. This is not autonomous tournament dynamics and not the full-order native `K` operation from Session65.

For the factorization `u a v`, retain

    L = tau(M_u x(empty)),       R = tau(f M_v).

Here the chronological product is `M_u=M_last ... M_first`. Then every marked replacement has the exact bounded output

    min(H(u b v),B+1) = tau(R M_b L).

This follows from the same semiring homomorphism. It is a direct useful merge of Session65's marked tape interface with the actual count-preserving transfer of THM-4880, rather than the tape observer `H(E(w))=3^wt(w)`.

The current scalar count does not suffice. Already at span one the two words `000` and `001` both have `H=1`. Change their common middle letter `0` to `1`: the outputs `010` and `011` have `H=1` and `H=3`. The same marked write therefore has different outcomes on the same scalar count. The right context is load-bearing.

This is an exact observer for a **supplied** local write, not a theorem that any desired count can be reached by such writes. Moving the head does not generally update `R` by matrix inversion, since matrices can be singular. Retain the word, a prefix/suffix deck, or a tree of block matrices and recompute the appropriate contexts. Entrywise-capped block matrices form a finite compositional monoid, so a balanced product tree supports a one-letter replacement in `O(log n)` dense matrix products (naively `O(d^3 log n)` saturated scalar operations), with `O(nd^2)` storage. This standard data-structure observation is not implemented here. Its storage grows with tape length; it supplies no finite-state simulation of an unbounded TM.

If a local edit temporarily leaves the `H<=B` region and later returns, the append-BFS pruning must not be reused as an edit-reachability rule. The two-sided observer remains sound when recomputed from the actual word.

## 5. Executed controls and results

`verify.py` reconstructs the profile matrices from the pinned source, checks nonnegative coefficients, and independently reconstructs adjacency masks and counts every emitted witness by subset/endpoint DP. It also runs 685 small-word graph/profile/cap checks and 220 canonical single-bit-edit checks against independent graph DP. Explicit hostiles cover scalar count, signed capping, and decrement microsteps. These are regression tests of the proof, not its all-length basis.

| Span | Budget | Capped states | Transitions | Longest shortest state word | Longest output witness | Missing positive odds |
|---|---:|---:|---:|---:|---:|---|
| 2 | 99 | 5,102 | 20,408 | 12 | 8 | 7, 21, 47 |
| 2 | 199 | 14,618 | 58,472 | 13 | 9 | 7, 21, 47 |
| 3 | 25 | 7,025 | 56,200 | 12 | 6 | 7, 21 |

The receipts include every shortest output witness as a list of letter bitmasks and a deterministic digest of the BFS state stream. Normal and optimized `B=99` executions are byte-identical. `--max-states` is a resource guard: hitting it raises an error and produces no completeness receipt. The theoretical guarantee does not guarantee a practical budget is enough.

Reproduce from any current directory:

    python3 -B /path/to/verify.py --span 2 --bound 99 --out span2_B99.json
    python3 -B /path/to/verify.py --span 2 --bound 199 --out span2_B199.json
    python3 -B /path/to/verify.py --span 3 --bound 25 --out span3_B25.json
    python3 -B -O /path/to/verify.py --span 2 --bound 99 --out span2_B99_optimized.json
    cmp /path/to/span2_B99.json /path/to/span2_B99_optimized.json

Read source hashes as byte provenance, not as independent mathematical verification. The source profile transfer is inherited and imported; the emitted tournament witnesses and bounded edit controls use the separate DP path. The independent auditor reproduced the span-two budget99 closure, independently closed span-three budget47, and checked the quotient, termination, shortest-witness and edit proofs. Its report is delivered separately; this packet claims only the executed ranges in its own receipts.
