# The exact path-observable quotient of a deletion–substitution grammar

**Status: PROVED deduction from the observer equivalence in [LIFTING.md](LIFTING.md); FINITE-EXACT refinement and substitution controls.** This is a scoped structural synthesis, not a new general reconstruction theorem or a claimed literature first.

## 1. Fix the operations before asking what the data determine

Work with labelled tournaments and named induced restrictions. Allow substitution: for a tournament Q on labelled slots and a nonempty tournament T_i in each slot, form Q[T_i], with all cross-block arcs prescribed by Q. Expressions are fixed independently of the input's hidden orientation. A context may use named restrictions, insert the input into a fixed slot of a fixed quotient, or use it as the quotient with fixed child tournaments. The output observable is H, the number of directed Hamiltonian vertex words.

Do **not** include orientation-sensitive endpoint tests, arbitrary adaptive source-dependent guards, or tests of ordinary automorphism groups under the same name. Those are additional observables/operations and can change the equivalence relation.

**Theorem (exact contextual quotient).** Two tournaments on the same labelled vertex set are indistinguishable by H in every expression of this full restriction–substitution grammar if and only if they have exactly the same labelled cyclic triples.

The reverse implication requires that all named triple restrictions are allowed. A smaller grammar inherits the forward implication but need not characterize C3: identity-only contexts observe merely total H.

## 2. The substitution congruence

Suppose C3(Q)=C3(Q') on a common slot set and C3(T_i)=C3(T'_i) on common named child sets. Then

    C3(Q[T_i]) = C3(Q'[T'_i]).

Proof: partition each three-element vertex set by child slots.

- All three in one child: use that child's C3 equality.
- Two in one child and one outside: the outside vertex beats both or loses to both, so the triple is never cyclic.
- Three distinct child slots: cyclicity is exactly the corresponding quotient triple's cyclicity.

These cases exhaust every triple. Named restriction plainly preserves C3 equality. Induction on an expression therefore preserves C3, and [LIFTING.md](LIFTING.md) turns this into equality of H, indeed of every labelled induced-subset path-cover count, at each output. Conversely, if a triple differs, restrict to its three named vertices: H is 3 on the cycle and 1 on the transitive triple. QED.

This does not assert a path-count computation is efficient: a compact description of an equivalence class need not make every function of that class easy to evaluate.

## 3. A persistent procedural counterexample, not only a snapshot collision

The all-height family in LIFTING.md shares one C3 hypergraph while its automorphism orders range over

    3^(3^h+j),  0 ≤ j ≤ (3^h−1)/2.

The theorem proves that no sequence of the allowed restrictions and substitutions, with H readouts, distinguishes these source choices. An arbitrarily large memory of such readouts does not restore the lost orientation information. This is stronger than showing that one finite list of statistics collides.

But the restriction matters: the labelled endpoint path matrix already distinguishes the four-vertex source/sink seeds. The full realization action in [FIBER_ACTION.md](FIBER_ACTION.md) restores the hidden choices and identifies ordinary symmetry as their stabilizer. The difference is therefore not “too little computational effort” but a precise change of retained observables.

This refines the older alignment-memory theme: first ask whether a summary is merely too small for the permitted future operations, or whether an exact quotient is intentionally invariant under every permitted operation. A richer named guard may break the latter quotient, but such a guard must be stated explicitly.

## 4. Finite counterexample search by partition refinement

For a finite operation-closed state universe, start with the partition by observed values. Refine each state by its current cell and the cells reached by every named operation (including an explicit invalid-action symbol for partial operations). Repeat until no cell splits. Induction proves that round r identifies exactly the observations obtainable by words of at most r operations; at stability, the resulting partition is the largest operation congruence contained in the initial observation partition. A split stores a distinguishing operation followed by a previously known shorter witness. This is a standard observational-equivalence construction, not a new general automata theorem.

For all labelled tournaments of orders zero through five, the public checker takes observation (n,H) and operations “delete current label i”, canonically renumbering the remaining increasing labels. The per-order cell counts are:

| refinement round | n=0 | 1 | 2 | 3 | 4 | 5 |
|---:|---:|---:|---:|---:|---:|---:|
| 0: (n,H) only | 1 | 1 | 1 | 2 | 3 | 7 |
| 1 | 1 | 1 | 1 | 2 | 11 | 107 |
| 2: stable | 1 | 1 | 1 | 2 | 11 | 208 |

The stable partition agrees with the actual labelled C3 sets, checked in both directions. One deletion is insufficient at order five: 107 cells still merge some of the 208 C3 structures. The all-order theorem explains why access to all triple restrictions eventually suffices, rather than extrapolating the finite table.

The first explicit scalar hostile in the checker's bit convention is (n,mask)=(4,2) versus (4,5). Their H values agree, but deleting original label 1 leaves H values 1 and 3. This is a named-label counterexample, not a claim that the two sources are nonisomorphic unmarked tournaments.

## 5. The common law with aliquot transport

For any observable q and operation F, an induced operation on q-values exists exactly when q∘F is constant on every q-fibre. This elementary descent criterion separates two outcomes in this session:

- Tournament C3 passes for restrictions and substitution, giving a genuine closed observational algebra.
- Scalar aliquot sum fails under dilation: s(27)=s(35)=13, but s(54)=66 and s(70)=74. The [aliquot observer theorem](ALIQUOT_OBSERVERS.md) supplies a different sufficient coordinate and a guarded inverse.

These are instances of one criterion, not an asserted isomorphism between tournaments and integers. The numbers 11,12,13 or 7,21 are not needed for the proof. To propose a stronger unifying theorem, specify the source objects, operations, quotient, and target predicate, then either prove descent or exhibit a fibre split.

## Reproduction and audit

Run `python3 check_observer_grammar.py` and `python3 -O check_observer_grammar.py` alongside `check_lifting.py`, which supplies the independently audited exact Hamiltonian-path oracle. The checker exhausts all 1,100 labelled states of orders 0..5, performs deletion partition refinement to stability, verifies exact agreement with C3, extracts a distinguishing deletion word, and checks 128 common-inflation cases through six output vertices. All validity tests use explicit exceptions.

A separate proof audit verified the three-case substitution argument and identified the necessary grammar-richness qualification for the reverse implication. The earlier independent all-subset census through six is recorded in LIFTING.md; it is not repeated by this smaller contextual checker.
