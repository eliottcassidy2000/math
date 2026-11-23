# A machine for every tournament tiling, and three different halting questions

**PROVED scoped / FINITE-EXACT.** This packet distinguishes the
finite program describing a tiling, an evolving machine stored in
tournaments, and the arithmetic search over all tournament orders.
A simulation or a shared logical form is not a reduction between RH
and Chan–Pak.

## 1. What a single tiling computes

Take a labelled tournament with a distinguished directed Hamilton
path as its frame. Its other directions are the tile bits. Compile
it to the Boolean language \(L_T\) of permutations whose consecutive
arcs point forward. For the associated counting program, its
states store a used set and last vertex, and each valid transition
adds one new vertex. Exactly \(H(T)\) accepting paths reach depth
\(n\); the program deciding any one input and the dynamic program
counting all inputs both terminate. Thus a tiling can be associated
to a machine without prescribing any dynamics on the metagraph.

This version retains far more than the path count: by the
[inherited support-parity reconstruction](../chanpak_support_birth_grammar_20261122/SUPPORT-RIGIDITY.md),
the directed Hamilton-permutation language determines even the
arcs never used in its accepting runs. This is faithful up to
alphabet relabelling. Converse reverses each accepted word, with
the appropriate alphabet relabelling for the frame. For the
converse-merged quotient, a small alternative orientation channel
is the list of directed triples up to word reversal. At each
vertex this records whether any two incident edges face in
opposite directions: exactly when the vertex is an allowable
middle in the triple. The line graph of the complete graph is
connected, so fixing one edge resolves all the others. The only
ambiguity without that choice is global converse. No claim is made
that undirecting just the long accepted words has this property.

**Exact streaming cost on a simple family.** Put a transitive
\(m\)-vertex block \(W\) between two special vertices in the
three-block cycle \(W\to b\to a\to W\), with \(m\ge1\).
The accepting words are exactly
\[
 \operatorname{inc}(S)\,b\,a\,\operatorname{inc}(W\setminus S)
 \quad(S\subseteq W),\qquad a\,\operatorname{inc}(W)\,b.
\]
Indeed \(b\to a\) is \(b\)'s only outgoing edge, so if \(b\)
precedes \(a\) they are consecutive. If \(a\) precedes \(b\),
they must be the first and last vertices. The number is
\(2^m+1\), but the **minimal complete DFA** has exactly
\(3\cdot2^m+m+2\) states. Its nonempty residual languages are:
\(2^m\) before \(b\) with both special letters left;
\(2^m\) immediately after \(b\), with just \(a\) left among
those two; \(2^m\) after \(ba\) and any sorted additional
letters, including the final empty suffix; and \(m+1\) along
the exceptional \(a\cdots b\) branch, with \(b\) left.
Different sets of remaining letters give different residuals,
and each listed type has a unique residual at each such set.
Add one dead state. Nondeterminism still needs at least \(2^m\)
states, using the pairs \((\operatorname{inc}(S)ba,
\operatorname{inc}(W\setminus S))\) as a fooling set. In contrast,
a two-way read-only checker can use logarithmic workspace by
rescanning the input for multiplicities and adjacent incidences.
A small formula for the **number** of accepted runs is not a
small streaming machine for their identities.

[Guarded macro machines](MACROS.md) make the other distinction
exact: ordinary Hamilton histories admit a small local
verifier with a sidecar, and the count mod any odd \(q\) has
a sharp finite automaton. Every residue occurs yet the same
language realizes a zero-density subset of absolute heights.

## 2. Machines that evolve inside an exact height fiber

[Padding machines](PADDING.md) give a complementary all-orders
construction. At every attained fixed \(h>1\), the left and
right singleton SCC gaps around one fixed residual core are
two independently unbounded counters. With explicit machine
control, universal counter computation lives at that same
\(H=h\), and **the first such height is 3**. In the same
gap-only model \(h=1\) is just a one-counter machine with
decidable halting. This requires insertion/deletion of singleton
components and marked control. It is not an undecidability
claim about original spectrum membership or about moving along
unmarked fixed-order complement lines.

The [actual complement-line comparison](COMPLEMENT.md) shows
where a bitwise tape ceases to be a metagraph tape. The raw
single-arc graph and the complement of all *non-frame* edges
are different operations. For sufficiently large order, a
bounded-width ordinal tape cannot take even one such global
complement step and stay in its low-width coding family.
Even two complement steps cannot implement one labelled arc
flip, irrespective of the deletion width. This does not rule
out protocols with larger overhead, re-framing or ancillary blocks.
For the *ordinary single-arc* graph, the [strong sentinel tape](STRONG-TAPE.md)
upgrades the old cube to a family of strong classes which remains
induced even after merging converse.

## 3. Which problem is the machine supposed to halt on?

- **One exact height \(h\).** This is decidable, including a
  negative answer. Strip the singleton SCCs. If the nontrivial
  components have orders \(n_i\) and counts \(h_i\), then
  \(h=\prod h_i\), and [Busch's bound](https://doi.org/10.37236/1141)
  gives \(h_i\ge3^{n_i/3}\): for \(n_i=3\) this is equality,
  and for \(n_i\ge4\) use \(5^{(n_i-1)/3}\ge3^{n_i/3}\).
  Hence a complete search needs at most \(3\lceil\log_3 h\rceil\)
  vertices and \(2^{O((\log h)^2)}\) time by subset path DP and
  orientation enumeration. The inherited
  [unit-seam decider](../chanpak_profinite_automata_20261122/BOUNDED-DECISION.md)
  gives an analogous order bound for a literal neighbouring pair.
  Nothing here claims these searches are practical at \(E_{61}\).
- **Chan–Pak.** Run the bounded decider on each allowed odd
  count, halting at the first missing one. With 7 and 21 already
  excluded, the assertion is exactly nonhalting of this specified
  machine. See the [inherited precise halting boundary](../session323_computability_and_native_descent_20261113/GUARD_CLOSURE_AND_HALTING.md).
  No indefinite “try bigger tournaments” loop is required inside
  each target. The unbounded quantifier is which target to test.
- **RH.** The [95-state Cryptids machine](https://github.com/google-deepmind/cryptids/tree/main/artifacts/riemann_95)
  comes with a kernel-checked blank-tape equivalence of halting
  to the negation of ordinary RH. Alpöge's
  [verified RH theorem](https://leanprover-community.github.io/blog/posts/riemann-hypothesis/)
  therefore implies this particular machine never halts. This
  is a mathematical composition of established statements; no
  merged Lean build is claimed. The [MIT verification record](https://math.mit.edu/~drew/grh/)
  distinguishes the ordinary statement and the checked Dirichlet
  generalization, not arbitrary automorphic GRH. A separate
  elementary checker uses [Lagarias's criterion](https://arxiv.org/abs/math/0008177)
  and rational enclosures for a strict violation, dovetailed over
  inputs and precisions. No equality test for arbitrary computable
  reals is needed.
- **An arbitrary supplied program.** It is the *uniform* family
  of program descriptions for which halting is undecidable.
  The RHS of a new mathematical equivalence must identify the
  actual program, initial state, overflow convention, and halting
  predicate. Classical Collatz is naturally an all-input
  totality question, not simply the blank-tape RH checker with
  a different label. Goldbach has a
  [separate certified checker](https://github.com/google-deepmind/cryptids/tree/main/artifacts/goldbach_21).

Using the [existing literal binary tournament tape](../session168_computational_structures_20261103/turing/NOTE.md),
the RH machine on \(L\) allocated cells uses \(480L\) vertices and
at most three arc flips per native step, while its entire slot
profile stays fixed. A replay here instantiates the official
transition table through its 3634-step initialized prefix
and counts the marked tournament bit flips without printing
the full adjacency matrix. The prefix visits 46 cells, so one
fixed allocation has 22,080 vertices, constant \(H=9^{4416}\),
and takes 8,593 arc flips. Same count at every checkpoint
does not mean the machine took no steps. This is a finite
simulation, not a re-verification of RH or a new proof of
nonhalting. For arbitrary time bounds one can allocate larger
intervals; hitting an artificial interval boundary is overflow
rather than halt.

The effective lesson from RH is not that a complicated-looking
tile should decode as a famous conjecture. A program can run
forever for a reason that is proved uniformly in its stages.
For Chan–Pak, the analogous proof obligation is a
**target-decreasing construction with a preserved ready-state
invariant**, or a uniform proof of the bounded decider's
positive outputs. A finite state chart that discards the
halting label or absolute height cannot do this. This is exactly
the gap between the [modular automata](../chanpak_profinite_automata_20261122/README.md),
the [guarded pure unit switch](../chanpak_guarded_ternary_packing_20261123/README.md),
and a genuine all-target constructor.

Replay with **python3 replay.py** and **python3 macro-replay.py**. The additional **python3 strong-tape-replay.py** checks the strong carrier and the short orientation channel.
The optional
**python3 rh-prefix.py PATH_TO_MACHINE_TXT** reads the official
Cryptids transition table linked above. The generated receipts
are mathematical computations with explicit scopes.
