# Exact-count repairs as nonnegative fibre connectivity

**Status:** PROVED (elementary results F1–F5 below); FINITE-EXACT (literal
tournament receipts); OPEN (extension supplying arbitrary Chan–Pak targets).
These are local result labels, not reserved canon IDs. The full Chan–Pak
conjecture remains open. All tournaments in the repair bank have counts
supported on the primes 3 and 5.

## Inheritance and the changed question

The closest proved mechanism is the SCC product rule and neutral singleton
padding. The earlier
[neutral-growth classification](../../fill_tree_matrix_20261113/multiplicity/PROOF.md)
shows why births alone stay in a fixed strong-component word. The hostile
example is the isolated H=37 class in
[METAGRAPH_RESETS.md](../../adaptive_state_induction_20261103/METAGRAPH_RESETS.md):
equal scalar counts need not be connected by the native count-preserving
moves. The corrected near miss is scalar transitivity: 47→73 and 73→117 do not
imply 47→117; see
[INDUCED_OBSTRUCTIONS.md](../../metagraph_primegame_20261102/INDUCED_OBSTRUCTIONS.md).
The least-used sidecar here is the *nonnegative feasibility domain* of a
factorization, rather than its integer relation lattice.

The outside connection is the distinction between a lattice basis and a
Markov basis in algebraic statistics. The latter connects every nonnegative
integer fibre, not just the ambient integer lattice. Diaconis–Sturmfels supply
the general framework; see the primary author text
[Lectures on Algebraic Statistics, Theorem 1.3.6](https://math.berkeley.edu/~bernd/owl.pdf).
The concrete result below is proved by an elementary normal form and does not
depend on a Gröbner-basis computation.

The question is now: **for a declared repair alphabet, how much temporary
neutral material makes every existing same-target witness reachable, and what
future response can be prescribed after that repair?** This separates target
existence from transport between witnesses of a target.

Concept board: SCC product; nonnegative fibre; temporary helper; future
two-path response; prime-target obstruction. A larger repair alphabet solves
the first transport problem completely, while the two-path response makes its
limitations exact.

## Objects and permitted operations

An adjacency row records out-neighbours: bit j of row i is 1 iff i→j.
Let H(D) count directed Hamiltonian paths, with H(empty)=1. Write AB for the
ordinal sum in which every vertex of A points to every vertex of B. All paths
visit A before B, so H(AB)=H(A)H(B).

Use these **five chosen atoms**, not an asserted canonical collection of five:

| Atom | Literal rows | Order | H | P₂ |
|---|---|---:|---:|---:|
| S | (0) | 1 | 1 | 2 |
| T | (2,4,1) | 3 | 3 | 12 |
| Q | (6,12,8,1) | 4 | 5 | 32 |
| U | (14,28,24,16,1) | 5 | 9 | 84 |
| V | (14,20,24,18,1) | 5 | 15 | 108 |

Here P₂(D)=Σ_{X⊆V(D)}H(D[X])H(D[V\X]); both bins may be empty.
All five atoms are strongly connected (including the singleton convention).
V is the displayed H=15 atom; no regularity is claimed. Atom counts are checked
both by permutation enumeration and by terminal-vertex subset dynamic
programming in [check.py](check.py).

The repair objects are ordinal words in these atoms. The permitted moves are:

1. swap two adjacent whole atoms, by reversing all arcs between them;
2. replace one of the bounded blocks below by the other, retaining its vertex
   set and all arcs to the rest of the ordinal word.

These are *batch rewiring operations*. They are not single neutral arc flips,
native K moves, induced embeddings, or a promise that individual intermediate
arc flips preserve H. All displayed tournament states preserve H and order.
Adjacent swaps let any selected multiset of atoms become contiguous. The
operation extends in any ordinal outside context by the product rule. It does
not automatically extend in an arbitrary cyclic homogeneous context.

Multiplicity vectors are x=(s,t,q,u,v)∈N⁵. The map and its lost information are

\[
 A x=(n,\alpha,\beta),\qquad
 A=\begin{pmatrix}1&3&4&5&5\\0&1&0&2&1\\0&0&1&0&1\end{pmatrix},
 \qquad H=3^\alpha5^\beta.
\]

The map preserves order and H. It forgets the ordering and labels of SCCs,
their root choices, and all higher path-cover responses. Adjacent block swaps
handle the lost order for this declared repair problem. Future contexts
require an additional response sidecar.

## F1. An all-order minimal three-move repair basis

**PROVED.** Every nonempty fibre {x∈N⁵:Ax=(n,α,β)} is connected by

\[
 B_1:TT\leftrightarrow SU,\qquad
 B_2:TQ\leftrightarrow S^2V,\qquad
 B_3:QU\leftrightarrow STV.                     \tag{F1}
\]

Together with adjacent block swaps, these connect all bank words of the same
order and H. Each of the three moves is indispensable for an all-fibre basis.

**Proof.** Fix n,α,β and write

\[
 t=\alpha-2u-v,\quad q=\beta-v,\quad
 s=n-3\alpha-4\beta+u+2v.                       \tag{1}
\]

Set v*=min(α,β). If v<v* and t≥1, apply B₂ in the displayed forward direction;
q≥1, so this is legal and increases v by one. If v<v* and t=0, then
2u=α−v>0 and q≥1. Apply B₃ forward; it decreases u by one, increases v and s
by one, and creates one T. It is again legal. Repeat until v=v*.

Then apply B₁ forward while t≥2. The unique resulting multiplicities are

\[
 v=v_*,\qquad u=\lfloor(\alpha-v_*)/2\rfloor,
\]

with the other coordinates given by (1). Every step has nonnegative
coordinates and preserves Ax. Every member reaches the same normal form,
proving connectivity.

Minimality follows from the two-point fibres (n,α,β)=(6,2,0),(7,1,1),(9,2,1).
Their respective pairs are {TT,SU}, {TQ,S²V}, {QU,STV}; no other move can connect
the corresponding pair. In fact any Markov basis for this matrix must contain
these three differences, up to sign. ∎

**Boundary and consequence.** The fibre is nonempty exactly when

\[
 n\ge 3\alpha+4\beta
       -\lfloor(\alpha-\min(\alpha,\beta))/2\rfloor
       -2\min(\alpha,\beta).                    \tag{2}
\]

Necessity follows because any feasible vector reaches the displayed normal
form. Sufficiency is its direct construction. This is a complete minimum-order
and connectivity classification for this bank, not for all tournaments.

## F2. One temporary singleton is both sufficient and necessary

**PROVED.** If the allowed replacements are only B₁ and B₂, every fixed-H,
fixed-order bank-word fibre becomes connected after adjoining one temporary
singleton. The singleton can be removed at the end. Zero temporary singletons
do not suffice uniformly.

**Proof.** B₃ is simulated by

\[
 S+(QU)\ \longleftrightarrow\ SUQ
 \ \longleftrightarrow\ TTQ
 \ \longleftrightarrow\ TSSV
 \ \longleftrightarrow\ S+(STV),                \tag{3}
\]

where the first and last arrows include adjacent block swaps, and the middle
two arrows use B₁ backwards and B₂ forwards. The extra S is present in the
start and finish. The helper can be reused for every B₃ in F1. At H=45,n=9,
the two states QU and STV have no available B₁ or B₂ move in either direction,
so zero helpers fail. ∎

If vertex identities matter, the temporary singleton can be returned under
its original name. In the SU replacement, place that name in one of the new T
blocks; swap the two T blocks if necessary; use that T in the TQ replacement;
assign that name to one of its two new singleton blocks. Move the resulting
singleton to the outside. Only the declared block rewrites and swaps are used.

The obstruction is genuine despite B₁,B₂ being a Z-basis of ker(A): every
integer kernel vector has the form

\[
 (\Delta u+2\Delta v,-2\Delta u-\Delta v,-\Delta v,
   \Delta u,\Delta v)=\Delta u B_1+\Delta v B_2.
\]

The failed implication is “integer generation implies legal paths inside the
nonnegative cone.” The third relation, or one temporary singleton, supplies
the missing boundary move. This is the concrete mathematical content of the
algebraic-statistics connection.

## F3. The complete future-response geometry of these fibres

**PROVED.** For a bank word with parameters in (1),

\[
 P_2=2^{n-\alpha+\beta}3^\alpha
          (7/6)^u(9/8)^v.                       \tag{4}
\]

The feasible (u,v) are precisely the integer points satisfying

\[
 u,v\ge0,\quad v\le\beta,\quad 2u+v\le\alpha,
 \quad u+2v\ge3\alpha+4\beta-n.                \tag{5}
\]

At fixed H and n, distinct multiplicities have distinct P₂. The three forward
moves multiply P₂ by 7/6, 9/8, and 27/28 respectively.

**Proof.** A two-bin colouring of an ordinal sum restricts independently to
each component, and each bin's paths concatenate in component order. Therefore
P₂(AB)=P₂(A)P₂(B). Substitute the atom table, obtaining

\[
 P_2=2^{s+2t+5q+2u+2v}3^{t+u+3v}7^u.
\]

Equation (1) yields (4), and nonnegativity gives (5). The exponent of 7
recovers u; the exponent of 3 then recovers v. The response ratios follow by
applying the move vectors. ∎

For a prescribed future response M, factor it over {2,3,7}, and let e₂,e₃,e₇
be its valuations. There is a bank word at (H,n) with P₂=M iff it has no other
prime factors and

\[
 u=e_7,\quad v=(e_3-\alpha+u)/2
\]

are nonnegative integers satisfying (5) and
e₂=n−α+β−u−3v. When it exists, the multiplicity vector is unique. F1 constructs
it from any other witness in the fibre; F2 does so using B₁,B₂ and one helper.
Thus this is an actual target-preserving, repeatable repair theorem for a
specified response, rather than just a connectivity assertion.

## F4. The consumer matters: exact repair versus context-safe replacement

Let C₁[D] add vertices b,a with D→b→a→D. This tournament is strongly connected
when D is nonempty. The established context calculation in
[METAGRAPH_RESETS.md](../../adaptive_state_induction_20261103/METAGRAPH_RESETS.md)
gives

\[
 H(C_1[D])=3H(D)+c_2(D)=H(D)+P_2(D),            \tag{6}
\]

where c₂ counts ordered pairs of nonempty paths covering D. The two empty-bin
terms contribute 2H(D) to P₂. The literal checks independently confirm (6).

**PROVED.** At fixed H and n, a desired next C₁ count R is attainable by these
repairs exactly when M=R−H passes the inversion test after (5). Every such
response can be reached from every other witness of this bank fibre. No
nontrivial change of multiplicities can preserve that next response.

The smallest indispensable third move gives the explicit boundary:

| Word | n | H | P₂ | H(C₁[word]) |
|---|---:|---:|---:|---:|
| QU | 9 | 45 | 2688 | 2733 |
| STV | 9 | 45 | 2592 | 2637 |

The full same-(H,n) fibre consists of those two multiplicities. Target 2735
under C₁ is impossible in this bank fibre. The successful repair changes future
capability; it cannot be substituted into an already fixed arbitrary strong
context while preserving that context's total H.

**Common-helper no-go.** If P₂(D)≠P₂(E), appending the same ordinal helper R
cannot make P₂(DR)=P₂(ER), because P₂(R)>0 and multiplication cancels. More
generally, the full cover polynomial

\[
 \pi_D(z)=\sum_k d_k(D)(z)_k
\]

(d_k is the number of unordered k-path covers) is multiplicative under ordinal
sum. Its equality is the existing full-profile context congruence. Since Z[z]
has no zero divisors, a common ordinal helper cannot hide a mismatch in that
polynomial either. F2 preserves the scalar target while deliberately changing
the response; it is not a catalyst for equality of all future responses.

## F5. What generalizes to every fixed finite atom bank

**PROVED.** For any finite list of actual tournament atoms, there is a finite
library of bounded block replacements that connects every ordinal-word fibre
of fixed order and H, after allowing adjacent block swaps. The bound depends
on the bank and is uniform over the number of atoms in a word.

**Proof.** Form a matrix whose first row is atom order and whose other rows
are valuations at the finitely many primes dividing the atom counts. Equality
of its images is exactly equality of order and H. Consider all nonzero pairs
(a,b)∈N^r×N^r with Aa=Ab and disjoint supports. The set of componentwise minimal
pairs is finite by Dickson's lemma. Use the corresponding replacements as the
library. For two equal-image vectors, remove their common coordinates. If a
nonzero difference remains, choose a minimal pair beneath the remaining
positive and negative parts. Apply that replacement; all coordinates stay
nonnegative and the L¹ distance to the target strictly decreases. Repeat.
Each relation is realized by literal ordinal blocks on the same number of
vertices, since order was a matrix row. ∎

This is the finite-bank instance of the Markov/Graver-basis mechanism. It does
not provide useful move sizes, a small helper budget for a preselected library,
or new target counts. Extra *multiplicative* response measurements can be
added as valuation rows. For example, adding P₂ to this five-atom matrix makes
the multiplicity map injective, so the remaining fibre transport consists
only of changing block order/labels. A finite bank supplies no new prime
factors of H; this alone prevents it from proving Chan–Pak.

## Exact verification and next decisive question

Reproduce with:

```text
python -B 05-knowledge/results/strategy_fibres_absorption_20261113/fibres/check.py
```

The [receipt](receipt.json) records all literal atoms and exchanged tournaments.
The finite universe is every nonnegative integer fibre with 0≤α,β≤12 and
0≤n≤60: 4,672 nonempty fibres, 35,082 states. Independent breadth-first
connectivity agrees with the constructive normal form. The two-move alphabet
fails on 252 fibres in that universe. The code also checks the order threshold,
response inversion, atom strong connectivity, permutation-vs-DP atom counts,
literal C₁-wrapper counts, and ordinal products. All validations use explicit
exceptions, so they also run under Python optimization. These finite checks
validate the implementation; the all-order assertions use the proofs above.

**What changed strategically:** ask for a *consumer-specific transport
theorem*. A repair may be allowed to alter P₂ precisely because a later step
needs a different P₂. Conversely, demanding all-context equivalence can make
the permitted fibre trivial. The required invariant must be chosen from the
actual next construction, not from a universal scalar slogan.

**The next hostile test was executed:** the absorption lane's distinct H=25
witnesses at the overlap of [23,25] and [25,31] are joined by a neutral
two-arc switch. It absorbs a singleton into a strong six-core, producing a
strong seven-core. That endpoint is precisely the orbit lane's M6² rank-one
reset, up to an explicit rooted isomorphism. Its equality fails under *every*
nonsingleton single-port module substitution, by sign-definite exact cover
coefficients; every two-port identical-transitive2 probe also fails. See
[H25-BRIDGE.md](H25-BRIDGE.md) for the positive connection, complete
single-port obstruction, and corrected next question. Scalar interval overlap
still does not itself glue witnesses; the literal two-arc move does.

The broader OPEN problem is a family of such consumer-specific repair banks
whose strong atoms supply the missing prime targets and whose required helper
budget is controlled. F1–F5 solve the finite-bank transport question, not that
prime-supply question or the 7,21 classification.
