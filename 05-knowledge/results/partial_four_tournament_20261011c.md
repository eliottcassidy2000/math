# The four-arc square is exactly the marked four-state tournament section

**Status: PROVED finite structural statements; FINITE-EXACT complete controls.**
The recognized four-state section is inherited. The new work here is its
recognition from the specified partial degrees and missing pairs, its complete
square-orientation comparison, and the explicit boundary between orientation
completion and the previous graph-colouring obstruction.

## 1. Recovering the exact object

A tournament on four vertices has six arcs. The user's four-arc object is a
**partial tournament**: the specified arcs are actual comparisons, while two
pairs remain unoriented. Missing arcs are not ties assigned an arbitrary winner.
The degrees here count only the specified outgoing arcs.

Suppose these degrees are $\{0,1,1,2\}$ and the missing pairs join the
degree-zero vertex to one degree-one vertex and the degree-two vertex to the
other degree-one vertex. The two missing pairs are a matching, so the underlying
four-edge graph is a square. The degree-two vertex is its source $a$, and the
degree-zero vertex is its sink $d$. They are adjacent: the stated missing pairs
exclude the opposite-source/sink case. Thus their edge is $a\to d$.

Let $b$ be the other neighbour of $a$ and $c$ the remaining vertex. The degree
conditions force

$$
 a\longrightarrow b\longrightarrow c\longrightarrow d,
 \qquad a\longrightarrow d, \qquad
 \text{missing pairs }ac,\ bd. \tag{1}
$$

Conversely (1) has exactly the stated degree and missing-pair pattern.
The frame $(a,b,c,d)$ is unique, so this partial orientation has no nontrivial
automorphism. It is acyclic and its Hamiltonian path is unique.

The closest proved mechanism is the $c=0$ half-section in
[four tournament states and recursive storage](tournament_recursive_four_state_20261004.md),
§1, recovered there from
[HYP-3244, half-tiling provenance](../hypotheses/HYP-3244-tiling-half-tiling-interlocking-recursions.md).
Its mathematical statement is already proved in the former packet; this note
does not promote the historical hypothesis as a new dependency. That packet
also preserves the MISTAKE-034 correction to the older Fourier half-band claim,
which is not used here.

The concept board is: partial comparison, marked completion bits, unlabelled
tournament type, binary cycle space, and actual colouring ports. Each next map
has a different information budget.

## 2. Exactly four completions, with a sharp degree-only hostile

Let $x=0$ mean $a\to c$ and $x=1$ mean $c\to a$; let $y=0$ mean $b\to d$
and $y=1$ mean $d\to b$. Completing (1) gives:

| $(x,y)$ | Type | Sorted full outdegrees |
|---|---|---|
| $(0,0)$ | $T$: transitive | $(0,1,2,3)$ |
| $(1,0)$ | $+$: a directed triangle over a sink | $(0,2,2,2)$ |
| $(0,1)$ | $-$: a source over a directed triangle | $(1,1,1,3)$ |
| $(1,1)$ | $S$: strongly connected | $(1,1,2,2)$ |

These are all four four-vertex tournament types. Indeed, a tournament with no
directed triangle is transitive. Given a directed triangle, the fourth vertex
either dominates it, is dominated by it, or has mixed incidence. Mixed incidence
gives a directed four-cycle. Its four possible diagonal orientations are
isomorphic under rotations of the cycle, so they give one strong type. The
companion independently verifies the list by all 64 labelled tournaments and
all 24 vertex permutations.

The missing-pair condition is essential. The opposite-source/sink orientation

$$
a\to b\to c,\qquad a\to d\to c
$$

has the same partial outdegrees $\{0,1,1,2\}$, but its missing pairs are $ac$
and $bd$, joining degrees $2$ to $0$ and $1$ to $1$. Its completions give
two transitive and two strong tournaments, with neither diamond type.

In fact every orientation of a square belongs to exactly one of these four
cases:

| Square orientation | Completions, with multiplicity | Count on one fixed labelled square |
|---|---|---:|
| Adjacent source and sink, as in (1) | $T,+,-,S$ | 8 |
| Opposite source and sink | $T,T,S,S$ | 4 |
| Directed four-cycle | $S,S,S,S$ | 2 |
| Alternating two sources and two sinks | $T,T,T,T$ | 2 |

To prove completeness, either all four partial degrees are one, giving the
directed cycle, or there is a source and sink. If there are two sources there
are two sinks and orientation is alternating. Otherwise the source and sink
are adjacent or opposite. The displayed completions follow by orienting the
two remaining pairs. Thus **among oriented squares**, having four distinct
completion types is equivalent to the user's specified pattern. No such
classification is asserted for arbitrary larger partial tournaments.

## 3. A valid two-port gadget, with its loss boundary

For this fixed partial object, the unlabelled completion type is a lossless
encoding of the two bits. The two permitted named flips act as

$$
\begin{array}{c|cc}
 & \text{flip }ac & \text{flip }bd\\ \hline
T&+&-\\
+&T&S\\
-&S&T\\
S&-&+
\end{array}
$$

They commute and square to the identity. Hence this is a genuine four-state
$\mathbb F_2^2$ interface, with a type-to-bit decoder, for contexts which retain
the partial source and those two ports.

The script exports:

- **frame(partial)**: authenticate the stated degree/missing-edge pattern and
  recover its unique frame;
- **complete(partial,(x,y))** and **read(partial,full)**: reconstruct and authenticate
  the same marked completion;
- **from_type(partial,type)** and **flip(partial,full,port)**: the exact section
  decoder and the two legal named updates.

This is more precise than storing an arbitrary four-vertex tournament type.
A relabelled transitive tournament has the same type but can reverse one of
the supplied partial arcs; **read** rejects it. External attachments therefore
still require actual vertex identification, or an explicit transport of the
frame. Type equality alone does not identify boundary ports.

Nor does this four-state quotient admit every possible free-arc flip. In the
inherited three-bit fixed-path cube, reversing $ac$ and $bd$ gives a strong
presentation, and reversing only $ad$ gives another strong presentation.
Reversing $ac$ again sends the former to the source-over-cycle type and the
latter to another strong presentation. The same type has two different next
types. This is the earlier predictive-memory hostile, and the repair remains
to retain the full presentation when the third flip is allowed.

Actual tournament converse, followed by path reversal to return to frame (1),
swaps $(x,y)$ and hence swaps the two diamond types while fixing $T,S$.
Complementing both bits sends $T$ to $S$; it is not converse.

No generic arithmetic certificate is supplied by these two bits. As in the
earlier recursive codec, any interpretation as stored proof data must retain
its parser, original source, guards and terminal proof.

## 4. What the preceding K4 cycle-basis hostile does and does not say

Use the labels $(a,b,c,d)=(0,1,2,3)$. In the undirected cycle space over
$\mathbb F_2$, let $C_{ijk}$ denote the triangle with those vertices. Then

$$
C_{012}\mathbin\triangle C_{023}
 = \{01,12,23,03\},\qquad
C_{012}\mathbin\triangle C_{013}\mathbin\triangle C_{023}
 = C_{123}. \tag{2}
$$

Thus the square underlying (1) is a sum of two of the previous star-basis
cycles. This is a literal shared edge-set calculation, not a claim that the
tournament completions solve the colouring test.

The preceding
[odd-cycle channel packet](plane_odd_cycle_channels_20261011b.md)
used channel labels $(1,0,0,0)$. Each original star-basis triangle is
non-monochromatic, yet the residual triangle $123$ is monochromatic. Equation
(2) explains the linear cycle-space relation, while the monochromatic predicate
is nonlinear and is not preserved by symmetric difference.

On the partial square, the monochromatic graph is only the path $1$--$2$--$3$,
so it has a binary lift. Adding the missing pair $02$ does not change this:
its endpoints have different channel labels. Adding pair $13$ creates the
monochromatic triangle $123$, **in either orientation**. All four completed
tournaments have underlying graph $K_4$, so all four give exactly this same
channel obstruction.

There are therefore three separate interfaces:

| Data | Question answered | Required retained information |
|---|---|---|
| Two missing-arc orientation bits | Which of four tournament types? | Partial orientation and port frame |
| Missing undirected edges present or absent | Which graph constraints exist? | Actual graph, or exact geometric realization |
| Equality/inequality of endpoint colours | Which boundary colouring extends? | Palette partition and complete lift data |

The geometric case is treated independently in
[unit rhombus interfaces](unit_rhombus_interfaces_20261011c.md).
In particular, an abstract tournament completion to $K_4$ is not an assertion
that both diagonals of a planar unit square or rhombus are unit edges. The
orientation table supplies neither a unit-distance realization nor a plane
six/seven colouring conclusion.

The useful reusable rule is to compose through the interface actually
preserved: tournament flips through marked orientation ports, colouring
extensions through the full boundary palette table, and unit-distance gluing
through exact coordinates. Switching between these interfaces requires an
explicit map; their small state counts do not supply one.

## 5. Exact universe and reproduction

Companion:
[partial_four_tournament_20261011c.py](../../04-computation/experiments/partial_four_tournament_20261011c.py).
It enumerates every choice of four of the six pairs with all orientations:
$\binom64 2^4=240$ partial tournaments. Among them are all 48 labelled
square orientations, of which exactly 24 have the user's pattern. It checks
all completions; all 24 relabellings of every accepted marked object and
completion; all 64 complete tournaments under all permutations; both named
updates; and the explicit class-loss, degree-loss and colouring hostiles.
Boolean/float label and bit aliases are rejected.

~~~text
python -B -X utf8 04-computation/experiments/partial_four_tournament_20261011c.py
python -B -O -X utf8 04-computation/experiments/partial_four_tournament_20261011c.py
~~~

Both runs perform 5,194 exact checks. Saved output:
[partial_four_tournament_20261011c.out](partial_four_tournament_20261011c.out).
No global graph or arithmetic coverage claim is inferred from this finite
classification.
