# Two clocks forced by nesting: an all-order reader dichotomy

The previous [physical-return session](../chanpak_physical_routes_20261115/DIRECTED_RETURNS.md)
asks for repeated positive integer shifts on an identical marked
interface. Background work observed that a stationary one-port
nest cannot do this. Here is a stronger obstruction: **even an
arbitrarily chosen or target-dependent schedule of single-port
nests along one chain, with a finite bank of readers, gives only
O(log X) distinct observed counts at most X.** The operation can
use unbounded skeleton size. No periodicity is assumed.

Throughout, S[D] replaces one marked vertex of S by a nonempty
tournament D and gives every vertex of D exactly the marked
vertex's old external arcs.

## Classify a reader by its marked vertex, not by its size

For a reader R with marked vertex v, precisely one case holds.

**No directed triangle through v.** Then v is a singleton strong
component of R. Indeed in any nonsingleton strong component,
if no out-neighbor of v beats an in-neighbor of v, the two sides
cannot return to v along a directed path. The triangle-free
marked case thus means v lies on no directed cycle. The ordered
strong components outside v concatenate their paths, and
substituting D for the singleton preserves that order. Therefore

    H(R[D]) = H(R) H(D).

**A directed triangle through v.** Write the directed triangle
v→b→c→v. Every ordered cover of D by two nonempty paths lifts
to a Hamiltonian path of the induced graph on D∪{b,c}: follow its
first path, then b,c, then its second path. Extend it by inserting
all other vertices in a fixed order using the standard tournament
insertion rule. The result is injective, since deleting those
vertices recovers the input path. This proves H(R[D])≥P₂(D).
Keeping D contiguous in any fixed Hamiltonian path of R also
gives H(R[D])≥H(D), for every reader. Thus in this case

    H(R[D]) ≥ Q₂(D)/3 ≥ 2^|D|/3,

where Q₂=P₂+2H counts optional ordered slots. Each of the 2^|D|
assignments of vertices to two slots contributes at least one
ordered cover, because every induced subtournament has a Hamiltonian
path (and an empty slot contributes one). The factor two here is
a count of assignments, not a recurrence step on the output.

## Follow any one nesting chain

Let D_(t+1)=S_t[D_t], with arbitrary marked S_t. Identity steps
leave D unchanged and other steps strictly increase its size.
If S_t is transitive then H(D_(t+1))=H(D_t): all other vertices
are just ordinal pads. Otherwise H(S_t) is odd and at least three.
Each Hamiltonian path of S_t lifts with D_t contiguous, independently
for each path of D_t. Consequently H(D_(t+1))≥3H(D_t).

For readers whose mark is in a directed triangle, any observation
at most X requires |D_t|≤log₂(3X). Only O(log X) distinct carriers
along the chain have that size, and there are only finitely many
readers. For the other readers, only the values H(D_t) matter.
The distinct such values grow by a factor of at least three, so
there are at most 1+floor(log₃ X) below X. This proves the upper
bound even if the chain visits infinitely many identities and
changes its choice of reader at every station. Quantitatively the
bound is at most the number of cyclic-mark readers times
floor(log₂(3X)), plus the number of singleton-mark readers times
1+floor(log₃ X).

The bound rules out any single infinite arithmetic progression
of positive constant shift produced by routing along one such
chain. It also rules out replacing that repeatable shift by a
fixed finite menu of bounded positive steps along the same chain.
It does **not** bound the union over all unrelated starting carriers
and all different schedules; that union can branch. Nor does it
cover moves that alter an already retained carrier, have multiple
live ports, or use an infinite menu of readers.

## What materially links the 2 and the 3 here

The two bases encode *different resources*. The 2 comes from free
assignment between two places of entry into the carrier. The 3
comes from the least odd count above one and is forced by a
nontransitive nesting step. Either growth rules out a bounded
absolute return, even when every congruence is programmable.
Contrast the prior [literal fixed-step arc routes](../chanpak_physical_routes_20261115/PORT_ROUTES.md):
they change arcs **outside** a fixed carrier rather than merely
nest that carrier again. A repeatable proof must evolve the carrier
through a different operation, or coordinate multiple live readers.
