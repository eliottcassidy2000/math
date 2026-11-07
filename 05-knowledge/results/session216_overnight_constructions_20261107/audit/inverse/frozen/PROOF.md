# Exact inverse membership for nonuniform apex-mask words

Status: PROVED specialization of positive weighted-automaton saturation;
FINITE-EXACT construction certificate. No eventual-spectrum theorem.

## Scope and inheritance

Fix a labeled tournament A on k apex vertices. A word w over {0,...,2^k-1}
builds a tournament T_A(w) by adjoining a transitive chain in word order.
Bit i of a letter is one exactly when apex i beats that chain vertex.
The transfer below is inherited from Session215's multi-exit packet (and its
credited earlier run transfers). Saturated finite membership is the general
technique already used in Session203/saturated-profile-membership. The present
contribution is its exact nonuniform mask-language specialization, an executable
inverse returning literal graphs, and a small-order certificate for the full
unrestricted two-apex word class. No generic automata or saturation priority
claim is made.

## Transfer and finite inverse theorem

For each permutation p of the k apexes, split a proposed Hamiltonian path into
k+1 gaps before, between, and after them. Chain vertices assigned to one gap
must occur in increasing chain order. Each gap has state 0 (empty), 1
(nonempty, latest vertex fails its successor boundary), or 2 (nonempty,
latest vertex meets its successor boundary). The last gap never has state 1.

Process the word from left to right. Assign its next vertex to one gap g.
An empty gap can be entered only if g=0 or its predecessor apex beats this
vertex. A nonempty gap can always be extended. Its new state is 2 exactly
when g=k or the vertex beats its successor apex; otherwise it is 1. Other
gap states stay unchanged. Every allowed assignment contributes one. The
initial vector has weight one at every permutation's all-empty state.
At the end accept states in which every nonempty gap has state 2 and every
empty interior gap has the required predecessor-to-successor apex arc.
Boundary empty gaps impose no condition. This is a bijection with Hamiltonian
paths: read the apex permutation and maximal intervening chain subsequences
from the path, and conversely concatenate them. Thus all transfer entries and
final coefficients are nonnegative integers.

Let d=k!*2*3^k be the number of potentially nonzero coordinates. Our simple
implementation retains k!*3^(k+1) coordinates, some permanently zero. For a
requested nonnegative target h, replace every coordinate x by min(x,h+1).
The set {0,...,h+1}, with capped sum and product, is a semiring, and this
capping map is a homomorphism from the nonnegative integers. Inductively,
every stored coordinate is precisely the cap of its genuine coordinate.
The capped final sum equals h if and only if the genuine count equals h.
A large hidden coordinate can later be discarded by a zero coefficient;
this is harmless. It is NEVER used by itself to reject a state.

Breadth-first search over these saturated states therefore terminates after
at most (h+2)^d states. If an accepting state is reached, parent pointers
recover a shortest mask word. If all reachable states are exhausted, h is
absent from this entire fixed-apex class at every chain length. An accepted
word has length at most (h+2)^d-1, by deletion of repeated states on a shortest
accepting run. This bound is intentionally enormous, not a practical
all-target algorithm. Elementary implementation complexity is O(2^k*(k+1)*d*
(h+2)^d) arithmetic steps and O(d*(h+2)^d) integer storage, with integers
bounded by h+1 (plus a temporary sum).

Optional exact pruning: if the current output already exceeds h, discard its
extensions. Hamiltonian-path counts are monotone under adding a vertex.
For each old Hamiltonian path, insert the new vertex at the start if it
beats the first vertex, at the end if the last beats it, or at an adjacent
transition old->new->old otherwise. Fix a deterministic choice. Removing
the inserted vertex recovers the old path, so these insertions are injective.
Iterate over new chain vertices. Testing output>h on capped coordinates is
exact: a summand exceeding h remains at least h+1, and otherwise every
summand is unchanged. This pruning never assumes every hidden coordinate
is bounded by h. Finite resource interruption reports INCOMPLETE, not ABSENT.

## Exact two-apex results

Set k=2 and A=(0->1). `inverse.py` exhausts 708 saturated states for h=7 and
4626 for h=21, giving all-order class absences. These are hostile controls,
not new ambient exclusions. It finds positive witnesses for 1,3,5,9,23.

`census.py` enumerates every word of length at most 10 with the justified
current-output pruning above and saves 4,998 literal witnesses. They realize
every positive odd h<=10,000 except 7 and 21 on at most 12 vertices. This is
not an extension of Session203's unrestricted 32-bit coverage; it is a
certificate inside the much narrower chain-plus-two-apex class. `audit.cpp`
checks the saved graphs independently using subset Hamiltonian-path DP.

The exit-budget probe is finite only. Among words of length<=10, three exits
at the cheaper apex suffice for every certified count. Exactly 8725,9153,
9643,9739 need more than two exits at both apexes in this finite universe.
Longer words are not excluded. The term exit here means an apex beating a
chain vertex, not an arbitrary boundary state.

## Digit probes and the remaining gap

`digits.py` tests constant suffix macros of length<=3, all seeds of length
4..6 having H>=23, and the two desired branches H->2H-1 and H->2H+1.
18 seed/macro-pair candidates satisfy both first steps; none survives all
second steps. `digits_boundary.py` permits prefix/suffix macros of length<=2
and duality after each action, with seeds of length4..5. Its96 first-step
pairs likewise all fail by depth2. These are explicitly bounded hostile
probes, not obstructions to larger macros, varying finite control, rebuilding,
or higher apex number. In particular, they do not contradict the finite
coverage certificate or imply a missing ambient count.

The inverse theorem decides each target but does not say every target is
accepted. The universal all-odd assertion, a uniform induction grammar, and
any eventual threshold remain OPEN.
