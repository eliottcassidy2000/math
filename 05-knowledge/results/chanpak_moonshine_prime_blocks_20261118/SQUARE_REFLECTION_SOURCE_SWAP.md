# A reflection-chain source moves the window while preserving every unmarked slot count

**Status:** PROVED uniform-module replacement principle and exact finite physical construction; no universal placement, minimum-order claim, novelty of scalar counts, or merged-metagraph adjacency is asserted.

## Recovered niche and the exact bridge

THM-4510 (`01-canon/theorems/THM-4510-sum-graphs-as-reflection-orbits-pythagorean-zigzags.md`) treats sum-graph edges as reflection matchings i↔s−i. Two reflections compose to translation by their target difference. In particular the two matchings with targets n,n+1 form a Hamiltonian chain on [n]. The square-sum work uses square targets, but the reflection-matching mechanism makes sense for arbitrary targets.

Here is a declared orientation map, rather than a numerical analogy: start with i→j for i<j, and reverse exactly those arcs whose endpoint sum lies in the target set. This gives an actual tournament. The undirected reflection chain and the directed Hamilton-path count are different objects; chain existence does not prescribe that count.

At n=5 with targets {5,6}, this map gives C_R with outgoing bitmasks

[6,16,26,19,1].

Compare the source C_E from the preceding equal-gain recontact experiment:

[16,1,3,7,14].

Both have H=9 and P_2=84. More strongly, **their full unmarked slot polynomials agree**:

P_q(C_R)=P_q(C_E)=q^5+6q^3+2q.                         (1)

They are not isomorphic, even up to converse. Their start-count multisets are respectively {1,1,1,3,3} and {1,1,1,2,4}; C_E has the latter multiset for its end counts as well. An isomorphism preserves starts, and an anti-isomorphism exchanges starts and ends.

### Proof of (1)

Directly enumerating cyclic orders gives exactly three directed triangles and one directed 5-cycle in each literal tournament. Five vertices cannot contain two disjoint directed odd cycles. The inherited odd-cycle/singleton slot identity P_q(T)=G_T(2q,q) therefore gives q^5+3·2q·q^2+2q, as claimed. The replay counts the cycles independently of its Hamilton-path dynamic program. This application of the slot identity does not assert that the two tournaments have equal induced or endpoint-marked slot data.

## Why equality of all slots survives arbitrary uniform contexts

**Uniform-module replacement lemma.** If nonempty tournaments A,B satisfy P_q(A)=P_q(B) for every nonnegative integer q, replacing an A-module by B inside any fixed tournament context with uniform external contacts preserves all unmarked slot counts of the resulting tournament.

Proof. Let Q_r(A) count ordered r-tuples of nonempty paths covering A. Since empty slots can be chosen in binomial(q,r) ways,

P_q(A)=sum_r binomial(q,r) Q_r(A).

The polynomials binomial(q,r) are linearly independent, so equality of all P_q gives equality of all Q_r. In any ordered k-path cover of the ambient tournament, cut each path whenever it crosses the module boundary. Fix the outside path fragments and the ordered locations of the module's r nonempty visits. Uniform external contacts make the validity of every entrance and exit independent of the internal endpoint labels. Filling those visits contributes exactly Q_r(A). Replacing it by Q_r(B) preserves every term, for every k. Empty exterior, unused slots and zero external vertices are included.

This explains both the freedom and the loss: unmarked slot data can be preserved in every uniform context, but a contact to one particular module vertex is not uniform and need not be preserved.

## A common height-45 host

Take A=TT_2, with masks [0,1], then a singleton b, then C. Orient A→b→C; add x with x→A,C and b→x. Number A as0,1, b as2, C as3..7, x as8.

Before any return is opened, C is a terminal module, and the other four vertices have slot polynomial q^4+4q^2. Hence both choices C_E,C_R give the identical full host polynomial

q^9+10q^7+26q^5+8q^3,

whose values at1 and2 are H=45 and P_2=2688.

Now reverse a single selected cross-seam arc, allowing a return from C into A. The resulting hosts are strong tournaments, and the contact is no longer uniform. Exact counts are:

| source | arc reversed, global labels | paths removed | paths added | resulting H |
|---|---|---:|---:|---:|
| C_E | 1→6 | 4 | 90 | 131 |
| C_E | 0→3 | 2 | 90 | 133 |
| C_R | 1→6 | 3 | 99 | 141 |
| C_R | 1→3 | 1 | 99 | 143 |

Thus changing the physical source, while keeping **every unmarked slot count of the unmodified host**, moves a two-count repair window from {131,133} to {141,143}. Both centers move by exactly10.

The old common gain90 becomes99 and the removed-count pair {4,2} becomes {3,1}. The net center shift is 9−(−1)=10. This is not inferred from matching H and P_2; the literal port-dependent seam counts are retained and verified.

## Exact accounting and matrices

For any original arc u→v define

J_T(u,v)=sum_{u∈S, v∉S} E_u(T[S]) S_v(T[V\S]).

Splitting a Hamilton path at that arc gives

H(T^{uv})=H(T)−J_T(u,v)+J_T(v,u).

The latter term counts paths gained after reversing the arc. A changed arc cannot lie inside either of the two disjoint induced parts, so these partial path counts are unchanged by reversal. Equal backward gains and forward usages differing by2 therefore give the table above.

Every full adjacency matrix is supplied as outgoing bitmasks in `source_swap.json`, together with start/end vectors, P_2, removed/gained counts and exact reachability checks. The two outputs in a row pair are **different single-arc reversals of a common host**; this is not a claim that they differ from each other by one arc, or that they form a genuine metagraph line.

## What this opens, and what remains

The square-reflection source is a constructive example of a neutral source replacement that is invisible to every uniform unmarked reader but becomes visible to a selective return. It is a physically verified route around the earlier terminal-factor obstruction.

The natural next selector is a library or recursive family of *equal-slot sources with different marked seam profiles*. For target N, one needs a member and ports satisfying H(host)+backward gain−forward usage=N. The current two-source result controls a finite center shift10 and a two-count window; it does not prove that these shifts concatenate, that all centers occur, or that continuation data match after the return. Those marked conditions, rather than equality of slot polynomials, are the remaining obligations.

## Reproduction and scope

Run `python source_swap_replay.py` or `python -O source_swap_replay.py`. Both produce byte-identical JSON. Checks are exact integer subset dynamic programming, literal cycle enumeration, and graph reachability. The two reflection targets here are 5 and6, not square numbers; their relevance to the prior square-sum work is the exact shared reflection-matching mechanism.
