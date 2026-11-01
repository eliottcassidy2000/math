# Two-apex continuation counts recover the actual word

**Status: PROVED by the gap bijection and exact finite identities; finite
replays are labelled separately. Prescribed-H completeness remains OPEN.**
No external literature-priority claim is made.

## 1. Inheritance and the new point

A two-apex tournament has marked vertices x->y and a transitive chain
v_1->...->v_m. A letter a_i in {0,1,2,3} has low bit 1 iff x->v_i and high
bit 1 iff y->v_i. Let H(w) count ordinary directed Hamiltonian vertex words,
not isomorphism classes, reversal classes or H/Aut.

The October20 feedback-word packet already supplies the exact gap automaton.
Session124 supplies its 36-state specialization, eight occupancy relations,
and scalar-compiler obstructions. The October30 square-sum/spectrum packet
already proves rational Hankel rank 17, the integer observation lattice, and
positive two-apex coverage through 99,999 except 7 and 21. Neither rank 17 nor
our smaller finite coverage recheck is new. In particular the 17-state
representation's lattice was expressly not its actual-word orbit.

Here the new conclusion is an explicit **digital readout from continuation
counts**, followed by an exact actual-orbit membership test. The readout
recovers both ear-neighborhood words, not only their Hamiltonian statistic.
It does not solve the existential problem of finding a state whose first
coordinate is a prescribed H.

## 2. Twelve bounded continuation probes

For a word w, write h_s=H(ws), and put

    S = (0,1,2,00,10,20,30,01,000,100,200,010).

Every probe appends at most three chain vertices. Define Q_x(w),Q_y(w) as the
ordinary binary integers whose chronological bits are respectively the low
and high bits of w, with the earliest bit most significant. Leading zeroes
are permitted. The empty word has codes (0,0).

**Theorem (exact readout).** For every finite word w,

    36 Q_x = -5h_0 -3h_1 +6h_2 +15h_00 -9h_10 -18h_20
             +6h_30 +9h_01 -4h_000 +6h_100 +6h_200 -9h_010,

    12 Q_y = -h_0 +3h_1 +3h_00 -3h_10 -3h_01 -2h_000 +3h_010.

These are exact integer identities. Their coefficients are not an
approximation or a fitted recurrence. They make no bound on |w| or H(w).

### Proof of the gap transfer used in the identity

Choose the order (a,b) in which a Hamiltonian path visits the two apices. It
has unique form P,a,Q,b,R, where each gap is an increasing subsequence and the
three gaps partition the chain vertices. For P,Q retain state 0=empty,
1=good last vertex,2=bad last vertex; good means that the last vertex points
to the apex on its right. For R retain 0=empty,1=nonempty. Q may start only
at an a-outneighbor, and R only at a b-outneighbor. Appending a chain vertex
to P, Q or R gives three additive alternatives, subject to these start guards.
This is 18 states per apex order, 36 in total.

At the end P must be empty or good. For order xy, Q may be empty or good;
for order yx, Q must be good because the marked edge is x->y, not y->x.
R has no final restriction. These conditions give a bijection with Hamiltonian
paths, including the empty chain, which has H=1.

Use column state v(w), transition matrices M_a, and output row c, so that
v(wa)=M_a v(w) and h_s=c M_s v(w), with later letters acting on the left.
The checker builds each suffix row explicitly; it does not confuse the
chronological suffix order with matrix multiplication order.

### Eight cross-order identities

Forget good/bad flags and retain only the occupancies of P,Q,R. An admissible
xy partition (P,Q,R) corresponds to the admissible yx partition (P,R,Q): the
first nonempty Q in the former must be an x-outneighbor and becomes the last
run in the latter, while the y-outneighbor condition moves in the other
direction. Terminal acceptance is deliberately not imposed on these states.

For each p,q,r in {0,1}, let E_pqr be total xy states with those occupancies
minus total yx states with occupancies (p,r,q). Thus E_pqr v(w)=0 for every
word, including empty. This is the inherited occupancy relation, not a new
symmetry of accepted Hamiltonian paths.

### Binary rows and exact certificate

Let q_x be the row summing all xy states with Q nonempty and R empty, with
no good/bad or P-occupancy restriction. Such a partial gap assignment chooses
the first Q vertex i with x->v_i. All earlier vertices go to P, and each later
vertex independently chooses P or Q. Hence its count is

    q_x v(w) = sum_i bit_x(a_i) 2^(m-i) = Q_x(w).

Likewise let q_y sum xy states with Q empty and R nonempty. Its count is
Q_y(w), by choosing the first R vertex. This also directly proves

    Q_r(wa)=2Q_r(w)+bit_r(a),     Q_r(empty)=0.

Let A_s=cM_s be each probe row. If n_x and n_y are the coefficient vectors
in the two displayed readout formulas, the following identities hold in
**all 36 coordinates**, not merely on a finite list of words:

    sum_s n_x(s) A_s -36q_x = -6(E_001+E_101)-12(E_010+E_110),
    sum_s n_y(s) A_s -12q_y = -6(E_001+E_101).

`check_word_recovery.py` reconstructs the integer rows from the stated gap
rules and verifies every coordinate of these two identities. It separately
verifies each occupancy relation's invariance under all four letters, its
zero initial value, and the binary-row recurrence. These finite algebraic
checks imply the displayed identities at **every word length** by induction.
No numerical rank decision, SymPy, discovery output, or search-success
assumption is a theorem dependency. This proves the theorem.

## 3. Word rigidity and exact orbit membership

Call w normalized if it is empty or its first letter is not 0. Delete all
leading 0s from an arbitrary word. Such chain vertices are universal sources,
including in every continuation T(ws), so deletion preserves all future H.

**Corollary 1 (rigidity).** For two words u,w, the following are equivalent:

1. H(us)=H(ws) for every continuation s;
2. H(us)=H(ws) for the twelve s in S;
3. u and w become the same word after deleting their leading 0s.

Indeed 1 implies 2. The readout gives equal code pairs under 2. Set
m=max(bit_length(Q_x),bit_length(Q_y)), pad the two binary strings to length
m, and interleave their bits as low/high bits of each letter. This recovers
the unique normalized word; when both codes vanish, it recovers empty.
Thus 2 implies 3. Universal-source deletion proves 3 implies 1.

There is no claim that the twelve probes are minimal. In particular, the
readout uses nonlinear binary decoding; it does not contradict the inherited
minimum dimension 17 for a **linear** representation with a linear output.

**Corollary 2 (exact membership algorithm).** Given any twelve integers z_s:

1. Compute the two displayed linear numerators. Reject unless they divide by
   36 and 12 respectively to nonnegative integers Q_x,Q_y.
2. Interleave their aligned binary digits to recover the unique candidate w.
3. Recompute all twelve H(ws). Accept iff they equal the supplied z_s.

This accepts exactly the actual twelve-probe profiles. On acceptance it
returns every representing word, namely 0^k w for k>=0. No enumeration of
candidate words is needed. The same algorithm decides membership in the
inherited 17-coordinate actual orbit: select the twelve coordinates, decode,
and compare **all seventeen** on replay. An actual 17-profile cannot have
another normalized representative.

If B is the maximum input integer bit length (taking B>=1), the extracted
codes have O(B) bits. The candidate has O(B) letters; the fixed-size gap
transfer requires O(B) arithmetic updates on O(B)-bit integers. Hence the
algorithm runs in O(B^2) elementary bit operations and O(B) space using
ordinary integer addition and constant-denominator arithmetic. This claim is
about fixed-dimensional state membership, not enumeration by target H.

The actual normalized orbit is therefore the nonnegative code-pair plane:
append letter a sends (Q_x,Q_y) to (2Q_x+bit_x(a),2Q_y+bit_y(a)). Its inverse
address is the pair of parities followed by componentwise halving. Every
nonroot state has a unique last letter and predecessor; the empty state's
0-loop is precisely source padding. The state-to-word inverse terminates.

## 4. Hostiles and the surviving scalar problem

Nonnegative integral codes are necessary but not sufficient without replay.
For example, 3 times the empty word's entire 17-profile is positive integral
and lies in the inherited observation lattice. Both extracted codes are 0,
so it decodes to empty, but its H coordinate is 3, not 1. Replay rejects it,
even though 3 itself is a genuine tournament count. The same test rejects
7 and 21 times the empty profile. These latter counts are globally forbidden
by inherited theorems, but no such theorem is needed for this rejection.

Thus the old lattice relaxation remains false; it now has a direct exact
membership repair. The target problem is still different:

    Given allowed odd h, find Q_x,Q_y>=0 whose decoded word has H=h.

Knowing h alone does not supply the twelve probes or the two codes. This
packet does not claim a decreasing scalar-H inverse, an interval induction,
or all odd values other than 7 and 21. The two-apex completeness strengthening
and the unrestricted tournament-spectrum conjecture both remain OPEN.

## 5. A finite arithmetic box for a scalar target

This is a separate **CITED-bound corollary**, not the missing completeness
argument. Busch proves H(C)^3 >= 5^(|C|-1) for every nontrivial strong
tournament C; strong-component counts multiply. See
[Busch, Theorem 1 and Corollary 1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf).

**Corollary 3.** If h is attained by a two-apex word, it is attained by one
whose spine length is at most

    L(h) = max {m>=0: 5^m <= h^3}.

Delete every chain vertex whose strong component is a singleton. This
preserves H, the two marked apices, and transitivity of the remaining chain.
Every remaining chain vertex belongs to a nontrivial strong component C_i;
each such component contains at least one apex, since the chain alone is
acyclic. If n_i is its order and a_i>=1 its number of apices, the remaining
spine length is

    m = sum_i(n_i-a_i) <= sum_i(n_i-1) <= log_5(h^3).

Singleton apices contribute no chain vertex and no Hamiltonian factor, so
no extra additive constant is needed. The empty sum handles h=1. The cyclic-core inequality is attained: w=100311 is the ordered sum of two
strong four-vertex one-apex tournaments, each with H=5. Its six chain
vertices are all cyclic, and L(25)=6. This is not a lower bound for the
minimum order of an H=25 witness; the shorter word1120 also has H=25.

Equivalently, h lies in the two-apex spectrum iff some pair

    0 <= Q_x,Q_y < 2^L(h)

has decoded count h. Thus there is a complete, terminating **all-order
membership decision for this restricted spectrum**, using at most
4^L(h) <= h^(log_5 64) candidates and polynomial-in-log(h) work per candidate.
This is pseudopolynomial in numeric h, not polynomial in its bit length.
`find_count` implements a depth-first version with the inherited Rédei oddness guard, monotone-H pruning and
leading-source normalization; `None` is an exclusion only from the two-apex
language. Its totality relies on the cited bound and the deletion argument,
not an arbitrarily selected maximum depth. The readout theorem supplies the
actual-word coordinates but is not needed to rediscover the known SCC bound.

The cutoff preserves scalar H, not the entire continuation profile. For
example H(3333)=1 but H(33330)=47, whereas H(empty)=H(0)=1. Deleting its
singleton chain vertices is legitimate for scalar search and illegitimate
for replacing an actual response state by the empty profile.

The box still grows with h. It does not prove the box contains a solution
for every allowed odd h, nor refute the unrestricted spectrum when it does
not. The expensive search is optional; the default verifier uses only small
explicit targets.

## 6. Finite verification, separate from the induction

Normal Python and `python -O` give identical receipts. The checker tests all
21,845 words of length at most seven and 24 seeded random words through
length 1,000, for 21,869 actual-word recovery tests. Independent adjacency
subset DP agrees with H on all 341 words through length four; all 341 also
pass singleton-chain-SCC deletion and the spine cutoff. Fifteen hostile
or malformed-input tests are rejected. Separate session audit reconstructs
the probe values by adjacency DP without importing this producer (4,437
literal continuation counts). It also independently proves the all-order
readout using a differently rebuilt transfer and a closed 28-dimensional
reachable span, rather than importing this packet's occupancy certificate.
The manuscript and its complexity/scope claims passed that independent
review.

An exploratory scan also recounted 19,876 distinct witnesses through length
ten and recovered a smaller interval through 24,105 except 7 and 21. This is
only a recheck, weaker than inherited coverage through 99,999, and is NOT a
new interval theorem or a dependency of the decoder proof.
