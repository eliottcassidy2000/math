# Exact dynamics of interacting regular tournament modules

**PROVED:** the module splitting and the all-word similarity classification below.
**FINITE-EXACT:** the recorded matrix controls and the two order 21 Hamiltonian counts.
**OPEN:** Chan–Pak, full tournament classification, and the proposed all-word decoder.

## 1. Inheritance and the new object

The inherited `THM-4857-the-gf2-automaton-of-a-tournament-classifies-the-merged-metagraph-exactly-through-order-five-circulants-are-additive-cellular-automata-and-paley-tournaments-give-quadratic-residue-codes.md`
uses the actual adjacency matrix over F2: a vertex updates to the parity of the
bits at vertices it beats. It proves that similarity loses tournament classes
from order 6, and recalls the Paley7 idempotent of rank4. We use those scoped
inputs without re-auditing the old census or making a new literature claim.

Session267's `topology/boundary_recursion_bridge.md` supplies an explicit
Paley7/C3 wreath-action model. Here the same physical substitution is studied
as a linear dynamical system, including **both** orders of the factors. The
new connection is the exact internal/quotient splitting, not a numerical
association between7,21 and periods.

Throughout Q[B] means replacing each vertex of Q by the same tournament B,
with all arcs between modules directed by Q. A tournament is regular if all
outdegrees equal (n−1)/2; its order is then odd. Write A_T for its F2 adjacency
matrix, 1 for the all-one vector, and W_T={x:1ᵀx=0}.

## 2. Splitting theorem

For regular odd-order B of order m and degree d_B modulo2,

    F2^(qm) = (F2^q tensor <1_B>) direct-sum (F2^q tensor W_B),
    A_(Q[B]) ~ (A_Q+d_B I_q) direct-sum (A_B|W_B)^[q].       (1)

The outer tournament Q need not be regular in (1). The superscript means q
copies, not a matrix power.

**Proof.** Since m is odd, <1_B> intersects W_B trivially and they span the
whole child space. Regularity gives A_B1=d_B1 and 1ᵀA_B=d_B1ᵀ, so both
subspaces are invariant. The actual adjacency matrix is

    A_(Q[B]) = I_q tensor A_B + A_Q tensor J_m.

J_m acts as identity on <1_B> and as zero on W_B. Restriction gives (1).
This is a direct-sum conjugacy, stronger than equality of characteristic
polynomials. It retains transient and periodic structure.

For regular factors T_1,...,T_l of odd orders n_i and degree parities d_i,
the iterated lexicographic product has the following complete block split:

- one global scalar eigenvalue d_1+...+d_l;
- for each i, the block A_(T_i)|W_(T_i)+(d_(i+1)+...+d_l)I,
  repeated n_1...n_(i−1) times.

This follows by induction from (1). The degree parity of an odd-order product
is the sum of factor degree parities. Internal modes are copied; deeper
global parity shifts the ancestor's internal dynamics. These are literal
interactions of state spaces, not independent copies without coupling.

## 3. A complete infinite two-letter dynamics classification

Take C=C3 and P=Paley7 with arcs i→i+s for s in {1,2,4} modulo7.
For a word w in {3,7}^l, let T_w be the lexicographic product of those factors
in the written outer-to-inner order. Both factor degrees are odd.

On W_C, A_C has the irreducible polynomial p(t)=t²+t+1, one two-dimensional
block. On W_P, the inherited idempotent has three zero eigendirections and
three one eigendirections. Both structures are unchanged up to similarity
by addition of I: p(t+1)=p(t), and the equally sized P eigenspaces swap.

Put

    N=product_i n_i,
    a(w)=sum_(i:n_i=3) product_(j<i) n_j,
    b(w)=3 sum_(i:n_i=7) product_(j<i) n_j,
    epsilon=l mod2.

Then the elementary-divisor multiplicities are exactly

    t: b+1−epsilon,
    t+1: b+epsilon,
    t²+t+1: a.                                      (2)

All blocks are semisimple, and N=1+2a+2b. Thus (2) classifies the F2 automaton
for **every word**, without constructing its N×N adjacency matrix. Two words
produce similar automata iff the three displayed multiplicities agree.
Computing them takes O(l) integer operations, versus storing N² bits.
This is classification of this infinite grammar's automata, not a proof that
similar automata identify their tournaments.

Let m0=b+1−epsilon, m1=b+epsilon. Every trajectory enters the periodic subspace
after at most one step. There are 2^m1 fixed points and

    (2^(m1+2a) − 2^m1)/3

cycles of length3. Each periodic point has exactly 2^m0−1 additional immediate
preimages in the transient layer. In particular every matrix obeys A^4=A.
The theorem includes the empty word (one vertex, A=0).

**Order matters.** C[P] on 21 vertices has (m0,m1,a)=(10,9,1), while P[C]
has (4,3,7). Their adjacency ranks are11 and17. Thus the ordinary size and
the multiset of module types are insufficient, even though all eventual
periods divide3. The same concrete graphs have, by the exact run-word
contraction implemented in `replay.py`,

    H(C[P]) = 122621733302733,
    H(P[C]) = 123156134945235.

The small-seed profiles are computed by subset dynamic programming; the
order 21 calculation uses512 and16384 quotient occurrence boxes respectively,
not the2^21 subset table of either finished graph. This is an exact
constructive example, not a claim that H is determined by (2).

## 4. An arbitrarily large family with only a three-step clock

For w=3^d, N=3^d and a=(N−1)/2, b=0. Hence

    A^3=I          when d is odd,
    A^3=I+J        when d is even.

For odd d the automaton is reversible, has 2 fixed points, and
(2^N−2)/3 three-cycles. For even d it kills exactly the global bit, has only
the zero fixed point, and has (2^(N−1)−1)/3 periodic three-cycles. There is
one transient preimage attached to every periodic point.

So a recursively growing tournament can have exponentially many vertices,
yet an elementary parity-alternating clock description. The inherited
THM-3213 non-P-recursiveness concerns a Hamiltonian counting sequence in its
specified balanced blowup regime; it is a different observable. A short
finite-field dynamics description never supplies an integer path-count
recurrence for free.

## 5. The connection contract and next question

- Source: a regular odd lexicographic module word with physical adjacency.
- Target: its parity-update linear automaton and exact block decomposition.
- Preserved: all transient/periodic trajectories up to linear conjugacy.
- Lost: individual vertex names, permutation conjugacy, Hamiltonian paths,
  and the full rooted or path-cover response used in the other lanes.
- Sidecar: the canonical module word and the appropriate response profile.
- Decisive controls: actual matrices through C3 depth5 (243 vertices), nine
  mixed words through 63 vertices, and two exactly evaluated order 21 graphs.

An inexpensive question survives: within the two-letter grammar, do N and
(a,b,epsilon) recover the factor word? Exhaustion through word length14
finds no collision, but this is recorded only as a **finite experiment**.
The theorem is (2), which is valid with or without an injective decoder.
Do not promote the experiment to a general isomorphism classification.

Run `python3 replay.py`; the emitted JSON is the committed `receipt.json`.
The checks use explicit exceptions and remain active under `python3 -O`.
