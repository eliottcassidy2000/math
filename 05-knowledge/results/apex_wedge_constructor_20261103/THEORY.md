# Constructing arbitrary even Hamiltonian edge jumps

**PROVED scoped.** A1--A4 below concern actual labelled tournaments, not a
weighted relaxation. They prescribe the change in the Hamiltonian path
count under one edge reversal. They do not prescribe either absolute count.

## A1 The contracted edge kernel

Let the base be the transitive tournament on 0,...,m-1, with i pointing to j
when i<j. Add named vertices a=m and b=m+1. For binary words w,u of length m,
let a point to i iff w_i=1 and b point to i iff u_i=1. Write T+ for a->b
and T- for b->a. H counts directed Hamiltonian vertex sequences.

If both words start with 1 and end with 0, both tournaments are strong:
each apex and the base contain the directed cycle
apex->0->...->m-1->apex. This is a sufficient guard used by our constructors.

For m>=1 define
\[
 K_{ij}=\begin{cases}
 1&j=i+1,\\
 2^{i-j-1}&i>j,\\
 0&\text{otherwise}.
 \end{cases}
\]
The number of paths in T+ using the consecutive block a,b is
\[
 F(w,u)=1-w_{m-1}+u_0+(1-w)^T K u.                 \tag{A1}
\]
For the empty base set F=1.

**Proof.** A path with this block partitions the base into an increasing
subsequence A before it and an increasing subsequence B after it. Empty B
contributes 1-w_last, and empty A contributes u_first. Otherwise put
i=max A,j=min B. If i<j, the partition is possible only for j=i+1 and is
unique. If i>j, all vertices strictly between j and i choose either side
independently, giving 2^(i-j-1). The two required boundary edges have
indicator (1-w_i)u_j. This proves the formula.

The other Hamilton paths use no a-b adjacency and remain unchanged by
reversal. Consequently
\[
 D(w,u):=H(T+)-H(T-)=F(w,u)-F(u,w).                \tag{A2}
\]
This is an exact derivative of the count with respect to a physical edge.

## A2 Balanced binary words

For positive integers P,R, write p for the ordinary binary expansion of P
and s for the reverse of the bitwise complement of the binary expansion
of R. Concatenation gives
\[
 w=p\,01\,s,\qquad u=p\,10\,s.
\]
These words satisfy the strong guard. Their jump is
\[
 D(w,u)=2\bigl(\lceil R/2\rceil-\lceil P/2\rceil\bigr).       \tag{A3}
\]

**Proof.** Put k=length(p), t=length(s). The integral alternating expansion
of A2 is proved in G1 of [GEOMETRY.md](GEOMETRY.md). Since u-w has entries
+1,-1 in positions k,k+1, its linear contribution to D/2 is 2^(t-1).
The prefix contribution is
\[
 -p_{k-1}-\sum_{i=0}^{k-2}p_i2^{k-i-2}=-\lceil P/2\rceil.
\]
The suffix contribution is
\[
 -s_0-\sum_{j=1}^{t-1}s_j2^{j-1}
 =\lceil R/2\rceil-2^{t-1}.
\]
These identities include one-bit prefixes and suffixes. Adding them proves
A3. Both even and odd P,R and unequal binary lengths are allowed.

This cancellation is why two apexes help: the earlier one-apex adjacent-swap
formula has only one of these charges. See
[feedback-word proofs](../feedback_word_hamilton_spectrum_20261020/proofs.md).

## A3 Every even jump has a logarithmic construction

For any integer d>=1 choose P=1,R=2d+1 in A2. Then T+ and T- are strong,
differ only in the named edge a-b, and
\[
 H(T+)-H(T-)=2d,\qquad |V|=\operatorname{bitlength}(d)+6.
\]
Indeed bitlength(2d+1)=bitlength(d)+1. Negative jumps follow by exchanging
the two endpoints of the move; zero follows from w=u with the strong guard.
Thus the edge-jump spectrum among strong tournaments is exactly the even
integers. The necessity of evenness here follows directly from A2/G1;
for arbitrary tournaments it also follows from the inherited oddness of H.

For example d=1 gives counts 133 and 131 on seven vertices. The construction
also works without enumeration for d=10^20+1. It specifies actual words and
all edges, not only an integer observation vector.

The order bound is a sufficient explicit bound; optimal order is not claimed.
This theorem does not turn a tournament with prescribed count h into one
with count h+2d. Its starting count is part of its constructed source.

## A4 Two native moves preserve the named apex ports

A native path complement K_P reverses every edge outside a supplied directed
Hamiltonian path P and keeps every edge of P. This is the native operation
used in the earlier [controller packet](../ternary_repair_clock_20261103/native.py);
the generic one-helper arc compiler is inherited.

For every pair of nonempty words with w_0=u_0=1, there is the following
uniform explicit specialization. Start at T+ on n=m+2 vertices. Add t=n,
which beats only a and loses to every other old vertex. Use
\[
 P=(a,b,0,\ldots,m-1,t),\qquad
 Q=(a,t,b,0,\ldots,m-1).
\]
P is Hamiltonian in the enlarged source. In K_P of that source, the edges
a->t and t->b were reversed from edges outside P, and the rest of Q lies
in P, so Q is Hamiltonian too. Apply K_Q and delete t.

The endpoint is exactly T- with every original name unchanged. To see this,
two complements cancel on edges belonging to both path supports or to
neither. The supports differ only at a-b, (m-1)-t, a-t, and t-b.
Only a-b has two old endpoints.

The vertex clock is [n,n+1,n+1,n+1,n]: birth, first move, second move,
deletion. The cost is exactly two native moves, one temporary birth and one
deletion. Both physical ports a,b survive. Intermediate H and strongness are
not asserted to remain constant. No reentry into a restricted birth alphabet
is inferred merely from this program.

## Exact finite controls

[verify.py](verify.py) compares subset DP, the inherited gap-partition
mechanism, and A1 on all 2,730 graphs from word lengths 0 through 5.
Permutation enumeration supplies 170 independent small graph counts.
It also replays the native program with a separate edge-matrix implementation,
tests all balanced P,R<=32, and supplies literal controls for the infinite
families in [RECURRENCES.md](RECURRENCES.md). These checks support the
all-order proofs above; their finite bounds are not theorem quantifiers.
