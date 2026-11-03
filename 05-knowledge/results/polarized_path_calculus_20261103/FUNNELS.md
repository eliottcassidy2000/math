# A complete unique-start grammar and its word compiler

**PROVED.** All tournaments are finite. H counts directed Hamilton vertex
words, H(empty)=1; s_T(r) counts those starting at r. A singleton is allowed
as a terminal core. A nontrivial strong tournament has at least three paths:
a directed Hamilton cycle supplies one from each vertex (the standard
Camion theorem). In particular every vertex of a strong core has s>=1.

## F1. Deleting the root gives an exact start-budget law

Let the strong components of T-r be U_1 -> ... -> U_k. Hamilton paths in
this ordinal sum traverse each component once in that order. Therefore

    s_T(r) = [product_(i=2)^k H(U_i)]
             sum_(v in U_1, r->v) s_(U_1)(v).              (1)

This is a bijection, including the edge from r to the first component.
Suppose T is strong and s_T(r) is one or two. The sum in (1) is positive.
Every U_i for i>=2 must be a singleton, since a nontrivial strong component
would multiply it by at least three. Write U=U_1 and let the remaining
vertices form a transitive tail v_1 -> ... -> v_t, dominated by U.
Let J be the out-neighbors of r inside U. Then

    sum_(v in J) s_U(v) = s_T(r).

For budget one, J={w} and s_U(w)=1. For budget two there are exactly two
possibilities: J={w} with s_U(w)=2, or J={w,z} with both starts equal one.
If t>0, strongness forces v_t -> r. Conversely these data produce a strong
tournament: r reaches U; U reaches every tail vertex; each tail vertex can
reach v_t and then r. For t=0, strongness holds when U is nontrivial strong
and J is a nonempty proper subset. The singleton t=0 output is transitive
and is excluded. These are the complete strongness guards.

### A constructive recognition algorithm

Equation (1) gives an O(n^3) algorithm, using O(n^2) strong-component
decomposition at each of at most n recursive levels. For any rooted
tournament, delete r; reject if a component after the first is nontrivial,
or if r has other than one out-neighbor w in the first component. Recurse
on that component at w. The singleton accepts. On acceptance concatenate
r, the returned core path, and the ordered singleton tail. This returns
the unique path exactly when s_T(r)=1; otherwise it rejects. It never
counts all Hamilton paths. The unique vertex sequence also gives a
canonical rooted-isomorphism labelling in this class. `recognize.py`
implements the inverse via Kosaraju decomposition. The algorithmic bound
uses an adjacency matrix; the reference code has Python's usual recursion
depth limit, which can be removed by an explicit stack for large inputs.

## F2. Fibonacci enumeration from the unique decomposition

Let u_n count rooted isomorphism classes of strong n-vertex tournaments
with s_T(r)=1, for n>=3. Equation (1) is a unique inverse:

* U is a singleton, t=n-2>=1: the last tail bit is forced and the others
  are free, giving 2^(n-3) possibilities.
* U is a smaller unique-start strong rooted core, of order m>=3. A zero
  tail costs one new root. A tail of length t>=1 costs t+1 vertices and
  contributes 2^(t-1) words.

No isomorphism overcount occurs. An isomorphism preserving r preserves the
ordered components of T-r, the unique entry w, and every tail vertex.
Conversely an isomorphism of these data extends to T. A rooted automorphism
also fixes the unique Hamilton path starting at r pointwise, so is trivial.
Consequently, as a formal power series,

    U(x) = x^3/(1-2x) + [x+x^2/(1-2x)] U(x)
         = x^3/(1-3x+x^2).

With Fibonacci numbers F_0=0,F_1=1, this gives u_n=F_(2n-4).
It is an all-size isomorphism-class enumeration, not an empirical fit.

For arbitrary, possibly nonstrong, rooted tournaments with s_T(r)=1,
the first component is either the rooted singleton or one of these strong
cores, followed by a transitive tail. Their class generating function is

    A(x) = (x+U(x))/(1-x) = x(1-2x)/(1-3x+x^2).

The counts are 1,1,2,5,13,34,89,... at orders 1,2,... . For any n, the number
of labelled rooted pairs (T,r) in either class is n! times its class count;
with a specified root label it is (n-1)! times that count.

This gives a precise weighted slice of tournament enumeration. If u(T)
is the number of unique-start vertices in an unrooted strong T, then

    u_n = sum_[T strong, |T|=n] u(T)/|Aut(T)|.

The automorphism action on these vertices is free. This is not A000568
itself, nor a new extension of its computed sequence.

### A smallest faithful fixed-order class carrier

Fix n>=3. Order the possible first core sizes as 1,3,4,...,n-1. For core
size m put t=n-m-1 and c_t=1 if t=0, otherwise c_t=2^(t-1). The branch
contains c_t u_m objects, using u_1=1 only for the terminal singleton.
Order its free tail words by their binary value b, then its core classes
by their recursive rank j. Its rank is

    offset_(n,m) + b u_m + j,                         (R)

where the offset sums all earlier branch masses. The unique decomposition
proves by induction that (R) bijects the classes with 0,...,u_n-1. Quotient
and remainder by u_m recover b and j. A decoder reconstructs every arc,
then names vertices in the unique path order. This gives a faithful
carrier with exactly u_n states, the smallest possible for distinguishing
these rooted classes at this fixed order.

With n retained as a separate size clock, a binary code needs exactly
ceil(log2 u_n) fixed bits, asymptotic to n log2((3+sqrt(5))/2)+O(1).
For example n=100 requires 135 bits. This bound is attained by the rank,
and its lower bound is simply injectivity into binary strings.

The carrier is closed under every strong one-port funnel: (R) updates the
rank directly from the old rank, tail length and word. The canonical old
vertex at position i moves to i+1, the new root is zero, and the tail
follows the core in its transitive order. Thus the construction also gives
explicit name transport. It is not a minimum carrier for the full
metagraph; arbitrary re-marking or complements can leave this family.
The polynomial carrier computes path counts; the rank carrier retains the
actual class. They solve different observation tasks.

## F3. Compile every binary tail, not a list of examples

Start with any rooted nonempty U at w. Adjoin a tail v_1,...,v_t, with
U -> tail, and a new root r. Set r->w, U-w->r. Write a binary word
epsilon_1...epsilon_t, where epsilon_i=1 means v_i->r and zero means r->v_i.
Let P,E,S be the path-bin polynomials defined in CALCULUS. Put

    L(m) = (m-1)^t + sum_(i:epsilon_i=1) m^(i-1)(m-1)^(t-i),
    D(m) = (m-1)^t + sum_(i:epsilon_i=0) (m-1)^(i-1)m^(t-i).

L counts assignments to m bins in which the last vertex of the first bin
beats r, or that bin is empty. D imposes the analogous condition that the
first vertex of the distinguished bin loses to r, or the bin is empty.

Let B(m) count assignments of the tail to two distinguished bins, before
and after r, and m-1 ordinary bins. The last before vertex must have bit1;
the first after vertex must have bit0. Either bin may be empty.

Here is a six-state weighted automaton for B. Its state is (b,a), where
b in {empty,good,bad} records whether the last before vertex is absent,
has bit1, or has bit0; a in {empty,started} records the after bin. Start
at (empty,empty). At each bit:

* an ordinary-bin choice keeps the state, with weight m-1;
* a before-bin choice replaces b by good or bad according to that bit;
* an after-bin choice is permitted if already started or the bit is zero,
  and sets a=started, with weight one.

Accept exactly b!=bad. Thus B is an integer polynomial and each query B(m)
costs O(t) arithmetic operations with six states. Six is a sufficient
state count; no minimality is asserted. Empty-word L=D=B=1.

The exact output polynomials, rooted at r, are

    P_new(m) = m[L(m+1) S(m+1) + B(m) P(m) - D(m) E(m)],
    E_new(m) = L(m) P(m) - (m-1)^t E(m),
    S_new(m) = m^t S(m) + D(m) P(m-1).                  (2)

Proof: choose the bin containing r and split it there. If some U vertices
follow r, their path must start at w, and the U-before path cannot end at w
because w is already used. This contributes L(m+1)S(m+1). If no U vertex
follows r, all core paths use m slots; the two tail parts contribute B(m).
The invalid case of an empty before-tail and a U path ending at w is
exactly D(m)E(m). Multiply by the m choices of root bin. Requiring r at
the end or start of the first bin gives the other two equations by the
same disjoint decomposition. This proves (2) for every tail word and core.

At m=1, write G=S(2), h=P(1), e=E(1). Then

    H_new = [1+sum_(i:epsilon_i=1) 2^(i-1)]G + B(1)h - D(1)e,
    s_new = s_U(w).

This preserves a unique start through arbitrary depth. If t>0 and the last
bit is one, e_new=h. For t=0, e_new=h-e. The full polynomials make the
operation recursively reusable, including mixed tail lengths and words.

## F4. The exact two-port correction

Permit r to beat a set J in U instead of just w. Define E_J=sum_(w in J) E_w
and S_J=sum_(w in J) S_w. Let C_J(m+1) count ordered m+1-bin covers of U
whose first bin ends in J and whose second bin starts in J. Both bins are
nonempty; the other m-1 may be empty. Their two endpoints are necessarily
distinct because the paths have disjoint supports.

The output laws become

    P_new(m)=m[L(m+1)S_J(m+1) - m^t C_J(m+1)
                +B(m)P(m)-D(m)E_J(m)],
    E_new(m)=L(m)P(m)-(m-1)^t E_J(m),
    S_new(m)=m^t S_J(m)+D(m)P(m-1).                    (3)

The extra subtraction removes exactly the paths with nonempty U portions
on both sides of r, no tail vertex before r, and an invalid predecessor
in J. After removing that forbidden before bin, tail vertices have m
unrestricted choices. For one port C_J=0, recovering (2).

Thus the two-start classification in F1 leads to an explicit joint
incidence coordinate. It does not justify summing two marginal one-port
controllers. A complete recursive calculus retaining two old ports would
need their further joint transport; that is the next construction target.
