# A bounded low-start interface for strong tournaments

**PROVED, independently audited.** The endpoint inequality, bounded-entry
decomposition, and actual-entry corollary have passed independent review.
This is an all-order structural reduction, not a proof of F(8,137), a finite
arithmetic automaton, or an eventual odd-spectrum theorem. No priority claim
is made for the endpoint inequality.

## 1. A quantitative endpoint lemma

Write s_x(T) for the number of directed Hamilton paths in T starting at x,
and H(T)=sum_x s_x(T). A singleton is strong and has H=s=1. Every strong
tournament of order at least three has a Hamilton cycle (Camion's theorem);
consequently every vertex starts a Hamilton path and H(T)>=|T|. This is the
same standard strong-tournament fact used by the independently audited
session203/rooted-profile-order-bound proof.

**Theorem 1.** If T is a nonempty strong tournament and S is a nonempty subset
of m vertices, then

    sum_{x in S} s_x(T) >= floor(m*m/4)+1.                 (1)

Consequently, for each integer k>=1,

    |{x in V(T): s_x(T)<=k}| <= 4k-1.                    (2)

In particular a strong tournament has at most 103 vertices with start count
at most26. This is a uniform bound on the eligible entry vertices, not a
bound on the order of the whole tournament.

### 1.1 One mixed-vertex extension

Let U be strong and add z with nonempty proper out-neighbor set A in U.
The child U+z is strong. Put p=s_z(U+z). Deleting the first vertex of a
root-z Hamilton path gives the exact identity

    p=sum_{a in A} s_a(U) >= |A|.                        (3)

Every old root-start path can be extended while preserving its first
vertex: if its last vertex lies outside A, append z; otherwise insert z
at its last transition from U-A into A. Such a transition exists unless
the path starts in A and never leaves A, which is impossible because U-A
is nonempty. Thus every old start count is nondecreasing.

There is a useful extra path when a Hamilton path begins and ends outside
A. Inserting z just before its first A-vertex and appending z produce two
distinct child paths. Different old paths remain distinct after deleting z.
For each x outside A,

    s_x(U+z)-s_x(U) >= number of old paths from x
                               whose last vertex is outside A.      (4)

Fix one Hamilton cycle of U. Cutting it just before x gives a path from x
to its predecessor on that cycle. Among any marked set S0 of m old
vertices, at most |A| have x in A and at most |A| have their cycle
predecessor in A. Therefore at least max(0,m-2|A|) marked vertices each gain
at least one path. If the new vertex is also marked, the total marked
start sum increases by at least

    p+max(0,m-2|A|) >= p+max(0,m-2p) >= ceil(m/2).        (5)

The last inequality uses only that p and m are integers. If the new vertex
is not marked, the marked sum is simply nondecreasing.

### 1.2 A strong growth sequence requiring no deletion theorem

Every strong T of order at least three contains a directed triangle.
Starting from such an induced triangle U, grow to T as follows.

If some outside vertex has both an in-neighbor and an out-neighbor in U,
add it, using the preceding mixed-vertex step.

Otherwise every outside vertex either dominates all of U or loses to all
of U. Call these sets A and D respectively. Strongness of T forces both
sets to be nonempty and an arc d->a with d in D,a in A: if either set is
empty or all arcs go A->D, a nontrivial ordered cut contradicts strongness.
Add this pair. Its induced orientation is

    U -> d -> a -> U.

The new tournament is strong. Every old path extends by appending d,a,
so all old start counts are nondecreasing. Moreover

    s_d(U+d+a)=s_a(U+d+a)=H(U)>=|U|.                     (6)

These identities follow from the forced root-d order d,a,U and root-a
order a,U,d. The procedure increases order and terminates at T.

### 1.3 Proof of (1)

Follow this growth sequence and mark precisely the vertices belonging to
the final S. Let m denote the number of marks already present. Maintain
sum of their start counts at least floor(m*m/4)+1 whenever m>0; for m=0
maintain only the sum0.

This holds in the initial triangle for m=1,2,3, whose marked sums are m.
A mixed extension adding one mark gains at least ceil(m/2) by (5), exactly

    floor((m+1)^2/4)-floor(m^2/4).

For m=0 its positive start count supplies the separate +1. An unmarked
mixed extension preserves the assertion.

For a two-vertex extension, let r in {0,1,2} be the number of newly marked
vertices. If r=0, monotonicity suffices. If r=1, the gain H(U)>=|U|>=m
is at least ceil(m/2), with positivity covering m=0. If r=2 and m>=1,
the gain 2H(U)>=2m>=m+1 equals or exceeds

    floor((m+2)^2/4)-floor(m^2/4)=m+1.

For m=0, the gain is at least6, exceeding the required2. This proves the
invariant and (1). The singleton boundary is immediate; there is no strong
order-two tournament.

For (2), apply (1) to all m vertices with start<=k. If m>=4k, then
floor(m*m/4)+1>mk (already strict at m=4k), contradicting the upper bound
sum starts<=mk. Hence m<=4k-1.

## 2. A complete bounded-budget SCC decomposition

Let (T,v) be strong, |T|>=3, with s_v(T)=k. Decompose T-v into its uniquely
ordered strong components

    C, D1, ..., Dt,

where all earlier components dominate later ones. Put

    A={a in C: v->a},  q=product_i H(Di),

with q=1 when t=0. Every root-v path first enters A and then traverses the
components in order, giving the exact formula

    k=q * sum_{a in A} s_a(C).                           (7)

Thus A is nonempty, q divides k, and q<=k. Conversely, choose any strong C,
a nonempty A subset of C, and an ordered list of strong tail blocks Di.
Add v beating exactly A in C and choose its incidences to the tail blocks
arbitrarily. If t=0, require A proper in C; if t>0, require some vertex of
Dt to beat v. These are exactly the conditions making the resulting T
strong. They ensure v reaches C and hence every tail block, while every
vertex reaches v through the final block (or directly through C when t=0).
The root-start count is exactly (7).

Theorem1 shows that every eligible member of A belongs to a set of at most

    4(k/q)-1 <= 4k-1                                   (8)

vertices of C. Indeed its individual start count is at most k/q. Every other
vertex of C beats the new root v. The entire nonuniform root interface with C is therefore bounded.

**Actual-entry corollary.** Applying Theorem1 to A itself, rather than to
all eligible vertices, gives

    floor(|A|^2/4)+1 <= k/q,  hence |A|^2 < 4k/q.       (9)

Consequently |A|<=10 whenever k<=26. This bounds the entry set actually
chosen at a junction; the larger eligible reservoir is still bounded by103.
The two bounds refer to different sets.

Each nontrivial tail component D has

    3 <= |D| <= H(D) <= q <= k.

There are at most floor(log_3 k) such components. Their total order is at
most q when at least one exists, because each H(D)>=3 and the sum of these
factors is at most their product. For fixed k, the nontrivial tail blocks
therefore come from a finite literal catalogue: all strong tournaments of
orders3 through k with H<=k. Arbitrarily many singleton tail blocks are
still allowed and cannot be deleted here: their incidences with v can change
B(T,v), so weak-tail profile deletion is not being applied inside T-v.

## 3. Where recursion descends and where it does not

Choose any a in A and regard (C,a) as a smaller rooted source.

- If q>1, then s_a(C)<=k/q<=k/3<k.
- If |A|>=2, then s_a(C)<=k/q-1<k.
- The sole non-descending case is q=1 and A={a}, with s_a(C)=k.

In that last case all tail components are singletons. The graph is obtained
literally from (C,a) by J (add a vertex beating only a and mark it), followed
by the tail vertices in order, each using O or Z while retaining the new
root. The singleton-C / rooted-edge boundary is handled literally, not by
applying the nontrivial-source singleton-Z formula.

Thus every strong rooted source with start at most K has a complete
budget-recursive construction with:

- at most4K-1 eligible entry vertices at every strong core;
- a finite catalogue of nontrivial tail blocks of order at mostK;
- same-budget chains exactly of the physical J/O/Z form;
- strictly smaller start budgets at every multi-entry or nontrivial-tail
  junction.

This is a complete physical decomposition, not yet a finite-state arithmetic
compiler. In particular, one must still transport the simultaneous marked
responses through the construction, and prove which integer readouts its
actual graphs attain. Bounded entry ports do not make arbitrary formal
nonnegative path-cover vectors physical.

## 4. Relevance and exact stopping point for F(8,137)

The staggered condition requires s in {1,2,4,5,...,25,26}, correlated with
B. The audited general profile order theorem makes every fixed target a
finite existence question. The present result separately bounds the
nonuniform initial-SCC interface by103 ports, regardless of the target B.
It also identifies every failure of the same-budget O/Z/J recursion as a
strictly lower-budget multi-entry or finite-tail junction.

This is a possible input to a complete inductive grammar or target-dependent
rebuilding method. It does not assert occupancy at any new B, does not turn
a finite scan into F(8,137), and does not prove a uniform finite list of
arithmetic states. F(8,137) and the eventual threshold remain OPEN.

## 5. Author controls

The standalone verifier enumerates all labeled tournaments through order6,
uses subset DP to count every root-start response, and checks (1) for every
subset by checking the sums of the m smallest responses. It separately
checks the strong growth construction on all strong inputs, and the marked
increment bound on every mixed one-vertex extension of strong inputs through
order5. The all-order arguments above, not the finite range, establish the
claimed scope. Normal and debug receipts must agree.
