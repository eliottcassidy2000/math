# Pure-gap continuation, the missing higher response, and composition

**Status: PROVED identities and scoped obstruction; FINITE-EXACT discovery
search. Arbitrary odd-target coverage remains OPEN.** The every-even-gap,
one-edge interpolation, target-preserving rebases, chronological caterpillars,
and congruence-exact square-root approximation are inherited from Session199.
The zero-defect classification in STRUCTURE.md is inherited from Session231.
This note does not reclassify those results as new.

## 1. Exact common-continuation boundary

For the Round3 pair, write the upper-minus-lower state as
`delta*(1,0,0,0)`, where delta=2U>0. Every common O/Z word preserves this
vector. Conversely, a common word containing J cannot finish with equal
three-coordinate rooted sidecars.

Indeed the first J sends `(1,0,0,0)` to `(1,1,1,1)`, the actual singleton
state. Every later letter acts on this difference by the same physical
transfer law as on a singleton construction. Its deletion count is at least
one: O and Z preserve that count, while J replaces it by the positive total
Hamilton-path count. Multiplication by delta therefore leaves a strictly
positive deletion difference. This proves the **iff**, for arbitrary finite
common O/Z/J words, rather than only showing preservation for O/Z.

It does not prohibit asymmetric continuation, rerooting, or non-word graph
operations.

## 2. A literal common deletion graph and its two cuts

Let S be any nonempty rooted tournament with root r and state `(D,E,U,d)`.
The two outputs are the actual words OJJJJOZ and OJOZJOZ. Delete their last
J-born root. The remaining graphs are isomorphic, by swapping the third and
fourth newborns and fixing every other vertex.

Here is a root-independent description of that common deletion graph. Start
with

    B0 = S ⇒ (t ⇒ C3[a→b→c→a]),

and reverse just the arc r→a, so a→r. Call the resulting graph B. Then the
common deletion graph is

    C = B ⇒ (x→y).

The lower cut is `{b,x}` and the upper cut is `{a,x}`. A cut records the
outneighbors of the removed root. The reversal is the only arc from the
four-vertex gadget back into S. This literal description follows by reading
the two newborn histories, not merely from equality of numerical coordinates.

Use the inherited seven-coordinate interface `(h,q,A,G,T,S0,E0)`:
h is the deletion graph's Hamilton-path count; q counts ordered two-path
covers, either path permitted empty; A and G are its start-selected and
end-selected two-path responses; T requires the first path to end in the cut
and the second to start in the cut; S0 and E0 are its full-path starting and
ending cut sums. All paths are directed and Hamiltonian on their own parts.

For the source define additionally

    Q = sum_X H(S[X]) H(S[V\X]),
    V = sum_X e_r(S[X]) H(S[V\X]),
    W = sum_(X,Y,Z partition V) s_r(S[X]) H(S[Y]) H(S[Z]).

Here s_r and e_r are zero if r is absent. The distinguished X path in W is
nonempty and starts at r; the other two paths are ordered and may be empty.
Thus W is the rooted ordered-three-path response. These definitions apply to
arbitrary rooted S, not only the response-selector family. In particular,
source reversal exchanges U and V, but W is not a function of the old four
coordinates asserted here.

**Theorem.** The upper-minus-lower deletion interfaces differ by

    (0,0,0,−12W,−2U,0,0).                         (1)

Their common coordinates include

    h = 3(D+2U),
    q = 4(24Q+42W−6V),
    A = 18D+44U,
    E0 = 0.                                      (2)

The individual T values are `4D+10U` and `4D+8U`. In particular h and q
are divisible by three. Equality of the small rooted sidecar does **not**
mean equality of the seven-coordinate interface: it misses the response W.

### Proof of the coefficient −12

First work on B before the transitive x,y suffix. In B0 the contributions
to G from endpoints a and b agree, by the cyclic automorphism of C3 fixing
t and every source vertex. Reversing r→a removes covers using r→a and adds
covers using a→r. Every changed cover uses that edge, so this is exhaustive.

For a removed cover, each of the two source portions is the initial segment
of its path. The source portion ending r is followed by a gadget path
starting a. In K=t⇒C3, the three possible a-starting portions are
`a`, `ab`, `abc`; their complements have respectively the unique paths
`tbc`, `tc`, `t`. If this is the distinguished first path, exactly one
cover ends at a and exactly one ends at b. If it is the other path, neither
of the two endpoint tests is satisfied. Thus the removed contribution to
G_a−G_b is zero, for each possible source cover.

For an added cover, the segment immediately after a is a nonempty source
path starting r. The two other source portions are precisely the initial
segments of the two cover paths, and may be empty. These three portions
are a unique object counted by W. Conversely, every such source object
can be joined to a gadget skeleton without dependence on its internal
arcs: the initial source portions avoid r and beat all gadget vertices;
the post-a portion may exit to t,b,c, while a is already used.

Replace the nonempty r-starting source portion by a token z. The skeleton
is a two-path cover of K+z using a→z, where z beats t,b,c. There are exactly
12 such skeletons whose first path ends b, and none ending a (a already has
its outgoing path edge to z). The finite list is checked by the independent
permutation-and-cut enumerator in verify_interfaces.py, not fitted from
source examples. Consequently G_a−G_b=−12W on B.

Passing to C=B⇒(x→y) and selecting x in both cuts adds the same quantity
to G, so the difference remains −12W.

### Proof of q and the other coordinates

The baseline two-path count is
`Q(B0)=Q(S) Q(t) Q(C3)=24Q`. For removed covers there are two choices of
which source path ends at r and three a-starting gadget skeletons, giving
`6V`. The same added-cover decomposition used above has **42** ordered
skeletons when no endpoint is prescribed, giving `42W`. Therefore

    Q(B)=24Q+42W−6V.

The x→y suffix has two-path count four, proving (2)'s q formula. The
42 and 12 constants have separate exhaustive finite skeleton certificates.

The already verified four-coordinate output matrices give h and A: the
root-deletion count is `3D+6U`, and rooted response is `21D+50U=h+A`.
Every full path of C ends at y, which is outside both cuts, so E0=0.
For the one-root extension the exact count is `h+A−T−E0`: paths ending at
the added root contribute h−E0, and paths continuing from it contribute
A−T. The known total-count difference 2U therefore gives ΔT=−2U and the
displayed individual T formulas.

Finally S0 is the same for the two cuts. No full path of C starts at x.
A full path of B starting a has forced gadget order
`a, [source path from r], t,b,c`; a full path starting b has forced order
`b,c,a, [source path from r], t`. Both are in bijection with full source
paths starting r. Thus ΔS0=0. Graph isomorphism supplies Δh=Δq=0, and the
small-sidecar identity supplies ΔA=ΔE0=0. This proves all seven coordinates
of (1).

## 3. Exact affine choice cubes, with a growing divisor

For physical nonempty interfaces the inherited ordinal law is

    (h,q,A,G,T,S0,E0)*(j,s,B,F,L,V0,W0)
     = (hj,qs,As+hB,Fq+jG,AF+Tj+hL+E0V0,S0j,hW0).

Close the deletion graph by adding an unselected sink v, then a root that
beats exactly the selected cut and loses to all other vertices, including v.
The count is `h+2A−T`. If every source root reaches all source vertices, the
closure of any finite ordinal chain of the modules above is strong: the
root reaches a selected vertex in the first module, hence that whole module
and all later modules, while every old vertex reaches the root through v.

Put k such modules in order, with arbitrary independent source parameters
U_i,W_i. The two choices in each module retain h_i,q_i,A_i,S_i,E_i and change
only G_i,T_i as in (1). Every output scalar is therefore exactly

    H(epsilon_1,...,epsilon_k)
      = H(0,...,0) + sum_i epsilon_i c_i,

    c_i = 2 h_suffix(i) [h_prefix(i) U_i + 6 A_prefix(i) W_i].       (3)

For an empty prefix use h_prefix=1,A_prefix=0, and for an empty suffix use
h_suffix=1. Formula (3) follows by changing a single module: multiplying on
the left changes ΔT to `A_prefix ΔG+h_prefix ΔT`, and multiplying on the
right then multiplies it by h_suffix. All coefficients are independent of
other choices. This is an actual affine physical cube, not just a gcd model.

**Divisibility consequence.** Every difference of two scalar values in this
k-module cube is divisible by

    2·3^(k−1).                                             (4)

Each h_i and q_i is divisible by three, by (2). Hence a prefix of j modules
has h divisible by 3^j and A divisible by 3^(j−1): use
`A(PQ)=A(P)q(Q)+h(P)A(Q)` inductively. In (3), the h_prefix U_i term and the
6 A_prefix W_i term, after the suffix product, are both divisible by
3^(k−1). The first-module case has A_prefix=0 and the same conclusion.
This proves (4). In particular, two or more of these modules cannot by
themselves furnish consecutive odd counts in a fixed cube.

This does not forbid a union of separately placed cubes, an asymmetric
rerooting scheme, or different modules with coprime deletion data. It is a
scoped obstruction, unlike an assertion that all ordinal compositions fail.

## 4. Approximation plus repair: the remaining positive bridge

The relevant inherited positive theorem is Session199
`session199/congruence-intervals/PROOF.md`, Theorem3. For each **fixed** m,
every sufficiently large odd target N has a constructed balanced strong C1
count H<=N with H≡N mod m and `N−H=O_m(sqrt(N))`. This is stronger than
modular density or existence of unbounded witnesses and must be retained.
It does not supply an embedding of the pure-gap compiler inside that H.

An exact repair would need a physical representation of the particular
approximant H containing a switchable module whose low state gives exactly
H, and whose retained source response equals half the required deficit
(after any known context multiplier). Equation (3) states the precise
context multiplier and the additional W response. Merely constructing an
independent pair with that deficit supplies neither the embedding nor the
low-state equality. A new module's baseline is not zero: for the response
selector its low total is `17g(U)+40U`, at least `40U+17`.

The modulus can accommodate (4) for any **fixed** number of modules, so the
cube divisor alone does not refute every approximation-and-repair plan.
Conversely, setting m=3^(k−1) with k growing with the target requires the
explicit m-dependent approximation constants; it is not the fixed-m theorem.
A successful proof must control the actual baseline and the correction
interval simultaneously. No physical baseline subtraction or all-deficit
repair is claimed here.

## 5. A positive bounded-error channel with coprime deletion

Dropping the opening O genuinely escapes the common factor three. For the
actual words JJJJOZ and JOZJOZ, the four-coordinate difference is

    (ΔD,ΔE,ΔU,Δd) = (3E+U,0,4E,0),
    common d = 3D−E+3U.                                (5)

The same common-deletion proof now uses C3 without the preceding t vertex.
The added-edge skeleton has14 ordered two-path covers, of which4 end at b;
the removed-edge endpoint difference is still zero. Therefore the exact
seven-coordinate difference is

    (0,0,4E,−4W,E−U,0,0),                              (6)

and the common h,q coordinates are

    h=3D−E+3U, q=4(12Q+14W−6V).

The ΔA coordinate follows from the full rooted response in (5), while
ΔT=ΔA−ΔD=E−U follows from the one-root closure identity. This is a
physical composition law, with its additional response retained.

### A shorter strong channel, and the cost of synchronizing endpoints

Before the final OZ, the pair JJJJ/JOZJ already has

    Δstate=(E+U,U,E,0),
    Δinterface=(0,0,E,−4W,0,0,−U).                    (5a)

Its common deletion graph is B, with singleton cuts {b} and {a}. Thus T=0
in both: the same selected vertex cannot end one disjoint path and start
the other. This also gives ΔE0=−U and ΔA=E directly from the rooted profile.
The q and h formulas above lose only the terminal factor four in q.

For any n>=1 take U=2n−1,E=1 and D=g(U)<=U, using the singleton at n=1.
Both four-letter outputs are strong and have gap2n, response discrepancy1,
equal deletion count coprime to three, and upper count

    5g(U)+7U−1 <= 24n−13.

The first pair is9/11. Its endpoint discrepancy is2n−1, **not one**. An
initial informal proposal incorrectly conflated those coordinates; direct
full-matrix subtraction repairs it before any theorem uses the proposal.
The strongness follows because J preserves root reachability and, once its
parent has at least two vertices and a reachable marked root, adds a new
root that is reached through an old nonroot. The O/Z branch is strong by
its final Z before the final J; singleton boundary steps are immediate.

Appending O to each graph preserves the scalar gap and yields
`Δstate=(2n,0,2,0)`. These graphs have a sink, so they are not strong.
Appending Z then restores strongness and gives `(2n+2,0,4,0)`. This explains
precisely how endpoint synchronization produces the six-letter channel,
rather than suggesting that the shorter pair already synchronizes it.

### Every even gap at least four with response error exactly four

For every integer n>=2, take the inherited canonical response source with
U=2n−3 and E=1; for n=2 use the singleton (U=E=1). Such a source has
D=g(U)<=U, and its root reaches every vertex. Both words above give strong
outputs with

    gap = 2n, endpoint difference = deletion difference = 0,
    response difference = 4,
    common deletion = 3g(U)+3U−1 ≡ 2 (mod 3),
    upper count = 17g(U)+21U−3 <= 76n−117.              (7)

Thus these are universal prescribed even gaps with a *constant*, rather
than gap-proportional, response discrepancy, and a deletion count coprime
to three. The n=2 pair is31/35; n=3 gives105/111. Gap two is available from
the pure-D compiler, but is not claimed from the E=1 formula.

After any common O/Z tail of length m, put

    a = sum_(i=0)^(m−1) [letter i is Z] 2^i.

The exact four-coordinate difference becomes

    (2n+4a,0,4·2^m,0).                                (8)

Every a=0,...,2^m−1 occurs, so at fixed n,m this gives every gap in that
step-four interval with the same response discrepancy and the same deletion
count coprime to three. The physical anchors also depend on the tail; (8)
is not arbitrary absolute placement.

### Exact contextual increment

For a module as in (6) between fixed physical contexts X and Y, the protected
closure count changes by

    h_X[4E(2q_Y−G_Y)+(U−E)h_Y] + 4W A_X h_Y.          (9)

For empty contexts use h=1,q=1,A=G=0. Formula (9) follows by multiplying the
seven-coordinate law before applying h+2A−T. In a chain, a flip in an
earlier module i and another in a later module j has the additional positive
mixed coefficient

    16 E_i W_j · (product h before i and after j)
                  · (product q strictly between i and j).

There are no higher-degree terms: A and G are affine in the choices, while
T contains their ordered pair products. This is an exact physical quadratic
choice model, not a pure affine cube. It supplies an explicit non-divisible-
by-three family that the preceding obstruction does not cover.

### Complementary switches cannot erase the response error

The proposed cancellation test has an exact answer for a fixed
ordinal chain of the E=1 modules, each allowed its own common O/Z tail.
Let tail i have length t_i. It preserves odd h_i and changes q_i to
`2^t_i q_i^base`, where q_i^base is even. The small rooted response difference
and the seven-coordinate A difference are both `4·2^t_i`.

Consequently the contribution of binary choice i to the final A is

    c_i = 4·2^t_i (product_(j<i) h_j)(product_(j>i) q_j).           (10)

These coefficients have strictly decreasing two-adic valuations. Indeed

    v2(c_i)−v2(c_(i+1)) = t_i+v2(q_(i+1)^base) >= 1.

The oddness of h_i follows directly from D,U odd and E=1. The two-path
count of any nonempty tournament is even, by pairing each ordered split
with its swapped split, so the right-hand valuation is positive.

For any two distinct binary choices, the rightmost changed coordinate has
the unique smallest two-adic valuation in their signed difference of (10).
That sum cannot vanish. Therefore **no two distinct choices in this fixed
chain have the same A**, and complementary switches cannot cancel the
response discrepancy. Common nonempty external contexts only multiply the
A difference by nonzero h_left and q_right, so they cannot repair it.

The same proof applies to the raw four-letter pair with arbitrary common O/Z tails: replace the leading4 in (10) by1 and measure t_i from that raw pair. Thus even the response-discrepancy-one strong channel cannot be canceled by complementary choices in a fixed ordinal chain.

The positive survivor is a faithful dyadic address of every choice, with an
explicit largest differing bit. The obstruction concerns exact response
cancellation in this particular ordinal-chain architecture. It says nothing
against scalar collisions, a separately placed chain, changing the source
modules along with their contexts, or general graph rerooting.

### Longer pure-channel search and a corrected quotient loss

The exploratory escape_divisor_probe.py exhausts the distinct O/Z/J
matrices at each length through14 and independently retains every matrix
having a terminal-Z witness. It tests unit coefficient differences
+2D,+2E,+2U,+2d with equal three-row sidecars. At length14 there are792415
distinct matrices; the accumulated terminal witnesses yield255 such pairs,
all +2U, and none has a deletion row not divisible by three. These are finite
matrix statements only, not an all-word rigidity theorem.

**Correction lineage.** The initial search kept only the first word for a
matrix and then filtered that representative by its final letter. The audit
observed JZ=OJ: choosing OJ loses the terminal-Z witness JZ. This made the
initial terminal-admissibility completeness claim invalid, though it did
not affect a displayed pair's verified identities. The repaired search
records Z(a) from every preceding-level matrix *before* representative
selection. The length14 run was repeated after this repair; the reported
unit-pair counts and no-escape outcome happened to remain unchanged. No
claim relies on the faulty representative filter.

## 6. Verification and scope

Run `python3 -I -B verify_interfaces.py`, then repeat under `-O`; compare
both outputs with INTERFACES_CHECK.json. The verifier independently counts
ordered gadget skeletons, literally recounts both seven-coordinate deletion
interfaces on65 source/root cases, verifies all common words through length8,
checks all choices in six fixed physical pure-gap module cubes, verifies16 prescribed bounded-error gaps and all their common tails through length6, and checks dyadic response injectivity through seven modules. It uses explicit
exceptions and has a failing `--fail` sentinel. The infinite statements
follow from the displayed combinatorial decompositions and inductions.

Live concept board: (i) pure-D small-sidecar transport; (ii) the higher root
response W; (iii) affine ordinal choices with actual baselines; (iv) inherited
congruence-exact sqrt approximation; (v) rerooting that escapes deletion
factor three at the expense of a new response coordinate. The new bridge
restores the missing coordinate and exact composition map. It does not
suppress the independently prescribed anchor, which is the unresolved step.
