# Mark a finite interaction graph, sum the unmarked world exactly

**PROVED in this note; prototype FINITE-EXACT controls.** This is an independent
combinatorial derivation of the fixed-depth coefficient closure in the sibling
[finite moment-jet compiler](../enumeration/PROOF.md). It is not a faster practical
counter, nor a proof of Chan–Pak. No literature-priority claim is made.

## 1. Inheritance and idea

Session278 `enumeration/PROOF.md` supplies the odd-cycle Burnside gas and its
signed elementary-symmetric coefficients C_k. It asks for low coefficients
without unrestricted partition enumeration. Session279's backbone replaces the
min kernel by gcd on every pair touching a chosen finite set B of odd lengths.
The present argument is a standard labelled-marking construction applied to
that gas. It gives a different explanation of why the moment jet closes:
a degree-k correction sees at most 2k exceptional cycle instances. All other
cycles remain an exactly summable external-field background.

The marked objects are cycle *instances*, not merely distinct lengths. Confusing
these two would lose repeated-length multiplicities and the m! normalization.

## 2. Formula

Fix B, not containing 1, and condition its multiplicities t_b. Let
N=n-sum(b t_b) and let w_B be its exact internal Burnside weight. All other
cycles have pair baseline min(a,b). A pair touching B has its exact gcd
baseline. For off-backbone lengths define

    q(a,b)=1-2^{gcd(a,b)-min(a,b)}.

For k>0, m ranges from 2 to 2k. Let G_{m,k} be all simple graphs on the labelled
vertex set [m] with k edges and no isolated vertices. Sum over ordered odd
length tuples A=(a_1,...,a_m), each a_i>=3, not in B, and sum A<=N. Define

    P_k(A)=sum_{G in G_{m,k}} product_{ij in E(G)} q(a_i,a_j),
    E_A=sum_i (a_i-1)/2 + sum_{i<j} min(a_i,a_j)
         +sum_i sum_{b in B} t_b gcd(a_i,b),
    W_A=2^{E_A}/product_i a_i.

Let R(s;A,t) be the total weight of all **unmarked** off-backbone cycles of
total length s, including their interactions with A and the backbone. Then

    C_k(n;B)=sum_t w_B sum_{m=2}^{2k} (1/m!)
               sum_A W_A P_k(A) R(N-sum A;A,t).             (1)

For k=0, use the empty marking: C_0=sum_t w_B R(N;empty,t).
Singleton marked vertices may be omitted since every q(1,a)=0 and isolated
vertices are forbidden.

### The factorial is exact even when lengths repeat

Fix an unmarked sector with off-backbone multiplicities u_a and a chosen set
of k distinct defect edges. Its incident vertices have multiplicities r_a and
total m. There are (u_a)_(r_a) ways to inject each labelled marked group of
length a into the original cycle instances. Summing all ordered length tuples,
all labelled incident graphs, and dividing by m! counts the original edge set
once. The weight changes from 1/product u_a! to
1/product (u_a-r_a)! precisely by these injections; factors 1/a for the marked
cycles supply the remaining cycle-index denominator. Thus the contribution is
W times the chosen product of q's. Summing chosen edge sets gives W e_k(q).
No automorphism-order guess or division by a graph stabilizer is used.

## 3. External-field size recurrence

For a fixed (A,t), initialize D(0)=1. Process odd b not in B in increasing
order. For each old size s and j new cycles of length b, update size s+bj by
multiplication with

    2^[j s + j(b-1)/2 + b binom(j,2) + j F_b] /(b^j j!),
    F_b=sum_i min(a_i,b)+sum_{a in B} t_a gcd(a,b).          (2)

The first term j s is exact because every previously processed unmarked cycle
has length at most b. Interactions with marked cycles are already fixed in F_b;
their lengths need not be earlier than b. The final coefficient at size r is
R(r;A,t). Marked and unmarked cycles may have the **same length**. Omitting that
case would destroy (1).

Equations (1)-(2) compute C_k with no unrestricted partition enumeration.
For fixed |B|=b and k, at most O(n^b) backbone choices and O(n^{2k}) ordered
length tuples occur. The number of labelled graphs is a k-dependent constant;
each size recurrence costs O(n^3) rational arithmetic operations (loose bound).
Hence a sufficient bound is

    O_{B,k}(n^{b+2k+3}) arithmetic operations.              (3)

The constants include graph enumeration and can be enormous. This is a fixed-
parameter polynomial statement, not a uniform polynomial bound when the desired
correction depth grows with n. Bit costs are not constant: exponents are O(n²),
centralizer denominators of recombined sectors divide n!, and intermediate
marked terms may include an additional m! denominator. Fixed-depth marked
rational arithmetic has polynomial bit length. No claim of a small bit-complexity exponent follows.

## 4. Independent coefficient controls and boundary

Run from repository root:

    python3 05-knowledge/results/session279_interval_tightening_20261112/marked/marked_graph.py

The producer uses graph markings and a size DP. The test harness separately
uses full odd partitions only as a reference. It compares C0,C1,C2 for every
n=0,...,18 in both the plain and B={3} models; additional C3 checks use n=14,18
plain and n=22 with B={3}. These include repeated-length stars, disjoint marked
edges, a triangle of marked interactions, singleton padding and empty sectors.
The k=3 first stars (3,3,3,5) and (5,5,5,7) retain four cycle instances despite
having only two length labels.

This solves the **fixed-depth coefficient computation** problem posed in
Session278 in a mathematically explicit way. It does not solve its all-size
rounding problem: the inherited fixed-depth absolute-error obstruction still
applies. An interval consumer must choose deeper corrections, more exact
backbone information, other certified bounds, or a residue as n grows.

## 5. How this differs from the sibling jet

The jet stores occupation moments up to a correction-dependent degree during
one length sweep. Formula (1) first exposes the finite graph that made those
moments necessary, then sums the background in its external field. They are two
representations of the same finite exceptional carrier, not independent claims
that a universal finite alphabet classifies tournaments. The practical jet avoids
much redundant marking work; the graph form makes the 2k-vertex limit and all
multiplicity factors visible. Both retain physical cycle incidence before
compressing the rest.

## 6. A complete no-reference producer at order64

With B={3,5,7,9,11}, every residual cycle has length at least13. At n=64
there are at most four such cycles and therefore at most six defect edges.
Equation (1) through C6 is the entire finite expansion, not an approximation.
Run:

    python3 05-knowledge/results/session279_interval_tightening_20261112/marked/marked_graph.py --order64

This producer computes all seven coefficients by graph marking and size DP,
without enumerating unrestricted partitions or reading a reference value.
`order64.json` records the result. An independent direct Burnside sum checks
it only after production. The two final correction coefficients are

    C5=183075581474758003459948593966307344384000/126293,
    C6=25039070082761392429033586688000/19.

This independently completes the same integer isolated by the sibling ground-
moment certificate. Its method differs: it computes the two higher correction
coefficients rather than retaining the ground-defect sidecar. The comparison
therefore concerns complementary exact representations, not identical work.

More generally, if the smallest off-backbone length is a, n permits at most
floor(n/a) residual cycles; depth binom(floor(n/a),2) is certainly exact. This
simple support bound gives finite termination at each n, but keeping it small
as n grows requires a growing backbone. It does not turn (3), whose constants
assume fixed B,k, into a polynomial-time all-order exact counter.
