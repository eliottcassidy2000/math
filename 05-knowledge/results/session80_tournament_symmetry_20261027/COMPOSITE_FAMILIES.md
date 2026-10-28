# From the aggregate prime bank to canonical composite-family data

**Status:** PROVED deductions from classical tournament modular decomposition;
FINITE-EXACT checks described below. No priority claim. This is an extension of
Session79's aggregate recovery, not an independent faster computation of A000568.

## Inheritance and question

Input: [the prior prime recovery](../session79_tournament_enumeration_20261027/PRIME_CYCLE_INDEX.md) supplies the
aggregate prime cycle index P, with prime orders at least three. The closest
proved mechanism is canonical strong-module substitution, not arbitrary
substitution trees. The canonical hostile is flattening transitive nodes or
forgetting the action on prime slots. The corrected near miss is the unsigned
cyclic-triple symmetry claim rejected by the inherited composition work.
The least-used sidecar pursued here is the full power-sum index, rather than
only its isotype-count specialization.

Concept board: canonical linear/prime nodes; weighted Adams operations; forward
filtration; maximum prime-depth; information lost by aggregate averaging.
Anchor: generate full composite indices from P. Niche: node-arity and depth
statistics. Wildcard: the thin extremal-depth boundary.

## 1. Universal canonical-tree equation

Work in the completion by vertex weight, wt(p_i)=i. All auxiliary marks have
vertex weight zero, and coefficients at each vertex weight are polynomials in
them. Let P_k be the degree-k part of P. Mark a prime node of arity k by u_k,
and a linear node with k ordered children by v_k. Let A be the weighted cycle
index of nonempty tournaments, S that of strong tournaments including the
singleton. Then

    S = p_1 + sum_(k>=3) u_k (P_k circ_w A),
    A = S + sum_(k>=2) v_k S^k.                         (1)

The operation circ_w applies to each inner cycle index the Adams operation

    psi_d(p_i)=p_(di),  psi_d(u_k)=u_k^d,  psi_d(v_k)=v_k^d.

The new root mark u_k is multiplied *after* composition; it is one root,
not d roots. The same applies to the new linear-root mark v_k.
For unmarked P_k, this definition simply means that every p_d in P_k is
replaced by psi_d(A). External scalar rational coefficients remain fixed.

Why (1) is canonical: maximal proper strong modules give either a prime root
whose children are arbitrary nonempty tournaments, or a transitive root whose
ordered children are strong tournaments. There are no adjacent linear nodes.
Thus linear-root arity is the number of maximal SCCs, never an arbitrary binary
parenthesization. A strong nonsingleton tournament has a prime root but need not
itself be prime. A singleton has no internal node.

Weighted fixed-point explanation: a root permutation d-cycle repeats one child
structure d times. Every vertex count and every additive child statistic is
therefore multiplied by d. The root itself occurs once. This proves the marking
rule directly, without assuming a weighted-species theorem as an extra input.

**Unique forward solution.** With A(0)=S(0)=0 and degree-one term p_1, (1)
has one solution. The degree-n root contribution uses only child degrees below n:
a prime has at least three nonempty children, and a linear node at least two.
Equivalently, start A=0 and iterate (1); each iteration fixes at least one more
vertex degree. Therefore P through degree N determines the entire marked A
through degree N, with no composite census. In the prime term its dependence
on an unknown inner degree-n term cannot start before degree n+2.

Setting all marks to 1 gives A=S/(1-S), hence the inherited identities. The
unmarked forward map is bijective between formal series P of valuation at least
3 and A with A=p_1+p_1^2+terms of degree >=3: its inverse is

    P = p_1/(1+p_1) - A^<−1>.

This is a formal-algebra bijection, not a claim that every rational P is a
realizable species. Nonnegative coefficients alone also do not prove species
realizability. Formal inverse plus the canonical-tree construction establish
exact invertibility at the aggregate level; they do not recover a prime atlas.

## 2. What composite-family outputs are now determined?

At marks 1, write S=p_1+P circ A. Then the disjoint root classes are

    singleton:                    p_1,
    prime tournament itself:      P,
    composite strong tournament:  P circ A - P,
    linear-root tournament:       A - S = S^2/(1-S).

Here 'composite' means order >=2 and not one of the designated primes; the
2-vertex tournament is included. Hence all composites have index A-p_1-P.
These are full cycle indices, not just ordinary counts. Specializing p_i=x^i
gives isotype counts; setting p_1=x and p_(i>1)=0 gives labelled counts divided
by n!. Coefficients times z_lambda give fixed points of a specified permutation
of type lambda within each root class and each marked tree stratum.

Equation (1) recovers joint distributions of all prime/linear arities. In every
nonzero monomial with n leaves, b_k prime k-nodes and l_k linear k-nodes,

    n-1 = sum_(k>=3)(k-1)b_k + sum_(k>=2)(k-1)l_k.           (2)

This is the tree edges identity (sum of outdegrees equals number of vertices
minus one). In particular 2b+l <= n-1, with equality iff every prime node has
arity 3 and every linear node arity 2. The checker coarsens to u_k=u and v_k=v
and verifies this necessary support bound.

Leaf colors are also admissible: replace the singleton p_1 by
p_1(c_1+...+c_r), extend psi_d(c_j)=c_j^d, and retain P. This enumerates
vertex-colored tournaments under color-preserving isomorphism. It does not
supply an unseen prime invariant such as self-converseness; marking such a
prime property requires the corresponding split prime cycle-index bank.

## 3. A hostile showing why ordinary marking fails

Inflate every slot of C3 by C3. The resulting nine-vertex tournament has four
prime nodes and no linear nodes. Its cycle index is

    u * ((u Z_C3)^3 + 2 psi_3(u Z_C3))/3.

Correct Adams sends u to u^3, so after unlabelled specialization the node weight
is exactly u^4. Holding u fixed instead produces

    (1/3)u^4 + (2/3)u^2,

which is not even an integral distribution of isotypes. This wrong expression
passes the unmarked u=1 check. Thus an unmarked round-trip alone is not an audit
of marked composition. `check_marked_composites.py` contains both controls.

Another hostile: P_3[A] alone counts all strong root-C3 inflations, not only
prime tournaments. A C3 with one slot replaced by a two-vertex transitive block
is strong and composite at order four. Equating 'prime root' with 'prime object'
would lose this entire family.

## 4. Maximum prime-depth, not the number of prime nodes

Define h(U) as the maximum number of prime nodes on a root-to-leaf path in the
canonical tree. Linear nodes count zero. Let H(S)=S+sum_(k>=2)v_k S^k. Then

    A_0 = H(p_1),
    A_h = H(p_1 + sum_(k>=3)u_k (P_k circ_w A_(h-1))) (h>=1). (3)

This is an exact depth filtration. A_0 consists precisely of transitive
tournaments (including the singleton), rather than only leaves. Formula (3)
works because every prime-root child has one fewer available prime level, while
an outer linear node does not consume a level. A_h-A_(h-1) gives exact depth h.
For n vertices, h<=floor((n-1)/2), so the filtration stabilizes by that level.
Depth is a maximum statistic, not an additive variable to be raised by Adams.

**Thin-boundary theorem.** For every integer h>=1:

* exactly one tournament on 2h+1 vertices has prime-depth h;
* exactly 4h-1 tournaments on 2h+2 vertices have prime-depth h.

Proof: a path of h prime nodes consumes at least 2h+1 leaves: at each prime node,
there are at least two off-path children, each requiring a leaf. Equality forces
all prime roots to be C3, every off-path child to be a singleton, and no linear
nodes. The chain is unique because C3 acts transitively on its three slots.

For one additional leaf, no extra off-path prime node is possible. No prime
node can have arity four, because P_4=0 (independently certified by the exhaustive
small census, also immediate from the recovered index). A prime arity >=5
requires at least two additional leaves. Thus there is exactly one binary
linear node inserted into the unique C3 chain. It can be:

1. above its root: 2 placements, according to which ordered child is the chain;
2. on one of its h-1 internal prime-to-prime edges: 2(h-1) placements;
3. at a side leaf of one of its h-1 non-bottom prime nodes: 2(h-1) placements;
4. at a leaf of its bottom C3: 1 placement modulo rotation.

Sum: 2+2(h-1)+2(h-1)+1=4h-1. These cases are nonisomorphic by the canonical
rooted tree and the ordered linear-node children. At a non-bottom C3 node, its
unique nonsingleton chain child fixes a quotient vertex, whose C3 stabilizer is
trivial; consequently the two side leaves are distinct. At the bottom all three
leaves are equivalent. This identifies exactly the boundary where the count
changes. The deficit-two extension is proved next; the all-defect extension follows below.

### Deficit two: a second exact polynomial boundary

**PROVED.** The number of tournaments of prime-depth exactly h on 2h+3 vertices
is 10 for h=1, and

    c_h = 8h^2 + 11h - 10          (h>=2).                  (4)

For a scalar proof, specialize marks to 1 and p_i=x^i. Write a_h(x) for (3),
delta_h=a_h-a_(h-1), and

    delta_h = x^(2h+1) + (4h-1)x^(2h+2) + c_h x^(2h+3)
              + O(x^(2h+4)).

For h>=3 the inner difference delta_(h-1) starts at degree 2h-1. Terms quadratic
in that difference in P start at degree 4h-1>2h+3. Nonidentity Adams differences
start at degree at least 6h-3>2h+3. Quotient sizes >=6 first contribute to the
linearized difference at degree (2h-1)+5>2h+3. Thus, to the required order, the
scalar step is linear with multiplier

    K(x) = (a^2 + 11a^4 + O(x^5)) / (1-s)^2,
    a = a_(h-2) = x+x^2+2x^3+O(x^4),
    s = x + P[a] = x+O(x^3).

The derivative coefficients 1 and 11 come from the identity terms a^3/3 in
P_3 and 11a^5/5 in P_5. Hence

    K(x) = x^2 + 4x^3 + 23x^4 + O(x^5),
    c_h = c_(h-1) + 4(4h-5) + 23 = c_(h-1)+16h+3.

The h=2 starting value c_2=44 follows by expanding two steps from
a_0=x/(1-x), using only

    P_3[a]=(a(x)^3+2a(x^3))/3,
    P_5[a]=(11a(x)^5+4a(x^5))/5.

Those two steps give [x^5,x^6,x^7](a_2-a_1)=(1,7,44).
Prime quotients of order >=6 cannot appear in a depth-two object of size seven
(the minimal leaf count with one such root is at least eight). The seed is
therefore an exact coefficient calculation, independently executed in both
checkers, not an extrapolated fit. Summing the recurrence proves (4).
The h=1 exception records where the nonlinear difference terms have not yet
moved above the observation window.

An independent supplemental check used ordinary rational series and only
the known P_3 and P_5, rather than the multivariate marked implementation. It
checked all three boundaries through h=20 (sizes up to 43), with -O giving the
same receipt. Higher prime arities are safely absent *for these boundary
coefficients*: a depth-h path plus one arity >=6 node requires at least 2h+4
leaves. They are not absent from the full order-43 tournament population.

### All fixed-excess diagonals become polynomials

**PROVED.** Fix r>=0, and let f_r(h) count tournaments on 2h+1+r vertices
with prime-depth exactly h. Put h_0=max(1,floor((r+2)/2)). Then f_r(h), for
integer h>=h_0, agrees with a polynomial in h of degree exactly r and leading
coefficient 4^r/r!. All these diagonal families are generated by one kernel.

Let a(x) be the nonempty tournament isotype series, B(x)=1+a(x), and define

    K(x)=B(x)^2 (partial P / partial p_1)(a(x),a(x^2),...),
    L(x)=K(x)/x^2=1+4x+23x^2+186x^3+2332x^4+... .

The partial derivative is of the *full* multivariate prime cycle index. It
includes nonidentity terms containing p_1; it is not the derivative of only
P's identity-permutation part.

For delta_h=a_h-a_(h-1), normalize

    D_h(x)=x^(-(2h+1)) delta_h(x).

Then for h>(r+2)/2,

    D_h = L D_(h-1)  modulo x^(r+1),
    D_h = L^(h-h_0) D_(h_0) modulo x^(r+1) (h>=h_0).       (5)

Proof of the cutoff: put m=2h-1, the valuation of delta_(h-1). A second-order
prime substitution term has valuation at least 2m+1=4h-1. A nonidentity odd
Adams substitution difference has valuation at least 3m=6h-3. Linearizing the
outer map S/(1-S) introduces errors only at still higher degree. Replacing the
base a_(h-2) by the stabilized a in the derivative also has error at least
4h-1: their difference has valuation m, and the derivative's lowest term is
a^2. Thus the complete linearization error starts at 4h-1. It is invisible
through degree 2h+1+r precisely when 2h>r+2. This proves (5).

Since L=1+4x+O(x^2), expanding its integer power modulo x^(r+1) via the binomial
theorem makes every coefficient a polynomial in h of degree <=r. Also
D_(h_0)(0)=1 by the unique minimum-size chain. Only choosing r copies of 4x
contributes to the degree-r leading term; it gives 4^r/r!, proving exact degree.

Only primes of order <=r+3 are needed to generate the *whole* r-th diagonal,
even as h grows without bound. Indeed a prime k-node on a maximum-depth path
uses excess at least k-3, while a prime node off that path uses excess at least
k-1. This supplies finite prime input for infinitely many composite orders.
This does not say the full order-n census is determined by a bounded prime bank.

Some exact binomial forms (t=h-h_0; binomial(t,0)=1):

| r | h_0 | coefficients of sum_j c_j binomial(t,j) |
|---|---:|---|
| 0 | 1 | 1 |
| 1 | 1 | 3, 4 |
| 2 | 2 | 44, 51, 16 |
| 3 | 2 | 334, 523, 296, 64 |
| 4 | 3 | 9626, 9991, 5561, 1808, 256 |
| 5 | 3 | 158146, 148557, 80591, 33500, 8704, 1024 |
| 6 | 4 | 7745707, 4663655, 1875218, 742027, 228192, 44800, 4096 |

The theorem is all-order; displayed coefficients are FINITE-EXACT seed outputs
and are checked against independent scalar iteration through h=15. These
seeds determine the stated polynomials exactly, rather than merely fitting
values with no structural degree guarantee.

### A bounded linear carrier at fixed excess

For fixed r work over Q in V_r=Q[x]/(x^(r+1)). Multiplication by L is an
(r+1)-dimensional closed transfer operator. Put q=L-1=4x+O(x^2).
The vectors 1,q,...,q^r form a basis: relative to 1,x,...,x^r their matrix is
triangular with diagonal 1,4,...,4^r, all nonzero. Multiplication by q sends
q^j to q^(j+1) and q^r to zero. Therefore M_L is unipotent with exactly one
Jordan block of size r+1; its minimal polynomial is (t-1)^(r+1).

Consequently every output sequence from the stabilized carrier satisfies
(r+1)-st forward difference zero. For the depth-r excess output f_r(h), the
proved nonzero leading term makes this order sharp: the generating function
sum_(m>=0) f_r(h_0+m) z^m has reduced denominator exactly (1-z)^(r+1).
The finite exceptional prefix before h_0 only changes its numerator after
combining the full generating function. Thus fixed excess admits a genuinely
bounded, operation-closed carrier, even while the unrestricted composite
inventory has no such bound supplied here. This statement is over Q; reducing
modulo 2 annihilates the crucial coefficient 4 and does not preserve the
single-block proof.

### A pointed-tournament expression for the same kernel

A second exact expression connects the depth transfer to vertex symmetry:

    K(x) = 1 - B(x)^2 / D(x),
    D(x) = (partial Z_T / partial p_1)|_(p_i=x^i).           (6)

Derivation: before scalar specialization, differentiate
A=(p_1+P circ A)/(1-p_1-P circ A) with respect to p_1. In a plethystic
substitution, psi_d(A) has no p_1 dependence for d>1. Consequently

    A_1 = (1+A)^2 (1+(P_1 circ A) A_1).

Since A_1=(Z_T)_1, rearrangement and specialization give (6).
D(x) is the ordinary generating series of derivative species T': its x^n
coefficient is the number of isomorphism classes of tournaments on n+1 vertices
with a distinguished vertex, equivalently the sum of vertex-orbit counts over
unpointed order-(n+1) tournaments. It is *not* B'(x), which weights each
unpointed type by n+1 instead of its number of vertex orbits.

For an odd partition lambda of n, let b_lambda=2^e(lambda)/z_lambda be the
classical Burnside coefficient. Differentiating in p_1 gives the independent
formula

    [x^n]D(x)=sum_(lambda odd, |lambda|=n) 2^ell(lambda) b_lambda.

Adding the distinguished fixed vertex creates one new pair-orbit with every
existing vertex cycle, which explains 2^ell(lambda). This formula computes
K directly from fixed-permutation data and independently checks the prime
substitution derivative. The diagonal checker matches all coefficients through
r=6, with L coefficients 1,4,23,186,2332,47020,1590832. This is a precise bridge
between aggregate vertex symmetry and the long-depth composite transfer law.

One subtle distinction: at order six, C3 followed by C3 has two prime nodes but
depth one. The exact b,l,h census has one (2,1,1) type and seven (2,1,2) types.
So prime-node count does not recover prime-depth.

### Extreme node density: branching versus depth

There is another infinite family generated by the single prime C3. Let C_b
count tournaments with b prime nodes on the minimum 2b+1 vertices. By the tree
identity every internal node is a C3 node and there are no linear nodes.
These canonical cyclic ternary trees have ordinary series

    C(z)=sum_(b>=0) C_b z^b,
    C(z)=1+z(C(z)^3+2C(z^3))/3.                             (7)

The z tracks prime nodes; the three repeated children in a quotient 3-cycle
explain z^3 inside C, another instance of weighted Adams. The first coefficients
are 1,1,1,2,5,11,28,74. Thus maximum prime-node density does not force maximum
prime-depth: the density family permits branching; only its chain has depth b.

With exactly one extra vertex, n=2b+2, absence of prime arity four forces one
binary linear node and otherwise C3 nodes. Let E_b count these tournaments.
The binary linear root gives C^2; a C3 root with one exceptional child containing
the linear node gives z C^2 E. A nonidentity quotient rotation cannot fix such
an assignment because it would replicate the exceptional node three times.
Consequently

    E(z)=C(z)^2/(1-z C(z)^2).                              (8)

Its first coefficients are 1,3,8,23,69,208,636. A supplemental calculation
computed both series through b=20 and compared the relevant coefficients with
the full marked cycle-index output through order15. The published
`check_density_symmetry.py` computes them through b=15 and compares the first
eight one-linear-node coefficients with independently typed-tree-enumerated
coefficients. These are aggregate family recursions, not a classified list of
representatives or automorphism groups. No literature novelty is claimed.

### Full automorphism-order distribution inside the C3 density family

This restricted family admits more information than a general averaged prime
bank: the sole prime C3 has a known actual action. Let u mark the integer
log_3|Aut(T)| and z mark prime nodes, including the singleton as the constant
term. The bivariate family series obeys

    C(z,u)=1+z[(C(z,u)^3-C(z^3,u^3))/3
                    +u C(z^3,u^3)].                      (9)

Proof: for an ordered triple of child types, its C3 stabilizer is trivial unless
all three child types are identical. Nonconstant triples have orbits of size
three, producing (C^3-C(z^3,u^3))/3. A constant triple has one orbit, gains a root
C3 factor, and triples its child's automorphism exponent, producing
u C(z^3,u^3). Child automorphism orders multiply. The diagonal correction refers
to *identical child isotypes*, not merely equal sizes or equal group orders.
That exact diagonal sum is C(z^3,u^3), which is why the aggregate generating
series suffices here. Setting u=1 recovers (7).

All groups in the family have order a power of three. More precisely,
log_3|Aut(T)| equals the number of prime-tree nodes whose three children have
the same isomorphism type. Every nontrivial finite tree has a bottom prime
node, and its three singleton children are identical. Hence no nontrivial
member is rigid. The unique chain with b>=1 prime nodes has only its bottom
node symmetric and therefore has Aut(T) isomorphic to C3 for every b.
Conversely any non-chain prime tree has at least two bottom internal nodes,
so its automorphism order is at least9. Thus, in this restricted family,
Aut order3 is equivalent to the maximum-depth chain.

The bound |Aut(T)|<=3^b is also immediate. Equality forces every node to have
three identical child trees, recursively forcing a complete ternary tree of
some height h. Thus equality is possible exactly at b=(3^h-1)/2, with one
isotype, and then |Aut(T)|=3^b. This is only a density-family claim, not an
all-tournament extremal-group theorem.

`check_density_symmetry.py` computes (9) through b=15, independently generates
canonical cyclic ternary trees through b=9 and compares every group-order
coefficient, then constructs their actual tournament adjacency matrices and
enumerates every vertex permutation through b=3 (n<=7). These actual group
orders are 3 at n3, 3 at n5, and 3 and9 for the two n7 types. Normal and -O
receipts are identical. Equation (9) does not contradict the general loss
ledger: it uses the explicit known action of the one allowed prime, and it
restricts to a family where the entire stabilizer test is equality of all
three children. A general prime palette need not admit such a scalar test.

## 5. Loss ledger and limits

Source -> target: the aggregate P goes by canonical weighted substitution to
aggregate A and its intrinsic tree-statistic strata. Preserved: cycle-type
averages, labelled fixed-point counts, ordinary isotype counts, marked arity
profiles, prime-depth counts. Destroyed: pairing of index summands with named
prime representatives, actual automorphism groups, individual composite
representatives, attachment/edge-flip incidence, anti-automorphisms.

This is an affirmative sufficiency statement for these aggregate outputs, not
an information-theoretic impossibility theorem about every invariant of the
one fixed tournament sequence. To output concrete composite representatives
one needs concrete prime representatives and their actions on slots; to output
stabilizer groups one needs those groups/actions plus child-type assignments.
The averaged cycle index alone is not a list of orbit stabilizers. In particular
its identity coefficient gives the sum of reciprocal group orders, not the
histogram of group orders.

## 6. Exact evidence and reproduction

Run with the Python standard library:

    python3 check_marked_composites.py 9 > check9.json
    python3 -O check_marked_composites.py 9 > check9-optimized.json
    cmp check9.json check9-optimized.json
    python3 check_marked_composites.py 12 > check12.json
    python3 check_marked_composites.py 9 --full > full-marked9.json
    python3 check_depth_diagonals.py
    python3 check_density_symmetry.py

The checker independently recovers the prime input from the odd-partition
Burnside formula, discards B as a construction input, and runs the forward
marked equation. It compares the entire unmarked output with B, not only counts.
An independent exhaustive oracle enumerates all labelled orientations through
six, quotients by all permutations, finds *all* modules directly, tests strong
modules by non-overlap, takes maximal proper strong modules, and extracts the
canonical tree. Every marked cycle-index coefficient and every depth cumulative
count agrees through six. It does not depend on a stored receipt or an unchecked
classification routine. All checks use RuntimeError, remaining active with -O.

[The separate tables](TABLES_THROUGH15.md) give prime-node, linear-node, and exact-depth marginals through order 15; the JSON contains joint node counts. `--full` exports every marked power-sum coefficient. Normal and optimized order-12 and order-15 receipts agree byte-for-byte. The order-15 run takes roughly eight minutes per mode here, versus about half a second for the default order nine; the small default is intentionally suitable for routine publication checks.

Through order nine, total isotype counts by maximum allowed prime-depth are:

| depth bound | n=5 | n=6 | n=7 | n=8 | n=9 |
|---|---:|---:|---:|---:|---:|
| 0 | 1 | 1 | 1 | 1 | 1 |
| 1 | 11 | 49 | 411 | 6535 | 187795 |
| 2 | 12 | 56 | 455 | 6869 | 191440 |
| 3 | 12 | 56 | 456 | 6880 | 191535 |
| 4 | 12 | 56 | 456 | 6880 | 191536 |

All rows beyond the exhaustive oracle range are exact consequences of the
proved recurrence, not brute-force censuses. The script also checks rational
coefficient positivity, integral marked isotype counts, the tree support bound,
and the deliberately wrong-Adams hostile.

## Structural sources and attribution boundary

The canonical decomposition input is Gallai's theorem as recalled in Theorem 16
of [Boussaïri–Chergui–Ille–Zaidi, 2018](https://arxiv.org/html/1805.04607v1).
The ordinary species substitution rule is Theorem 2.5 of
[Gainer-Dewar–Gessel, 2014](https://emis.de/ft/31891).
These primary sources were opened during this session. The weighted marking,
forward filtration, and extremal-depth formulas above are proved here as
consequences; no claim they are new in the literature is made. The supplied
small prime-bank values are independently rechecked by the script's census.

## Independent audit

An independent audit validates the proofs and performs an independent sparse scalar check through h=12,r=4, plus a Burnside-pointed kernel check through r=6. It confirms the hostile premature-transfer mismatch 45 versus44 at h=2,r=2 and the pointed order-three count4 versus the incorrect vertex-size weight6. Normal and optimized audit outputs agree.


The published density checker computes the one-linear-node recurrence through
b=15 and compares its first eight coefficients with independently typed-tree-enumerated
coefficients; it does not rerun that separate one-linear-node oracle. The published diagonal checker covers h<=15 and excess r<=6; the separate
h=20 thin-boundary check is supplemental evidence, not a hidden input. All public
checkers are self-contained and read no receipt files.
