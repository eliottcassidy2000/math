# Fixed-apex reset algebras and periodic constructions with polynomial counters

**PROVED, with standalone exact certificates. Session 192, 2026-11-05.**

A genuinely different use of the repository's exterior-transfer and Lie-algebra work is to keep the entire attachment operator, rather than guess a scalar update for H. The result is an explicit common tensor flag for every fixed-apex attachment, an all-order characteristic polynomial for every mixed block, and a constructive periodic-word calculus. Mixing letters creates polynomial counters that a single constant run cannot produce. Two-apex word `120` attains all three possible positive pole orders and starts an infinite strong family at H=13.

These results do **not** prove the global `{7,21}` spectrum conjecture, an interval of attainable values, or a recursively reusable ternary digit grammar. No literature-priority claim is made.

## 1. Inheritance, sources, and why this is a new operation

Closest mechanism: Session 191's [two-apex grammar](../../session191_tournament_growth_20261104/two_apex/TWO_APEX.md) gives the 36-state nonnegative two-apex transfer and proves that a **single repeated letter** has bases 1,2,3 with no polynomial factor. Its [sparsity theorem](../../session191_tournament_growth_20261104/two_apex/SPARSITY.md) excludes a fixed finite collection of constant-run templates as a universal construction. Here the repeated object is an arbitrary nonconstant **block**, so the number of letter changes grows without bound.

Canonical hostile: Session 191's H=13 example is outside the entire one-apex class, but the two-apex word `120` realizes it. Corrected near miss: the eigenvalues of each generator do not normally determine the eigenvalues of their products. We prove the common flag that makes this inference valid here, and retain Jordan data that eigenvalues omit.

Two less recently used repository sources supplied concrete questions:

- `01-canon/theorems/THM-3186-full-exterior-continuant-path-convolution-and-cancellation-wall.md`: retain invariant-line extensions and their transfer amplitudes; eigenvalues alone can forget the projected observable. Here the actual transfer is derived from Hamilton-path run assignments, not imported from factorial moments.
- `01-canon/theorems/THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge.md`: test a possible constructive Stern shear mechanism at the operator level. The exact common flag gives a decisive boundary for this candidate, while producing useful Jordan counters instead.
- `session35-lie/RESEARCH_NOTE.md`: dimensions do not specify a bracket. The Lie statement below uses an actual common representation and commutator, not a matching dimension.

Live concept board: endpoint resets; tensor flags; positive integer spectra; nilpotent extension data; periodic strong constructions; Stern shear hostiles.

## 2. Arbitrarily many fixed apexes

Fix a tournament A on r>=1 named apexes, and a word w over `{0,1}^r`. Add a transitive chain x_1 -> ... -> x_m. The i-th bit of its letter says that apex i points to that chain vertex. Denote this actual tournament T_A(w).

For each of the r! possible orders p_1,...,p_r of the apexes in a Hamilton path, distribute each chain name among the r+1 runs

    L_0, p_1, L_1, p_2, ..., p_r, L_r.

Each run is increasing. For each run L_i, i<r, retain a three-dimensional state space with basis E,G,B: empty; nonempty with its last name pointing to the following apex; nonempty with the opposite last edge. The last run has basis E,O: empty/nonempty. This gives one space

    V = (Q^3)^(tensor r) tensor Q^2

of dimension 2*3^r per apex order. A bad last edge must be retained because a later chain name can repair it.

For a current letter c, assigning its chain name to run i is a rank-one operator R_i(c) on that factor, tensor the identities on all other factors. In a bounded run it sends G and B to the same current good/bad target. It sends E to that target if the incoming apex edge is permitted, and otherwise to zero. L_0 has no incoming condition. In the final run it sends O to O and E to either O or zero according to the incoming apex edge.

Thus the letter transfer is the concrete tensor sum

    M_c = R_0(c) + ... + R_r(c).                       (1)

Operators on different factors commute. Each nonzero local reset is idempotent. First-name constraints are imposed at birth; last-name constraints and consecutive-apex edges across empty runs are imposed by the accepting covector. These are precisely all adjacency conditions in the displayed Hamilton word. Summing over the r! apex orders therefore gives H(T_A(w)) exactly. This is a bijection, not a fitted recurrence.

## 3. One common tensor flag

In each bounded factor choose the ordered basis

    E, B, D=G-B,                                      (2)

and in the final factor choose E,O. Give them levels 0,1,2 and 0,1 respectively. Every reset has the following form in (2):

    E -> alpha*(B + delta D),
    B -> B + delta D,
    D -> 0,                                          (3)

where alpha is its first-name permission and delta says its current last-name edge is good. In the final factor,

    E -> alpha O,    O -> O.

All matrix entries move a tensor coordinate only upward in these orders. The diagonal contribution of any local reset is 1 exactly on the coordinate B, or O in the last factor, and 0 otherwise. Consequently the diagonal of M_c at a tensor state s is

    j(s) = number of coordinates occupying their level-1 basis vector. (4)

Crucially, (4) is independent of c. This proves simultaneous triangularizability explicitly over Q; no Lie theorem or algebraic closure is needed.

### Theorem 1: every mixed block has the same spectrum at fixed length

Let v be any nonempty block of length L, and let B_v=M_(v_L)...M_(v_1). On one apex-order space,

    char(B_v;X) = product_(j=0)^(r+1) (X-j^L)^m_(r,j),
    m_(r,j) = [z^j](2+z)^r(1+z).                     (5)

On the full r!-order transfer space each multiplicity is multiplied by r!. The statement is independent of the apex tournament A as well as the letters of v. A enters the accepting covector, not the transfer matrices.

**Proof.** Each diagonal entry multiplies to j(s)^L. There are two zero-diagonal choices and one unit-diagonal choice in each bounded factor, and one of each in the final factor. Their multiplicity generating polynomial is exactly (2+z)^r(1+z). Taking the r! direct summands proves the full statement.

For r=2 this specializes to

    X^8 (X-1)^16 (X-2^L)^10 (X-3^L)^2.               (6)

This is a constructive transfer theorem, not a claim that the actual tournament adjacency matrices have these eigenvalues.

## 4. Exact Jordan bounds and the periodic construction theorem

The common tensor flag has more information than a bare triangular matrix. Let rank(s) be the sum of its coordinate levels. If j(s)=j>=1, then

    j <= rank(s) <= 2r+2-j,
    rank(s) = j + 2*(number of D factors).             (7)

Hence an increasing tensor-coordinate path can visit at most r+2-j **distinct** states having the same positive diagonal j. For j=0 it can visit at most r+1 states: the final coordinate is empty and only the r bounded coordinates can move E->D.

For completeness, this bounds Jordan blocks without assuming that intermediate states have equal diagonal. Write B_v=D+N with D diagonal and N strictly triangular in any linear extension of the tensor order. Expand its resolvent as the finite series

    (zI-B_v)^(-1)
      = sum_(a>=0) ((zI-D)^(-1)N)^a (zI-D)^(-1).      (8)

Each term is a sum over increasing tensor-coordinate paths. A pole at j^L can occur once for each visited state with that diagonal. The bounds just established therefore bound the pole order in every matrix entry. Equivalently,

    B_v^(r+1) product_(j=1)^(r+1)
      (B_v-j^L I)^(r+2-j) = 0.                      (9)

### Theorem 2: arbitrary repeated blocks have controlled polynomial counters

Fix any prefix u, nonempty block v of length L, and suffix w. There are rational polynomials P_1,...,P_(r+1) with

    H(T_A(u v^k w))
      = sum_(j=1)^(r+1) P_j(k) j^(Lk),               (10)
    deg P_j <= r+1-j,

for every k>=r+1. The possible finite zero-eigenvalue transient is confined to smaller k. Boundary cases can be checked by the same finite transfer; no assertion that they always vanish is needed.

**Proof.** Equation (9) bounds each generalized positive eigenspace, so powers contribute a polynomial of degree at most r+1-j times (j^L)^k. The zero generalized eigenspace vanishes by k>=r+1. Fixed prefix and suffix only change the initial vector and final covector. All matrices and distinct eigenvalues are rational, so the polynomial coefficients are rational.

For two apexes this is

    P_1(k) + P_2(k) 2^(Lk) + C 3^(Lk),
    deg P_1<=2, deg P_2<=1.                          (11)

A constant-letter block can have much smaller Jordan blocks, exactly as in S191. Replacing it by a mixed block is an actual way to turn on additional counters; it does not introduce new exponential bases.

## 5. Three explicit all-order strong families

Use the S191 convention a->b, with bits 0 and 1 of each letter recording a->x and b->x. The following rational generating functions start at k=0 with the empty-chain tournament, whose H is 1. Every tournament at k>=1 is strong.

### A. A linear counter and a 4^k clock

For w_k=(20)^k,

    H(T(w_k)) = (10*4^k + 6k - 1)/9,                 (12)
    sum_(k>=0) H(T(w_k)) t^k
      = (1-t-2t^2)/((1-t)^2(1-4t)).

The first values are 1,5,19,73,287,1141,... . This is a one-apex subfamily inside the two-apex notation; no new feedback-number claim is made. It isolates the polynomial counter in the simplest form.

### B. A mixed exponential linear counter

For w_k=(12)^k,

    sum_(k>=0) H(T(w_k)) t^k
      = (1-14t+57t^2-76t^3)
        /((1-t)^2(1-4t)^2(1-9t)),                   (13)

or explicitly,

    H(T(w_k)) = 19/54 + 4k/9
              + (316/675 + 8k/45)*4^k
              + (9/50)*9^k.

The first values are 1,5,29,197,1485,12021,101949,... . This realizes a nonzero k*4^k term, which cannot occur in any fixed-prefix/fixed-suffix **constant-letter** clock from S191.

### C. Simultaneous sharpness of every two-apex Jordan bound

For w_k=(120)^k,

    sum_(k>=0) H(T(w_k)) t^k
      = (1-33t+255t^2-615t^3-284t^4)
        /((1-t)^3(1-8t)^2(1-27t)).                  (14)

The first values are

    1,13,225,4917,121381,3166885,84462961,2270871133,... .

The numerator is nonzero at 1,1/8,1/27. Thus the denominator is reduced and all three bounds in (11) are attained: a genuine quadratic counter, a genuine k*8^k term, and a nonzero 27^k term. For example the coefficient of k^2 is 13/49. The first member is exactly the H=13 hostile that crosses the all-order one-apex boundary.

**All-order certification of (12)-(14).** Let q(X) be the product of the displayed positive characteristic factors, with their indicated powers. The verifier constructs each exact 18-dimensional block in both apex orders and applies q(B_v) to the initial empty-state vector. In each case the resulting vector is killed by B_v. Their total accepting value is also zero. Therefore the scalar q recurrence holds at k=0 and at every later k, not merely for a finite prefix. The finitely computed initial values determine the numerators uniquely. These are exact matrix identities constituting finite certificates of all-order formulas.

**Strongness.** Each nonempty displayed word starts with a nonzero letter, ends with a letter other than 3, and contains a 2. The letter 2 supplies b->x->a, putting both apexes in one strong component. That component reaches the first chain name, and the last chain name reaches it. The increasing chain then lies entirely in that component.

## 6. What the algebra rules out, and what it supplies instead

The Lie algebra generated by all letter matrices is solvable: every element is triangular in the explicit common basis, and every commutator has zero diagonal. This is an actual bracket statement, not an analogy based on a dimension.

It rules out a specific attractive Stern construction. The usual rational-enumeration shears

    U = [[1,1],[0,1]],    V = [[1,0],[1,1]]

have UV eigenvalues (3+sqrt(5))/2 and (3-sqrt(5))/2. No two fixed attachment blocks, even of different positive lengths, can induce these two actions on a common invariant two-dimensional subquotient: their concatenation has only eigenvalues j^(L_U+L_V), all integers, by Theorem 1. This excludes a faithful standard-Stern linear carrier inside the fixed-apex transfer. It does **not** exclude nonlinear target coding, noninvariant coordinates, changing the apex set, or an unbounded-state construction.

The productive substitute is explicit: the nilpotent extensions, rather than a semisimple SL2 action, provide polynomial counters. Equations (12)-(14) are the first three physical output families. Theorem 2 gives an algorithmic compiler for any chosen periodic block: build its finite transfer, obtain the bounded recurrence, compute enough initials, and certify the boundary transient. Searching mixed blocks for a prescribed numerator or nilpotent coupling is a genuinely different design problem from fitting single constant-run slopes.

## 7. Loss ledger and exact scope

- Source: a fixed-apex tournament presentation with all boundary states.
- Map: Hamilton run assignment -> tensor reset operator -> common triangular flag.
- Preserved: exact H when initial/accepting vectors are kept, mixed-block characteristic polynomial, and Jordan bounds.
- Lost by taking only the spectrum: essentially all attachment order and every scalar path-count distinction. The length-two blocks `00` and `12` have identical full characteristic polynomial (6), but H=1 and H=5.
- Needed sidecar: actual off-diagonal transfer entries and the initial/accepting pair; (14) requires the nonzero pole numerators.
- Cheapest decisive test: compare those two blocks, and then compare the simple-pole single-letter law with the double/triple poles in (13)-(14).

No density or arbitrary-target theorem follows from these families. A fixed periodic family remains sparse; an unbounded catalogue of blocks is a different problem. The global `{7,21}` conjecture, constructive tail coverage, and a universal closed marked grammar remain open.

## 8. Reproduction

Run `python3 verify.py` and `python3 -O verify.py`. The script is standalone, standard-library-only, deterministic, and uses explicit exception gates rather than removable asserts. It checks:

- 30 local exact tensor-conjugacy identities for r=1,2,3,4;
- every block through length 4 at r=1, length 3 at r=2, and length 2 at r=3;
- the entire matrix annihilator (9) on every basis vector for the first two of those block universes;
- 341 two-apex words independently by full vertex-subset DP;
- all-order residual certificates for the three generating functions;
- 12 actual strong-family tournaments by vertex DP, including the order-14 member of (14);
- all three nonvanishing pole numerators in the sharp family.

Normal and optimized runs produce identical stdout and byte-identical `receipt.json`. The general theorems rest on the explicit flag and resolvent proofs, not on the bounded verification universe.
