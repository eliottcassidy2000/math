# Colored score-energy and rooted-forest laws on the fixed-path cube

Status: PROVED by the elementary derivations below; FINITE-EXACT controls are
in `verification.json` and `wiring_verification.json`. A separate independent
review is not claimed. No literature-priority claim. No Git mutations were
made in this lane.

## 1. Definitions and exact agreement with repository conventions

Fix the directed path `1 -> 2 -> ... -> n`, with `n>=2`. The free pairs are
`E={ (i,j): 1<=i<j<=n, j-i>=2 }`, of cardinality

    D=(n-1)(n-2)/2.

Let `r(i)=n+1-i`. Reflection sends the forward bit at `(i,j)` to the forward
bit at `(r(j),r(i))`. This is converse followed by label reversal, with no
bit negation. It is exactly sigma in
`01-canon/theorems/THM-643-blue-black-line-parity-structure.md`, not a new
interpretation of blue. A presentation is blue precisely when it is fixed
by this reflection.

The number of fixed free edges and reflection orbits are

    f=floor((n-1)/2),
    o=(D+f)/2=floor((n-1)^2/4).

For `n=2m`, the reflection-fixed pair `(m,m+1)` is a path edge, leaving
`m-1` fixed free edges. For `n=2m+1`, all `m` pairs crossing the central
vertex are free. Thus the total and blue presentation counts are

    M=2^D,                 B=2^o.

The all-free-edge complement C commutes with reflection. For `n>=3`, C
and C followed by reflection have no fixed presentations: the latter
would toggle the reflection-fixed pair `(1,n)`. Consequently:

    blue complement lines                    B/2,
    literal black complement lines           (M-B)/2,
    fundamental paired-black edges           (M-B)/4.

The last division by two is essential after pairing the reflected black
lines. It is the convention of THM-849 section 2 and THM-4731 statement 3,
respectively in `THM-849-n8-black-self-line-obstruction.md` and
`THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md`.
At n=2, C is the identity and the line formulas do not apply. At n=3,
every presentation is blue, so there is no conditional black distribution.

## 2. Exact colored score-energy moments

Let A be the tournament adjacency matrix, let

    d=(A-A^T)1,        W=||d||^2,
    mu=n^2-3n+4,
    delta=0 if n is even, and n-3 if n is odd.

The uniform all-presentation and blue-presentation expectations are

    E_all W = mu,
    E_blue W = mu-delta.

Equivalently the blue mean is `n^2-3n+4` for even n and `n^2-4n+7` for
odd n. For `n>=4`, subtraction gives

    E_black W = mu + B delta/(M-B).

These are means over presentations, not means over unweighted tournament
isomorphism classes.

Proof. Put `b=e_1-e_n` and, for each free pair e=(i,j), put
`v_e=e_i-e_j`. Write its forward sign as `x_e in {-1,1}`. Then

    d=b+sum_e x_e v_e.

For independent free signs, cross terms average to zero and
`E W=||b||^2+sum_e ||v_e||^2=2+2D=mu`.

On the blue cube there is one independent sign x_O per reflection orbit,
so `E_blue W=2+sum_O ||sum_(e in O) v_e||^2`. A fixed edge contributes
2. A two-edge orbit with disjoint endpoints contributes 4. Only an odd
n can give an orbit of two distinct edges sharing a vertex: they are
`(i,c),(c,r(i))`, with `c=(n+1)/2` and `i<c-1`. Their incidence vectors
have inner product -1, so each such orbit contributes 2 rather than 4.
There are `(n-3)/2` of them. The total decrease is n-3. This is the exact
source of the parity split; all other orbit covariances are silent.

## 3. Degree-two forest coefficients and their old cyclic-triangle meaning

Use the outdegree Laplacian `L=diag(s)-A`, and let

    F_T(t)=det(I+tL),                 N=binom(n,2).

The directed matrix-forest interpretation is the one in the companion
Hamiltonian-size-biased forest note: every nonroot has one outgoing arc,
and every component points toward its root. Since tournaments have no
oppositely oriented pair,

    F_2(T)=[t^2]F_T=(N^2-sum_i s_i^2)/2
          =(4N^2-n(n-1)^2-W)/8.

Therefore, setting

    q_n=(n-2)(n^2(n-1)+2)/8,

we obtain

    E_all F_2=q_n,
    E_blue F_2=q_n+delta/8,
    E_black F_2=q_n-B delta/(8(M-B))       (n>=4).

This statistic is not a new individual tournament invariant. If c3(T)
counts cyclic triangles, the familiar source-count argument gives

    F_2(T)=c3(T)+binom(N,2)-binom(n,3).

Thus the degree-two forest coordinate is exactly the old cyclic-triangle
axis plus a constant. The content here is its exact color-resolved
presentation-weighted mass. Higher forest coefficients can retain more,
as the explicit equal-energy nonloop in section 7 demonstrates.

## 4. Exact blue forest polynomial as a compatible-orbit sum

Let F range over all spanning directed rooted forests on [n]. Declare it
compatible if:

1. it uses no backward path arc `i+1 -> i`;
2. whenever it uses free edges from one reflection orbit, their signs
   relative to increasing endpoint order agree.

Let q(F) be the number of free reflection orbits touched by F. Then

    E_blue F_T(t) = sum_(compatible F) t^|F| 2^(-q(F)).       (1)

Equivalently, its total over all B blue presentations is

    sum_(blue T) F_T(t) = sum_(compatible F) t^|F| 2^(o-q(F)). (2)

Proof. A compatible forest pins one bit in each touched orbit; its
probability is exactly `2^(-q(F))`. An inconsistent demand pins one orbit
to both signs and has probability zero. Sum the matrix-forest expansion.
The ordinary all-presentation counterpart weights a path-compatible forest
by `2^(-r(F))`, where r(F) is its number of free edges, without the orbit
consistency requirement.

This is an exact finite combinatorial formula, not a claimed elementary
closed product. It retains the missing correlation coordinate. In
particular one cannot replace the blue arc weights by their marginal
means and then take a determinant: all free marginal means remain 1/2,
but at n=5 the resulting answer is already wrong at t^2 by 1/4.

For comparison, the previously derived all-presentation mean is

    P_n(t)=sum_(l=0)^floor((n-1)/2)
      binom(n-1-l,l) (1+nt/2)^(n-1-2l) (t^2/4)^l.

Writing Q_n(t) for (1), the exact conditional black mean is

    (M P_n(t)-B Q_n(t))/(M-B),                 n>=4.

The sign of Q_n-P_n coefficientwise is not established at all n. Small
positive examples are evidence, not a universal positivity theorem.

## 5. Reflection gives one additional exact coefficient

Let J be the all-ones matrix, R the vertex-reversal permutation matrix,
and `K=2L-nI+J=diag(d)-(A-A^T)`. For every blue tournament,

    K1=0,          RKR=-K.

Thus K descends to the quotient by the constant vector and is similar
there to its negative. If `a=1+nt/2`, the forest polynomial has the form

    F_T(t)=a^(n-1)+c_1(T)a^(n-3)t^2
             +c_2(T)a^(n-5)t^4+...,

where terms stop when their a-exponent is nonnegative, and

    c_1(T)=(n(n-1)-W(T))/8.

Indeed `tr K^2=W-n(n-1)`, and the second characteristic coefficient is
`-tr K^2/2`; the factor t/2 supplies the additional factor 1/4. This
argument uses characteristic polynomials, not diagonalizability.
Consequently, for `n>=4`,

    [t^3](Q_n-P_n)=delta*n*(n-3)/16.

At even n both the t^2 and t^3 means agree, but equality need not persist:

    Q_5-P_5 = t^2/4 + 5t^3/4 + 27t^4/16,
    Q_6-P_6 = (5/16)t^4(1+3t).

So equal colored score-energy means do not imply equal forest laws. The
full n=2..8 polynomials are recorded in verification.json.

## 6. Intrinsic merged-node inventory constraints

For an ordinary tournament class, h=H/|Aut| is its fixed-path fiber size.
Give a self-converse merged node mass w=h, and a distinct converse pair
mass w=2h. Let b be its blue fiber size (zero on a non-self-converse pair).
Its blue degree is b, and its fundamental paired-black degree is

    k=(w-b)/2,

with loops counted twice. Write Fbar=F on a self-converse node and
`Fbar=(F_T+F_(T^op))/2` on a distinct converse pair. This averaging is
required for the full polynomial; it is unnecessary for W and F_2,
which are converse invariant. Then

    sum_nodes b Fbar = B Q_n,
    sum_nodes k Fbar = (M P_n-B Q_n)/2,

and in particular

    sum_nodes b W = B(mu-delta),
    sum_nodes k W = ((M-B)mu+B delta)/2,
    sum_nodes b F_2 = B(q_n+delta/8),
    sum_nodes k F_2 = ((M-B)q_n-B delta/8)/2.

These are necessary colored-incidence inventory constraints. They are not
sufficient realizability conditions and do not reconstruct the endpoint
matching. The following inherited wiring laws retain an additional
coordinate rather than being inferred from these averages.

## 7. Inherited wiring laws, and a minimal forest refinement witness

The exact endpoint flux and odd/even parity laws below are already proved
in THM-785 `cyclic-triangle-flow-on-blue-black-metagraph.md`, sections 2-3,
and THM-790 `blue-parity-law-proved.md`, sections 1-4. The no-odd-blue-loop
and even endpoint necessary condition are THM-648
`blue-selfloops-only-at-even-n.md`. They are reproduced here to identify
the closest inherited mechanism, not promoted as new results.

For any fixed-path presentation, using the decomposition in section 2,

    d(T)=b+X,     d(CT)=b-X,
    W(T)-W(CT)=4(X_1-X_n)=8(s_1-s_n-1),
    c3(CT)-c3(T)=s_1-s_n-1.

For blue presentations, `s_1+s_n=n-1`, hence

    c3(CT)-c3(T)=2s_1-n.

Therefore the odd-order blue graph is canonically bipartite by cyclic-
triangle parity and has no loops, even after converse merging. At even
order blue edges preserve cyclic-triangle parity, and any blue loop must
satisfy `s_1=n/2`, `s_n=n/2-1`. This is necessary, not sufficient. The
score-lattice version `(W-epsilon*n)/8 mod 2`, with epsilon=1 for even n
and 0 for odd n, is the same parity up to an n-dependent constant, not a
new invariant or an asserted Arf invariant.

A minimal even-order hostile is n=6. Its outgoing neighbors, using labels
1,...,6, are

    1: 2,3,4       2: 3           3: 4,6
    4: 2,5,6       5: 1,2,3,6     6: 1,2.

This presentation is blue. Its free-orbit bits, in zero-based coordinates
and the order below, are `1,1,0,0,0,0`:

    {(0,2),(3,5)}, {(0,3),(2,5)}, {(0,4),(1,5)},
    {(0,5)}, {(1,3),(2,4)}, {(1,4)}.

It and its free-edge complement have the same score multiset
`{1,2,2,3,3,4}`, W=22, c3=6, and the same first three forest coefficients.
Nevertheless,

    F_T  = 1+15t+91t^2+279t^3+433t^4+273t^5,
    F_CT = 1+15t+91t^2+279t^3+432t^4+270t^5.

The distinct forest polynomials certify nonisomorphism. Both tournaments
are self-converse, so they remain different after converse merging. This
is a concrete place where the full forest coordinate adds information
lost by its degree-two shadow. Direct permutation checking independently
confirms nonisomorphism. It is minimal: odd orders have no zero blue flux;
at n=4, the only two zero-flux blue presentations are loop endpoints. At
n=6 there are 24 zero-flux presentations but only four loop endpoints,
matching the inherited THM-648 count without using it as an input.

## 8. Reproduction, controls, and stopping boundary

Run:

    python3 session31-blue-forests/check_blue_forests.py
    python3 session31-blue-forests/check_wiring.py

The scripts use only the Python standard library. The first exhaustively
checks all-presentation energy means through n=7, blue means through n=9,
full blue forest polynomials through n=8, and full all-presentation
polynomials through n=6. An independent functional-parent enumeration of
all rooted directed forests through n=6 verifies the compatible-orbit
formula coefficient by coefficient, rather than reusing the determinant
calculation. Bareiss determinants at t=-2,1,3 check the Newton-identity
polynomial routine on selected hostile and endpoint masks. Klein orbits
are directly counted through n=7.

The second checks exact endpoint, triangle, energy, and parity fluxes on
all presentations through n=7 and all blue presentations through n=9,
then independently checks the n=4,6 zero-flux/loop pools by permutations.
The adjacency and polynomial witness is persisted in its JSON receipt.

Scope boundary: no new closed product for the full blue forest mean is
claimed; no full blue-graph wiring determination follows from these
moments; no universal forest coefficient monotonicity is inferred from
the finite rows; the low-order forest statistic is explicitly identified
with the old triangle axis. The successful new bridge retains reflected
edge-orbit compatibility and supplies exact color-resolved mass checks,
plus a higher-forest separator of a tied blue line.
