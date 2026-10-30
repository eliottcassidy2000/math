# Signed-graph equations: the complete integral obstruction

**PROVED elementary extension of the inherited homogeneous signed-incidence
mechanism.** The earlier [odd-order arithmetic note, §4](../../odd_prime_function_cycle_20261018d.md)
already gives unsigned homogeneous kernels. This note retains inhomogeneous
edge data, arbitrary signs and the integer cokernel. No priority claim is made.

Let a finite graph have n vertices and e edges; parallel edges and loops are
allowed. Give each edge an orientation v->w and sign s_e in {+1,-1}. Its equation is

    x_v-s_e x_w=b_e.

The matrix B has one such row per edge. A positive loop gives a zero row; a
negative loop gives 2x_v=b_e. A cycle is balanced when its sign product is +1.
A connected component is balanced when every cycle is balanced; an isolated
vertex is balanced. Let c be the total number of components, b the number of
balanced components, and o=c-b the unbalanced count.

## 1. Complete Smith form

The nonzero Smith entries of B are

    n-c copies of 1 and o copies of 2.

There are b free kernel directions and e-n+b free cokernel directions. Thus

    ker_Z B = Z^b,
    coker_Z B = Z^(e-n+b) plus (Z/2)^o.

The statement includes the zero-row/zero-column boundaries. Over Z/M,

    |ker B| = M^b gcd(M,2)^o,
    |coker B| = M^(e-n+b) gcd(M,2)^o.                       (1)

Each nonempty inhomogeneous solution fiber has the kernel size in (1).
Cokernel cardinality is not by itself a right-hand-side solvability test.

**Proof.** In each component choose a spanning tree and a root. Tree equations
express every x_v as alpha_v t+beta_v, with alpha_v=±1 and beta_v an integer
linear combination of the tree right-hand sides. This elimination uses only
unit pivots, yielding n_component-1 unit Smith entries. Every non-tree equation
then has coefficient 0 or ±2 on t. Its coefficient is zero precisely when the
corresponding fundamental cycle is balanced. If all are zero there is a free
t-variable. Otherwise their coefficient column has gcd 2, giving one additional
Smith entry 2; all other rows vanish after integer row elimination. A connected
component has a negative cycle iff some fundamental cycle is negative, since
cycle signs multiply under symmetric difference. Summing components proves
the claims and the modular counts.

## 2. Constructive inhomogeneous solvability

After tree elimination, write each remaining equation as a_e t=d_e, where
a_e in {0,+2,-2}. Keep the actual d_e; merely knowing the cycle signs is not
enough.

Over Z:

- Every zero-coefficient equation requires d_e=0.
- In a balanced component this is sufficient, and t is free.
- In an unbalanced component choose one nonzero a_0 and normalize it to +2 by
  changing both signs. Its d_0 must be even. Every other nonzero equation,
  similarly normalized to +2, must have exactly the same right-hand side d_0.
  Then t=d_0/2 is the unique solution, and tree reconstruction gives all x_v.

Over Z/M replace equalities by congruences modulo M; require gcd(M,2)|d_0.
The unbalanced component then has gcd(M,2) solutions, not necessarily one.
These conditions are both necessary and sufficient, because they explicitly
solve every transformed row of the original matrix.

For all-negative edges, every cycle has product (-1)^length. Hence the
balanced components are exactly the ordinary bipartite components. An odd
cycle changes one free component parameter into one factor-two obstruction.
Extra cycles still impose independent compatibility on the edge data; one
odd cycle does not guarantee arbitrary b is soluble.

## 3. Hostiles and safe algorithmic use

The triangle equations x0+x1=1, x1+x2=0, x2+x0=0 have no integer solution,
but a unique rational solution and a unique solution modulo every odd M.
The even square with edge data (1,0,0,0) violates its alternating-sum condition;
it has no solution even modulo 3. Thus odd target modulus removes 2-torsion,
not the free cycle-space obstruction.

A spanning forest gives an O(n+e) arithmetic-operation algorithm to form the
constraints and reconstruct a solution. Bit complexity also depends on integer
coefficient sizes. This is an exact linear-equation solver, not an algorithm
for unrestricted graph coloring, Ramsey numbers, or nonlinear Boolean circuits.
Prescribed vertex values can be added as equations x_v=p_v and tested by the
same tree-coordinate method; they must not be discarded when passing to a cover.

The determinant/SNF statement must not be replaced by a modular Gram test:
B^T B can acquire an extra kernel in finite characteristic even when B does not.
The inherited K4 modulo3 hostile remains applicable.

## 4. Reproduction

The companion `check_signed_integral.py` independently computes determinantal
divisors by gcds of literal minors, and compares a tree-coordinate solvability
test with brute-force modular fibers. It covers all signed simple graphs on
one through four vertices, positive/negative loops, parallel-edge conflicts,
and isolated vertices: 766 matrices and 134,000 right-hand sides. Normal and
optimized Python outputs agree with `receipt.json`. This does not replace the
uniform spanning-forest proof or certify a nonlinear coloring/search problem.
