# Reusable sufficient criterion for exactly n squared HP blockers

## Theorem

Let T be a tournament of order n=4lambda+3>=11, put r=(n-1)/2=2lambda+1 and m=n(n-1)/2. Assume:

1. T is doubly regular: every indegree/outdegree is r and every distinct vertex pair has lambda common in-neighbors.
2. A subgroup G of Aut(T) acts transitively on directed arcs.
3. The bipartite support graph with tail copies, head copies, and edges u_tail--v_head for the tournament arcs u->v is connected.
4. The nonempty family of directed Hamiltonian-cycle incidence vectors has affine dimension m-2n+1 over Q. It suffices to exhibit one base cycle and m-2n+1 cycle differences with a nonzero square minor, even merely modulo a prime.
5. For some tournament arc u->v there exists a Hamiltonian path starting at u and ending at v.

Then beta_HP(T)=n-1, and the minimum arc-deletion sets destroying all Hamiltonian paths are EXACTLY:

- all incoming arcs at each unordered pair of vertices, making two sources;
- all outgoing arcs at each unordered pair of vertices, making two sinks;
- all arcs incident with each single vertex, isolating it.

These are disjoint families of sizes binomial(n,2), binomial(n,2), and n. Thus the total is n^2.

This is a sufficient criterion, not a theorem that its hypotheses hold for every Paley or every regular tournament. No integral lattice saturation is assumed in condition4.

## Proof

Every HP-blocker B meets every directed Hamiltonian cycle in at least two arcs: deleting zero or one edge from a directed spanning cycle leaves a directed spanning path. Averaging its intersection with the orbit of any cycle under the arc-transitive group gives mean |B|/r, so |B|>=2r. Isolating a vertex costs2r, hence beta=2r. For a minimum blocker the mean is exactly2; all intersections in each orbit are therefore exactly2. Starting with each cycle separately proves that EVERY Hamiltonian cycle meets B in exactly two arcs.

The connected bipartite support makes the row-column potential space b_uv=a_u-d_v have rational dimension2n-1. Every such vector is orthogonal to cycle differences. By condition4 the cycle-difference annihilator has exactly that dimension. The binary incidence vector of B therefore has a rational potential representation.

Fix one potential to0. Propagate along a support spanning tree using the integer equations a_u-d_v=b_uv. All potentials become integers; existence of the rational solution ensures consistency. This establishes integrality without claiming a unimodular cycle minor or a saturated cycle-difference lattice.

Shift all potentials by their global minimum so that their minimum is0 and maximum is M. At level k=1,...,M put S_k={u:a_u>=k}, T_k={v:d_v>=k}. Nonnegativity of b forces N^-(T_k) contained in S_k. The level cut D_k=S_k -> (V\T_k) has size r*delta_k, where delta_k=|S_k|-|T_k|. Each arc occurs in exactly a_u-d_v cuts, hence in at most one. Thus B is their disjoint union.

All active delta_k are positive integers. Regularity gives nonnegativity; a zero-deficiency compatible cut would have no edges crossing its bipartite partition in either direction, contradicting support connectivity unless both sides were full or empty. The normalized active levels are neither. Since |B|=2r, sum delta_k=2: there is either one deficiency-two cut or two nested deficiency-one cuts.

For T of size t, double regularity and Cauchy-Schwarz give

    |N^-(T)| >= r^2*t/(lambda*t+r-lambda).

For deficiency one, subtracting the desired upper bound t+1 after clearing denominators gives

    r^2*t-(t+1)*(lambda*t+r-lambda)
    =r*(r-2)+lambda*(t-1)*(2r-1-t)>0

for1<=t<=2r-1. Only t=0 or t=n-1 remain, so every deficiency-one cut is one full outgoing row or one full incoming column.

Two nested row cuts coincide and overlap; two nested column cuts do likewise. Thus the two-level alternative consists of one incoming column at w and one outgoing row at v. Disjointness requires either v=w, giving a star, or w->v. In the latter case condition5 and arc transitivity supply a spanning path from w to v. It avoids the incoming column of its first vertex and outgoing row of its last, so the mixed deletion is not a blocker.

For deficiency two, the analogous difference is

    r^2*t-(t+2)*(lambda*t+r-lambda)
    =r*(r-3)+lambda*(t-1)*(2r-2-t)>0

for1<=t<=2r-2, since r>=5. Only t=0 or t=n-2 remain, giving two outgoing rows or two incoming columns. Those Hall cuts and all stars really do destroy every HP. Their residual degree profiles distinguish the three families, completing the classification and count. QED.

## P11 verification and the P7 exception

For P11 the conditions are all checked directly: n=11,r=5,lambda=2; the square-affine group has55 elements and is arc-transitive; the support has a22-vertex spanning tree;35 listed Hamiltonian cycles provide34 independent differences (the older exact minor is-2; here it is99 mod101); and the forward endpoint path is

    0,3,4,2,5,6,7,10,8,9,1.

Thus P11 has exactly121 minimum ten-arc blockers.

At P7, r=3, the deficiency-two comparison vanishes at t=1 and4 rather than staying strictly positive. Those two equality cases contribute the seven common-out and seven common-in Hall families. The same verified rank/potential proof gives21+21+7+7+7=63. This is a separate small-order exception, not an application of the n>=11 count formula.

A complete small-order modular probe finds cycle-difference ranks8/8 for P7 and34/34 for P11 over F2/F101. Therefore the previously selected P11 determinant-2 does not signal an actual mod2 rank loss. No torsion conclusion is inferred from either that minor or this rank comparison.
