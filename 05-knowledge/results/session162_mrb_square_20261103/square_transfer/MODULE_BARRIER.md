# A run-number barrier and a nonuniform escape

## Prime-modulus module barrier

For a nonempty tournament B, let C_k(B) count ordered covers of its vertices by k nonempty directed paths, and F_B(r) count ordered r-slot path covers allowing empty slots. Then

  F_B(r)=∑_{k=1}^r binom(r,k) C_k(B).

The symmetric group S_k acts freely on the k-tuples of paths, so k! divides C_k(B). Let p be prime. If F_B(r)≡0 modulo p for all 1≤r<p, integer binomial inversion gives C_k(B)≡0 for 1≤k<p; and p divides C_k(B) automatically for k≥p. Hence every C_k(B) is divisible by p.

Now let A be any tournament with a vertex v replaced by B, meaning each outside vertex sees every vertex of B in the same direction. Every Hamilton path of A[v←B] decomposes uniquely into its maximal B-runs, interspersed with single outside vertices. For each k, fix an order of the outside vertices and k insertion gaps. Whether their adjacencies are legal depends only on A, the order, and the selected gaps, not on the endpoints of the internal B paths. Thus there is a nonnegative integer a_k(A,v), independent of B, such that

  H(A[v←B])=∑_{k≥1} a_k(A,v) C_k(B).

It follows that p divides H(A[v←B]) for every outer tournament A and every substitution vertex v, regardless of the number of runs. This proves the no-go for all one-hole modular contexts, not merely contexts analyzed using one or two runs. The argument even holds for a composite modulus m by requiring the vanishing of F_B(r) for 1≤r<κ(m), where κ(m) is the least k with m dividing k!.

For p=3 it suffices that H(B)=F_B(1)≡0 and F_B(2)≡0, since C₂=F_B(2)−2H(B). The full ordinal-sum blocks B built from a square-sum Hamilton path on [n], n≥4, satisfy both conditions by PRIME_OBSTRUCTION.md. Therefore arbitrary outer substitutions, including those allowing three or more B-runs, cannot escape 3-divisibility. The elementary divisibility statement allows the small exceptions 3 as a prime count and 21 as a count 7p; to exclude 21 one may invoke the separate S157 gap theorem. It makes no claim about tournaments in which the outside vertices see different B vertices differently.

## A one-vertex nonuniform escape

The square-sum edge 9→7 has 7+9=16. Let B be the triangle 7→8→9→7. Add vertex 10 with 10→7, 10→8, and 9→10. Vertex 10 sees B non-uniformly, so B is no longer a module. The resulting four-vertex X has two cyclic triangles, H(X)=5, F_X(2)=32. It is isomorphic to S158's carrier X: map its labels 7,8,9,10 to S158 labels 2,3,0,1 respectively. Hence appending n transitive vertices and applying the strong wrapper G yields

  H(G(X⊕T_n))=5+2^(n+5).

For n=2 this gives an eight-vertex strong tournament with exactly 133=7·19 Hamilton paths; for n=8 it gives 8197=7·1171. These counts were already realized by S158. The contribution here is the explicit mechanism: one non-uniform attachment to a real square-sum edge defeats a barrier that every one-hole modular substitution obeys. Only the edge 9→7 is claimed to have square-sum endpoints; the other edges of X and the wrapper are not constrained that way.

## Controls

module_runs_check.py verifies C_k(B)=(3,6,6) for the triangle and compares the full run expansion to an independent endpoint-subset Hamilton-path count for every outer tournament of order at most 5. It verifies the nonuniform X has (H,F₂)=(5,32) and directly counts the eight-vertex 133 witness. Run under Python or Python -O. These finite computations are controls, not proofs of the all-order theorem.
