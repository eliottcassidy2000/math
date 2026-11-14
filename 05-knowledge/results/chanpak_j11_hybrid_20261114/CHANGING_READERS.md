# Exact changing-reader windows and a port-state return

These are physical hybrid families. They are not claims that the bare Wave-20 eleven-matrix orbit is cofinite.

## Reflective carry supplies the source pair

For an apex with outgoing binary word W over a transitive spine, the rooted state is
`(H,e,o)=(g(W),ceil(W/2),1+floor(W/2))`.
Here e and o sum Hamiltonian counts of the complement over even-/odd-order paths starting at the marked root, and
`g(0)=1, g(2u+1)=g(u), g(2u)=g(u)+u+(u mod 2)`.
The e/o formula follows by partitioning extensions of a path from the apex among subsets of the later spine vertices. Except for the choice where the last spine vertex is the first step, even and odd continuations pair by toggling that last vertex. The recurrence follows by appending a terminal contact bit.

Fix r≥2, an odd m=2z+1, and put A=m·2^(r−1).
Consider W0=m·2^r−2 and W1=m·2^r+1. Exact substitution gives
`H0=A+g(m−1)`, and
`H1=A+g(m)−m+1`.
Indeed W0=2(A−1) and g(A−1)=g(m−1) by appending ones; W1=2A+1 and g(A)=g(m)+(2^(r−1)−1)m+1 by appending zeros. Therefore the source difference is
`(H1−H0,e1−e0,o1−o0)=(-3z−p(z),2,1)`.
Geometrically the two low-bit limits are ...1110 and ...0001, which both reflect to the real binary number 0.1000...=0.0111.... The equal leading slope is an exact carry identity, not a numerical coincidence.

## A changing aperture cancels the jump

A port-chain seam appends a reader B, with all old-to-new arcs pointing forward except for the arrow from the reader's entry to the old root, and then selects a port in B. In the port coordinates above, the first transfer row is `(eta,2q_o,2q_e)`, where eta=H(B) and q_e,q_o are the analogous path sums ending at the entry. This is THM-5070's physical transfer. Its count difference here is `-eta(3z+p(z))+4q_o+2q_e`.

1. For any w≥1, choose the transitive reader of order w+1 entering at its sink. It has eta=1, q_e=q_o=2^(w−1). Choose z=2^w−1. The reader changes the phase by +3·2^w, canceling the source difference −3·2^w+2 down to exactly +2. The two final absolute counts are
   `(4^(w+1)−1)·2^(r−1)+{1,3}`.
   Given a target, try offsets 1 and 3, remove all powers of two from the positive difference, then require that odd part plus one to be 4^(w+1), with r−1≥1 and w≥1. This decodes this whole unbounded two-parameter family without enumeration.

2. For k≡3 mod 18, k≥3, let y=27·2^k−1 and choose as reader the converse of the one-apex tournament of address y, entering at its apex. Then eta=g(26)=19, since y is formed by appending ones to 26, and q_e=q_o=(y+1)/2. Since 2^18≡1 (mod 19) and 27·2^3≡7 (mod 19), z=(y−6)/19 is an odd integer. With m=2z+1, the jump cancels to +2 as above. The exact pair is
   `m(19m+14)·2^(r−1)+19g(z)+12+{0,2}`.
   The first pair at k=3,r=2 is 20815,20817. To decode a target try the allowed k up to log₂ N, then test the exact remaining power-of-two quotient. Both the core and the marked response can change with k.

These families are sparse; in the first family only O(log² X) distinct totals are possible below X. The result is exact physical selection, not universal placement.

## A two-phase matrix pump

The original Wave-20 state (H,E,L,h) changes to the HYP-9640 port coordinates by `e=(L−E)/2, o=(L+E)/2`. Its J becomes the singleton seam transfer
`J'(H,e,o)=(H+2e,o,H+e)` exactly.

At w=1 the source pair uses m=3. Let the transitive two-vertex reader enter at vertex 0 (the sink of the arc 1→0), and move the port to vertex 1. It has matrix C below. Write matrices with their three rows:
- C=(1,2,2;1,1,0;1,2,1).
- S=(5,10,8;4,5,4;5,8,5), a four-vertex reader with row-major adjacency bits 0001100011000110, entry 1, port 2.
- R=(1,4,4;2,3,2;2,4,3), a transitive three-vertex reader with row-major bits 000100110, entry 0, port 2.

The physical windows after C have states `(5A+1,2A+2,4A+1)` and `(5A+3,2A,4A+2)`, where A=3·2^(r−1). Their difference is v=(2,−2,1). Direct multiplication gives S v=−v; J'v=−w for w=(2,−1,0); and Rw=−w. Consequently for **every** r≥2 and a,b≥0, the two baselines after C, then S^a, J', and R^b differ in count by exactly ±2. During this process the order changes, the marked port changes, and the baseline grows without losing the odd phase choice. S and R are physical extra letters from HYP-9640's alphabet; they are not asserted to be Wave-20's eleven.

Only O(log³ X) such parameter triples give counts below X: C starts at order-of-magnitude 2^r, S multiplies H by at least 5 at each step, and every two repeats of R multiply H by at least 9. This is an exact multi-stage return, not all-target coverage.
