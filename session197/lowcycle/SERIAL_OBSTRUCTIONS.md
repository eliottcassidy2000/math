# Infinitely many minimal obstructions to one-or-two-cycle components

Session 197. **All-order explicit graph family and an infinite forbidden-induced-subgraph obstruction.** This packet does not classify all one- or two-cycle tournaments and does not settle the Hamiltonian-path spectrum.

## 1. The class and the conclusion

Let C(T) count directed Hamiltonian cycles modulo rotation. Let L_2 be the class of tournaments in which every nontrivial strongly connected component has C<=2. This includes all weak ordered joins of strong one/two-cycle tournaments and singleton components.

**Theorem.** For every k>=0 there is an explicit tournament F_k of order k+8 such that

- F_k is strong and C(F_k)=4;
- every proper induced subtournament of F_k belongs to L_2;
- its Hamiltonian-path counts begin

      105,187,349,641,1177,2167,...,

  and obey an exact tribonacci recurrence;
- its feedback vertex number is at least `2+floor(k/3)`;
- its minimum backward span is exactly two.

Consequently L_2 has infinitely many minimal forbidden induced tournaments, and cannot be characterized by any finite list of forbidden induced tournaments. The F_k form an induced-subtournament antichain. In particular, no fixed order of local forbidden tests can decide this low-cycle component class.

This obstruction concerns a finite forbidden-induced-subgraph strategy. It does not exclude a recursive gluing classification, a finite-state language for some separately specified ordering, or an arithmetic construction of Hamiltonian-path values.

## 2. Cycle monotonicity in a strong extension

We first justify that L_2 is hereditary.

**Lemma.** If U is an induced strong subtournament of a strong tournament T, then C(U)<=C(T).

**Proof.** Extend each directed Hamiltonian cycle on U while preserving the cyclic order of its old named vertices. For a current directed cycle D, a vertex v outside D that has both an in-neighbor and an out-neighbor on D can be inserted into some arc x->y with x->v->y.

If no outside vertex is mixed toward D, partition the outside vertices into A, which lose to all D vertices, and B, which beat all D vertices. Strongness of T forces both sets to be nonempty and some arc a->b with a in A and b in B. Otherwise one of D,A,B is separated by a directed cut. Replace any fixed arc x->y of D by x->a->b->y. This is again a directed cycle.

Choose all available vertices, arcs and ties deterministically, and repeat until every vertex is included. Deleting the inserted vertex names from the final cyclic order recovers the original cycle on U, so this extension is injective. QED.

An induced strong subtournament of any tournament lies in one of its strong components. Applying the lemma there shows that every induced strong component of a member of L_2 still has at most two cycles. Hence L_2 is hereditary.

## 3. The serial-seam construction

Take the rooted four-tournament

    A=[2,8,3,5],    root a=2,

and its converse

    B=[12,5,8,2],   root b=2.

They are strong, each has C=1 and H=5, and

    s_a(A)=2,       e_b(B)=2.                       (1)

Every proper induced subtournament containing the designated root has at most three vertices, and therefore has at most one Hamiltonian path starting, respectively ending, at that root. This last fact is immediate for transitive tournaments and directed triangles, the only possibilities at orders at most three up to ordered components.

For k>=0 place the blocks in the order

    A, z_1, z_2, ..., z_k, B.

Initially orient every edge between blocks forward. Reverse only the adjacent spine edges

    b -> z_k -> ... -> z_1 -> a,                    (2)

with b->a when k=0. Preserve all arcs inside A and B. This defines an actual tournament F_k; the singleton blocks are named vertices, not formal parameters.

It is also the ordinary one-seam operation from `ROOTED.md`: start at (A,a), append k singleton gadgets while moving the root each time, then append (B,b). Thus the gluing mechanism is inherited; the induced-minimality analysis below is the new-in-session use of it.

## 4. Every spanning cycle is forced across the spine

Across each cut between consecutive displayed blocks, exactly one arc points backward: the corresponding edge of (2). A directed cycle that has vertices on both sides of that cut must use that edge. In a Hamiltonian cycle of F_k, all spine arcs are therefore forced and occur consecutively as the directed path in (2).

The rest of the cycle is exactly an A-path starting at a, followed by a B-path ending at b. The forward edge joining these paths is legal; both endpoint blocks have four vertices, so it cannot coincide with the reversed seam. Hence

    C(F_k)=s_a(A)e_b(B)=2*2=4.                     (3)

This constructs Hamiltonian cycles and proves strongness at the same time.

## 5. Why every proper induced graph is low-cycle

Let U be a vertex set whose induced tournament is strong. If it meets more than one displayed block, it must meet a consecutive interval of the blocks. Moreover it must contain every spine vertex needed to cross each cut in that interval. A missing intermediate singleton or a missing relevant endpoint root leaves no backward arc across that cut, contradicting strongness.

There are only four cases for such a strong induced graph.

1. **Inside A or B.** It has at most one Hamiltonian cycle.
2. **Inside the singleton spine.** If its order is at least three, it has all adjacent arcs backward and all other arcs forward. Every adjacent cut forces its backward arc into a spanning cycle, leaving exactly one way to close from the first vertex to the last. Thus C=1. Smaller strong cases are singletons and have C=0.
3. **One endpoint block and some spine vertices, but not the other endpoint block.** The forced-cut argument gives C equal to the appropriate starting or ending path response of the selected endpoint subtournament, which is at most two by (1) and its proper-subset observation. The two-singleton boundary case cannot be strong and causes no exception.
4. **Both endpoint blocks.** All k spine vertices and both roots are required. The same forced-cut argument gives

       C(F_k[U]) = s_a(A[U intersect A])
                    * e_b(B[U intersect B]),       (4)

   except for a possible graph of order two when k=0, which is not strong. If U is proper, at least one selected endpoint block is a proper subtournament. Its response is at most one; the other response is at most two. Thus (4) is at most two.

These cases cover every strong induced subtournament. Therefore every proper induced subtournament, even if weak, has only strong components with at most two Hamiltonian cycles. It belongs to L_2.

If L_2 had a finite forbidden induced list, each F_k would have to contain a member of that list. It cannot contain one properly, because all its proper induced subtournaments are in L_2. Thus every F_k itself would have to appear in the list. Their distinct orders give a contradiction.

## 6. Exact path growth and unbounded feedback vertex number

For a rooted input write the inherited state x=(h,e,R). The singleton seam has matrix

    S=[[1,-1,1],[1,-1,0],[1,0,1]].

For A at its specified root, x_A=(5,1,9). For B at its specified root, the scalar output row is (5,-1,9). Therefore

    H(F_k)=(5,-1,9) S^k (5,1,9)^T.                (5)

The exact identity `S^3=S^2+S+I` proves

    H(F_(k+3))=H(F_(k+2))+H(F_(k+1))+H(F_k),

with initial values 105,187,349. The formula is not fitted from those initial values; it follows from the returning physical transfer and a finite exact matrix identity.

A contains the triangle on vertices 0,1,3, as does B after taking the converse. The k singleton spine vertices contain floor(k/3) additional vertex-disjoint triangles, obtained from consecutive disjoint triples of spine positions. All these triangles are mutually disjoint. Every feedback vertex set must meet each one, giving the claimed lower bound `2+floor(k/3)`. In particular this family leaves every fixed-apex class while maintaining exactly four Hamiltonian cycles.

### The antichain already lies at exact backward span two

Order A by its local vertices `[3,0,1,2]`, then list the spine in order, then order B by its local vertices `[2,1,0,3]`. Both four-vertex blocks have maximum backward-edge distance two in these orders. Their roots are respectively last and first, so every reversed seam in (2) is adjacent. Every other cross-block arc is forward. Hence every F_k has a span-two presentation, uniformly in k, including k=0.

No F_k has a span-one presentation. In any span-one tournament, a forward adjacent pair separates an ordered join, because every other edge across that cut is forced forward. Strongness therefore requires all adjacent edges to be backward. Those are the only backward arcs across their respective cuts, so every Hamiltonian cycle must use all of them; the remaining forward edge from the first vertex to the last closes exactly one cycle. Since C(F_k)=4, span one is impossible. Thus the minimum span is exactly two.

This gives an infinite forbidden-induced antichain inside a fixed narrow span class, even though the graphs have unbounded feedback vertex number. Finite frontier transfer and a finite forbidden-induced-subgraph basis are different kinds of finiteness.

## 7. Verification and source boundary

The graph construction and path transfer are based on the explicitly cited S195 one-seam law in `ROOTED.md`. The cycle-cut proof, minimality argument and infinite-antichain consequence are proved directly above. No unverified global low-cycle gluing assertion, asymptotic estimate or external order bound is used.

Run

    python3 verify_serial_obstructions.py > serial.normal.stdout
    python3 -O verify_serial_obstructions.py > serial.optimized.stdout
    cmp serial.normal.stdout serial.optimized.stdout

The standard-library checker constructs F_0 through F_5, of orders8 through13, and counts Hamiltonian cycles on **every vertex subset** using independently anchored subset dynamic programs. It checks that each full graph has four cycles and every proper induced graph has at most two. It separately checks full counts by direct cycle DFS, the actual disjoint triangle packing, the explicit span-two ordering, the endpoint profiles, the Hamiltonian-path counts and the exact matrix recurrence. The finite checks are implementation controls; the cut argument proves all k and all proper induced subgraphs.
