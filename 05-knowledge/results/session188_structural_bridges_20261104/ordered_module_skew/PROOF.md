# Positive ordered-cover skew inside an arbitrary exterior module

**PROVED:** if a tournament R has a module M=A->B, then the aggregate ordered-cover skew from B back to A is nonnegative, and is strictly positive when |M|>=3. This strengthens the marked-pair quotient criterion beyond exteriors with a global ordered cut. Consequently every five-vertex quotient inflated by nontrivial children, at least one of which is a pair, is separated from its quotient converse with all child interiors unchanged. In particular all strong prime five-vertex quotients are covered under this child-size hypothesis.

## 1. Aggregate cover notation

Let A,B be nonempty disjoint sets whose union M is a module of a tournament R, and suppose every edge between them goes A->B. Their internal tournaments are arbitrary. Let C_AB count ordered spanning covers (P,Q) of R by two nonempty directed paths with P ending in A and Q starting in B. Define C_BA analogously. Endpoints stay physically named until the indicated sums are taken. Put K=R/M, the tournament contracting M to one vertex.

Then

    C_BA-C_AB
      >= H(K)[q_2(A)H(B)+H(A)q_2(B)+q_2(A)q_2(B)].       (1)

In particular C_BA>=C_AB, with strict inequality whenever |M|>=3. The statement holds in arbitrary strong or non-strong R, and M may equal R.

## 2. A support-preserving injection

For a cover counted by C_AB, let P=P_out P_A, where P_A is the maximal terminal M-run. Because P ends in A and no edge goes from B to A, this run lies wholly in A. Likewise write Q=Q_B Q_out, where its maximal initial M-run lies wholly in B.

Map this cover to

    (P_out Q_B, P_A Q_out).                              (2)

The new cover ends its first path in B and starts its second in A. Both exchanged runs are nonempty. If P_out is nonempty its last vertex lies outside M, and if Q_out is nonempty its first vertex lies outside M. Module homogeneity validates the new joins. All support is preserved, no name is repeated, and the paths remain disjoint.

The new terminal M-run of the first path is exactly Q_B and its new initial M-run of the second is exactly P_A. Therefore the operation is invertible onto the subset of C_BA covers whose terminal first M-run is purely in B and whose initial second M-run is purely in A. Other visits to M earlier in the first path or later in the second path cause no difficulty. This gives C_BA>=C_AB without discarding any covers.

## 3. Quantified covers outside the image

Choose a Hamilton path of K. Split it at the contracted vertex into an exterior prefix and suffix. Make the first cover path from the prefix followed by a path inside M ending in B; make the second from a path inside M starting in A followed by the suffix.

There are three disjoint choices outside the image of (2):

- Split A into two ordered nonempty paths A1,A2 and choose a Hamilton path of B. Use A1 followed by B in the first slot and A2 in the second. The first terminal M-run is mixed. This gives q_2(A)H(B) choices.
- Choose a Hamilton path of A and split B into two ordered nonempty paths B1,B2. Use B1 in the first slot and A followed by B2 in the second. The second initial M-run is mixed. This gives H(A)q_2(B) choices.
- Split both sides into two ordered nonempty paths, and use A1 followed by B1 in the first slot and A2 followed by B2 in the second. Both marked M-runs are mixed. This gives q_2(A)q_2(B) choices.

All joins between A and B are directed, as are the exterior joins by modularity. Empty exterior prefix or suffix is allowed; both cover paths remain nonempty. The three classes are disjoint by their pure/mixed marked-run types. Within each class the resulting cover determines its two M-paths and the quotient Hamilton path, so the counting map is injective. Multiplying by H(K) proves (1).

Every nonempty tournament has positive H, and q_2(D)>0 when |D|>=2. Thus (1) is strictly positive if at least one of A,B contains two names. For two singleton sides, the lower bound is zero and the swapping map is a bijection; strictness is not claimed. When M=R, the formula reduces to the exact ordered-cut identity in the adjacent marked-pair packet.

## 4. Transfer to quotient converse using a pair child

Let T=Q[A_i] and U=Q^op[A_i], with unchanged named child interiors. Suppose A_z={p,q}, p->q, and put R=T-A_z. The adjacent marked-pair theorem proves that for disjoint unions A,B of remaining full child blocks,

    [a_q a_A c_p c_B](F_T-F_U)=C_BA(R)-C_AB(R).            (3)

If their union M is an ordered module A->B in R with |M|>=3, (1) makes this coefficient strictly positive. Hence one of its sixteen actual finite-difference queries separates T and U. The grouping is uniform on each other child, but is deliberately nonuniform on the pair child. Exact quartic data, not a cubic truncation, drives this separation.

This criterion also applies when the deleted quotient Q-z remains strongly connected, unlike the earlier global-cut corollary.

## 5. All five-vertex quotients under a pair-child hypothesis

**Theorem.** Let Q have five vertices. Inflate every quotient vertex by a nontrivial tournament, with at least one child exactly a pair. Then Q[A_i] and Q^op[A_i] have different full named two-port responses.

Choose a pair child at z. The quotient Q-z has four vertices.

If Q-z is not strong, split its linearly ordered SCCs into a nonempty prefix and suffix. Their inflated unions A,B form an ordered module M=R. They contain at least four physical names in total, so (3) is nonzero.

If Q-z is strong, it contains a directed triangle a->b->c->a. Let d be the fourth vertex. It cannot beat the whole triangle or lose to the whole triangle, since that would make d a source or sink SCC. If d beats only a, then {c,d} is a homogeneous pair: both beat a and both lose to b. If d beats a,b and loses to c, then {a,d} is a homogeneous pair: both beat b and both lose to c. Cyclically renaming the triangle covers every mixed case. Thus every strong four-vertex tournament has a homogeneous pair i->j.

Inflating that pair produces a module M=A_i->A_j in R. Both children are nontrivial, so |M|>=4. Apply (3) and (1). This proves the theorem for every five-vertex quotient, in particular every strong prime one, without assuming its isomorphism type.

The hypothesis that all children are nontrivial ensures the required ordered module has at least three physical names in the strong-deletion case. The theorem does not claim that singleton-child cases are collisions; many are already separated by exposed-anchor results. Likewise the special pair child is essential to this proof, not asserted necessary for separation.

## 6. An all-order cyclic strong-prime family

Let r>=2 and let Q_r be the tournament on residues modulo 2r+1 with

    u->v iff v-u modulo (2r+1) lies in {1,...,r}.

These quotients are strong and prime. Strong connectivity follows from the directed cycle of successive residues. For primeness, suppose a module contains distinct u,v, oriented u->v. Write v=u+d with 1<=d<=r and put w=u+r+1. Then w->u while v->w, so modularity forces w into the module. The pair u,w has opposite incidences to every other vertex: the outgoing sets on the remaining names are disjoint and exhaust those names. Thus every other vertex must also belong to the module.

For any z, deleting z leaves the homogeneous pair

    i=z+r,  j=z+r+1,  with i->j.

After rotating z to zero, the names 1,...,r-1 beat both i,j, and both i,j beat r+2,...,2r. These are all remaining names. This proves the pair assertion without a finite classification assumption.

**Corollary.** Inflate Q_r by arbitrary nonempty child tournaments. If the child at some z is a pair and |A_(z+r)|+|A_(z+r+1)|>=3, quotient converse with all child interiors unchanged has a strictly positive grouped quartic of the form (3). In particular, any nontrivial-child inflation of any Q_r with at least one pair child is separated.

The exterior module is exactly A_(z+r)->A_(z+r+1); all other child interiors and sizes are unrestricted. The conclusion holds at unbounded quotient order and unbounded sizes of the other children. It does not claim that every strong prime quotient belongs to this cyclic family.

## 7. Scope and verification

This closes a concrete part of the strong-prime quotient frontier. It does not settle arbitrary child tournaments all of order at least three, higher prime quotients outside the cyclic family and lacking the indicated ordered module after deletion, arbitrary unequal C3 profiles, or arbitrary strong outer ambient lifting.

The checker verifies the literal injection on all labelled tournaments through order four, every eligible ordered module and both marked cover slots, plus a larger reentry control. It tests the four-vertex classification and the resulting block selection for all five-vertex quotients and all choices of deleted root. It retains sixteen-query quartic controls for representatives of every strong prime five-vertex quotient and all five root positions, using nontrivial children. It also checks the cyclic quotient structure through order fifteen and actual larger-quotient query batteries. Exact finite classifications here support, but are not substituted for, the elementary all-order injection and four-vertex argument above.
