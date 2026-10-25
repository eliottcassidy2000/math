# Independent review of the protected-port boundary theorem

2026-10-25. Reviewer: the counting/recurrence agent, independent of the
port theorem's author. **PASS:** the stronger theorem in
[NEXT_BOUNDARY.md](NEXT_BOUNDARY.md) follows for its full stated range,
using the literal exceptional repair in
[GENERAL_TAIL.md, sections 2--3](../../trijection_hourglass_ports_20261025/endpoint/GENERAL_TAIL.md).
No repair to the statement or proof was required.

The proved scope is: a directed, closed, blue marked Hamiltonian word
P=(0,...,n-1), n>=9, incoming interior port set I with min I>=5, and
{n-4,n-2} contained in I. Every such source has three directed black K
moves to the exact labelled target K_P(T). Other ports and all remaining
arcs remain arbitrary. This is not an all-port theorem or a proof of
unrestricted black-core connectivity.

## Proof audit

**The 4->2 branch.** In the first word L1=(0,2,3,...,z,1), vertex2 is
at position1 and vertex4 at position3, so its actual reversal sends them
to z and n-3. Original blue reflection sends the arc a->2 to (n-3)->z.
The guard 4->2 therefore gives the direct contradiction
T(2,4)=0 versus T(n-3,z)=1. It uses no regularity bound, no cyclic
classification, and no restriction on the other endpoint ports.

After that move, a is the first vertex of L1, whose only retained
incident path arc is a->2. Its outdegree is k+1. At z, one incoming and
one outgoing path arc are retained, so its new indegree is k. The second
word's reversal exchanges a and z and cannot be anti. The final reversal
product with r_P has two disjoint four-cycles; n>=9 ensures their
supports are disjoint. Such a permutation cannot be a tournament
automorphism: an even-order automorphism has a nontrivial involution
power, and any transposed pair would reverse its own arc. The inherited
legal joins and path-edge XOR therefore complete this entire branch.

**The 2->4 branch.** The Q words remove vertex3, insert a immediately
before the first port h, then restore3. The new joins use precisely
2->4, z->3, a->(h-1), and h->a. These follow from the chord guard,
the reflected protected port n-4, and h=min I. Since h>=5, h-1 is not
the removed vertex3. The protected port c=n-2 gives the first witness
T(c,a)=1 versus T(3,2)=0. After Q1, both relevant nonpath pairs have
been complemented, giving the second witness S(z,b)=0 versus
S(3,2)=1. Neither argument imposes a middle-orientation condition.

The third word is the insertion word Q=(1,...,h-1,0,h,...,z). The
rotation lemma in [README.md, section 2](README.md) computes the actual
product r_Q r_P, including overlap. In the disjoint case its a-orbit
is {0,...,h-1}. In the even-overlap case its a-orbit also lies below h
and has n/2 vertices. In either orbit, a has exactly one outneighbour
in B, since N^+_B(a)={b,z} union I. Regularity of an automorphism orbit
would respectively force h=3 or n=6, both excluded. Odd overlap gives
a transitive automorphism and forces n=2k+3 from deg^+_B(a)=k+1.

**The exceptional size is fully divided.** At n=2k+3, the k-1 interior
ports all lie in [h,n-2], so h<=k+3. The values h<=k+1 are already
in the disjoint case. For h=k+2, one has n=2h-1; a is the central
position of Q and its reversal fixes a and exchanges b,z. It would
equate B(a,b)=1, a retained P edge, with B(z,a)=0, the complemented
closing edge. This is a direct anti-reversal contradiction, independent
of the other ports.

For h=k+3 the containing interval [h,n-2] has exactly k-1 vertices.
It is therefore filled: I={k+3,...,2k+1}. This is exactly the inherited
consecutive tail, not just a port set of the same size. Inclusion of
n-4=2k-1 forces k>=4. If Q is blue, transitivity and the complete
a-row determine every labelled arc of B, giving precisely the cyclic
source stated in GENERAL_TAIL.md. No family of unexamined sources is
hidden at this boundary.

**The inherited repair is applicable as written.** For the exceptional
source, the three alternative words (A,w,u,v,s,Z), (A,w,u,Z,v,s),
(A,u,Z,v,w,s) start at the original T; they are not appended to a
partially executed Q sequence. Their nonpath joins have cyclic
differences2,k, and k-2 outside the defining set D and hence the
claimed source orientations after complementing. Their anti-reversal
witnesses are respectively a->v versus (k-1)->z, a->w versus w->s
after the first move, and a->3 versus z->s at B. The last difference
is k+5 in D exactly throughout k>=4, including the equality endpoint
k+5=2k+1 at k=4. The nonempty Z block and all word positions are valid
throughout that range. Thus the inherited proof repairs the exact
remaining source at every possible k.

Blackness may be tested at B for the last word because its unordered
path-edge set is invariant under its positional reversal; complementing
the other edges preserves the truth or failure of the anti-reversal
condition. The literal path-edge XOR identifies B separately from this
blackness test. No isomorphism quotient or clone contraction occurs.

## Boundary and finite verification

The one-hole corollary follows without an additional estimate: for k>=5,
its interval contains n-4 and n-2, its minimum is n-k-1, and n>=k+6
is precisely h>=5. At n=2k+3 its minimum is k+2, so the new fixed-a
witness applies and the consecutive-tail exception cannot occur. The
k=4 one-hole family lacks n-4 and is correctly excluded from this theorem.

At h=4, after the first two Q words the pair {0,3} has been complemented
twice. Its original direction is0->3, since3 is not a port. The proposed
last word needs3->0 and is illegal. This is a real failure of these words,
not a disconnection conclusion. The retained n=15,I={4,11,13} cell
confirms that this boundary is feasible.

Both maintained reproduction commands passed during this review:

    python boundary_check.py
    python -O boundary_check.py

They certified all2,044 protected-port affine guard cells for n=9..17,
all660 declared one-hole cells for k=5..24 and n=k+6..2k+7, and the
conditional consecutive-tail exception for k=4..40. The protected cells
contain four feasible conditional exceptions. Each certified cell also
has two literal model replays. The parity-equation source system is
satisfiable before the anti-reversal test; contradiction after that test
certifies every source completion in that cell. The checker does not
discard a source because its proposed repair fails.

An additional independent ephemeral Python check used no affine solver
or imported port implementation. It directly constructed the cyclic B
from its difference set for all77 values k=4..80, formed T by the literal
K_P complement, verified the complete tail row and blue P, and replayed
all three inherited repair words using arc sets. Every join, actual
anti-reversal failure, and exact final labelled target passed. Separate
literal position checks covered the new L1 witness for n=9..200 and
the fixed-a endpoint witness for k=3..80. These controls also passed
under optimized Python with explicit exceptions.

The finite ranges corroborate the formulas. The unbounded conclusion
comes from the case split and the all-k exceptional repair audited above.
This strengthens sufficient scope; it does not retract a previously
false theorem, eliminate the h=4 word obstruction, or close the remaining
black-core connectivity problem.
