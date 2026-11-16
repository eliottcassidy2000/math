# The 3→5→7 divisor ladder as a refinement about a converse axis

The archived divisor balance counts nontrivial proper divisors F,
proper distinct prime divisors U, and nontrivial proper squarefree
divisors S. For distinct primes p,q,r and N=p²qr, the
balance is F=U+S=3+7=10. See
[the earlier divisor-balance theorem](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/arithmetic_braids_20260917_divisors.md).
Here the 3→5 is a **different count**, about observable states of these
ten divisors. The earlier
[arithmetic chart and metagraph correspondence](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/trijection_reset_20261020/arithmetic.md)
uses exponent triples (a,b,c), where a ranges from 0 through 2
and b,c from 0 through 1, excluding (0,0,0) and (2,1,1).
The two involutions are C(d)=N/d and Q, exchanging q,r;
Q is the chirality flip in the archived dictionary. Their product
CQ is independent of both on the ten divisors.

Regard x=b-c as a three-state observer. The three fibers, from
x=-1 to x=+1, have sizes 3,4,3. The central fiber is
{p,p²,qr,pqr} and is fixed pointwise by Q. As the user
suggested, each outer fiber can split transversely, leaving the central
state alone. More precisely, map the middle to (0,0) and the
outside to (x,y)=(b-c,2|a-1|-1).
The new states are the center and the four points (±1,±1).
The aligned outer pairs are {pq}↔{pr} and {q,p²q}↔{r,p²r}.

Both C and Q induce the same horizontal reflection, preserving
the transverse coordinate. Their product CQ is not visible.
The split is unique under the following constraints: preserve the
coarse chirality; collapse the center; have two cells per outer row;
and require both C,Q to descend. Indeed Q identifies the two
outer partitions as the same partition of a=0,1,2. The action of
C forces that partition to be invariant under a↦2-a.
The singleton must therefore be {1}. In particular there is no
extra vertical reflection lifting to the divisors: the two outer
fibers that vertical reflection would have to swap have sizes 1 and 2.

If the center stays one point, the **minimal faithful** next chart
has seven points. Retain y=a-1 on each outer side instead. Then
Q flips x, C flips x and y, and CQ flips y.
The exact information ladder is 7→5 by forgetting the transverse
sign, then 5→3 by forgetting the whole transverse coordinate.
Under Q its three states are one fixed point plus one free pair,
five states one fixed point plus two free pairs, and seven states
one fixed point plus three free pairs. This is different from
quotienting the original ten divisors by complement: that also
produces five objects, but Q then fixes three and exchanges just
one pair. This check is useful for comparing any proposed literal
metagraph model: an equal state count does not certify an equal
involution or an adjacency-preserving map.

Reproduce the small checks with python3 check_divisor.py.
