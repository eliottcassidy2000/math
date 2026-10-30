# Independent audit: the AP2 row deletes at least four gates, sharply

**PROVED local statement under the stated legal-state hypotheses.** This audits
an improvement to the AP2 row in the pinned provisional
[elimination candidate](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/blob/4a61f024a4993592d97b49ca7d00979029b7fece/submission/elimination.md).
It does not certify that source's complete circuit theorem, disperser construction,
or benchmark status. No literature-priority claim is made.

## Hypotheses retained

There is one output and an acyclic B2 gate part, with all dead gates removed.
No nonlinear gate feeds the affine system. The chosen gate G has parameter
parents p,w from distinct marked pairs (p,r),(w,h). Each of those four
parameters has exactly one child, a gate, and the partner children Hr,Hh are
different from G. The AP2 restriction fixes r,h safely and w to the killing
value for G. It is legal on the product domain and leaves a nonconstant output.
The source ensures the last fact by a'=a-1>d, not by ignoring terminal outputs.

Each old gate with a fixed pin is deletable after literal/constant propagation.
Counting distinct old gates avoids double-counting common descendants. Unary
rerouting may remove a gate even when its replacement remains nonconstant.

## The four-gate proof

### Distinct partner children

Suppose Hr != Hh. The three distinct gates G,Hr,Hh are deletable: G is
constant, and Hr,Hh each receive their fixed partner. Take a directed path
from G to the output. Such a path exists because dead gates were exhausted.
While that path stays in S={G,Hr,Hh}, every path gate becomes constant.
Indeed Hr or Hh has one pin equal to its fixed partner, and its other pin is
its preceding constant path gate. A path ending inside S would make the
output constant, contrary to the hypothesis. Therefore it first exits S at
a fourth gate K. K has a constant input from its predecessor and is deletable.
The path argument is about old graph incidence; early deletion/rerouting of
an object already contributes its deletion and cannot erase the count.

### Shared partner child

Suppose Hr=Hh=H. Its two pins are r,h, so H becomes constant; G also becomes
constant. Neither is the output. Neither is a child of the other, by their
specified parameter pins. Both have at least one child. If their child sets
have a union of size at least two, those two distinct children are deletable,
giving four gates together with G,H. Otherwise their common unique child is K.
Both its pins are G,H, so K is constant. K cannot be the output, and a child L
of K is a fourth distinct deletable gate. Acyclicity makes L distinct from the
three earlier gates. This covers full sharing, not just tree-like circuits.

Consequently AP2 has Delta s<=-4 with the same Delta I<=0, Delta q=-2,
Delta a=-1 as before. The proof uses no Boolean table assumption on descendants.

## General constant-antichain deletion lemma

Let r gates in a live single-output B2 DAG be pairwise ancestry-incomparable.
Suppose a legal restriction makes their outputs constant and the final output
remains nonconstant. Then constant/unary propagation deletes at least 2r distinct
old gates, including the r seeds. Equations may precede this DAG, but no gate may
feed an equation; gate-ancestor liveness and the nonconstant-output condition are
required.

Proof. Form C from the seed gates by closing under gates whose two input parents
already lie in C. Every gate in C is constant. Let B be all gate children of C
outside C. Every gate in B has a constant pin, hence is deletable. A seed receives
no edge from C: every C vertex descends from a seed, and such an edge would make
two seeds comparable (or create a cycle). Every other C vertex receives at most
two edges from C; every B vertex receives at most one distinct edge from C, or it
would belong to C. The output is outside C, so liveness supplies an outgoing edge
from every C vertex. Select one such edge per vertex. Counting selected incoming
edges gives

    |C| <= 2(|C|-r)+|B|, hence |C|+|B| >= 2r.

This counts gates, not paths or pins, so arbitrary shared descendants are allowed.
Further constant propagation from one-pin-forcing truth tables only adds possible
deletions and is not needed for the argument.

The factor two is sharp: take independent AND seeds A_i=a_i b_i, boundaries
B_i=A_i XOR x_i, and an AND tree on the B_i as output. Set every a_i to zero.
Exactly the r seeds and r boundaries disappear, and the surviving AND tree is
minimum for conjunction of the r independent x_i (for r=1 the output is x_1).
The ancestry condition is necessary. In the live chain G1=ab, G2=G1*c,
K=G2 XOR x, setting a=0 makes both comparable seeds G1,G2 constant but deletes
only three gates, with nonconstant output x. Thus counting arbitrary constant
seeds as though they were an antichain is false.

`check_antichain.py` enumerates every distinct-pin B2 DAG topology through five
gates over two external input vertices, with one output and no dead gates. For
all eligible antichains whose formal constant closure excludes the output it
checks the incidence inequality. It examines 343 live topologies and 1,544
eligible seed cases, and separately verifies the comparable-seed hostile. This
is a finite structural control, not a proof that arbitrary forced seed values
are realizable under a particular restricted domain. The uniform incidence
proof above supplies the all-size conclusion.

## Sharpness with a minimum core and a reachable loop-eligible state

For any M>=4 use independent input bits

    p,w,r,h,z,t,u, x_1,y_1,...,x_M,y_M.

Let IP_M=sum_i x_i y_i. Define the core circuit

    G=pw, Hr=rz, Hh=ht, K=G+u,
    J=Hr*Hh, O=J+K, Q=O+IP_M.

Here + is XOR and juxtaposition is AND. Compute IP_M with M AND gates and
M-1 XOR gates. The core has 2M+6 gates on 2M+7 input bits. Mark the pairs
(p,r),(w,h), each excluding the corner (1,0). Both safe partner values are one.
All paired parameters have one child, the partner children are distinct, and
there is no earlier deletion, affine move or preliminary paired row. G may
be chosen first in topological order.

The AP2 restriction r=h=1,w=0 deletes exactly G,Hr,Hh,K and leaves

    Q = zt + u + IP_M.

No further D1-D4/M simplification applies. The residual has 2M+2 gates and
2M+3 essential input coordinates; p is unused but retained in the parameter
map. The fan-in-two connected-ancestor bound gives 2M+2 gates as a lower bound,
so no extra semantic simplification of this restricted function can save a
fifth gate.

The pre-restriction core is minimum as well, even when agreement is required
only on the paired product domain. Every one of its 2M+7 coordinates has a
legal single-coordinate sensitivity witness. For p use r=w=h=1,z=t=0; for
w use h=p=r=1,z=t=0. For r use p=w=0,h=z=t=1, and for h use p=w=0,r=z=t=1.
For z or t set the other three factors r,h,t or r,h,z to one and p=w=0.
For u use zero controls; for x_i or y_i set its partner to one and all other
controls and IP coordinates to zero. Unspecified coordinates are zero. Thus
any agreeing circuit must depend on all 2M+7 coordinates and uses at least
2M+6 gates. Restricted input dependence is proved by legal same-fiber witnesses,
not inferred from unrestricted truth tables.

To make this state genuinely reachable from q=0, prefix two generators

    B=p(1+r), D=w(1+h), C=Q+B, E=C+D,

and use E as output. Place B,D before the core in a topological order. GP on
B deletes B,C and creates the first marked pair. GP on D then deletes D,E
and creates the second pair. The core Q remains, exactly as above. Each GP
has one child and exactly two deletions. No preliminary paired row intervenes.
The full function is Q+B+D. This is an explicit circuit construction, not an
abstract inventory walk.

Every constant affine flat of the full function has dimension at most M+7.
Indeed intersect its direction with the kernel of the seven control-coordinate
projections. This leaves a constant affine flat of IP_M of dimension at least
D-7. Its direction is isotropic for the nondegenerate IP polar form, so its
dimension is at most M. Therefore D<=M+7. We may take d=M+7.
The safe-core dimensions along the trace are

    2M+7 -> 2M+6 -> 2M+5 -> 2M+4.

For M>=4 the first three are at least d+2, and the final one exceeds d.
Thus both GP rows and the sharp AP2 row satisfy the actual loop guard, and the
constant-output exclusion is justified. This is not a small-d asymptotic
hardness example: d is large, and no new circuit lower-bound record follows.

## Independent reproducibility

`check_ap2.py` imports no producer code. It constructs the entire B2 DAG and
performs literal constant/literal propagation for all three steps. For M=4 it
checks all 32,768 assignments for the relevant domain restrictions, confirms
all paired-domain sensitivity witnesses, and confirms the exact gate counts
18 -> 16 -> 14 -> 10. It checks distinct child sets and loop parameters.
Run normal and optimized Python; both use explicit exceptions and must return
the same JSON. Finite checks corroborate the displayed general proof, not the
converse. They are not a Lean proof or external human review.
