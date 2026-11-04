# Fanout surplus is consumed at constants and transferred at aliases

**Status: PROVED for the explicit ordinary DAG operations below.** This is
a sequential refinement of `PROOF.md`, not a global circuit lower bound or
a certification of an affine elimination algorithm.

## 1. The conserved resource and permitted operations

For a finite acyclic gate DAG D, let d_D(g) count distinct gate children,
not pins or output references. Define the gate-only potential

    Phi(D) = sum_g max(d_D(g)-1,0).

Primary inputs, external boundary values, and constants are not gates in this
sum. Allow only the following operations, with no new gates:

- replace a certified constant gate by a constant;
- replace a gate that has become a literal of a predecessor value by that
  literal, redirecting every use and absorbing signs into user truth tables;
- delete dead gates;
- substitute constants or input literals for primary inputs, followed by
  the same simplifications.

Repeated pins are allowed and are counted only once in distinct child sets.
This statement concerns the DAG, so it does not independently authorize a
restriction, certify semantic constancy, or make an affine coordinate change
free. Any new-gate resynthesis lies outside this potential claim unless its
potential increase is separately paid.

## 2. Literal rerouting transfers rather than creates fanout

**Lemma.** Each permitted operation weakly decreases Phi.

Constant and dead-gate deletion only remove a gate and incoming/outgoing
edges, so cannot increase any other gate's child degree.

For a literal reroute g -> p, first suppose p is a gate and g has at least
one gate child. Write D=d(p)>=1, d=d(g)>=1, and

    t = |(children(p) minus {g}) intersect children(g)|.

After rerouting, p's degree is exactly `D-1+d-t`. The removed contribution
of g is d-1, while p's contribution increases by exactly d-1-t. Thus their
combined potential drops by t. Any other gate parent q of g loses the child
g, which can only decrease its own potential. In a two-distinct-parent gate,
that additional decrease is exactly one if q is a gate with degree at least
two, and zero otherwise. No other child degree changes.

If p is a primary input or boundary value, its increased fanout is not in Phi;
g's whole d-1 contribution disappears. If g has no gate child (in particular
if it is the output of a live DAG), d=0, and p loses
its edge to g while gaining no gate child, which again cannot increase Phi.
Repeated-parent cases are handled by the same distinct-child counting.
Dead cleanup after rerouting can only decrease the potential further. QED.

In particular, if p was g's sole gate parent with no shared users and the
other pin was a constant, the full d-1 branch contribution moves from g
to p without being consumed. This is the exact boundary between an alias
transfer and a genuine constant deletion.

## 3. A no-reuse theorem for constant-region surplus

Consider a sequence of restrictions and permitted simplifications. At step j,
take a current live DAG and a certified constant region C_j as in `PROOF.md`,
excluding the output, and put

    E_j = sum_(g in C_j) (d_current(g)-1).

Assume the step actually removes all gates of C_j and its constant/unary
boundary and introduces no new gate. Then

    Phi(D_j) - Phi(D_(j+1)) >= E_j,
    sum_j E_j <= Phi(D_0) - Phi(D_final) <= Phi(D_0).  (1)

**Proof.** Delete the C_j gates in topological order before performing any
boundary reroutes or dead cleanup. When each C_j gate is removed, all its
children in C_j occur later and its children outside C_j still exist. Its
child degree is therefore the degree measured before the step. Deleting it
consumes its entire d-1 contribution. Changes at other parents only decrease
Phi. The subsequent reroutes and cleanup do not increase Phi by the lemma.
Sum over the steps. QED.

Thus the same gate-fanout resource cannot pay E twice in consecutive
constant-region certificates. It can move through arbitrarily many unary
boundaries before the step that finally consumes it. Counts must be taken
in the current DAG, not from historical degrees.

For a live one-output binary circuit with s>=1 gates and n_used primary inputs,
all of which have positive degree, `Phi = e_gate_gate-(s-1)`. Since the
total number of distinct gate input edges is at most 2s and there are at
least n_used primary-input edges,

    Phi <= s-n_used+1.

This optional bound is in ordinary input coordinates; it is not a replacement
for the source algorithm's separately maintained parameter dimension.

## 4. Sharp hostile: arbitrarily many boundary charges reuse one branch point

Here is a precise counterexample to the tempting stronger assertion that the
sum of fanout excesses of *all deleted unary boundary gates* is also bounded
by initial Phi, or provides additional distinct-gate surplus at every step.

For k>=1, use independent bits a,b,y,z and a_i,b_i for 1<=i<=k. Build

    G_0 = a AND b,
    A_i = a_i AND b_i,
    G_i = G_(i-1) XOR A_i,            i=1,...,k,
    L = G_k XOR y,
    R = G_k XOR z,
    O = L AND R.

There are 2k+4 gates and 2k+4 used inputs. Every gate is live. G_k is the
only gate with more than one child, and it has exactly two, so Phi=1.
The output remains nonconstant in every restriction below.

Fix a_k=0, then a_(k-1)=0, down to a_1=0. At step i:

- the constant region is C={A_i}, with r=1, f=E=0;
- its boundary is B={G_i};
- A_i and G_i are exactly the two deleted gates;
- G_i has two children at this moment, but its branch contribution transfers
  to G_(i-1), which becomes their shared parent;
- Phi remains exactly one.

The sum of the contemporaneous boundary excesses `d(G_i)-1` is k, although
initial Phi is one and the actual deletion count is exactly 2k. Therefore
adding one presumed extra deletion for each such boundary would claim 3k
deletions where only 2k occur. The first invalid implication is treating a
boundary gate's fanout as consumed when literal rerouting transfers it.

After these k steps, fix a=0. Now C={G_0}, its two boundary gates are L,R,
and the surviving circuit is the single gate y AND z. This final step has
r=1,f=0,E=1 and exactly three deletions. Phi drops from one to zero.

Across the complete trace, the original surplus identity is exact:

    total deletions = 2(k+1)+1 = 2k+3,
    sum E_j = 1 = Phi(initial)-Phi(final).

The family is a literal DAG counterexample at arbitrarily large k. It makes
no minimum-circuit or small-affine-flat claim. Its role is to separate a
valid constant-region charge from an invalid repeated boundary charge.

## 5. Reproduction and consequence

`verify_sequential.py` exhausts literal reroutes in every live distinct-pin
binary DAG through five gates on two external inputs, checking nonincrease
and the exact transfer law before dead cleanup. It then constructs the sharp
ladder for k=1,...,32, verifies each contemporaneous C/B set and potential,
and checks the Boolean functions exhaustively at k=1,2,3.

Run normally and with Python optimization; the two receipts must agree.
These are finite controls for the displayed all-k proof. The original
global obligation remains: one must force actual consumed surplus in the
relevant circuit class. A large historical count of unary branch points is
not such a proof.
