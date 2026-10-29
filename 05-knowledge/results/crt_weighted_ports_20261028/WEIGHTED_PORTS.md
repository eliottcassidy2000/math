# Weighted shared ports without uniformizing their sources

**PROVED scoped / FINITE-EXACT; same-agent audit.** This extends the inherited
[joint-port moments](../lucas_prime_reflection_towers_20261028/PORT_MOMENTS.md)
and Session93's [complete tuple invariant](../session93_reversal_lucas_clocks_20261029/ports/PORT_ORBITS_AND_LCA.md).
Incoming Session94 already gives tag-invariant weights, full cycle indices,
unordered subsets and positive depth templates. Here the additional target
is **arbitrary vertex weights and actual joint source mass**, which need not
be invariant under the group. No new unrestricted A000568 formula is claimed.

## 1. Arbitrary weighted fixed-set distributions

Use A_k,C_k,alpha from [LINEAR_REPAIR](LINEAR_REPAIR.md). For a uniformly
sampled group element g and arbitrary real vertex weights w, set
F_g(w)=sum_(v fixed by g) w(v). Let A_k(w;t) and C_k(w;t) be its formal
exponential moment series for the respective groups. For a singleton,
both are exp(tw). Split w into its three child restrictions w_0,w_1,w_2.
Writing beta for child digit negation gives the exact recurrences

    A_k(w;t)=2/3 + (1/3) product_(b=0)^2 A_(k-1)(w_b;t),       (W1)
    C_k(w;t)=C_(k-1)(w_0;t)
               A_(k-1)(w_1+w_2 composed with beta;t).        (W2)

Proof. A nonzero root rotation fixes no vertices. Conditional on zero
rotation, the three child automorphisms are independent uniform elements.
For C_k, the zero child and paired child factor are independent, but the
same b acts on the two other children, conjugated by beta on one side.
Its fixed vertex j contributes w_1(j)+w_2(beta(j)). This proves both formulas.
They also hold with formal multivariate weights or signed scalar weights.

Setting all weights to one recovers the critical population X_k and the
centralizer population Y_k, with E X_k=1, Var X_k=2k, E Y_k=2k+1 and
Var Y_k=4k(k-1). These are inherited specializations, not new observations.

**Weighted Burnside boundary.** For any weight mu on ordered tuples,

    (1/|G|) sum_g sum_(x fixed by g) mu(x)
       = sum_(orbits O) [sum_(x in O) mu(x)]/|O|.             (W3)

Each x is fixed by |G|/|O| group elements, which proves the formula.
Only for an orbit-invariant weight is this a sum of one common weight per
orbit. General weighted moments give orbit-*averaged* mass, not the original
source mass of each orbit. This distinction is essential for arithmetic weights.

## 2. An exact joint-orbit mass recurrence

For an ordered named-port tuple p=(p_1,...,p_r), repetitions allowed, and
separate source functions w_i on vertices, define

    M_G(p;w_1,...,w_r)
      = sum_(y in G.p, distinct tuples) product_i w_i(y_i).   (W4)

The empty tuple has mass1; at depth zero the mass is product_i w_i(0).
These bases also handle empty blocks correctly. Partition the port indices
into I_b according to their first digit b, and let p_i' be their suffix.

For A_k, sum over the root shift s in F3, and take a product of three child
masses. In child b, use ports (p_i')_(i in I_b) and weights
j↦w_i((b+s,j)). The three root-shift images are disjoint for a nonempty
tuple, and each child group acts independently. Therefore

    M_(A_k)=sum_s product_b
       M_(A_(k-1))((p_i')_(I_b); (w_i(b+s,·))_(I_b)).        (W5)

For C_k, multiply the zero-child C_(k-1) mass by **one** A_(k-1) mass
formed from *all* ports in blocks1 and2. For a side1 port retain p_i' and
j↦w_i(1,j). For a side2 port use beta(p_i') and j↦w_i(2,beta(j)). This is
one diagonal orbit, not the product of two orbit masses. The direct-product
factorization of C_k proves the rule. All-one weights recover orbit size.

This preserves source mass exactly: summing (W4) over the tuple orbits gives
product_i sum_v w_i(v). It does not produce a full transition operator under
arbitrary later native moves. A correlated source law mu(y_1,...,y_r) needs
that joint tensor, or a valid sum of product tensors, rather than just its
marginals. Moving labels by g requires pushing the weights by g as well.

**Smallest hostile.** In paired three-vertex blocks, the action is
(j,l)↦(j+a,l-a). The laws uniform on j+l=0 and on j+l=1 have identical
one-port marginals, but occupy disjoint joint orbits. Separate marginal
weights cannot distinguish them. This is the inherited mod3 charge in a
source-weighted form; no independence assumption repairs it.

## 3. Arithmetic weights from the actual Berggren tree

Use Euclidean parameters (m,n) starting at(1,2), with branches

    M(m,n)=(n,2n+m), L(m,n)=(n,2n-m), R(m,n)=(m,2m+n).

Digit0 denotes M,1 denotes L,2 denotes R. The triple is
Psi(m,n)=(n²-m²,2mn,n²+m²). Thus every named vertex at depth k carries a
specific primitive triple. Choose weights from its odd leg, even leg or
hypotenuse and use (W4)-(W5). The tests compare the recurrence with literal
joint orbit enumeration for every ordered tuple of up to three ports under
C_1,C_2,C_3. They retain all original source weight and unequal heights.

The arithmetic mirror has an exact guard. The inherited THM-3357,
[Berggren Walsh/parent circuit](../../../01-canon/theorems/THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit.md),
gives H=[[-1,1],[1,1]], with H²=2I and

    HL=RH, HM=MH, HR=LH,       Psi(Hv)=2 swap_legs(Psi(v)).   (W6)

Therefore digit negation on words is intertwined with H only when the
root changes from(1,2) to H(1,2)=(1,3). At the fixed primitive root,
the exchanged first-level L/R hypotenuses are13 and17, so this is not a
height-preserving symmetry there. Applying H twice scales parameters by2
and triples by4; primitive normalization must retain that scale sidecar.
The code checks (W6) for all words of depth at most5.

This is a map of branch addresses with an arithmetic frame and scale. The
branch matrices are not elements of order three, and the C3 action on tree
digits is not asserted to be one global Berggren matrix. The weights make
the difference explicit rather than hiding it in an equinumerosity.

## 4. What this adds to the concurrent cycle carrier

Session94's [cycle-index recurrence](../session94_recursive_ports_crt_clocks_20261029/ports/WEIGHTED_CYCLE_CARRIERS.md)
is needed for unordered attachments: fixed-point moments alone lose moving
cycles. Its regular C9 versus C3×C3 hostile has identical named-port moments
at every order, but10 versus12 three-subset orbits. Our recurrence answers
a different query: the mass of a specified **joint named-port orbit** under
arbitrary supplied product weights. Neither result subsumes the other.

The retained tuple/LCA data provide the shared interface. Its source weights,
joint alignment, immutable side tags, and the chosen operation menu must
survive each later refinement. The computation gives3651 active checks;
normal and optimized replay must agree. A polynomial cost bound independent
of the supplied weight table is not asserted.
