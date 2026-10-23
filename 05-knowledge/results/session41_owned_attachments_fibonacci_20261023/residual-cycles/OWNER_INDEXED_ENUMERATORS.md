# Owner-indexed attachment enumerators

**PROVED, all finite forests and supplied tree decompositions.** This is a
standard finite-state contraction mechanism applied to the shared-frame
carrier of README.md. It gives a useful exact generating function without
pretending that a large cyclic support is a forest.

## 1. The count must retain the shared owner until contraction

A local attachment state a at cell i belongs to a finite domain D_i.
Allowed pairs R_ij retain the actual shared charts and owners. Assign
integer vector observations p_i(a) and p_ij(a,b), and let X^p denote the
corresponding Laurent monomial. Negative exponents are harmless because
all state spaces are finite.

For a forest, root each component. Define

    Z_i(a) = X^(p_i(a)) times
             product over children j of
             [sum over b with (a,b) in R_ij of
              X^(p_ij(a,b)) Z_j(b)].

At a leaf the empty product is one. Then

    Z(X) = product over roots r of [sum over a in D_r of Z_r(a)]

is exactly the generating polynomial of all globally compatible assignments,
with a term for each assignment and exponent its total observation.

**Proof.** With the root state fixed, removing the root separates child
subtrees, and their only common data are the already fixed root state.
The allowed edge pair accounts for the boundary compatibility. By induction
each child message counts every compatible subtree assignment once. Product
joins independent child subtrees; the outer sum forgets the root only after
all its incidences have been checked. Components are independent. QED.

In particular Z(1) counts compatible states. If Z is nonzero, the observed
wiring is unique exactly when Z has one distinct exponent, regardless of
its coefficient. Zero Z means infeasible, not unique by an empty-set
convention. The actual observed wiring may have many frame realizations.

## 2. A fast certificate without expanding the polynomial

For each scalar observation coordinate separately, replace every sum in
the recursion by minimum or maximum, product by addition, and monomial by
its exponent. Infeasible options carry +infinity for minima or -infinity
for maxima. The same induction gives the exact minimum and maximum over
all compatible assignments. A nonempty fiber has a constant vector
observation iff every coordinate's minimum equals its maximum.

The work is proportional to the number of retained allowed pairs and local
states per observation coordinate, with ordinary arithmetic-cost caveats.
Store the messages and predecessor choices: checking each recurrence gives
a finite certificate; differing extrema give two actual compatible
assignments. This is an efficient all-size positive class, not a complexity
claim for unrestricted cyclic supports.

## 3. Why multiplication recycles an owner

Two local switches have state t in {0,1}. Independently their switch-count
polynomials are 1+z and 1+w. If the port is shared, the constraint is that
the **same** t is used twice. Correct gluing yields

    sum_t z^t w^t = 1+zw,

whereas ordinary multiplication yields 1+z+w+zw and introduces the two
incompatible mixed choices. Setting w=z in the product gives 1+2z+z^2,
not the correct 1+z^2. Variable identification after forgetting the owner
cannot restore the lost compatibility.

Equivalently retain a polynomial with an owner index until its coefficient
arrays are contracted; for this equality constraint the operation is a
coefficientwise match (a Hadamard-type contraction), not multiplication.
This is exactly the joint-switch hostile in README.md, now seen through
its enumeration rather than its residual kernel.

The map to a count-only polynomial loses which owner contributed a term.
The restoration sidecar is the common owner index. The two-state example
is the cheapest decisive hostile; equal cardinalities do not repair it.

## 4. Extension to actual separator decompositions

For a supplied tree decomposition, keep a table on every full separator
assignment. A bag is processed by multiplying its own factors and child
messages, then summing over variables not in the parent separator. Assign
each local factor to exactly one bag. The running-intersection property
(each variable's containing bags form a connected subtree) guarantees that
all occurrences of a shared owner are identified before it is forgotten.
The same induction proves exactness, including observation enumerators.

Thus a boundary compression is valid when its retained boundary relation
and observation data determine every interaction with the rest of the
graph. Marginal cardinality alone does not. The incoming 64-state theorem
is a special relation compression for a particular three-boundary affine
carrier; it does not supply a three-boundary decomposition of the actual
order-31 K8 triple support. A single bag containing all eight cells is a
correct but potentially expensive fallback. A smaller decomposition must
be exhibited and checked, not inferred from the small boundary alphabet.

## 5. Counting frames is not counting tournament classes

For the frame-comparison carrier, Z(1) counts marked cell-preserving
isomorphisms, not isomorphism classes. If an isomorphism f:A->B exists,
composition f h gives a bijection from the retained automorphism group of
A to all such isomorphisms. Thus the count is zero or |Aut(A)|. This also
explains why two empty parallel-section summaries can conceal frame counts
three and one in the order-31 triple mechanism.

For a family of candidates, any division by a group order requires a
specified common action, a complete invariant inventory, and the correct
weights. In particular, for a finite G-invariant set X, the valid formula
is sum over x in X of |Stab_G(x)| divided by |G|; varying stabilizers do
not invalidate that formula. Arbitrary comparison-frame totals are not
automatically that stabilizer sum. The observation and orbit quotient
must be declared separately. Likewise, a polynomial over
abstract compatible states is not automatically the native metagraph
Hamiltonian-port enumerator: the actual path inventory and all target and
converse constraints must be included first.

## Controls

The companion check enumerates every binary pair-relation assignment on a
four-vertex star (16^3=4096 instances), computes the owner-indexed polynomial
by the recurrence, and compares every coefficient with brute-force global
assignments. It also verifies the two-switch contraction hostile. This is
an independent finite check of the proof's common-owner bookkeeping.

## 6. Zero-length slots: elimination is not reconnection

The Fibonacci-color proposal supplies a useful exact analogy, not an
identification with tournaments. Consider binary hard-core variables on
A--Z--B, where adjacent selected variables are forbidden. If the zero slot
Z is required to be unselected, its exact eliminated factor on (A,B) is

    sum over z=0 of 1[not(A and z)] 1[not(z and B)] = 1.

All four assignments of (A,B) remain possible. Deleting Z and inserting the
edge A--B instead adds a new prohibition and loses the (1,1) assignment.
There is no cycle here: even on a tree, deleting a zero-weight vertex is
not automatically constraint-preserving elimination.

If instead Z is a zero-weight state allowed to be selected or unselected,
its eliminated multiplicity table on 00,01,10,11 is (2,1,1,1), not (1,1,1,0).
Thus the numerical weight zero does not decide its state-space convention.
One must say whether Z is forced absent or whether both zero-length choices
are counted. In both cases the correct contraction retains the induced
factor, not a guessed replacement edge.

## Coupling can also rescue uniqueness

In the two-switch fiber of README.md, X00+Y01 is identically one. Switching
X alone changes this observation by minus one; switching Y alone changes
it by plus one. Neither is allowed, and the allowed simultaneous switch
has zero total change. Thus an uncoupled zero-cycle-cost potential can fail
while the coupled observation is unique. In the same fiber X00 alone is
not unique. The observer matters as much as the compatible move bank.
