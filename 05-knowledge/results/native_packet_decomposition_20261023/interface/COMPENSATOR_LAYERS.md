# A single active compensator gives translated charge layers

**PROVED, scoped interface; FINITE-EXACT native realization.** Assume the
independent binary U/V supply and complete physical-current criterion of
`../algebra/BOUNDED_CIRCUITS_AND_PRECEDENCE.md`. All full records, normalized
resources, pins, converse/inverse partners and kernel lifts remain supplied.
The group H below may include exact Z coordinates as well as finite factors.

Let O be the old source slots. Initially and finally exactly k of them use V,
with source-charge sum c in H. Let P be a disjoint supply of genuine helper
slots, each initially and finally at U. Permit at most one helper at V at a
time. Write lambda_i for the authenticated source charge of slot i.

The complete state space under this one-active restriction is the disjoint
union of the following layers:

    L_0 = {I subset O : |I|=k,   sum_I lambda_i=c},
    L_h = {I subset O : |I|=k-1, sum_I lambda_i=c-lambda_h}, h in P.

In L_h the helper h is at V and all other helpers are at U. In L_0 they are
all at U. Thus the helper supplies a **translated charge level**, together
with one unit of target margin; it is not an extra copy of an old path.

Between states (h,I) and (g,J), treating 0 as the no-helper label, the actual
changed-slot count is exactly

    |I symmetric_difference J| + d(h,g),

where d=0 if h=g, d=1 if exactly one label is0, and d=2 otherwise.
Replace each changed-slot indicator by its supplied path-clock cost for the
exact work metric. These formulas include switching one active helper for
another in the SAME atomic event; both old/new records are not simultaneously
occupied after arrival.

**Theorem.** The layered graph with this metric and the full eligible sets
L_h is isomorphic to the original one-active-helper assembly graph. In
particular, a path of maximum metric b is equivalent to a legal route using
at most b changed source slots per event under this restriction.

Proof. Actual target margin forces the stated cardinalities. The authenticated
additive current criterion forces the translated charge equations. Independence
supplies every remaining pin/resource condition. These steps are reversible,
so the state correspondence is bijective. Comparing old source bits and the
one-hot helper vector gives the displayed distance exactly. Retaining the
full-record sidecar reconstructs every native arrival. QED.

The unused helper slots are part of the supplied inventory throughout; a
layer name is not a license to allocate the same normalized port again. If
shared constraints couple sources, replace each L_h by its actual compatible
records and retain the coupling factors. The unconstrained subset description
then need not apply.

## Physical row accounting

If an exact character X has X_U=X_V=t, and the original old-source V sum is S,
then at a state with old sum S_I and a helper h the U-row changes are

    old rows:  2[(S_I-S)-t(|I|-k)],
    helper:    2[X_h-t].

They sum to0 when S_I+X_h=S and |I|=k-1. The V-row changes are their negatives;
every source-row change is zero because the two target values agree. The
factor2 is the original paired pure-black conductance. The diagonal degree
term t cannot be dropped when reporting a partial physical-row current. For
finite character coordinates the same identity is interpreted in the required
modulus with the appropriate paired factor. A scalar total does not identify
the owner at which cancellation occurs.

## Native consequence and what remains

For the original sixteen-source endpoints, k=8 and S=1016. The new
`../native/README.md` authenticates twelve possible helper sources. Every
single translated layer fails to connect the two endpoints with budget6.
The two layers belonging to sources1532 and1923 succeed together, using only
one helper at V at a time. Thus the source-pair interaction is geometric in
this finite layered graph, rather than a claim that either helper supplies
enough scalar charge by itself.

This interface turns the next question into a precise one: which actual
translated layers connect the components left by small legal moves? A proof
that sufficiently many suitable native layers always exist would be a genuine
all-size advance. The present theorem represents the problem exactly; the
finite eighteen-source construction solves one authenticated instance.
