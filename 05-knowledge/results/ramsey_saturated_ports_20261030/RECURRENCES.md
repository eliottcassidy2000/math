# Exact growth with marked neighborhoods

**PROVED counting identities, classical in form.** These connect the
Ramsey attachment problem to the repository's tournament enumeration
work without identifying different categories of objects.

Fix targets `(k_1,...,k_r)`. For an avoiding coloring C on n named
vertices, let L(C) be the set of functions `f:V(C)->{1,...,r}` such that
the color-i graph on `f^-1(i)` has no `K_(k_i-1)`, for every i.

Adding a new vertex with incident colors f is legal **if and only if**
`f in L(C)`. An old forbidden clique was already absent, and a new
forbidden clique consists of the new vertex and exactly such a clique
in one of its ports. The predicate uses the induced colored subgraphs,
not just the cardinalities of the ports.

If A_n counts labelled avoiding colorings, then

    A_(n+1) = sum_(C labelled on n vertices) |L(C)|.

Writing M_n for the sum of `1/|Aut(C)|` over unlabelled C gives

    M_(n+1) = 1/(n+1) sum_[C] |L(C)|/|Aut(C)|.

For unweighted orbit counts there is a distinct identity:

    sum_[C on n] |Aut(C) \ L(C)|
      = sum_[D on n+1] |Aut(D) \ V(D)|.

An attachment orbit is a child coloring with its added vertex marked.
Deleting the mark does not yield one copy of each child; it yields one
copy per vertex orbit. This is the same mark-retention issue that
appears in canonical tournament generation, now with an r-color edge
alphabet and different legal-attachment predicates. Color permutations
are not quotiented in these formulas.

`recurrences.py` checks complete labelled universes through n=5 for
targets `(3,3,3)`, `(3,3,5)`, `(3,3,6)`, and both orbit identities for
growth from n=3 to4 using all vertex permutations. The two latter
universes first differ at n=5, by the single all-blue coloring.
This finite audit is a correctness control, not a competitive Ramsey
enumerator or an improvement to the known bounds.

The relevant next compression is operation-dependent. For the
17-vertex `(3,6)` ports, seven unrooted graphs yield67 rooted types.
Two simultaneously selected vertices require pair orbits, rather than
an independent product of single-root orbits. A future transfer matrix
must retain those joint incidences until an equivalence under every
allowed extension has been proved.
