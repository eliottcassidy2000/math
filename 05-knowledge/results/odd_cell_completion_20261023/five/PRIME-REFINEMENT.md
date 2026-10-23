# Prime regular rotation groups force proper equitable partitions to split

23 October 2026. **PROVED**, with exact finite controls. This strengthens the
five-cell one-owner propagation theorem to every prime cell whose internal
automorphism group is exactly a regular C_p. It does not include the Paley
seven-cell, whose internal group has order21, or general prime tournaments
without a regular cyclic automorphism group.

## Theorem

Let T be a tournament on p vertices, p an odd prime, and suppose
Aut(T)=C_p acts regularly. Then its only equitable partitions are the
indiscrete partition and the partition into singletons. Consequently any
proper initial colouring becomes discrete under directed equitable
refinement, and one actual owner suffices to individualize T.

**Proof by rational invariant spaces.** Index the vertices by F_p so the
generator is translation. There is a connection set S in F_p^* with exactly
one of x,-x present and adjacency A(i,j)=1[j-i in S]. Write zeta for a
primitive p-th root of unity. A has constant eigenvalue (p-1)/2 on the
constant vector; its other eigenvalues are

    lambda_k = sum_(s in S) zeta^(ks),  k=1,...,p-1.

These p-1 numbers are distinct. If lambda_k=lambda_l, the integer polynomial

    sum_(t in kS) X^t - sum_(t in lS) X^t

vanishes at zeta. Its degree is at most p-1 and its constant coefficient is
zero. The minimal polynomial 1+X+...+X^(p-1) has nonzero constant
coefficient, so the displayed polynomial must be identically zero. Thus
kS=lS. Multiplication by l/k is an automorphism fixing vertex0. Since a
regular C_p has no nonidentity point stabilizer, l/k=1, so k=l.

The Galois conjugates of lambda_1 therefore have size p-1, giving an
irreducible minimal polynomial of degree p-1. The rational zero-sum space
U={v in Q^p: sum v_i=0} has dimension p-1 and is invariant under A because
the tournament is regular. Hence U has no nonzero proper rational
A-invariant subspace: for any nonzero u in U, its cyclic span has minimal
polynomial dividing that irreducible polynomial, so it has dimension p-1.

For an equitable partition, the rational span W of its cell-indicator
vectors is A-invariant and contains the all-one vector. If the partition
has at least two cells, subtracting a suitable rational constant from one
cell indicator produces a nonzero vector in W intersect U. Invariance
forces U to be contained in W; together with the constant vector this
gives W=Q^p. Its dimension is the number of cells, so every cell is a
singleton. QED.

This proof uses only the cyclotomic polynomial and the actual automorphism
group. It does not infer rigidity from a generic simple spectrum. In
particular, the rationality and irreducibility of the nonconstant spectral
piece are load-bearing.

## Consequence for attachments at every total size

Take a supplied named equitable collection of such prime cells, with
matching corresponding cells on a target when chart comparison is needed.
A nonuniform equitable block between two cells of different prime sizes
cannot exist: if their sizes are p,q and row/column sums d,e, then
p d=q e, so coprimality forces a zero or full block. Thus each connected
nonuniform support component in this carrier has a common prime size.

One actual owner makes its root cell discrete by the theorem. Once a
neighbor's source-side owners are discrete, a nonuniform block gives a
proper colouring of its target cell. The theorem makes that target
discrete too. Induction proves that **one owner per connected nonuniform
component suffices for ordinary refinement at every total tournament
order**, not only for regular-five cells.

The exact zero-or-one external-owner saving in the main five-cell note
also extends in principle if full block canonicalization is retained.
Fix an internal isotype template and one generator of its C_p, carrying
that generator through every comparison chart. This transported generator
is independent of the selected chart because C_p is abelian. Arbitrarily
changing the generator independently on an incident block is not allowed.
For any rigid block its C_p x C_p rotation orbit has no stabilizer, so
the lexicographically least cyclic matrix realization gives a unique
actual root pair. For a symmetric block B(i,j)=1[j-s i in S], the proper
nonempty target subset S has trivial translation stabilizer, so its least
cyclic binary word gives a unique actual selected target owner in each
row, producing a bijection of slope s. Nonunit cycle gain then selects a
unique fixed owner. This is an exact finite algorithm using the supplied
actual blocks, not a bounded-size motif observer or fast enumeration of
unrestricted tournaments.

At five vertices the main note supplies much smaller actual pair flags,
rather than requiring general block canonicalization. The existence of
the prime rule does not make the seven-cell theorem redundant: it handles
all regular-seven types, including Aut=C3 and Aut=C7 semidirect C3.

## Exact controls and boundaries

`check_prime_refinement.py` tests every two-colour split up to swapping
the colours in every trivial-multiplier Cayley connection set for primes
3,5,7, and in the standard cyclic R11 and R13 representatives. This is
finite corroboration of the proof, not its justification. The p11/p13
controls do not enumerate all connection sets or all tournaments.

The standard cyclic R_p has automorphism group exactly C_p: for an arc
i->i+d, 1<=d<=(p-1)/2, the common-outneighbor count is (p-1)/2-d. Hence
the directed short successor i->i+1 is intrinsic, and its automorphisms
are just rotations. This is the mechanism already recorded in the
incoming prime-cell Burnside packet.

The spectral obstruction also explains why the field/set distinction
matters. A nontrivial multiplier stabilizing S makes several lambda_k
coincide, so the irreducibility argument no longer applies. This is an
actual information loss, not merely a numerical coincidence of primes.

    python -B 05-knowledge/results/odd_cell_completion_20261023/five/check_prime_refinement.py
    python -O -B 05-knowledge/results/odd_cell_completion_20261023/five/check_prime_refinement.py
