# A finite parity quotient of every transitive-blowup family

2026-11-13, Codex. **PROVED IN NOTES**, derived from the inherited
run expansion and excursion-permutation lemma. This strengthens the
single-module halving obstruction in a different direction.

**Theorem.** Fix any tournament A and any set of its vertices to be
replaced by nonempty transitive tournaments of lengths m_1,...,m_t.
Modulo4, its Hamilton count depends only on which lengths equal1 and
which are at least2. No restriction is placed on the induced tournament
of module ports, the order of A, or the remaining lengths.

**Proof for one port.** Write H(A[r<-TT_m])=sum_j a_j K_j(TT_m).
Here K_j(TT_m)=j!S(m,j). For j>=4 this is divisible by4. For j=3,
a_3 is even by freely swapping the two internal excursions, so
a_3 K_3 is also divisible by4. Finally K_1=1 and K_2=2^m-2, which is2
modulo4 whenever m>=2. Therefore

    H(A[r<-TT_m]) = a_1+2a_2 (mod4),   m>=2.          (1)

For several ports, fix all lengths but one and regard all other
expanded blocks as part of the finite exterior frame. Apply(1), then
repeat for each port. This proves the theorem even when the module
ports contain directed cycles. It does not require a finite simultaneous
run expansion for all ports. QED.

For a frame with t chosen ports this is a finite response table on2^t
subsets, rather than an infinite table on positive length vectors.
The retained predicate is H modulo4. Exact counts, endpoint information
and most of the graph structure are lost. It is an exact quotient for
this one invariant, not a classification of the full tournament metagraph.

**Singleton resource law (PROVED).** Keep a fixed frame and fixed marked
ports, and only increase their transitive lengths. Every exact3H+4 step
must change at least one length from1 to at least2. Thus a trajectory
has at most t such steps, where t is the number of initially singleton
ports. Indeed, for odd H, 3H+4 is the opposite odd residue modulo4.
If no singleton status changes, the theorem forbids that change.
Each singleton can be consumed at most once under increasing lengths.
The law permits arbitrarily many3H+2 steps because those preserve H
modulo4. Rechoosing modules or changing the core lies outside its scope
and requires a new physical preparation argument.

**Corollary (PROVED).** No fixed finite frame with finitely many
transitive module lengths, each eventually either1 or at least2, can
realize

    Q_m = 9*3^m + 9*2^m - 6

for all sufficiently large m. In particular, no choice of positive
integer affine length functions q_i m+c_i with fixed coefficients can
do so. This holds with any finite number of modules, not only one.

Indeed Q_m is3 modulo4 at even m>=2 and1 modulo4 at odd m>=3, whereas
the theorem makes the proposed blowup's residue eventually constant.
The affine-length case allows q_i=0 with positive fixed length and
q_i>0 with eventual positive lengths.

The boundary matters. A length that alternates between1 and2 can switch
the residue class. Different frames for different target residues are
also allowed. Thus the corollary does not exclude a piecewise two-module
constructor, a nontransitive growing module, or an unbounded core
selector. Restricting to even or odd m removes this particular residue
obstruction, although the single-module excursion argument still rules
out an entire arithmetic progression in that more restrictive grammar.

This further sharpens the Collatz test: **adding a second steadily
growing transitive module alone cannot supply the missing first halving**
of the new plus-branch family. A candidate must change a module's
singleton status, its internal type, or the fixed frame, and must supply
the exact count after doing so.

`clone_parity.py` independently checks every subset of doubled vertices
of all ordinary frames through5, increasing each doubled vertex's length
to3 in turn. It also checks simultaneous blowups of the cyclic three-port
frame. The proof above, rather than that finite check, supplies all orders
and arbitrary module lengths.
