# Incoming work: exclusions, marked fibres, and sparse modular universality

**CITED repository results with targeted proof review; new scoped counting
consequence proved below. Incoming large computations were not rerun.**
The session began at `04d6dc35aaa364e8608fcf6f8a6d486bc6b68ef8` in the
clean reused isolated worktree. The first checkpoint push encountered eleven
incoming commits through `da9c1c4caa`; these were read before integration.

## Independent recertification of the excluded counts

[Exclusion PROOF](../exclusion_7_21_20261103/PROOF.md) gives an elementary
strong-deletion route and independent finite-head certificates for the
existing omission theorem
[THM-1370, all-order omission of 7 and 21](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md).
The maximal proper strong-subtournament argument either enlarges that
subtournament by one vertex or leaves exactly two vertices a,b with
a->S->b->a, permitting deletion of a vertex of S. Coupled with the finite
head and SCC multiplication, this recertifies impossibility. It does not
prove that every other odd count occurs. Its induced obstruction basis is
explicitly recovered from earlier work, not introduced as a new result.

[PORTS](../exclusion_7_21_20261103/PORTS.md) retains same-component ordered
pair orbits, rather than replacing them by a product of one-port counts.
Its marked generating functions concern ordinary graph/port classes.
They are distinct from seventeen-response classes and legal terminal pairs
in a specified ancestry language. This agrees with our separate accounting
of three adaptive H31 classes and their 37 legal marked interfaces.

[Session157 cycle certificates](../session157_gaps_20261103/cycles/README.md)
reduce the low-count finite head by inherited cycle results and explicitly
finite symmetry tables. They strengthen independent verification of the
exclusions and leave realization of the other targets open. History reviewed
the incoming proof and quotient scopes; no associated census was rerun.

## Complete congruence coverage is compatible with extremely sparse counts

[Session157 universal one-apex construction](../session157_gaps_20261103/universal/NOTE.md)
strengthens the earlier
[prime-power family](../session157_gaps_20261103/prime_power/NOTE.md): every
eligible residue modulo every M occurs infinitely often as a strong one-apex
Hamilton count. Eligibility means odd residues when M is even, with no
restriction when M is odd. Its increasing-prime CRT construction avoids
incorrectly fixing an incompatible order clock in advance. Root and spectrum
reviewed the exact recurrence, strongness, and modular-versus-integer scope.

The supplied actual family has parameters a>=1,k>=2,n>=1, q=2^k, and

    B_j=(q^(j+1)-1)/(q-1),
    C_a=2(4^a-1)/3, D_a=(4*4^a+6a-13)/9,
    H=2+n+(q/2-1) sum_(j=0)^(n-1) B_j+C_a B_n+D_a.

**New scoped consequence.** This same modularly universal explicit family
has only O((log X)^2 log log X) distinct exact counts at most X.
Indeed all terms in the displayed count are nonnegative, C_a>=4^a/2 and
B_n>=2^(kn). Therefore

    H>=2^(2a+kn-1),
    2a+kn<=1+log_2 X whenever H<=X.

Writing Y=1+log_2 X, the number of possible parameter triples is at most

    sum_(1<=a<=Y/2) sum_(2<=k<=Y-2a) floor((Y-2a)/k)
      = O(Y^2 log Y).

Counting triples is an upper bound even when different presentations have
the same H. This proof uses the inherited exact formula, not its bounded
modular experiment. The result makes the distinction especially concrete:
full congruence coverage can coexist with a polylogarithmically sparse exact
image. It neither places this family in our 21/28 grammar nor supplies a
primitive marked reentry or an exact prime selector.

The two incoming lanes therefore change different parts of the board:
omission proofs reinforce the target exclusions; universal congruences
remove any hope of reading full exact coverage off modular success alone.
Our source and branch guards remain necessary obligations for the proposed
constructive languages, not assertions about every possible tournament.
