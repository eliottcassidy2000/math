# Incoming work: exclusions, marked fibres, and sparse modular universality

**CITED repository results with targeted proof review; new scoped counting
consequence proved below. Incoming large computations were not rerun.**
The session began at `04d6dc35aaa364e8608fcf6f8a6d486bc6b68ef8` in the
clean reused isolated worktree. The first checkpoint push encountered eleven
incoming commits through `da9c1c4caa`; these were read before integration.
A second publication race brought seven further commits through
`fa0a643acebf35322cc42d8051a18c26bd7028a3`, reviewed below. The computational
inputs to our frozen jobs did not change in either integration.

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

## Further constructive operations and a literal cover-profile map

[Session158 cyclic blocks](../session158_creative_20261103/cyclic_blocks/NOTE.md)
and [carrier](../session158_creative_20261103/carrier/NOTE.md) consume ordered
nonempty path-cover counts C_k(T). These are obtained exactly from our
eighteen polynomials:

    Q_T(z)=sum_s z^(|s|) P_s(z),
    C_k(T)=k! [z^k] Q_T(z).

Every cover has |s| marked-endpoint components and the anonymous components
counted by P_s; the factor k! orders its components. Thus this aggregation
preserves the entire input needed for uniform cyclic-block gluing. It loses
the paired retained-port labels, so the eighteen polynomials must remain
available when a later attachment distinguishes those ports. This is an
explicit source/target map, not an analogy from similar dimensions.

The [two-bin extension](../session158_creative_20261103/carrier/extension.md)
uses f(T)=C_2(T)+2H(T), allowing either labelled path slot to be empty.
In our rooted notation, f(T)=H(T)+sum_v F_v(T): each cover with a nonempty
first slot has a unique endpoint v, and the remaining empty-first-slot
case contributes H(T). This f is different from one fixed-root F.
The pair (H,f) multiplies under ordinal sums. A cyclic wrapper with two
singleton blocks gives H_new=H+f and is strong, even when its input block
is weak. The supplied actual family has counts 5+2^(n+5), n>=0, including
133=7*19 and8197=7*1171. It proves neither prime-value infinitude nor entry
in our adaptive source grammar. It supplies a different operation to test
at the sparse balanced boundary.

[The reflection sandwich](../session158_creative_20261103/reflection/NOTE.md)
uses an actual source A and Hamilton path p, its native partner X=K_p(A),
and the endpoint a of p. Its portal variable r is exactly our rooted F(X,a),
so its odd reflected native endpoint has count

    H_new=F^2+2hF+(h-s)^2.

The supplied91/187 hostile retains h,s but changes F. This construction
uses two copies, a new central vertex and a specified reflected path; it
changes H and yields a strong self-converse endpoint. It is outside the
H-neutral physical-port selector and the balanced-renewal output language.
No two-port entry or exact-target selector follows without a further map.

[The binary bank](../session158_creative_20261103/binary_bank/NOTE.md)
provides an exact inverse for its restricted strong transitive-deletion
family and a47X/120+O(sqrt(X)) upper bound on distinct counts. This does not
contradict the stronger polylogarithmic bound for the smaller Session157
modular family. Its [FRACTRAN program](../session158_creative_20261103/binary_bank/FRACTRAN.md)
is a terminating forward evaluator specialized to a supplied word, not
an all-target inverse or a universal primegame realization theorem.

[Square-overlap gadgets](../session158_creative_20261103/square_overlap/NOTE.md)
give an actual two-feedback escape from the one-apex restriction, including
H13, and a path-classified 2^A3^B response. The square-sum realization guard
is separate from that unrestricted recurrence; it cannot be omitted to
claim an all-B square-sum family. This continues the original prime-group/
square-sum lane through actual boundary attachments, without identifying
prime-group labels with the present tournament source states.

History, spectrum and toothpick reviewed these operation scopes. Root
checked the literal cover aggregation and the distinction between f and F.
The incoming finite computations were not rerun, and none changes the open
status of the full prime and7p lanes.
