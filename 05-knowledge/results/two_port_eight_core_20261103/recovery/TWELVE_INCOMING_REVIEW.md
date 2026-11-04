# Independent review of incoming twelve-vertex principles

**ACCEPTED: T1–T4 in
[two_port_twelve_principles_20261103/PROOF.md](../../two_port_twelve_principles_20261103/PROOF.md).**
No mathematical correction was required. The incoming proof and executable
were first reviewed directly from commit `f642c25f32`, without changing the
then-dirty worktree. After the clean checkpoint and rebase onto `7a64b40972`,
their local bytes were checked to be identical to those reviewed blobs.
This review separates the all-order arguments from their finite controls.

## T1: the full modulo-four descent is valid

For two nonempty modules on the same number of named vertices, the ordered
k-path-cover difference is divisible by 2k!: its unordered-partition formula
is k! times a sum of differences of odd products of induced Hamilton counts.
The same-size hypothesis supplies the same partition universe. Every k>=2
term in the integer uniform-run expansion consequently disappears modulo
four. The remaining coefficient is the odd H(Q), which is a unit modulo
four. No exact integer cancellation or primitive-content assumption is used.

The port step preserves this argument. Give each new port the same incidences
to all old exterior vertices as the contracted module vertex z. The old
module together with the new ports is then a module in the actual extended
tournament. Its internal attachments remain arbitrary, as does the fixed
tournament orientation among the ports. The quotient is still Q. Thus an
equality of complete ambient m-port responses forces equality of the core
m-port responses modulo four, for every m>=0. This is a necessary modular
descent for arbitrary same-size module pairs, not an exact response quotient.

## T2 and T3: paired paths and support parity are retained

The odd-cycle congruence imported from THM-466(i) is valid. Its maintained
correction retracts the reverse claim that residues recover individual cycle
collection counts; it does not affect H=1+2o modulo four used here.

The definition of delta divides the even integer difference of the two
Hamilton counts before reducing modulo two. It does not divide in F2.
Finite differences commute with that integer division, so its Boolean
coefficients have the stated meaning.

The decomposition of an odd cycle into its zero-port, one-port, adjacent-port,
and separated-port cases has the required multiplicities. In the last case,
the labelled ports canonically distinguish the two nonempty old paths:
the path after x and the path after y. Squarefree multiplication keeps their
supports disjoint and their endpoint pairings intact. There is no extra
factorial and no substitution of endpoint marginals for the paired paths.
The algebraic expansion using e=Us, G+G^t=ss^t, U^2=1, and p^2=q^2=0 gives
the displayed cubic expression. L_even and L_odd are linear support-parity
functionals; treating them as evaluation of the squarefree variables at one
would invalidate the proof, and the incoming proof does not do so.

For the shared-label coefficient a_u a_v c_v, u!=v, the surviving paired
path in X is exactly its singleton entry G_vv=z_v. The alternative endpoint
assignment repeats a vertex across disjoint supports and vanishes. Removing
v yields the endpoint-parity sum on even supports of T-v. At even n, the
superset count cancels precisely the supports of size n-2, leaving the even
sizes at most n-4. In particular the twelve-vertex sizes are 2,4,6,8; the
empty support contributes zero. The empty sums at n=2 and n=4 are consistent.
The second monomial a_v c_u c_v follows from the same paired-path symmetry.

T3 also has the correct sign and type. Each ordered-cover gradient is
independent of the mutual endpoint edge because its endpoints lie in different
paths. Flipping that edge therefore changes H by the reverse gradient minus
the forward gradient. Removing the changed edge from each odd cycle created
or destroyed gives a simple path between those endpoints with an odd number
of vertices. Such a path does not use the mutual edge. After integer division
by two and reduction modulo two, the signs disappear. This is an equivalent
path-parity implementation of the old C-skew modulo-four condition, not an
additional exclusion of its survivor set.

## T4: the 1,753-query count is for exact core/converse comparison

Reversing an extended Hamilton word and exchanging the port names gives
F_(T^op)(a,c)=F_T(1-c,1-a). Every quartic term comes from words in which
both ports have two old neighbors, so those words do not use the port-port
edge. Exchanging the ports preserves that word family, making its homogeneous
quartic part symmetric in the two masks. The two degree-four parts cancel
over the integers, before any residue reduction.

Both pure quadratic families also cancel. With a=0, x is a universal source,
leaving the one-port response. With c=0, y is not a universal source; the
separate specialized squarefree calculation in the proof is necessary and
correct. In both cases the quadratic coefficient for two distinct old labels
is `-(C_uv+C_vu)`, which is invariant under converse. The constant contrast
is zero by Hamilton-path reversal.

The remaining possible supports are zero, 2n singletons, n^2 mixed pairs,
and 2n*binomial(n,2) mixed triples. Triples may share an old label across the
two ports; they still contain three distinct attachment variables and are
included. Their total is 1+2n+n^3, equal to 1,753 at n=12. Omitted pure-pair
values follow from the singleton contrasts because their quadratic and
constant contrasts vanish. Ordinary Boolean finite differences then recover
every remaining coefficient.

This compares a core with its entire converse. It is not a 1,753-query
theorem for arbitrary tournament pairs, nor an exact comparison theorem for
reversing only one module inside an unchanged larger ambient tournament.
T1 transfers the resulting core distinctions only modulo four. Each counted
assignment requires one Hamilton evaluation for each of the two cores.

## Replayed finite universe and hostile controls

The incoming checker was copied from its reviewed git blob to a temporary
directory outside the repository and executed there normally. Its parsed
output matched the committed receipt exactly. The local normal and optimized
receipts also agree. The declared universe is all 1,099 labelled cores of
orders one through five, plus the specified order-six and order-eight controls;
complete two-port truth tables through order four; and the stated finite
same-size module comparisons in all labelled quotients through order three.

The replay confirms:

- 42,582 response comparisons, including 612 literal-permutation comparisons;
- 21,386 ordered path-skew checks and 21,386 shared-label cubic entries;
- 22,976 actual Hamilton queries used to check selected cubic coefficients;
- 1,089 actual ambient module comparisons;
- 77 literal odd-cycle congruence controls;
- 16,932 exact Boolean coefficients from the complete small tables,
  including 2,040 nonzero exact cubic coefficients.

The order-six hostile has H=45, symmetric C modulo four, and twelve nonzero
shared-label cubic entries. It fails endpoint parity. The order-eight hostile
has H=249, balanced endpoint parity, and twenty-eight nonzero cubic entries;
it fails C symmetry modulo four. These show failures of the respective
single-condition shortcuts. They do not prove the cubic test logically
independent of the conjunction of the two older modular conditions, and
the incoming proof correctly makes no such assertion.

## Independent transfer to the complete twelve-vertex data

A separate calculation in this review used the previously audited complete
integer response tables for the named twelve-vertex core with rows

```
(3706,680,3603,1508,3914,1300,294,1397,519,1256,2370,1002)
```

and its converse. It computed all 6,241 contrast coefficients by finite
differences, confirmed every quartic and pure-quadratic cancellation, counted
the 1,753 permitted query supports, and independently checked all 132 ordered
shared-label coefficients against integer subset endpoint counts. Sixty of
those shared-label residues are nonzero. The integer contrast has 24 nonzero
linear, 132 nonzero mixed-quadratic, and 1,584 nonzero cubic coefficients.

The separately implemented [converse_transfer.py](converse_transfer.py) was
also source-reviewed. Its Boolean predecessor recurrences compute start
parities and odd-support directed-path parities without importing the incoming
checker. Its [normal receipt](converse_transfer.json) and
[optimized receipt](converse_transfer.optimized.json) agree. They additionally
verify reconstruction of all 6,241 table entries from the permitted contrast
coefficients, both shared-label monomial families (264 checks), and all 132
ordered odd-support path-skew identities. Those counts and its sixty nonzero
shared-label residues agree with the initial independent calculation.

This is one specified twelve-vertex core and its converse. Neither the
incoming small controls nor this transfer is a complete order-twelve census.
Zero of the cubic residue polynomial remains only necessary for an ambient
hidden converse; it does not certify an integer full-response collision.
Universal named two-port reconstruction and exclusion of every prime
twelve-vertex module remain open.

## Reviewed source identity

SHA256 values at this review:

| Object | SHA256 |
|---|---|
| Incoming PROOF.md | `06811fe45366352f924793be7dadddacef287607a380efc92043a1066f3056c0` |
| Incoming check.py | `8ca59e0834cc20c9f3ea390fcf1078bc2e7e052ee5776e7c5087ab1e0cd9437d` |
| Incoming check.json and check.optimized.json, LF-normalized | `075260cde03b22429acb1cf3309139c192c1c0a176aa4f851cf957a8de841bd0` |
| Same two receipts, current worktree bytes | `9679b3552057030d9cdb46784f707d3d089a4e53465e661139a250e3b38132b1` |
| Existing twelve-core response/normal.json | `425a80fe1c311461295e6c51ef089bec93cdb22a647a9d340d0e399d21ca6a88` |
| Independent converse_transfer.py | `1c51d0e8700ba981278d81d45059833a20a3660d71b7c267426c7797e33d8187` |
| Maintained THM-466 source | `2964ced6d038c6dd9f2b798cbebbf3d26459d051d7ff1316eceb1dcaf5f1521c` |

The incoming receipt byte difference from the initially reviewed commit is
only CRLF versus LF; their parsed contents are identical. No incoming proof,
checker, receipt, or previously frozen packet file was edited for this review.
