# Unbounded port transport, a strong tripling channel, and the remaining digit entry

2026-11-03. **PROVED scoped operator grammar, all-order zero-letter channel,
and conditional coverage criterion; FINITE-EXACT controls and obstructions.**
Reusable closure under both nonzero digits and the assertion that 7 and 21
are the only missing positive odd Hamiltonian path counts remain **OPEN**.

The main advance is an exact unbounded grammar for the information needed by
the repairs, together with an actual indefinitely reusable strong tripling
construction. These are different achievements: the first transports counts
and defects under legal operations; the second supplies one arithmetic letter.
Neither yet supplies a closed construction for both H->3H+2 and H->3H+4.

## 1. The state now keeps unbounded component depth

[transport/PROOF.md](transport/PROOF.md) records spanning directed path covers
by the **paired endpoints** of components at currently active named ports.
Endpoints outside that boundary become `*`. A polynomial variable z counts
components with both endpoints forgotten. At p active ports there are at most

    g_0=1, g_1=4, g_p=4g_(p-1)+2(p-1)g_(p-2)

polynomial coordinates: 1, 4, 18, 88, 460, 2544, ... . Their degree grows with the
graph. This retains the unbounded path-cover information that a fixed list of
endpoint totals would discard.

The grammar has exact integral rules for:

- births with arbitrary active-port incidences and one uniform orientation
  toward all forgotten vertices;
- reversal of an edge between active ports;
- forgetting an active port and taking the converse.

Named edge reversal is `I+J_reverse-J_forward`, where J joins two distinct
path components. All named joins commute, but mixed products must be retained.
Uniform births use polynomial derivatives to count choices of anonymous
components. Forgetting can increase the z degree. These operations are closed
on the whole polynomial family, not on a fixed component cutoff.

Actual vertex names, source words and native guards remain separate
certificates. A forgotten name cannot later be distinguished without a richer
reconstruction. Strongness is not inferred from an aggregate polynomial.
The coordinate bound is not a polynomial-time algorithm for initializing
the state of an arbitrary tournament or solving global spectrum coverage.

The [phase and hostile-control readouts](transport/PHASE_CONTROLS.md) retain
the named joint ports and recover the integer inputs to both F49/F343
decoders. Fresh small checks reproduce the actual 170807 gauges and the
literal 81/27 hostile's distinct responses to the same cut. Changed graphs
require fresh decoding; field coincidences do not identify physical sources.

## 2. The changing inner defect has exact transport

The inherited H512425->1537279 branch changes its actual B14 suffix as follows:

|Stage|H(B)|E_32(B,2)|Inner defect|
|---|---:|---:|---:|
|Before repair|46073|46073|0|
|Reverse4->33|40243|40243|0|
|Reverse32->22|73259|75961|2702|
|Reverse33->31|80319|83279|2960|

The full third mixed response of H and E32 is 396. It comes from a genuine
five-component input layer; discarding that layer gives the wrong answer.
The corresponding third mixed response of the defect cancels, and the proof
explains this cancellation using its endpoint constraints.

The [positive quadratic potential](DEFECT_POTENTIAL.md) makes the mechanism
explicit. For fixed deletion D=T-r and incoming set I at r,

    b_r(I)=sum_(i in I) start_i(D)
           +sum_({i,j} subset I) (J_D(i,j)+J_D(j,i)).

At the failed collar, both start terms vanish and the pair term is 890+2070.
This is a local supermodular potential on the root's incoming set, not a
decreasing global rank when the deletion graph also changes.

There is also a smaller closed **defect-only** subgrammar: with two inports
and a protected preceding vertex, changes to one port's star act affinely on
a fixed bank of pair responses. For the actual collar, the entire 1024-choice
star bank has defect at least 614. The full count still needs additional
coordinates, and this restricted bank does not repair the collar to zero.

Independently, [neutral/README.md](neutral/README.md) exhausts all 125672
suffix changes using at most three edges, all outer marks, and every
nonshrinking transitive-prefix length. Of these, 4324 strong suffixes regain
some zero-ready mark; **none preserves the full 1537279 count**. Larger repairs
and other attachment types remain possible. This justifies carrying the
defect rather than assuming a free local reset.

## 3. A genuine unbounded strong zero-letter construction

[construction/RESULT.md](construction/RESULT.md) retains an old mark w and
an auxiliary a. Define h=H(B), F=E_w(B,2), and two nonnegative entry defects
sigma and tau, determined by exact two- and three-slot responses.
Two literal births followed by moving the auxiliary to the newborn give

    (h,F,sigma,tau) -> (3h+2sigma, 3F+2tau, 0, 0).

Consequently the zero-defect class is closed under `(h,F)->(3h,3F)` for
arbitrarily many steps. Every step keeps every old vertex and edge; strong
input produces strong intermediate and final graphs. In the outer zero-ready
attachment this multiplies H by 3 while retaining the old root star and inport.

Both nonzero entry conditions have actual small witnesses:

|Entry defect|Marked suffix response|After the same two-birth operation|
|---|---|---|
|(sigma,tau)=(1,0)|(h,F)=(5,7)|(17,21)|
|(sigma,tau)=(2,0)|(h,F)=(9,21)|(31,63)|

These give all-clock +2 and +4 families respectively, on different input
states. They are not an iterable two-branch family. The inherited reset
obstruction, after taking converse, proves that a C2 image of a strong base
of order at least 4 has **no** auxiliary with start defect 1 or 2. Re-marking
alone therefore cannot rearm this same operator for either next digit.
The order 3 exception is explicitly retained.

The actual failed B14 suffix is outside even the zero-letter chart:
it has no start-ready marking (minimum start defect 616), as well as no
end-ready marking (minimum 2960). Entry requires an actual change of source,
not a switch between equivalent scalar descriptions.

## 4. Both uniform digits from one infinite entry family

The late incoming merge led to a concrete positive result. For every t >= 1,
the same strong source with transitive prefix TT_t and an all-ported TT3
suffix has

    H_t = -5+7*2^t.

Replacing the terminal suffix vertex by a directed triangle with that old
vertex as its only port gives a TT2 -> C3 suffix and

    H_t^(2) = -13+21*2^t = 3H_t+2.

Two strong births retain every old name and edge. The incoming alternating
TT5 port pattern gives the other output

    H_t^(4) = -11+21*2^t = 3H_t+4.

The physical root defects are 6 at the source, 16 at the +2 output, and 12
at the +4 output. See the [literal entry construction](transport/TWO_DIGIT_ENTRY.md).
This gives both letters on the same infinite source family, but neither
output has been proved to return to a class admitting both letters.
Native Hamilton-complement guards remain separate from adjacency strongness.

The [construction bridge](construction/INCOMING_BRIDGE.md) identifies the +2
operation as C2 with entry defects (1,0). The +4 operation changes an old
root incidence and escapes C2's nonnegative endpoint correction. Equal H
values alone would have hidden this distinction. It also proves a separate
rigidity result: a nontrivial strong zero-ready port with exactly one
Hamilton path starting there forces H=2^m+1; only 3 -> 9 triples within that
specified boundary.

A new [transitive-chart theorem](TRANSITIVE_CHART.md) explains why the triangle
matters. In a transitive suffix the port mask f determines its joint defect
j(f). Even allowing arbitrary new port masks and independent fixed output
clock shifts, a uniform return has the unique possible digit

    d(f) = 3j(f)-j(3f)-2.

Thus both digits cannot stay in those transitive charts. The new +2 output
leaves that class in the smallest possible strong component. This is a
proved boundary for that operation language, not an unrestricted obstruction.

The [singleton/triangle component inverse](COMPONENT_LANGUAGE.md) then tests
the first larger candidate language completely, with no arbitrary word-length
cutoff. It accepts both first digits, only 22 and 42 at level two, and only
422 at level three, at the common prefix clock. The missing second-level
intercepts cannot be supplied by any word in that atom language. These are
typed acceptance results, not physical return programs or globally forbidden
Hamilton counts. The next concrete extension is the full port bank of the
strong four-vertex atom, whose response ratios break the integral mixed-radix
restriction behind this exact inverse.

## 5. Global coverage has an exact arithmetic target

[coverage](coverage/COVERAGE.md) separates source-compatible construction from
scalar arithmetic. On positive odd integers the three maps

    h -> 3h, 3h+2, 3h+4

form a ternary tree under the address m=(h-1)/2. Removing the forbidden values
7 and 21 leaves exactly six components, rooted at

    1,23,25,63,65,67.

A grammar accepting actual representatives for these entries and supplying
all admissible three-letter transitions, with its guards and source types
preserved, would prove the desired coverage. This is a **CONDITIONAL proof
criterion**, not a discharged acceptance lemma. Known multiplication-based
prime/7p reductions remain available with their own source requirements.

Even an indefinitely closed +2/+4 tree from 170807 alone would be too sparse.
Finite unions of such fixed-seed binary trees have sublinear power growth
and miss infinitely many primes and 7p values. The incoming binary-word reader
has a different, unbounded output language and is not ruled out by this bound.

## 6. Inheritance, incoming work, and the next proof obligation

The anchor was the actual failed inner collar; the niche was a growing
path-cover interface; the wildcard was the third arithmetic letter and its
entry type. The live concept board was: named source / paired components /
polynomial depth / positive defect / arithmetic letter / source acceptance.

The main recovered mechanisms are the older active-port pairing grammar,
[THM-3121 path-cover substitution](../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md),
[THM-3134 endpoint jets](../../../01-canon/theorems/THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform.md), the weighted-cut
defect, and the dual two-birth reset obstruction. Every connection in the
component notes states which operations its coordinates support and what
extra guards remain. No numerical coincidence is used as a theorem transfer.

Incoming Session181 independently supplies the unbounded attachment
`M(A;B,W)` with outer defect `H(A) J_W`, and a distinct 17-vertex +4 example
whose outer defect is 3690. Concurrent component-reader work supplies three
finite branching rounds with a complete word inverse. The final incoming
checkpoint refutes unlimited closure of its fixed 219-key bank at 13835447,
then repairs that actual miss by seed enlargement. Acceptance for a growing
recursive selector remains open. These are complementary interfaces; equal
H values do not permit identifying their named sources.

The later [Session182 integration](INCOMING.md) adds actual continuations of
that distinct defect-carrying source and the infinite clock-uniform +4 family
used above, with changing port membership. Its seven-coordinate
ordered-component law is a useful small quotient of our carrier for that
specific operation. It does not supply the missing two-branch closure.

The clearest sufficient next obligation for the C2 route is a **named digit
preparation operation**: preserve h and F while producing an active auxiliary
with `(sigma,tau)=(1,0)` or `(2,0)`, then apply the proved two-birth map.
Such preparations must return to an explicitly defined source class and work
at arbitrary depth. More general multiport grammars may change the attachment
as well, so this sufficient criterion is not asserted necessary for every
possible proof. The polynomial transport now supplies exact semantics for
testing those operations; the preparations and global entry coverage remain
the unresolved mathematical content. See [AUDIT.md](AUDIT.md) for the
independent controls and their exact scopes.
