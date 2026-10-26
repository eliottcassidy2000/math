# Integration audit: intrinsic tape bits, native clocks, and signed defects

2026-10-26. **ACCEPTED in the stated scopes.** This independently reviews
the incoming Session65 checkpoint `bcf9b658eb` and connects it to this
packet. The incoming results and their existing audits retain their own
provenance; no novelty is claimed here. Elementary proofs were read and
checked independently. The additional small controls below are
**FINITE-EXACT**. Supersingular genus identities and Ogg's complete
classification remain **CITED arithmetic**, not new proofs or new censuses.

Incoming fixed sources:

- [Checkpoint README](../session65_cryptid_tournament_interfaces_20261026/README.md).
- [Machine macros and word encoding](../session65_cryptid_tournament_interfaces_20261026/session65-cryptids/REPORT.md), §§2–3;
  [its independent audit](../session65_cryptid_tournament_interfaces_20261026/session65-audit-cryptids/REPORT.md).
- [General signed-pair polynomial](../session65_cryptid_tournament_interfaces_20261026/session65-fifteens/GENERAL_PERMUTATION.md);
  [Frobenius specialization](../session65_cryptid_tournament_interfaces_20261026/session65-fifteens/REPORT.md), §§2–3;
  [its independent audit](../session65_cryptid_tournament_interfaces_20261026/session65-audit-fifteens/AUDIT.md).

The frozen incoming proof drafts contain historical audit-pending wording;
the final linked PASS audits resolve that state. This review does not alter
their sealed files. Closest proved mechanism: order-join metric uncrossing.
Hostile: equal path counts with different next native macros. Corrected
near miss: parity of a vertex permutation is not its minimum orientation
defect. Least-used sidecars: native traversal length and the exact pair
action. Board: word order / clone mergers / edit metric / clock / signed
pair cycles / Frobenius. Anchor is faithful arithmetic state; niche is a
sharp clock obstruction; wildcard is a quantitative action defect.

## 1. Bits can be intrinsic even when zero blocks merge

Let `G0=TT3`, and let `G1=C3` differ from G0 by reversing its long arc.
For a finite binary word w, put `E(w)=G_(w1) join ... join G_(wm)`, with
every earlier block beating every later block.

The canonical strong-component order recovers the word without bit colors.
Each one is a C3 component; each zero is three singleton components.
Every maximal intervening run of 3r singleton components decodes as r
zeros. No claim that an original zero block remains a maximal transitive
module is needed. A write can split or merge those modules while this
decoder remains exact. For example, the pair of vertices {2,3}, using
zero-based labels, is homogeneous in `E(00)=TT6`; in E(01), vertex 5
distinguishes the two. The original cuts are recoverable from triples of
singleton components, not from an assumed unchanged clone partition.

For words u,v of the same length,

    d_iso(E(u),E(v)) = Hamming(u,v).

This is inherited from [THM-2183, order-join metric product](../../../01-canon/theorems/THM-2183-order-join-is-an-exact-tournament-metric-product.md).
Its proof was inspected: swapping the images of two crossing vertices
strictly decreases the mismatch count, while no third-vertex pair cost
increases. Hence every minimizing bijection respects each equal-size
ordered cut. The factors here all have order three, so iteration applies.
The distance uses isomorphism only, without identifying complements.
[THM-2195](../../../01-canon/theorems/THM-2195-transitive-quotients-exactly-control-universal-substitution-products.md)
explicitly blocks generalizing this product rule to arbitrary cyclic
substitution quotients.

This adds something concrete to our [colored configuration carrier](machines/README.md):
the tape bits themselves can be stored in uncolored tournament structure.
The window origin, head cell, machine state, and transition table still
have to be retained. Neither representation removes unbounded tape size.
A native write becomes one selected internal arc flip; the head marker
updates separately. This is a controlled rewrite, not the tournament
metagraph's K_P operation or a black-detour certificate.

The exact scalar identities are

    c3(E(w)) = #ones(w),       H(E(w)) = 3^#ones(w).

No directed triangle crosses an order-join cut, and a Hamiltonian path
must traverse the factors in order. These prove both formulas. Thus
E(01) and E(10) both have `(c3,H)=(1,3)` but have distance two. In the
literal bouncer+bell table `1RB0LA_0RC1LA_1RD0RA_0LB1RB`, with head zero
and state A, the first input reaches tape 11, head zero in two steps;
the second reaches tape 00, head -1 in one step. Consequently `(n,H,c3)`
does not determine this next macro. This is a precise obstruction to
using those scalars as a deterministic state for this operator. It is
not an impossibility theorem for every A000568 recurrence.

## 2. One changed arc can have an unbounded native clock

Use exactly our irregular-bell table
`1RB0RC_1LC1RA_1RA1LD_0LC0LA`. Begin in C on the rightmost cell of
`0(01)^k`. Each literal pair of instructions C1,D0 crosses one 01 pair
leftward without changing either cell. After k such pairs, C0 changes
only the initial zero to one and enters A one cell right. Thus the
native duration is `2k+1`, while exactly one tape bit changes.

Encode both tapes in the same window, adding a blank padding cell for
the final head when k=0. The literal tournaments differ in one arc.
Their triangle counts differ by one, so they are nonisomorphic; their
unlabelled edit distance is therefore exactly one. For every k>=0 this
gives the same distance but an arbitrarily large clock. No finite bound
depending solely on edit distance can control all these section macros.

The quantifier is important: these are legal finite configurations of the
stated macro section. Reachability of every such configuration from the
designated blank tape is not proved or assumed. A bound retaining window
length, sweep length, or an additional blank-orbit invariant is not
excluded. The lesson for arithmetic excursion compression is specific:
paying an edit or defect count does not pay elapsed native time unless
a separate clock identity or bound survives the quotient.

## 3. Full signed pair actions measure more than permutation parity

For a vertex permutation sigma, decompose its action on unordered pairs
into orbits O. Let L_O be the orbit length and epsilon_O the sign of its
endpoint return: +1 when both endpoints return individually, -1 when
they swap. If D_sigma counts failed orientation transports, then

    sum_T z^D_sigma(T)
      = product_O ((1+z)^L_O + epsilon_O*(1-z)^L_O).

The proof is sound. On one pair-orbit, the parity of the failure bits
must equal its endpoint twist. Every pattern with that parity admits
exactly two edge assignments, chosen by the initial edge orientation.
Independent pair-orbits use disjoint edge variables, giving the product.
Here D_sigma measures the prescribed labelled action. It is not
`d_iso(T,sigma T)`: that unlabelled distance is always zero because sigma
only relabels T. The tape isometry and the action defect have different
types, even though both count orientation disagreements before quotienting.

Exactly one negative pair-orbit occurs per even vertex cycle: its
antipodal pairs. Endpoints in different cycles cannot swap, and in one
cycle only antipodal endpoints can swap on return. Consequently, if e
is the number of even vertex cycles and q the number of pair-orbits,

    min_T D_sigma(T)=e,
    number of minimizers = 2^q product_(even cycles C)(|C|/2).

Every defect has parity e. The minimum is not determined by permutation
sign: `(12)(34)` is even but has minimum two and polynomial
`16z^2+32z^4+16z^6`. A 4-cycle has polynomial `8z+48z^3+8z^5` and
minimum one. This quantifies the familiar odd-order automorphism
boundary; it does not introduce a new automorphism theorem or prove
that a defect can be repaired by permitted black metagraph moves.

For supersingular j-values S_p, Frobenius is an involution, with r_p
fixed points among s_p points. It therefore has `(s_p-r_p)/2`
transpositions. Combining the elementary minimum formula with
[Aricheta's cited genus identity](https://arxiv.org/abs/1809.07421) gives

    min_T D_Frobenius(T) = (s_p-r_p)/2 = g(X_0^+(p)).

The minimization is over **all** tournament orientations of S_p. It
does not select an arithmetic comparison relation on j-values. The
incoming independent arithmetic audit identifies at p=37 the three
points `{8,3+10t,3-10t}`, with `t^2=2`, and Frobenius interchanges the
last two. Its pair-action polynomial is `4z+4z^3`, giving minimum one.
Our additional test below independently checks that permutation's
polynomial; it does not repeat the field/Hasse calculation. The full
Ogg prime classification is inherited CITED arithmetic, not established
by either this control or a finite prime scan.

This connects the [historical atlas](atlas/HISTORY_ATLAS.md)'s genus
defect to a precise optimization over binary relations. The arithmetic
input remains in the action and genus theorem. Forgetting those inputs
would lose the mechanism, even though its scalar minimum survives.

## 4. Bounded independent controls and exact limits

[check_incoming_s65.py](check_incoming_s65.py) imports only Python's
standard library, no incoming or existing repository engine. It checks:

- both irregular C-section exits for k=0,...,50: 102 native runs stopped
  on the first return to A, complete tape/head/state/clock and one-arc edit;
- all 63 binary words of lengths zero through five, after reversal of
  vertex labels, by a separately implemented SCC decoder;
- the clone-pair split at `00 -> 01`, and the complete 6! isomorphism
  comparison/path counts for the 01/10 scalar hostile;
- eight chosen permutation cycle types through order six, enumerating
  all 34,001 corresponding labelled tournaments and comparing full
  defect histograms, minima, parity, and minimizer counts to the formula.

This small universe is deliberately not a repeat of the incoming large
macro audit or its every-conjugacy-type census. The all-order claims use
the proofs above. Reproduce from the repository root:

    python 05-knowledge/results/root27_machine_bridges_20261026/check_incoming_s65.py
    python -O 05-knowledge/results/root27_machine_bridges_20261026/check_incoming_s65.py

Both modes pass with identical output. No new large census, morphology
classification, blank-orbit reachability theorem, formal-system build,
Collatz proof, or unrestricted black-core connectivity proof is claimed.
