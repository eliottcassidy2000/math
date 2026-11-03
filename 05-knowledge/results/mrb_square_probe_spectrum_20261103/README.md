# MRB phase filters and a square-sum probe bank for actual tournament sources

**PROVED scoped synthesis / FINITE-EXACT controls. Same-agent proof review,
with separate computational implementations. The full odd Hamilton-path
spectrum remains OPEN.** No new numerical coverage record or priority claim.

This packet connects two older mechanisms to the current source-selection
problem. A square-sum bank implements every one- and two-port birth on the
same scaled source. Those probes determine the *entire* one-birth count
response, with one pair probe recoverable for free. Period-four MRB filters
give a nonzero, rigorously signed test for the two named port phases of the
actual H31 branching source. Neither statement supplies recursive coverage.

## 1. Inheritance and the live board

Base snapshot: `d79221bc31`. Counts throughout are raw directed Hamilton-path
counts H, never H/Aut or square-compatible path counts.

- Closest proved mechanisms: [square-scaled simultaneous pair extension](../induced_square_pair_extension_20261019/SCALED_PAIR_EXTENSION_AND_TWO_DEGENERATE.md),
  and **THM-4104**, [selected-order-eleven-strong-ear-solid-interval](../../../01-canon/theorems/THM-4104-selected-order-eleven-strong-ear-solid-interval.md),
  §1's exact quadratic ear response. The pair-extension theorem already
  embeds all finite two-degenerate graphs. Neither theorem is new here.
- Current actual-source input: [BRANCH B4](../source_selection_closure_20261103/BRANCH.md)
  gives one named H31 source with branches 93,95,97 and marked reentries.
  Its complete next-level proper-cut closure fails. We transport its literal
  witnesses, not a scalar recurrence pretending to choose its children.
- Canonical hostile: [attachment theorem](../prime_lane_ports_five_groups_20261102/attachment.md)
  already shows that opposite named triangles have the same full one-birth
  response but different two-birth responses. This is recovered and retested.
- Corrected near miss: the earlier [MRB contour audit](../session31_mirror_mrb_20261022/mrb/MRB_CONTOUR_AUDIT.md)
  rejects the proposed exact equality D=q^4. The stronger
  [MRB research note, Theorems 1–2](../session144_jacobian_eth_mrb_20261103/mrb/RESEARCH.md)
  preserves finite prefixes and all half-plane algebraic asymptotics while
  changing the alternating functional, and classifies periodic filters by
  their mean. Similarity of summaries does not determine a response.
- Least-used relevant sidecar: the common square scale and the separation
  between the *vertex-name order*, numeric labels, and the *query clock*.
  [Shared square clocks](../attachment_lattices_square_sum_20261023/square_sum/SHARED_CLOCKS.md)
  supplies the older warning against replacing indexed clocks by residues.

| Lane | Object and operation | Preserved target / decisive boundary |
|---|---|---|
| Anchor | Actual H31 source and its two-port births | All named arcs, raw H, physical ports; one-level branching is not closure |
| Niche | Common square-scaled pair bank | Full old square/non-square relation; original numeric labels change once |
| Niche | Quadratic ear field | All one-birth counts from sparse probes; not a recursive state |
| Wildcard | Period-four MRB filters | Exact signed port-phase distinction; no MRB closed form or coverage implication |
| Hostile | Triple anchors 1,2,3 | Three-neighbor extension can fail after every square scaling |

The global {7,21} exclusions are inherited. The remaining prime/7p
realization reduction is **THM-4094**,
[hamiltonian-matching-deficit-and-two-prime-lane-completeness](../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md).
No finite test here changes that open obligation.

## 2. One common arithmetic bank for every small cut

**Proposition P1 (square-controlled birth bank).** Let b_1,...,b_n be distinct
positive integers, H=max b_i, and let T be an actual tournament on these
named vertices. There is one common integer-square scale L and distinct new
labels X_S for every S subset of V with |S|<=2 such that:

1. Lb_i+X_S is a square exactly for i in S;
2. all X_S exceed all Lb_i, and no sum of two distinct X_S is a square;
3. scaling preserves every old square and nonsquare pair;
4. attaching X_S by i->X_S iff its sum with Lb_i is square gives exactly
   the ordinary tournament birth with incoming old neighborhood S.

The actual old arcs of T are retained. For arbitrary T they remain an
explicit sidecar: the numbers alone are not asserted to encode arbitrary
old orientations. If a supplied total order already encodes T by
`i<j: i->j iff b_i+b_j is square`, the entire one-birth chart obeys the same
rule. This is an intrinsic pairwise observable with an explicit order gauge,
not a claim that an undirected square graph is itself a tournament.

**Construction and proof.** Add a fictitious anchor b_*=H+1. Let lambda be
the least power of four at least 32(H+2), and L=4lambda. For a two-set S
use its labels a<b; for singleton S={i}, use a=b_i,b=H+1. Set

    X_S=(lambda(b-a)-1)^2-4lambda a.                 (1)

Apply the inherited simultaneous pair theorem to the augmented seed, then
delete the fictitious anchor. It gives every claimed incidence, collision
exclusion and old-edge preservation. All these X_S are 1 modulo4. For the
empty set take X_empty=LH+2: its sums with old labels are 2 modulo4, and
with all other candidates are 3 modulo4. It too exceeds all old labels.
Thus all candidates form an independent square-sum bank. All labels are
below the inherited loose bound 65536(H+1)^4. QED.

This adds the simultaneous singleton/empty alternatives and the actual
source-preserving birth interpretation to the older pair construction.
The alternatives may coexist as arithmetic probes; H(T+X_S) uses only the
selected candidate, not the full bank as additional tournament vertices.

**Sequential compiler.** Every finite birth word in which a newborn loses
to at most two old vertices has a faithful square-sum encoding with the
same named source and command history. Start with label1 and apply P1
inductively. For fewer rescalings, with current maximum M use

    X=(M+1)^2                  for no neighbor,
    X=(M+1)^2-b_i              for the sole neighbor i.   (2)

These exceed M. Any unintended square would differ from (M+1)^2 by a
nonzero integer of absolute value at most M, less than either adjacent
square gap. Both moves retain all existing numeric labels. A two-neighbor
move may rescale them; the scale is recorded. Finite construction does not
give a fixed infinite integer embedding under repeated rescaling.

## 3. A whole ear cube from the bank, with one probe omitted

Use incoming-set convention R_T(S)=H(T+x), where i->x iff i in S. This is
the complement of THM-4104's convention. Let s_i,e_i be source start/end
counts, and Q_ij count old permutations with i immediately before j and
every other step legal. The exposed step may have either orientation. Then

    R_T(S)=H + sum_(i in S) l_i - sum_({i,j} subset S) K_ij,
    l_i=e_i-s_i+sum_j Q_ij,       K_ij=Q_ij+Q_ji.        (3)

This is inherited: removing the newborn from a Hamilton path gives either
an endpoint path or a two-path cover, not necessarily an old Hamilton path.
Simply counting insertion slots in old Hamilton paths would omit those
covers. The two constant cuts give R(empty)=R(V)=H and hence

    sum_i l_i = sum_(i<j) K_ij.                         (4)

**Proposition P2 (arithmetic observation and one free pair).** Values on
the P1 bank determine every R_T(S). Explicitly,

    l_i=R({i})-H,
    K_ij=R({i})+R({j})-R({i,j})-H.                     (5)

For n>=2, omit any one pair {a,b} from the queried bank and recover its
K_ab from (4). Thus H, all n singleton probes and all but one pair probe
suffice: n(n+1)/2 scalar values in total. If H is already known, only
n(n+1)/2-1 new-vertex queries are required. This is a linear identity, not
an average-coverage assertion. Equivalently, for k=|S|,

    R(S)=sum_(pairs in S) R(pair)
          -(k-2)sum_(i in S)R({i})
          +(k-1)(k-2)H/2.                              (6)

Formula (6) also holds for k=0,1. Negative coefficients recover the count
of an explicitly defined actual cut; they do not turn a formal difference
of unrelated tournaments into one source.

The real space of unrestricted quadratic fields (3) satisfying (4) has
dimension n(n+1)/2. The omitted-pair bank is an invertible coordinate system
on that space. The exact finite audit finds the actual tournament row span
has that full dimension at n=4 and5 (10 and15), but ranks at n=1,2,3 are
only1,2,4. No all-n rank theorem or nonlinear/adaptive query lower bound is
claimed. The finite full-rank certificates are checked over Q separately.

The three-anchor obstruction is sharp information about P1's scope:
1,2,3 admit no common square neighbor after *any* square scaling, by the
inherited rational-right-triangle-area-one reduction. Nevertheless (6)
computes the three-port tournament count. Arithmetic observation needs only
pair probes; faithful arithmetic realization of a specified three-port
birth needs additional information or a different placement.

## 4. Actual shared labels for the H31 branching source

The inherited source B has outgoing rows

    (42,4,1,22,7,94,31),        ports (5,3),        H=31.

One faithful increasing square-sum labeling is

    (1,3,13,195,38221,1460921283,2134290996583287373).

After the single common scale L=295147905179352825856, choose the candidates
from P1. The four port masks have these counts:

| Incoming old set | empty | {5} | {3} | {5,3} |
|---|---:|---:|---:|---:|
| R | 31 | 53 | 97 | 95 |

The 95 and97 branches are single births. The 93 branch first uses {5}, then
uses its newborn7 as the sole incoming neighbor of newborn8, applying (2)
without further rescaling. Consequently all three complete branches retain
the *same seven old numeric labels*, not just isomorphic old graphs.
Their terminal physical ports are respectively (8,3), (5,7), (7,3), as in
the inherited native command words. Literal labels and all endpoint rows
are in `receipt.json`; their largest labels have 508,193,254 bits.

The values and marked reentries were already known; the common arithmetic
chart is the new transport. The 93 branch really needs two births here:
the complete 128-cut table of B does not contain93. The old next-level
proper-cut failure through H295 remains in force. Arithmetic relabeling
does not change which abstract continuation graphs exist.

## 5. An MRB filter that provably sees the missing port phase

Put f(n)=n^(1/n)-1. On the repeating four-query clock
`empty,{5},{3},{5,3}`, let

    alpha=(-1,+1,-1,+1),    beta=(-1,-1,+1,+1),
    gamma=alpha*beta=(+1,-1,-1,+1).

These three filters span the mean-zero period-four space. In particular
alpha_n=(-1)^n, so C_alpha=sum alpha_n f(n) is the usual MRB constant.
Define C_beta and C_gamma analogously. Their convergence follows from the
inherited periodic-mean theorem, or bounded partial sums and eventual
decrease of f. Linear independence of these *filters* is not an assertion
of arithmetic independence of the resulting constants.

For the actual four responses above the Walsh decomposition is exactly

    R=69+5alpha+27beta-6gamma.                           (7)

The interaction probe R00-R10-R01+R11=-24 is a mean-zero filter that also
annihilates both singleton terms. Symmetrizing under simultaneous port-bit
reversal retains 69-6gamma and loses the two directional terms. Swapping
the two named ports retains the mean and interaction but exchanges 5,27.

**Proposition P3 (a nonzero source-specific analytic selector).** Define

    Psi(B;5,3)=sum_(n>=1)(R_mask(n)-69)f(n).

The source and query clock are fixed. Then

    Psi(B;5,3)=5C_alpha+27C_beta-6C_gamma,
    Psi(B;5,3)-Psi(B;3,5)=22(C_beta-C_alpha)<-319/500.    (8)

**Proof of the strict sign, without an MRB decimal evaluation.** Pair terms:

    C_beta-C_alpha=2 sum_(j>=0)[f(4j+3)-f(4j+2)].        (9)

This difference series is absolutely convergent by eventual monotonicity
and telescoping domination. Every term after j=0 is negative, since
log(x)/x decreases for x>e. Integer-power-certified rational bounds for
the six roots at2,3,6,7,10,11 give

    2 sum_(j=0)^2[f(4j+3)-f(4j+2)] < -29/1000.

The remaining tail is negative. Multiplication by22 proves (8). The exact
root bounds are in the receipt and independently checked by integer powers.

There is an all-source version, which also gives the exact limit of this
selector. For any T and distinct ports a,b, center its four responses by
their mean and define Psi on the same clock. Directly subtracting the two
port orders gives

    Psi(T;a,b)-Psi(T;b,a)
      = [R_T({b})-R_T({a})](C_beta-C_alpha)/2.           (10)

Thus the MRB filter separates these port orders **if and only if** their
singleton birth counts differ, with the opposite sign. The analytic
functional carries precisely this existing directional information; it
does not discover a further hidden source invariant. In particular it
cannot resolve the identical-table triangle hostile below.

This is an exact connection through the same finite character basis and a
genuine nonzero distinction of the two port phases. It is not an MRB closed
form, a new irrationality theorem, or a source-coverage proof. The query
clock n is not the potentially huge square-sum vertex label. Replacing one
clock with the other would change the functional and lose this statement.

## 6. What remains to turn this into spectrum coverage

The arithmetic bank represents two-port choices faithfully, and the signed
filter detects a specific lost orientation. Neither makes one-step data
recursively closed. We replay the inherited minimal hostile:

    T=(2,4,1), U=(4,1,2),
    R_T=R_U=(3,5,5,5,5,5,5,3).

Births with incoming named sets {0}, then {1}, give13 from T and11 from U.
Therefore even all P1 probes and every linear MRB filter of that one-step
table miss a circulation that a later legal operation sees. Full graph data
or the appropriate joint continuation observer must remain available.

Likewise, a nontrivial native chord complement cannot be implemented just
by common square scaling with the same order gauge: scaling leaves every
old arithmetic orientation unchanged. A reset needs actual arc data and a
new arithmetic chart, or an explicitly proved chart-changing operation.

**Next decisive target.** Couple the existing seventeen-query two-port
carrier to actual neutral-native source resets, and certify that its chosen
children reenter a closed source family. Use P1/P2 to test all port cuts of
each supplied source efficiently; use the actual graph and joint response
to audit reentry. Try a common chart for the child *family*, rather than
mistaking one source's three exits for induction. A proposed MRB scalar
selector must first survive the triangle hostile and native return loops.

Reproduce from this directory:

```sh
python3 verify.py > receipt.json
python3 -O verify.py > /tmp/mrb-square-probe-optimized.json
cmp receipt.json /tmp/mrb-square-probe-optimized.json
python3 independent_audit.py > audit.json
```

The universe is all1,099 labelled sources of orders1..5, all33,866 one-birth
cuts, and all637 subsets of1..10 of sizes1..5 for the square-bank check.
The finite triple-anchor hostile uses1<=s<=12,1<=x<=10000; its universal
exclusion is inherited, not inferred from this scan. The all-size statements
P1–P3 are proved above. See `AUDIT.md` for scope and independent paths.
