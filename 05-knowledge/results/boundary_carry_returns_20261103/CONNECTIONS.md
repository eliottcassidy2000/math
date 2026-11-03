# Connection ledger: the arithmetic clock and its retained boundary

The connections below are maps or tests with explicit limitations. Equal
integers in unrelated classifications are not used as proof dependencies.

## 31,33 and 12,14: the same projective boundary mechanism

For a finite field F_q, the projective line has q+1 points. Multiplication
by a primitive field element fixes 0 and infinity and cycles through the
other q-1 points. Multiplication by a suitable element of the quadratic
extension, acting on one-dimensional F_q subspaces, instead cycles through
all q+1 points. The quotient is E*/F_q*, and its missing coordinate is the
scalar multiplier in F_q*.

Here are explicit, exactly checked instances; the field constructions and
power cycles in `arithmetic.py` suffice without a literature assumption:

* F_32=F_2[z]/(z^5+z^2+1). The element z has order 31; its 31 powers give
  every nonzero element and certify that this quotient is a field. Adjoin
  u with u^2+u+1=0; the polynomial has no root in F_32. The element
  eta=z^2+z+u has order 1023=31*33. Its first scalar power is eta^33,
  so its projective order is 33. The split clock is 31 and the nonsplit
  clock is 33; the lift retains a 31-valued scalar phase.
* F_169=F_13[phi], phi^2=phi+1. Its discriminant is 5, nonsquare modulo 13.
  The scalar 2 has order 12. The element eta=2+phi has norm 5 and order 56,
  eta^14=5, so its projective order is 14. Its scalar fibre has order 4.

Thus **32+-1 and 13+-1 have the same two-fixed-points versus no-fixed-points
mechanism**. This is stronger than a list of recurring numbers. The count
difference two comes from the split action's two distinguished endpoints.

The choice of generator matters. The tempting z+u has projective order
11, not 33: its linear order is 341. And phi itself has order 28 but
phi^7=8, so its projective order is 7, not 14. It acts as two 7-cycles on
the 14 projective points. Adding 2 changes the norm from -1 to 5 and reaches
the full 14-cycle. These are hostile controls against extrapolating an
ambient group order to a particular generator. The rejected z+u guess is
recorded as a draft error, not a retracted published theorem.

**Source -> target:** multiplication on E* -> its action on E*/F_q*.
**Preserved:** projective position and exact quotient period.
**Lost:** the accumulated scalar; restore a 31-phase or 4-phase sidecar in
these examples. **Cheapest test:** compute the first scalar power, not only
the first power equal to 1. In particular a universal one-hidden-bit model
fails even for the golden example, which needs four scalar states.

This does not implement a tournament operation. Even an explicit enumeration
of a tournament's 31 or 33 Hamilton paths would only provide a set bijection;
it would not prove that projective multiplication is a native move, respects
the selected pins, or preserves the next birth response. BC2 is a legal
compiler for supplied named targets, but cannot supply a target from a field
address alone. The arithmetic connection motivates the boundary bookkeeping,
not a shortcut around the outstanding realization theorem.

## Recovered old and concurrent mechanisms

### An arithmetic address for the carry

The inherited [six-root return forest](../triple_reentry_interfaces_20261103/README.md)
has a particularly useful normalization. Write L=(H+1)/2 for an odd raw
count H. The three desired count operations become

    3H       -> 3L-1,
    3H+2     -> 3L,
    3H+4     -> 3L+1.

Thus the arithmetic address is balanced ternary. The forbidden raw counts
7 and 21 have addresses 4 and 11. Removing those vertices leaves roots
at addresses 1,12,13,32,33,34, corresponding to raw counts
1,23,25,63,65,67. Here 12 and 13 are addresses, not even Hamilton counts;
address 14 corresponds to raw H=27 and has parent 5, whereas addresses
12,13 have the removed parent 4. A sliding triple such as 13+-1 therefore
crosses a genuine parent boundary. Its members need not have one source.

For the two-offset sublanguage at our raw H=31 seed, use instead D=H+1.
Each branch bit b=0 or 1 acts as D->3D+2b, starting at D=32. After ell steps,

    D_ell / 3^ell = 32 + sum_(i=1)^ell 2b_i/3^i.

This is an exact integer clock with a finite ternary fractional address.
The limiting address set is the translated ternary Cantor set in [32,33].
The four newly realized length-two addresses are D=288,290,294,296.
Beyond length two the formula specifies requested counts, not proved
realizations from the current source. The inverse digit word also does not
recover a tournament's boundary relation: several named sources can have
the same H and different legal continuations. That is precisely the missing
sidecar our matching and path-cover records retain.

Even complete realization of this binary sublanguage from one fixed seed
would give only 2^ell counts over a range of size proportional to 3^ell.
It would not prove the full odd spectrum. The reusable multiplicative
branch and all required source seeds remain essential. The 31/33
refinement pair in BC4 is a different, literal graph comparison; its
changing gap must not be replaced by this address interval.

| Prior source | Map / preserved predicate | Information that must survive | Effect on this session |
|---|---|---|---|
| [THM-1865, score and inflation](../../../01-canon/theorems/THM-1865-h-not-score-determined-and-inflation-neutrality.md), July | Same score multiset -> different actual path counts; source/sink insertion gives a path bijection | Named arcs and endpoint role | Produced the literal 31/33 pair; BC4 shows its gap changes sign under the same refinement. No old random maximality claim is inherited. |
| [THM-3121, run substitution](../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md) | Actual module paths -> ordered skeletons | Run counts, factorials, attachment endpoints | BC1 identifies their exact top degree by a matching; BC3 certifies infinite identities from coefficients. |
| [H31 source switching](../source_selection_closure_20261103/BRANCH.md) | A neutral six-state component changes an available ternary branch | Actual pins (5,1), named path, return permutation | Recovered the smaller source instead of extending only the 1457 bank. Neutrality alone did not imply the same clone response. |
| [Six concurrent clone returns](../adaptive_four_returns_20261103/LIFTING.md) | Common predecessor/successor permits two-complement lifting | Full intermediate path, clone address, separate refinement clock | Extended the one-neighbour sufficient condition to a source with three run modes; retained rather than re-proved its common-frame mechanism. |
| [H701 ready obstruction](../return701_algebra_20261103/RETURN701.md) | Complete inverse source grammars test ready defects | Class existence is different from existence in a ready interface | Chose a broader, explicitly charged helper language; no claim that re-marking alone repairs readiness. |
| [Marked native compiler](../adaptive_metagraph_reformulation_20261103/METAGRAPH.md), §5 | One-helper two-complement word reverses one arc | Named cut, path guards, separate birth/delete clocks | BC2 adds a common-frame guarantee for every clone size, at a larger explicit helper budget. |
| [CM trace carrier](../cursor_cm_phase_20261028/TRACE_CARRIER.md) | Sampling a quadratic unit -> a trace quotient | Smith kernel and actual sampled exponent | Suggested testing first scalar return instead of only total order. The older a=5, discriminant-29 example modulo 13 is distinct from the present golden discriminant-5 field. |
| [Five ready-tail families](../ready_return_pentagon_20261103/PENTAGON.md) | Actual tail operations form a five-vertex pentagon; prefix doubling commutes | Current tail and physical root | Its 13*2^k+3 and 14*2^k+3 neighbours differ by 2^k, not a constant. This is the same kind of changing-gap warning as BC4, with a different precise map. |
| [Square-sum nice-pair recursion](../ninja_pairs_square_sum_number_system_20261031.md) | Paired endpoint interface -> a recursively admissible pair | Both endpoints, parity and block order | The useful analogy is proving closure of the boundary class first. Scalar square-sum Hamiltonicity or a tournament count does not itself give that closure. |
| [Mersenne clocks and F21](../mersenne_switch_parity_f21_compression_20261006.md) | Literal source arithmetic and finite-depth CRT periods | 2-adic and 3-adic digits, marked source | Keeps the 2^m and 3^m modes separate. Neither finite periodicity nor our return tree proves Collatz coverage. |
| [Scheduled boundary carry](../clock_boundary_simulation_20261103/PROOFS.md) | One-site history -> logical state plus an observed lift bit | The other deck bits, gate schedule and integer clock | A faithful computer can be hidden in a constant-profile bank, but ordinal attachment multiplies H by powers of 29. It is not a free prime-count controller. |

Every link was read for its stated mechanism; historical numerical summaries
were not promoted over corrected current results. The numerical search was
targeted to these actual mechanisms, not a claim that every occurrence of
31,33,12,14 or 5 in the repository has been exhausted.

## The reusable move and a sharper open question

### Incoming work integrated after the first checkpoint

The clean-tree pull added eight commits through `2d137c00ea`. The most
relevant sources were read at their actual proof level:

* [S170 guard resets](../session170_holonomy_homotopy_20261103/guards/NOTE.md)
  exhibits the two literal words DCD and D C^14 D. Their Hamilton endpoint
  transfers have rank one, while the full path-slot polynomial transfer is
  injective. Its value plus first derivative recovers the input value for
  these two reset words. The long reset has outer vectors (47,112,0) and
  (47,112,47), and repeating it multiplies by 14753. These exact matrix
  identities are distinct from the finite-field 14-cycle above. The source's
  all-gap exclusion uses an external Tribonacci zero theorem; this packet
  neither re-proves nor needs that exclusion.
* [S169's collision ledger](../session169_runner_jacobian_20261103/collision/META_LEDGER.md)
  derives an exact K+1 finite-horizon budget for its nonuniform head from
  the positive new exponential mode at each step. Comparing this with BC1
  led to **BC7's new d+K bound**: a matching cap pays for root returns and
  each permitted path component pays for another run. The two theorems
  concern different operators but share a proved triangular-rank mechanism.
* [S170 flat covers](../session170_holonomy_homotopy_20261103/monodromy/NOTE.md),
  §§3--5, gives an actual labelled tournament torus whose deck return is
  winding parity. Mixed-coordinate squares preserve the cover, whereas
  filling a coordinate triangle would erase its sheet bit. This sharpens
  the boundary-carry interpretation without changing the price H=29^N of
  attaching that clock bank. Its literal edges are arc flips; they are not
  silently counted as our native Hamilton-path complements.
* The [701 update](../return701_algebra_20261103/ALGEBRA95.md) supplies an
  exact degree-24 annihilator for a positive weighted sum of its persistent
  generators. It controls those matrix responses and their radical layers,
  not the existence of a positive word giving our desired affine offsets.
  Its incoming-note warning about an unpaired complement is retained: such
  a move changes a transitive block and cannot use common-frame cancellation.

The revised board therefore separates three costs: matching-limited clone
moments, full path-slot memory at a reset, and actual named reachability.
BC2 handles the third only for a supplied target; BC7 handles the first
only in a specified linear query space. Neither is substituted for the
remaining positive source-return theorem.

A final clean-tree pull added S171 through `942626a6e4`. Its
[one-vertex read head](../session171_knots_information_20261103/information/NOTE.md)
stores nonisomorphic five-vertex alternatives with the same full profile
and H=9, yet the same specified attachment reads H=17 versus 19. Its
[update barrier](../session171_knots_information_20261103/information/UPDATE_BARRIER.md)
shows that a constant-profile update can leave the promised error-correcting
code family between valid endpoints. This makes our boundary precise:
**BC7's moment-space dimension is not a sufficient controller state**, even
if many unpointed cover moments are supplied. Attachment context and the
completion clock still matter. Our native repairs already retain them and
do not claim that their intermediate tournaments remain two-cover sources.

The companion [winding-code theorem](../session171_knots_information_20261103/winding_codes/NOTE.md)
prices a hidden physical return at exactly three times code distance in
its arc-flip carrier. That is a different clock from a native complement;
BC2's explicit binom(L-1,2) conversion retains this distinction. These
incoming results strengthen the stopping conditions of the present proof,
without replacing its source-specific realization steps.

The next close-out pull through `37e6cd63a4` brought the most direct
convergence: [six-output COMMON_GRAMMAR](../six_output_return_20261103/COMMON_GRAMMAR.md)
independently proves the matching boundary, credits this packet's earlier
checkpoint, and identifies a cycle guardian at matching rank one. Its
arbitrary-child naturality applies to our rank-two programs as BC8. The
exact correction cH(S) isolates a selected-path realization problem; our
C3 and H=5 child controls check it on all six supplied steps. Its
[RETURN_BARRIER](../six_output_return_20261103/RETURN_BARRIER.md) also proves
that no finite fixed dyadic-ray bank closes under ternary offsets, and
that exact integer carry recovers word height and arithmetic LCA. Our
new d+K dimension theorem is not offered as a way around that obstruction:
linear query dimension and the number of realizable source states differ.

[The simultaneous 2107 return](../digit_cycle_ternary_20261103/RETURN2107.md)
advances that separate marked lane to an actual 6325 source. It does not
supersede our all-m family result or supply its missing third-level rules.
The next work should compare these actual source relations and their
cycle covers, rather than repeat the already completed 701/2107 searches.

The common operation across these lanes is:

1. identify a scalar clock or a quotient period;
2. recover the smallest boundary data used by the next actual operation;
3. prove the realization/transport map on that enlarged state;
4. separate target existence from the cost of reaching a supplied target.

In this packet the boundary has a combinatorial meaning: a cover of the
O->I graph, not merely a metaphorical carry bit. The new question is whether
the two-cover class admits a compact **positive** continuation grammar after
retaining its pinned ordered-cover responses. Formal linear combinations of
responses are insufficient if they do not correspond to legal programs.

The direct hostile is now available at every size: one-birth continuation
from the 95 family cannot reach either correct second-step function at the
same base order. A successful grammar must charge the extra permanent
vertices or change its representation. The complete second-level repair
shows that this is an interface/size obstruction, not a counterexample to
the scalar spectrum conjecture.

**A precise next hypothesis, OPEN and deliberately restricted.** Starting
from the marked A_m family, can representatives be chosen for every finite
word in {+2,+4} so that each step adds at most three permanent vertices,
changes at most four old exterior arcs, leaves the old clone star unchanged,
and has the requested 3H+c identity for every m? All words through length
two have the stated witnesses. BC2 would compile such an endpoint rule into
a uniformly bounded number of legal native events per step (with its stated
quadratic arc-clock cost). The first new test is both continuations of the
actual 295 family, not an arbitrary graph with H=295. This hypothesis is
stronger than scalar existence for this family and weaker than full spectrum
coverage; neither direction of that distinction is suppressed.

For the four-old-arc bound on the 97->293 branch, choose program 1 of
`next_0.json`. BC3's displayed program 0 uses six old arc changes but only
two native events; the two clocks need not rank programs the same way.
BC6's boundary table refers to the displayed program 0, not this alternative
target. Both named targets and all their path guards remain in the packet.
