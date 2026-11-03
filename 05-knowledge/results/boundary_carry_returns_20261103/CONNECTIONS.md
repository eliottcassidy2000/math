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
