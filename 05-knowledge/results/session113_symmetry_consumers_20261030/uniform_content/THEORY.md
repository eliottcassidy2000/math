# A finite obstruction to the universal phase-gated epoch envelope

**Status: exact, independently reproducible counterexample to the envelope target.**
This is NOT a counterexample to the sharp Collatz content inequality, nor to the
fixed-precision quotient certificates at E≤8. The current upper-envelope
consumer cannot support an all-layer potential, irrespective of how it is
constructed or lifted.

Source lineage: [the inherited budget](../../layer_lifts_cyclic_intersections_20261030/COLLATZ.md), [the E6–8 certificates and first obstruction](../collatz/THEORY.md), and [independent boundary audit](../audit_collatz/AUDIT.md).

## 1. Explicit positive valuation word

Use chronological concatenation and the original full affine coefficients.
Set

    R=(5,1^8,5), W=R 1^233 R, U=W 1^8 W.

Here `1^j` means j repeated letters1. Word U has length514 and binary clock546.
All its prefixes have square multiplier modulo5: letters1 and5 both preserve
this phase. At its eight non-one letters the defect heights are

    (1,4,1,8,1,4,1,11).

The intervening runs of ones have lengths

    (8,233,8,8,8,233,8).

The valuation κ(5)=1, so replacing the envelope's `max(h,κ)` by the tighter
`h` changes nothing here.

## 2. Exact excess and failed universal target

With λ∞=5/8, the untruncated envelope excess is exactly

    C∞ = (1+4+1+8+1+4+1+11) − (5/8)*40
         − (5/8)*5^(-1)*32 − (5/8)*5^(-4)*466
         − (5/8)*5^(-8)*8
       = 479371/312500 > 0.

This is precisely J−B for the inherited phase-gated, discounted-clock
budget. Thus the prior OPEN hypothesis J≤B for every positive word is false.
At E=10, only the last height is truncated, while the finite baseline is
smaller by (5/8)5^-10. Therefore

    C10 = C∞ − 1 + (5/8)5^-10 * 546
        = 1043012/1953125 > 0.

A putative E10 potential satisfying every entry inequality, every transition
inequality, and `max V≤0` would imply C10≤0 along this very word, contradiction.
This rules out all potentials for these weights, not merely the previously
computed lift normalization or a particular producer.

For comparison, exact full-coefficient arithmetic gives

    S∞(U)=S10(U)=317 < (5/8)*546=341.25.

The endpoint has defect height11 but actual content height1. The enclosure
loses return timing, despite preserving square phase and defect height.

## 3. Why the obstruction appears: resonant doubling

For any block X with normalized defect z_X and multiplier m_X, ones scale
z by3/2. Consequently

    z_(X 1^L X) = z_X * (1 + m_X (3/2)^L).

This identity retains the actual entry and does not reset source or clock.
For R, m_R=3^10/2^18 and h_R=4. The exact valuations

    v5(3^9+2^13)=3,
    v5(2^251+3^243)=4,
    v5(2^277+3^261)=3

give h_R=1+3=4, h_W=4+4=8, and h_U=8+3=11.
The low defect heights at proper non-one phases are unaffected when the
second copy is entered with strictly higher defect height; the ultrametric
inequality explains the listed eight heights. The checker also independently
iterates integer (P,Q,B) throughout the original word.

The two copies of W each have only a small negative budget at precision8;
waiting eight ones at height8 costs only 1/78125. Resonant doubling then
adds three units of endpoint height. The high-height clock discount makes
that wait extremely cheap, while the envelope charges the new height without
requiring multiplier identity at the same precision. This is the specific
mechanism that breaks the proposed universal induction.

## 4. Consumer and scope ledger

- Source: formal finite positive valuation words, with exact P,Q,B and clock.
- Quotient: square phase plus z modulo5^E; closed and valid.
- Preserved: truncated defect height and phase gate, hence an upper envelope.
- Lost: fine multiplier phase and the simultaneous identity-return predicate.
- Hostile probe: U above, positive envelope but content strictly under budget.
- Conclusion: sharp finite certificates remain valid; their current envelope
  cannot lift universally. The all-layer content bound remains OPEN.
- Repair direction, NOT a theorem: retain more multiplier phase or tighten
  finite-epoch counting, rather than trying a new potential for unchanged
  weights. The universal exact-return state lower bound still does not bar
  such one-sided, consumer-specific refinements.

## Reproduce

    python check_obstruction.py

Only standard-library integers and Fractions are used; `obstruction.json`
contains E8, E9, E10, E11, and untruncated endpoint/accounting receipts.
`witness.cpp` is an exploratory producer which extracts a shortest tight path
attaining the E8 maximum. It is unnecessary for verifying the counterexample.
No full E10 state enumeration is required.

## 5. A bounded low-phase repair of the interface (not universal feasibility)

Fix a small precision r≥1, independently of the target truncation E≥r. At an
active epoch start retain `m mod5^r`, in addition to defect height h≤E. Write
`T_e=2*5^(e−1)` and let d_e be the unique integer in `[0,T_e−1]` such that

    m (3/2)^d_e = 1 mod5^e.

For an epoch with L following ones, all counts at e≤min(r,h) are exact:

    n_e = max(0, 1+floor((L−d_e)/T_e)).

For e>r, the actual first-return offset d_e reduces to d_r modulo T_r, so
`d_e≥d_r`. In particular, there are no high returns if L<d_r. Otherwise
`n_e≤1+(L−d_r)/T_e`. Summation gives the rigorous upper bound

    F_r(m,h,L) = sum_(e≤min(r,h)) n_e
      + [h>r and L≥d_r] * ((h−r)+(λ_h−λ_r)(L−d_r)).

Inactive epochs contribute0. For h≤r the second term is absent. This is an
all-layer, finite-epoch **upper-bound formula**, not a proof that its sum is
always at most λ_E A. No difficult phase has been assumed to align.

At r=1, only m modulo5 is newly retained: the active values1 and4 determine
whether the first possible return occurs at d_1=0 or1. This is just one
additional bit on active states. It already repairs both known obstruction
words:

| word | refined upper bound F_1 summed | 5A/8 |
|---|---:|---:|
| U in this packet | 26448046/78125 =338.5349888 | 341.25 |
| sibling's E9 obstruction U3 | 198661/1250 =158.9288 | 169.375 |

At r=2 the final height11 of U contributes only its true level1 return,
because its multiplier is not1 modulo25. The resulting total upper bound
is5094453/15625=326.044992, still above the exact317.

**Bounded online interface.** One can implement F_r without retaining an
unbounded run length: retain `m mod5^r`, `z mod5^E`, and one pending flag.
At each epoch start reset the flag until its first low-phase return occurs.
Count low-layer returns exactly at each step; at that first possible high
return charge h−r, clear the flag, and on subsequent ones charge the tail
rate λ_h−λ_r. No tail-rate charge occurs while the flag is pending.
A new non-one letter ends the old epoch and resets this logic. Thus at most
`2*phi(5^r)*5^E` states suffice (an overcount allowing unreachable flags).
For r=1 this is8*5^E, still linear in5^E rather than quadratic. The old
phase-gated quotient had2*5^E states. This is a precise consumer-specific
repair candidate; existence of a universal potential for it is OPEN.

`check_phase_refinement.py` evaluates this bound using exact rational
arithmetic. It does not treat success on these two hostile probes as an
all-word theorem. The finite exact-return-observer lower bound is unchanged:
this still overcounts high returns and cannot answer all identity queries.
