# Audit, corrections, and handoff

2026-11-03. Same-agent proof review with independent computational paths.
No subagent or external referee was used. The inherited seed/cover/compiler
is credited to S167/S168; no new literature-priority claim is made.

## Exact universes and independent paths

- All64 one-vertex cuts on each of the six specified seeds are examined by
  subset Hamilton-path DP. Exactly16 cuts separate all six. The selected
  cut2 is independently checked by all7! literal vertex orders, on all six
  extensions. The rooted formula is checked by subset-cover convolution.
- All logical states for widths1 and3 are classified by depth and exact
  orbit period. All36 and46,656 physical states, respectively, are exhausted
  without prefiltering. The independent matrix audit checks all fixed-space
  dimensions at widths1,3,9,27, rather than inferring them from a sample.
- Impulse logical/physical clocks are checked at widths9,27,81. Uniform
  source classification comes from CB1's proof, not those samples.
- Logical clock orders are computed by polynomial primary-factor orders at
  widths1 through12 and16. A separate dense-matrix implementation verifies
  return and every prime-divisor minimality test. Both verify cyclicity;
  one-site observability is also checked for every width1 through20 and27.
- Nonlinear physical-boundary histories exhaust every logical source and
  both observed deck bits at widths1,2,3,4. Remote deck bits are quotiented
  only after their invariance is proved. Counts come independently from
  literal gate simulation and dense logical matrices plus the carry formula.
- Both scripts use explicit exceptions and remain active with `python3 -O`.
  The normal and optimized generated JSON must be byte-identical.

## Proof audit and hostile controls

1. **Actual update recovered.** The initially available S167 I+S rule has
   order2m, not8m. During this session origin/main supplied S168's exact
   two-track R rule and literal hexagon. The exploratory alternative
   `(x,y)->(Sx+y,x)` was not canonized or substituted for R.
2. **All source periods.** CB1 first forces a multiple of8 using the
   nonzero leading coefficient, then proves exact nilpotent order. It does
   not assume all nonzero sources have maximum period. A constant source
   is the period8 hostile. Zero has4^m distinct physical fixed lifts.
3. **Physical return.** The cancellation uses a state-specific half-turn
   and the specified linear gate schedule, not merely an even period.
   Widths2,4,8 double the impulse period. Arbitrary-width logical orders
   do not imply no doubling.
4. **Readout type.** CB2 retains a physical label and adds one new vertex.
   It does not recover a deck bit from an unmarked isomorphism class.
   The probe changes the profile; querying a copy preserves the source.
   The simple product multiplier is asserted for the ordinal-sum carrier,
   not automatically for the separate strong wrapper.
5. **Sharp observation count.** Linear rank alone would not give a lower
   bound for six-valued measurements. CB3 supplies the actual carry formula
   and exact shorter-history classes. It also proves which remote bits
   remain invisible at every time. Changing the observation alphabet can
   change this equivalence relation.
6. **Every width.** CB4 treats the inseparable spatial multiplicities and
   the critical quadratic separately. The repeated factor occurs at
   widths divisible by8 after removing powers of3. No diagonalization of
   a nonsemisimple ring is silently used.
7. **Refinement and time.** CB5 uses explicit compatible projections, not
   independently selected cycles. CB6 keeps gate attempts, nontrivial flips,
   parallel layers and rounds distinct. A finite profile fibre supplies no
   natural-number descent certificate for this cyclic operation alphabet.

## Corrected draft implication

During checking, the draft expression H+K_q for the sole-in-neighbor birth
failed on the very first seed: it predicts112 instead of109. The omitted
quantity is s_q, the old Hamilton paths starting at q, which cannot follow
the new vertex. The exact formula is H-s_q+K_q. A singleton is the minimal
hostile: H=K_q=s_q=1, while the extended chain has H=1 rather than2.
The correction is in CB2 and the shared mistakes ledger. No previously
published theorem is retracted. An audit-output variable-name collision was
also repaired before filing: the birth-count field now holds the six actual
path counts, not the last history-census list.

## Handoff

**Completed:** all-state depth/lift theorem; explicit one-vertex six-state
reader; sharp logical and physical boundary-history length; exact invisible
deck subgroup; arbitrary-width logical period prescription; formal-limit
nonreturn; integer schedule accounting; older-work connection map.

**Do not reprove:** the S168 8m impulse theorem, its literal hexagon lift,
the S167 full-profile product, or the older golden F9 order8. The separate
joint-port packet already proves finite31- and27-dimensional observers
for its specified attachment languages. Generalize its operation alphabet
before treating observer closure as an open question.

**Next decisive tests:** (i) two fixed spatial probes: optimize the exact
joint observation rank and carry refinement; (ii) a formula for the impulse
deck parity at general widths, tested first at2,4,5,7,8,11; (iii) apply
boundary-response comparisons to the actual spectrum return programs,
retaining their source, ports and costs. A clock theorem by itself does
not close their universal marked-reentry gap.

The small board and four helper questions are in [CONNECTIONS](CONNECTIONS.md).
The reproduction commands and source-pinning scope are in [README](README.md).
