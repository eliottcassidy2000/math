# Signed address chart-change sheaf for the three recursion modes

The repo already knew the formulas:

```text
A+B+C-D-E-F+G
A+B-C
A+B-C+D-E-F+G
```

The deep-dive correction is that the similarity is not "these sign strings
look alike."  The similarity is that all three are charted inclusion-exclusion
operators, and the chart determines what the letters mean.

The full `A+B+C-D-E-F+G` formula is the THM-442 full-staircase `B3` Mobius
face.  The even `A+B-C` formula is the THM-549/550 complement-folded `B2`
edge: pronic, fixed-line aware, and compatible with HYP-3232's pure
Eisenstein scale fold.  The odd formula is the subtle one.  In prompt order it
is `A+B-C+D-E-F+G`; in the corrected Venn geometry it is
`A+B+D-C-E-F+G`, with corners `A,D,B`.  The scalar cancellation between `C`
and `D` is not simplification.  It is the warning that a corner and an overlap
occupy the same size stratum and must be kept apart.

That is the recursive pattern that was still underdescribed:

```text
signed recurrence = local chart + destroyed-coordinate sidecar.
```

This sits naturally on top of the latest route.  HYP-3231 says the route is
scale-normal: normalize, expose what scale cannot kill, attach the sidecar,
and recurse.  HYP-3232 says the three modes concentrate at the apex:
Mobius/moment order, Eisenstein/modulus fold, and Legendre/speed ratio.  The
new HYP-3233 layer says the proof packet must also know which signed-address
chart it is in before it is allowed to move between those modes.

The other related concepts found in the repo now line up as chart changes:

- HYP-2899: divisor Mobius `phi=mu*id` is the exact-period denominator chart,
  so `Div(D)` has to stay beside the Boolean/character chart.
- HYP-2901/2902: the corrected odd Legendre Venn is a local address sheaf, not
  a global reuse of the same letters.
- HYP-2681: the cube-root `S_omega/P_omega` coordinates are phase sidecars for
  the seven packets; the pair-tax sign can reverse actual residual sign.
- HYP-2704: survival middle mass is a tail-depth quotient of the seven-packet
  chart, not the chart itself.
- THM-553/HYP-2690: tile addresses `(beta,tau)` make the half-tiling fold a
  layer clock rather than a raw count.
- HYP-3004: odd/even and sign are recursion modes only when their preserved
  predicate and destroyed coordinate are declared.
- HYP-3230/HYP-3232: the cap kernel is the speed-ratio chart; it is stable
  under modulus scaling until the apex break exposes the antipode-side debt.

The practical proof rule is:

```text
packet = scale orbit
       + signed chart
       + local slot basis
       + exact-period denominator packet
       + fixed-line/apex sidecar
       + moment-depth target
       + terminal discharge or named debt.
```

This prevents a recurring failure mode in the archive: scalarizing a sign word
and then trying to read cap risk, Hamiltonian-path behavior, or witness
denominator behavior from it.  The sign word is only the shadow.  The proof
object is the chart-change packet.

Next useful audit: take a small HYP-2963/HYP-3083 packet sample and add
`signed_address_chart`, `local_slot_basis`, `cancelled_same_size_slots`, and
`apex_break_defect`.  The first theorem target is a chart-change legality
lemma: every move full -> even -> odd -> survival/cap either preserves the
LRC predicate, reconstructs the lost coordinate, or emits named debt.

-> HYP-3233, HYP-3232, HYP-3231, HYP-3230, HYP-3216, HYP-3004, HYP-2902,
HYP-2901, HYP-2899, HYP-2704, HYP-2685, HYP-2681, THM-553, THM-550, THM-549,
THM-442, LTI-331, LTT-231, T1331, OPEN-Q-108.
