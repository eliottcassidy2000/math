# A small board of clocks, boundaries and niche connections

2026-11-03. The board is **clock valuation / fixed probe / deck bit /
integer schedule / compatible limit / native return**. This is a targeted
inheritance and synthesis pass, not an assertion that all historical uses
of these topics have been reread. Proof references below use both path and
local theorem label; none relies on a bare, potentially colliding THM ID.

## Inheritance and portfolio

- **Anchor:** [S168 cellular](../session168_computational_structures_20261103/cellular/NOTE.md)
  gives the precise 8*3^a rule and its real tournament gates. Extend it to
  every source and identify a useful boundary observer: CB1–CB3 in
  [PROOFS](PROOFS.md).
- **Niche:** the older golden-field reader and ramified-prime clocks explain
  why a period needs both a residue-field order and a nilpotent depth. CB4
  turns this into an exact prescription at every spatial modulus.
- **Wildcard:** the inverse limit has no nonzero periodic state although
  every finite truncation is periodic (CB5). This is a hostile test for
  finite-modulus proof strategies, with a precise compatible source.
- **Canonical hostile:** widths2 and4 double the physical period; a ternary
  three-cycle cannot close after three labelled bit flips. Also, remote
  deck-bit changes are permanently invisible to the fixed boundary probe.
- **Corrected near miss:** adding a vertex with sole in-neighbor q does not
  give H+K_q. Paths beginning at the new vertex cannot next visit q; the
  exact count is H-s_q+K_q. The correction is tested and logged.
- **Least-used sidecar:** the operation-specific observation alphabet. A
  port fixed in physical labels is different from one transported under
  the deck involution, even when their unmarked base graphs are identical.

## 1. Golden arithmetic: the same residue clock, a different thickening

The [October 3 golden prime clock](../golden_prime_clocks_20261003.md),
sections1–2, retains `(a,b)` for a+b*phi and proves that phi has order8
modulo3, order20 modulo5, and the nilpotent obstruction to arbitrary
division. For the current constant spatial mode, M(1) has polynomial
T^2+T-1, so it is multiplication by **-phi**, not phi, in F9. This is an
explicit companion-field identification; both elements have order8.

Here u=S-1 supplies spatial nilpotence. At m=3^a its nilpotence index is m,
and the surviving source depth d fixes the additional power of3. In the
older golden modulus5 example the nilpotence arises from ramification of
the quadratic number ring. The preserved predicate is an element's exact
order; the origin of the nilpotent coordinate differs. Keeping only the
residue field destroys the extra period in both cases.

CB4 adds a second spatial mechanism: the rational map
`r(T)=T-1-T^(-1)` has critical points T^2=-1. Their images have order8, so
widths divisible by8 acquire one extra factor3 in the logical period.
This is a ramification calculation, not another deck-bit effect.

**Cheap test:** widths4 and8 have prime-to-three order40; the latter has
logical order120. Their physical impulse periods are80 and240. Field
ramification and labelled-cover doubling multiply separately.

## 2. Joint tournament clocks: scalar recurrence is not boundary recurrence

The incoming [joint-clock packet](../sheet_return_compiler_20261103/JOINT_CLOCK.md)
has two cubic sectors, a constant endpoint imbalance, and a seven-coordinate
integer realization for a fixed rooted observer. Its F49 and F343 growth
components give period2736 modulo7; its parity inverse guard remains
necessary over integers. The associated [sheet construction](../sheet_return_compiler_20261103/SHEETS_AND_SINGER.md)
retains the scale lost in a projective quotient.

Our clock has a particularly explicit cyclic observable: one site and its
first2m iterates span the whole dual. The map is v to its boundary history.
It preserves every logical source, loses all but the observed cell's deck
bit, and needs named cells plus the schedule. The analogue in the native
return work is to ask which rooted observables have a closed, finite
future-observation space under the **actual permitted operations**.

This suggests a test, not a solution: take a proposed scalar/phase return
state, attach the permitted rooted probes, and compare their next responses
on same-scalar sources. A mismatch identifies the missing coordinate before
trying to prove a universal budget. Our one-vertex switch reader is a
worked finite instance; it does not supply reentry for every ternary H output.

## 3. Cayley–Dickson transport: fibre size and lap count are different

The newly merged [S168 signed transport clock](../session168_computational_structures_20261103/clocks/NOTE.md),
Theorems2–4, proves that a fixed-level signed-unit transport has form
`(s,x)->(s+f(x),x+a)`, hence exponent at most4. A 26*3^a base clock can have
one-, two-, or four-lap lifts regardless of how large the unit register is.

The tournament switch has2^(2m) deck sheets but exponent2, so a logical
orbit can at most double. CB1's half-turn cancels the accumulated deck
translation entirely at every source in the 3-power-width family. The map
in both examples is projection of a retained register; the preserved
predicate is base evolution, while the lost quantity is the one-lap action
on the initial sheet. Sheet count does not determine that action's order.

**Boundary:** a six-state machine computing a basis multiplication sign is
not a summary of arbitrary accumulated all-sheet transports. Likewise the
one-vertex probe reads one switch, not an entire clock in a single query.

## 4. Older LCA and Zeckendorf work: a boundary must answer future queries

The [October 25 RMQ/LCA atlas](../session56_rmq_lca_connections_20261025/CONNECTIONS.md),
cards1,6,12–15, distinguishes query structure from value, p-adic prefix
depth from integer magnitude, and Fibonacci colour from its carry cocycle.
It also records how the same ancestry can have linear and exponential
integer addresses.

The current u-adic valuation is a common-prefix depth of polynomial
coefficients, just as v_p(x-y) is a least-significant-digit LCA depth. The
exact map is coefficient expansion at z=1, not ordinary base3 digits of
the Hamilton count. It preserves the period stratum and discards the
phase and residue-field direction. CB2's boundary readout restores a
different part of the state; it cannot infer phase from valuation alone.

This gives a useful layered representation: **depth chooses the possible
clock length; source and phase choose the trajectory; a port chooses which
trajectory data can be queried.** A Fibonacci carry label or a switch deck
bit matters only when an admitted subsequent operation can expose it.
CB3 gives both directions of that statement for one explicit clock/probe
alphabet, rather than only an information-loss example.

## 5. July odometer work: retained clocks can coexist with an autonomous obstruction

[THM-2693, odometer-skew-product-three-event-escape-and-uniform-delayed-depth-four-nilpotence](../../../01-canon/theorems/THM-2693-odometer-skew-product-three-event-escape-and-uniform-delayed-depth-four-nilpotence.md)
proves a surviving three-event interval, yet an autonomous terminal
coordinate obeys y'={13y} for every integer lift and makes the stated
delayed carrier empty at four events. The clock/carry repair does not steer
that coordinate. Its scope does not settle LRC(14).

The current source-to-boundary map has an equally exact steering boundary:
remote deck bits never affect the observed cell under any number of the
specified rounds. No improvement in time resolution overcomes it. A new
allowed port query does. This is a structural analogy between autonomous
factors, not an identification of the LRC dynamics with this cellular rule.

**Next decisive test in a new problem:** compute which source coordinates
actually influence the claimed terminal predicate before enlarging the
clock register. A longer clock alone may leave the obstruction untouched.

## 6. Fixed twists versus fresh twists, and the compatible-limit hostile

The [October 25 twist mass theorem](../session63_integral_observers_root_arrangements_20261025/twists/RESULT.md),
Theorem B, proves that a fixed twist g changes over degree-r extensions by
the power map g->g^r. Freshly retwisting at each field does not describe one
compatible family. Its square-extension point counts explicitly separate
these choices.

Similarly, choosing a convenient cycle independently at each modulus is
not a compatible source. CB5 goes further: even the **same compatible
impulse** can return at every finite level with unbounded return times,
while never returning in the inverse limit. The bonding maps and an
integer-time bound are indispensable. This is relevant as a hostile for
Collatz modular controllers, not an embedding of a nonconvergent integer
Collatz orbit.

## Next questions, with cheap decisive tests

1. **Expanded native boundary closure (OPEN).** The
   [joint-port packet](../joint_port_geometry_20261103/README.md) already
   proves minimal31- and27-dimensional observers for two specified
   attachment languages. Retain those results. The new question is closure
   after adding actual variable-order path complements and return repairs,
   beyond those languages. Start with the actual195/197/199 sources and
   their retained port6 from
   [NATIVE_BRANCHES](../six_roots_hidden_bits_20261103/NATIVE_BRANCHES.md).
   Same H is insufficient; compare the first permissible probe responses.
2. **All-width physical clock (OPEN here).** CB4 computes the logical order
   at every width. Find an arithmetic characterization of the one-lap deck
   parity for the impulse. Widths2,4,8 double;5,7 do not. The mechanism
   cannot simply be the uniform half-turn: the 5-site factors have root
   orders8 and164, with different powers of two.
3. **Boundary history compression (resolved here, then extended).** The
   initial helper question asked whether the six-valued reader could beat
   2m queries. The carry identity CB3(9) proves it cannot in the worst case:
   the last observation only partially separates the next trit. The exact
   shorter-history count is10*3^(L-1). The next OPEN question is to optimize
   two fixed probe locations, retaining their unknown initial deck bits:
   compute the joint observation rank and then the carry refinements.
4. **Refinement as a legal move (OPEN).** Logical folding to a shorter ring
   is explicit, but a port-preserving physical compiler for it must account
   for the destroyed source data, inverse guards, and costs. Try two sites
   sharing one folded target; test injectivity before calling it reversible.

The current packet does not fill the spectrum's universal marked-reentry
gap. It supplies a sharper test of a proposed interface and a complete
small example in which the arithmetic clock, hidden bits and probe are
all explicit.
