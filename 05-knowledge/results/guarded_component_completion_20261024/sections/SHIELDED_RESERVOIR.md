# Shielded finite reservoirs behind arbitrarily large components

Status: **PROVED**, conditional on the explicitly supplied finite sections
and exact full-record compatibility. The finite40-source instance is
authenticated in the prior native supply packet.

## The boundary theorem

Let Q be any finite abelian charge group. Partition a finite actual native
inventory into core slots B, controls C and reservoir slots R. Retain every
slot's complete allowed records, normalized resource capacities, source and
target owners, actual clocks and all overlap predicates. Each reservoir
slot has two fixed actual records, one for U and one for V; a bit that stays
unchanged retains the same full record, including its kernel/lift marking.
This binary full-record hypothesis is essential for the reset bounds below.
Suppose:

1. At the core states being considered, there is a chosen feasible control
   configuration c on C with exactly kappa V slots. Its internal component
   size and predicate structure are unrestricted.
2. The reservoir has a complete fixed-rho charge section: for every q in Q
   at least one rho-subset of R has charge q, with actual legal records.
3. The reservoir is **shielded**: every complete reservoir configuration
   used by the section or prescribed as an endpoint is jointly compatible
   with every permitted core/control configuration. This is a full
   configuration predicate, including arbitrary higher-arity factors, not
   a test of isolated reservoir records. In
   particular, no untested owner equality, local partner arc, source or
   inverse-target resource conflict is hidden across that boundary.
4. The full physical target-margin requirement is the core occupancy plus
   kappa+rho, and its required Q-charge is specified. The relevant native
   paired-code theorem is applicable, with the actual degree terms retained.

Then every such feasible core/control configuration extends to every
required total Q-charge. Choose the reservoir section at

    required charge - charge(core) - charge(controls).

This is a direct proof: the section supplies exactly rho V slots and the
missing charge; shielding supplies all omitted cross-boundary predicates.
No independent copies of R have been made. A single actual bank is reset
between successive arrivals, and credit is not reused inside an arrival.

For a transition between two chosen configurations, let b and c count the
changed core/control slots. The reservoir reset changes at most
2 min(rho, |R|-rho) slots. Hence the support bound is

    b + c + 2 min(rho, |R|-rho).

For actual clock cost, add the actual core and control costs to the sum of
the largest2 min(rho, |R|-rho) eligible reservoir slot-change costs. Overlap
may reduce the realized cost but is not silently subtracted. The theorem
also permits prescribed legal actual reservoir endpoint configurations
when they have occupancy rho and the required charge and satisfy shielding;
use them at the corresponding endpoint instead of the default section.
A repeated core/control state with the same required total charge can be
assigned the same section, providing literal endpoint return.

If control occupancy varies with the boundary state, this statement does
**not** apply using one fixed-rho section. One needs a complete reservoir
section separately at every required occupancy rho=r-kappa. Charge
surjectivity does not propagate automatically between adjacent occupancies:
the four distinct labels of F2^2 cover the group at occupancies1 and3 but
miss zero at occupancy2. See the exact K4/no-resource-reuse hostile in
`../../compatible_charged_helpers_20261023/controller/PADDING_PROBE.md`.

## Actual40 plus balanced controls

The complete40/k3 bank from
`../../compatible_charged_helpers_20261023/supply/selected_supply.json`
is a subset of the actual53 bank. Its complement consists of exactly

    368,556,753,1369,1952,2066,2135,2397,2439,2483,2638,2948,3187.

Select an even number t<=12 of these actual outside slots as controls,
with exactly t/2 V records in every boundary state, leaving all unused
outside slots fixed U. Once shielding and all full-record overlap predicates
have been authenticated, the helper occupancy is3+t/2 and the core two-slot
support bound is8+t. This numerical statement does not assert that every
choice of control predicates is physically realizable.

The companion `../native/NATIVE_COMPONENTS.md` provides the actual
realization: ten controls share an owner pair with core176, while two
exceptional controls share another pair with core1532. Five same/five
opposite, and one same/one opposite, give exactly six control V slots at
all four root states. Its explicit source and partner embeddings certify
shielding of the40 private reservoir cells. This realizes two shared-owner
components, with occupancy9 and support at most20 (the ten-control version
has occupancy8 and support18). The companion receipt retains the precise
actual clock bounds; the proof here does not replace them by bit counts.

The construction's scale freedom is in the component predicate, not in
free resources: an arbitrarily large feasible control component can be
hidden behind such a boundary only after its actual slot inventory and
constant occupancy have been supplied. No all-size finite-supply claim is
made. If a small factor graph defeats the section, as the XOR tree in
`PAIR_RELATIONS.md` does, either strengthen the boundary summary or supply a
new authenticated helper. The relative-block source2125 provides an example
of the latter without increasing occupancy.

## Preserved and forgotten information

The map is from feasible full-record boundary states to (occupancy, Q-charge),
followed by a section into the shielded actual reservoir. It preserves the
specified target, aggregate electrical equations, indexed or relative pins,
kernel lifts, inverse/converse resources and actual clock sidecars. The
charge summary forgets internal component arrangement; universal shielding
is precisely the hypothesis that permits that forgetting. Without shielding
the conclusion is false, even if every separate two-slot constraint is
individually harmless. The independent pair census and native XOR-tree
obstruction give an exact finite failure boundary.
