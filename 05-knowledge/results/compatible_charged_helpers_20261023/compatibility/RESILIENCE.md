# Native charge banks under moving helper restrictions

Status: FINITE-EXACT native subset coverage; PROVED conditional compatibility
corollaries for the supplied banks. Full records remain in
`../supply/selected_supply.json`; no additional helper path is invented.

## 1. What is counted

The inherited helper records have targets U=1965,V=2234, equal chosen
per-slot original phase, fixed literal common owner clocks, complete pinned
U/V domains and globally disjoint normalized option resources. Their paired
charge group is Q=Z/128 plus (Z/2)^4. The inherited supply checker authenticates
these claims; `resilience.py` pins its complete source sidecar by LF-normalized
SHA256 and checks the resource/phase fields again.

For a bank B with exactly r helpers at V, define

    W_q = {S subset B : |S|=r and sum_(s in S) label(s)=q}.

Every other helper still occupies its U record. Forcing helper h to U means
requiring h not in S; it does NOT delete a slot or leave a resource vacant.
Forcing h to V means h in S. Charge values are encoded as
16*(X sum modulo128) plus the four-bit parity xor.

These tests concern fixed-choice restrictions on ACTUAL supplied records.
An arbitrary geometric pin, shared-owner equation or resource constraint is
not automatically equivalent to such a restriction. The application must
first compile those constraints, as in `NEUTRAL_FIBRE_LIFTING.md`.

## 2. Exact one-bit gates and minimal blocking sets

Forcing h to U destroys charge q iff h belongs to the intersection of every
member of W_q. Forcing h to V destroys q iff h belongs to none of them.
More generally a forced-U set E destroys q iff E meets EVERY member of W_q;
this is its hypergraph hitting condition. Witness count alone does not decide
resilience: many witnesses can share a single mandatory helper.

The checker enumerates every fixed-cardinality subset in each bank and obtains:

| Helpers | V occupancy | Minimum witnesses per charge | Helpers individually safe when forced U | Helpers individually safe when forced V |
|---:|---:|---:|---:|---:|
|17|8|4|9|1|
|18|6|1|1|0|
|22|5|3|19|0|
|32|4|5|32|0|
|53|3|2|53|0|
|40|3|1|0|0|

The unique forced-V-safe helper in the seventeen-helper bank is source1493;
forcing that SAME helper to U is also safe. This does not make every helper
bit safe in that bank.

The universal forced-U tolerance is EXACTLY one for the original32/k4 and
53/k3 banks. Explicit two-erasure adversaries are:

- 53/k3: charge305 is lost if sources776 and3139 are both forced U. Its complete
  witness list is{776,1182,2638},{1391,1948,3139},{1400,3111,3139},
  {2439,2728,3139}.
- 32/k4: charge247 is lost if sources862 and2616 are both forced U. All six
  witnesses contain at least one of those sources.
- 40/k3 is fragile under EVERY single forced-U helper: the surviving charge
  count ranges from2034 to2046. For example charge10 has just the witnesses
  {776,1182,2218} and{1182,2728,2752}, so forcing1182 to U destroys it.

For each charge, the checker determines the minimum hitting size exactly
when it is1,2 or3. It otherwise reports only a lower bound4 and an explicit
greedy upper-bound witness. Completeness through3 follows by branching on an
element of the first witness, then an element of the first still-unhit
witness, then testing the intersection of the remaining witnesses. Every
hitting set of size at most3 must appear in that branching. The routine is
independently compared with literal subset search on all2,223 nonempty uniform
hypergraphs on at most five vertices. No larger hitting minimum is asserted.

## 3. Frozen robust subbanks

Deleting redundant helpers while retaining every single-U gate gives these
two authenticated subbanks. They are inclusion-minimal for that property under
deletion, not proved globally smallest.

The48-helper bank, with V occupancy3, is exactly:

    368,556,650,676,753,760,776,789,811,827,862,921,
    1182,1248,1254,1369,1391,1395,1400,1493,1574,1596,
    1615,1624,1650,1702,1770,1882,1948,1952,2066,2135,
    2218,2283,2319,2384,2393,2439,2483,2616,2638,2677,
    2728,2752,2948,2983,3139,3187.

The28-helper bank, with V occupancy4, is exactly:

    368,556,650,676,776,827,862,921,1248,1254,1369,1391,
    1395,1493,1596,1615,1702,1770,1882,1952,2135,2218,
    2283,2319,2397,2616,2677,2728.

The48 bank has all2,048 charges after ANY one of its48 helpers is forced U,
and with no forced choice. The28 bank has the analogous29 gates. Every gate
is replayed through a complete2,048-entry section and saved by its section
hash. Complete reconstructible certificates are retained, rather than trusting
the deletion search.

For compact storage, each q keeps at most r+1 witnesses whose intersection is
empty. Such a subfamily always exists: begin with one r-set and repeatedly
choose a witness deleting at least one element of the current intersection.
At most r additions are needed. Any single forced-U helper is then avoided
by one saved witness. The verifier checks each chosen mask's exact cardinality,
actual source labels, charge, bank membership and forced-U exclusion for every
gate. Masks index the ORIGINAL corresponding bank's `source_ids`, not a
silently renumbered source set.

Both original two-erasure adversaries remain inside these reduced banks.
Thus their universal single-U tolerance is still exactly one.

## 4. Increasing occupancy gives arbitrary one-bit resilience

Use the SAME53-helper pool, or the SAME frozen48-helper pool, but place exactly
four helpers at V. Direct enumeration gives every charge under EVERY single
forced-U restriction and EVERY single forced-V restriction.

There is also a useful proof of the forced-V half: fix h at V and request
charge q. Choose a triple avoiding h with charge q-label(h), using the
single-erasure triple theorem. Adding h gives the required quadruple. The
53-pool forced-U half also follows from its32/k4 subbank, which survives any
single erasure inside that subbank and is unaffected by erasures outside it.
The48-pool forced-U half is separately checked by complete enumeration.

For both pools, the bounded hitting search finds NO forced-U blocker of size
at most3 for any charge. Thus up to THREE helpers may simultaneously be
forced U while retaining every charge. This is a lower bound on that one-sided
resilience; minimum blockers of size4 or larger were not searched.

**Incoming improvement, full53 pool only:** the independent frozen certificate
in `../../session47_resilient_native_wiring_20261024/audit-helper-resilience/`
(git27a7118d83) gives at least NINE pairwise-disjoint quadruple witnesses for
every charge. Any eight forced-U owners can hit at most eight witnesses, so
this full53/k4 pool survives all eight such restrictions. It also supplies an
inclusion-minimal46/k3 one-U-resilient bank. Our48-pool bounded hitting result
and actual shared-owner certificate remain separate. The larger one-sided
radius does not improve the arbitrary signed-bit obstruction below.

Arbitrary two-bit resilience is impossible at occupancy4. Requiring two named
helpers at V leaves only C(h-2,2) possibilities, namely1,035 for h=48 and1,275
for h=53, fewer than the2,048 charges. The single arbitrary-bit guarantee is
therefore sharp in the number of arbitrary forced bits. There are also
explicit mixed U/V failures: translate the triple obstruction for
{776,3139} by requiring3139 at V and776 at U. The remaining triple must avoid
both sources, so charge305+label(3139) is unavailable. Exact encoded charges
and witnesses are in the receipt.

## 5. Conditional dynamic-compatibility controller

Apply `../controller/COMPLETE_CHARGE_BANK.md` to the old authenticated28-source
core. At each core arrival, compile its actual interactions with the helper
records. Suppose they permit ALL fixed-occupancy helper subsets except for
the following declared restrictions, with no other changing joint predicate:

1. For the48/k3 bank, at most one helper may be forced U. The forbidden helper
   may be DIFFERENT at every core arrival. Every required charge still has a
   compatible section. A balanced two-core-slot move plus its helper reset
   changes at most8 slots and costs at most50 actual owner-clock occurrences.
2. For the28/k4 bank, the same single forced-U guarantee gives support10 and
   owner-clock bound56.
3. For the48/k4 bank, one helper may be forced to EITHER U OR V, with both its
   identity and required value changing at each arrival. The support bound is
   10 and the owner-clock bound is64. Alternatively, up to three helpers may
   all be forced U at each arrival, under the same bounds.

Proof. At each core arrival choose the certified helper section for its
required charge and current restriction. Consecutive size-r helper subsets
differ in at most2r slots, regardless of whether the forbidden helper changes.
The controller's actual row-current and fixed-margin equations are preserved.
The supplied option resources remain globally disjoint, and the compiled
restriction guarantees every arrival's retained shared pins. The two largest
old core path distances sum to8. The relevant largest helper distances sum
to42,48,56, respectively, giving50,56,64. These are slot-clock OCCURRENCES,
not unmarked graph-edge distances. Supplied feasible full endpoints may be
joined to their chosen sections by a helper-only reset of cost at most2r.
Any requested final helper return must itself be compatible with the final
core records; it is not promised against an inconsistent final pin.

All48 (or28) helper slots remain occupied throughout. A forced-U helper uses
its genuine U path and normalized resource. A charge-dependent choice of a
different witness and a later return to a released option are allowed;
simultaneous resource reuse, double counting at an arrival, or releasing an
old shared pin is not. The guard is the full inherited paired dyadic bank;
an additional exact unscaled integer-X equation would require its own section
test and is not inferred from coverage of Q.

These are native-charge sufficient conditions for the controller's
state-dependent projection-completeness criterion. They do not assert that
arbitrary geometric compatibility constraints have these simple forms.

## 6. Reproduction and complete data

Run from the repository root:

    python 05-knowledge/results/compatible_charged_helpers_20261023/compatibility/resilience.py
    python -O 05-knowledge/results/compatible_charged_helpers_20261023/compatibility/resilience.py

`resilience.json` is the readable summary. `resilience.full.json.gz` preserves
the complete per-charge hitting certificates, all single-pin missing-charge
lists and compact robust sections; its raw hash is in the summary. Generation
uses a deterministic gzip header with mtime0. Verification decompresses and
compares the mathematical certificate, avoiding dependence on a platform's
compression-library byte choices. `--write` regenerates both files.

The full native paths, common pins, kernel/converse resources and source labels
remain linked by exact source hash to `../supply/selected_supply.json`. This
resilience calculation changes none of those records.
