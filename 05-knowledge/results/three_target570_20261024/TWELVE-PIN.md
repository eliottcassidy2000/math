# Late integration: unrestricted occupancy with twelve moving pins

Status: **PROVED finite-supply extension**, using the incoming authenticated
native frames at commit2bc5d106db and the independently checked graded bank.

The incoming
`../guarded_component_completion_20261024/native/NATIVE_COMPONENTS.md`
constructs ten controls on source176 and two on1532. Their total V occupancy
is always six; their total switching costs are61 and11 respectively. It
keeps2948 at U and fixes the protected forty-helper reservoir occupancy at
three, giving a20-source/120-clock controller on that declared family.

Here retain EXACTLY its twelve indexed-owner equations and immutable frames,
keep2948 at U, and remove the separate core/helper occupancy restrictions.
This is a different pin family from `SHARED-PIN.md`: in particular368 has a
different relation to176. We do not conjoin incompatible pin equations or
change a pin during a route. Add the same four fresh570 sources with their
own shared anchor in new owner names. In the incoming host2393 has a private
cell, verified from its full embedding, so its extra570 option creates no
new overlapping-pair condition. All source and partner gluing remains valid.

## Separate a constant population from its changing charge

Remove the twelve controlled bits and the fixed2948 bit from the binary
variables. There remain M=68 or67 free bits (according to2393's W status):
the20-source bank R and F=M-20 other bits. Let k be the number of these free
bits at V. The actual binary V population is k+6, not k.

Write z for the176 bit and t for the1532 bit. The full controlled assignment
C(z,t) has six V records and its actual source charge ell(C(z,t)). Thus the
reservoir section at an outside subset J and phase z must pay

    q - ell(J) - ell(C(z,t(J))),

while containing176 exactly when z=1. The conditional7..13 table covers
every such residual charge. The controlled sources remain occupied and their
cost is counted whenever either root changes.

## Constructive route and budgets

For7<=k<=M-7, first enter reservoir grades7..13 while holding z fixed.
Then normalize to grade7 when k<=F+7 and to grade13 otherwise. This always
leaves a valid outside cardinality, because F>=47. Changes of a free outside
1532 bit also flip its two controls; other outside changes are ordinary
individual changes. At the normalized grade, switch z to zero and choose
the exact new conditional bank record, paying all ten controls. The
normalized outside sets now connect by Johnson exchanges at grade7 or13.
Reverse the second endpoint's normalization to recover all its full records.

The estimates hold for every bank representative:

* Entering the band costs at most63 from the previous argument, plus11 if
 1532 flips: at most74 clocks and16 sources.
* A consecutive-grade move costs at most45+11=56 clocks and22 sources.
* A normalized outside exchange costs at most28+14+11=53 clocks and18 sources.
* Switching z at grade7 or13 changes at most14 bank bits plus ten controls,
 costing at most28+61=89 clocks and24 sources.

No extra phase change is hidden in a Johnson step. Each normalization takes
at most eight events including its phase reset. Outside exchanges need at
most24. Thus the binary route uses at most40 events. The optional four-source
W exchange gives at most41 events for the full three-target fibre.

For k<7, any two free-bit configurations differ on at most12 bits, costing
at most84 clocks. At most all twelve controls also change, costing72.
The full simultaneous endpoint event therefore costs at most24 sources and
156 clocks. For M-k<7 complement the free bits and all controlled bits while
leaving2948 at U. Every control equation survives, and exactly six controls
remain selected, so the same bound applies. This complementation is used for
the distance estimate, not inserted as an unguarded event.

Consequently **every nonempty complete three-target fibre with these twelve
pins and2948 fixed U is connected at24 sources/156 clocks per event**.
The interior k range has the stronger24/89 bound. The result frees helper
occupancies; it is not a cost improvement over the incoming20/120 theorem
on its smaller fixed-occupancy family. The W-factor proof is unchanged because
no W-eligible source is one of the twelve controls or the fixed2948 source.

## Weighted lift and checks

Use visible state(J,z), keep all full compatible bank records, and use the
actual weights of C(z,t(J)) in the marginal. All cross-fibre choices in the
interior routing graph obey24/89. Therefore the conditional-refresh identity
in `WEIGHTS.md` transports a lazy visible chain with **no extra gap loss**.
Boundary full resets may require156 clocks; they are not asserted to fit89.

`twelve.py` loads and hashes the exact incoming embeddings and guard equations,
tests every total free grade and feasible bank occupancy, and evaluates every
indexed-owner equality in global owner coordinates at every arrival. All four
root phases occur. Its2037 endpoint tests contain23769 events;309672 literal
pin checks pass identically in normal/optimized runs. Universal connectivity
uses the proof above, not extrapolation from those tests. Kernel units,
converse clocks and original paired dyadic currents retain the base audit.
