# LRC / Unit-Distance Channel Ledger S624

The point of this bridge session was not to force a decorative analogy between
LRC and unit distance, but to test a shared method: start with a state-local
frontier, deliberately damage one retained channel at a time, and ask which
damage changes the decision.

The LRC side already has the clean state-local shadow from S599w-x: survival
bitmasks over `Z/(2n-1)`. The recent n=14 work also warns that the shadow alone
is too thin. Owner, carry, pinch, lift/CRT, and Cprime windows are not metadata;
they are part of what makes the quotient proof-relevant.

The unit-distance side has the matching state-local shadow from S617:
frontier-gain edge counting on a Moser carrier. S622 and S623 then show how to
break it productively: starve width, drop directions, cap gains, damage
compactness bias, or weaken canonicalization. The losses name the side
channels.

S624 put those two lessons into one small bridge table. On the LRC side, it
used a witness-orbit jackknife over the `C=2n-1` survival masks. On the
unit-distance side, it recomputed a lightweight Moser direction/gain impairment
and leaned on the stored S622/S623 atlas for the heavier target-10/14 facts.
The shared language is a damage ledger over channels.

The LRC result is pleasantly sharp. At `n=14`, the AP row, repo `Vstar`, and
doubled AP are all blocked on `C=27`, and their orbit damage profiles are the
same: omit the gcd-1 orbit and release 18 shallow residues, omit gcd-3 and
release 6 shallow residues, omit gcd-9 and release only 2 residues, but each of
those two is killed four times in the intact ledger. The redundancy prices are
therefore `18`, `6`, and `8`; gcd-9 beats gcd-3 once redundancy is retained.
Small mass, high redundancy: that is exactly the kind of side channel that a
scalar survivor count can underprice.

The unit-distance result rhymes with it. At Moser target `9`, width `30`, the
full shell reaches `18` edges. Direction-drop loss is `{0:3, 2:6}`: three
antipodal directions are replaceable and six are edge-critical. The replaceable
directions have usage `0` in the full cluster; each critical direction has
usage `3`, loss `2`, and usage-loss price `6`. Gain ceilings `1,2,3,4` give
`9,15,16,18` edges. Again the important feature is not only the scalar edge
total; it is which direction/gain packets the construction cannot afford to
forget.

Tournament Analysis used channel/proof-obligation vertices, not runners, unit
edges, or points. The pairwise observable was whether a channel preserves
frontier decisions under impairment. The two switches were proof relevance and
scaling value. Both tournaments were transitive with one Hamiltonian path, but
they flipped `19/36` pairwise edges. That is the practical moral of the
session: the fastest frontier quotient is not automatically the proof quotient.

The next bridge test should move from witness-orbit deletion to lift/carry
deletion: keep the same `C=27` scalar row, but impair owner route, carry
cocycle, or pinch status and measure which false floor rows reappear. On the
unit-distance side, the matching test is the side-channel-complete `21`-core
ledger from HYP-2196: impair direction support, gain packets, canonical orbit
class, deletion resilience, and obstruction labels one at a time.
