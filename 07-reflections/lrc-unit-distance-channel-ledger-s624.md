# LRC / Unit-Distance Channel Ledger S624

This is the reservation note for the overnight bridge session. The point is not
to force a decorative analogy between LRC and unit distance, but to test a
shared method: start with a state-local frontier, deliberately damage one
retained channel at a time, and ask which damage changes the decision.

The LRC side already has the clean state-local shadow from S599w-x: survival
bitmasks over `Z/(2n-1)`. The recent n=14 work also warns that the shadow alone
is too thin. Owner, carry, pinch, lift/CRT, and Cprime windows are not metadata;
they are part of what makes the quotient proof-relevant.

The unit-distance side has the matching state-local shadow from S617:
frontier-gain edge counting on a Moser carrier. S622 and S623 then show how to
break it productively: starve width, drop directions, cap gains, damage
compactness bias, or weaken canonicalization. The losses name the side
channels.

S624 will try to put those two lessons into one small bridge table. On the LRC
side, the first impairment should be a witness-orbit jackknife over the
`C=2n-1` survival masks. On the unit-distance side, the first comparison should
reuse the Moser direction/gain/width impairments rather than inventing another
large search. The shared language is a damage ledger over channels.

Tournament Analysis should use channel/proof-obligation vertices, not runners,
unit edges, or points. The pairwise observable is whether a channel preserves
frontier decisions under impairment. The switch/gauge can prioritize either
proof relevance or scaling value; any edge flips between those gauges are part
of the result.

The honest risk is that this bridge only restates the parent packets. The way
to avoid that is to compute an LRC jackknife table that looks structurally like
the unit-distance direction jackknife: which omitted witness orbit unlocks a
false survivor, and at which `n` or row family?
