# Unrestricted routing also survives the literal moving176/368 pin

Status: **PROVED finite-supply extension**, using **FINITE-EXACT** conditional
sections and the incoming authenticated two-cell gluing.

The base theorem retains the inherited fixed pins. We now ADD the genuinely
changing equation from
`../compatible_charged_helpers_20261023/supply/shared_owner.md`:

    path_176[2] = rho(path_368[1]),
    rho=(8,9,10,2,6,11,12,13).

On the actual U/V records its four cases give owner pairs (2,6),(2,2),(6,6),
(6,2), so compatibility is exactly b_176=1-b_368. The two source cells share
the coherent oriented pair6->2; it is a non-path pair for every supplied
option, hence both local partner matrices reverse it. `guarded.py` checks
these literal equations. The incoming fourteen-owner local-host realization
is retained; no claim about a global Hamiltonian path on the combined host is
made. The other selected records remain in the same compatible independent
context. Arbitrary additional cross-cell constraints are outside this result.

The20-source reservoir contains176;368 is outside. At each outside arrival
the required176 bit is known. Remove176 from the bank and enumerate all2^19
subsets. The deletion bank covers all2048 charges at every grade6 through13.
Its grades5 and14 cover only2045. Therefore BOTH values of176 have complete
conditional charge coverage at every TOTAL reservoir grade7 through13.
This retains the shared pin as an equation, rather than ignoring it and
checking charge afterwards. The missing grade5 charges explain why the
unconditional5..15 band cannot simply be reused with the new pin.

Repeat `THEOREM.md`'s normalization with band[7,13]. For7<=k<=n-7, entering
the band from r<7 changes at most r+7 reservoir slots and7-r outside slots,
costing at most

    2(r+7)+7(7-r)=63-5r <=63.

The upper boundary follows by complement. Consecutive-grade events still
cost at most20 slots/45 clocks. Horizontal grade7 exchanges cost at most
16 slots/42 clocks. A within-fibre refresh costs at most20 slots/40 clocks.
Every conditional table used by these steps is complete, regardless of
whether the outside change flips368 and hence the shared owner.

At k<=6, retain the fact that exactly one of176 and368 is selected. A direct
DP on the other n-2 actual sources, followed by adjoining either required
source, gives the COMPLETE attainable charges and a cheapest witness. For
both n=80 and81 the maximum minimum witness costs at grades1,...,6 are

    6,13,20,13,12,13.

Thus an arbitrary endpoint of grade k can reach that witness within
7k plus the displayed cost, at most55 clocks. The number of changed sources
is at most12. Complementation handles k>=n-6; it preserves the opposite-bit
equation. The exact DP is over every actual source, not a sampled charge table.

Consequently every nonempty three-target fibre WITH this moving-owner
equation is connected at **20 changed sources and63 owner clocks per event**.
There are at most44 binary events (seven normalization events at each end
and at most30 outside swaps), and at most45 including the optional W event.
No helper occupancy is held fixed. All original source records and fixed
pins remain present. This extends the incoming fixed-occupancy10/64 theorem
to unrestricted occupancies; it does not improve its smaller support bound
on that restricted state space.

The conditional refresh theorem in `WEIGHTS.md` applies to this same projected
graph with band7..13 and the actual compatible reservoir fibres. On the
interior grades it retains the prescribed measure with at most20 slots/63
clocks per transition. Boundary fibres still require their own sparse
weighted graph; complete unconstrained refresh is not inferred from a route.

Kernel lifts remain one capacity unit. Converse reverses both paths and
changes the displayed equality to clocks5 and6, exactly as in the incoming
certificate. A changed shared owner is synchronized at each accepted arrival;
unaccepted one-slot intermediates are not introduced.
