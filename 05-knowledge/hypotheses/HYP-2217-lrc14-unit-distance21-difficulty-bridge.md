# HYP-2217: LRC n=14 and Unit-Distance n=21 Share a 27-Quantum Section/Bulk Obstruction

**Status:** OPEN synthesis, supported by S641 and parent results HYP-2164,
HYP-2166, HYP-2167, HYP-2197, THM-408, HYP-2206, HYP-2215, HYP-2216, and
the monad S6 exhaustive `n=9` H-spectrum cross-signal.

## Claim

The connection between the difficulty of LRC `n=14` and unit-distance `n=21`
is not the raw numbers `14` and `21`, and not a literal tournament `H=21`
realization.  The real common object is:

```text
27-quantum + retained section + bulk/lift side channel.
```

On the LRC side, the quantum is the THM-401 shell modulus

```text
C = 2n - 1 = 27 = 3^3.
```

On the unit-distance side, the same `27` is the edge increment per full Moser
slab in THM-408:

```text
E(P_m^+) = 27m + 6,
E(P_m^-) = 27m + 3.
```

The exact unit-distance `n=21` row is `P_2^-`, hence

```text
|P_2^-| = 21,  E(P_2^-) = 2*27 + 3 = 57,
spine = 20,  bulk = 37 = C_hex(3).
```

Thus both problems become hard at a retained-side-channel seam.  The easy
section is already visible:

- LRC: the unit-shell / `Res_27` section and the classified AP, `V*`,
  `2*AP` base quotient;
- unit distance: the unit Hamiltonian spine through Moser slabs.

The remaining difficulty is the small side channel that the scalar quotient
forgets:

- LRC: lift/CRT carry, owner, and Cprime endpoint data;
- unit distance: bulk, direction support, endpoint-compatible ears, embedding,
  and obstruction labels.

## S641 Evidence

`04-computation/lrc14_ud21_difficulty_bridge_s641.py` records the comparison
and stores `05-knowledge/results/lrc14_ud21_difficulty_bridge_s641.out`.

For LRC `n=14`, the script recomputes the AP `C=27` orbit packets:

| gcd class | orbit size | shell reps | released residues | depth histogram | redundancy price |
|-----------|------------|------------|-------------------|-----------------|------------------|
| `1` | `18` | `9` | `18` | `{1:18}` | `18` |
| `3` | `6` | `3` | `6` | `{1:6}` | `6` |
| `9` | `2` | `1` | `2` | `{4:2}` | `8` |

The tiny gcd-9 channel has only two residues but depth four, so it is small
mass with high proof relevance.  This is the LRC prototype for a side channel
that scalar counts underprice.

For the Moser side, S641 verifies the named THM-408 rows:

```text
P_1^+ : vertices=14, unit_edges=33, spine=13, bulk=20
P_2^- : vertices=21, unit_edges=57, spine=20, bulk=37
P_2^+ : vertices=22, unit_edges=60, spine=21, bulk=39
```

The exact `n=21` row is therefore the second minus slab:

```text
57 = 2*27 + 3 = 20 spine + 37 bulk.
```

## H-Spectrum Cross-Signal

While this bridge was being written, monad S6 pushed an exhaustive
isomorph-free tournament H-spectrum at `n=9`:

```text
iso classes = 191536
distinct H values = 1520
H range = [1, 3357]
low odd gaps <=609 = [7, 21]
```

All the previous `n=8` high gaps above `609` fill in at `n=9`, but `7` and
`21` remain the only low gaps.  This does not make unit-distance `n=21` a
literal realization of forbidden tournament `H=21`; in fact the bridge above
explicitly avoids that interpretation.  It does say the visible `21` is not
just decoration: it is a durable obstruction value in the tournament channel,
while the LRC/unit-distance difficulty comparison is controlled by the
`27`-quantum section/bulk ledger.

## Dictionary

```text
LRC C=27 shell quotient       <-> Moser full-slab edge quantum 27
LRC unit-shell section        <-> unit Hamiltonian spine section
LRC nonunit/carry bulk        <-> unit-distance tile/bulk and ears
lift/CRT conservativity       <-> endpoint-compatible ear conservativity
gcd-9 tiny high-depth channel <-> direction/gain packets with small mass
```

This dictionary says why the two problems feel similarly difficult: both have
a clean quotient that solves most of the visible structure, followed by a small
but load-bearing residual channel.  In both settings, widening the raw search
is less valuable than proving the retained channel is conservative.

## Tournament Analysis

S641 uses proof routes as vertices.  The bridge tournament is transitive with
one Hamiltonian path:

```text
27_quantum_section_bulk_ledger
lift_ear_conservativity_transfer
side_channel_jackknife_damage_table
same_scalar_twin_search
centered_hex_vs_gcd9_microchannel
raw_14_21_number_match
```

The ranking says the transfer should not be raw `14/21` numerology.  The first
usable route is a channel-complete ledger that keeps the `27` quantum, the
traceable/clock section, and the bulk/lift side channel simultaneously.

## Next Tests

1. Extend S624 from witness-orbit deletion to lift/carry deletion: impair owner
   route, carry cocycle, and pinch status over the same `C=27` rows and measure
   which false floor rows return.
2. Extend the unit-distance `21`-core ledger by endpoint-compatible ears:
   delete direction support, bulk class, spine endpoint, and gain packet
   channels one at a time.
3. Build scalar twins: same `27` quantum but different retained side channel;
   same `57` unit-edge count but different spine/bulk/direction packet.
4. Try a common conservativity statement:

```text
If the retained section is fixed and all local channel jackknives are stable,
then the only possible frontier failures lie in explicitly named lift/ear
cocycles.
```

## See Also

`04-computation/lrc14_ud21_difficulty_bridge_s641.py`;
`05-knowledge/results/lrc14_ud21_difficulty_bridge_s641.out`;
`07-reflections/lrc14-unit-distance21-difficulty-bridge-s641.md`;
HYP-2216; HYP-2215; HYP-2206; HYP-2197; HYP-2188; HYP-2167; HYP-2166;
HYP-2164; HYP-2200; THM-408; THM-407; THM-401; THM-115.
