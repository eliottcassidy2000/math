# q3/q11 collar–endpoint bridge: exact scratch report

Status: **FINITE-EXACT SCRATCH**. This is not a canon promotion, not a proof
dependency, and not an LRC(14) conclusion. It composes the proved THM-2825
collar, the proved THM-2847 q3/q11 endpoint horn, and the proved THM-2851 carry
triangle. THM-2857 is used only as a pinned **proof-candidate comparison
fixture**.

## Exact universe and inheritance

The audit reconstructs all 193 semantic cells and all 587 rooted collar paths.
Its focused universe is the 42-cell q3/q11 intersection (114 roots), including
the 20-cell q7 horn (52 roots). Hostile controls include every non-distinguished
root, the 22 non-horn intersection cells, short paths, both endpoint frames,
the h=9 and h=4 legs, the step-90 boundary, two finite-field embeddings, and
normal versus `python -O`.

In every one of the 42 intersection cells there is exactly one labelled
index-zero root whose second collar iterate is

```
I = M2 = (142004992589460, 142005019034340).
```

For all 20 horn labels these roots collapse to one physical ladder

```
R  = (142004190428100, 142004216872980)
M1 = (142004591508780, 142004617953660)
M2 = I.
```

The semantic word is live/dead/live. At `R`, E3 alone is missing; at `M1` and
`M2`, all six physical factors hold. All 20 horn paths survive at least 14
steps. Thus the collar supplies a genuine labelled source section, but its
birth is locked to E3 repair.

## Positive physical first-carry reference

The source endpoint mask at `R` is empty. The masks at `M1` and `M2`, and all
three target-frame masks, are the same 81-point Cartesian rectangle. Under
`iota(v,w)=v+13w`, cyclic successor is `Y` away from `v=12` and `ZY` at the
wrap; the high carry is the vertical endpoint translation `Z=(0,1)`.

On exactly 12 long horn labels,

```
s in {0,3,12},  t in {5,6,9,10},
```

the physical common-tangent move from collar step 2 to step 68 is `+66h` and
acts on both endpoint frames by `(0,8)=Z^8`. It preserves:

- all six physical factors in both frames;
- the source and target carrier delta-zero masks;
- the live semantic value;
- the ancestry chamber `(14,8)`; and
- the literal ancestry contributor sets, of sizes `(966606,28534)`, with
  symmetric differences `(0,0)` and digest
  `15c804c7cea9f61feab3b641eccdc035d937142b446d1cc14e059210eb1534fd`.

Since `gcd(8,13)=1`, `Z^8` is a faithful generator of the first-carry
coordinate. The labelled rank is 12, while the physical interval quotient has
rank 1. Carry faithfulness is therefore no longer the missing bridge.

The full plateau is slightly asymmetric and this is a useful boundary:
the source-frame `(0,8)` translate persists for steps 68 through 90, while the
target-frame translate persists only through step 89. All factors, carrier,
and ancestry chamber persist at all 276 labelled occurrences. Semantics
alternate: 144 even-step occurrences are live/live and 132 odd-step
occurrences are live/dead. Step 90 is live in the source frame but loses the
target-frame endpoint translate in all 12 labels. Hence step 68 is the first
semantic-preserving two-frame reference and step 88 is the last such reference
on this plateau.

## Character comparison

For the `+66h` move, exact cyclotomic powers are

```
X0*R*Delta/(N/13) = 396 = 6 mod 13,
Y0*R*Delta/(N/13) = 1254 = 6 mod 13.
```

The signs of the actual currents matter. Exact evaluation in both pinned
finite-field embeddings gives:

- source `x_sweep` phase `omega^7 = omega^(-6)`;
- target `endpoint_sum` phase `omega^6`.

Under the THM-2857 candidate reparameterization `r=10b`, the physical shift
`b -> b+8` gives `Delta r=2`, so its centered character 3 predicts
`omega^(3*2)=omega^6`. The target current matches this candidate character,
while the source current carries the inverse character; their paired phase
cancels. This is an exact one-sided character match, not the missing
semilinear clutch, and THM-2857 remains only a comparison fixture here.

## Exact obstruction to calling this a q-arrow

The pulled full endpoint masks on q0, q3, q7, q11 have masses
`81,90,0,81`. The q0 and q11 masks are not translates in either endpoint
frame. Their convolutional difference has rank 168 and augmentation zero:
only the trivial character dies. The q7 mask is empty and its boundary from
q11 has rank 169 and augmentation `-81`.

Therefore `(0,8)` is an **endpoint-address arrow**, not an allocation-q arrow.
The physical first-carry reference lives inside the E3 block, whereas the q7
horn requires the off-diagonal state `(source present, E3 absent)`. Current
endpoint allocation annihilates the whole q7 address plane, not merely the two
origin columns seen by the transverse minor.

The complete directed-path scan reinforces the typing. Among h=8, h=9, h=4
translations on both coordinate axes, only vertical h=8 occurs. It occurs
only on the distinguished roots (6624 common-intersection pairs, 3312 horn
pairs). No h=9 leg, h=4 leg, horizontal leg, or 8–9–4 triangle occurs.

## Typed bridge contract and remaining gate

- **Source:** the 12 labelled step-2 horn sections over the single physical
  interval `I`.
- **Target:** their step-68 endpoint-address sections.
- **Map:** physical translation by `+66h`, inducing `Z^8` on both endpoint
  masks.
- **Preserved:** six-factor support, carrier, semantic truth, ancestry
  chamber, literal ancestry, and faithful first-carry coordinate.
- **Not supplied:** any map from endpoint addresses to q3/q11/q7 allocation,
  any E3-to-complement transporter, or any macro semantic-word intertwiner.
- **Destroyed by the present q7 constructor:** the entire endpoint address
  plane.
- **Needed sidecar:** a lawful q-coupling that retains the source section while
  moving into the E3-absent block, equivalently a new q7 endpoint support or an
  E3/complement transporter compatible with `Z^8`.
- **Cheapest decisive next test:** enumerate admissible endpoint
  pattern/owner recharts on these same 12 labels and reject any candidate
  unless q7 becomes nonempty, its coefficient current (not just support)
  commutes with the physical `Z^8` reference in both frames, and the q3/q11
  transverse minor survives.

## Reproduction

```bash
PYTHONHASHSEED=0 python3 \
  .scratch/lrc_q3q11_collar_bridge_20260728/cell_path_scout.py
PYTHONHASHSEED=0 python3 \
  .scratch/lrc_q3q11_collar_bridge_20260728/endpoint_carry_intertwiner_audit.py
PYTHONHASHSEED=0 python3 -O \
  .scratch/lrc_q3q11_collar_bridge_20260728/endpoint_carry_intertwiner_audit.py
```

The stored normal and optimized endpoint-audit outputs are byte-identical.
