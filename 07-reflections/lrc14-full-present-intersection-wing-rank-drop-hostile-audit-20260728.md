# Hostile audit of the root-zero Mayer--Vietoris wing candidate

**Independent audit, 2026-07-28.  Verdict: REJECT the provisional THM-2751
candidate in commit `7918b065db7`.**  The all-label two-clock rank-one law
survives.  The candidate wing gain does not.

## First false implication

The natural one-sided source and target sheets are unclocked in the physical
`c1`-present coordinate.  THM-2749 freezes that coordinate at `e=1`.  The
candidate silently identifies the resulting fixed-`e` two-sided section with
the full intersection of the unclocked sheets:

```text
THM-2749 fixed-e=1 section = A intersect B.              (false)
```

The correct common object is the disjoint union over every physical-present
clock,

```text
M=A intersect B = disjoint-union_(e in F_7) M_e,
M_e=A_e intersect B_e,                                   (1)
```

after pulling the target sheet into the source coordinate.  All cross-clock
pieces `A_e intersect B_f`, `e != f`, are empty.  Thus the failure is not a
hidden cross-clock correction: it is the omission of the nonzero `e=2,3`
same-clock pieces.

## Object-level checks

The audit constructed the seven exact weighted interval objects before
integration and checked:

- the physical-present sections are pairwise disjoint;
- `disjoint-union_e A_e=A` and `disjoint-union_e B_e=B` as literal weighted
  objects;
- `disjoint-union_e M_e` equals the direct literal `A intersect B` object at
  both endpoints;
- every `e != f` cross-clock intersection is empty; and
- direct literal restrictions to `A minus M` and `B minus M` integrate to the
  same wing coefficients as subtraction of the full common amplitude.

These are equality-of-object checks, not only equality of total masses.

## Exact decomposition at `(s,t)=(0,3)`

Put

```text
C  =339633525654239542165440,
D  =750593782703678965571520,
E  =719200126392878704654080,
Es =722054095148406001101120,
C' =345341652135823400016960,
D' =756301720214733558465600,
Et =724908063903933297548160.
```

The physical-present row amplitudes are

```text
A_e=(0,C, D, Es,0,0,0),
M_e=(0,C, D, E, 0,0,0),
B_e=(0,C',D',Et,0,0,0).                                 (2)
```

Consequently

```text
A=1812281403506324508838080,
M=1809427434750797212391040=C+D+E,
B=1826551436254490256030720,                            (3)

L=A-M= 2853968755527296447040,
R=B-M=17124001503693043639680.                          (4)
```

The literal wing rows are

```text
L_e=(0,0,0,2853968755527296447040,0,0,0),

R_e=(0,
     5708126481583857851520,
     5707937511054592894080,
     5707937511054592894080,
     0,0,0).                                             (5)
```

In particular, subtracting only `C=M_1` from either unclocked sheet is not a
Mayer--Vietoris wing decomposition.

## Root-normalized boundary

After division by the inherited content `26`, reduction modulo `13`, and the
declared endpoint root normalizations, the component rows are

```text
                   source root 12       target root 1
common M_e         (0,4,10,8,0,0,0)    (0,9,3,5,0,0,0)
true wing           (0,0, 0,12,0,0,0)   (0,9,2,2,0,0,0). (6)
```

Augmenting in the physical-present coordinate gives

```text
common: 9 -> 4,       gain 4/9=12=-1,
natural: 8 -> 4,      gain 4/8=7,
wing:   12 -> 0.                                      (7)
```

Equivalently, `v_13(L)=1` whereas `v_13(R)=2`.  The target wing is nonzero
before augmentation, but

```text
9+2+2=13=0 mod 13.                                    (8)
```

This cancellation is the mechanism of the extra target valuation.  After the
delayed-clock vector `(0,1,1,1,1,1,1)` is reduced modulo `Phi_7`, the source
wing profile is `(1,0,0,0,0,0)` and the target wing profile is zero, with
determinants `1` and `0`.  There is therefore no normalized wing unit and no
wing gain `2`.

## Strongest survivor

The correct exact statement is a full-physical-present Mayer--Vietoris
decomposition with a target augmentation rank drop.  The common raw amplitude
still transports exactly and retains normalized gain `-1`; the unsplit
one-sided sheets retain normalized gain `7`.  The wing datum is a
three-component present-clock vector whose target augmentation vanishes, not a
scalar clutch current.

The fixed-`e=1` arithmetic in the rejected candidate may be retained only as a
local subchart comparison.  It is not the complement of the full common
intersection and cannot support a global Mayer--Vietoris or external `C3`
interpretation.

The independent two-clock audit also confirms, for all `81*7*7=3969` cells,

```text
B_(e,j)(s,t)=a_e(s,t) 1_(j != 0).                       (9)
```

That separation supplies a physical-present coefficient row and a delayed
support mask.  It does not supply a lawful order-three action on the present
rows, an external semantic arm, a same-support fixed reference, or an endpoint
current.  In particular, the tempting multiplier `e -> 2e` is not an
established physical operation on this carrier, while `e=0` and `j=0` are
empty structural faces rather than fourth reference states.

## Reproduction and status

The full object unions, cross-clock hostile, amplitudes, valuations, residues,
profiles, and determinant checks are now included in
`04-computation/lrc14_root_zero_full_target_semantic_clutch_20260728.py`.
Run

```bash
python3 04-computation/lrc14_root_zero_full_target_semantic_clutch_20260728.py
python3 -O 04-computation/lrc14_root_zero_full_target_semantic_clutch_20260728.py
```

and compare with
`05-knowledge/results/lrc14_root_zero_full_target_semantic_clutch_20260728.out`.
This audit rejects only the provisional THM-2751 gain-`2` interpretation; it
does not weaken the corrected THM-2754 two-clock theorem.
