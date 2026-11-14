# A nonnegative seam defect and exact small-gap reader obstruction

**Proved:** the exposed finish response cannot cancel a small-gap
contact cost, even if that response is allowed to vary with the
source. A signed version of the ordinal transfer hides this fact.
The restriction is on jointly specified absolute coordinates, not on
the scalar Chan–Pak spectrum.

## The disjoint finite supports

A puncture B is a nonempty tournament, with one external apex v.
Let U be the old vertices beaten by v. Write h=H(B), K=H(B+v),
and P for the number of ordered pairs of complementary paths whose
right path is nonempty and starts in U. The left path may be empty
(and then has exactly one choice). This is the unsigned start response
of the [common-apex gluing law](../chanpak_movable_reader_returns_20261114/INDEPENDENT_RETURNS.md).
Use R_+=h+P. Let S count full old paths ending in U and let T
count full old paths starting outside U (nonnegative conventions).

Write E for ordered pairs of complementary **nonempty** paths with
left path ending in U and right path starting in U. Exactly these
pairs in P fail insertion of the apex: the left endpoint would
need to beat the apex. Exactly S full old paths cannot place the
apex last. Partitioning literal augmented paths at the apex gives

```
L = R_+ − K = E + S.
```

This uses individual labelled paths, not cancellation of a
signed count. Likewise M=f−h−Q counts complementary pairs whose
first path is nonempty and ends in U; the second may be empty.

For C → D (all arcs oriented from C to D) with the same apex,
split the pairs counted by E according to the block holding the
last vertex of the left path and the first vertex of the right path:

| Left ends | Right starts | Exact number |
|---|---|---:|
| C | C | h_D E_C |
| D | D | h_C E_D |
| D | C | P_C M_D |
| C | D | S_C(h_D−T_D) |

Each term is a Cartesian product of choices. For example in the
third row the right portion in C and left portion in D must be
nonempty, but the other two portions may be empty. The ordinal
orientation supplies all needed edges. Since S_CD=h_C S_D,
this proves the nonnegative identity

```
L_CD = h_D E_C + h_C L_D + P_C M_D + S_C(h_D−T_D).
```

In particular the negative term −S_C T_D in the signed expression
for L_CD can at most consume the S_C part of h_D L_C. It cannot
offset P_C M_D. This is why the field-with-one-element viewpoint
helps here: retain the actual sets and the cancellation vanishes.

## Exact budget on a right reader

If D has r apex outcontacts, each nonempty subset of those vertices
admits a path, as does its complement (with the empty convention).
So M_D ≥ 2^r−1. For r≥1 equality holds exactly when D is transitive
and the contacts form an initial segment of its order. Indeed a
reverse cross arc from outside the contacts into them gives an
extra two-vertex subset with a path ending in the contacts. With
no reverse cross arc, the subset of all contacts contributes
H(U)H(D−U), which must equal one. A tournament has count one
exactly when transitive: a directed triangle contributes a positive
extra term in the odd-cycle packing expansion; acyclic tournaments
have a unique path.

Thus for prescribed d=R_+−H, any source C with P_C>0 permits only

```
r ≤ floor(log2(1 + d/P_C)).
```

If 0≤d<P_C, the right reader has **no** such contacts. It can have
arbitrary internal geometry, but all of its vertices beat the
apex. The transfer reduces to

```
d = h_D E_C,                 R_+ = h_C h_D + P_C f_D.
```

For fixed C and d>0 this gives an exact decision procedure.
Require E_C>0, divide d by E_C for h_D, divide
R_+−h_C h_D by P_C for f_D, and reject unless both are positive
integers. Accept iff (h_D,f_D) is a genuine signature of a
nonempty tournament. That last check is finite since
f_D=Σ_A H(D[A])H(D−A)≥2^|D|. For d=0, require E_C=0, and test
the finitely many positive h_D≤R_+/h_C instead. These tests are
sufficient as well as necessary; no residue lift substitutes for
the absolute coordinates. The threshold is sharp: take C to be
the transitive word 1 followed by n≥1 zeroes (P_C=2^n), and D
the word 10. Then E_C=S_C=L_D=0, M_D=1, so d=P_C with one
reader outcontact.

## Published base-five bank, and a literal escape by changing source

For the [joint bank](JOINT_FIVE.md) before its optional reader,
S_C=0, h_C=5^k, P_C even and at least 40,960 for k≥1.
Its seed has L=2, and each digit gives L_new≥5L_old by the
identity. Thus E_C≥2·5^k. For an odd prime p with
2p<40,960, prescribing d=2p forces h_D=1: the only other
positive odd divisor p would require E_C=2. The nonempty
reader has even f_D by swapping complementary subsets, so
R_+≡1 (mod 4). All prescribed pairs with d=2018 and
R_+≡3 (mod 4) are therefore out for this bank, even over
**every** ordinal right reader with arbitrary contacts.

The exclusion is not a universal tournament obstruction. Change
the source to a 12-vertex transitive puncture with contact word
100101101011 from first to last (1 means beaten by the apex).
Its exact (h_C,P_C,E_C,S_C) is (1,2411,2018,1).
Append one singleton beating the apex: (h_D,f_D)=(1,2).
This gives the strong 14-vertex augmented tournament
(H,R_+)=(2805,4823), with d=2018 below P_C. Changing the source
works here; retaining it and relying on a cancellation at the
reader cannot work. Neither this pair nor the reader obstruction
settles absolute coverage of H alone.

The initial arbitrary-reader calculation was found in a background
pass; the disjoint-support equation strengthens it to arbitrary
S_C, gives the full contact budget and exhibits the changed-source
escape. Run **python3 replay_ordinal_defect.py** for the literal
checks, including 5,476 small joins and the 14-vertex witness.
