# A source that requires changed old edges, and two literal native escapes

**FINITE-EXACT** source/operation certificates; **PROVED** all-order reduction
using the cited strong-order bound. No reusable two-offset return is claimed.

## E1. The chain has source geometry

Starting with the unique strong four-vertex H5 class, complete two-ear
enumeration gives two strong six-vertex H17 classes. Continuing gives six
strong eight-vertex H53 classes. Exactly two of those six have strong two-ear
extensions to H161; they both reach the same ten-vertex class. Neither gives
a strong H163 extension. The two successful H53 sources are the converse pair
T_(3,1),T_(1,3) in WINGS.md. Every relabelling permutation and incidence mask
along the successful chain is saved in `branches.json`.

Thus W5 explains the stages17,53,161 by the actual parameter changes
(1,1)->(3,1)->(3,3). The preceding5->17 construction is also saved literally.
All four displayed stages are strong. A temporary intermediate after the
first of a two-birth pair need not be strong.

The final canonical source is

    U=(2,16,3,7,45,271,63,127,735,255),   H(U)=161.        (E1)

## E2. No induced supertournament of U has H485 or H487

The complete two-birth cut calculation finds neither target among any of the
2^10*2^11=2097152 labelled extensions, without a strongness filter. The replay
prunes a first birth only when its H already exceeds487: H is nondecreasing
under induced extension, since each old Hamilton path has an insertion of a
new vertex and deletion recovers the old path. This is a proved exhaustive
filter. One-birth extensions are also excluded, since a universal extra
vertex would preserve their H and make a two-birth witness.

[Busch2006, Corollary1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf)
gives minimum strong H at order13 equal to625, and a larger minimum at every
subsequent order. Thus no strong extension with either desired H exists at
any order. This applies the published theorem; the strong minimum is not a
new result of this packet.

Nor can arbitrarily many weak exterior vertices help. In any extension, all
of the strong U lies in one strong component C. Insertion monotonicity gives
H(C)>=161, and H(C) divides the full H by the SCC product law. The only divisor
of485=5*97 or487 at least161 is the target itself (487 is prime, checked by
trial division through22). Hence C would be one of the impossible strong
extensions. This proves the unrestricted induced-extension statement.

This is a source obstruction, not a missing spectrum value. Both values have
the explicit strong realizations in the next paragraph.

## E3. One actual path complement plus one birth realizes each target

For a directed Hamilton path p, the native complement keeps its consecutive
arcs and reverses every other arc. The newborn mask lists old vertices that
point to the new vertex. All old labels below refer to E1.

| Target | Selected legal path p in U | H after complement | Birth mask |
|---|---|---:|---:|
| 485=3*161+2 | (3,2,0,1,4,5,8,9,7,6) | 25 | 153 |
| 487=3*161+4 | (2,1,4,5,8,9,7,6,3,0) | 169 | 22 |

Both final tournaments have eleven vertices and are strong. The complete
intermediate and final rows are in `native_escape_161.json`. The mask153 has
in-neighbours {0,3,4,7}; mask22 has {1,2,4}. Guarded complement replay, subset
counting, and an independent depth-first Hamilton-path enumeration agree.
These are native complements; no claim about their black/blue colour is made.

The physical old root0 is retained. Its F49/charge triple changes from
(1,6,1) to (1,6,3) on the485 branch and to (6,5,6) on the487 branch.
Thus even the branch preserving alpha changes the charge g; neither is
silently called phase-neutral. The complete integer observers are in the receipt.

For contrast, none of U's161 legal single native complements preserves H161.
Among all45 one-arc flips and990 two-arc flips, respectively four and eight
preserve161, and every such target is isomorphic to U. Thus those small
H-neutral local searches do not supply a reset. They do not exclude larger
edge changes or multiple native moves.

There is an all-size continuation after each escape. The H487 target has
F_r(2)=H at root2. The converse of the H485 target has the same zero guard at
root9. The inherited two-birth reset sends any such rooted source to a strong
zero-guard source with H'=3H. It therefore realizes strong counts

    485*3^j and 487*3^j,  j>=0, at order11+2j.

This zero continuation is credited to the prior reset/Session163 mechanism.
It does not keep adding offsets2 or4. The new positive step here is the
actual native escape into those two infinite continuations.
The incoming canonical return calculus further supplies multipliers3 and5,
with exact orders, as proved and tested in [INCOMING.md](INCOMING.md).

## E4. Same H and F49 phase do not determine the next attachment language

The bad H53 source, in a named frame, is

    B=(250,40,51,52,2,16,62,126).

The good source is T_(3,1) of W5 in natural block order. At named root7 their
golden observers are respectively

    (53,225,226) and (53,204,205).

Both have the F49 chart (a,b,g)=(2,2,2). Nevertheless B has no two-birth
extension to either161 or163, while the good source has a strong161 extension.
This finite hostile is reproduced from both literal relations in `receipt.json`.

Because the entire three-coordinate observers agree modulo7, every word of
the five inherited Q4 matrices continues to agree modulo7 as well. Their
integer difference (0,21,21) is retained. This does not claim equality of
their integer Q4 responses, joint tensors, or legality predicates. It shows
why a correct F49 transport coordinate is a sidecar, not a replacement for
the source and its jointly observed ports.

Next targets are now concrete: lift the two displayed native paths to a
parameterized family with endpoint positions retained; or find a coupled
wing operation that changes the constant in W2 while satisfying its integer
realization guards. Neither task is settled by the finite escapes.
