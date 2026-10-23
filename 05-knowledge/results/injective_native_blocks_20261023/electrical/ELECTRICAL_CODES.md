# Exact additive codes for paired native trades

23 October 2026. **PROVED:** the finite-group compression argument and the
consequences of the authenticated native coordinates below. **FINITE-EXACT:**
the original order-eight dyadic filtration, code maps, short charge relations
and literal image certificates. **OPEN:** a complete native wiring or move
classification. Electrical admissibility does not allocate actual ports.

## 1. Inheritance and the changed question

The closest proved mechanism is the
[Session44 image/annihilator quotient](../../session44_guarded_native_moves_20261023/electrical-duality/UNIVERSAL_ELECTRICAL_GUARDS.md).
The canonical hostile is the
[two scarce-path rectangles](../../session44_guarded_native_moves_20261023/rectangle-bank/THEOREM.md):
the chosen primitive character accepts amplitude one, the full modulo-eight
bank requires amplitude two, and distinct-port capacity forbids it.
The corrected near miss is replacing an additive charge quotient by a
short label that merely separates individual vertices. The least-used
sidecar is the full integer image witness for a *sum* of source charges.

The anchor is simultaneous cancellation across several source pairs.
The niche is multiplication by the actual paired-black coefficient two.
The wildcard is keeping the entire native dyadic tower rather than only
modulo eight. The concept board is: additive charge, paired multiplicity,
source circuit, frozen per-slot phase, literal capacity, and image witness.

Here `L` is the ORIGINAL 3,528-node converse-merged native Laplacian, with
root 3,527. A selected ordinary-source port together with its actual
converse partner changes black conductance by two. Let `K` be this graph's
critical group, and `C_m=K/mK`. For a zero-sum integer charge `d`, write
`[d]_m` for its class. The universal electrical test is

    d.x=0 for every grounded harmonic x modm
       iff [d]_m=0 iff d in im(L modm).

The new gate needs **`2d`**, not `d`. Its exact effective charge group is

    C_m / C_m[2], isomorphic to 2 C_m.                  (1)

This is where the genuine factor two removes most coordinates. It does
not permit division by two in the original modular equations.

## 2. Authenticate the small group before using it

The checker replays the original dyadic producer against every native
edge. It reconstructs the same owner-marked lifts and filtration, not just
the three exceptional vertex aliases. The successive kernel dimensions
at levels two, four and eight are

    3363, 5, 1.

Thus `|C_8|=2^3369` and `|C_8[2]|=2^3363`, so the group (1) has exactly
64 elements. The complete two-primary group, with its original basis
data retained, is

    K_(2) = (Z/2)^3358 + (Z/4)^4 + Z/256.              (2)

Consequently `2 K_(2)` has 2,048 elements and is isomorphic to
`Z/128 + (Z/2)^4`. These invariant factors alone do not tell us the code of
an actual owner. The following authenticated potentials supply that map.

Let `Y0,...,Y4` be the five original `radical_lifts_mod4`; let `Z` be
`eighth_radical_lifts_mod8[0]`; and let `X` be the saved primitive
modulo-256 lift. They are grounded and verified by literal multiplication:
`LYj=0 mod4`, `LZ=0 mod8`, `LX=0 mod256`. All coordinates below are in the
original owner indexing. For every native vertex `v`, the receipt verifies

    Z_v = Y0_v+Y4_v mod2,
    X_v = Z_v+2(Y0_v+Y3_v) mod4.                      (3)

Define the two owner codes

    sigma(v) = (Z_v mod4, Y1_v,Y2_v,Y3_v,Y4_v mod2),
    Sigma(v) = (X_v mod128, Y1_v,Y2_v,Y3_v,Y4_v mod2). (4)

For charges, extend these maps **additively over the integers**. The root
has code zero. All 64 possible `sigma` values occur at actual vertices.
For `Sigma`, the actual vertices have 1,695 distinct values, but their
additive subgroup is all 2,048 values: the five root-relative vertex
charges at nodes `1,2,4,5,8` already generate it. The difference between
1,695 vertex values and 2,048 group elements is part of the representation.

**Exact native code theorem.** For every zero-sum integer charge `d`,

    2d in im(L mod8)   iff sigma(d)=0;                (5)

and

    2d in im(L mod2^h) for EVERY h>=1
                      iff Sigma(d)=0.               (6)

**Proof.** Pairing `2d` with `2Yj` and `Z` modulo eight yields precisely
the five binary sums and the quaternary sum in (4); relation (3) removes
the redundant `Y0` coordinate. Therefore a zero charge in (1) has zero
`sigma`. The latter map is onto a group of size 64, exactly the size of
(1), so this induced surjection is an isomorphism. This proves (5).

For modulus 256 use the authentic harmonics `64Yj` and `X` instead.
They produce `Sigma`, with the same redundant binary coordinate removed.
The generated image has size 2,048, exactly `|2 K_(2)|` by (2), so the
same argument proves the modulus-256 version. Since the full two-primary
exponent is 256, it is equivalent to every dyadic modulus. QED.

This is a cardinality-and-surjectivity proof on additive charge groups.
It does not infer generation of all harmonic characters merely from
separating vertex pairs. Odd-primary characters and preservation at all
integer moduli are NOT settled by (6).

## 3. Executable gate and its precise interface

For distinct source owners `S` and target owners `T`, with `S` disjoint
from `T`, let `M` be signed new-minus-old paired-port incidence. Require
its actual integer row and column sums to be zero. The associated
Laplacian change has off-diagonal blocks `-2M,-2M^t` and zero diagonal
blocks. Therefore universal preservation is equivalent to:

    every source row has zero additive target-code sum;
    every target column has zero additive source-code sum.     (7)

Use `sigma` for the modulo-eight bank, or `Sigma` for all dyadic banks.
This is an immediate application of (5),(6) to each literal row of the
Laplacian change. The general transport and completion theorem is in
[the sibling two-sided lift packet](../dyadic/TRANSPORT_LIFTS.md).

`check.py` exports `add_code` and `gate_trade`. Pass the saved `codes`
with cyclic modulus four, or `all_dyadic.codes` with cyclic modulus 128.
The routine rejects overlapping owner blocks, repeated labels, bad shapes
and nonzero integer margins. It returns ALL row and column syndromes,
not merely one total checksum. Its arithmetic cost is linear in the
number of matrix entries times five small coordinates.

| Item | Exact interface |
|---|---|
| Source | Signed table of actual owner-labelled paired-port transitions |
| Target | Five-coordinate additive syndromes on every source row and target column |
| Preserved | The specified full original harmonic bank remains harmonic |
| Lost | Actual paths, index pins, per-slot phases, unused capacity and global assignment |
| Needed sidecar | Authenticated literal before/after records, chosen per-slot phase, injectivity and converse partners |
| Cheapest hostile | A zero source-total checksum with a nonzero target-column syndrome |

Passing (7) means inclusion of the ORIGINAL harmonic bank in the perturbed
one; it does not prove equality of kernels or identify a new native wiring.

## 4. Why the original four sources cannot pay unit-capacity debt together

In order `(1238,2242,2246,2315)`, their modulo-eight effective codes are

    (3;0000), (0;0010), (3;0011), (2;0101).

For an integer balanced charge with coefficients `(a,b,c,d)`, its code
vanishes exactly when

    a,b,c,d are all even, and a+c=0 mod4.              (8)

To see this directly, the three active binary coordinates force `d`, then
`c`, then `b` even; balance forces `a` even. The remaining cyclic equation
is `3a+3c+2d=0 mod4`, which gives the second condition. This proves (8)
for ALL integer coefficients, beyond the checker's bounded controls.

Hence there is no nonzero relation using only coefficients `0,+1,-1`.
Any universally modulo-eight-safe bipartite trade supported on these four
sources has **every entry even**, irrespective of how many target pairs
or rectangles are combined. Under a unit net capacity at every entry,
only zero incidence change remains. Actual literal capacities are still
needed to justify such unit bounds; a generic black-edge weight does not.

The smallest nonzero balanced integer relations, up to sign, have l1 norm
four:

    2(e_1238-e_2246),       2(e_2242-e_2315).

For the full dyadic code the exact condition strengthens to

    a,b,c,d even;  d=-a-b-c;
    13a-66b-19c=0 mod128.                              (9)

The smallest l1 norm is eight, attained only up to sign by
`4(e_1238-e_2246)`. The finite minimality check exhausts every balanced
integer vector with l1 norm at most eight; (9) supplies the general kernel.
These are source-column statements, not sufficient target-row tests.

## 5. Expanded source circuits pay the full dyadic debt

The native search lane supplied three charges:

    w1 = e_1477+e_2837-e_1425-e_1921,
    w2 = e_1477+e_1841-e_1394-e_2185,
    v  = e_1965-e_2234.

All three have zero `Sigma`. In fact their `X` pairings vanish over the
integers:

    152+223=158+217,
    152+99=113+138,
    X_1965=X_2234=146.

The difference

    w6=w1-w2=e_2837+e_1394+e_2185-e_1425-e_1921-e_1841

also vanishes. Every nonempty proper balanced signed subblock of each
displayed source circuit has nonzero `sigma`. Thus supports four, four
and six are sharp for these declared binary before/after menus, already
at modulus eight. This is not a claim of minimum support among every
possible move using arbitrary other ports or sources.

The two-target trade `M=w_i v^t` therefore passes (7) for every dyadic
modulus. Its native paths, retained endpoints and arcs, chosen per-slot
phases, one-port capacities and strict injectivity are authenticated in
[the search packet](../search/README.md). Those data are not inferred from
the charge identity. In particular, the available outside inventory must
be respected if the packet is embedded in a larger assignment.

### Full image certificates, independent of code completeness

`image_witnesses.json` saves full grounded vectors `z1,zv,z2` satisfying

    L z1 = 2w1 + 256j1,
    L zv = 2v  + 256jv,
    L z2 = 2w2 + 256j2,                               (10)

with the complete integer correction vectors saved sparsely. The checker
reconstructs all 3,528 integer residuals from the entire original edge
file. It also checks the derived witness `z6=z1-z2 mod256` for `2w6`.
Thus the selected trades' modulo-256 image claims do not depend on the
generator, a claimed Smith basis or our code-completeness theorem. The
known exponent in (2) upgrades them to every dyadic modulus; it does not
turn these congruence solutions into integral solutions of `Lz=2w`.

## 6. Universal per-slot guards still forbid these retargetings

There is no contradiction with Session44's target-injectivity result.
That result tests the unmultiplied scalar phase in EACH actual port for
EVERY original character. Here we preserve the chosen `X` phase in each
slot, but only the OWNED AGGREGATE currents universally.

An exact separator makes the difference visible. Every original incident
conductance at target 1965 is even; its weighted degree is 926. Hence
`h=4e_1965` is a genuine grounded modulo-eight harmonic character.
It distinguishes targets 1965 and 2234 by four. A single retargeted port
therefore changes its `h` phase, while the paired-black contribution
multiplies that difference by two and vanishes modulo eight. The remaining
target-coordinate debt is paid by the source circuit, as (7) requires.

The checker saves the exact current `Lh/8`. The positive packet is thus
**chosen-character per-slot, universal-dyadic aggregate**, never universal
per-slot. Literal target ownership has not been identified or discarded.

## Reproduction and boundary

```text
python -B 05-knowledge/results/injective_native_blocks_20261023/electrical/check.py
python -O -B 05-knowledge/results/injective_native_blocks_20261023/electrical/check.py
```

Optional image regeneration:

```text
python -B 05-knowledge/results/injective_native_blocks_20261023/electrical/solve_images.py
```

The generator adapts Session44's exact bit-plane elimination. The checker
validates the result by direct edge multiplication. Original text inputs
are pinned with LF-normalized SHA256; the compressed edge export uses its
raw hash. The checker also exhausts 486 balanced three-by-three integer
trade shadows on six small composite-ring carriers, testing the full
image criterion against literal harmonic banks.

The decisive next obligation is native, not another scalar congruence:
find enough compatible distinct actual ports to compose these circuits
inside a retained global inventory. Electrical completion and literal
resource completion are different parts of that task.
