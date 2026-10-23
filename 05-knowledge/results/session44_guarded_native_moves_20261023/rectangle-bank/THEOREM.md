# Full harmonic-bank obstruction for native black rectangles

23 October2026. **PROVED:** exact all-bank image criterion and the general
fixed-source-pair valuation obstruction. **FINITE-EXACT:** three original
order-eight charge orders, two native rectangle capacities, and their
universal modulo-eight exclusion. Complete original edge multiplication
verifies both the separating characters and explicit image witnesses.
No completeness assumption about a proposed harmonic bank or Smith basis
is used for these three order certificates.

## 1. The actual inherited rectangles and the changed question

Incoming commit `a201515762059411a745360467afff224f7db4c1`,
`native_guard_moves_farey_20261023/native_moves/NATIVE_GUARDED_MOVES.md`,
authenticates selected native observer rectangles on sources

    S=1238, R=2246

and target pairs `(881,2762)` and `(2650,2843)`. Source and target owners
are converse-merged original native tournament classes. The actual black
conductances on `(S,U),(S,V),(R,U),(R,V)` are respectively
`(4,2,4,4)` and `(4,2,4,6)`; all four blue entries are zero in each case.
A selected ordinary-source observer **plus its actual converse partner**
contributes black conductance2. Thus amplitude a has coefficient `c=2a`,
not a. No cancellation of this nonunit is allowed modulo8.

Those rectangles preserve one saved primitive character modulo8, with
higher paired dyadic limits16 and64. The present question is different:
do they preserve **every** original harmonic character modulo8? They do
not at any nonzero distinct-port amplitude allowed by those menus.

This is a typed extension of the existing annihilator/rowwise criterion in
`session44.../electrical-duality/UNIVERSAL_ELECTRICAL_GUARDS.md`, itself an
application of the original electrical image/kernel duality. The classical
module identity is not claimed as new. The new native increment is an
explicit exact order4 certificate and a capacity obstruction that was
invisible to the selected high-order character.

## 2. Exact all-bank criterion and amplitude

Let L be a symmetric integer Laplacian on a connected graph, reduced
modulo m. Let `H_m=ker(L mod m)`. Grounding one root gives the same tests
for zero-sum charges. The perfect finite character pairing gives

    (ker L)^perp = im(L mod m).                                  (1)

For four distinct vertices set `u=e_S-e_R`, `v=e_U-e_V`. The switch
removing coefficient c from edges `(S,U),(R,V)` and adding it to
`(S,V),(R,U)` has

    DeltaL = c(uv^t+vu^t).                                       (2)

Its action on a character x, in order `(S,R,U,V)`, is

    c(v.x, -v.x, u.x, -u.x).                                    (3)

**Theorem.** Every original character remains harmonic after the switch
iff both `cu` and `cv` belong to `im(L mod m)`.

Here preservation means `H_m(L) is contained in H_m(L+DeltaL)`: it does
not assert equality of the old and new kernels or critical quotients.

**Proof.** Disjoint supports and the unit coefficients of u and v make
`DeltaL x=0` equivalent to `c(v.x)=c(u.x)=0`. Quantify over the full kernel
and use (1) for each zero-sum charge. This proves both directions. QED.

In the finite charge quotient `C_m = D^0/(LZ^V+mD^0)`, let the orders of
[u] and[v] be `d_u,d_v`, and put `d=lcm(d_u,d_v)`. The exact scalar
condition is `d divides c`. For paired black coefficient `c=2a`,

    d/gcd(d,2) divides a.                                        (4)

This is only an **electrical** criterion. Actual source/target owners,
Hamiltonian paths, literal flags, converse pairing, nonnegative
conductances, distinct-port capacity and global matching remain separate
obligations. Zero Laplacian debt never supplies those missing objects.

## 3. Three exact charge orders, without a claimed complete bank

For the actual original3528-node Laplacian, all three differences

    d0=e_1238-e_2246;
    d1=e_881-e_2762;
    d2=e_2650-e_2843

have exact order **four** in `C_8`.

**Lower bounds from actual harmonic characters.** Let Y0,Y1 be entries0,1
of `radical_lifts_mod4` in the original electrical receipt. Literal
multiplication by the entire native Laplacian verifies `L(2Yi)=0 mod8`.
The differences of these actual mod8 characters are

| Character | d0 pairing | d1 pairing | d2 pairing |
|---|---:|---:|---:|
|2Y0|2|2|0|
|2Y1|0|6|2|

Thus `2d_i` pairs nontrivially with a harmonic character for every i,
so it is not in the image. These are genuine order-four separators.
The character2Y0 has1757 nonzero coordinates; it is **not** described as
sparse. The simple sparse characters `4e_j`, which are harmonic at all six
named vertices, have order two and are automatically blind to the even
black coefficient. Retaining the factor2 explains that failed shortcut.

**Upper bounds from explicit image witnesses.** `image_witnesses.json`
contains three full grounded vectors `z_i in {0,...,7}^3528`, each checked
by direct integer multiplication against every original edge:

    Lz_i = 4d_i + 8w_i.                                         (5)

The full integral correction `w_i` is also saved and independently checked.
Hence `4d_i` lies in the image. Combining (5) with the lower bounds proves
exact order4 without assuming that any saved harmonic vectors generate
the full bank. The witness generator used exact dyadic row elimination;
its unit-pivot counts were164 modulo8,3358 modulo4, and4 modulo2, with one
wholly unpivoted coordinate; lower-modulus pivot variables retain their
higher-bit freedoms. The **certificate checker does not trust those ranks or
the elimination algorithm**: equation(5) alone proves each upper bound.

Consequently both rectangles preserve the full original mod8 bank
precisely when their amplitude is even. Their minimum positive electrical
amplitude is two, despite passing the selected character at amplitude one.
For the explicit separator2Y0, the amplitude-one defects are respectively

    (4,-4,4,-4),     (0,0,4,-4),

both nonzero modulo8. In contrast the old selected character produces
`(-112,112,64,-64)` and `(128,-128,64,-64)`, both zero modulo8.
A character with a larger ambient order need not dominate other characters
of smaller order on a given charge.

## 4. Actual distinct-port capacity excludes the next electrical amplitude

The checker independently enumerates all409 actual Hamiltonian paths of
source1238 from its canonical matrix and the original forward-bit owner
map. Its complete target fibres have sizes

    target881:2; target2762:1;
    target2650:2; target2843:1.

The scarce paths are exactly

    to2762: (7,1,0,5,2,6,4,3);
    to2843: (7,0,6,4,1,5,3,2).

Therefore either orientation of either selected rectangle has maximum
absolute amplitude one when its observers must be distinct: amplitude two
would need two distinct actual source1238 ports on a side with only one.
Aggregate black-edge capacity does not create a second path. The incoming
paired/converse type and actual endpoint/arc restrictions are retained;
further restrictions can only shrink these already complete target fibres.

Combining this capacity with the exact electrical condition `2 divides a`
leaves only `a=0` for each incoming menu. Since the target pairs are disjoint,
the two-rectangle combination also has only the zero net amplitude under
these bounds. This obstruction is already in the **universal mod8 bank**,
not merely in the one-character modulus256 layer considered previously.

## 5. Many targets: a source-character gap imposes coordinatewise divisibility

Fix distinct sources S,R and `u=e_S-e_R`. Let v be any integer zero-sum
vector supported away from S,R; its entries are **net target amplitudes**
after adding all signed rectangle contributions. For any k,b>=0, consider

    DeltaL = 2^b(uv^t+vu^t) modulo2^k.                           (6)

Suppose a genuine original harmonic character x has
`u.x=2^s h mod2^k`, with h odd and `0<=s<k`. At every target coordinate j,

    (DeltaL x)_j = 2^b v_j (u.x).

Therefore a necessary condition for preserving even this one character is

    2^max(0,k-b-s) divides every v_j.                            (7)

Disjoint source/target support is load-bearing. For example modulo8, if
`u=v=(1,-1)` and `x=(0,4)`, then `2u(u.x)=0` although the coefficients
of v are odd. The source and target contributions coincide and cancel;
this falls outside (6)'s disjoint-support hypothesis.

This is exact for these target-coordinate equations; it is **not** claimed
sufficient for the source equations. The proof is ordinary integer
2-adic divisibility, retaining the nonunit coefficient2^b.

Here the authenticated character2Y0 has source gap2, so `k=3,b=1,s=1`.
Every net target amplitude must consequently be even, whatever target
pairs were chosen and whatever cancellations occurred between rectangles.
If the actual residual menu supplies `|v_j|<=1` at every target, then v=0
and `DeltaL=0`. Thus no incidence-changing universally mod8-safe trade
exists on this fixed source pair under those unit net capacities.

This does not force every named rectangle count to vanish: different
rectangles may cancel to v=0. Nontrivial reassignments of actual
presentations with unchanged incidence are allowed by this necessary
condition. Conversely, larger capacities can permit even nonzero v, and
additional source pairs can create new cancellations. Those cases are not
excluded by (7).

For completeness there is a precise all-bank criterion for general v.
Let `g=gcd_j |v_j|`, with g=0 when v=0. Because the source rows are
`+/-2^b v^t` and the target rows are `2^b v_j u^t`, the criterion is

    2^b v in im L, and 2^b g u in im L.                          (8)

The second condition is equivalent to all the target-row conditions by
integer Bezout, and the first is the source-row obligation. For primitive
rectangle v, g=1 and (8) reduces to section2. For v=0 both requirements
hold. This identifies exactly what the divisibility condition alone loses.

An original-atlas **algebraic hostile** makes this limitation visible.
Let Z be the saved primitive mod8 lift, whose difference `Z0-Z5=7 mod8`,
and set `v=2(e0-e5)` with the same source pair. All net target amplitudes
are even, but the source defect is `2(v.Z)=4 mod8`. This is a signed
matrix/character counterexample, not an assertion that those target
choices form a legal native rectangle.

## 6. Scope, provenance and reproduction

The source pin file records original atlas nodes, owner map, full438216-row
edge export, original electrical receipt, and incoming native rectangle
receipt at commit `a201515762059411a745360467afff224f7db4c1`. The earlier
coupled-moves and holonomy packets remain unchanged.

    python3 -B check.py --repo /path/to/math --out receipt.json
    python3 -O -B check.py --repo /path/to/math --out receipt.optimized.json
    cmp receipt.json receipt.optimized.json

The optional witness regeneration is deterministic and uses no external
linear-algebra package:

    python3 -B solve_images.py --repo /path/to/math --out image_witnesses.json

The fast checker validates the produced vectors directly, rechecks actual
characters and native capacities, tests the full-bank criterion on explicit
small composite rings, and exhausts the stated valuation controls. No
floating calculation, optimization-sensitive assertion or inferred
harmonic-bank completeness is used.

These results concern selected native packet capacities and necessary or
exact electrical tests on their incidence shadows. They neither rewrite
the canonical all-path atlas nor prove that every electrically admissible
candidate is an actual native wiring. An extension to a pre-used global
inventory still needs a compatible injective assignment of all affected
literal ports, including the unchanged remainder.
