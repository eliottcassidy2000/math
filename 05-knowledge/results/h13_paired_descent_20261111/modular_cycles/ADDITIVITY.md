# Root amalgamation: actual additive switches and a larger exact order

**PROVED IN NOTE:** the full response of an ordered common-root amalgam;
an additive switch theorem with its actual baseline; preservation of the
positive source order by every right amalgam; a sufficient two-sided
order on seven coordinates; and a specific obstruction
to repaying the 211/213 reset by one converse switch. **FINITE-EXACT:**
27 named literal response controls, three four-corner constructions,
and twelve two-sided reset controls.
**OPEN:** selection of coupled switches that restores arbitrary target
counts. This extends [PROOF.md](PROOF.md), not the known spectrum range.

## 1. The operation and the coordinates it needs

Given rooted tournaments A and B, identify their roots r and orient every
edge from a nonroot vertex of A to a nonroot vertex of B. Call the result
`A ≺_r B`. Keep all original edges incident to r. This is an actual
tournament, and the operation is associative. A rooted singleton is its
identity. There is no choice of an unknown count at this step.

For a rooted tournament T use the following seven quantities:

|Coordinate|Meaning|
|---|---|
|H|Hamilton paths of T|
|E, S|Hamilton paths ending, respectively starting, at r|
|U|Sum of `H(T−vertices(P))` over all simple paths P starting at r|
|V|The same sum over paths ending at r|
|d|`H(T−r)`|
|K|`Σ_X H((T−r)[X]) H((T−r)[Xᶜ])`|

Empty paths are allowed in the two factors defining K, and H(empty)=1.
The paths defining U and V include the singleton r. Thus K is an ordered
two-path-cover count of the deletion tournament. The ordinary root state
is `(H,E,U,d)`; V and S retain the opposite-end information required here.
All seven quantities are determined by the literal graph and root.

Write lower-case a-coordinates for A and capital B-coordinates for B.
To avoid confusing H with the object, the main formula is written with
subscripts:

    H(A ≺_r B)
      = d_B(H_A−E_A) + d_A(H_B−S_B) + E_A S_B
        + (U_A−d_A)(V_B−d_B).                       (1)

When both factors have a nonroot vertex, the other coordinates are

    E = d_A E_B,              S = d_B S_A,
    U = (U_A−d_A)K_B + d_A U_B,
    V = (V_B−d_B)K_A + d_B V_A,
    d = d_A d_B,              K = K_A K_B.          (2)

The singleton is handled directly as the identity; in particular the
endpoint formulas must not be applied across an empty nonroot part.

**Proof of (1).** Split a Hamilton path immediately before and after r.
On either side of r, its nonroot vertices from A precede those from B.
There are four disjoint cases, according to whether the pre-root B part
and the post-root A part are empty. If only the former is empty, the
contribution is `(H_A−E_A)d_B`. If only the latter is empty, it is
`d_A(H_B−S_B)`. If both are empty, it is `E_A S_B`. If neither is empty,
the root-start portion in A and the root-end portion in B give
`(U_A−d_A)(V_B−d_B)`. All remaining crossing edges point forward and
require no further condition. These cases also prove (1) with singleton
factors.

For U, if the root-start path contains an A vertex, its A contribution
is U_A−d_A, and the B vertices may be assigned to its continuation or
the complementary path in K_B ways. If it contains no A vertex, the
contribution is d_A U_B. This proves U; reversing all arcs proves V.
Deleting r gives the ordinal sum of the two deletion tournaments.
Its one-path and ordered two-path counts multiply. The endpoint
formulas follow from the last and first nonempty parts. QED.

## 2. A genuine independent-switch construction

Suppose `A_0,A_1` have equal E,U,d and

    H(A_1)−H(A_0)=w_A.

For **every fixed B**, equations (1)-(2) give a pure-H switch

    state(A_1 ≺_r B)−state(A_0 ≺_r B)
        = (d_B w_A,0,0,0).                         (3)

This transports the source's real baseline and root, not merely the
abstract existence of a gap.

Now suppose `B_0,B_1` also have equal E,U,d, and have equal K. Set

    w_B = H(B_1)−H(B_0),
    σ = S(B_1)−S(B_0),    ν = V(B_1)−V(B_0).

The four actual graphs `A_i ≺_r B_j`, i,j∈{0,1}, have a common E,U,d
sidecar, and their counts are exactly

    N + i d_B w_A
      + j [d_A w_B + (E_A−d_A)σ + (U_A−d_A)ν],    (4)

where N is the value of (1) at `A_0,B_0`. In particular the mixed
second difference is zero. The extra K equality is what keeps U fixed;
it must not be discarded merely because the smaller root states agree.
If the A pair also has equal K, all four deletion K values agree too.

**Proof.** Subtract (1) in each direction. The coefficients of the
B switch involve only the fixed A sidecars, so no interaction term
survives. Formula (2) proves the common E,U,d assertion. QED.

If the B pair is instead the converse of a pure-H pair, its S,V,d agree.
The scalar rectangle simplifies to

    N + i d_B w_A + j d_A w_B.                     (5)

Its E and U need not agree across j. Thus (5) is an exact additive count
construction, but does not alone give a portable pure-H rectangle.

More generally, any finite ordered chain of switch pairs with equal
E,U,d,K within each pair gives an actual additive count cube. The common
root states have E,U,d,K independent of all choices, and H is the
baseline plus a sum of independent switch increments. This follows by
induction from (1)-(2). Explicitly, a switch at one position with prefix
P and suffix R has increment

    d_R [d_P w + (E_P−d_P)σ + (U_P−d_P)ν].         (6)

The prefix sidecars are choice-independent. For an empty prefix use
the singleton values d_P=E_P=U_P=1; an empty suffix has d_R=1.
Formula (6) exposes the actual transfer coefficient of every switch.
The increments need not all be positive in an arbitrary naming of the
two local alternatives.

## 3. A concrete pentagon rectangle

The H13 pair has core `(18,8,3,5,14)` and cuts6/3, with states

    (H,E,S,U,V,d,K) = (31,7,4,45,60,13,104),
                      (33,7,5,45,63,13,104).

Use on its right the cyclic pentagon core `(6,12,24,17,3)` and cuts3/5.
That pair has states

    (41,9,6,57,78,15,120),
    (43,9,6,57,78,15,120).

The right pair has the same forward **and** opposite-end sidecars.
The baseline supplied by (1) is

    15(31−7)+13(41−6)+7·6+(45−13)(78−15)=2873.

The four counts are therefore

|Left choice / right choice|41 side|43 side|
|---|---:|---:|
|31 side|2873|2899|
|33 side|2903|2929|

Their common forward sidecar is `(E,U,d)=(117,4581,195)` and their
common deletion two-cover count is12480. The two independent increments
are26 and30, and activating both adds56. These are literal tournaments
on11 vertices, not a formal sum of separately realized gaps.

The two H25 sources from [PROOF.md](PROOF.md) have equal deletion K608
at their211/213 realization, although their deletion graphs are different.
With H13 on the left, the same construction gives the pure-H rectangle

    [[16435,16359], [16485,16409]].                  (7)

Here the right increment is−76: σ=1 and ν=−3 give
`13·2+(7−13)·1+(45−13)(−3)=−76`.
This is a useful hostile control against assuming that a positive local
H switch remains positive in every right context. Its converse instead
gives the scalar rectangle `[[11001,11027],[11051,11077]]`, with the
expected independent increments26 and50, but differing endpoint sidecars.

The [construction continuation](../construction/PROOF.md) proves that the
H25 charts supply every positive even pure-H width, and that their
deletion K values agree for every fixed pair of flank addresses. Thus
the theorem applies to those wider switches as well. It does not identify
the resulting baseline with an arbitrary target.

## 4. Every right amalgam preserves the exact positive source order

The root-dynamics coordinates from [the order theorem](../dynamics/PROOF.md)
are

    q=(h,u,v,t)=(H,U,U−d,H+U−E).

Fix a nontrivial B. Formula (1) becomes

    h' = S_B h + (H_B−d_B)u
         + (S_B+V_B−H_B−d_B)v + (d_B−S_B)t.
    u' = U_B u + (K_B−U_B)v.
    v' = (U_B−d_B)u + (K_B−U_B+d_B)v.
    t' = h' + (U_B−E_B)u + (K_B−U_B+E_B)v.        (8)

Every coefficient in (8) is nonnegative. Here are the needed elementary
counting inequalities, including the potentially non-obvious one:

    S_B <= d_B <= H_B,
    d_B <= U_B <= K_B,
    U_B >= H_B >= E_B,
    V_B >= H_B+d_B−S_B.                            (9)

The endpoint bounds follow by deleting the root from a full endpoint
path. The insertion property of a path in a tournament gives H_B>=d_B:
each path of B−r admits an insertion of r, and deletion recovers it.
The root-alone path gives U_B>=d_B. Removing r from a path counted by
U_B injects it into the ordered two-path covers counted by K_B, giving
U_B<=K_B. Splitting any full Hamilton path immediately before r injects
it into U_B, giving U_B>=H_B.

Finally V_B contains the root-alone covers (d_B of them) and the covers
obtained by splitting a Hamilton path immediately after r (H_B of them).
The intersection consists exactly of Hamilton paths starting at r,
and has size S_B. Their union proves the last inequality in (9).

Consequently componentwise comparison of q is preserved by every
right amalgam. Together with the previously proved positive O/Z/J maps,
this establishes:

**Exact order theorem for the enlarged grammar.** For two actual rooted
sources, q dominance is equivalent to domination of H after every common
finite program using O, Z, J and arbitrary fixed right-root amalgams.

Sufficiency follows by multiplying the nonnegative response matrices.
Necessity follows from the O/Z witness words already in that grammar,
as proved in the cited root-dynamics note. This includes every finite
right gadget, not just the low-order controls in the accompanying script.

For the211 replacement, direct subtraction gives the exact cost

    H(old ≺_r B)−H(new ≺_r B)
      =90V_B−12H_B−78d_B.                          (10)

Since V_B>=max(H_B,d_B), this is nonnegative. It vanishes exactly when
V_B=H_B=d_B, equivalently the root of B is a global source. Indeed any
incoming edge x→r supplies an additional non-singleton contribution to
V_B beyond the d_B root-alone covers. Thus fixed-gadget amalgamation
preserves the reset's decrease; it does not silently repay its debt.

## 5. A sufficient order for arbitrary fixed contexts on both sides

There is a stronger sufficient comparison when the opposite-end and K
coordinates are retained. Suppose two sources have old minus new

    ΔH, ΔU, ΔV, ΔK >=0,
    ΔE, ΔS, Δd <=0.                               (11)

Both left and right fixed-root amalgamation preserve every sign in (11).
This is a **sufficient** cone, with no necessity claim.

For a fixed right factor B, the coefficients of E_A and d_A in (1) are
`S_B−d_B<=0` and `H_B+d_B−S_B−V_B<=0`; those of H_A and U_A are
nonnegative. For a fixed left factor A, the coefficients of S_B and d_B
are `E_A−d_A<=0` and `H_A+d_A−E_A−U_A<=0`, while the H_B and V_B
coefficients are nonnegative. The latter inequality follows by the
root-reversed version of the union argument in (9): `U_A>=H_A+d_A−E_A`.
Thus ΔH remains nonnegative in both placements. Formula (2) preserves
the other signs; the potentially negative coefficients are U−K and
V−K, both nonpositive by deleting the root from the marked path.

The actual reset's full states, in `(H,E,S,U,V,d,K)` order, are

    old211: (211,13,4,347,476,13,832),
    new211: (211,25,10,269,436,25,608),
    old213: (213,13,5,347,479,13,832),
    new213: (213,25,11,269,433,25,608).              (12)

Both satisfy (11), with V savings40/46 and K saving224. Therefore every
fixed ordered common-root context assembled on either side preserves the
decrease. In particular, moving the same gadget to the opposite side
does not repair the count debt. For a single left gadget P the lower
branch's debt is exactly

    52U_P−12H_P+6E_P−46d_P
      >=40(H_P−d_P)+46(d_P−E_P)>=0.                (13)

It vanishes exactly when P's root is a global sink: equality forces
H_P=d_P=E_P and U_P=d_P, and any outgoing root edge contributes to U_P
beyond its root-alone covers. The upper branch has52 replaced by58 and
46 replaced by52, with the same equality condition.

This order applies to fixed amalgamation contexts. It is not asserted
to survive every asymmetric birth or root change. The construction's
native asymmetric paired birth, which moves a difference between H and
E and can reverse a pair's order, remains a distinct possible escape.

## 6. One tempting compensation is impossible even with unbounded width

There is a stronger statement for a natural use of (5). Take **any**
pure-H pair B−,B+ of positive width w. Write its common coordinates as
E,U,d and its lower count as h. Compare

    old211 ≺_r converse(B−),
    new211 ≺_r converse(B+).

The difference of their actual counts is

    90U−12h−78d−25w.                               (14)

The common U is at least the upper Hamilton count h+w, while h>=d.
Therefore (14) is at least

    78(h−d)+65w > 0.                               (15)

No choice of that pair's width or baseline can make these two constructions
meet. The obstruction uses the same source whose gap is being used;
it is stronger than observing that an unrelated gap lacks a baseline.

More generally, suppose a same-H reset changes
`d_new−d_old=E_new−E_old=a>=0` and `U_old−U_new=b>=0`.
The corresponding difference after this converse-switch correction is

    (a+b)U−a h−b d−d_new w
      >= b(h−d)+(b−d_old)w.                        (16)

In particular b>d_old forbids this correction for every positive width.
For our reset b78>d_old13. This does not forbid a different left/right
placement, a chain whose genuine increments solve its actual debt,
or a construction changing the boundary coordinates before the join.
Equations (4)-(6) supply the appropriate coupled equations for those
next proposals; independent arithmetic gap addition is unnecessary.

## 7. Reproduction and scope

```bash
python3 -B 05-knowledge/results/h13_paired_descent_20261111/modular_cycles/amalgam.py
```

The script builds literal graphs, computes all seven coordinates by
subset path enumeration, checks (1)-(2) and the q response on27 named
ordered pairs, and checks all four corners of the C5, H25, and converse-H25
rectangles. Twelve additional literal controls place transitive endpoint
roots and a cyclic root on either side of both reset branches. The results,
source rows, exact roots, and full seven-coordinate
states are retained in `amalgam_receipt.json`. The universal results use
the path splitting proofs above, not an extrapolation from those checks.

The recovered modular common-quotient idea led first to paired compatible
birth cuts in [PROOF.md](PROOF.md), and now to explicit common-root
amalgamation. The connection preserves a shared boundary response while
allowing two independent internal lifts. Its required sidecar is now K
as well as both endpoint responses. The missing target-selection law
has become an equation among actual additive increments and an actual
baseline, rather than a claim that every even numerical gap can be used
wherever it is wanted.
