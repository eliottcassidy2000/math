---
id: THM-4431
title: "Colored lattice basis and three-direction LRC network closure"
status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED
source: overnight-hexagon-sep05 second research wave
---

# THM-4431 -- Colored lattice basis and three-direction LRC network closure

**PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The [full proof and independent raw-carrier/row/literal controls](../../05-knowledge/results/lrc14_colored_basis_three_ray_overnight_hexagon_sep05.md)
are part of this theorem.

## Colored extremal replacement

Let L be a rank-two lattice, H any subgroup, and B=-B a convex subset
of its real span. If (L minus H) intersect B spans the plane, it contains
a basis of L. No boundedness or closedness assumption on B is needed.
Choose an admissible pair of minimum positive lattice determinant.
A centered nonzero fundamental-parallelogram representative either is live
and reduces the determinant directly, or belongs to H and its signed
complement with a live endpoint stays live in B and reduces the determinant.
The forbidden subgroup is retained throughout; a dead short vector alone
is not a replacement.

For primitive sorted distinct positive odd ternary-unit speeds w=(a,b,c),
let Lambda be the **complete** THM-4414 raw-carrier dictionary. Its ambient
colored lattice has index three in ker_Z(w), and its dead subgroup has
index three in that lattice. Every spanning Lambda contains u,v with
u cross v=+/-3w.

## Exactly three directions and the network bound

If Lambda has exactly three primitive unoriented rays, suitable live-basis
coordinates make their directions exactly one of

```text
+/-{u,v,u+v},             +/-{u,v,v-2u}.
```

The circuit coefficients are (1,1,1) or (1,1,2), respectively; the second
type cannot be deleted as a cosmetic variant of the first. All raw multiples
remain in the physical network sum.

For exactly three rays, every projection satisfies E_i<6/77 at every
height. Each ray norm is at least seven, and pair products are at least
3c/2. They give the sufficient bound

```text
E_i < max((36/49)*sqrt(2/(3c))+12/(7c), 12/343+4/c).
```

This closes c>=99. The full head 1<=a<b<c<99 contains5,409 eligible
triples and1,791 three-ray triples (1,107 first type;684 second type);
independent raw formulas, Beatty rows and native interval controls close
all of them strictly. The head maximum18/301 at(5,37,43) is not asserted
to be the all-height sharp constant.

For a complete colored carrier set with no three distinct collinear points,
parity zero is absent (otherwise -C,C/2,C are three live collinear points).
Each of the two live colors uses at most one point in each other parity
class, so its size is at most six, sharply at(23,29,37).
The same argument gives2(2^d-1) for rank d and an index-three dead subgroup;
no higher-rank sharpness claim is made. This is a **cap** bound, not a bound
of six on arbitrary raw dictionaries or on multiple points of one ray.

The colored-basis/circuit mechanism was independently obtained in the
concurrent empty-core session; its contribution is credited in the full
proof and is not counted twice. The later independently audited
[universal network proof](../../05-knowledge/results/lrc14_global_slope_empty_core_certificate_sep06.md)
subsumes the three-ray target, not this structural classification. Its
THM-4434 namespace is still RESERVED at the integration check and is not
a proved dependency here.
Neither local network theorem proves chart entry, synchronization or LRC(14).
