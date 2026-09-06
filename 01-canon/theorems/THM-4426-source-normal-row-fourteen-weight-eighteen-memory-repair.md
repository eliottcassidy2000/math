---
id: THM-4426
title: "Source-normal row-fourteen weight-eighteen memory repair"
status: >
  PROVED FINITE-ROW REPAIR RELATIVE TO THM-4410/4415 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. Restoring the complete omitted weight-eighteen pair
  z=[p^9]H and h=[p^6*y^2]H turns THM-4415's unpaid row-fourteen bracket
  class into a rationally split conic with a Q-defined section over every
  characteristic-zero base point. On the exact boundary Phi=eta=0,
  alpha11=1, the simultaneous bracket-and-depth locus is a Q-rational G_m;
  every one of its points has terminal fibre A^10. The earlier two quadratic
  conjugate points are the h=0 slice of this rational curve. This is a partial
  source result, not full B2 membership, termination, chart entry, a Keller
  pair, JC(2), or DC(2).
source: root + row14_incoming_bridge + row14_memory_referee / JC2 continuation session, 2026-09-05
depends_on:
  - THM-4410-source-normal-least-weight-twenty-row-thirteen-affine-continuation
  - THM-4415-source-normal-row-fourteen-same-row-response-rank-obstruction
  - THM-4315-source-normal-student-stein-row-nine-high-contact-collapse
related:
  - THM-4403-source-normal-two-channel-weight-eighteen-row-twelve-affine-continuation
  - THM-4368-diagonal-boundary-valuation-triangular-address-and-simplex-stream-rank
script: 04-computation/jc2_source_normal_row14_weight18_memory_independent_referee_thm4426.py
output: 05-knowledge/results/jc2_source_normal_row14_weight18_memory_independent_referee_thm4426.out
script_sha256: 99071f872ccd1d3a25541599cf9e866b439d7919b1c7a2b7cc9b18e3ccd662ed
output_sha256: 8048778247f23d56b70ce8112576305c071cd8521bb77cb970adebcb1717af11
companion_script: 04-computation/jc2_source_normal_row14_weight18_two_memory_rational_conic_companion_thm4426.py
companion_output: 05-knowledge/results/jc2_source_normal_row14_weight18_two_memory_rational_conic_companion_thm4426.out
companion_script_sha256: c647bdf2c894d60413b041197cff446884191bf1469372e6b4a6a55ae96fd52a
companion_output_sha256: c9ac320b1adcdc16bfb282478188dbdbceabf343f04820cf9d633724bea5b0e7
polarization_script: 04-computation/jc2_source_normal_row14_weight18_two_memory_axis_diagonal_independent_audit_thm4426.py
polarization_output: 05-knowledge/results/jc2_source_normal_row14_weight18_two_memory_axis_diagonal_independent_audit_thm4426.out
polarization_script_sha256: 07c8be7da9da295005b5d209433c04f2b1b53c93541b686362186f2a1bb1b616
polarization_output_sha256: af99b4f5ee0c69875542e340e75323219fcddbf36960c429424bfd1718f1f4bb
hash_basis: raw LF bytes
audit: >
  PASS. A clean-room certificate imports only the audited THM-4308/4315 row
  operators, reconstructs rows nine through fourteen without importing the
  discovery scout, and performs 287 exact characteristic-zero checks. It
  verifies every prefix rank and graph, the literal THM-4415 J14 hash at
  z=0, the global quadratic coefficients, all fifteen boundary bracket
  coefficients, the complete 78-equation projected-depth remainder modulo
  Q(z), the nonzero depth pivot, and a rational hostile. Normal, optimized,
  and fixed-hash-seed runs byte-match the frozen LF output. A second
  clean-room companion restores the other weight-eighteen coordinate and
  reconstructs rows four through fourteen. Its 239 exact checks verify the
  global split-conic normal form and section, degeneration divisor, boundary
  G_m parameterization, all fifteen bracket and 78 depth coefficients, and
  the terminal pivot. Normal, optimized, and fixed-hash-seed runs byte-match
  its frozen LF output. An independently built axis/diagonal polarization
  audit checks the two mixed coefficients and all slice hashes in ten further
  gates, again with three byte-matching execution modes.
---

# THM-4426 -- Source-normal row-fourteen weight-eighteen memory repair

**PROVED FINITE-ROW REPAIR RELATIVE TO THM-4410/4415 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED.** The complete omitted weight-eighteen pair repairs
the same-row bracket obstruction over the ground field, and an exact rational
one-parameter boundary family also survives projected depth. The source is
still partial and the calculation stops at row fourteen. In particular,
`JC(2)` and `DC(2)` remain **OPEN**.

## 1. The restored memory coordinate

Work in the characteristic-zero source-normal chart of THM-4410. Its selected
tail through row thirteen retained

```text
lambda*p^3*y^4 + nu*y^6 + kappa*p*y^6.                 (1)
```

The complete residual-weight-eighteen face also contains `p^9` and
`p^6*y^2`. Write the two omitted coefficients as

```text
z=[p^9]H,       h=[p^6*y^2]H.                          (2)
```

At row fourteen, retain the paid same-row channel

```text
rho*p^2*y^6.                                           (3)
```

The resulting partial source is

```text
H_14 + z*p^9 + h*p^6*y^2 + lambda*p^3*y^4 + nu*y^6
     + kappa*p*y^6 + rho*p^2*y^6.                      (4)
```

Rows nine through thirteen have the same bracket and depth rank profile as
the inherited continuation. Their source projection is now the global
affine six-space

```text
(Phi,eta,alpha11,c51,z,h),                              (5)
```

and the row-thirteen terminal fibre remains `A^9`. Every solved graph has a
constant denominator. Setting `h=0` recovers the one-coordinate slice; the
complete two-coordinate consumer is computed after that slice below.

## 2. The unpaid class is a universal quadratic

After the row-thirteen graphs are substituted, the row-fourteen bracket
matrix has shape `15 x 9`, rank nine, and two active source conditions. The
first condition still has THM-4415's nonzero constant `rho` pivot

```text
1536794921806853409824256000000000.                    (6)
```

The second, in the primitive normalization agreeing literally with
THM-4415 at `z=0`, is

```text
J_18(z)=J_14 + L*z + A*z^2,                            (7)

A=643145476450643616480975000000,                      (8)

L=
-211041194570570137414961250000*Phi^2
+1125125025346225769232150000000*Phi*alpha11
+127589692647756992119425000000*Phi*c51
+373391507810486191201875000000*Phi*eta
+127589692647756992119425000000*alpha11*eta
+562562512673112884616075000000*eta^2
+63406745144628813691929829048320000.                  (9)
```

Here `J_14` is exactly the frozen unpaid quartic in THM-4415; the verifier
checks its expression hash, not merely its values at sample points. Since
`A` is a nonzero integer, `(7)` has a root over an algebraic closure at every
base point `(Phi,eta,alpha11,c51)`. Equation `(6)` then solves the other
condition. Therefore the obstruction proved in THM-4415 is a same-row
obstruction, not an obstruction after the lower-row memory `(2)` is restored.

Over a general characteristic-zero ground field, a root of `(7)` need not
belong to that field. The preceding solvability statement is over its
algebraic closure, or over an extension containing a chosen root.

## 3. The complete face has a ground-field bracket section

Restoring `h` changes the same primitive unpaid class to

```text
J_18(z,h)=J_14+L_z*z+L_h*h
            +A_zz*z^2+A_zh*z*h+A_hh*h^2,               (9a)

A_zz=643145476450643616480975000000,
A_zh=40291481888765365932450000000,
A_hh=626898609743995317881250000.                       (9b)
```

Here `L_z=L` from `(9)`, while

```text
L_h=
-7786670517605627944380000000*Phi^2
+35276293010813403389400000000*Phi*alpha11
+3970357861711970346581250000*Phi*c51
+12089412187921687824843750000*Phi*eta
+3970357861711970346581250000*alpha11*eta
+17638146505406701694700000000*eta^2
+2093570097805563863475514245120000.                   (9c)
```

The key new fact is not merely that this is quadratic. Its homogeneous
part splits over `Q`. Put

```text
u=801h+27826z,
v=46174707h+1363633922z,
k=16949646393750000.                                    (9d)
```

The coordinate determinant is `-192586625460`, and exact row reduction gives

```text
J_18=k*u*v+alpha*u+beta*v+c_0,                          (9e)
```

where `alpha`, `beta`, and `c_0` are polynomials over `Q` in
`(Phi,eta,alpha11,c51)` with only constant denominators. With

```text
U=u+beta/k,       V=v+alpha/k,
gamma=c_0-alpha*beta/k,                                 (9f)
```

the conic is exactly `kUV+gamma=0`. Therefore, for an auxiliary
`tau!=0`,

```text
U=tau,       V=-gamma/(k*tau)                           (9g)
```

is a rational parameterization. At `tau=1`, the recovered `h`, `z`, and the
paid `rho` graph have only rational constant denominators. Thus every point
of the four-dimensional characteristic-zero base has a ground-field
row-fourteen bracket solution; no quadratic extension is needed after the
complete weight-eighteen face is restored.

The exact primitive `gamma` is an irreducible quadratic over `Q` and is
frozen in the companion output. Off `gamma=0`, `(9g)` covers the full affine
conic. On `gamma=0`, the fibre is the union `U=0` and `V=0`; the displayed
chart covers only `V=0, U!=0` and misses the other component. This is the
precise degeneration boundary, not a failure of the `tau=1` section.

## 4. The one-coordinate boundary slice

First set `h=0` and impose the exact boundary

```text
Phi=eta=0,       alpha11=1,                             (10)
```

while keeping `(c51,z)` free. The unpaid bracket equation becomes

```text
D*c51+A0*z^2+B0*z+C0=0,                                (11)

D =45455611113114234086880000000,
A0=10049148069541306507515234375,
B0=990730392884825213936403578880000,
C0=10771463517407623861844320362236985344.             (12)
```

It solves `c51` with a constant denominator. After `(11)` and the paid
`rho` equation are imposed, the row-fourteen projected-depth matrices have

```text
P2:       150 x 353, rank 108,
P3:       165 x 564, rank 129,
joint:     78 x 15,  rank   5, pivots 10,11,12,13,14.  (13)
```

There is one remaining source equation, namely

```text
Q(z)=A0*z^2+B0*z+C1=0,                                 (14)

C1=10771463883409470380030782972892985344.             (15)
```

Its discriminant is

```text
Delta=
548570569425327391286894008779398603921084860366624411440000000000
     =(4141464600600000)^2
       *31983397604359613083640070627484254.            (16)
```

The squarefree factorization of the second factor is

```text
2*3*311374477*20767412909*824342968318013.             (17)
```

Thus `(14)` has two distinct real conjugate roots in

```text
K=Q(sqrt(31983397604359613083640070627484254)).         (18)
```

At both roots, the combined bracket and depth ideal fixes

```text
c51=1087/135.                                          (19)
```

The paid equation fixes `rho` in the same quadratic field. All fifteen
bracket coefficients and all 78 projected-depth coefficients vanish
coefficientwise modulo `Q(z)`. The depth pivot minor is exactly `1/32`, so it
does not vanish at either root; the terminal fibre over each point is exactly

```text
A^10_K.                                                (20)
```

This proves existence over the explicit field `(18)` and over every
characteristic-zero field containing it. It does **not** assert such points
on the line `h=0` over every characteristic-zero ground field. The full face
does have ground-field points, as follows.

## 5. The full boundary is a rational `G_m`

Keep `h` and retain the boundary `(10)`. The bracket equation remains linear
in `c51`. After its graph and the paid `rho` graph are imposed, the projected
depth rank profile is still `(13)` and the unique remaining condition is

```text
Q(z,h)=
 39181163108999707367578125*h^2/4
+629554404511958842694531250*h*z
+32712032778211935366804910080000*h
+10049148069541306507515234375*z^2
+990730392884825213936403578880000*z
+10771463883409470380030782972892985344=0.             (20a)
```

Use the same rational linear forms

```text
u=801h+27826z,       v=46174707h+1363633922z.           (20b)
```

Exact completion of the product turns `(20a)` into

```text
(u+U_0)(v+V_0)=R,                                      (20c)

U_0=85855050266495746048/37533020625,
V_0=5869475532385201397235712/262731144375,
R=4114203158584849143084737646705784371355648
  /394443738426269109375.                              (20d)
```

In particular `R!=0`, so the entire affine curve is isomorphic over `Q` to
`G_m`. For every nonzero parameter `s` in a characteristic-zero field, set

```text
u=s-U_0,       v=R/s-V_0,                              (20e)
```

and recover `(h,z)` through the invertible map `(20b)`. Conversely, every
affine point of `(20a)` arises from one such `s`; `s=0` and `s=infinity` are
the two projective points at infinity. Along this family the exact joint
ideal fixes

```text
c51=1087/135,                                          (20f)
```

and gives `rho` as a rational function with denominator proportional to
`s^2`. Substitution kills all fifteen bracket coefficients and all 78 depth
coefficients individually. The depth pivot is identically `1/32`, hence the
terminal fibre at every point of the family is `A^10`.

Choosing any rational `s!=0` gives a rational survivor. Thus the quadratic
extension in `(18)` is a property of the line section `h=0`, not an
irrationality obstruction on the complete solution component.

## 6. A hostile and the exact boundary of the claim

The rational choice `h=0,z=1` is a decisive hostile to the shortcut “bracket
solvability implies projected depth.” Both bracket equations can be solved
rationally there, but one raw depth coefficient remains

```text
-468367592341369707599392429738569553
 / 3074292539051204237760000000,                       (21)
```

so the projected-depth condition fails. The integer displayed in summaries
is the primitive numerator in `(21)`, not the literal rational residue.

The primitive value `Q(1,0)` in the companion normalization is
`10772454623851503274786025883987099719`; the two residues differ by the
nonzero factor `-1/70708728398177697468480000000` and represent the same
cokernel class.

The theorem proves only row-fourteen statements for the partial source `(4)`.
It does not compute the global row-fourteen depth cover away from boundary
`(10)`, establish full `B_2` membership, continue to row fifteen, prove
polynomial termination or chart entry, or construct a Keller pair.

## 7. Triangular addresses retain an intercept sidecar

For a source monomial `p^a*y^b`, define

```text
ell=2a+3b,     n0=a+2b,     N=a+b,
s=ceil(ell/2), u=n0-s+1,    v=N,
Addr=T(u+v-2)+u.                                        (22)
```

The weight-eighteen face has consecutive natural addresses

```text
p^9, p^6*y^2, p^3*y^4, y^6  ->  37,38,39,40.          (23)
```

The inherited tail kept only addresses 39 and 40. Address 37 repairs the
one-coordinate obstruction `(7)`; restoring adjacent address 38 completes
the weight-eighteen face and exposes the rational conics `(9e)` and `(20c)`.
This is a useful scheduling signal but not a lossless state code. The
weight-nineteen monomial `p^8*y` also has bare address 37 while having a
different intercept, row alignment, and parity.
Accordingly the consumer-complete label is at least

```text
(ell,Addr,coefficient),                                (24)
```

not `Addr` alone. This is an exact collision warning, not a theorem relating
the Jacobian predicate to a separate tournament or lonely-runner predicate.

## 8. Next exact consumers

The immediate consumers are now sharply ordered:

1. form the row-fifteen relative response problem over the exact boundary
   `G_m`, retaining `s` rather than specializing to one point;
2. compute the global projected-depth locus over the four-dimensional base
   conic `(9e)`, with special attention to the irreducible degeneration
   divisor `gamma=0` and its two components;
3. restore the next omitted weight/intercept packets only after checking
   which row-fifteen cokernel classes actually consume them.

## 9. Reproduction

Run from the repository root:

```powershell
python -B 04-computation/jc2_source_normal_row14_weight18_memory_independent_referee_thm4426.py
python -B -O 04-computation/jc2_source_normal_row14_weight18_memory_independent_referee_thm4426.py
python -B 04-computation/jc2_source_normal_row14_weight18_two_memory_rational_conic_companion_thm4426.py
python -B -O 04-computation/jc2_source_normal_row14_weight18_two_memory_rational_conic_companion_thm4426.py
python -B 04-computation/jc2_source_normal_row14_weight18_two_memory_axis_diagonal_independent_audit_thm4426.py
python -B -O 04-computation/jc2_source_normal_row14_weight18_two_memory_axis_diagonal_independent_audit_thm4426.py
```

The frozen LF output records the complete rank ledger, exact quadratic data,
field quantifier, hostile residue, and 287 optimization-live checks. The
companion adds 239 exact optimization-live checks of the two-memory conics,
their parameterizations, and every substituted bracket/depth coefficient;
the independent polarization audit adds ten.
