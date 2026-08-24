# Planar JC: rational moving-root decic, mixed Kummer support, and anchored infinity response

**Session status (2026-08-23):** mathematical research synthesis.  This
session promotes THM-3918 and THM-3919.  Subsequent incoming work promotes
THM-3920--3924: affine address and radial-chart closure; the audited six-
address degeneration; the affine-plane boundary class-basis theorem; the
fixed four-ray obstruction; and the index-five ramification-class
obstruction.  None settles the planar Jacobian conjecture.  `JC(2)` and
`DC(2)` remain **OPEN**.

## 1. Inheritance pass

The assigned anchor was the planar Jacobian conjecture, with the current
degree-three finite-normalization programme as its closest live mechanism.
The inheritance packet was:

- **closest proved mechanism:** THM-3913, whose moving-triple-root binary
  cubic simultaneously gives a normal globally nonmonogenic `S3` order, an
  irreducible one-place decic, and a nonzero quadratic-resolvent three-class;
- **canonical hostile example:** the same THM-3913 decic has elliptic
  normalization, so its discriminant branch is not polynomially uniruled and
  its cubic completion has no plane atlas;
- **corrected near miss:** THM-3914 identifies exactly when a pure degree-ten
  boundary residue globalizes under a prime-to-three finite-remainder
  hypothesis, but does not say that the actual THM-3913 class is pure
  boundary;
- **least-used sidecar:** the complete removed-divisor lattice, including
  infinitely-near origin components and local Kummer inertia rather than only
  the visible split-boundary Gram matrix.

The live concept board had five entries:

1. rational degeneration of the THM-3913 elliptic incidence curve;
2. exceptional-lattice gluing and the actual Kummer residue vector;
3. equality colors and first nonzero response jets;
4. nonproperness, polynomial uniruledness, and collapsed valuations;
5. an honest pairwise response carrier inherited from tournament work.

The Anchor / Niche / Wildcard portfolio was therefore:

- **Anchor:** rationalize the moving-root decic without losing normality,
  nonmonogenicity, one-place confluence, or the resolvent three-class;
- **Niche:** compute the full THM-3913 removed lattice and decide whether its
  Kummer class is actually pure boundary;
- **Wildcard:** determine whether the elliptic-to-rational transition is
  visible in a tournament, or whether a richer discrete carrier is forced.

All three lanes produced exact structure.

## 2. Incoming work as signal rather than conclusion

Three concurrent arrivals changed the useful comparison class.

### 2.1. THM-3915: a different rational decic

THM-3915 constructs a rational one-place decic with complete singular packet

```text
ordinary eightfold + two A2 + two A5,
```

an explicit rational normal nonmonogenic cubic maximal order, a mixed
three-class, and maximal-etale Euler characteristic `14`.  This proves that
rationality, normality, nonmonogenicity, one-place confluence, and global
three-torsion can coexist.  It does not provide a plane atlas.

The useful signal was not to copy its formulas.  It was to separate five
invoices that had looked entangled in THM-3913:

```text
branch rationality;
normal cubic order;
existence of a global three-class;
its support in the full removed lattice;
plane/Keller realizability.                               (1)
```

THM-3918 pays the first three by a sparse deformation inside the THM-3913
family, with a different singular packet but no rational cubic-field chart.
THM-3919 pays the fourth only for the original elliptic fibre; the rational
fibre's full removed lattice remains open.

### 2.2. THM-3916: rationality is not the Keller invoice

THM-3916 proves a general obstruction: a divisorial valuation of
`k(x,y)` with positive-genus residue field cannot be positive on both entries
of a scalar-Jacobian polynomial pair.  In the THM-3915 rational chart, the
common-zero divisor is genus two, so that particular completion contains no
Keller plane open.

THM-3916 alone does **not** automatically close THM-3918.  The rational curve
in THM-3918 is the discriminant normalization.  The curve in THM-3916 is a
divisor collapsed by both target coordinates in a chosen rational model.
They are different typed objects.  A successor must either construct a
rational chart for the THM-3918 cubic field and inventory its common-zero
divisors, or find a different chart-independent obstruction.  THM-3920 later
supplies that different obstruction through affine normalization addresses,
without proving cubic-field rationality.

### 2.3. THM-3917: genus zero escapes one gate and meets the next

During the audit, THM-3917 was proved and independently promoted.  Its
quintic-parameter depressed cubic makes the collapsed divisor rational, so
it genuinely escapes the positive-genus obstruction of THM-3916.  The
success is not a Keller chart: a rational ramification component has six
distinct normalization addresses forced through one point of the finite
normalization.  An irreducible boundary curve of an affine-plane open in a
normal surface is unibranch, so the six-branch collision closes that model.

This incoming result first sharpened the THM-3918 successor test into a
two-stage address audit:

```text
collapsed residue genus > 0  -> THM-3916 obstruction;
collapsed residue genus = 0  -> count boundary normalization branches.
```

Genus zero is therefore a necessary escape from one obstruction, not a
positive atlas certificate.  THM-3920 is now **PROVED + VERIFIED-EXACT +
INDEPENDENTLY HOSTILE-AUDITED**: boundary unibranchness and finite flatness
give the universal address cap `r<=d`; in degree three it excludes every
affine packet here with more than three addresses, including THM-3918.
Its character/Mason/third-Veronese/Kummer case split further closes the whole
irreducible radial chart `A=F/4,C=sF/4`.  THM-3921 is also proved and audited:
all algebraic invoices persist, but the same six addresses become four smooth
plus two cuspidal origin branches.  THM-3917/3920 close that model.

Two later arrivals expose the actual-completion invoice that the resolvent
could not see.  THM-3922 proves that the prime boundary divisors of any normal
affine completion containing `A2` form a primitive `Z`-basis of its actual
class group.  THM-3924 computes for THM-3921's cubic completion

```text
Cl(B)=Z<g>,                    [E_ram]=5g.
```

The mandatory ramification prime is nonprimitive, independently closing the
model.  The integer five also equals the `A^5` power-index exponent and the
ramification pole-valuation gap `6-1`.  This is an exact bridge in this order,
not a general index-to-class theorem and not a transfer from the quadratic
resolvent.

THM-3923 closes the fixed THM-3808/3853/3855 inverse-discriminant grammar:
four distinct tangent rays give four normalization addresses, but a cubic
Keller completion permits at most three.  Formal lifting remains correct as
algebra; it lands in the wrong local branch packet.  The live local designs
must collide the quartic rays to `3+1`, `2+2`, or `2+1+1`, make the quartic
jet vanish, or leave the common-zero packet.

## 3. Anchor result: the elliptic incidence has a rational boundary fibre

Start from the one-parameter binary cubic

```text
ell=AU+CV,
Phi_r=ell^3+C ell U^2+r C ell V^2+A^2V^3.                 (2)
```

The repeated-root incidence is birational to

```text
18y^2=(w+r)(w^2-r^2-6).                                  (3)
```

For generic `r`, this is genus one.  At `r^2=-6`, it becomes

```text
18y^2=w^2(w+r),                                           (4)
```

a nodal cubic with rational normalization.  Pulling the normalization back
through the incidence calculation gives the polynomial decic map

```text
L=162t^6-18r t^4+6t^2-r,
G=(18t^2-r)L,
A=tG,
C=-(3t^2-r/6)G.                                          (5)
```

The target coordinates have degrees `(9,10)`.  Their exact discriminant is
an absolutely irreducible decic with degree-ten row `-27A^10` and derivative
`-4` at its sole infinity point.  Thus its projective normalization is
`P1`, and its affine normalization is `A1`.

The sparse cubic order survives the specialization:

- the coefficients stay primitive but share the proper common-zero ideal
  `(A,C)`, forcing global nonmonogenicity;
- the `C=0` Kummer specialization proves generic cubic irreducibility;
- the irreducible discriminant has height-one exponent one, so the
  finite-free Delone--Faddeev order is maximal in codimension one and normal;
- the generic Galois group is `S3`;
- the split degree-ten boundary still forces constant units, and the
  connected cyclic cubic layer gives nonzero `Cl(Q)[3]`.

This is a real movement of the frontier: the exact elliptic obstruction of
THM-3913 is gone inside a close sparse deformation, while its order and
resolvent invoices remain.

## 4. The singular packet is complete, not inferred

The map `(5)` has eight reduced normalization addresses over the origin.
Its only two nonimmersive addresses satisfy

```text
t^2=-r/18,                                                (6)
```

and are ordinary `A2` cusps.  Equal radial slope pairs are controlled by the
involution

```text
iota(t)=-r/(18t).                                         (7)
```

Among the eight origin addresses, only the pair with `t^2=r/18` is paired by
`iota`.  Exact collision numerators vanish to orders six and seven in the
two target coordinates, giving intersection multiplicity seven.  All other
origin pairs are transverse.  Hence

```text
delta_origin=binom(8,2)+(7-1)=34.                         (8)
```

Together with the two cusp contributions,

```text
34+1+1=36=p_a(decic).                                     (9)
```

The arithmetic-genus ledger leaves no room for a hidden singularity.  This
is stronger than a successful parametrization or a finite singular-ideal
search: it identifies why completeness holds.

The comparison with THM-3915 is informative:

```text
THM-3915: 28 + 1 + 1 + 3 + 3 = 36;
THM-3918: 34 + 1 + 1         = 36.                        (10)
```

The same genus budget is paid by different contact allocations.  Rational
decic is therefore not a single singularity type; the contact ledger is a
state coordinate.

## 5. Niche result: the THM-3913 Kummer class is maximally mixed

The actual THM-3913 removed lattice is

```text
R=K_10 orthogonal_sum L_origin orthogonal_sum A2(-)^4.    (11)
```

The rank-eight origin block has determinant `12` and Smith form

```text
(1,1,1,1,1,1,2,6).                                       (12)
```

Its unique nonzero order-three residue has square `2/3`.  Each external
`A2(-)` block contributes a nonzero local residue of square `1/3`, and the
boundary `K_10` direction contributes a nonzero residue of integral square.
Thus

```text
A_R(3)=Z/9 direct-sum (Z/3)^5.                            (13)
```

Local triple-root charts prove that all four cusp coordinates of the actual
cyclic cover are nonzero.  The valuation-one Eisenstein perturbation at
infinity proves the boundary coordinate is nonzero.  Isotropy then forces
the origin coordinate:

```text
4*(1/3)+2/3=0 modulo Z.                                   (14)
```

Therefore the class uses every finite block as well as the boundary.  It is
not a pure boundary class and cannot be read from the `K_10` Gram matrix.
This is an exact realization of the repository's saturation warning:
pairwise data inside a sublattice do not determine its embedding in the
ambient Picard group.

THM-3914 is preserved, not retracted.  Its prime-to-three remainder
hypothesis fails because

```text
det(L_origin)*det(A2)^4=12*3^4.                           (15)
```

The apparent collision was a hypothesis mismatch, and the full lattice
identifies the missing coordinate.

## 6. Wildcard result: rationalization is a boundary tie, not a color collision

The radial squareclass of the family is

```text
h_r(z)=(r^2+6)z^4+2rz^2+1.                               (16)
```

With `q=z^2` and `eta^2=-6`, its binary homogenization factors as

```text
H_r(q,s)=(s+(r+eta)q)(s+(r-eta)q).                       (17)
```

The two color roots have constant discriminant `-24`; they never collide.
The genus drop happens because one color becomes the infinity factor `s`.
Introduce vectors

```text
I=(0,1),                 v_+=(r+eta,1),
v_-=(r-eta,1),           omega(p,q)=det(p,q).             (18)
```

Then

```text
omega(v_+,v_-)=2eta,
omega(I,v_+)omega(I,v_-)=r^2+6.                          (19)
```

The unanchored full gain is constant and cannot distinguish elliptic from
rational fibres.  At the critical parameter one anchored edge vanishes
simply, with derivative `-1`, and the next `q` coefficient is nonzero.  The
faithful carrier is therefore

```text
anchored valued graph + first tie response,               (20)
```

not a tournament.  Forcing a direction at the tie deletes the event being
measured.

This recovers the useful content of the older tournament programme without
forcing an irrelevant tournament object onto the geometry.  The native
binary observable is alternating, but the target predicate lives at its
zero and in its next response.

## 7. Natural-number encodings and their decoders

The user-supplied odd-square convention suggests treating an admissible
rank as the natural number itself rather than continually carrying its
embedding in a sparse set.  That is productive here, provided the decoder is
kept.  The following are honest natural coordinates:

| object | natural coordinate | decoder that must remain |
|---|---:|---|
| rationalizing anchor edge | `1` | which conjugate anchored gain vanishes and its derivative |
| lost radial `q` shell | `1` | `q=z^2`, so the visible `z`-degree drop is two |
| exceptional origin pair | `7` | the two normalization addresses and their target-coordinate jets |
| origin singularity | `34` | eight smooth branches, one contact-seven pair |
| total decic budget | `36` | degree ten and normalization genus zero |
| origin Kummer residue | `2` mod `3` | the origin block and its discriminant quadratic form |
| four cusp residues | `4` | four distinct local `A2` blocks, not an unlabeled count alone |

The rank is structural when the next operation pulls back to a law on it:
tie order adds under products, intersection length is a local algebra length,
and discriminant residues add in `Q/Z`.  It is only a scheduler when the
decoder is discarded.  In particular:

- the transitive ordinal order of response times does not retain their gaps;
- the number of normalization addresses does not retain contact depth;
- the boundary determinant does not retain saturation;
- the number of nonzero residue blocks does not retain their quadratic-form
  squares.

This is the planar-JC version of “regard the `n`th odd square as `n`”: strip
away a redundant ambient embedding, but do not strip away the operation law
or the coordinates needed to recover the predicate.

## 8. Typed connection contracts

| source | target/map | preserved predicate | destroyed information | required sidecar | cheapest decisive test |
|---|---|---|---|---|---|
| radial squareclass `(16)` | anchored gain graph `(18)` | elliptic versus rational fibre | color labels over `k`, higher target geometry | infinity anchor and first tie response | factor `(17)`, check discriminant `-24` |
| normalization map `(5)` | address/contact graph | complete delta ledger | individual branch equations after scalarization | conductor delta and intersection-length matrix | collision factorization plus genus budget |
| removed lattice `(11)` | discriminant residue packet | possible and actual three-class support | effectivity and full ambient geometry | local inertia and ambient saturation | Smith form plus local `uv=s^3` charts |
| affine branch normalization | finite-flat completion | address cap `r<=d` | contact orders and projective-infinity packets | affine target point and branch labels | THM-3920 fibre-length injection |
| four-ray tangent cone | cubic branch normalization | per-component address count | reducible global grouping | reduced irreducible component label | THM-3923 `4>3` gate |
| radial cubic field | same-field plane atlas via finite normalization | derivative divisors and common fibre | nonradial order structures | common radial factor and global `z` | THM-3920 character/Mason/Veronese/Kummer split |
| actual normal completion | affine-plane open via Weil localization | primitive free boundary basis in `Cl(X)` | effectivity and etaleness | actual completion, not resolvent | THM-3922 class-basis gate |
| `A^5` power index | ramification boundary class | integer five in this completion | a general index-to-class law | two A=0 primes and pole valuations | THM-3924 `[E]=5g` |
| natural response data | ordinal scheduler | chronological order only | gaps, ties, labels, coefficients | update law and decoder | compare two states with same rank order |

No row licenses a transfer beyond its preserved predicate.  In particular,
the class-group three-torsion does not imply a plane atlas, rational branch
normalization does not imply a rational collapsed divisor, and a valued
graph is not automatically a tournament.

## 9. Connections to earlier arithmetic and dynamics lanes

There is a useful but limited analogy with the earlier `abc`/IUT discussion.
The discriminant, index, boundary, cusp, and origin terms behave like a
local-to-global conductor budget: omitting one place can make a divisibility
permission look global when isotropy actually requires compensating local
coordinates.  The precise result here is geometric lattice arithmetic over
an algebraically closed field, not an `abc` inequality and not an
inter-universal Teichmuller statement.  No height comparison or arithmetic
field descent has been constructed, so the analogy remains a task generator.

The same restraint applies to LRC and Mahler response work.  Natural depth,
first tie response, and retained sidecars are common research operations,
but the objects and targets differ.  The transferable mechanism is:

```text
find the first vanishing native response,
retain its natural depth,
then carry the first nonzero transverse coefficient.      (21)
```

It is not a theorem transporting runner gaps or Mahler orbits into a
Jacobian proof.

## 10. Generated successor tasks

After THM-3920--3924, the strongest next tasks, in order of diagnostic value, are:

1. **Break both closed grammars.**  Construct a normal nonmonogenic binary-
   cubic order without the common radial factor/global `z` and without the
   fixed four-ray tangent packet.  Start with `3+1`, `2+2`, `2+1+1`, a
   vanishing quartic jet, or a unit-ideal packet.
2. **THM-3918 independent arithmetic.**  Decide cubic-field rationality and
   compute the full removed lattice.  Resolve the eight-branch origin with
   its contact-seven pair and the two external `A2` cusps; compute the full
   Smith/discriminant packet and actual cyclic-cover residue vector.  These
   questions no longer control Keller entry, because THM-3920 closes it.
3. **Euler tariff as a comparison invariant.**  Compute the compact-support
   Euler characteristic of the THM-3918 maximal etale locus and, if it is not
   one, the exact tariff a hypothetical plane open would have had to delete.
4. **Actual-completion class invoice.**  For the next order, compute
   `Cl(X)` itself—not merely a resolvent group—and test whether every mandatory
   ramification divisor extends to a primitive boundary basis.  Then combine
   units, Euler tariff, address cap, and nonproperness valuations.
5. **Compare the rational decics.**  Test whether THM-3915 and THM-3918 are
   related by `GL_2` source change, target automorphism, or a shared cubic
   resolvent family.  Use singular/contact packets as the first obstruction.
6. **Compare the THM-3917 escape.**  THM-3921 identifies its six addresses
   with the complete origin packet, proves the genus gap stays one, and
   retains the actual cyclic-layer class; THM-3924 also finds `[E]=5g` in the
   actual cubic completion.  Search outside both obstructions rather than for
   another radial genus-zero tuning.
7. **Response automaton.**  For general sparse tangent deformations, store
   the infinity anchor, conjugate norm gains, tie order, and first surviving
   shell.  Classify which states preserve branch genus and which merely
   detect degree loss.
8. **Positive-characteristic hostile probes.**  Identify exactly which steps
   fail in characteristics `2` and `3` and whether tame inertia is the only
   lost coordinate.  This is a boundary audit, not a route to characteristic
   zero by itself.

The highest-value immediate probe is Task 1.  THM-3920 shows that another
radial parameter search cannot succeed: the derivative-attachment grammar,
not the chosen coefficients, is the obstruction.  The next object must alter
the representation before another expensive plane-atlas computation.

## 11. Reproduction packet

```powershell
python -B 04-computation/jc2_rational_moving_triple_root_decic_thm3918.py
python -B -O 04-computation/jc2_rational_moving_triple_root_decic_thm3918.py

python -B 04-computation/jc2_thm3913_mixed_kummer_lattice_thm3919.py
python -B -O 04-computation/jc2_thm3913_mixed_kummer_lattice_thm3919.py

python -B 04-computation/jc2_affine_plane_boundary_unibranch_depressed_cubic_chart_thm3920.py
python -B -O 04-computation/jc2_affine_plane_boundary_unibranch_depressed_cubic_chart_thm3920.py

python -B 04-computation/jc2_quintic_decic_degeneration_order_thm3921.py
python -B -O 04-computation/jc2_quintic_decic_degeneration_order_thm3921.py

python -B 04-computation/jc2_binary_cubic_four_ray_tangent_address_thm3923.py
python -B -O 04-computation/jc2_binary_cubic_four_ray_tangent_address_thm3923.py

python -B 04-computation/jc2_decic_cubic_index_five_ramification_class_thm3924.py
python -B -O 04-computation/jc2_decic_cubic_index_five_ramification_class_thm3924.py
```

The companions have respectively `45`, `48`, `97`, `69`, `28`, and `20`
active gates.
All LF-normalized normal and optimized streams match their frozen LF outputs.  The lattice
theorem also survived an independent reconstruction of the log resolution, local cusp
inertia, infinity Eisenstein packet, and missing-divisor census.

The durable frontier statement is narrow but real: the moving-root family
now contains a rational one-place normal nonmonogenic `S3` decic with a
nonzero resolvent three-class, and the original elliptic fibre's three-class
is known to be maximally mixed.  THM-3920 closes the radial family, THM-3923
closes the fixed four-ray formal grammar, and THM-3922/3924 convert actual
completion classes into a primitive-boundary obstruction.  The live Keller-
design gap is a genuinely nonradial, at-most-three-ray or quartic-vanishing
order whose mandatory boundary classes form a primitive basis; `JC(2)`
remains open.
