# The D5 `H^1` map reaches the marked Kummer normal slice and stops before flux

**Status: historical synthesis with a FINITE-EXACT sidecar, independently
audited and promoted only in the repaired scope of
[THM-3496](../01-canon/theorems/THM-3496-marked-graph-kummer-degree-square-and-finite-coefficient-frobenius-flux-extinction.md);
this reflection is not a truth source and gives no LRC(14), JC(2), or
physical-current result.**  Reproduce the original finite claims and
polynomial identities with
`04-computation/d5_marked_kummer_normal_slice_frobenius_flux_no_go_20260816.py`.

## Independent audit correction -- 2026-08-16

The marked graph/deck map, normal Kummer square, additive no-go, and
Frobenius telescope survive.  Three scope repairs are load-bearing.

1. The equality
   `H^1_et(K((lambda));mu_13)=Z/13[kappa_lambda]` uses the stated
   algebraically closed characteristic-zero residue field.  It is false over
   a general characteristic-zero field: `[2]` is an extra valuation-zero
   class over `Q((lambda))`.
2. Degree and orientation naturality leave all twelve nonzero line
   isomorphisms.  The exponent-one marking selects one normalized map; there
   is no unmarked canonicity or maximality theorem.
3. Equations (43)--(44) prove that `[1]` becomes exact over `Z/13^r` for
   every finite `r`, and that ordinary `13`-adic completion kills it.  They do
   **not** prove a universal no-go for derived completion, `lim^1`, or every
   construction called a Bockstein tower.  All stronger Bockstein wording
   below is superseded by this finite-coefficient statement.

## 1. Inheritance pass and verdict

The closest general mechanism is
[THM-2622](../01-canon/theorems/THM-2622-affine-torsor-holonomy-fixed-section-spectrum-and-v4-c13-dictionary.md):
after a linear local system is fixed, an affine cyclic bundle is classified
by its return translation.  The constant-coefficient specialization in
[THM-3487](../01-canon/theorems/THM-3487-two-twenty-four-state-fibonacci-bundles-cycle-type-obstruction.md)
is the cellular rule

```text
edge translations / vertex gauges  =  H^1(circle;V),
class                               =  seam sum.             (1)
```

The canonical hostile is
[THM-3354](../01-canon/theorems/THM-3354-inequivalent-h1-carriers-and-typed-obstruction-cospan.md):
the actual LRC, sporadic-monodromy, quartic-Kummer, and Hamiltonian-response
objects have different sites and coefficients, and every direct additive
coefficient map from the odd LRC class to the characteristic-zero response is
zero.  The corrected near miss is
[THM-3450](../01-canon/theorems/THM-3450-marked-d5-carrier-isomorphism-and-full-germ-margin-obstruction.md):
it gives a strong marked characteristic-zero *carrier* isomorphism, but not
an `H^1` map.  The least-used relevant sidecar is the connecting morphism in
[THM-3406](../01-canon/theorems/THM-3406-affine-modification-power-jets-and-principal-part-transgression.md),
together with the full root-arm persistence of
[THM-3412](../01-canon/theorems/THM-3412-hamiltonian-principal-part-differential-and-prufer-torsion-arms.md).

The resulting verdict has a positive and a negative half.

1. After selecting one JC response divisor and its oriented normal
   meridian, there is an explicit marked homomorphism

   ```text
   Phi_lambda:H^1_graph(C7;F13)
       -> <kappa_lambda> subset H^1_et(K((lambda));mu_13).   (2)
   ```

   It sends LRC chart holonomy to the Kummer torsor
   `y^13=lambda`, respects graph/deck gauges, and intertwines cyclic cover
   degree with normal ramification degree.

2. Equation (2) does **not** extend to the additive Hamiltonian flux.  The
   direct characteristic-zero homomorphism is forced to be zero.  Pole order
   is not additive.  If one instead matches coefficients by reducing the JC
   complex modulo thirteen, the smallest one-root model

   ```text
   P=x+x^2 z                                                   (3)
   ```

   acquires a Frobenius first integral which kills the
   principal-part transgression exactly.

Thus the most explicit lawful D5 map currently available is a **marked
divisor-normal Kummer map**, not a word-current-to-Hamiltonian-flux map.

## 2. Which LRC object is actually an `H^1` class?

Several live LRC objects use current or cocycle language, but only one is the
source of (2).

| object | actual type | why it is not automatically the source of an `H^1` map |
|---|---|---|
| THM-2461 temporal blocker word | a directed prescribed-clock atom-to-word coupling | composition is temporal and root-typed; it is not an invertible local system |
| THM-2334 relation current `A_rho(r)` / `A(q)` | an absolutely convergent complex measure on relation-address orbits and its finite pushforward | addition is amplitude addition; no graph coboundary quotient is supplied |
| THM-2337 first jet `beta mod 13` | a factor-coloured divided-difference sidecar | the atomic gauge changes it while fixing the exact relation address |
| THM-2512 doubly-centred interaction | a rational `7 x 13` ANOVA/Fourier carrier | THM-3450 maps its primitive characteristic-zero line, not a finite `H^1` class |
| THM-2542 chart transition `g` | an `F13`-valued one-cocycle on the seven-chart nerve | this is the one actual graph-`H^1` source |

In particular,
[THM-2337](../01-canon/theorems/THM-2337-expiration-word-residue-invisibility-and-first-bockstein-sidecar.md)
has the exact gauge

```text
(u,beta,v,m) -> (u-13^k gamma,beta+gamma,v,m),          (4)
```

which fixes the full relation address but changes `beta mod 13`.  Therefore
a formula built from the first word jet is not well defined on the unmarked
relation current.  Likewise, THM-2334's current is a measure, not a class in
the graph group below.

The source of (2) is precisely the spatial chart cocycle constructed by
[THM-2542](../01-canon/theorems/THM-2542-seven-chart-cech-holonomy-and-c91-arrival-obstruction.md),
which inherits the word/current machinery but still lacks its physical
semantic-arrival two-cell.  No arrow from a full THM-2334 or THM-2512
physical current to this chart class is claimed here.

## 3. The source complex, gauge, and deck transgression

Orient the seven-chart nerve and put

```text
C_L^0=F13^(V(C7)),       C_L^1=F13^(E(C7)),
(delta f)_i=f_(i+1)-f_i.                                  (5)
```

There are no two-cells, so every edge cochain is a cocycle.  Vertex gauge is

```text
g -> g+delta f.                                           (6)
```

The seam functional

```text
s(g)=sum_(i in C7) g_i                                    (7)
```

kills every coboundary and induces an isomorphism

```text
H^1_graph(C7;F13) -> F13.                                 (8)
```

For THM-2542's constant chart step `g_a=(a,...,a)`,

```text
s(g_a)=7a !=0                 when a!=0.                   (9)
```

Now pull back to the degree-thirteen cover `C91 -> C7`.  The pullback class
is zero because its total sum is `13s(g)=0`.  If `h` is a primitive, the
deck generator `tau(j)=j+7` has defect

```text
chi_g(tau)=tau^*h-h=s(g).                                 (10)
```

Changing `h` by a constant does not change `chi_g`.  This is exactly the
isomorphism

```text
T_L:H^1_graph(C7;F13) -> H^1_group(C13;F13)               (11)
```

proved in
[THM-3431](../01-canon/theorems/THM-3431-d5-secondary-h1-descent-defects-and-valuation-persistence.md).
It identifies the source base graph, coefficient local system, cocycle, gauge
action, cover groupoid, and the distinguished class without treating the
length `91` carrier as a `Z/91` coefficient.

## 4. The smallest coefficient-compatible JC target

Use a selected one-root JC response chart as in THM-3422/3427 and THM-3431,

```text
P=ax+b+c(x-alpha)^e z^d,
lambda=(P-(a alpha+b))/a.                                (12)
```

After base change to an algebraically closed characteristic-zero field, take
the punctured formal normal slice of the response divisor `lambda=0`:

```text
U_lambda=Spec K((lambda)).                                (13)
```

The Kummer sequence gives

```text
H^1_et(U_lambda;mu_13)
  =K((lambda))^*/K((lambda))^(13)
  =Z/13 * kappa_lambda,                                  (14)
```

where

```text
kappa_lambda=[y^13=lambda].                              (15)
```

The equality in (14) is valuation modulo thirteen: because `K` is
algebraically closed of characteristic zero, every constant unit has a
thirteenth root and formal Hensel lifting gives a root of every principal
unit; `lambda` has valuation one.  Over a general residue field the unit
quotient need not vanish (for example `[2]` over `Q((lambda))`).
Consequently, under the stated field hypothesis, an orientation-preserving change of uniformizer
does not change `kappa_lambda`; reversing the meridian sends it to
`-kappa_lambda`.  Changing the sheet coordinate `y` by a deck element is a
torsor gauge and also leaves the class unchanged.

This is a genuine JC **divisor cocycle**.  It is not yet a Hamiltonian flux.
It remembers the oriented meridian and ramification exponent modulo thirteen,
but not the polynomial derivation, root multiplicity, number of collided
root arms, integral observer, or principal-part depth.

## 5. The explicit marked map

Evaluation of (11) at the chosen deck generator and Kummer realization give

```text
Phi_lambda
  = (t |-> t kappa_lambda) after ev_tau after T_L,         (16)

Phi_lambda([g])=s(g) kappa_lambda.                        (17)
```

Thus the canonical class maps as

```text
[g_a] -> 7a kappa_lambda.                                (18)
```

Every map in (16) is an isomorphism between one-dimensional `F13` groups.
Before the exponent-one normalization there are twelve nonzero scalar
choices, all compatible with degree and orientation.  Thus the markings
select the displayed map but prove no unmarked canonicity or maximality.
At the finite exponent-cochain level, choose one target meridian seam
`e_*`.  Then

```text
F^0=0,
F^1(c)=s(c)e_*,
F^1(delta f)=0.                                          (19)
```

Moving the target seam changes `F^1(c)` by a target coboundary.  Realizing
`e_*` as the monodromy exponent of (15) gives (17).  This is the concrete
finite slice checked by the companion.

The word **marked** is load-bearing.  The data are:

```text
source orientation and deck generator tau;
selected JC root-value divisor lambda=0;
target meridian orientation;
identification of exponent one with y^13=lambda.          (20)
```

Without (20), the nonzero identifications form an `F13^*` torsor.  Also,
(17) is not induced by a morphism from the seven-chart nerve to the JC
normal slice.  It is a correspondence through their marked cyclic quotient.
This is more than a shared dimension, because it carries the gauges and the
degree action below, but it is not a natural transformation between the
original geometric sites.

## 6. The nontrivial commuting square

Let

```text
p_k:C_(7k) -> C7                                        (21)
```

be the degree-`k` cyclic cover, and let

```text
r_k:Spec K((t)) -> Spec K((lambda)),       lambda=t^k   (22)
```

be the degree-`k` normal ramification.  If `kappa_t` is the target Kummer
generator, then

```text
s(p_k^*g)=k s(g),
r_k^*kappa_lambda=k kappa_t.                            (23)
```

Therefore the following marked square commutes:

```text
H^1(C7;F13)          --Phi_lambda-->  <kappa_lambda>
     | p_k^*                                | r_k^*
     v                                      v
H^1(C_(7k);F13)      --Phi_t------->  <kappa_t>.          (24)
```

Both paths send `[g]` to `k s(g)kappa_t`.  In particular,

```text
k=13: both classes die;
k=14: both return to their original value.              (25)
```

This upgrades THM-3431's common death-bar grammar to an actual marked map on
the cyclic normal-slice operation.  The boundary is equally important:
`r_k` is a normal ramification operation, not composition of Keller maps or
a polynomial mate.  Equation (24) transfers a divisor-meridian class and no
JC injectivity predicate.

## 7. Why the Kummer map cannot be promoted to flux

The characteristic-zero Hamiltonian object is the cohomology of the
two-term complex

```text
K_P=[R --D_P--> R],
H^1(K_P)=C_P=R/D_P(R).                                  (26)
```

For `P=f(x)+g(x)z`, THM-3406 also has

```text
S=R[g^-1],       Pi=S/R,
delta:H^0([Pi --Dbar--> Pi]) -> H^1(K_P),
delta([h])=[D(h)].                                      (27)
```

THM-3412 computes the whole filtered source of (27): one Prüfer arm for each
geometric root of `g`, with root multiplicity and collided-root copies still
visible.  THM-3431 separately embeds a selected cyclic observer as

```text
xi_q=[lambda^-q]
  in H^1_((lambda))(K[lambda]).                         (28)
```

There are three independent obstructions to continuing (17).

### 7.1 Additive coefficient order

The source of (17) has exponent thirteen, while `C_P` and (28) over
characteristic zero are vector spaces over a characteristic-zero field.
Hence

```text
Hom_Add(F13,C_P)=0,
Hom_Add(C_P,F13)=0.                                    (29)
```

The same argument applies to the additive local-cohomology module.  Kummer
works only because its coefficient sheaf is `mu_13`, not because the old
additive no-go disappeared.

### 7.2 Pole order is not additive

One might try to send the monomial principal part `xi_q` to
`q kappa_lambda`.  This is a valuation operation, not a homomorphism.  The
classes `xi_q` and `-xi_q` have the same pole order, so that rule sends both
to `q kappa_lambda`, whereas an additive homomorphism must send the second to
`-q kappa_lambda`.  For `q` nonzero modulo thirteen these disagree.  Also,

```text
q and q+13 have the same Kummer image                         (30)
```

but different annihilator depth and different finite stages of a THM-3412
Prüfer arm.  Multiple roots with the same `f`-value supply another collision:
one normal meridian, several persistent arms.

### 7.3 Matching coefficients by reduction kills transgression

The next section gives the minimal exact witness.  It is stronger than (29):
even when both sides are made `F13`-linear, the principal-part class survives
but its connecting flux dies.

## 8. The minimal Frobenius flux no-go

Work first over `Z` with

```text
P=x+x^2 z,
g=x^2,
D=(1+2xz) partial_z-x^2 partial_x.                     (31)
```

The Bezout row

```text
A=1-2xz,             C=4z^2,
A P_x+Cg=1                                             (32)
```

gives

```text
h=-A/g=-x^-2+2z/x,
D(h)=6z.                                                (33)
```

Thus `eta=[h]` is the THM-3406 principal-part cycle and its transgression is
`mu=[6z]`.  Multiplication by `P` gives

```text
P h=-x^-1+z+2xz^2,
[P eta]=[-x^-1] in Pi,
D(-x^-1)=-1.                                           (34)
```

Consequently

```text
delta(P eta)=[-1] in C_P.                              (35)
```

Over characteristic zero this class is nonzero.  There is a direct
one-line recurrence proof.  Give `x,z` weights `1,-1`; then `D` raises weight
by one.  A polynomial primitive of `1` could only use the weight-`-1`
monomials

```text
x^j z^(j+1).                                           (36)
```

The coefficient equations force `c_j=(-1)^j`.  A finite sum ending at `N`
then leaves terminal coefficient `(N+2)c_N`, which cannot vanish in
characteristic zero.

Now define

```text
Q_13=sum_(j=0)^11 (-1)^j x^j z^(j+1).                 (37)
```

Exact telescoping gives the integral identity

```text
D(Q_13)=1-13(xz)^12.                                  (38)
```

Therefore over `F13`,

```text
D(Q_13)=1.                                             (39)
```

The Hamiltonian flux class `[1]` vanishes after the coefficient change that
would be needed to receive the LRC `F13` class.  More structurally,

```text
Q_13-x^-1=-P^12/x^13             over F13,             (40)
```

because `1+(xz)^13=(1+xz)^13`.  The right side is a localized Frobenius
constant: `D(P)=0` and `D(x^-13)=0` in characteristic thirteen.  Its image in
`Pi` is `[-x^-1]`, so the principal-part cycle in (34) now comes from
`H^0(S,D)` and exactness forces its connecting class to be zero.  This is the
specific mechanism behind the characteristic-`p` warning in THM-3406.

The phenomenon is not a one-prime accident.  For every odd prime `p`, put

```text
Q_p=sum_(j=0)^(p-2) (-1)^j x^j z^(j+1).               (41)
```

Then

```text
D(Q_p)=1-p(xz)^(p-1),
Q_p-x^-1=-P^(p-1)/x^p                    over F_p.      (42)
```

The exact all-finite-coefficient statement is stronger than reduction modulo
thirteen alone.  For every `n>=1`, the analogous finite telescope gives

```text
[1]=(-1)^n(n+1)[(xz)^n]        in coker(D over Z).     (43)
```

Taking `n=13^r-1` yields

```text
[1]=13^r[(xz)^(13^r-1)]         for every r>=1.        (44)
```

The natural integral response class is infinitely thirteen-divisible.  More
precisely, reducing (44) modulo `13^r` makes `[1]` exact for every finite
`r`, while the characteristic-zero flux remains nonzero.  This proves loss
under every finite coefficient quotient and under ordinary completion; it
does not adjudicate derived completion, `lim^1`, or arbitrary higher
Bockstein packages.

## 9. Connection and loss ledger

| field | exact content |
|---|---|
| source site | oriented seven-chart nerve `C7`, followed by its degree-13 deck group |
| source coefficients/class | constant `F13`, `[g]`, `s(g)=chi_g(tau)` |
| target site | punctured formal normal slice of a selected JC response divisor |
| target coefficients/class | `mu_13`, `kappa_lambda=[y^13=lambda]` |
| map | `Phi_lambda([g])=s(g)kappa_lambda` |
| gauge actions | source vertex coboundary; cover-primitive constant; target Kummer deck gauge; simultaneous orientation reversal |
| preserved | zero/nonzero seam class, exponent mod 13, cyclic cover/normal-ramification action, degree-13 death and degree-14 return |
| destroyed | THM-2334 amplitudes, word jet, target mask, positivity, ancestry, owner/clock data, root multiplicity, collided-arm count, principal-part depth, additive flux, polynomial mate predicate |
| source sidecar still missing | a physical word/current-to-chart chain map or semantic two-cell |
| target sidecar still missing | a coefficient-compatible derived realization carrying Kummer meridian data to the integral Hamiltonian complex |
| cheapest flux hostile | `P=x+x^2z`, equations (37)--(40) |
| status | marked divisor-normal `H^1` map proved; flux extension blocked |

This is exactly the kind of boundary THM-2622 and THM-3487 warn about.  A
seam class is complete only after the linear coefficient system is fixed.
Here the Kummer and LRC systems can be synchronized, but the Hamiltonian
complex is a different realization functor.

## 10. Relation to the marked mode-one carrier

The present map and THM-3450 are complementary rather than composable:

```text
LRC chart H^1(F13) ----Phi_lambda----> divisor Kummer H^1(mu13)
       ^                                      |
       | physical semantic arrival OPEN      | flux realization BLOCKED
       |                                      v
THM-2334/2512 current ----THM-3450----> marked char-zero mode-one carrier.
                                                               (45)
```

THM-3450 preserves the rational doubly-centred interaction and a marked
characteristic-zero amplitude.  It necessarily discards the mod-thirteen
`H^1` class.  Equation (17) preserves the mod-thirteen seam and cyclic degree
action.  It necessarily discards the amplitude and all Hamiltonian response
depth.  Neither route supplies the missing vertical arrows in (45).

## 11. The sharpened D5 frontier

A future D5 theorem should no longer ask for a bare homomorphism between two
abelian `H^1` groups.  The exact residual is a derived, filtered
correspondence with four obligations.

1. **Physical source realization.**  Construct a chain-level map from a
   THM-2334/2337 or THM-2512 word/current packet to the THM-2542 chart
   cocycle, invariant under the atomic gauge and retaining target mask,
   common ancestry, owner, and clock.
2. **Nonadditive divisor-to-response realization.**  Replace pole order by a
   filtered or Picard-groupoid object which can retain both the Kummer
   meridian and the full principal-part arm.  It cannot be an additive map
   from the local-cohomology module: the sign hostile already forbids that.
3. **Integral/Frobenius control.**  Any coefficient change must explain
   (38)--(44).  Every finite coefficient reduction modulo `13^r` kills the
   canonical unit flux in the minimal model.  A derived-complete or `lim^1`
   replacement remains unadjudicated.
4. **Root and operation naturality.**  Retain multiplicities and separate
   collided roots, and distinguish the proved normal-ramification square
   from composition of Keller maps or physical time.

The most plausible replacement category is therefore not another ordinary
`H^1` group.  It is a filtered cyclic torsor carrying

```text
(meridian Kummer class, principal-part arm, connecting flux,
 root label/multiplicity, integral divisibility).                    (46)
```

Equation (24) supplies its degree-action skeleton.  Equations (29)--(44)
prove the additive characteristic-zero no-go and finite-coefficient
extinction.  They do not exclude a genuinely derived or nonadditive filtered
realization carrying the missing extension data.

## 12. Comparison with the incoming endpoint and level-three divisor work

Two concurrent packages landed while this note was being checked.

The FINITE-EXACT U_full ancestry gate in
[`lrc_endpoint_ufull_frozen_five_common_ancestry_gate_20260816.md`](lrc_endpoint_ufull_frozen_five_common_ancestry_gate_20260816.md)
shows that the frozen endpoint residues admit a synthetic nonnegative Boolean
realization and nonzero role-graph products, but the endpoint worker has
already multiplied separate marginals.  Its first missing coordinate is a
shared ancestry key before either endpoint sum.  This does not provide the
left vertical arrow in (45); it sharpens obligation 1 to an atom-level map
`J_ell(omega_L,omega_R)` before marginalization.  In particular, the positive
scalar realization cannot be used to declare the THM-2542 chart class a
physical word current.  The concurrent synthesis identifies the exact loss as

```text
ker(epsilon_L) tensor ker(epsilon_R),                  (46a)
```

the universal marginal fiber-product kernel.  A physical source realization
must select or annihilate this whole mixed-Haar kernel; appending one more
scalar to the marginals cannot recover it.

The incoming independently proof-audited fixed-map theorem
[THM-3495](../01-canon/theorems/THM-3495-level-three-sporadic-keller-norm-divisor-and-three-component-nonproperness.md)
proves a new prime divisor `J` with the exact normalization

```text
N(H)=J/(2^35 L^7).                                        (47)
```

On the complement of `LJ=0`, the constant factor is a Kummer-trivial unit
under the algebraically closed field gate, and the associated `mu_13` class
has two-component divisor residue

```text
(res_J,res_L)=(1,-7) in F13^2.                           (48)
```

Equation (48) is a useful **target hostile** for (17): the single
marked normal circle must choose a component and therefore forgets the other
residue.  The appearance of `7` in both (9) and (48) supplies no map by
itself.  The source seven is the number of chart overlaps; the target seven
is an escaping-sheet pole order.  No common base, operation, or predicate
currently identifies them.  A lawful two-divisor refinement would have to
name the map to the ordered meridian pair and show why its image is the line
spanned by `(1,-7)`, while still confronting the additive flux obstruction.

Thus the incoming work reinforces both boundaries rather than filling either
missing arrow: the LRC side needs pre-marginal ancestry coupling, and the JC
side now offers a richer Kummer divisor vector whose relation to flux remains
unconstructed.

## 13. Exact sidecar and scope

Run

```text
python -B 04-computation/d5_marked_kummer_normal_slice_frobenius_flux_no_go_20260816.py
python -B -O 04-computation/d5_marked_kummer_normal_slice_frobenius_flux_no_go_20260816.py
```

The dependency-free original companion checks:

- rank six of the `C7` coboundary and the one-dimensional seam quotient;
- the explicit `C91` primitive and all primitive-constant gauges;
- surjectivity onto all thirteen Kummer exponents;
- orientation reversal and the cover/ramification square through degree 200;
- the pole-order sign hostile and the depth `1/14` collision;
- `D(h)=6z`, `[P h]=[-x^-1]`, and the connecting representative `[-1]`;
- the integral, mod-thirteen, and localized-Frobenius identities
  (38)--(40);
- the general odd-prime telescope at thirteen primes;
- the integral depth formula through 200 and the `13^r` instances for
  `r=1,2,3`; and
- the characteristic-zero terminal recurrence.

THM-3496's independent SymPy/incidence companion separately audits the full
typing, the residue-field hostile, all twelve normalization scalars, and the
finite-coefficient-versus-derived boundary.

The calculation constructs no physical LRC current, no semantic arrival
two-cell, no nonzero Hamiltonian cross-map, no polynomial mate, no Keller
inverse, no scalar-row exclusion, and no result on LRC(14), JC(2), or DC(2).
