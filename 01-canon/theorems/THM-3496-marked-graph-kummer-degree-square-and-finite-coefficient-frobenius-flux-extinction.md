---
id: THM-3496
title: "Marked graph--Kummer degree square and finite-coefficient Frobenius flux extinction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED, WITH SCOPE REPAIRS.
  After fixing the source orientation and deck generator, an algebraically
  closed characteristic-zero residue field, an oriented target uniformizer,
  and the exponent-one Kummer torsor, seam sum gives a normalized isomorphism
  from the LRC seven-chart graph H1 line to the target normal Kummer H1 line.
  Cyclic graph pullback and normal ramification by degree k commute; k=13
  kills and k=14 restores the class.  For P=x+x^2z, the characteristic-zero
  unit response is nonzero but becomes exact over Z/13^r for every finite r.
  This is a marked divisor-normal statement, not a physical LRC current,
  Keller-composition, Hamiltonian-flux realization, derived-completion no-go,
  LRC(14), JC(2), or DC(2) theorem.
source: codex2 independent D5 audit, 2026-08-16
audit: >
  independent finite-field incidence and C91 prefix integration; Kummer unit
  hypothesis and Q((lambda)) hostile; normalization-torsor and degree-square
  audit; independent SymPy Hamiltonian differentiation, odd-prime Frobenius,
  finite-coefficient, terminal-recurrence, normal/-O/stored, hash, type, and
  scope audit
depends_on:
  - THM-2542-seven-chart-cech-holonomy-and-c91-arrival-obstruction
  - THM-3431-d5-secondary-h1-descent-defects-and-valuation-persistence
  - THM-3406-affine-modification-power-jets-and-principal-part-transgression
related:
  - THM-3354-inequivalent-h1-carriers-and-typed-obstruction-cospan
  - THM-3412-hamiltonian-principal-part-differential-and-prufer-torsion-arms
  - THM-3450-marked-d5-carrier-isomorphism-and-full-germ-margin-obstruction
  - HYP-9031-d5-h1-dictionary-lrc-word-current-vs-jc-flux
script: 04-computation/d5_marked_graph_kummer_finite_coefficient_flux_thm3496_independent_audit.py
output: 05-knowledge/results/d5_marked_graph_kummer_finite_coefficient_flux_thm3496_independent_audit.out
script_sha256: 75a231dd17d2bc847c01700c6c4665b25ac841dc07e52fdfe7433d8cc4385acf
output_sha256: e2aad9e9298afad36ae00abeb4d10dcb92a733902f75a8608891bda991c093e6
semantic_sha256: 0c5756987947c93c271b7db5443a1618e483269c490b6710c30b9056562be1bb
hash_basis: LF-normalized bytes
---

# THM-3496 -- marked graph--Kummer degree square and finite-coefficient Frobenius flux extinction

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED, WITH SCOPE REPAIRS.**

## 1. Exact statement and markings

Orient the seven-cycle and write

```text
C_L^0=F_13^7,             C_L^1=F_13^7,
(delta f)_i=f_(i+1)-f_i.                                  (1)
```

Let `K` be an algebraically closed field of characteristic zero, put

```text
E=K((lambda)),             U_lambda=Spec(E),               (2)
```

and fix all four markings

```text
source orientation and the deck generator tau(j)=j+7;
the selected target divisor and its uniformizer lambda;
the target meridian orientation;
the exponent-one torsor y^13=lambda.                       (3)
```

Then the seam functional

```text
s(g)=sum_(i in C_7)g_i                                    (4)
```

and the Kummer generator

```text
kappa_lambda=[y^13=lambda]
  in H^1_et(U_lambda;mu_13)                                (5)
```

define a normalized isomorphism

```text
Phi_lambda:H^1_graph(C_7;F_13)
  -> H^1_et(U_lambda;mu_13),
Phi_lambda([g])=s(g)kappa_lambda.                          (6)
```

For the degree-`k` cyclic graph cover `p_k:C_(7k)->C_7` and the
degree-`k` normal ramification

```text
r_k:Spec K((t))->Spec K((lambda)),       lambda=t^k,       (7)
```

the marked square commutes:

```text
H^1(C_7;F_13)       --Phi_lambda--> H^1_et(K((lambda));mu_13)
      | p_k^*                              | r_k^*
      v                                    v
H^1(C_(7k);F_13)    --Phi_t------> H^1_et(K((t));mu_13).   (8)
```

Both paths multiply the marked exponent by `k`.  Degree thirteen kills the
class and degree fourteen restores it.

For the integral Hamiltonian model

```text
P=x+x^2z,
D=(1+2xz)partial_z-x^2partial_x,                          (9)
```

the characteristic-zero class `[1]` is nonzero in

```text
C_Z=Z[x,z]/D(Z[x,z]),                                    (10)
```

but its image is zero in the response cokernel over `Z/13^r` for every
finite `r>=1`.  Equivalently, `[1]` lies in `13^r C_Z` for every `r`.

The quantifiers in (2)--(3) and the word **normal** in (7) are load-bearing.
There is no maximality or unmarked-canonicity claim, and the finite-coefficient
statement makes no universal claim about derived completion, `lim^1`, or all
possible constructions called Bockstein towers.

## 2. The graph quotient and the C91 deck defect

There are no two-cells in the cycle graph, so every edge cochain is a
cocycle.  The incidence matrix in (1) has rank six.  Every coboundary has
seam zero by telescoping.  Conversely, if `s(g)=0`, prefix integration around
the seven-cycle gives `g=delta f`.  Hence

```text
H^1_graph(C_7;F_13)=C_L^1/delta C_L^0
  --s (isomorphism)--> F_13.                              (11)
```

Pull `g` back to the degree-thirteen cover `C_91`.  Its total seam is
`13s(g)=0`, so choose `h` with

```text
delta h=pi^*g.                                            (12)
```

The primitive is unique up to a constant.  With the forward action convention

```text
(tau^*h)(j)=h(j+7),                                       (13)
```

its deck defect is constant and equals

```text
chi_g(tau)=tau^*h-h=s(g).                                 (14)
```

Replacing `g` by `g+delta f` replaces `h` by `h+pi^*f`; the added term is
deck invariant.  Thus (14) is gauge independent and recovers THM-3431's
isomorphism

```text
H^1_graph(C_7;F_13) -> H^1_group(C_13;F_13).             (15)
```

For THM-2542's constant chart cochain `g_a=(a,...,a)`, equations (4) and
(14) give

```text
s(g_a)=chi_(g_a)(tau)=7a.                                (16)
```

This is still an `F_13` class.  The carrier length `91` does not turn it into
a mixed `Z/91` coefficient class or a physical endpoint current.

## 3. The Kummer unit gate and its sharp hostile

The Kummer sequence applies because `13` is invertible in `E`:

```text
H^1_et(E;mu_13)=E^*/E^(13).                              (17)
```

The valuation decomposition is

```text
E^*=lambda^Z times K[[lambda]]^*.                        (18)
```

Every unit in (18) has a thirteenth root.  Its constant residue has one
because `K` is algebraically closed.  After dividing by that root, the
remaining principal unit has a unique root by formal Hensel lifting: the
derivative of `T^13-u` at the residue root is a unit in characteristic zero.
Therefore valuation modulo thirteen gives

```text
E^*/E^(13) = Z/13 * kappa_lambda.                        (19)
```

The residue-field assumption cannot be suppressed.  Over `Q((lambda))`, the
unit `2` has valuation zero but is not a thirteenth power: the exponent of
the prime `2` in a rational thirteenth power is divisible by thirteen.  Thus
`[2]` is an additional Kummer class.  Over such a field, the line generated
by `[lambda]` still exists, but (19) is false and changing the uniformizer by
an arbitrary unit need not preserve that line generator.

Under (2), an orientation-preserving change `lambda'=u lambda` preserves
`kappa_lambda` because `[u]=0`; orientation reversal sends it to its inverse.
Changing `y` by a deck element changes a torsor representative, not its class.

## 4. What marks and naturality do -- and do not -- select

Equations (11) and (19) identify two one-dimensional `F_13` groups.  Before
the exponent-one normalization in (3), there are exactly twelve nonzero
isomorphisms between them.  Every one commutes with degree multiplication and
orientation reversal.  Consequently neither the degree square nor orientation
naturality selects a scalar.

The last marking in (3) selects the scalar-one map (6): a source class with
seam exponent `c` is sent to the torsor exponent `c`.  This proves uniqueness
**after normalization**, not maximality among cross-domain correspondences.
Moreover, (6) is not induced by a morphism from the chart nerve to the JC
normal slice.  It is a marked correspondence through their cyclic exponent
lines.

## 5. The degree square

The pullback cochain on `C_(7k)` repeats the seven source edges `k` times, so

```text
s(p_k^*g)=k s(g).                                         (20)
```

On the Kummer side,

```text
r_k^*kappa_lambda=[t^k]=k kappa_t.                       (21)
```

Equations (20)--(21) prove (8).  In particular,

```text
k=13:  both images are zero;
k=14:  both images equal the original marked class.       (22)
```

The target operation in (7) is ramified base change in one chosen normal
coordinate.  It is not composition of Keller maps, multiplication of generic
fiber degrees, a polynomial mate, or physical time.  Thus (8) transfers only
the divisor-meridian exponent.

## 6. Why this Kummer map does not realize additive flux

For a characteristic-zero Hamiltonian response complex

```text
K_P=[R --D_P--> R],              C_P=R/D_P(R),            (23)
```

the additive group of `C_P` is a vector space over a characteristic-zero
field.  It is torsion-free and divisible.  The source and target of (6) have
exponent thirteen.  Therefore

```text
Hom_Add(F_13,C_P)=0,             Hom_Add(C_P,F_13)=0.     (24)
```

The first equality follows from torsion-freeness.  For the second, write any
`v` as `13(v/13)` and use that the target has exponent thirteen.

Pole order cannot repair (24).  The principal parts `xi_q` and `-xi_q` have
the same pole order, whereas an additive map must send their classes to
opposites.  Also `q` and `q+13` have the same Kummer exponent while retaining
different principal-part depths.  Kummer records the meridian exponent
modulo thirteen; it does not carry the additive response module.

## 7. The exact Frobenius hostile

For (9), take the Bezout row

```text
A=1-2xz,                 C=4z^2,
A P_x+C x^2=1.                                             (25)
```

The localized class represented by

```text
h=-A/x^2=-x^(-2)+2z/x                                    (26)
```

satisfies

```text
D(h)=6z,
P h=-x^(-1)+z+2xz^2,
D(-x^(-1))=-1.                                            (27)
```

After base change to `Q`, the principal-part quotient has
`[P h]=[-x^(-1)]`.  THM-3406's connecting map therefore sends this cycle to
the unit response class `[-1]`.  The identities themselves are integral, so
we now study their class in the integral cokernel (10).

Put `y=xz` and

```text
m_j=x^j z^(j+1).                                          (28)
```

Direct differentiation gives the recurrence

```text
D(m_j)=(j+1)y^j+(j+2)y^(j+1).                            (29)
```

The derivation raises the weight `wt(x)=1`, `wt(z)=-1` by one.  A polynomial
primitive of `1` could therefore be reduced to a finite weight-minus-one sum
`sum_(j=0)^N c_j m_j`.  Equating successive powers of `y` forces

```text
c_j=(-1)^j,                                               (30)
```

but leaves the terminal coefficient `(N+2)(-1)^N`, which is nonzero in
characteristic zero.  Hence `[1]!=0` in (10).

For every `n>=1`, define

```text
Q_n=sum_(j=0)^(n-1)(-1)^j x^j z^(j+1).                   (31)
```

Telescoping (29) gives the integral identity

```text
D(Q_n)=1+(-1)^(n-1)(n+1)(xz)^n.                          (32)
```

Thus in `C_Z`,

```text
[1]=(-1)^n(n+1)[(xz)^n].                                 (33)
```

Taking `n=13^r-1` yields, for every `r>=1`,

```text
D(Q_(13^r-1))=1-13^r(xz)^(13^r-1),
[1]=13^r[(xz)^(13^r-1)].                                 (34)
```

After coefficient reduction modulo `13^r`, equation (34) becomes

```text
D(Q_(13^r-1))=1.                                         (35)
```

This proves both infinite thirteen-divisibility of the nonzero integral
class and its exact extinction in every finite coefficient quotient.

For any odd prime `p`, the specialization `n=p-1` gives

```text
Q_p=sum_(j=0)^(p-2)(-1)^j x^j z^(j+1),
D(Q_p)=1-p(xz)^(p-1).                                    (36)
```

Over `F_p`, finite geometric summation and Frobenius give

```text
Q_p-x^(-1)=-P^(p-1)/x^p.                                 (37)
```

The right side is a localized first integral: `D(P)=0`, and
`D(x^(-p))=0` in characteristic `p`.  Its image in the principal-part
quotient is `[-x^(-1)]`, so the connecting unit flux vanishes after this
coefficient change.  At `p=13`, this is the required hostile to extending
(6) to the Hamiltonian response.

## 8. The finite-coefficient conclusion is not a universal Bockstein no-go

Let `M=C_Z`.  Equation (34) says

```text
[1] belongs to intersection_(r>=1) 13^r M.               (38)
```

Consequently the ordinary completion map

```text
M -> inverse_limit_r M/13^r M                            (39)
```

kills `[1]`; equivalently, its image is zero at every finite coefficient
level.  The explicit primitives in (35) prove the same statement directly in
the reduced response complexes.

Nothing here proves that every derived `13`-complete enhancement, every
`lim^1` extension, or every spectral sequence called a Bockstein tower loses
all information about the integral class.  Indeed, (38) says precisely that
the integral cokernel is not separated at this element, so derived extension
data are the natural place where information could remain.  The stronger
wording in the original reflection is repaired to the finite statement
(35)--(39).

## 9. Connection and loss ledger

| field | exact content |
|---|---|
| source | oriented `C_7` graph cocycles modulo vertex coboundaries |
| source sidecar | degree-thirteen `C_91` primitive and marked deck generator |
| target | `mu_13` Kummer line of an oriented punctured formal normal slice |
| map | `Phi_lambda([g])=s(g)kappa_lambda` after exponent-one normalization |
| preserved | zero/nonzero seam, exponent mod thirteen, graph-cover/normal-ramification degree action |
| destroyed | word amplitudes, positivity, target masks, ancestry, owner/clock, root multiplicity, arm depth, additive flux, Keller predicate |
| sharp field hostile | `Q((lambda))`, where `[2]` is an extra unit Kummer class |
| sharp flux hostile | `P=x+x^2z`, equations (25)--(37) |
| finite conclusion | `[1]=0` over `Z/13^r` for every finite `r` |
| unproved | physical source realization, additive/derived target realization, Keller composition, LRC(14), JC(2), DC(2) |

In particular, no map is constructed from THM-2334's relation current,
THM-2337's gauge-dependent word jet, or THM-2512's doubly-centred interaction
to the chart class in (11).  No common ancestry atom or semantic arrival
two-cell is supplied.  On the target side, normal ramification of a selected
divisor is not composition in the Keller degree monoid and says nothing about
injectivity or polynomial inverses.

## 10. Independent exact audit

The companion uses a route independent of the original unnumbered sidecar:
SymPy finite-field matrices and differentiation rather than its custom sparse
arithmetic, plus direct prefix integration on all seven basis cochains of
`C_91`.  It checks:

- rank six and all `49` basis/vertex-gauge deck-defect cells;
- the constant defect `7`, all scalar degree squares through `k=500`, and the
  twelve unmarked nonzero isomorphisms;
- a formal thirteenth root of a nontrivial principal unit through order eight
  and the `Q((lambda))` unit hostile;
- the Bezout row and all three principal-part identities in (27);
- the weight recurrence and terminal obstruction through `201` terms;
- (36)--(37) at thirteen odd primes; and
- direct coefficient extinction modulo `13` and `13^2`, with (32) supplying
  the all-`r` proof.

Reproduce with

```text
python -B 04-computation/d5_marked_graph_kummer_finite_coefficient_flux_thm3496_independent_audit.py
python -B -O 04-computation/d5_marked_graph_kummer_finite_coefficient_flux_thm3496_independent_audit.py
```

Both runs byte-match the stored transcript.  **QED.**
