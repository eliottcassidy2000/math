# The five old and three new branches are one indivisible `3:5` packet

**Status:** PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED as
THM-2713.  The result upgrades THM-2704's generic integrality to a uniform
geometric-integrality theorem on the whole nonzero-`eta` split even-Faber
prime-23 family.  It does not exclude rational specializations, restore the
odd seeds, or prove the split branch or `JC(2)`.

Companion:

```text
04-computation/jc2_split_prime23_component_divisor_budget_20260728.py
05-knowledge/results/jc2_split_prime23_component_divisor_budget_20260728.out
script SHA-256 81cb01e4dff454fc18417b3cfa2136fcab97901317d2b8d2df805a5a0c677350
output SHA-256 483754f20a492b3a77892cfea2f6c06612459d4366dc52647ce3d07cd4b9468b
```

## 1. Inheritance pass and the missing object

The closest proved mechanism is THM-2704.  It gives the exact weighted
complete intersection

```text
F2=0,                    zeta F1^4=eta t^23,
deg(F2,F1)=(6,5),        eta!=0,                      (1)
```

with five fixed `(4,23)` cusps, three fixed smooth branches, arithmetic genus
`254`, and generic normalization genus `89`.

The canonical hostile is not a bad generic specialization.  It is a special
member that becomes reducible, nonreduced, or gains exactly `89` extra delta
units.  THM-2636/2692 handled the nonsplit curve by enumerating fixed-factor
types and killing their finite Hensel carriers.  Here the fourth-root Kummer
cover is globally trivial, and the old eliminant has degree `23`, so copying
that enumeration looked much larger.

The corrected near miss is the degree-six `C2` order-raising route.  The exact
Keller differential row still sees only the first odd tooth; it supplies no
specialization-stable lattice that kills the eleven omitted odd seeds.  That
problem is orthogonal to the even-Faber curve and must not be smuggled into
the present proof.

The least-used sidecar was not another Hensel coefficient.  It was the degree
of the two *sections* whose local orders THM-2704 had already computed:

```text
F1 is a section of O(5),             zeta is a section of O(3). (2)
```

Their degree budgets point in opposite directions.  That opposition is the
missing global carrier.

The live concept board was

```text
five old branches O; three new branches N; the weight-one pencil [h:t];
the O(5) section F1; the O(3) section zeta; the order-23 norm coupling;
the rational q-map and the omitted third flux.
```

## 2. Hostile audit of the pencil before counting anything

The component argument needs more than the local word `4^5,1^3`.  Each of
the following gates was checked explicitly.

### Ambient and line-bundle gate

The two quotient-singular points of `P(1,1,2,3)` are avoided uniformly,
because `F2` has the fixed nonzero `v^3` and `zeta^2` values there.  Thus
`O(1)` restricts to an honest line bundle on every member.

### No common hypersurface and no embedded component

The exact `t=0` support is finite in `h=1`, and the already-frozen corner
Groebner calculation gives no point at `h=t=0`.  If `F2` and the degree-23
equation shared a hypersurface, that projective surface would meet `t=0` in a
curve, a contradiction.  The two equations are therefore a regular sequence.
The curve is lci/Cohen--Macaulay, pure of dimension one, and has no embedded
associated point, even when the degree-six surface is special or singular
away from the fixed fibre.

### No vertical component

The sections `h,t` have no common zero, so they define

```text
pi:C -> P1,                         O_C(1)=pi^*O_P1(1). (3)
```

If `pi` were constant on a component, the right side would be trivial there;
the left side is ample of positive degree.  Hence every component maps
nontrivially and surjectively to the `t`-line.  This simultaneously excludes
an affine constant-`t` component and a component hidden entirely at infinity.
It does not use generic infinity transversality.

### Local-domain and multiplicity gate

At a new `G3` point the local equation is

```text
zeta=unit*t^23,               (ord t,ord zeta,ord F1)=(1,23,0). (4)
```

It is regular.  At an old `L5` point it is

```text
F1^4=unit*t^23,               (ord t,ord F1,ord zeta)=(4,23,0). (5)
```

The final zero in (5) belongs to `zeta`; `zeta` is a unit there.  Since
`gcd(4,23)=1`, (5) is one reduced analytic branch, not four branches being
silently grouped.  Every fixed point therefore belongs to one and only one
global component.

After the divisor budget leaves one topological component, the regular new
points make it generically reduced.  The lci is `S1` and has no embedded
prime, so generic reducedness upgrades to reducedness.  This is why the
conclusion is geometric integrality rather than merely connected support.

## 3. The opposing divisor budget

Take the normalization of one geometric component.  Suppose it contains
`r` old and `s` new branches.  Because the map (3) is surjective, it contains
at least one fixed branch, and its degree is

```text
d=4r+s.                                                (6)
```

Neither `F1` nor `zeta` can vanish identically: equation (1) would then kill
the nonzero generic section `eta*t^23`.  Their visible zeros give

```text
23r <= 5d,                       23s <= 3d.            (7)
```

Substitution of (6) turns (7) into

```text
3r <= 5s,                        5s <= 3r.             (8)
```

Hence `3r=5s`.  The global bounds `r<=5,s<=3` leave only

```text
(r,s,d)=(5,3,23).                                      (9)
```

This is substantially sharper than degree counting alone.  Before the two
section budgets, `4r+s` realizes every integer from `1` through `23`; after
them, no proper degree survives.  Every component would have to contain all
eight local domains, so there is exactly one.

The structural slogan is:

> weight three prevents a component from taking too many new branches;
> weight five prevents it from taking too many old branches; order twenty
> three makes the two inequalities meet only at the full `3:5` packet.

This explains the otherwise surprising primes.  It is not numerology: the
Diophantine equality is literally `3r=5s`, and the only available packet has
three new and five old places.

## 4. The norm form and why order 23 really is the first coupling

Write

```text
f1=a*zeta+b,
f2=A*zeta^2+B*zeta+C0,                  A=15944049,
D=A*b^2-B*a*b+C0*a^2.                                 (10)
```

The companion proves in an abstract polynomial ring that

```text
Res_zeta(f2,zeta*f1^4-X)=C0*D^4-X*T+A^5*X^2.          (11)
```

Here `D=Res_zeta(f2,f1)` and `T` is the scaled quadratic trace.  On the
prime-23 curve,

```text
R=C0*D^4-eta*t^23*T+A^5*eta^2*t^46.                  (12)
```

The exact supports are

```text
deg_(t,v)D=(10,5), terms=60;
deg_(t,v)T=(23,11), terms=1220.                       (13)
```

At the fixed section,

```text
C0=-672G3,             D=-16071601392L5,
gcd(T,G3L5)=1.                                         (14)
```

Thus the `G3*L5^4` product is uncoupled through order `22`, and the first
coupling at order `23` is a unit on every branch.  This is the correct exact
Hensel instrument, but the section-degree proof is stronger for components:
it kills every proper partition without enumerating 127 complementary
subsets.  Formula (12) should now be spent on the rational perfect-power
residual, not on repeating the component classification.

## 5. Equality gives global divisors, not just inequalities

For the unique component, both inequalities in (7) are equalities.  Therefore
there are no unrecorded zeros:

```text
div_0(F1)=23O,                   div_0(zeta)=23N.      (15)
```

The split-sheet reconstruction is rational and, up to a nonzero constant,

```text
q=-t^5/F1.                                             (16)
```

Since `div_0(t)=4O+N`,

```text
div(q)=5N-3O.                                          (17)
```

So `q` has exactly three order-five zeros and five order-three poles, and no
other zero or pole anywhere on the normalization.  The Kummer cover is not
merely trivial; its trivializing function has a completely frozen divisor.

Uniform lci adjunction gives `p_a=254`.  The five fixed cusps spend `165`, so

```text
g=89-Delta_extra.                                      (18)
```

The reducible/nonflat escape from the old singularity-budget statement is
gone.  A rational member must carry exactly `89` additional delta units.

## 6. A rational survivor is a simultaneous `3,5,23` perfect-power map

If the normalization is `P1`, pull back `O(1)` as `O_P1(23)`.  Let `alpha`
and `beta` be the squarefree binary forms cutting out `N` and `O`; their
degrees are `3` and `5`.  Equations (15) and the `t`-fibre force

```text
t=tau*alpha*beta^4,
F1=phi*beta^23,
zeta=sigma*alpha^23.                                  (19)
```

Let `H` be the other weight-one section and `V` the weight-two section.  Then

```text
(h,t,v,zeta)=(H,tau*alpha*beta^4,V,sigma*alpha^23),
deg(H,V)=(23,46),
gcd(H,alpha*beta)=1,
sigma*phi^4=eta*tau^23.                               (20)
```

Affinely, and suppressing nonzero constants,

```text
q=alpha^5/beta^3,
t=alpha*beta^4/H,
zeta=alpha^23/H^3,
f1=beta^23/H^5.                                       (21)
```

The first curve equation has become tautological, but the identities

```text
F2(H,t,V,zeta)=0,                 F1(H,t,V,zeta)=phi*beta^23 (22)
```

remain highly overdetermined.  This is the exact survivor, not a claimed
contradiction.

There is a second useful view.  The map `q` has degree `15`; its visible
zero/pole fibres contribute

```text
3*(5-1)+5*(3-1)=22
```

of the genus-zero total `28`, leaving only six other ramification units.
The degree-23 map `t` leaves `29`.  The `q`-Hurwitz space is therefore the
smaller residual object.  Its logarithmic derivative has the familiar form

```text
q'/q=5 alpha'/alpha-3 beta'/beta,                    (23)
```

whose numerator has degree at most six after the leading terms cancel.  The
six-unit invoice is built into the `3:5` perfect-power form rather than being
an independent coincidence.

## 7. Relation to the older routes

- **THM-2636/2692:** their root/pair Hensel atlases were necessary because the
  nonsplit degree-five eliminant had no opposing section budget.  In the
  prime-23 split curve, the `O(3)` and `O(5)` sections collapse all component
  partitions before a factor-field census.
- **THM-2704:** generic genus `89` remains correct, but generic integrality was
  weaker than the actual statement.  The new proof makes integrality uniform
  while leaving genus drop genuinely special.
- **Degree-six `C2` stopping certificate:** the odd-seed order-raising repair
  still fails on the first canonical grade.  Uniform integrality of the
  even-Faber curve supplies no missing divisibility for `E3,E1` and does not
  restore the eleven odd channels.
- **Tournament or graph language:** no intrinsic binary relation is needed.
  The correct carrier is a pair of section divisors with opposite degree
  inequalities.  Forcing a tournament would discard the multiplicities
  `23`, weights `3,5`, and the line bundle that make the proof work.

## 8. Next decisive work

The frontier is now narrower and more algebraic.

1. Substitute (20) into the exact `F1,F2` rows and reduce successively modulo
   `alpha`, `beta`, and `H`.  Retain all scalar units and the third flux.
2. Use the norm form (12) only at its first coupling order `23`; test whether
   (22) forces a forbidden common factor, a degree overflow in `V`, or a
   finite-free quotient killed by one good prime.
3. Treat the `q`-map as a degree-15 Hurwitz object with profiles `5^3` and
   `3^5` and only six remaining ramification units.  Classify the possible
   six-unit monodromy packets before returning to the much larger `t`-map.
4. If a hostile perfect-power solution survives the first two fluxes, spend
   the third flux there.  Do not call failure of the inherited `C2` lattice a
   reason to omit this exact even-Faber residual.

## 9. Reproduction

```bash
python3 04-computation/jc2_split_prime23_component_divisor_budget_20260728.py
python3 -O 04-computation/jc2_split_prime23_component_divisor_budget_20260728.py
```

Both modes byte-match the frozen output.  All identities, gcds, support
counts, divisor-budget cases, and exponent checks are exact.
