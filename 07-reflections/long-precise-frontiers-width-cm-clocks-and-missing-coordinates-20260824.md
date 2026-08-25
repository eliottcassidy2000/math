# Long precise frontiers: width, CM owners, clocks, and missing coordinates

**Session synthesis, 2026-08-24.** Current truth is owned by the linked canon,
not by this reflection. The session promoted four results and one `FINITE-EXACT`
sidecar. None of LRC(14), planar `JC(2)`, Euclidean four-dimensional Kakeya,
the all-order tournament inequality, or any Rule 30 prize is claimed solved.
Sun's `2-4-6-8` conjecture was already refuted; its least hole and exceptional
set remain open.

| lane | new status | exact advance | live boundary |
|---|---|---|---|
| tournaments | [THM-4051](../01-canon/theorems/THM-4051-tournament-order-seven-strong-base-exact-frontier.md), `FINITE-EXACT` | complete order-seven strong base; sharp ratio `27/8` | all orders; Paley family is now mandatory hostile |
| LRC(14) | [THM-4052](../01-canon/theorems/THM-4052-lrc14-affine-component-width-escape-cones.md), `PROVED` relative to inherited lower-dimensional input | exact `d=2` width escape and unbounded `d=2,3,4` cones | strict complementary component-intersection cones |
| planar Jacobian | [THM-4053](../01-canon/theorems/THM-4053-jc2-live-max-eight-trichotomy-and-eisenstein-survivor.md), `PROVED` on one live seam | complete weight-eight lower model; two strata excluded | three collision walls and target-compatible `p*y^2` survivor |
| sixty/dyadic clocks | [THM-4055](../01-canon/theorems/THM-4055-sixty-dyadic-response-fibre-law.md), universal `PROVED`; fixed bank `FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED` | sharp fibre-product law; `ell_29(t+60)=1-ell_29(t)` | moving observers and temporal claims |
| Sun clock audit | [finite artifact](../05-knowledge/results/sun_2468_clock_height_firewall_20260824.out), `FINITE-EXACT` | positive target factors at `60,420,27720`; exact nonzero `F_61` atlas | lift height, nearest-square terminal, global leastness |

## Inheritance pass

The session began by naming a mechanism, hostile, repaired near miss, and
least-used sidecar for each live lane.

| lane | closest proved mechanism | canonical hostile | corrected near miss | least-used sidecar |
|---|---|---|---|---|
| LRC affine boundary | THM-4041 gives exact `d=2` components; THM-4032/4030 give exact `d=3,4` defect tests | fixed endpoint banks fail at `(1,55,56)` and `(26,31,57)` | a safe midpoint is not a component crossing | closed pack-safe component versus open spoiled component |
| planar `JC(2)` | THM-4045 forces residual weight at least eight | THM-4017's special `p^4` degeneration does not attach to the generic map | a `C_60` Hasse observer loses exactly the forced second boundary jet | unique positive-genus component and its specialization degree owner |
| Rule 30 | THM-4047 gives exact fixed-column affine monodromy | the center is the moving diagonal `c_t=ell_t(t)` | fixed-tail periodicity is not temporal periodicity | dyadic phase torsor and onset certificate |
| AP/Fibonacci/triangular | THM-4038/4035 give pointed `C_60` addresses | phases `0,15` collide in the scalar Fibonacci/triangular pair | periodic template was conflated with periodic evaluated value | affine height and full pointed state |
| Kakeya `R^4` | THM-4035 supplies exact finite broad/narrow controls | the same clock gives a twisted-cubic broad spine and a planar narrow spine | a direction atlas was nearly treated as a Kakeya set | basepoint, line, multiplicity, shading, plane field, and scale |
| Sun `2-4-6-8` | THM-4027 gives universal modular solubility | THM-4026 is an empty integer fibre despite that local support | mean growth and one CRT-class census do not give global leastness | role-labelled lift quotient and exact square terminal |
| tournament strong base | THM-1950 reduces to strong tournaments | order-seven Paley tournament minimizes the complete ratio census | sampled `4.22` and monotone-room language were false | complete orbit quotient and high-symmetry extremals |

## Live concept board

The board was kept at seven objects. New ideas were compared against every
row before being promoted.

| object | representation | invariant | operation | symmetry/quotient | scale loss |
|---|---|---|---|---|---|
| affine LRC row | pack `H` plus typed exceptions | safe/spoiled component width | labelled lift and phase motion | divide common affine factor | residue data forgets component placement |
| weight-eight Keller seam | four Newton faces `L,D,V,R` | genus and target-compatible abelian factor | specialization and clutching | toric normalization | lower hull forgets generic map degree |
| fixed Rule 30 column | least dyadic word plus onset | affine monodromy `(M,C)` | time translation | `t mod 60` quotient | fixed depth does not reach moving depth |
| AP tail | sixty rational laws `L_r(n)` | owner moments and exact deficit | `n -> n+60` | phase address | phase forgets unbounded height |
| Sun fibre | role-labelled binomial tuple | exact multiplicity / triangular square | residue projection and bounded lifting | CRT phase | every finite modulus forgets lift height |
| four-dimensional directions | `P^3(F_61)` vectors | three-/four-fold wedge ranks | attach affine lines and change scale | projectivization | direction forgets placement and overlap |
| strong tournament | skew matrix modulo relabelling | `H,E_even,E_odd` | deletion/root response | isomorphism and converse | sampling forgets thin symmetric extrema |

The common lesson is not that these objects are secretly the same. It is that
each useful quotient has a different fibre, and the target usually lives in
that fibre.

## Exact advance I: complete order seven before extrapolating

THM-4051 exhausts all `456` order-seven tournament isomorphism classes,
exactly `353` of which are strong. Every strong class satisfies

```text
H(C) >= max(E_even(C),E_odd(C;1)),
```

strictly. The minimum additive slack is `21`, at representative `38` with
`(H,E_even,E_odd)=(25,2,4)`. The minimum ratio is instead

```text
H/max(E_even,E_odd)=27/8,
```

at the Paley representative `85298`, with `(189,8,56)`. Two independent
paths—subset DP/Bareiss inversion and literal `7!` permutation/Pfaffian
enumeration—agree.

This was a useful wildcard because it corrected a method-level error relevant
everywhere else in the session: random or incomplete search is least reliable
near algebraically symmetric extrema. The surviving problem is the all-order
strong base. A serious next attack should compute the Paley family
symbolically or asymptotically before proposing any growing margin.

## Exact advance II: LRC width is orthogonal to residue firewalls

On the physical `d=2` affine boundary, write the eleven-pack maximum as `M`
and the odd exceptions as

```text
alpha=ga < beta=gb,        gcd(a,b)=1.
```

THM-4052 proves that when `a+b>7`, the largest fully spoiled component has
exact length

```text
W=min(2/(7 beta),(alpha+beta-7g)/(7 alpha beta)).
```

The inherited eleven-runner phase gives a **closed** pack-safe arc of length
`1/(42M)`. Since spoiled components are open, equality already forces escape.
Thus the row is lonely when

```text
1/(42M) >= W,
```

and any hypothetical failure must lie in the strict wedge

```text
a+b>7,
beta<12M,
alpha beta<6M(alpha+beta-7g).
```

The same topological operation gives the coarser closures

```text
d=3: E>=11M,
d=4: 3E>=44M.
```

The new physical rank-eleven row with

```text
P=713721382004055345
```

meets forbidden residues `11,23 mod 56`, so THM-4049's residue firewall does
not apply, while `beta=P/3>12M` closes it by width. Residue exclusion and
component-width escape are therefore genuinely orthogonal operations, not
two encodings of the same family.

The sharp next object is not another scalar inequality. It is the literal
intersection of each connected pack-safe component with each affine spoiled
component inside the strict wedge, retaining endpoints and lift labels.

## Exact advance III: the live planar-Jacobian weight-eight trichotomy

On the live reduced `(2,3)`, `b=d=0` seam, THM-4053 inventories the complete
weight-eight residual

```text
H=lambda*p+alpha*p^2+epsilon*p^3+kappa*y^2
  +phi*p^2*y+delta*p^4+theta*p*y^2,
```

with the two inherited coefficient rows. After `Q=sigma^24`, all four lower
face functions are integral and every face has multiplicity one. The exact
support split gives three nonresonant cases.

1. `delta!=0, theta=0`: off
   `Delta_D=phi^2-4*kappa*delta=0`, the sole positive-genus factor is either a
   `j=1728` elliptic curve or absent; it cannot own positive degree to the
   good `j=0` target.
2. `delta*theta!=0`: off `delta+theta=0`, the sole positive-genus factor is
   Bolza with `Jac~E_8000^2`; it has no nonconstant map to the target `E_0`.
3. `delta=0, theta!=0`: off
   `Delta_V=phi^2-4*epsilon*theta=0`, the unique positive-genus component is
   target-compatible `j=0`. This is the only nonresonant survivor.

If the survivor owns generic degree `n`, then

```text
n=a^2-ab+b^2.
```

Every prime `2 mod 3` must therefore occur in `n` to even valuation. The
generic source genus is eight, so its ramification divisor has degree `14`.
The elliptic specialization is unramified; locating those fourteen generic
units among the rational, node, and boundary charts remains open.

The exact live frontier is now four objects, not the vague phrase “weight
eight”: the repeated `Delta_D` wall, the cuspidal `Delta_V` wall, the
`delta+theta=0` attachment collision, and the smooth theta-only survivor's
generic degree. Computing that degree is the cheapest potentially decisive
move because the Eisenstein norm filter can exclude it without another
global classification.

## Exact advance IV: one sixty clock, many incompatible consumers

### AP, Fibonacci, and triangular states

For the consecutive seven-sector AP cover, THM-4038 gives

```text
D(n+1)=(1/7) sum_(c=0)^5 A_c(r)/(n-c),
r=n mod 60.
```

The phase template is periodic, but the evaluated value is not. In fact, for
every `n>=12`,

```text
D(n+61)-D(n+1)
 =-(60/7) sum_(c=0)^5 A_c(r)/((n-c)(n+60-c)) < 0.    (1)
```

This is the cheapest hostile to deleting height from the AP state.

The full pointed Fibonacci state modulo `10` and the pointed triangular state
`(r mod 30,T_r mod 30)` are also `60`-cycles, but for different reasons. The
scalar shadows collide, and even the scalar pair leaves twelve doubletons.
They can relabel the AP phase; they do not derive its rational owners.

The terminology matters: the persistent AP owners here are rational
mechanical/Christoffel words. Calling the evaluated AP deficit an “exact
60-periodic Sturmian sequence” collapses both the rational-owner mechanism
and the unbounded height.

THM-4042 adds a second guardrail. The prime-sector clocks begin

```text
Pi_5=6, Pi_7=60, Pi_11=420, Pi_13=27720,
```

but this initial agreement with a familiar lcm staircase is not a general
law: already `Pi_17=120120`, not `lcm(1,...,16)=720720`.

### Dyadic Rule 30 response phases

THM-4055 couples `C_60` to a response of least period `p=2^e`. With

```text
g=gcd(60,p),       m=p/g,
```

the combined clock is the fibre product

```text
C_lcm(60,p) ~= C_60 x_(C_g) C_p.
```

Every fixed sixty phase hides exactly

```text
m=2^max(e-2,0)
```

response-phase states. The sharp lossless coordinate is

```text
q=(t-(t mod 60))/60 mod m.
```

It is a torsor coordinate, not a canonical product factor. A least-`p`
response factors through `C_60` exactly when `p<=4`.

Applied to THM-4047's certified fixed bank, only `82` of `100001` columns
factor through the sixty clock. The periods `8,16,32` leave one, two, and
three phase bits. The first exact hostile is stronger than a finite mismatch:

```text
ell_29(t+60)=1-ell_29(t)       for every t>=90.        (2)
```

Both companions compare the full closed `r<=29` state at times `90,94,98`,
certifying a repeat after eight steps and ruling out period four. Equation
`(2)` cannot be
pushed to the center because `c_t=ell_t(t)` changes columns and the fixed-bank
onset bounds do not uniformly reach the diagonal.

### Polynomial clocks are a third type

Under THM-4044's field hypotheses, the depth-`k` root-of-unity Hasse observer
has kernel

```text
((P^60-1)^k).
```

On the planar-Jacobian pure residual lane, its first alias is
`P^2(P^60-1)^k`, of degree `60k+2`, and the missing coordinate is the second
boundary Hasse jet. This is neither the AP height fibre nor Rule 30's dyadic
torsor. Equal clock syntax does not imply equal kernels.

## Four-dimensional Kakeya: finite four-wise transversality survives

THM-4035 uses `phi=44 in F_61`, of order `60`. Three independent clocks give

```text
(a,b,c) |-> [1:phi^a:phi^b:phi^c],
```

an exact chart of `216000` nonzero-coordinate directions in
`P^3(F_61)`. The remaining `14764` directions are boundary charts.

One clock already supplies the clean broad/narrow hostile pair

```text
B(r)=[1:t:t^2:t^3],       N(r)=[1:t:1:t],       t=phi^r.
```

All `C(60,4)=487635` four-minors of `B` are nonzero by Vandermonde; all of
those of `N` vanish because it lies in a two-dimensional vector subspace (a
projective line). Phase injectivity and cyclic equivariance survive both
representations. Four-fold transversality does not.

This is a finite algebraic laboratory, not a Kakeya set. It supplies no
basepoint, line, tube, shading, overlap multiplicity, moving three-plane,
two-ends condition, or induction on scale. The cheapest lawful next test is
to add the boundary direction charts and one line per direction, then compare
concurrent, owner-derived, and shuffled basepoints. A Euclidean attack begins
only after the resulting family carries multiscale nonconcentration and
broad/narrow data simultaneously.

## Sun's counterexample: the exact clocks are positive hostiles

THM-4026 proves that

```text
N=896315812331399
```

has no canonical-domain representation

```text
C(w,2)+C(x,4)+C(y,6)+C(z,8)=N.
```

Its controls are `a(N-1)=89`, `a(N)=0`, `a(N+1)=67`. THM-4027 nevertheless
proves every target residue modulo every modulus is represented. THM-4040
makes `N` least and unique only inside the one class
`459490 mod 1062347` through `1,001,999,999,999,999`; the global leastness
interval remains open.

The new finite sidecar checks the tempting clocks as **output moduli**. Their
true role-labelled binomial periods and target local factors are

| modulus `q` | role periods for degrees `2,4,6,8` | `N mod q` | exact local factor |
|---:|---|---:|---:|
| 60 | `(120,720,3600,7200)` | 59 | `2486/3375` |
| 420 | `(840,5040,25200,352800)` | 299 | `773146/1157625` |
| 27720 | `(55440,332640,1663200,23284800)` | 13319 | `19117792/46690875` |

All are positive. They are not the AP clocks with the same numerals. At
`p=61`, the four nonzero input roles have image sizes `(31,24,26,24)`, already
`S_2+S_4=F_61`, and exactly `212846` of `60^4=12960000` nonzero phase tuples
map to `N mod 61=21`. The counterexample is locally ordinary on this clock.

What each pure periodic phase projection used here destroys is the lift

```text
input = residue + period*height,
```

together with the coupled anisotropic bound and final square test. For fixed
`(x,y,z)`, put

```text
m=N-C(x,4)-C(y,6)-C(z,8)>=1,
D=8m+1.
```

The canonical triangular role exists exactly when `D` is an odd square at
least `9`, equivalently `D=(2w-1)^2` for some `w>=2`. The next experiment
should therefore retain each `p=61` phase tuple on THM-4026's existing
survivor stream, attach the lift quotient, and record the nearest admissible
square gap. It has hostile controls `89/0/67`. If that conditional lift
statistic does not separate the target, enlarging a bare clock should stop.

## Connection ledger

Every proposed bridge was typed by source, target, map, preserved predicate,
loss, sidecar, and cheapest decisive test.

| source | target | map | preserved predicate | destroyed information | needed sidecar | cheapest decisive test |
|---|---|---|---|---|---|---|
| LRC pack-safe phase | affine full row | `d` labelled lifts (two/three/four) | pack clearance | exception component placement | component endpoints and lift label | compare closed safe width with open spoiled width |
| AP height-phase state | `C_60` | `n -> n mod 60` | exact rational template | height and evaluated value | `n` plus owner coefficients | strict same-phase drop `(1)` |
| pointed Fibonacci/triangular state | AP phase | cycle conjugacy | successor address | AP evaluator and cause | full state, then height/owners | phases `0,15` scalar collision |
| `C_60` plus dyadic tail | response phase | CRT fibre product | common phase mod `gcd` | `m` dyadic phases | torsor `q` | `ell_29` shift-sixty complement |
| fixed Rule 30 columns | center | `c_t=ell_t(t)` | exact value with full diagonal data | uniform onset | onset/current certificate | prove a sufficient onset bound `<=r`; current theorem does not |
| `C_60` | polynomial residual | Hasse evaluation at roots | class modulo `((P^60-1)^k)` | ideal multiples and second boundary jet | degree cap or second jet | alias `P^2(P^60-1)^k` |
| one finite clock | broad/narrow directions | `B` versus `N` | phase injectivity | representation-dependent wedge rank | exponent/weight vector | all `487635` minors |
| three clocks | nonzero chart of `P^3(F_61)` | torus coordinates | nonzero directions | boundary charts and every line placement | charts, basepoints, multiplicity, scale | `216000` versus `230764` directions |
| Sun input tuple | output residue | binomial sum modulo `q` | local existence/multiplicity | integral lift and exact equality | role heights plus square terminal | `89/0/67` phase-conditioned lift audit |
| weight-eight generic map | special-fibre components | specialization | pullback-degree conservation | which component owns degree and how it attaches | degree owner and clutch | compute `n`, apply Eisenstein norm sieve |
| complete tournament quotient | incomplete sample | select representatives | values on sampled rows | guarantee of the global finite extremum | full quotient plus Paley control | complete 456-class census |

## A useful taxonomy of missing coordinates

The session exposed six distinct losses which should not all be called
“phase loss.”

1. **Torsor loss:** `C_60` hides a dyadic phase coordinate `q`.
2. **Height loss:** AP and Sun residue phases hide unbounded affine lifts.
3. **Representation loss:** the same direction clock can be broad or planar.
4. **Placement loss:** directions do not specify lines or their overlaps.
5. **Attachment loss:** a special-fibre genus ledger does not determine the
   generic degree owner without clutch data.
6. **Extremum loss:** a sample can miss a high-symmetry minimizer even when
   every computed inequality is correct.

The successful research move was to compute the actual fibre before adding a
sidecar. The AP/Sun, Rule 30, Kakeya, and polynomial-observer
counterindications show the conditions a future reusable meta-pattern would
need to retain.

## Ranked next attacks

1. **LRC anchor:** inside THM-4052's strict `d=2,3,4` cones, compute actual
   pack-safe/spoiled component intersections with endpoints and lift labels.
   Width alone is exhausted at equality; topology and placement remain.
2. **Planar-JC niche:** compute the theta-only survivor's generic degree and
   stable models on `Delta_D=0`, `Delta_V=0`, and `delta+theta=0`. Test the
   degree first against `a^2-ab+b^2` and the ramification budget `14`.
3. **Sun height lane:** add `p=61` lift quotients and nearest-square gaps to
   the exact survivor stream, with neighbours as positive controls. Do not
   enlarge the modulus before this test.
4. **Rule 30 moving lane:** combine THM-4050's marked cylinder with THM-4048's
   observer variation/carry identities to study temporal discrepancy. The
   fixed-column clock is now a hostile, not a proposed solution.
5. **Kakeya wildcard:** complete boundary direction charts and attach actual
   affine lines; then introduce multiplicity, two-ends, shadings, and nested
   scales before comparing with four-dimensional broad/narrow arguments.
6. **Tournament wildcard:** derive the Paley-family energy and Hamiltonian
   asymptotics, or find the next complete-order extremum. Any margin conjecture
   must survive the `27/8` order-seven control.
7. **Prime-sector AP lane:** study the arithmetic/distribution of THM-4042's
   explicit `Pi_P` and determine which owner-word modes and prime-adic
   stabilizers survive sparse deletion. The clock formula is known; the
   sparse consumer response is not.

## Stopping boundary

This session produced exact closures, obstructions, and experiments across
seven problems, but no syntax-only reduction among them. A clock is an
address; a determinant is a broadness predicate; a special-fibre component is
a possible degree owner; a modular solution is a local lift class; a width is
a component bound. Each becomes useful only after the consumer-specific
fibre is retained.
