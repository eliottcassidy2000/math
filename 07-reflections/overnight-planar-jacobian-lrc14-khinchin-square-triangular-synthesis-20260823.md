# Overnight synthesis: planar Jacobian, LRC(14), Khinchin, and the square-triangular Pell spine

**Status.** **PROVED** only where this document points to promoted theorem
files; **CITED** for the literature statements routed through
[`CORE-PAPERS.md`](../05-knowledge/reference/CORE-PAPERS.md);
**FINITE-EXACT** for the named companions; **SYNTHESIS** for the connection
contracts; **OPEN** for LRC(14), JC(2), and every claimed bridge between them.

The two main problems remain open.  The overnight sift nevertheless produced
four exact results and one useful negative architecture:

1. [THM-3742](../01-canon/theorems/THM-3742-square-triangular-pell-mod13-central-sign-projective-cycle.md)
   proves a complete signed Pell clock on the twelve nonzero norm fibres over
   `F_13`, and pinpoints the quotient that destroys the central sign.
2. [THM-3743](../01-canon/theorems/THM-3743-lonely-runner-polyhedron-khinchin-flatness-relation-reduction.md)
   turns any hypothetical LRC(14) counterexample into an `l1`-at-most-`356`
   Graver relation, improving the existing total coefficient cap `367` by
   `11`.
3. [THM-3744](../01-canon/theorems/THM-3744-pell-prefix-loneliness-constant-carry-exact-formula.md)
   solves the lonely-runner maximum of **every** initial Pell prefix, with a
   unique phase, exact carries, and square-triangular factorizations.
4. [THM-3745](../01-canon/theorems/THM-3745-monomial-plane-branch-conductor-triangular-pell-selector.md)
   proves that `k[F(b),bF(b)]` has triangular normalization defect and
   conductor `F^(m-1)`, with square defects selected by a Pell equation.
5. The Khinchin/continuant audit proves that digit means, denominators, and
   even equal digit products do not carry either an LRC owner/tie packet or a
   polynomial Cohn word.  The shared survivor is an **ordered Euclidean word
   plus a target-specific cocycle**.

Nothing here is a disguised proof of LRC(14) or the planar Jacobian conjecture.

## 1. Portfolio and inheritance pass

### Anchor — LRC(14)

The closest proved frontier mechanism is the THM-3701--3718 semantic packet:
mass, endpoint tariff, two detectors, adaptive atom, excluded-target literal,
and the still-missing canonical owner/root/word transport.  The canonical
hostile is the AP/far-AP pair: `{1,...,13}` and `{1,...,12,5460}` can share a
tiny relation while having radically different global geometry.  The
corrected near miss is the deep rational-defect clock: a scalar successor
marginal can be constant while the hidden detector remains nonzero.  The
least-used sidecar was the translate of a lattice slice, not another relation
count.

### Niche — planar JC through Cohn words

The closest proved mechanism is
[THM-3736](../01-canon/theorems/THM-3736-automorphic-cohn-complete-constant-sl2-polynomial-exposure-classification.md),
which completely classifies the constant `SL_2` exposure lane, followed by
[THM-3740](../01-canon/theorems/THM-3740-automorphic-cohn-one-variable-right-shear-binomial-tower-classification.md),
which closes the one-variable right-shear tower.  The canonical hostile is a
constant determinant-one matrix word whose rows are not gradients of a
polynomial pair.  The corrected near miss is “matrix determinant one implies
Jacobian one”: curl/integrability and polynomial exposure remain separate
obligations.  The least-used sidecar was cyclic multiplier holonomy.

### Wildcard — Khinchin and classic polygonal numbers

“Kinchin's content” was split rather than guessed into one object:

- Khinchin's almost-everywhere continued-fraction digit theorem;
- ordered finite continuants and Euclidean words;
- Khinchin's geometry-of-numbers Flatness Theorem; and
- unrelated namesakes (Khintchine inequalities, Wiener--Khinchin), which were
  excluded.

The canonical hostile is the dual rational expansion `[0;a_1,...,a_r] =
[0;a_1,...,a_r-1,1]`: the same rational has two finite digit lists.  The
least-used sidecar was the carry sequence.  That sidecar became the main
positive signal in THM-3744.

## 2. The live concept board and how it changed

The session maintained seven live concepts:

1. signed Pell clocks and their central involution;
2. quotient zonotopes and integral width covectors;
3. ordered continued-fraction words plus midpoint/owner cocycles;
4. Cohn matrix words plus curl and multiplier holonomy;
5. LRC physical root/owner/phase/tie/temporal semantics;
6. Pell-prefix safe phases with a constant carry;
7. collision-ring conductors and triangular colengths.

The important board updates were:

- The mod-13 clock did **not** connect directly to the THM-3713 owner packet.
  It instead exposed two inequivalent quotients: ordinary projectivization
  gives two 7-cycles, while stereographic coordinates retain the antipode and
  give one 14-cycle.
- Flatness did not produce sparse support, but the exact width formula made an
  `l1`-minimal direction a Graver element.  That split the lane into a finite
  pair-ratio atlas and a genuine higher-support partition branch.
- The Pell prefix, initially only a numerical positive control, became an
  all-length theorem when the safe interval collapsed inductively.
- The Jacobian lane rejected every constant Pell-matrix transplant.  Its live
  object is now a **longer word with at least two interacting polynomial
  factors**, not another constant continuant identity.
- Squares and triangular numbers were useful when they factored an actual
  phase numerator/denominator or conductor length.  Bare congruence overlap
  never carried a target predicate.

## 3. Four Khinchin mechanisms that must not be conflated

### 3.1 Metric digit content — CITED

For almost every real number, the geometric mean of its first `n` continued-
fraction digits tends to Khinchin's constant.  The pinned source and exact
scope are in
[`CORE-PAPERS.md`](../05-knowledge/reference/CORE-PAPERS.md#continued-fractions-and-khinchin-content).
This is an almost-everywhere asymptotic scalar.  It does not canonically assign
a finite mean to a rational and does not recover a digit order.

### 3.2 Ordered continuants — FINITE-EXACT

The companion
[`jc_lrc_khinchin_continuant_sidecar_probe_20260823.py`](../04-computation/jc_lrc_khinchin_continuant_sidecar_probe_20260823.py)
checks `34,310` exact gates.  Its decisive hostiles are:

- `[0;14]=[0;13,1]`, so a rational does not have one unqualified finite digit
  list;
- for odd `m`, `[0;1,2m]` and `[0;2,m]` have the same length, digit product,
  denominator, and final midpoint phase, but different tie counts;
- among denominators at most `200`, there are `1,222` mixed-tie signature
  groups after scalarization;
- `43/182=[0;4,4,3,3]` and `55/182=[0;3,3,4,4]` have the same digit product
  and denominator but different ordered words; and
- equal products of digits can give different constant `SL_2` words.

Thus the lawful common object is the ordered word.  LRC additionally needs a
Christoffel/midpoint/owner cocycle; JC additionally needs polynomial entries,
curl-free rows, and cyclic multiplier holonomy.

### 3.3 Khinchin flatness — CITED + PROVED ALGEBRA

This is a geometry-of-numbers theorem, historically sharing a name but not an
object with the digit constant.  THM-3743 uses the projected lonely-runner box

```text
K=pi([1/14,13/14]^13) in R^13/Rv,
Lambda=pi(Z^13).
```

Its integral dual is exactly the speed-relation lattice and its width is

```text
Lambda*={a in Z^13:a.v=0},
width_a(K)=(6/7)||a||_1.
```

The explicit general flatness estimate gives

```text
Flt(12)<=60 sqrt(26),
counterexample => ||a||_1<=70 sqrt(26)<357,
counterexample => ||a||_1<=356.
```

Choosing a shortest relation makes it a Graver element.  Support two gives
exactly `19,314` unordered reduced speed ratios with numerator plus denominator
at most `356`; support at least three is a real multiway branch, as witnessed
by `(1,-2,1)` on `(3,4,5)`, which beats every pair relation.

### 3.4 The Pell/silver word — PROVED

The phase limit in THM-3744 is

```text
alpha=1-1/sqrt(2)=[0;3,2,2,2,...].
```

This is a deterministic quadratic orbit of measure zero, not Khinchin-typical
data.  Its value comes from the ordered digit-2 recurrence and an affine carry
of one.  It is therefore a useful exact control for the distinction between
metric digit statistics and word semantics.

## 4. The mod-13 square/triangular Pell clock

THM-3742 begins with the integral matrices

```text
S_0=[[1,2],[1,1]],        A=S_0^2=[[3,4],[2,3]],
Q(x,s)=x^2-2s^2.
```

The two integral orbits give the complete square-triangular families
`T_n=q^2` and `T_N=2T_m`.  Modulo `13`, in `(x,q)` coordinates,

```text
M=[[3,8],[1,3]],          x^2-8q^2=e.
```

Each of the twelve nonzero norm fibres has fourteen points and is one
`M`-orbit; `M^7=-I` and `M^14=I`.  This is the exact fourteen-state signal.

The critical distinction is:

- direct vector projectivization `[x:q]` identifies antipodes and leaves two
  7-cycles (square-norm and nonsquare-norm slopes);
- stereographic `theta=q/(x+1)` is a bijection of the conic with
  `P^1(F_13)`, retains the antipode, and conjugates the action to one 14-cycle.

Its half-period is the antipodal involution `theta -> 5/theta`.  The reflection
`r -> -1-r` formally resembles the midpoint reflection in THM-778, but no map
assigns that state to a physical LRC runner owner.  The attractive THM-3713
offset match is therefore **wrong-torsor**, not a hidden LRC certificate.

Raw linear observers have an exact loss law.  For `L=ax+bq`,

```text
#L^(-1)(y)=1-chi(y^2-e(a^2-5b^2));
|image L|=8 or 7 according to the square character of e(a^2-5b^2).
```

Every scalar observer forgets some signed clock state.  This is the mechanism
behind the failed residue bridge, not a bad choice of one particular linear
functional.

## 5. The Pell-prefix theorem: the positive signal became structure

Let `P_0=0,P_1=1,P_(n+2)=2P_(n+1)+P_n` and

```text
A_n=(P_n-P_(n-1)+1)/2.
```

THM-3744 proves, for every `N>=1`,

```text
M({P_1,...,P_N})=A_N/(P_N+1).
```

The unique maximizer on `[0,1/2]` is that same fraction.  At the maximizer,

```text
floor(P_i t)=A_i-1,
y_(i+2)=2y_(i+1)+y_i-1,       y_i={P_i t}.
```

This `+1` carry is the missing coordinate that a scalar phase alone would
erase.  The exact distance packet is

```text
||P_i t||=
 [D-|P_(i-1)+(-1)^i P_(N-i)|]/(2D),      D=P_N+1.
```

The all-length upper proof is a genuine interval-collapse induction, not a
finite pattern extrapolation.  The exact companion has `161,773` gates through
`N=200` and independently exhausts all lower-envelope peak candidates through
`N=9`.

At `N=13`,

```text
M(S_13)=99/338,
distance numerators=
99,140,157,164,167,168,169,168,167,164,157,140,99.
```

The family is highly safe: `99/338` is more than four times `1/14`.  It is an
exact owner-labelled laboratory, not a near-counterexample.

## 6. Why squares and triangular numbers mattered here

The classical square-triangular equation

```text
T_((x-1)/2)=(s/2)^2  <=>  x^2-2s^2=1
```

is the integral Pell spine.  With

```text
a=P_(2k-1),  s=P_(2k),  x=s+a,  b=P_(2k+1),
```

one has

```text
ab=s^2+1,
A_(4k-1)=2a^2,        M(S_(4k-1))=a/x,
A_(4k+1)=x^2,         M(S_(4k+1))=x/(2b).
```

At the cannonball row `(a,s,x,b)=(29,70,99,169)`,

```text
29*169=70^2+1,
T_49=35^2,
M(S_11)=29/99,
M(S_13)=99/338.
```

This is a lawful connection because the square identities factor the exact
loneliness numerator and denominator.  By contrast, the mod-13 scalar sets

```text
squares={0,1,3,4,9,10,12},
triangular={0,1,2,3,6,8,10}
```

and their missing union values `{5,7,11}` do not identify an LRC owner or a JC
polynomial factor.  Congruence overlap is only an address unless accompanied
by a predicate-preserving map.

## 7. Planar Jacobian archaeology: what survives the sift

The planar Jacobian conjecture `JC(2)` remains **OPEN**.  The newest proved
Cohn interface is sharp enough to reject the tempting Pell shortcut:

- THM-3736 completely classifies constant `SL_2` polynomial exposure.  The
  constant Pell matrix supplies no new row-integrable polynomial pair.
- THM-3740 classifies the one-variable right-shear binomial tower.  A genuine
  open cell must contain at least two interacting polynomial factors.
- The continuant scalar probe preserves neither matrix order nor curl.
  Equal digit products can give different `SL_2` words.
- The cyclic multiplier product is necessary but not sufficient.  The exact
  `(2,1/2)` control closes its multiplier holonomy, whereas `(-2,1/2)` does
  not; neither fact alone constructs a Keller pair.

The correct JC packet for a candidate word is

```text
(ordered elementary factors,
 polynomial variable/exponent exposure,
 both row curls,
 constant determinant,
 cyclic multiplier holonomy,
 collision/injectivity witness,
 degree and infinity sidecars).
```

The cheapest live experiment is therefore a length-three or length-four word
with **two variable factors**, audited first for row integrability and only
then for Keller/collision geometry.  Another census of constant Pell powers is
saturated by THM-3736.

## 8. Typed connection contracts

| Source | Target | Map | Preserved | Destroyed | Needed sidecar / cheapest test |
|---|---|---|---|---|---|
| Pell norm conic over `F_13` | 14-state clock | stereographic `q/(x+1)` | antipode, full orbit | integral height | physical runner owner; test THM-3713 detector packet |
| Pell conic | slope residues | `[x:q]` | norm character | central sign, phase | signed lift; direct quotient is already hostile |
| LRC quotient box | speed relations | minimum integral width covector | counterexample implies resonance | support, owner, slice translate | active facets + intercept; traverse rank-11 stars |
| support-two Graver relation | Euclidean word | reduced speed ratio | exact pair arithmetic | global multiway cancellations | THM-778 owner/midpoint cocycle |
| finite CF digits | scalar digit content | product/geometric mean | coarse size | order, rational presentation, ties | keep canonical word and both endpoint conventions |
| ordered constant word | planar JC candidate | Cohn product | determinant-one matrix | polynomial exposure and curl | two variable factors + integrability audit |
| Pell prefix | LRC phase packet | `t=A_N/(P_N+1)` | all runners, order, carries | generality/extremality | perturb one speed; `169->170` fails |
| square-triangular Pell row | phase factors | `(a,s,x,b)` | exact numerator/denominator | arbitrary LRC geometry | demand a factorization of the target statistic |
| monomial-plane normalization | Pell conic | `delta=T_(m-1)=q^2` | scalar defect | `F`, branches, conductor | retain the full module/conductor packet |
| degree-nine conductor | hypothetical JC degree cell | homogenize `H`; `36 -> (H^8,H^12)` | degree/common-power/root directions | Keller pair, lower terms, Jacobian one | common-power pair has Jacobian zero; tangent cone fails component equalization |

## 9. Conductor theorem and the triangular/Pell selector

A separate JC arithmetic probe suggested the subring

```text
k[F(b), bF(b)] subset k[b]
```

for a monic polynomial `F` of degree `m`.  THM-3745 now proves, in every
characteristic,

```text
B/A=direct_sum_(i=1)^(m-1) k[X]/(X^i),
delta=m(m-1)/2=T_(m-1),
conductor=F^(m-1)k[b].                                (32)
```

If `F` is separable, the origin is geometrically an ordinary `m`-fold point;
the algebraic formulas do not require separability.  The nonsquarefree hostile
`F=b^2` is the cusp `k[b^2,b^3]`: it keeps `(32)` but loses the ordinary-node
label.

The square rows are now a lawful connection, not numerology:

```text
delta=q^2 <=> (2m-1)^2-8q^2=1,
m=1,2,9,50,289,... .                                  (33)
```

This map preserves only the scalar normalization defect.  It forgets `F`,
branch labels and incidence, the conductor ideal, and multiplication.

At `m=9`, `delta=36` and `length(B/conductor)=72`.  The cited sub-`125` JC
classification leaves only the hypothetical degree pair `(72,108)`.  After
homogenizing `H(X,Y)=X^9F(Y/X)`, the appended pair `H^8,H^12` has those total
degrees.  But `H^12` is not selected by the conductor, and

```text
(H^12)^2=(H^8)^3,
```

`36,72,108` resonance is an exact leading-form shadow, not a Keller map.
For separable `F`, the squarefree tangent cone `H` cannot enter even the
log-canonical gate: `J(H,W)=lambda H` would force all component values of `W`
to equal at the common origin, hence `W-W(0)=HV` and `J(H,V)=lambda`, which
contradicts `grad H(0)=0`.  This direct argument is the component-equalization
boundary isolated in THM-3770.  It rules out this homogeneous entry, not an
arbitrary lower-term Keller completion.
Likewise THM-3734's nine cyclotomic components are an axis plus eight disjoint
hyperbolas, whereas the separable degree-nine conductor curve glues all nine
branches at one point.  Component count survives; incidence does not.

## 10. Stopping reasons for attractive failed bridges

- **Khinchin constant -> LRC.** Stops at rational presentation ambiguity and
  tie collisions.  Strongest survivor: ordered CF word plus LRC cocycle.
- **Khinchin constant -> JC.** Stops at noncommutative word order and curl.
  Strongest survivor: ordered polynomial continuant plus holonomy.
- **Pell mod 13 -> THM-3713.** Stops at the wrong torsor: direct slopes erase
  antipodes; stereographic states have no physical owner map.
- **Short relation -> loneliness.** Stops at AP/far-AP and the safe AP
  hostile.  Strongest survivor: relation plus active slice/facet/translate
  packet.
- **Square/triangular residues -> either conjecture.** Stops because residue
  membership is an unordered scalar image.  Strongest survivor: exact
  factorization of a target statistic, as in THM-3744/3745.
- **Constant Pell matrix -> JC.** Stops at the complete THM-3736
  classification.  Strongest survivor: test interacting variable factors.
- **Pell--Chebyshev/radial profiles -> JC.** THM-3765 closes the normalized
  three-consecutive-charge ansatz, including Pell--Chebyshev profiles as
  hostile controls; independently audited THM-3771 closes the separate cubic
  radial-carrier dressing by its unequal zero-fibre addresses.  THM-3772 is
  only a **PROVISIONAL proof candidate under audit**, but its exact evidence
  says that varying both flanks still collapses to a one-product quotient.
  The triangular-conductor Pell orbit supplies no missing mate.  Strongest
  survivor: leave that quotient or add a response channel that distinguishes
  vertical components.
- **Pell prefix -> hard LRC frontier.** Stops because `M(S_13)=99/338`, far
  from `1/14`.  Strongest survivor: exact carry/owner laboratory.

## 11. Next decisive pulls

1. Traverse the existing THM-2052 rank-eleven star templates and ask whether
   their relation spaces contain an `l1`-at-most-`356` Graver vector.  Outside
   that span, THM-3743 forces the mixed rank-twelve cofactor terminal.
2. For each surviving short relation, retain the minimizing flatness facets,
   slice intercept, induced lattice, and two-torsion centre.  Test whether this
   packet determines any THM-3718 owner/root/word class.
3. Use THM-3744's constant-carry profile as a positive control for proposed
   owner-aware LRC transducers.  Any transducer that cannot reconstruct its
   symmetric `N=13` packet is too lossy for the frontier.
4. In JC, enumerate only the first unsaturated Cohn cells with two interacting
   polynomial factors.  Do not re-enter the audited THM-3765/3771 radial
   cells; keep THM-3772 conditional until independent promotion, then move
   beyond its one-product quadratic quotient or add a component-address
   response channel.  Gate in the order: exposure, both curls, determinant,
   holonomy, then collision/infinity.
5. Starting from proved `(32)`, ask whether its conductor filtration can enter
   an actual bivariate Keller collision ring.  First compute the actual-target
   linearized Jacobian image and its cokernel, as in THM-3737's retained-value
   hyperplane; do not substitute a SAGBI semigroup conductor for an unknown
   global conductor.  Reject any degree-only shadow before lower terms, both
   curls, and Jacobian one are supplied.
6. For the mod-13 clock, test every candidate physical owner map against the
   central sign and THM-3713's two detectors before doing any larger modular
   census.

The synthesis is therefore not “Pell solves both problems.”  The durable
lesson is narrower and more useful: **ordered arithmetic becomes a carrier
only after the relevant carry, owner, slice, curl, or holonomy coordinate is
retained.**
