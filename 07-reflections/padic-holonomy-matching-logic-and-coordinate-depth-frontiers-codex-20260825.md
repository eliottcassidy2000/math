# P-adic holonomy, matching-logic sort flow, and coordinate-depth frontiers

**Session date:** 2026-08-25
**Status:** external p-adic theorem `PREPRINT CLAIM / UNDER SPECIALIST AUDIT`;
external interval replay `FINITE-EXACT`; matching-logic source `CITED / ARXIV
V1`; the referenced THM-4088--4091 and THM-4093--4096 are `PROVED` in their
exact scopes. None of these statuses proves a new p-adic irrationality
statement in repository canon.

## 1. Outcome first

The two supplied links are independent.

1. The GitHub repository is Christopher D. Long's August 24 research draft
   claiming 22 individual Kubota--Leopoldt irrationalities. Its numerical
   certificate replays exactly, but it certifies only the terminal formulas;
   several global arithmetic-geometric interfaces still need specialist audit.
2. `arXiv:2608.13306v1` is Xiaohong Chen and Grigore Rosu's paper on basic
   matching logic. It contains no p-adic-zeta material. Treating the two links
   as one bibliographic object would be a source-identity error.
3. The p-adic formulas expose a sharp structural boundary. The interior
   arithmetic optimizer exists exactly for `p<11`, explaining the levels
   `p=2,3,5,7`. Even after granting `p=13,s=3` an idealized zero-codimension
   one-power Hasse saving, its complete margin remains strictly negative:

   ```text
   ideal tau >= 613/288,       ideal margin < -37/72;
   actual inf tau = 515/224,   actual margin < -67/56.
   ```

   Thus better one-power Hasse codimension alone cannot reach `zeta_13(3)`.
4. Integral coordinate changes preserve cumulative LCM clearing at every
   depth and literal coefficientwise depth at `e=1`, but literal depth fails
   universally for every `e>=2` at the minimal output degree three. This is
   THM-4091.
5. The matching-logic paper's three-sort incompleteness witness compresses to
   exactly two sorts and one unary symbol. The semantic implication survives,
   while the hypothesis's output sort cannot feed the conclusion sort in the
   calculus. This is THM-4090, independently refereed against every primitive
   rule.
6. A bare order tournament cannot detect rational/algebraic/transcendental
   type (THM-4088). The lawful rational-edge survivor is narrower: unit
   endpoint ratios give diagonal similarity, and a tournament's first p-adic
   determinant-zeta tangent is controlled exactly by triangle-count
   divisibility (THM-4093).
7. THM-4094 replaces a selected Hamiltonian-path matching by the full deletion
   incidence and recovers the exact missing deficit. It also reduces the open
   global Hamiltonian-path spectrum conjecture to two strong prime lanes.
8. THM-4095 proves exact generated-field transport for every fixed finite LRC
   observer, while optimized primitive-pair margins form a closed countable
   gapped set, with unique accumulation point `1/6` and literal gap
   `(0,1/15)`. Bare order tournaments see neither fact; their determinant
   graph zeta is identically one because they are acyclic.
9. THM-4096 repairs the character branch in the Cover14 regularization lane:
   the first ray point has a unique positive-side rational affine optimizer,
   while every later point has a sign obstruction. This is not a sum across
   p-adic fields.

## 2. Inheritance pass and live concept board

### Anchor: named p-adic-zeta irrationality

- **Closest proved mechanism:** the published multi-value linear-form results
  routed in the August 21 p-adic-zeta reflection.
- **Canonical hostile:** `there exists an irrational coordinate` does not name
  which coordinate; a terminal positive margin does not verify the geometric
  construction that produced it.
- **Corrected near miss:** the earlier scalar-gauge route lost modular weights.
  The new draft explicitly works on a de Rham frame torsor at large primes and
  claims only one small-prime determinant saving.
- **Least-used sidecar:** the character-DFT / simultaneous Pade / Casoratian /
  Smith selector for the six `p=13` reflection orbits.

### Niche: denominator transport

- **Closest proved mechanism:** THM-4056's exact denominator/LCM compiler.
- **Canonical hostile:** an LCM address is not a nonzero small integer linear
  form, and literal `n^e` depth is not invariant under reparametrization.
- **Corrected near miss:** the smallest counterexample is
  `R=q^2/2^e`, `h=f+f^2`, output degree `3`.
- **Least-used sidecar:** the full scaled composition matrix, not only its
  diagonal denominator labels.

### Wildcard: logical localization and sort flow

- **Closest cited mechanism:** Chen--Rosu semantic localization to a
  backward-closed core and its doubled countermodel.
- **Canonical hostile:** logical localization is not localization at a prime;
  a logical double cover is not a geometric cover.
- **New positive object:** the two-sort witness of THM-4090.
- **Least-used sidecar:** a directed sort-consumer graph on proof obligations.
  It is not a tournament because most pairs have no intrinsic comparison.

The five live concepts used throughout the session were:

| Concept | Native representation | Preserved predicate | Operation | Information destroyed | Cheapest decisive test |
|---|---|---|---|---|---|
| p-adic arithmetic packet | weighted de Rham/torsor rows | local lattice and slope cost | Hasse/Frobenius transport | weights under scalarization | verify each transport lemma on its declared bundle |
| denominator staircase | coefficient lattice `D_e` | integrality after clearing | `R -> R(h)` | literal index depth when `e>=2` | hostile `(q^2/2^e,f+f^2)` |
| LCM clock | divisibility address | finite cumulative clearing | truncate at `N` | decay and nonvanishing | integer-form gate of THM-4057 |
| logical theory | result-sort-labelled formulas | semantic totality / derivability | proof rules | hypothesis lineage under sort compression | primitive-rule sort audit |
| `p=13` selector | six labelled Hurwitz/reflection orbits | named trivial-character coordinate | DFT/Pade/Smith selection | row label under unlabelled rank | trivial-character minor with height and decay |

## 3. External p-adic source audit

### 3.1 Source maturity and claimed theorem

The audited snapshot is
[`b46a177`](https://github.com/octonion/p-adic-zeta-irrationality/commit/b46a1770901551961710e155d775aae7c5ea39e7).
It contains three commits, all dated August 24, 2026. The TeX calls itself a
"product-digit repaired research draft," but the history begins with an
initial import and supplies no public correction diff. No separately indexed
paper was located. The README asks for independent specialist review.

The draft claims irrationality of

```text
zeta_2(s), odd 3<=s<=29;        14 values,
zeta_3(s), odd 3<=s<=11;         5 values,
zeta_5(3), zeta_5(5), zeta_7(3); 3 values.
```

The claimed architecture is genuinely different from the scalar route
rejected in the earlier reflection:

- below a cutoff, a Hasse-invariant source kernel claims one determinant-level
  saving by rank-nullity;
- above it, divided Frobenius and Taylor no-backflow act on the algebraic de
  Rham frame torsor, without scalarizing nonzero weights;
- modular Jensen energy and a Bost slope inequality turn a positive terminal
  margin into a rationality contradiction.

The elementary Dwork identity, Taylor coefficient no-backflow calculation,
Jensen integral, and numerical formula replay survived this audit. The
load-bearing unresolved review gates are global BGG/coarse descent,
single-unipotent deepest-row projection, genuine-source Hasse codimension and
saturated lifts, the full-product pole-digit Cartier argument, transfer of
scalar CDT height asymptotics to one fixed vector bundle, exact continuation
radii at levels five and seven, and compatibility of the local lattices with
the global Bost filtration.

Therefore the 22 headlines remain **external preprint claims** here. They do
not promote the August 21 singleton `OPEN` ledger to proved canon.

### 3.2 What the certificate does and does not certify

The repository's standard-library script uses exact `Fraction` arithmetic and
integer fixed-point intervals. On Python 3.13.14 it reproduced the stored
output line for line. Independent 120-digit checks enclosed its `pi`, `log`,
and `arccos` values; all 22 lower endpoints remained positive at multiple
precisions. The smallest stated margin is the `(p,s)=(5,5)` cell,

```text
margin > 0.1317993568270168...
```

The LF-byte manifest matches:

```text
script  1408c9092d0c34917253c8f520db56853112c58640d15c98d58b76a61a0478f5
output  9e1d8f74eb198c7e7b0c7420a1e8021d7bad4d2bc8014cc939b6b353111c0da0
```

Three verifier hardening gaps matter:

1. input and positivity gates use Python `assert`, which disappears under
   `python -O`, while the final `CERTIFIED: True` line is unconditional;
2. no `.gitattributes` fixes LF checkout, so the advertised raw hashes change
   under default Windows CRLF conversion;
3. the manifest binds the script and output only, not the TeX/PDF or the claim
   that the code matches every proof formula.

These defects do not change the signs of the immutable 22 shipped cases. They
do prevent treating the script as a hardened proof checker for arbitrary
inputs, and the calculation cannot certify the geometry or slope theorem.

## 4. THM-4089: exact optimizer boundary and `p=13` no-go

The draft defines, for `m=s+1`,

```text
K=floor(s/xi),
I=(s+1/2)H_K-K xi+[s+1-(K+1)xi]_+^2/(2(K+1)),
c_p=(p+1)/12,
J_p(xi)=integral_0^xi (1-c_p x)_+ dx,
S=s^2 xi-(s-1)J_p(xi),
tau=2(S+sI)/(s+1)^2.                                    (1)
```

In the unsaturated first chamber

```text
1<xi<(s+1)/s,        K=s-1,        c_p xi<1,             (2)
```

the second derivative of `S+sI` is positive and its stationary point is

```text
xi* = 12(s(s+1)-1)/(12s^2+(s-1)(p+1)).                   (3)
```

The exact subtraction

```text
xi*-1=(s-1)(11-p)/(12s^2+(s-1)(p+1))                    (4)
```

makes the boundary structural: `xi*>1` exactly when `p<11`; `p=11` is the
wall; for `p>11` the stationary point is outside the admissible interval.
Among primes, the interior list is exactly `2,3,5,7`. This is an optimizer
boundary, not a genus-zero classification; `X_0(13)` is also genus zero.

For the next named frontier `p=13,s=3`, make the arithmetic artificially more
favorable by replacing `J_13(xi)` with its absolute upper bound `xi`. This is
the best possible one-power Hasse cost. Complete piecewise minimization on
`1<xi<4` gives

```text
xi=11/9,       I=305/108,       S=77/9,
S+3I=613/36,  tau=613/288.                              (5)
```

The analytic margin is

```text
M=3 Lambda_13(Y)-4 tau-C_13(Y),
Lambda_13(Y)=log(13)-2 pi Y,      C_13(Y)>=0.             (6)
```

For every `0<Y<1/13`, equations (5)--(6) give

```text
M < 3 log(13)-613/72.                                    (7)
```

The first six positive terms of the exponential series give

```text
exp(8/3) > 49621/3645 > 13,
```

so `log(13)<8/3`. Hence the rigorous global bound is

```text
boxed: M < 8-613/72 = -37/72 < 0.                        (8)
```

This is stronger logically than a numerical optimizer: it excludes every
continuation parameter `Y` in the stated domain. For orientation, the
one-layer analytic stationary point is

```text
13Y=sqrt(87)/16,       sqrt(1-(13Y)^2)=13/16,
```

and the exact margin there is

```text
3 log(13)-613/72-(48pi/169)acos(sqrt(87)/16)
             -(3pi sqrt(87))/208
= -2.087947241221729472... .                              (9)
```

For the actual formal `p=13` integrand, `c_13=7/6` and
`J_13(xi)=3/7` for every `xi>1`. Its arithmetic cost is strictly increasing
across all five chambers and has the unattained boundary infimum

```text
inf_(xi>1) tau_(13,3)(xi)=515/224,                       (10)
```

so the corresponding rational margin bound is `M<-67/56`. At the same
analytic optimizer, the exact supremum replaces `613/72` in (9) by `515/56`
and is approximately `-2.770486923761`. The idealized calculation remains
useful because it proves that improving the same one-power Hasse codimension
all the way to its absolute best still cannot change the sign.

The consequence is scoped but decisive: within this margin architecture, a
better estimate for the same one-power Hasse codimension cannot prove
`zeta_13(3)` irrational. A viable continuation must provide more than one
determinant-level saving, enlarge the analytic continuation width, or change
the template/energy functional.

Independent exact derivative brackets and concave tangent bounds prove the
global next-cell sign changes

```text
p=2:  s=29 positive, s=31 negative;
p=3:  s=11 positive, s=13 negative;
p=5:  s=5  positive, s=7  negative;
p=7:  s=3  positive, s=5  negative.
```

These four negative maxima and the 22 shipped positive cells are now promoted
in THM-4089. They prove a one-step stopping boundary, not a classification of
all higher odd weights.

## 5. THM-4091: the exact coordinate-change depth boundary

Let

```text
D_e={R(q)=sum a_k q^k: a_0 in Z and k^e a_k in Z for k>=1}.
```

For every `h in f Z[[f]]`:

1. `L_N^e R(h) mod f^(N+1)` is integral for every `R in D_e`, because every
   contributing `k<=N` divides `L_N`.
2. `D_1` itself is preserved, by the derivative identity

   ```text
   (n/k)[f^n]h^k=[f^(n-1)]h^(k-1)h' in Z.
   ```

3. Every `D_e`, `e>=2`, fails universal preservation at the first possible
   output degree:

   ```text
   R=q^2/2^e, h=f+f^2,
   3^e[f^3]R(h)=3^e/2^(e-1) not in Z.
   ```

4. For a fixed `h`, preservation is equivalent to integrality of every scaled
   matrix entry

   ```text
   (n/k)^e [f^n]h^k.
   ```

This types a warning in the external draft and a broader repo frontier. The
`6! binom(n,6)` clearing in the Sturmian/AP 60-phase law, THM-4056's LCM clock,
and any Apéry-style denominator staircase may be transported cumulatively.
A literal coefficientwise depth at exponent at least two requires the full
composition-matrix sidecar. Neither form supplies decay or nonvanishing.

## 6. THM-4090: the exact one-sort/two-sort logic boundary

For sorts `b,a`, sole symbol `f:b->a`, set

```text
Gamma={forall x:b forall y:b. f(x and y)},
phi=forall x:b. x.                                        (10)
```

If a `Gamma`-model had distinct `r,s in M_b`, then `x and y` could denote the
empty set, and pointwise `f(empty)=empty` could not be total in the nonempty
carrier `M_a`. Thus `M_b` is a singleton, so `Gamma models phi`; a singleton
model proves satisfiability.

Proof-theoretically `Gamma` has result sort `a`, while `phi` has sort `b`.
The only symbol edge is `b->a`. A primitive-rule induction shows that a
`b`-sorted derived line cannot depend on the `a`-sorted hypothesis: framing is
the only premise-bearing sort-changing rule, and it only travels `b->a` here.
If `Gamma proved phi`, the empty theory would prove `phi`. But `phi` is invalid
on a two-element `b` carrier, contradicting soundness.

Chen--Rosu's cited global completeness theorem rules out a one-sort witness in
the exact basic, definedness-free, fixpoint-free language with ordinary
set-valued symbols and no nominals or set variables. Therefore the sort count
two is minimal in that scope. This does not assert completeness or
incompleteness uniformly for every two-sort signature.

## 7. What legitimately connects the two supplied sources

There is no object-level theorem bridge. Four tempting identifications are
false:

```text
logical localization        != localization/completion at p,
logical C x {0,1} cover      != geometric or 2-adic cover,
least-set mu fixpoint        != integer interval fixed point,
"hybrid" matching logic      != hybrid arithmetic holonomy.
```

The shared method is a typed representation-loss audit:

- **source:** a structured object with sort, weight, lattice, or carrier data;
- **target:** a compressed one-sort or scalar representation;
- **preserved:** selected semantics or numerical values;
- **lost:** which source coordinate can feed which consumer;
- **needed sidecar:** sort-flow, torus weight, de Rham lattice, or denominator
  matrix;
- **decisive test:** a faithful-translation/transport theorem plus a hostile
  input whose conclusion depends on the lost coordinate.

The matching-logic non-edge `a !-> b` and the p-adic refusal to scalarize
weighted rows are two independent examples of that method. Similar syntax is
not the connection; consumer-faithful transport is.

The two concurrent tournament lanes are now promoted with those declarations.
THM-4088 proves that order and p-adic-quality rankings collapse to transitive
tournaments and lose arithmetic type. THM-4093 proves that unit vertex-ratio
weights are literal diagonal similarity for determinant/Bowen--Lanford zeta;
for a tournament its first p-adic tangent has depth `3v_p(x)` exactly when `p`
does not divide the directed-triangle count. This graph zeta is not the
Kubota--Leopoldt function. Dependency and transport graphs remain sparse, so
forcing their missing pairs into a tournament would still invent data.

## 8. Updated frontier priorities

1. **Specialist p-adic audit:** mechanize one complete claimed cell, preferably
   `(p,s)=(5,5)`, from source lattice through the final slope contradiction.
   Rechecking the positive decimal is no longer the bottleneck.
2. **`p=13` selector route:** stop optimizing one-power Hasse codimension.
   Compute character-DFT simultaneous Pade/Casoratian minors, Smith saturation,
   cyclotomic norm height, and prime-primary cost. Success must name the trivial
   character and provide either an additional determinant saving or new width.
3. **Verifier hardening:** replace `assert` with unconditional checks, make
   certification conditional on all gates, bind TeX/formula metadata, add
   `.gitattributes`, and run normal plus optimized CI.
4. **Formal finite-certificate export:** the finite rational interval layer is
   a plausible target for fixpoint-free proof export. Matching-logic
   completeness would concern that finite encoding only, not analytic
   continuation or irrationality.
5. **LCM/coordinate transfer:** use THM-4091's exact matrix criterion before
   transferring the 60-phase AP law, Apéry denominators, or Duffin--Schaeffer
   clocks through a nonlinear coordinate. Preserve decay/nonvanishing as
   separate proof obligations.
6. **Decorated tournament tangent:** THM-4093 makes `p|c_3` the exact first
   cancellation wall. Any stronger arithmetic observer should retain the
   next nonzero determinant coefficient, its p-adic unit, and vertex-height
   data; the orientation or triangle count alone remains type-blind by
   THM-4088.

## 9. Late incoming frontiers and their exact interfaces

### 9.1 Full incidence, not a chosen matching

THM-4094 shows that for a tournament `T`, vertex `v`, and every Hamiltonian
path `P` of `T-v`, the number `a_v(P)` of legal insertions satisfies

```text
H(T)-H(T-v)=sum_P(a_v(P)-1)+O_v,                         (23)
```

where `O_v` counts paths whose deletion shortcut fails. Selecting one legal
extension of each `P` proves only monotonicity; it deletes both excess
extensions and orphans. The directed triangle is the minimal hostile: after
choosing one edge, the transitive triple and `C_3` have the same selected
matching but deficits zero and two.

This is the combinatorial analogue of the sort/weight/lattice sidecar rule,
with a sharper conclusion. The missing coordinate is not “another sort” in
the abstract; it is the unmatched right witnesses and their full fiber
profile. The theorem's separate reduction of the open H-spectrum conjecture
to strong realization of `p` and `7p` lanes uses SCC multiplicativity, the
proved omissions `7,21`, and the finite-prefix strong carry atoms
`49,63,343`; it does not follow from matching logic.

There is nevertheless a precise logical interface. Each fixed finite
tournament has a faithful two-sort object/path incidence, but finitary
first-order compactness says that a theory whose fibers are finite in every
model has a uniform fiber bound. Hence one theory cannot isolate all unbounded
finite path fibers while excluding infinite ones under ordinary all-structure
semantics. This compactness boundary is different from THM-4090's sort-flow
obstruction and points toward external finite-model or standard-code semantics
for a mechanized H-spectrum project.

### 9.2 Fixed-observer field transport versus optimized-family gaps

For fixed finite integer speeds `S` and rational `beta`, THM-4095 proves

```text
Q(F_S(t)-beta)=Q(t),                                      (24)
```

because an active branch writes `t=(k+epsilon F_S(t))/v`. Thus rational,
algebraic-irrational, and transcendental field fibers are each transported
exactly and are dense inside every strict witness interval. This strengthens
THM-4088's density statement before the order quotient is taken.

Changing the quantifier changes the answer. After optimizing over time and
then varying primitive pairs, the complete two-runner `1/3`-margin spectrum is

```text
{1/6} union {1/6-1/(2q):q>=3 odd},                       (25)
```

with open gap `(0,1/15)`. A fixed observer has dense typed values; the family
of optimized observers is closed, countable, gapped, and accumulates only at
`1/6`.
Every fixed observer also has a maximizing witness with denominator at most
twice its largest speed. The destroyed coordinate is the order of
“fix `S` / vary `t` / optimize / vary `S`,” not arithmetic type itself.

The exact THM-4093 connection is negative but structural. Increasing-order
prefix tournaments have strictly upper-triangular adjacency matrices, hence
`det(I-uA)=1` and graph zeta `1`; they land in THM-4093's `c_3=0`
infinite-depth branch. The order quotient has killed every cycle before any
unit vertex-ratio gauge is applied.

### 9.3 The lawful p-adic/LRC shadow is rational and character-labelled

THM-4096 studies the different Cover14 ray

```text
V_m={1,...,12,182m}.                                      (26)
```

After a rational normalization `C_m` of THM-2057's witness phase, it compares
`C_m` with the correctly typed twisted vertices

```text
L_p(-1,omega_p^2)=(p-1)/12.                              (27)
```

The affine moment equation has a nonnegative solution only at `m=1`; every
`m>=2` requires a negative moment and lies below `zeta(-1)`. At `m=1`, the
level-prime `{2,7}` optimizer puts all allowed mass at seven. Two Dedekind
observers collide there and separate strictly from `m=2` onward.

This is a lawful connection because every vertex in `(27)` is first mapped to
the same rational affine line and the character is retained. It deliberately
forgets p-adic topology. It is unrelated to THM-4089's external variables
`(p,s,xi,Y)`, despite both results containing an optimizer and a next-case
obstruction. Neither theorem identifies an `f_14` residual or proves LRC(14).

MISTAKE-506 supplies a useful branch-sensitive control. At `p=7,s=-5`, the
interpolation index six is divisible by `p-1`, so the trivial branch has the
finite value `2801/42` with valuation `-1`; at `s=-1`, THM-4096 must instead
name the `omega^2`-twisted rational vertex. In both cases negative valuation
means a denominator, not a pole.

### 9.4 Uniform ear averaging turns one-defect words into all-order growth

THM-4111 supplies the cleanest positive matching-to-tournament transfer in
this session. For a tournament `T` of order `n`, delete the new vertex `x`
from a Hamiltonian path of the cut ear `T+x_S`. If `x` is internal, the result
is a base ordering with either zero bad adjacencies and a marked site, or one
bad adjacency marked at its unique failure. Summing over the free cut bits
gives exactly

```text
sum_S H(T+x_S)=2^(n-2)((n+3)H(T)+F_1(T)),                (28)
```

where `F_1(T)` counts one-defect orderings. Thus THM-4099's first squarefree
defect layer is not merely analogous to an error term: it is the exact surplus
in the Boolean cut trace.

The connection ledger is explicit. The source is the labelled cut-ear fibre;
the target is the zero/one-defect word layer; deletion of `x` together with
its marked position is the map. It preserves total path multiplicity and the
unique repair site, but summing over `S` destroys every individual cut value
and the `(w,h)` dispersion from THM-4097. Consequently every full-cut recursive
one-witness-per-value bank has

```text
M_(n+1)>=(n+3)M_n/4,                                     (29)
```

and hence unbounded maxima under every representative choice. It does not
prove that those maxima lie in, or extend, a solid interval. The new H-spectrum
frontier is a variance/small-ball or local-step theorem for the cut quadratic,
not another proof that the right tail grows.

This also gives a disciplined negative p-adic comparison. Equation `(28)` is
an exact finite trace, while THM-4091's coefficientwise depths are a filtration.
Both warn that an aggregate can retain enough information for one conclusion
while erasing the coordinate needed for the next. There is no p-adic zeta or
irrationality consequence of `(28)`.

### 9.5 Newton totals and coefficient depths are parallel lossy quotients

THM-4103 closes a different finite response problem on the smooth,
nonresonant, theta-only JC survivor: its ramification packet is
`{1,2,2,3,3,3,7}`, its Riemann--Hurwitz total is fourteen, and its necessary
cover degrees are `{7,12,21}`. This is a response atlas, not a Jacobian
counterexample.

The exact interface with THM-4091 is methodological. Newton-polygon totals and
the Riemann--Hurwitz degree remember aggregate slope/ramification mass; the
labelled Puiseux edges remember which coefficient created it. Likewise an LCM
depth remembers the worst p-adic denominator, while the coefficientwise depth
vector remembers where it occurs. The source-to-target map in both cases is a
maximum or sum over labelled coordinates. It preserves a necessary obstruction
and destroys the attachment data needed for sufficiency. No step transports a
p-adic irrationality certificate to JC(2), which remains **OPEN**.

### 9.6 Residual-component ancestry is a two-shell 2-adic filtration

THM-4100 lowers the arbitrary-parity AP8-plus-three cutoff to `93` under the
exact survival tariff

```text
1/b+min(4/c,2/c+4/d)<=7L.                                (30)
```

The improvement comes from retaining whether the third danger family shares
an ancestral component with the second; the two branches `4/c` and
`2/c+4/d` are not interchangeable scalar estimates. Odd/even dilation makes
this a genuine two-shell `v_2` filtration: the map retains parity ancestry and
component overlap, but destroys higher residue labels and physical entry
times. The analogy to coefficientwise p-adic depth is exact at this two-level
bookkeeping scale only. It proves neither arbitrary-core closure nor LRC(14).

### 9.7 Reciprocal powers separate exact gauges from cosmetic tournaments

THM-4105--4108 provide four useful firewalls. Primitive reciprocal-exponent
commutators characterize one physical abelian-group orbit and quantify near
arrival, but do not certify loneliness. A mixed-parity pair's optimized margin
and owner imbalance reconstruct its primitive labelled ratio; cross-`v_2`
summaries on a spanning tree reconstruct a primitive row that is not all odd.
The decoder loses the common time, maximizer locations, and gcd sheets, so it
is arithmetic reconstruction rather than synchronization.

Pairwise gcd normalization of `a^b` versus `b^a` does create a genuine
tie-free tournament: an edge reverses exactly at divisibility or the ratio
`2:3`, and triangle holonomy records incompatible edge gauges. This object is
not THM-4093's diagonal endpoint gauge, whose closed-walk ratios cancel. More
decisively, AP13 and an infinite loose family have the same complete labelled
exponent tournament. The tournament therefore loses LRC margin, phase and
arrival, and THM-4111's Hamiltonian ear average has no automatic role in it.

Finally, THM-4108's reciprocal-power separation and gauge obstruction are
unconditional, while radical saturation is explicitly **CONDITIONAL ON ABC**
and ineffective. Radical support is the zero/nonzero projection of signed
valuation depth; lift dependence, AP13's zero boundary, and an unbounded
pair-sum/determinant ratio prevent transfer to LRC or to p-adic irrationality.
The useful common lesson is to retain the labelled valuation tree before
passing to radical support, not that ABC supplies a missing LRC theorem.

### 9.8 Spanning-tree saturation retains valuation depth, not safety

THM-4110 supplies the missing synchronization invariant for a sparse graph of
reciprocal phase equations. For a primitive row `v` and connected graph
`Gamma`, the circle-phase quotient has exact order

```text
I_Gamma(v)=gcd_(T spanning tree of Gamma)
             product_i v_i^(deg_T(i)-1).                  (31)
```

Primewise, its depth is the minimum spanning-tree energy

```text
v_p(I_Gamma)=min_T sum_i(deg_T(i)-1)v_p(v_i).             (32)
```

Saturation is therefore equivalent to the `p`-free vertices inducing a
connected dominating set for every prime. On AP13 the unequal-`v_2` graph has
53 edges and quotient `(Z/2Z)^6`: 64 total sheets, one physical and 63 extra.
Six same-shell anchors remove precisely that tariff. Physical orbit membership
still says nothing about `1/14` clearance.

This is a genuine interface with THM-4091. The source is a labelled valuation
vector plus the full edge-incidence lattice; the target is the finite phase
quotient. The spanning-tree minimum preserves prime locations and depths.
Radical support, one LCM, or THM-4107's edge orientation cannot recover it.
Conversely, THM-4100 begins with an actual time, so its commutators already
vanish; its missing coordinate is continuous survivor-component ancestry, not
another torsion sheet.

THM-4109 sharpens the AP7 selected-pair certificate to floors `197,232,268`
for gaps `4,8,12`. Its modulus `980=2^2*5*7^2` is an endpoint/split branch
atlas, not p-adic analysis. A composition with THM-4100 would need the length
of a specific surviving component: total positive measure and an endpoint
clock do not provide that sidecar.

There is no direct THM-4110-to-THM-4111 theorem. The shared research move is
nevertheless precise: expand the full fibre before selecting. Equation `(31)`
takes a gcd over every spanning tree; THM-4111 expands every nonconstant cut
before retaining one witness per scalar value. Choosing one tree or cut early
destroys the invariant. This is a reusable method, not an LRC/H-spectrum
identification and not an irrationality result.

## 10. Reproduction

```bash
python -B 04-computation/hybrid_padic_zeta_margin_frontier_thm4089.py
python -B -O 04-computation/hybrid_padic_zeta_margin_frontier_thm4089.py
python -B 04-computation/hybrid_padic_margin_boundary_thm4089.py
python -B -O 04-computation/hybrid_padic_margin_boundary_thm4089.py
python -B 04-computation/matching_logic_two_sort_obstruction_thm4090.py
python -B -O 04-computation/matching_logic_two_sort_obstruction_thm4090.py
python -B 04-computation/integral_coordinate_change_depth_thm4091.py
python -B -O 04-computation/integral_coordinate_change_depth_thm4091.py
python -B .scratch/gauge_p13_referee_20260825/independent_audit.py
python -B -O .scratch/gauge_p13_referee_20260825/independent_audit.py
python -B 04-computation/tournament_uniform_ear_average_growth_thm4111.py
python -B -O 04-computation/tournament_uniform_ear_average_growth_thm4111.py
python -B 04-computation/tournament_uniform_ear_average_growth_thm4111_independent_audit.py
python -B -O 04-computation/tournament_uniform_ear_average_growth_thm4111_independent_audit.py
python -B 04-computation/lrc_ap7_even_gap_overlap_atlas_thm4109.py
python -B -O 04-computation/lrc_ap7_even_gap_overlap_atlas_thm4109.py
python -B 04-computation/lrc_ap7_even_gap_overlap_atlas_thm4109_independent_audit.py
python -B -O 04-computation/lrc_ap7_even_gap_overlap_atlas_thm4109_independent_audit.py
python -B 04-computation/sparse_reciprocal_phase_graph_thm4110.py
python -B -O 04-computation/sparse_reciprocal_phase_graph_thm4110.py
```

The external certificate must be replayed from its pinned source commit; its
result is a source audit, not an in-repo dependency of THM-4089--4091.
