# The apex-cubic proof is a boundary-state compiler with a sharp defect budget

**Research synthesis, 2026-08-25.** The external theorem in
[arXiv:2608.22870v1](https://arxiv.org/abs/2608.22870) is a one-day-old,
computer-assisted preprint claim and remains **CITED / UNDER AUDIT** here.
The internally proved consequences are THM-4113 and THM-4116. HYP-9080 and
all cross-problem proposals below remain **OPEN** unless a stronger status is
displayed explicitly.

The decisive abstraction is not “snarks behave like difficult primes.” It is

```text
bounded local defects
  -> a planar core with ordered interfaces
  -> a native boundary-move safety game
  -> a finite reducibility atlas
  -> a global surplus beating all local defect tariffs
  -> collision-complete quotient descent.                         (1)
```

That compiler overlaps strongly with several repo mechanisms, but each
transfer needs its own state carrier. Cyclic order, component identity,
boundary face, rooted signs, phase sheets, and collision partitions are not
decorations; they are the coordinates on which the proof works.

## 1. Inheritance pass and live concept board

The required inheritance pass is:

- **closest proved mechanism:**
  [THM-4116](../01-canon/theorems/THM-4116-boundary-state-gluing-and-ap-odd-shell-tree-synchronizers.md),
  where ordered boundary extension vectors glue by an exact dot product;
- **canonical hostile:** the five-boundary Petersen carrier in
  [THM-4113](../01-canon/theorems/THM-4113-maximal-noncrossing-half-kempe-atlas.md):
  disjointness allows 15 pair-pair edges, while cyclic planarity removes a
  crossing `C5` and leaves only 10;
- **corrected near miss:** MISTAKE-507. Cubic three-edge-coloring is a
  nowhere-zero 4-flow over `(Z/2Z)^2`, not a nowhere-zero `Z/3Z`-flow;
- **least-used relevant sidecar:**
  [THM-3173](../01-canon/theorems/THM-3173-six-state-free-factor-actions-and-pointed-frame-cube.md),
  which distinguishes a genuine `C2*C3` action from its `S3`, `C6`, and
  eighteen-state fibre-product quotients.

The live concept board is:

```text
D  defect budget:       global surplus versus local exceptional tariffs;
B  boundary response:   ordered extension vectors/tensors, not totals;
K  move groupoid:       simultaneous subset switches of Kempe components;
Q  quotient atlas:      every lawful collision image, ordered by size;
G  gluing operation:    contraction/dot product over the full interface;
U  unavoidability:      a global argument forcing one local certificate. (2)
```

Every useful connection below touches at least two of these six concepts.
The portfolio is **Anchor / Niche / Wildcard**:

- **Anchor:** audit and sharpen the apex-cubic compiler itself;
- **Niche:** convert `H>=disc` for tournaments into an exact deletion charge;
- **Wildcard:** expose the literal modular/Farey action on closed triangular
  dart systems and determine precisely where digons and boundaries break it.

## 2. What the preprint claims, and what it does not

The paper claims that every 2-connected apex cubic graph is
three-edge-colorable. Combined with the Robertson--Seymour--Thomas reduction
and the published doublecross case, this claims the final apex piece of
Tutte's **cubic** Petersen-minor-free three-edge-coloring conjecture. It does
not prove Tutte's full 4-flow conjecture for arbitrary graphs.

For a smallest counterexample, the proof deletes an apex-incident edge or the
apex. The resulting planar subcubic graph has two or three degree-two
vertices. Its dual is a pseudo-triangulation with two or three digons, no
vertex incident with two digons, and controlled short separators. The proof
then combines:

1. multi-boundary islands and singleton half-Kempe chains;
2. semi-D and semi-C reducibility;
3. `915` configurations and `35` discharging-rule families;
4. `12,391` possible bad cartwheels;
5. `109,501` generated islands, all reported semi-reducible; and
6. recursive enumeration of every free homomorphic image of each
   configuration.

The authors provide detailed pseudocode and public code/data. Their own
reproducibility instructions nevertheless state that the headline metrics
are not by themselves a correctness guarantee; some standard reducer logic
comes from cited implementations, configuration order is load-bearing, and
the v1 paper pins no immutable repository commits. This session reconstructed
the interfaces and independently checked the smaller half-Kempe generator,
but did not replay the three large external obligations. See the dedicated
[source audit](../05-knowledge/reference/CORE-PAPERS-APEX-CUBIC-2026-08-25.md).

### 2.1 How much is inherited, and how much overlaps this repo

A numerical percentage would be artificial, but the dependency split is
clear.

| Layer | Provenance | Assessment |
|---|---|---|
| minimal counterexample, Kempe switching, reducible configurations, discharging, cartwheels | Four Color Theorem / doublecross lineage | inherited proof skeleton |
| darts, flexible/free homomorphisms, collision-image recursion | recent cited configuration-homomorphism machinery | substantially adapted |
| two/three primal degree-two defects as dual digons | apex case | load-bearing new input |
| singleton half-Kempe chains and multi-boundary islands | apex case | load-bearing new state space |
| three new rule families, defect allowance, dynamic refinements, 915-config universe | apex case | load-bearing finite engineering |
| exact `2--3` half-Kempe law and asymptotics | this repo, THM-4113 | new independent sharpening |
| ordered response-vector gluing and AP odd-shell tree law | this repo, THM-4116 | new independent theorem |
| tournament deletion charge | this repo, HYP-9080 | exact proposal, global step open |

So the global strategy is recognizably inherited, while nearly every place
where the apex destroys planarity's clean cubic interface requires new
machinery. The direct overlap with prior repo canon is concentrated at the
**interfaces**—matching carriers, response vectors, component quotients,
torsion sheets, and descent—not in a pre-existing theorem that implies the
paper's coloring claim. No repo result bypasses the external computer checks.

## 3. The sharp four-digon wall

Let a closed plane pseudo-triangulation have `q` digon faces and all other
faces triangular. If it has `V` vertices and `E` edges, then

```text
2E=2q+3(F-q),             V-E+F=2,
E=3V-6+q.                                                   (3)
```

For the paper's initial charge `T_0(v)=10(6-d(v))`, equation `(3)` gives the
exact total

```text
sum_v T_0(v)=60V-20E=120-20q.                              (4)
```

The paper's five one-digon neighborhood cases have charge ceilings

```text
18, 12, 16, 8, 12,                                        (5)
```

so `18` is the exact maximum of the proved bounds. The current scalar
argument has enough surplus exactly when

```text
18q < 120-20q  iff  q<=3.                                  (6)
```

This is a genuine method boundary. At `q=4`, the global surplus is `40` but
four independent allowances may total `72`. A four-digon extension using
only the same additive budget would need a uniform tariff `<10`, or else a
new negative interaction/cluster term proving that the four worst local
cases cannot coexist. Merely enlarging the configuration list does not alter
this arithmetic unless it changes the local tariff or the global charge.

This suggests the correct next discharging search:

```text
source:       four marked digon neighborhoods;
target:       a cluster-interaction charge, not four independent charges;
preserve:     Euler total and every legal rule transfer;
destroyed by scalarization:
              overlap, separation, and incompatible worst cases;
needed sidecar:
              digon incidence graph plus distance-two overlap type;
cheapest test:
              optimize the exact maximum total tariff over all legal
              four-digon radius-two clusters.                         (7)
```

## 4. The boundary reducer is a greatest-fixed-point safety game

The paper's semi-consistency definition has a clean finite fixed-point form.
This is an exact reformulation, not a new assumption.

Let `D_R` be the ternary color words on an ordered multi-boundary `R`. For a
word `phi` and an unordered color pair `p={x,y}`, let `A(phi,p)` be the
admissible planar semi-matchings which partition the positions colored `x`
or `y`. For `M in A(phi,p)`, let

```text
Orb(phi,M)={all words obtained by switching x,y on any subset M' of M}. (8)
```

For a bad-state set `B` and `S subseteq B`, define

```text
Phi_B(S)={phi in S:
  for every color pair p there exists M in A(phi,p)
  with Orb(phi,M) subseteq S}.                              (9)
```

`Phi_B` is monotone and deflationary. Starting from `S_0=B`, set

```text
S_(i+1)=Phi_B(S_i).                                        (10)
```

Because `B` is finite, `(10)` stabilizes. Its limit `nu Phi_B` is exactly
the largest semi-consistent subset of `B`:

- a fixed point satisfies the definition of semi-consistency;
- any semi-consistent `C subseteq B` lies in every `S_i` by induction and
  hence in the limit.

Thus semi-D-reducibility is exactly `nu Phi_B=empty`. Semi-C-reducibility by
a deletable set `F` is exactly

```text
(nu Phi_B) intersect C_(I dotdiv F)=empty.                 (11)
```

The quantifier order in `(9)` matters. The adversary selects a color pair;
the exterior chooses one admissible semi-matching; then every subset of its
components must remain safe. This is a finite safety game with a
hyperedge-valued response, not a pairwise tournament. Orienting which of two
boundary words is “better” destroys the simultaneous-subset operation.

This fixed-point form is portable. For another problem one needs a finite
bad-state carrier, its native moves, and a proved lift of every move back to
the original object. A pruning algorithm without the lift is only a search
heuristic.

## 5. THM-4113: the half-Kempe topology atlas is algebraic

[THM-4113](../01-canon/theorems/THM-4113-maximal-noncrossing-half-kempe-atlas.md)
proves that Appendix B.1.2 produces exactly the maximal noncrossing partial
matchings on a cyclic boundary. An unmatched point is a singleton
half-Kempe chain. Rooting a boundary gap gives a bijection with either

1. one noncrossing partition into blocks of sizes two and three; or
2. one root singleton plus an ordered pair of such partitions.

If `u` marks singletons, the bivariate generating functions are

```text
C=1+x^2 C^2+u x^3 C^3,
M=C+uxC^2.                                                (12)
```

Consequences proved and exactly audited in the theorem include:

```text
|H_n| = 1,1,1,3,4,10,20,42,98,210,492,1122,2607,...,
x^4M^3-2x^2M^2+(1+x^2)M-(1+x)=0,
|H_n|=(-1)^(n+3) A321197(n+2).                            (13)
```

The implicit schema also gives a proved asymptotic

```text
|H_n| ~ 1.295313565966669... * 2.610718613276039...^n
                                  * n^(-3/2).              (13a)
```

The growth base is the positive root selected by
`4alpha^3+alpha^2-18alpha-31=0`. Thus the new generator removes redundant
search but the topology output itself remains exponentially large.

For boundary sizes `n_1,...,n_r`, the topology count is the literal product
`product_i |H_(n_i)|`; boundary count is therefore a state-space coordinate.

At five boundary points, the coarse compatibility graph on endpoint pairs is
the Petersen graph `K(5,2)`, equivalently THM-261's `A4` positive-root
orthogonality graph. It has 15 disjoint-pair edges. Cyclic planarity removes
the five crossing diagonal-diagonal edges, a `C5`, leaving the 10 atlas
states. This is the cheapest exact firewall against two recurring errors:

- unrestricted perfect matching/hafnian logic is not planar half-Kempe
  logic; and
- a compatibility graph, even a distinguished Petersen graph, does not
  recover cyclic order.

The direct `2--3`-partition generator can replace the paper's
generate-and-prune topology routine. It does not count ternary boundary
colorings or compute the fixed point `(10)`.

There is a further algebra firewall. After rooting, these diagrams sit inside
a Motzkin link-state basis, but maximal states are not a module for the usual
local Temperley--Lieb/Motzkin reconnection action. Already at four points,
the maximal state with chord `1--3` and singletons `2,4` is sent by the first
local reconnection to chord `1--2` with adjacent singletons `3,4`, which is
nonmaximal. The correct algebra is state-dependent: a chosen semi-matching
gives one Boolean cube of subset switches, and semi-consistency takes an
AND--OR greatest fixed point over those cubes. No Temperley--Lieb spectral
shortcut transfers without a new invariant-subspace theorem.

## 6. THM-4116: the actual snark interface is a response-vector pairing

[THM-4116](../01-canon/theorems/THM-4116-boundary-state-gluing-and-ap-odd-shell-tree-synchronizers.md)
proves the exact gluing law independently of the preprint. If a graph is cut
along ordered edges and `f_X(sigma)` and `f_Y(sigma)` count proper labelled
three-edge-coloring extensions of the two sides with boundary word `sigma`,
then

```text
#Col_3(G)=sum_sigma f_X(sigma)f_Y(sigma).                   (14)
```

The Petersen two-vertex cut has disjoint support vectors, so the dot product
is zero. The analogous `K4` cut has dot product six, one coloring modulo the
free global `S3` color gauge. This is a literal obstruction certificate,
unlike a scalar “snarkness” score.

The same operation, applied to the already-proved AP odd-shell phase quotient
of THM-4110, gives an exact incoming result. If `F` is the graph of added
odd--odd phase constraints on the `q` odd speeds of `AP_(2q-1)`, then

```text
|K_(Gamma_F)(v)/P_v|=2^(c(F)-1).                          (15)
```

The phase quotient is saturated exactly when `F` is connected. Minimal
saturators are the labelled spanning trees, so there are `q^(q-2)` of them.
For AP13 this explains why six extra constraints are necessary but why an
arbitrary six are insufficient: a cycle can waste one constraint while an
odd speed remains isolated.

Equation `(15)` is the session's cleanest cross-frontier success because the
source and target share an actual operation—interface constraint gluing.
Even here it proves only phase synchronization. It supplies no safe time and
does not prove LRC(14).

## 7. Collision-complete quotient descent

The paper's `allHomImages` routine handles a configuration whose darts may be
identified in the target. It does not test one preferred quotient. Whenever
a map is noninjective, the collision factors through a free homomorphic image
with strictly fewer darts, and the recursion explores every resulting image.
The dart count is therefore a well-founded descent measure.

The abstract compiler is:

```text
input object K
  -> test injective realizations of K
  -> enumerate every lawful first collision relation
  -> form each quotient K/~ while preserving typed incidences
  -> recurse on strictly smaller quotients
  -> test the target predicate on every terminal image.              (16)
```

This is a factorization-complete quotient atlas. It preserves dart head,
edge reversal, non-null successor/predecessor incidences, degree data, and
the relevant rotation system. An adjacency-only graph quotient is too coarse:
it forgets face rotation and boundary nil pointers.

The closest proved repo analogues are:

- [THM-3042](../01-canon/theorems/THM-3042-subdirect-graph-order-common-quotient-and-singleton-owner-criterion.md):
  two projected structures glue through one explicit common quotient;
- [THM-4067](../01-canon/theorems/THM-4067-seminormal-period-kernel-and-figure-eight-completeness-obstruction.md):
  cycle periods are complete exactly when the ambient realization equals the
  full graph equalizer; and
- THM-4110/4116: the finite quotient also needs its Smith/component sidecar.

The shared lesson is not “quotient everything.” It is “enumerate every
allowed quotient, retain the target predicate's missing coordinates, and
prove descent.”

## 8. A literal modular action, with three stopping boundaries

On a closed dart system let

```text
alpha=rev,           sigma=succ,           F=alpha sigma.   (17)
```

`alpha` reverses an edge and `F` walks a face. On a closed triangle-only
triangulation,

```text
alpha^2=1,           F^3=1,                               (18)
```

so the darts carry a permutation action of

```text
C2*C3 ~= PSL_2(Z).                                        (19)
```

The action need not be faithful. Vertex rotations are the cycles of
`sigma=alpha F`, analogous to cusp-width cycles. This is the first genuinely
literal bridge from the paper to the repo's Farey/Stern--Brocot work: the
binary edge reversal and ternary face rotation are the two modular free
factors.

There are three sharp stopping boundaries:

1. With both digons and triangles, `F` has cycles of lengths two and three;
   globally only `F^6=1` is automatic. The coarse action is through
   `C2*C6`, with a mandatory face-type `2/3` sidecar.
2. Triangle darts alone need not be stable under reversal, so deleting the
   digon darts does not automatically leave a modular subaction.
3. In outer extensions, `succ` or `pred` can be nil and the face axiom may
   fail. This is a partial transformation system until a closed completion
   is chosen.

[THM-3173](../01-canon/theorems/THM-3173-six-state-free-factor-actions-and-pointed-frame-cube.md)
is the exact guardrail: a six-state `S3` or `C6` quotient does not reconstruct
Farey ancestry or the free-product action. The snark darts give a new literal
carrier, but no map from their modular orbits to an LRC time, a Jacobian
cofactor, or an arithmetic owner has been constructed.

## 9. Exact tournament proposal: deletion-slack unavoidability

The apex proof separates **local reducibility** from **global
unavoidability**. [HYP-9080](../05-knowledge/hypotheses/HYP-9080-tournament-deletion-slack-local-unavoidability.md)
turns this architecture into an exact native tournament identity.

For a tournament `C`, vertex `x`, `T=C-x`, and slack

```text
S(W)=H(W)-disc(W),                                         (20)
```

let `u_x` be the actual rooted incident sign vector of `x`, and let
`Delta H_x=H(C)-H(T)`. THM-3729 and THM-4094 give

```text
chi_x(C)
 =2 Delta H_x-E_odd(T;u_x)+disc(T)
 =2[S(C)-S(C-x)].                                          (21)
```

The open local-unavoidability conjecture is

```text
for every C, some x has chi_x(C)>=0.                       (22)
```

It would prove `H>=disc` by deletion induction. The averaged strengthening is

```text
1/2 sum_x chi_x(C)
 =|C|S(C)-sum_x S(C-x)>=0.                                 (23)
```

Both `(22)` and `(23)` are **OPEN**. Exact computation found no all-negative
row among all labelled tournaments through order six or all 456 order-seven
isomorphism classes. Pointwise monotonicity is false: a source over `C3` has
one charge `-2`, although its other three deletions are positive. This is the
minimal hostile and explains why the desired statement must be existential.

OCF and the principal-minor expansion sharpen `(23)` exactly. If `Gamma`
ranges over vertex-disjoint directed odd-cycle collections with support
`U(Gamma)`, and `p_A=Pf(K[A])^2`, then

```text
|C|S(C)-sum_xS(C-x)
 =sum_Gamma |U(Gamma)|2^|Gamma|
  -2^(1-|C|)sum_(A even)(2|A|-|C|)p_A.                    (23a)
```

Only even supports larger than half the vertices are adverse. Positivity of
the OCF atoms does not compare the two support families, so `(23a)` reframes
ALU without proving it.

There is also an exact ear tariff. For base `T`, let `d=disc(T)`,
`R=(I-K^2)^(-1)`, take THM-4104's cut weight `w` and zero-sum field `h`, and
put `wtilde=w+2dR`, `a_0=disc(T+sink)-d`. Then

```text
chi_(x_S)/2=cut_wtilde(S)+h(S)-a_0,
chi_(x_(V-S))/2=cut_wtilde(S)-h(S)-a_0.                    (23b)
```

The open strong tariff is `cut_wtilde(S)>=a_0+|h(S)|` for every nonconstant
cut. A new exact companion checked all 2,031,616 nonconstant ears over every
labelled order-six base and all 57,456 ears over the inherited 456-class
order-seven bank: minimum `chi=0`, no negative case; the minimum combined
edge weight remains `1/2`. Constants over `C3` retain `chi=-2`, while all six
mixed ears have `chi=4`. These are finite facts, not an all-order proof.

Incoming THM-4115 supplies the exact Walsh transform of this carrier. If
`Y_S=chi_(x_S)/2` and `Wtilde=sum_(i<j)wtilde_ij`, then

```text
E Y=Wtilde/2-a_0,
Var(Y)=1/4(sum_i h_i^2+sum_(i<j)wtilde_ij^2).              (23c)
```

This restores the full quadratic coefficient norm, but not the labelled
support floor: `min(Y_S,Y_(V-S))=cut_wtilde(S)-a_0-|h(S)|`.
The all-order bound `Y_S>=-d(n-1)/2` does give a support-sensitive maximum
inequality, sharp on the `C3` base, but it selects a favorable extension
orientation. LU and ALU prescribe charges coming from different deletion
bases. Thus mean/variance growth and local unavoidability remain distinct
problems.

THM-4118 adds the right kind of sidecar: the exact singleton/pair response
lattice plus labelled flip/exchange components whose unit-step images are
solid intervals. Its proof transfers to the combined integer quadratic `Y`,
with a new gcd `delta_Y`; its Hamiltonian conclusion `d_T=2` does not. This
organizes tariff hostiles by actual cut-state reachability instead of a
scalar moment, but still leaves the sign of each component minimum open.

The typed comparison is:

```text
snark local certificate:       reducible boundary island;
tournament local certificate:  chi_x(C)>=0;
snark global forcing:           Euler discharging;
tournament global forcing:      LU or averaged LU, still missing;
mandatory tournament sidecar:   u_x plus full insertion/orphan response;
forbidden scalarization:        H, disc, or a selected matching alone.    (24)
```

THM-4099's mixed insertion polynomial and THM-4104's full `Start/End/Q` ear
response are the promising expansion coordinates. THM-4111 is the firewall:
equal aggregate ear means can hide different response images, so mean growth
cannot replace unavoidability.

## 10. Overlap ledger with earlier and incoming work

| Repo object | Exact overlap | Information that still fails | Status |
|---|---|---|---|
| THM-261 | Petersen graph is the coarse five-boundary pair carrier | cyclic crossing `C5`, face, colors | PROVED |
| THM-2290 | endpoints and contraction family are load-bearing | hafnian includes crossing/unrestricted pairings | PROVED |
| THM-857 | exact component recursion and certified pruning | LRC component is not a Kempe component | PROVED / finite |
| THM-3990 | repair must be componentwise | no Laplacian map to Kempe chains | PROVED |
| THM-4009 | a counterexample has bounded relation defects | no boundary response or physical arrival | CITED + PROVED algebra |
| THM-4099 | full mixed insertion boundary is compositional | proper faces lose mixed coefficients | PROVED |
| THM-4104 | complete one-ear response is a finite cut tensor | selected image is not a global census | PROVED / finite |
| THM-4111 | exact zero-mode/mean formula | mean loses image and reducibility | PROVED |
| THM-4112 | parent-gap component-span envelopes iterate on scale-separated LRC families | exact component identities, owners, physical origin, and arbitrary cores remain | PROVED under explicit hypotheses; relative to THM-2061/2066/2072 |
| THM-4114/4115/4118 | exact ear quadratic, moments, response lattice, and unit-state intervals | combined-charge component minima and cross-family domination remain | PROVED / independently audited |
| THM-4117/4119 | a physical row and an infinite lonely residue family miss all three known supplier shapes | full primitive support and marked unit/core origin are lost | PROVED; physical scope only for THM-4117 member |
| THM-4110/4116 | exact phase quotient and component repair law | synchronized phase is not LRC safety | PROVED |
| THM-3042/4067 | common quotient/equalizer controls gluing | ambient realization and torsion can survive | PROVED |
| THM-3173 | exact `C2*C3` quotient firewalls | cardinality six is not a modular action | PROVED |
| THM-4085 | injective characteristic addresses give a triangular bijection | repeated addresses have response `(1,3,3,1)` | PROVED |
| THM-2230 | a full Jacobian response has an exact target-shear fibre | response quotient does not imply affine owner | PROVED |

Two incoming promotions materially sharpen the comparison. THM-4116 turns the
older AP13 “six missing constraints” observation into the complete component
law `(15)`. THM-4112 proves an arbitrary-depth parent-gap span-envelope
recursion under its separation hypotheses, including ratio-two depth through
six, every finite adjacent ratio at least `12/5`, and explicit eleven- and
thirteen-speed families. It does not recover exact component identities,
multiplicities, owners, or physical origin. Neither promotion covers arbitrary
LRC(14) cores. THM-4114 is now
**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED**; its OCF presence-cube
positivity and opposite-ear directed-cut curvature are genuine inputs, but
they do not compare the OCF support moment with the adverse Pfaffian support
moment in `(11c)`, so they still do not prove ALU or HYP-9080.

THM-4117 is the mandatory physical-origin hostile. Its physical row

```text
(1,4,6,8,10,12,14,15,16,18,22,2^45,3*2^45)
```

is `1/14`-lonely with clearance `2/19`, yet no common dilation or primitive
normalization enters THM-4112's `AP7+4`, `AP8+5`, or `D0+6` supplier classes;
the forced AP7 quotient first misses the marked unit `1`. Rank, `11+2` type,
pair type `(1,3)`, mod-56 residues, parity, and phase-sheet data therefore do
not determine supplier entry. Failure of a sufficient THM-4112 certificate
is neither a cover nor an irreducible boundary state.

THM-4119 proves this is not an isolated arithmetic accident. For

```text
S_t={1,4,6,8,10,12,14,15,16,18,22,t,3t},       t>22,
```

fourteen of the nineteen residues of `t` give the exact common-phase
clearance `2/19`, while every `t>22` misses all three supplier shapes under
common dilation and primitive normalization. The canonical `t=2^45` member
is the THM-4117 physical row; the other parameters are not claimed to come
from the finite physical producer. This separates two issues that a boundary
atlas must retain: arithmetic loneliness and physical-origin realizability.

## 11. LRC(14): a concrete boundary-atlas programme

The paper suggests a sharper LRC compiler than another scalar sieve.

```text
source:       a bounded THM-4009 short-relation packet plus circular
              strict-safe components;
target:       a finite interface state indexed by component, owner, clock,
              residue, and phase sheet;
map:          restrict a row to a separator and record its exact extension
              response under every lawful alias/collision quotient;
preserved:    nonempty circular component and physical time compatibility;
destroyed by current coarse maps:
              literal primitive support, marked unit/core origin, exact
              component identity, endpoint owner, phase, and Smith torsion;
needed sidecar:
              THM-4100 ancestry + THM-4110/4116 sheet component + exact
              boundary response rather than a mean;
cheapest decisive tests:
              the THM-4100 equal-shadow/different-component hostile and
              AP13's 64 phase sheets.                                  (25)
```

The analogue of `allHomImages` is to enumerate every coordinate alias forced
by a short relation, factor each noninjective state through a smaller
interface, and retain the full circular response. THM-4112 supplies a real
positive control for the span-envelope coordinate: on its scale-separated
families the parent-gap recursion can be iterated to arbitrary prescribed
depth, and the proved ratio-two and `12/5` regimes show that this is not merely
a finite toy. The analogue of reducibility is that every surviving exterior state
extends to a `1/14`-safe physical cell. The analogue of unavoidability is still
absent for arbitrary cores: bounded relations do not yet force a reducible
interface type, and the scale-separation hypotheses cannot be silently erased.

THM-4116 closes the phase-synchronizer subproblem on odd AP shells: connect
the odd-shell constraint graph. THM-4112 closes the parent-gap span-envelope
recursion under its declared separation hypotheses. Their conjunction still
does not close exact component identity, physical origin, owner arrival,
alias-complete response, or safe-cell intersection. THM-4117 proves this is a
real obstruction even for a lonely physical row. The snark quotient is sound
because it retains exact extension states and quantifies over every admissible
colored semi-matching—not because topology reconstructs exterior origin. An
LRC atlas must likewise retain every physical-origin fibre's full response;
a catalogue of known suppliers is not an unavoidable atlas.

## 12. Rule 30: collision quotients after characteristic injectivity

THM-4085 proves that distinct left-characteristic addresses give a triangular
bijection and hence independent Haar outputs. The first repeated-address
hostile has joint counts

```text
(00,01,10,11)=(1,3,3,1).                                  (26)
```

This is exactly the point at which an injective chart must be replaced by a
collision atlas. The snark-inspired task is:

```text
source:       a finite set of spacetime observations;
quotient:     the partition identifying equal characteristic addresses;
local data:   within-block Boolean response, beginning with `(26)`;
operation:    merge collision blocks and factor through fewer addresses;
target:       entropy/correlation or distinguished-seed cylinder predicate;
sidecar:      right-tail variables and chronological order.              (27)
```

This is an architectural transfer, not a Kempe action. The Rule 30 local map
is irreversible and directed in time, while Kempe switches are involutive.
The named-seed prizes also cannot inherit Haar conclusions.

## 13. Planar Jacobian conjecture: quotient every collision, retain every branch

For planar JC the promising analogy is not discharging. It is the
factorization-complete collision atlas `(16)` combined with THM-2230's exact
response quotient.

```text
source:       labelled generic inverse branches / punctured target fibre;
target:       target-response classes modulo the proved shear kernel;
candidate quotient:
              every collision partition of labelled branch or puncture data;
must preserve:
              multiplicity, ramification, monodromy action, target value,
              Jacobian response, and affine-source regularity;
destroyed by degree/genus/discriminant alone:
              branch owner and which sheets collide or escape;
cheapest test:
              enumerate the first seven-puncture collision partitions and
              compare response fibres before and after every fold.       (28)
```

The descent measure would be the number of labelled branches/punctures, just
as dart count decreases in `allHomImages`. But a graph quotient cannot be
substituted for a branched-cover quotient. Ramification multiplicity and
escape at infinity are mandatory; MISTAKE-319 already records that a drop in
finite fibre count can be sheet escape rather than collision. No planar JC
case follows from the current analogy.

## 14. Arithmetic, ABC/IUT, Mahler `3/2`, and torsion: useful negatives

Several earlier arithmetic analogies become clearer after the exact defect
budget and flow-group correction.

1. **ABC/IUT.** “Height is global surplus and radical is local defect” is at
   most a heuristic schema. ABC supplies no planar overlap bound or finite
   local configuration atlas. It may enter only through an actual primitive
   three-term relation already produced by the LRC relation compiler; it
   cannot replace `(7)` or `(25)`.
2. **Mahler `3/2`.** A finite symbolic chart is useful only if carries,
   valuation depth, and termination lift through every state quotient. The
   snark compiler suggests enumerating collided address/fibre quotients, but
   no finite atlas or Z-number decision is obtained.
3. **Binary versus ternary torsion.** Kempe changes for a selected color pair
   are exponent-two switches, and cubic edge-color flows use `(Z/2)^2`.
   They cannot additively encode THM-4110's independent `Z/3` or higher Smith
   sheets. A shared number `3` is not a homomorphism.
4. **Hafnian/Pfaffian.** Boundary planarity produces Motzkin/Temperley--Lieb
   style noncrossing states, while THM-2290's endpoint-selected kernel sums
   over unrestricted perfect matchings. The missing crossing monomial occurs
   already on four cyclic endpoints.

These negatives are progress: they identify which sidecar a real bridge
would have to construct.

## 15. Generated high-value tasks

### Anchor

1. Archive immutable commits of all four external verifier repositories and
   reproduce the `747`, `12,391`, and `109,501` ledgers with an implementation
   independent of both the authors' and the AI-generated code.
2. Replace Appendix B.1.2 by THM-4113's direct `2--3`-partition generator,
   then verify that the full colored semi-consistency fixed point is unchanged.
3. Solve the four-digon optimization `(7)`: either push the realizable tariff
   below `10`, find a negative cluster interaction, or freeze a legal hostile
   showing that this discharging scheme cannot extend.
4. Build a face-labelled response tensor whose indices retain boundary
   component, cyclic order, singleton profile, and color pair; test exact
   contraction under island gluing.

### Niche

5. Enumerate order-eight tournament isomorphism classes for HYP-9080 and
   seek a double-counting proof of `(23)` from the full `Start/End/Q`, orphan,
   and rooted odd-Pfaffian responses.
6. Expand `chi_x` into local ear-cut atoms. Classify equality and first
   negative summands on strong tournaments before attempting a universal
   charge rule.
7. For the even-graph/five-cycle lane, audit the proposed weighted deletion
   identity

   ```text
   sum_v(F_n(H)-F_(n-1)(H-v))
     =sum_(k=3)^6 k c_k^-(H)-n(A_n-A_(n-1)),                (29)
   ```

   with `A_n=n^4-13n^3+63n^2-132n+100`; decide whether its
   finite-positive behavior through `n=8` has a structural double count.
   This is an **OPEN proposal**, not canon.

### Wildcard

8. On the LRC `11+2` packets, build the collision-complete alias atlas in
   `(25)` and retain exact component response, odd-shell tree connectivity,
   owners, clocks, and Smith sheets.
9. Implement the Rule 30 repeated-address quotient atlas `(27)` through the
   first nontrivial collision partitions and measure which within-block
   responses determine entropy beyond the distinct-address floor.
10. Implement the planar-JC seven-puncture collision atlas `(28)` with
    multiplicity and infinity-escape sidecars; stop immediately if two folds
    have the same quotient graph but different target-response fibres.
11. Classify closed invariant dart subsets of the mixed digon/triangle maps.
    Determine when the coarse `C2*C6` action contains a genuine
    `C2*C3` subaction and record its Farey cusp-width data.

## 16. Bottom line

The preprint's main claim is substantial but externally computer-assisted and
not yet internally replayed. The session nevertheless produced three durable
advances:

1. **PROVED:** the half-Kempe topology generator is the maximal noncrossing
   partial-matching class with exact `2--3` partition law (THM-4113);
2. **PROVED incoming:** boundary extension vectors glue by a dot product, and
   the AP odd-shell phase repair is exactly graph connectivity (THM-4116);
3. **OPEN but exact:** tournament `H>=disc` admits the deletion charge
   `chi_x=2(S(C)-S(C-x))` and a precise unavoidability conjecture (HYP-9080).

The most surprising literal connection is the closed-dart
`C2*C3 ~= PSL_2(Z)` action. The most useful general connection is the
greatest-fixed-point boundary game plus collision-complete quotient descent.
The sharpest obstruction is the four-digon wall: the current local tariff is
not remotely strong enough beyond three defects. Across LRC14, Rule 30, and
planar JC, the actionable proposal is therefore the same but fully typed:

> keep the complete boundary response, enumerate every lawful collision
> quotient, attach the sidecar lost by the quotient, and prove a global
> argument that forces one locally reducible state.

Nothing in this synthesis proves LRC(14), planar JC, a Rule 30 prize, Mahler's
`3/2` problem, ABC, or the full Tutte 4-flow conjecture.
