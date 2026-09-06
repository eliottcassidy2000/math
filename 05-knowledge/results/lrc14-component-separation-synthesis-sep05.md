# Sharp pair arithmetic and the scale-three `44/13` tail closure

**Status: PROVED ANALYTICALLY + FINITE-EXACT arithmetic obligations +
INDEPENDENTLY AUDITED.** The independent `nc2_seed` referee and root reviewed
the all-height proof, its quantifiers, clipping, residue universes, and strict
tail boundary. No scarce theorem identifier is used in this research note.

For positive distinct odd integers `a<b<c`, each nonzero modulo three, the
degree-zero certificate of
[THM-4409 — third-sheet component-network certificate](../../01-canon/theorems/THM-4409-lrc14-third-sheet-component-network-certificate.md)
has the following all-height properties:

1. Every contact graph is a star forest, and its maximum flow is the sum of
   the edgewise minimum component lengths.
2. The six ordered distinct-sheet pair intersections have total mass
   `T_ab<=12/77` and total component count `N_ab<=24b/13`. Both constants
   are sharp, at reduced pairs `(1,11)` and `(11,13)`, respectively.
3. If `c/b>=44/13`, the fixed pair `(a,b)` already gives a strict certificate
   `U_net(a,b;0,0)<6/77`.

The first assertion independently recovers a result published concurrently in
[Sparse interval transport removes the all-height flow obstruction](synthesis_20260905_lrc_sparse_transport.md).
It is included below as a self-contained prerequisite, not claimed as a
distinct new contribution. The sharp pair arithmetic and the third assertion
extend that report: an infinite scale-separated region is now closed. The
remaining comparable speeds, arbitrary chart entry, and synchronization stay
open. **LRC(14) remains OPEN.**

## Inheritance and connection ledger

The closest proved mechanism is THM-4409's component-capacity relaxation.
The canonical equality hostile is `(1,5,11)`; the information-loss hostile is
`(1,19,79)`, whose physical mass is `108/10507` while its best degree-zero
capacity is `8/553`. The corrected near miss is
[THM-4396 — finite dual and exact-pair hybrid certificate](../../01-canon/theorems/THM-4396-lrc14-finite-dual-exact-pair-hybrid-certificate.md):
dropping the third incidence prevents equality at every finite pair degree.
The relevant sidecar is the separation *between* components, inherited from
the fastest tooth in a pair. Its role is distinct from the nested-owner
Hunter star in
[THM-1283 — terminal endpoint transfer and gcd tax](../../01-canon/theorems/THM-1283-terminal-endpoint-transfer-and-gcd-tax.md).

The live concept board is: contact capacities; tooth separation; shifted
Bernoulli pair arithmetic; component count; equality versus clipped crossings.
The map is from labelled physical interval intersections to a bipartite graph
with component lengths. It preserves an upper bound on physical triple mass,
but loses overlap endpoints. Separation makes the optimization exact for that
quotient, without making the quotient exact for physical mass. The hostile
`(1,19,79)` continues to distinguish those two statements.

The final freshness pass found concurrent commit `566677ae1d` and its
[sparse-transport report](synthesis_20260905_lrc_sparse_transport.md). It
already proves the separated-star identity, its sharp general gap boundary,
and an exact raw-roof expression for each network capacity. These overlap
with the first part of this work and are credited as concurrent independent
work. Its reserved namespace `THM-4414` is an empty stub and is not used as a
proved dependency. The report explicitly leaves its universal inequality
`min_ij sum_C K_ij(C)<=6/77` open; the pair bounds and `44/13` region below
are not present there. Searches for the exact sharp pair statements and
`44/13` found no earlier match. This is a repository comparison, not a
literature priority claim.

## 1. A separated-interval lemma

Suppose disjoint interval families `A,B` have respective maximum lengths
`L_A,L_B`, and distinct intervals within each family are separated by gaps
at least `2L_A,2L_B`, respectively. Join two intervals when they overlap in
positive length, and put length capacities on vertices.

Assume `L_A<=L_B`. Every `A` interval meets at most one `B` interval, since
its length is smaller than the gap between two `B` intervals. Thus every
nontrivial connected component is a star centered at a `B` interval.
If a center touches `n>=2` leaves, it spans the `n-1` intervening gaps;
positive contact gives center length strictly greater than `2(n-1)L_A`.
The sum of the leaf lengths is at most `n L_A<=2(n-1)L_A`, so assigning
each edge its entire leaf length is feasible. For a one-edge star assign the
minimum endpoint capacity. Every individual edge load is bounded above by
this minimum, so the assignment simultaneously attains all edgewise bounds:

```text
Cap_G(1) = sum_(I,J in E(G)) min(|I|,|J|).              (1)
```

The argument allows end pieces clipped by a fixed interval boundary.
It does not hold for arbitrary nonconstant weights: length separation need
not imply separation of weighted capacities.

For a danger comb of radius `0<lambda<=1/6`, a speed-`v` tooth has length
`L_v=2lambda/v` and adjacent gaps `(1-2lambda)/v>=2L_v`.
For a pair with largest speed `b`, each `b` tooth intersects at most one
tooth of the other comb. Hence the pair components have length at most
`L_b` and inherit gaps at least `2L_b`. Apply the lemma to this pair family
and the third comb. Therefore (1) is valid for arbitrary integer speeds,
arbitrary shifts, and every radius `lambda<=1/6`; no parity, primitivity,
or mod-three condition is needed for this structural assertion. At the LRC
radius `1/14`, the available gap is in fact `6L_v`.

The gap factor two is sharp for this abstract lemma. For any `0<=gamma<2`,
take unit leaves `[0,1]` and `[1+gamma,2+gamma]` and center
`[1-epsilon,1+gamma+epsilon]`, where
`0<epsilon<min(1,(2-gamma)/2)`. The center length `M=gamma+2epsilon` is
strictly between zero and two. Both contacts have positive length, and the
edge-minimum sum `2 min(1,M)` exceeds the center capacity `M`. In particular
gap `3/2` is refuted by center `[9/10,13/5]`, whose capacity is `17/10`.
These are abstract interval hostiles, not counterexamples in the fixed
three-sheet LRC family.

## 2. Exact ordered-pair count

Now fix radius `1/14` and the LRC sheets

```text
D_(v,s)={x in R/Z : ||v(x+s/3)||<1/14}, s in F_3.
g=gcd(a,b), p=a/g, q=b/g, 1<=p<q.
```

The six pairs have distinct sheets. Their components avoid the cut `x=0`:
at zero, only sheet zero of a ternary-unit speed is dangerous, so two distinct
assigned sheets exclude an entire neighborhood of zero. In particular the
usual `[0,1]` implementation does not split a pair component. A third-sheet
tooth at zero is still clipped into two intervals, exactly as in THM-4409;
its table has `68` third pieces rather than `6*11` for pair `(1,5)`.

For reduced coprime `p,q`, the differences of tooth centers run once through
the lattice of spacing `1/(pq)` on the circle. A distinct-sheet shift has
relative phase `pq/3`, whose fractional part is `1/3` or `2/3`. Reflection
identifies those two cases. Two teeth meet exactly when their relative center
distance is less than the sum of their radii. Therefore each ordered pair has
exactly `g n` components, where

```text
r=(p+q)/14,
n=#{z in Z : |z+1/3|<r}
 =floor(r+1/3)+floor(r-1/3)+1,
N_ab=6g n.                                             (2)
```

There is no endpoint equality: `r` has denominator dividing seven since
`p+q` is even, and `r+/-1/3` is not an integer.

Writing `r=t/7`, the seven values of `n-2r` modulo `t->t+7` are

```text
0, -2/7, -4/7, 1/7, -1/7, 4/7, 2/7.                  (3)
```

Thus `n<=(p+q+4)/7`. Since distinct odd `p,q` obey `p<=q-2`, for `q>=13`

```text
N_ab/b=6n/q <=12/7+12/(7q) <=24/13.                   (4)
```

The only eligible reduced denominators below thirteen are `q=5,7,11`.
Their largest `6n/q` values are respectively `6/5,12/7,12/11`, all strictly
smaller. Equality in (4) forces and is attained by `(p,q)=(11,13)`.

## 3. Sharp ordered-pair mass

The intersection of two intervals of radii `1/(14p),1/(14q)` and center
distance `|z+1/3|/(pq)` has length

```text
(1/(pq)) [min(p/7, (p+q)/14-|z+1/3|)]_+.
```

Dilation `y=gx` sends sheet label `s` to `g*s mod3`. Because `g` is a
ternary unit, this permutes the six ordered distinct-sheet assignments.
It repeats the pattern `g` times and divides each length by `g`.
Consequently the total mass of the six pair intersections is exactly

```text
T_ab=(6/(pq)) sum_(z in Z)
       [min(p/7, (p+q)/14-|z+1/3|)]_+.                 (5)
```

The old piecewise-quadratic/Bernoulli representation of periodic interval
overlap becomes useful after retaining the distinct-sheet phase `1/3`.
Put `B_2(x)={x}^2-{x}+1/6`. The elementary periodized-tent identity is

```text
sum_z (u-|z+theta|)_+
 =u^2+B_2(theta)
  -(B_2(theta+u)+B_2(theta-u))/2.                       (6)
```

It follows directly by comparing the piecewise-linear derivatives in `u`
and the value zero at `u=0`. Subtract (6) at
`u=(p+q)/14` and `u=(q-p)/14`; their squared difference is `pq/49`.
Equation (5) becomes

```text
T_ab=6/49+d(p,q)/(pq),
d=3[B_2(1/3+(p-q)/14)+B_2(1/3-(p-q)/14)
    -B_2(1/3+(p+q)/14)-B_2(1/3-(p+q)/14)].             (7)
```

The 49 odd residue pairs `(p,q) mod14` give `d<=26/49`. This is a finite
exact rational calculation; the maximum occurs at residues `(3,11)` and
`(11,3)`. Thus, if `pq>=16`,

```text
T_ab <=6/49+26/(49*16) <12/77.                         (8)
```

The entire remaining universe consists of the four pairs

| reduced pair | exact total pair mass |
|---|---:|
| `(1,5)` | `4/35` |
| `(1,7)` | `6/49` |
| `(1,11)` | `12/77` |
| `(1,13)` | `12/91` |

This proves `T_ab<=12/77`, with equality exactly for reduced pair `(1,11)`.

## 4. The infinite `c/b>=44/13` closure

Let `U_ab` be THM-4409's degree-zero certificate selecting pair `(a,b)`.
By (1), it is exactly the sum of all contact-edge minimum lengths.

If `b<c<6b`, a pair component has length at most `1/(7b)`, strictly smaller
than the gap `6/(7c)` between two `c` teeth. Hence it touches at most one
third tooth. Each edge contributes at most `1/(7c)`, giving

```text
U_ab <=N_ab/(7c) <=24b/(91c).                          (9)
```

If `c/b>=44/13`, the right side is at most `6/77`. Equality in that ratio
would imply `c=44h,b=13h`; odd `b,c` make this impossible. Thus the bound
is strict throughout this branch.

For `c>=6b`, use a component `I` of length `ell`. A `c` tooth touches it
only if its center belongs to the interval obtained by expanding `I` by
`1/(14c)` on each side. Centers have spacing `1/c`, so the number of touched
teeth is at most `c ell+8/7`. This also bounds the cut-at-zero implementation:
the pair components avoid zero, so splitting a third tooth cannot create an
additional contact with the same pair component. Summing edge lengths yields

```text
U_ab <=T_ab/7+8N_ab/(49c)
     <=12/539+(192/637)(b/c)
     <=12/539+32/637
      =508/7007 <546/7007=6/77.                        (10)
```

Equations (9)--(10) prove the claimed infinite tail closure. The surviving
local obligation now lies entirely in `c/b<44/13`. This is a ratio bound,
not a finite height bound: primitive triples in that region still have
unbounded height.

The concurrent sparse-transport report identifies `U_ab` with
`sum_(C in Omega_w) K_ab(C)`, where `Omega_w` is the complete owner-live raw
carrier support and `K_ab` retains the pair overlap roof. Thus (9)--(10)
also close an infinite region of its new raw-roof inequality with the fixed
global selector `(a,b)`. The full carrier support remains load-bearing;
this statement does not replace it by all carriers with merely a positive
pair roof.

## 5. Exact verification and remaining scope

The standalone verifier is
[lrc14_component_separation_synthesis_sep05.py](../../04-computation/lrc14_component_separation_synthesis_sep05.py).
It imports no repository mathematics. One path uses cleared integer physical
intervals and all six sheet permutations; another uses the reduced pair roof
and Bernoulli formulas. Feasible edge-minimum loads and the identical
individual-edge upper bounds certify every network optimum without trusting
a flow package. Checks remain enabled under `python -O`.

The explicit universe is every primitive distinct positive odd ternary-unit
triple through height `79`: `2,910` rows, of which `233` lie in the new
all-height tail. It rederives all three pair capacities and the physical mass
by all three groupings: `2,909` strict selected certificates and the single
equality `(1,5,11)`. The canonical clipped hostile remains strictly slack.
Additional controls include `(11,13,43)`, `(11,13,47)`, `(11,13,79)`,
`(1,11,331)`, and `(7,11,997)`. The exact source/output are frozen together;
this finite replay is a control on the analytic argument, not a substitute
for its infinite quantifiers.

Reproduction:

```bash
python3 -B 04-computation/lrc14_component_separation_synthesis_sep05.py
python3 -B -O 04-computation/lrc14_component_separation_synthesis_sep05.py
```

Normal and optimized output streams are byte-identical and record `216,039`
explicit checks. Frozen raw-byte SHA256:

```text
source 3164de8aeab6ae2376fd057e86460fd26bc71574100669299d0e2c372ecc72dc
output e52cc962ce15528a276d5bfa84019b10a622ff069fb487c43d54cf86fbb7d8c3
```

The bounds here do not imply that the contact quotient reconstructs physical
mass; do not extend the edge-minimum identity to arbitrary weights; and do
not provide the remaining scale-three bound for `c/b<44/13`. Entry,
synchronization, and LRC(14) remain **OPEN**.
