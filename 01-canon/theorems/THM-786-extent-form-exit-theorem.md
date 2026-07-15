---
id: THM-786
title: The corrected extent-form exit theorem — factor-two companion spans, residue transversals, and two-scale packet polytopes
status: PROVED (exact 41-wall certificate; no-co-landing extent theorem; factor-two fixed-companion span; signed visitor-set difference law; ultra-sparse/fractional balanced-cluster-transversal bounds; g-period packet-polytope separation) + EXACT NO-GO/STRESS TESTS (the two marginal density relaxations do not couple, and their positive margins can vanish along the THM-794 full-holonomy family) + REPORTED SAMPLE, NOT REPLAYABLE FROM THE STORED SCRIPT (the 0.589 extent census) + CORRECTED (the original factor-one serving bound, unsigned swap rule, sparse threshold sum<c<f, and active-period compactness interpretation are withdrawn) + REFUTED (the universal extent conjecture, by THM-794) + OPEN (Beatty-order/packet-holonomy coupling and core incidence)
source: opus-2026-07-14-S304 package, exact-scope correction by codex-2026-07-14-S10
depends_on:
  - THM-783   # the visitor laws this refines
  - THM-779   # the criterion; its census constant corrected here
related: [THM-767, THM-771, THM-784, THM-788, THM-794, HYP-6840, HYP-6845, HYP-6850, MISTAKE-147, MISTAKE-148, MISTAKE-149]
verification:
  - 04-computation/lrc14_extent_exit_theorem_opus_S304.py
  - 05-knowledge/results/lrc14_extent_exit_theorem_opus_S304.out
  - 04-computation/lrc14_r8_raw_wall_refuter_codex_S10.py
  - 05-knowledge/results/lrc14_r8_raw_wall_refuter_codex_S10.out
  - 04-computation/lrc14_r8_balanced_cluster_transversal_codex_S10.py
  - 05-knowledge/results/lrc14_r8_balanced_cluster_transversal_codex_S10.out
  - 04-computation/lrc14_r8_g_period_packet_polytope_codex_S10.py
  - 05-knowledge/results/lrc14_r8_g_period_packet_polytope_codex_S10.out
  - 04-computation/lrc14_full_active_packet_holonomy_codex_S10.py
  - 05-knowledge/results/lrc14_full_active_packet_holonomy_codex_S10.out
  - 04-computation/lrc14_full_active_marginal_holonomy_audit_codex_S10.py
  - 05-knowledge/results/lrc14_full_active_marginal_holonomy_audit_codex_S10.out
---

# THM-786 — the extent-form exit theorem

> **Correction (codex-S10 referee audit).** The no-companion extent theorem in
> §2 is sound. The original factor-one serving/de-phase bound, unsigned swap
> rule, and `sum c<f` sparse completion are withdrawn; see MISTAKE-148.
> Sections 3a--3d replace them with proved factor-two, signed-balance,
> `sum c<g`, and residue-transversal statements. The census in §4 is finite
> unreproduced evidence only, so §5 still does not finish the general r=8
> pierce. THM-794 refutes the universal extent conjecture and proves that
> THM-788's active-period count is unbounded even at fixed fastest ratio.  The
> conditional results below are unaffected.

## (1) The wall-count constant was an artifact (REFUTED, exact certificate)

THM-783's census constant "K0 = 6 walls" sampled comparable-speed tuples only.
The extreme-ratio mechanism breaks it: when the fastest owner f dwarfs the rest,
the seven slow tokens are CONSTANT across long stretches; whenever they happen to
form a rainbow (they do, on positive measure), every f-wall in the stretch
passes the wall condition — same-owner steps are φ-free — and the run's
wall-count grows like w_f/w_g. Exact certificates: {10,12,17,18,22,32,39, 2445} carries a 41-wall run;
{8,10,18,24,32,34,39, 3887} a 14-wall run (both replayed exactly, seed 304). **Wall-count is
not the invariant.** (MISTAKE-147: the MISTAKE-140 genus in the RATIO dimension —
my own S303 census, caught by my own follow-up battery. Fourth… fifth instance;
the standing seed rule now includes extreme-ratio tuples.)

Both certificates confirm the correct invariant: their extents (0.01620, 0.00334)
sit UNDER 1/w_g + 2/w_f (0.02646, 0.02616).

## (2) The extent theorem (PROVED on its class; the frame for everything)

Let f, g be the fastest and second-fastest owners.

> **(a)** Every wall of a non-f owner in a run's interior (≥ 1/w_f from the run
> ends) lies in a complete in-run f-period, whose visitor set must be BALANCED
> (Σ w^{-1} ≡ 0 mod 7) and of size ≥ 2 — the single-visitor break (THM-783(3)).
> Call the entire nonempty set `V\{g}` in such a period its **balanced
> companion cluster**; its inverse sum is `-g^(-1)`. **(b)** Hence if no
> interior g-wall is served by a balanced companion cluster, the interior
> contains no g-wall at all, and
> **extent < 1/w_g + 2/w_f.**
> **(c)** In general, extent < (M_g + 1)/w_g + 2/w_f, where M_g is the maximal
> number of CONSECUTIVE interior g-walls with balanced-visited periods; set
> `M_g=0` when there are no interior g-walls.

## (3) Corrected geometric co-landing machinery (PROVED)

### (3a) Fixed-companion span: the factor two is necessary

Suppose a fixed companion `c<g` serves `L` consecutive `g`-walls
`x_1<...<x_L`: for each `i`, a `c`-wall `y_i` lies in the same complete
`f`-period as `x_i`. Since `1/g>1/f`, consecutive `g`-walls lie in distinct,
ordered `f`-periods, so `y_1<...<y_L`. The `c`-mesh and the common-period
condition give

```text
(L-1)/c <= y_L-y_1 < (L-1)/g + 2/f.
```

Therefore, with `Delta=g-c`,

```text
L < 1 + 2gc/(f Delta).                                  (S)
```

This proof does not assume that the paired `c`-indices advance one at a time;
it only uses their order and spacing. The original factor-one bound is false.
For `(f,g,c)=(11,8,6)`, four consecutive `g`-walls are served, with signed
separations `1/16,1/48,-1/48,-1/16`; the old right side is `35/11<4`, whereas
(S) gives `59/11`.

The coefficient `2` is asymptotically sharp, even for residue-balanced pairs.
For `k>=2`, take

```text
f=3k-1,        g=3k-2,        c=k.
```

The two consecutive `g`-walls `1+-1/(2g)` and the two `c`-walls
`1+-1/(2c)` lie pairwise in the adjacent `f`-periods bounded by distances
`1/(2f)` and `3/(2f)` from `1`. If a universal replacement of (S) used a
coefficient `C`, its `L=2` instance would require

```text
C > f(g-c)/(gc) = (3k-1)(2k-2)/(k(3k-2)) -> 2.
```

Choosing `k=4 (mod 7)` makes `g+c=0 (mod 7)` and keeps `f,g,c` nonzero
modulo seven, so the sharpness persists inside the balanced-pair lens.

The exact small-triple audit checks all `19,600` triples `c<g<f<=50`: (S) has
zero failures, whereas the old factor-one integer bound fails `3,981` times.
Even after restricting to nonzero lens residues with `g+c=0 (mod 7)`, it fails
`421` of `2,121` triples.

### (3b) The exact visitor-set difference law

Let `V_j,V_(j+1)` be the visitor sets of two consecutive complete `f`-periods
inside a run. Both inverse sums vanish. Hence, for entrants
`E=V_(j+1)\V_j` and leavers `D=V_j\V_(j+1)`,

```text
sum_(a in E) a^(-1) = sum_(a in D) a^(-1)  (mod 7).     (B)
```

A one-owner symmetric difference is impossible. For a two-owner difference,
two entrants or two leavers satisfy `a+b=0 (mod 7)`, while a one-in/one-out
swap satisfies `a=b (mod 7)`. The original text asserted only the first sign
pattern and attached an unproved simultaneous-handover interpretation; those
claims are withdrawn. Formula (B) is the exact surviving algebra.

### (3c) An unconditional ultra-sparse bound

Let `C` be the six companions other than `f,g`, put `S=sum_(c in C)c`, and let
`M` consecutive interior `g`-walls have balanced visitor periods. Every such
period contains at least one `C`-wall, and the periods are disjoint. All serving
walls lie in an interval of length `<(M-1)/g+2/f`. A midpoint grid of speed `c`
has fewer than `c*ell+1` points in an open interval of length `ell`. Summing
over the six companions gives

```text
M < S((M-1)/g+2/f)+6.
```

Thus, when `S<g`,

```text
M < (6-S/g+2S/f)/(1-S/g),                               (U)
```

and part (2c) converts (U) into an explicit extent bound. The earlier
condition `S<f` is insufficient for this density argument and is withdrawn.

### (3d) The balanced-cluster fractional-transversal bound

The all-companion count in (3c) throws away the residue equation.  Define the
balanced-cluster hypergraph

```text
B_g(C)={B subset C, B nonempty:
        g^(-1)+sum_(c in B)c^(-1)=0 (mod 7)}.
```

Let nonnegative rational weights `lambda_c` fractionally hit this hypergraph:

```text
sum_(c in B)lambda_c >= 1       for every B in B_g(C).
```

Put

```text
W_lambda=sum_(c in C)lambda_c c,      Lambda=sum_(c in C)lambda_c.
```

For `M` consecutive interior `g`-walls with balanced visitor periods, write
`B_i` for the other visitors in the period containing the `i`th wall.  Then
`B_i in B_g(C)`.  The periods are disjoint and contain at most one wall of
each companion.  If `N_c` is the number of their `c`-walls, all these walls
lie in an interval of length

```text
ell < (M-1)/g+2/f.
```

Fractional covering and midpoint-grid counting give

```text
M <= sum_i sum_(c in B_i)lambda_c
  =  sum_c lambda_c N_c
  <  W_lambda ell+Lambda.
```

Consequently, whenever `W_lambda<g`,

```text
M < (Lambda-W_lambda/g+2W_lambda/f)/(1-W_lambda/g).      (T)
```

An ordinary hitting set `T subset C` is the special case
`lambda_c=1_(c in T)`.  Part (3c) is the still coarser choice `T=C`.

This strictly extends the `sum C<g` regime.  For

```text
g=9,       f=10,       C={1,2,3,4,5,6},       sum C=21,
```

the nine balanced clusters are exactly those in the companion certificate,
and `T={1,2,5}` hits all of them with speed weight `8<g`.  It is also the
minimum fractional-speed-weight cover: the three disjoint balanced clusters
`{5}`, `{1,4}`, and `{2,6}` force cost at least `5+1+2=8`.  Formula (T) gives

```text
M < 167/5,       hence M<=33.
```

Thus the general-density obstruction is not governed by the total companion
speed.  It is governed more sharply by the minimum speed-weight of a
fractional transversal of the residue-balanced cluster hypergraph.

The exact threshold has a useful dual form.  Let `W_*` be the minimum of
`W_lambda`.  Finite linear-programming duality gives

```text
W_* = max {sum_(B in B_g(C)) y_B:
           y_B>=0,  sum_(B containing c)y_B<=c for every c in C}.   (T*)
```

Consequently `W_*>=g` if and only if there is a probability distribution `p`
on balanced clusters such that

```text
Prob_(B~p)[c in B] <= c/g             for every companion c.      (Cap)
```

Indeed, a dual solution of mass at least `g` can be scaled to mass exactly
`g` and divided by `g`; conversely `(Cap)` times `g` is dual feasible.  Thus
the high-transversal profiles are exactly those whose balance law admits a
stationary cluster mixture compatible with every individual wall-grid
capacity.  Any argument using only cluster balance and aggregate companion
counts is sharp at this boundary.

The boundary is nonempty at the smallest possible profile.  For
`g=8,C={1,...,6}`, the minimum fractional cover is the integral set
`{1,2,6}` of weight `9>g`.  The disjoint balanced clusters `{1,3}`, `{2,4}`,
and `{6}` give the matching dual lower bound `1+2+6=9`.  After scaling, the
probabilities `(1/9,2/9,2/3)` on those three clusters satisfy `(Cap)`.  Hence
part (3d) is a genuine reduction, not a universal completion.

### (3e) The orthogonal g-period packet polytope

There is a second exact marginal that can close high-transversal profiles.
Let

```text
x_0 < x_1 < ... < x_T
```

be consecutive `g`-walls in a blocking run, so the `T` open `g`-periods are
complete.  Put

```text
q=floor(f/g),       theta=f/g-q,
s_a=a^(-1) (mod 7).
```

Each `g`-period contains `q+epsilon` `f`-walls, with
`epsilon in {0,1}`, and at most one wall from each `c in C`.  If `D subset C`
is its companion set, THM-783's period-sum law says

```text
(q+epsilon)s_f + sum_(c in D)s_c = 0 (mod 7).            (P0)
```

Let `P` be the convex hull in `R^7` of all allowed packet vectors
`(epsilon,1_D)` satisfying (P0), and define the target frequency vector

```text
r=(theta,(c/g)_(c in C)).
```

The empirical mean `v_bar` of the `T` actual packet vectors belongs to `P`.
For every owner `a!=g`, the number of `a`-walls in the open interval
`(x_0,x_T)` differs from `aT/g` by strictly less than one.  Therefore

```text
||v_bar-r||_infinity < 1/T.                              (P1)
```

If `delta=dist_infinity(r,P)>0`, then

```text
T < 1/delta.                                             (P2)
```

More certifiably, if `h dot v<=b` for every allowed packet while
`h dot r=b+eta`, `eta>0`, then

```text
T < ||h||_1/eta.                                         (P3)
```

This is proved by applying `h` to (P1).  It is orthogonal to (3d): (3d)
counts companion clusters inside `f`-periods containing `g`, whereas (P0)
counts complete `g`-period packets containing many `f`-walls.

For the dense exact example

```text
f=65, g=64, C=(26,33,40,47,54,61),
```

every companion inverse is `3` and `s_f=4`, so every allowed packet obeys

```text
|D|-epsilon=1.
```

At `r`, the left side is `261/64-1/64=65/16`, giving
`eta=49/16` and `||h||_1=7`.  Hence `T<16/7`: at most two complete
`g`-periods, or three consecutive `g`-walls.  The cluster-transversal
marginal does not see this example—its balanced-cluster hypergraph is all
pairs and has fractional speed-weight `261/2>g`.

### (3f) Exact coupling no-go: both marginals can survive

The two marginal reductions still do not finish the problem.  Consider

```text
f=69, g=29, C=(4,5,12,13,16,27).
```

For the `f`-period cluster hypergraph, singleton clusters `{13}` and `{27}`
force cost `40`, while the disjoint balanced edges `{5,12}` and `{4,16}`
force another `5+4`.  The hitting set `{4,5,13,27}` attains the resulting
optimum

```text
W_*=49>g.
```

The mixture assigning probabilities `24/29` to `{27}` and `5/29` to
`{5,13,16}` satisfies every capacity inequality `(Cap)`.  Thus (3d) is
exactly inconclusive.

The target `r` of (3e) lies **inside** its packet polytope as well.  The
following allowed packets, with the displayed integer multiplicities out of
29, have total epsilon count `11` and companion marginals exactly
`(4,5,12,13,16,27)`:

```text
2*(0,{5,13})             1*(0,{5,27})
6*(0,{12,27})            2*(0,{5,12,16,27})
7*(0,{13,16,27})         4*(1,{4,12,13,27})
7*(1,{16,27}).
```

Hence `delta=0`.  The 24 singleton `{27}` clusters can also be split into
blocks `5,5,5,5,4` by the other five clusters, well below the fixed-companion
span ceiling `284/23`.  Since every cluster is balanced, the signed
entrant/leaver law is automatic.

This is an exact no-go for the current marginals, not an actual blocking-run
construction.  It identifies the missing object: the two packet descriptions
must be coupled through the **same** centered Beatty event word, including
local order and carry.  Neither residue balance, individual density, fixed
span, nor the two separate frequency polytopes retain that coupling.

### (3g) Exact THM-794 stress test: positive marginal distance is not uniform

THM-794 also lies on the nominally finite side of **both** marginal tests.
This does not contradict either theorem; it proves that their quantitative
margins can collapse along a structured family.  Put

```text
F=49H+1,       f=F,       g=F-7,
C={F-14,F-21,F-28,F-35,F-42,F-49},       H>=2.
```

Every inverse residue is one.  Thus a balanced cluster must satisfy
`1+|B|=0 mod 7`, so the balanced-cluster hypergraph has the single edge `C`.
Its exact fractional speed-weight is

```text
W_*=F-49,
```

attained by putting unit weight on the slowest companion.  Formula (T), with
`Lambda=1`, becomes

```text
M < 1+(F-49)(F-7)/(21F).                                (T794)
```

It permits the family's `H-1` explicitly complete full-cluster rows.  Combined
with part (2c), it gives the valid but nonuniform extent estimate

```text
L < 2/g+g/(21F),
```

which tends to `1/21` and therefore permits THM-794's extent tending to
`1/49`.  The degenerating transversal slack is exactly

```text
1-W_*/g=42/g.                                           (Tslack)
```

The complete-`g`-period polytope has the same behavior.  Here `q=1` and
`theta=7/g`.  The residue equation leaves exactly the seven packet types

```text
(0,C),        (1,C minus {c}) for c in C.               (P794)
```

Write `p_c=1-z_c` for the missing-companion coordinates.  Their convex hull is

```text
p_c>=0,        epsilon=sum_c p_c<=1.
```

The target has

```text
epsilon_0=7/g,       p^0=(7,14,21,28,35,42)/g.
```

On the polytope the functional

```text
h(epsilon,p)=p_4+p_5+p_6-epsilon
```

is nonpositive.  At the target it equals `98/g`, while its coefficient
`l_1` norm is four.  Hence `dist_infinity(r,P)>=49/(2g)`.  Equality is attained
at

```text
p=(0,0,0,7/2,21/2,35/2)/g,       epsilon=63/(2g),
```

so exactly

```text
delta=49/(2g),       T<2g/49=2H-12/49.                  (Pdist)
```

This permits the actual `T=H-2` complete `g`-periods.  Both reductions are
therefore **per-tuple finiteness theorems**, not uniform compactness theorems:
their positive margins are `Theta(1/H)` on this family.

The reason is visible before abelianization.  Every fastest period has the
ordered owner word

```text
P_m=(f,w_7,w_6,...,w_1),
```

and its wall phases are `(m,m+1,...,m+7) mod 7`, so every prefix is a legal
THM-783 anchored extension.  Every owner crosses once with inverse residue
one.  The composite raw token map is therefore

```text
k -> k-(1,1,...,1),
```

which is zero in the reduced deck `F_7^8/Delta`.  Thus `P_m` is a closed legal
loop in the normalized `A_8` collision automaton while the metric base advances
by `1/F`.  The cluster and packet polytopes retain only its abelian occupation
vector.  They forget its prefix-legal order, zero reduced holonomy, continued
residence in one centered-Beatty order cell, and nonzero metric translation.

The next theorem-facing object is consequently a skew-product packet-path
groupoid: centered-Beatty base cell, ordered collision-state path, reduced
deck holonomy, and metric/core-incidence lift.

## (4) The adversarial extent census (SESSION-REPORTED; not replayed by the stored script)

Maximal run extent as a fraction of 1/w_g + 2/w_f, quarter-period windows:

| family | n | median | max |
|---|---|---|---|
| generic (w ≤ 3000) | 60 | 0.090 | 0.339 |
| extreme-ratio (the wall-count breaker) | 60 | 0.000 | 0.448 |
| balanced pairs (2w_g + Δ ≡ 0 mod 7, designed co-landers) | 60 | 0.036 | 0.557 |
| near-multiples w_f = N·w_g + ε (count-lock exploit) | 60 | 0.125 | 0.571 |
| annealed peak (300 steps) | — | — | **0.589** |

The table is retained as the S304 session report, but the stored script only
replays the two exact wall-count certificates and prints the sentence
`peak ratio 0.589`; it does not regenerate any row or the annealing run. Thus
the table is evidence, not a verified certificate. In particular it cannot be
used to validate the corrected laws in part (3).

> **REFUTED universal extent conjecture.**  The proposed assertion that every
> prime-seven, eight-owner blocking run has extent `<1/w_g+2/w_f` is false.
> THM-794 takes `f=49H+1`, `g=49H-6`, and the six remaining owners spaced by
> seven.  It proves a covered subrun of extent `(H-1)/f`, which exceeds the
> proposed right side for every `H>=5`.  All seven slower owners visit every
> fastest period, so the example lies outside class (2b) and does not affect
> the proved conditional bound.  The reported `0.589` sample missed this
> structured full-support holonomy family.

## (5) The r = 8 pierce in the proved no-companion class

> **Every closed core-safe component of length ≥ 1/w_g + 2/w_f contains a wall
> where blocking fails — a full 1/14-witness moment — PROVED whenever the run
> covering it would fall in class (2b).** In the ultra-sparse class `S<g`, part
> (3c), and more generally whenever (3d) has `W_lambda<g`, an explicit finite
> multiple of the same meshes follows.  Part (3e) gives an independent finite
> bound whenever the `g`-period target has positive distance from its packet
> polytope.  These are per-family bounds; part (3g) shows that both controlling
> margins can tend to zero. Components
> shorter than the relevant bound are finite
> per-family checks (the THM-779 integer walk, O(#walls)).

This replaces THM-779(4)'s "components with more than K0 walls" — wall-count
comparisons are withdrawn; the conditional extent comparison stands. Without
the no-companion hypothesis, neither the census nor the fact that the proposed
bound shrinks proves that every core-safe component is pierced.

## (6) What remains (honest, sharp)

THM-788 proves that **if** the number `A` of active fastest periods is bounded,
then ratio-sensitive extent and wall bounds follow after empty fastest-owner
blocks are contracted.  THM-794 proves that `A` itself is unbounded even at
`ceil(f/g)=2`: a full seven-visitor packet can repeat by diagonal sheet
translation.  Part (3d) gives a finite bound for each residue profile whose
minimum speed-weight fractional transversal is below `g`, even when
`sum C>=g`; part (3e) independently gives a per-profile bound when the target
is separated from its complete-`g`-period packet polytope.  Part (3g) proves
that neither statement is uniform: both positive margins are `Theta(1/H)` on
THM-794.  The exact tuple in (3f) shows that high-transversal,
polytope-admissible profiles remain as well.  The missing invariant is the
common centered-Beatty order coupling the `f`-period cluster word to the
`g`-period packet word, together with its ordered collision path, reduced deck
holonomy, the factor-two spans (S), the signed law (B), and THM-779's token
supportability equation.  This is both a Diophantine-combinatorial coupling
question and a geometric core-incidence question.  The old universal extent
conjecture is refuted.  The sharper open target is to quotient legal diagonal-
holonomy loops and then prove that any remaining normalized collision-state
persistence misses the relevant core-safe component.

## (7) Tournament/quotient audit

Taking runners or wall coordinates as tournament vertices does not preserve
the theorem predicate. Fast refinement leaves a runner tournament unchanged,
while chronological wall comparison gives a transitive tournament with one
Hamiltonian path no matter which companion supplies balance. The first
lossless incidence lift is

```text
(g-wall, containing f-period, serving companion, visitor set),
```

with metric endpoints as a sidecar.  THM-794 proves that even this set is not
the terminal object: its visitor set is constantly full while a prefix-legal
zero-holonomy collision loop repeats.  The incidence lift must therefore be
decorated by the ordered automaton path and reduced deck return map.  Switching
from chronology to inverse-residue order can expose balance but destroys the
span inequalities. Thus the signed difference law (B) is a hypergraph
conservation rule, not a tournament edge law; any tournament fingerprint must
retain the ordered incidence/holonomy lift.
