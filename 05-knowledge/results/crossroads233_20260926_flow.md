# Two cutoffs force a gap between lower density and ordinary density

Status: **PROVED + INDEPENDENTLY AUDITED by the root lane; exact computations
independently reproduced. Section 6 was also independently audited by the
geometry lane.** This concerns modified pairing maps
with global two-step descent. Collatz and global longer-horizon peak-price
gluing remain open.

## 1. Inheritance and the decisive target

[THM-4491](../../01-canon/theorems/THM-4491-pairing-two-step-tree-extra-density.md)
proves that the minimum lower natural density of global two-step pairings
is a convergent tree-cost series `alpha`. It constructs a member attaining
that lower density by gluing conditionally optimal prefixes at rapidly
separated cutoffs. Its natural-density attainment question remained open.

The least-used sidecar is the conditional-prefix cost bound
`O(Y(1+log X))`: it makes gluing cheap when the next cutoff `X` greatly
exceeds the already fixed cutoff `Y`. It does not make gluing cheap at all
intermediate scales. The present result tests precisely that loss.

The canonical hostile control is now very small. The clauses

```
epsilon_2+epsilon_3=1, epsilon_3<=epsilon_5
```

give `C(2)=0` and `C(5)=1`, but the minimum of the *sum of both prefix
counts under one assignment* is two. The first optimum wants bit two
zero; the second wants it one. Separately minimized objectives cannot be
added as a jointly attainable objective.

The live concept board is: shared pair bits; costs at multiple cutoffs;
exact residue-dependent fringe costs; lower versus ordinary density; and
the horizon-four boundary of the two-step tree. This is an intrinsic
compatibility relation between prefix optima. It is symmetric, so a
conflict graph and weighted optimization are the appropriate objects;
no tournament orientation is introduced.

## 2. Exact common-assignment objective

For pair `P_i={2i-1,2i}`, the bit `epsilon_i` gives the map

```
0: 2i-1 -> 3i-1,  2i -> i;
1: 2i-1 -> i-1,   2i -> 3i.
```

Let `F={i:epsilon_i=1}` and `A_F(X)=|F intersect [1,X]|`. Global `P_2`
means every source `n>=3` descends below itself within two steps. The
inherited exact clauses, for source-pair indices at least two, are

```
epsilon_(2k)+epsilon_(3k)=1;
epsilon_(3k+1)<=epsilon_(2k+1)<=epsilon_(3k+2).
```

They form a rooted tree on indices at least two, and every feasible prefix
extends to a global member. Put `r=2/3` and define

```
C(X)=min_(global P_2 F) A_F(X),
N(X)=floor(r^2 X),
G(X)=min_(global P_2 F) [A_F(X)+A_F(N(X))].
```

`G` is an exact tree optimization with vertex weight two through `N(X)`
and weight one from there through `X`. Every actual flip is counted once
in each of the two specified prefix counts. There is no private-source
multiplicity or inferred ownership convention.

The first strict incompatibility in this declared cutoff family is
`X=5,N=2`: `G(5)=2>C(5)+C(2)=1`. At the user's corrected cutoff 233,
`N=103`, `G(233)=100`, while the separate optima are 68 and 29. The value
233 selects an experimental universe here; no special arithmetic property
of 233 enters the argument.

## 3. The weighted fringe series

For a full tree of depth `h`, assign weight one to its bottom two levels
and weight two to all earlier levels. Equivalently the root weight is

```
w_h=1 for h=0,1;       w_h=2 for h>=2,
```

and each child carries the same pattern at depth `h-1`. Let `F_h^b(i)` be
the weighted subtree optimum with root bit `b`; write `Delta_h=F_h^1-F_h^0`
and let `tau_h` be the unrestricted root-subtree optimum minus the sum of
the unrestricted child-subtree optima. Initially `Delta_0=1,tau_0=0`.
For `h>=1`, using the child differences `x,y` at depth `h-1`,

```
i even, child 3i/2:
    Delta_h(i)=w_h-x;
    tau_h(i)=min(max(x,0),w_h).

i odd, children (3i-1)/2, (3i+1)/2:
    Delta_h(i)=w_h+min(x,0)+max(y,0);
    tau_h(i)=min(max(-x,0),w_h+max(y,0)).
```

These formulas follow by fixing the root bit in the complement and
monotonicity clauses. In particular,

```
|Delta_h|<=2(h+1),       0<=tau_h<=2h.
```

The residue `i mod 2^h` determines each value. Put

```
B_h=sum_(i mod 2^h) tau_h(i),
gamma=sum_(h>=1) B_h/3^(h+1),
beta=(9/13)gamma.
```

**Proposition 1.** `G(X)/X -> gamma`. Consequently
`G(X)/(X+N(X)) -> beta`.

**Proof.** As in THM-4491, the local nonnegative costs telescope the entire
finite-tree optimum. Every descendant at depth `j` below `i` satisfies

```
(3/2)^j(i-1)+1 <= v <= (3/2)^j(i+1)-1.
```

For fixed `h`, roots in the macroscopic band
`(r^(h+1)X,r^h X]` therefore have their full depth-`h` tree, apart from a
bounded number of endpoints. A descendant at depth `j` in this band lies
below `N(X)=r^2 X+O(1)` exactly when `j<=h-2`, again away from bounded
endpoint errors. Thus its weight pattern is precisely the one above.
The band has asymptotic length `(1/3)r^h X`; its residue average is
`B_h/2^h`. The limiting contribution is `B_h/3^(h+1)`.

For a general tree cut at `X`, weights are at most two. The root-bit
difference and local cost are bounded by twice the depth plus a constant.
Hence roots below `epsilon X` contribute at most
`O(epsilon X(1+log(1/epsilon)))+o(X)`, using the same descendant bound.
This uniformly controls the omitted small-index region. The fixed fringe
bands exhaust the normalized optimum and prove convergence.

The nonnegative series has the explicit tail estimate

```
0<=gamma-sum_(h=1)^H B_h/3^(h+1)
    <=2(H+3)(2/3)^(H+1).
```

The factor `9/13` comes from `(X+N(X))/X -> 1+4/9=13/9`; omitting it
would compare two different density normalizations.

## 4. Natural-density attainment of alpha is impossible

Define

```
ell=liminf A_F(X)/X,       u=limsup A_F(X)/X.
```

**Proposition 2.** Every global `P_2` member satisfies

```
4u+9ell >= 13 beta,
9u+4ell >= 13 beta,
u >= beta.
```

For the first inequality, take a sequence of the *larger* cutoffs `X`
where `A_F(X)/X -> ell`. The exact inequality

```
A_F(X)+A_F(N(X)) >= G(X)
```

then gives `ell+(4/9)u>=gamma`. For the second, choose the smaller cutoff
along a liminf sequence and round the larger one by the factor `9/4`.
For the last inequality, both normalized prefix counts have limsup at most
`u`. These arguments do not assume that the density exists.

In particular, any member with natural density `delta` obeys
`delta>=beta`. Twenty exact coefficients give the rigorous lower bound

```
beta >= b = 4538724347/15109399071
          = 0.30039079156439336... .
```

THM-4491 gives the rigorous upper bound

```
alpha <= U = 3089623223/10460353203
            = 0.29536509552219564... .
```

Since `b>U`, **no global two-step pairing has natural density alpha**.
This refutes the previously open attainment target; it leaves the theorem
that alpha is an attained *lower* natural density completely intact.
Already twelve weighted coefficients suffice for separation:

```
beta >= 684037/2302911 = 0.29703145280039045... > U.
```

Moreover, every member attaining lower density `ell=alpha` must obey

```
u >= (13/4)b-(9/4)U = 13417603/43046721
                     = 0.3116986076593383...,
u-ell >= (13/4)(b-U) = 170854306/10460353203
                       = 0.016333512137142698... .
```

Thus the sparse-scale gluing in THM-4491 necessarily produces substantial
density oscillation. The earlier use of widely separated cutoffs was a
load-bearing feature, not a dispensable convenience.

These are lower bounds on natural density and upper-density tradeoffs.
They do not identify the minimum natural or upper density. No construction
attaining `beta`, or attaining the displayed rational lower bound, is claimed.

## 5. Exact controls and the remaining connection

The weighted coefficient sequence through depth twenty is

```
1,2,6,12,30,60,124,250,498,982,1996,3958,
8008,16062,32318,64712,129682,259672,519980,1041548.
```

Reproduce with
`python3 04-computation/experiments/crossroads233_20260926_flow.py`, or
with `-O`. The [retained output](crossroads233_20260926_flow.out) agrees in
both modes. The program checks 32,766 complete Boolean assignments through
cutoff 15 against the weighted DP, reconstructs valid optimum assignments
at the reported larger cutoffs, and checks all 508 old/new full-tree
residue cases through depth seven by a separate recursive optimization.
It recomputes the inherited unweighted upper bound and compares the
separation and oscillation certificates using exact fractions.

The independent root audit separately reconstructed the first nine
weighted coefficients, the small binary cases, and all final fractions.
The finite ratios are not assumed monotone: for example the joint ratio
at large cutoff 22,500 is 0.3008 and at 225,000 is about 0.300486. The
proof uses nonnegative fringe terms and an exact asymptotic decomposition.

The source of the connection is a family of individually optimal prefix
assignments. Its target is a single shared global assignment. The map is
the joint weighted tree optimization; it preserves every actual pair bit
and both cutoff objectives. Separately minimizing loses their compatibility.
The missing coordinate is therefore not another scalar marginal: it is
the common assignment across scales. The test at cutoffs two and five
isolates that loss, and the fringe series makes it occur at positive scale.

The same issue is relevant to private-to-global pairing certificates in
HYP-9140, but the present theorem does not transport directly to longer
horizons. The inherited global `P_4` counterexample to the two-step clause
shows why the exact tree must be rebuilt with its legal path alternatives.

## 6. Finite positive cutoff kernels strengthen the ordinary-density bound

**PROVED; independently audited by the root and geometry lanes.** The obstruction is not
specific to two cutoffs. Let `a_0,...,a_K` be a finite list of nonnegative
rational coefficients, not all zero, and put

```
J_a(X)=min_(global P_2 F) sum_(k=0)^K a_k A_F(floor(r^k X)),
W=sum_(k=0)^K a_k,
Lambda=sum_(k=0)^K a_k r^k,             r=2/3.
```

Define a full-depth tree with root weight
`w_h=sum_(k<=min(h,K)) a_k`, and each child with the same pattern at
depth `h-1`. Use the difference and toll formulas in Section 3 with these
weights, starting from `Delta_0=a_0,tau_0=0`. Write

```
B_h(a)=sum_(i mod 2^h) tau_h(i),
gamma_a=sum_(h>=1) B_h(a)/3^(h+1),
beta_a=gamma_a/Lambda.
```

**Proposition 3 (finite cutoff kernel).** The exact joint optimum satisfies
`J_a(X)/X -> gamma_a`. Every global `P_2` member has upper natural density
at least `beta_a`; if its natural density exists, that density is at least
`beta_a`. In particular, every finite truncation gives the rigorous bound

```
upper natural density >= (1/Lambda) sum_(h=1)^H B_h(a)/3^(h+1).
```

**Proof.** A vertex receives weight `sum_k a_k 1{i<=floor(r^k X)}`.
For roots in the depth-`h` fringe band, a descendant at depth `j` lies
below the `k`th cutoff precisely when `j+k<=h`, away from the bounded
endpoint exceptions already controlled in Section 3. Thus its weight is
`w_(h-j)`, which proves the stated full-tree pattern. The residue average
again contributes `B_h(a)/3^(h+1)`.

The bound needed for this extension is the *sum* `W`, not the largest
individual coefficient. All vertex weights lie in `[0,W]`. Induction on
the displayed recurrences gives

```
|Delta_h|<=W(h+1),       0<=tau_h<=Wh  (h>=1).
```

The same bound applies to a pruned tree with its actual depth, even though
its vertex weights need not have the full fringe pattern. For `i>=2`,
the descendant inequality bounds this depth by
`log_(3/2)((X-1)/(i-1))`. The roots below `epsilon X` consequently
contribute at most
`O(W epsilon X(1+log(1/epsilon)))+O(W log X)`. Fixed fringe bands and
then `epsilon -> 0` prove the exact limiting series. Its explicit tail is

```
0<=gamma_a-sum_(h=1)^H B_h(a)/3^(h+1)
    <= W(H+3)(2/3)^(H+1).
```

Finally every global assignment has its weighted cutoff sum at least
`J_a(X)`. If its upper density is `u`, that sum divided by `X` has limsup
at most `Lambda u`, so `u>=beta_a`. Existing natural density is the
special case when the lower and upper densities agree. Nonnegative tolls
justify truncating for a lower bound; no empirical convergence assumption
is used.

**Eight adjacent scales.** The bounded comparison tested equal-mass
kernels with adjacent cutoff exponents `0,...,m-1`, using integer
coefficients `a_j=3^j 2^(m-1-j)`. Each cutoff contributes the same mass
`a_j r^j=2^(m-1)` to the density normalization. Among the tested cases
`2<=m<=8` at depth sixteen, eight scales gave the strongest bound. At
depth twenty the retained kernel is

```
a=[128,192,288,432,648,972,1458,2187],
W=6305,                 Lambda=1024.
```

Exact integer recurrence gives

```
sum_(h=1)^20 B_h(a) 3^(20-h) = 3286041555236,
B_20(a) = 3259399538.
```

Therefore every global two-step pairing has

```
upper natural density >= 821510388809/2677850419968
                         = 0.30677978974599224... .
```

The same lower bound applies to its natural density whenever that density
exists. It strengthens the two-cutoff ordinary-density lower bound. The
two-cutoff argument still gives the stronger, conditional upper-density
bound `0.3116986076...` for members attaining lower density `alpha`.

The preceding comparison with six cutoffs at exponents `0,2,...,10`
gave `4915373767757/16067102519808=0.30592782747836345...` at depth twenty.
Adjacent scales capture the complement-edge relation that skipping a
level can miss. These are the best certificates in the declared bounded
comparison, not a claimed optimum over all positive kernels, all numbers
of cutoffs, or all truncation depths. Kernel search stops here.

The retained script prints the seven adjacent-scale depth-sixteen bounds,
the six-scale spacing-two comparison, and all twenty coefficients of the
final certificate. It checks its direct fixed-root subtree optimization
against the difference recurrence for every residue through depth seven
in each kernel run, and checks the final kernel's weighted prefix optimum
against 4,094 exhaustive binary assignments through cutoff twelve.
These add finite controls to the general fringe proof; they are not a
substitute for its scale and quantifier argument. The first nine final
coefficients were separately reproduced by the geometry lane using only
fixed-root costs: `128,640,2416,7792,25440,80650,249609,600912,1320201`.
The root's independent `crossroads233_20260926_scale_audit.py` computes
both fixed-root costs directly through depth twenty, without importing the
difference or toll recurrence, and obtains
`sum_(i mod 2^20) min_b F_20^b(i)=3286041555236`. This is the same aggregate
numerator because each child residue occurs exactly three times in a full
period: if `T_h=sum_i min_b F_h^b(i)`, then
`T_h=B_h(a)+3T_(h-1)` and `T_0=0`.

The mechanism is a finite positive functional separating incompatible
prefix optima. It preserves the actual shared assignment across all
cutoffs, while the density normalizer remembers how much each cutoff
counts. It proves neither the optimal natural density nor a Collatz
descent certificate, and does not transport the two-step tree clauses to
longer horizons.
