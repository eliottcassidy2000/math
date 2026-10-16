# A sharp arithmetic envelope for source-clock refinement energy

2026-10-15. **PROVED:** an exact interval-variance envelope, nested upper and
lower certificates under arbitrary interval refinement, and a square-root
mesh reduction of the dyadic energy problem. **FINITE-EXACT:** all source
labels through 2^29 grounded independently; 29 complete octaves checked.
**OPEN:** a uniform bound on the actual coarse means. Neither universal
Collatz convergence nor unconditional nonvanishing of the amplitude is proved.

## Inheritance and live concepts

The anchor is the integer-clock energy problem. The nearest proved mechanism
is the block Gram/ANOVA compiler in
[the repaired break-clock note](collatz_break_clock_repair_20261013d.md),
especially (5a)--(5c). The canonical hostile is its pair of equal event words
with energies 417 and 81. The corrected near miss is substituting an expected
unweighted break statistic for a source-window bound. The least-used sidecar
is the zero jump at every multiple of three: it improves the elementary
Lipschitz variance coefficient from 1/12 to 1/27.

Concept board: source ports and incoming current; integer lengths; coarse
means; positive refinement innovations; projective scale; prime-power lift
delay. The niche is the all-depth split/inert unit theorem in
[the companion](unit_lift_clock_20261015a.md). The wildcard is the
[clock-preserving duality and block/glue comparison](clock_duality_glue_20261015a.md).
LRC(14), the planar Jacobian problem and plane colouring are not advanced by
these estimates.

## 1. The sharp interval bound

Use the ROOT-supported shortcut parity label t(n) in {-1,0,1}, with zero
for a source that does not reach 1. This definition assumes no universal
convergence. More generally everything here
holds for any real sequence |t(n)|<=1. Define

    W(n)=0,-1,+1 for n=0,1,2 mod 3,
    J(1)=0, J(n+1)-J(n)=W(n)(t(n+1)-t(n))/2.

Thus |Delta J(n)|<=1 and Delta J(3m)=0. For I=[a,a+L), put

    c_a(i)=#{0<=h<i: a+h != 0 mod 3},  0<=i<L,
    V_*(a,L)=sum_i c_a(i)^2 - (sum_i c_a(i))^2/L.

**Theorem M1 (PROVED, all a>=1 and L>=1).**

    D_I(J):=sum_(n in I)(J(n)-mean_I J)^2 <= V_*(a,L)
                                                  <= L^3/27+L/3.       (1)

The first inequality is sharp in the stated class. It does not assume
independence, a gap distribution, Collatz convergence, or a particular
incoming current.

**Proof.** For i<j, telescoping the allowed jumps gives
|J(a+j)-J(a+i)|<=c_a(j)-c_a(i). Apply
`D_I=(1/L) sum_(i<j)(J(a+j)-J(a+i))^2`. The periodic binary labels
t(n)=+1 for n=2 mod3 and -1 otherwise have jumps -1,-1,0 on the
allowed phases and attain equality. They are a hostile/positive control for
this relaxation, not asserted to satisfy the Collatz relation.

For the second inequality, write L=3m+u. The three offsets in
`c_a(3j+s)=2j+offset_r(s)` are

| r=a mod3 | offsets |
|---|---|
| 0 | (0,0,1) |
| 1 | (0,1,2) |
| 2 | (0,1,1) |

Summing squares gives the following exact formulas (omit L=0).

| r | u | V_* |
|---|---|---|
| 0 or 2 | 0 | m(3m^2-1)/3 |
| 1 | 0 | m(m^2+1) |
| 0 or 1 | 1 | m(3m^3+4m^2+3m+1)/(3m+1) |
| 2 | 1 | m^2(m+1) |
| 0 | 2 | m(m+1)^2 |
| 1 or 2 | 2 | (m+1)(3m^3+5m^2+4m+1)/(3m+2) |

Subtracting these from L^3/27+L/3 gives, in table order,
`4m/3`, `0`, `(54m^2+39m+10)/(27(3m+1))`,
`2(18m+5)/27`, `2(18m+13)/27`,
and `(54m^2+69m+25)/(27(3m+2))`, all nonnegative. In particular
the coefficient 1/27 and the linear correction L/3 are jointly attained
when a=1 mod3 and L=0 mod3. ∎

## 2. Every refinement improves both bounds

For a partition P of the same source interval [X,2X), retain its actual
block sums M_I=sum_I J and lengths L_I. Define

    C(P)=sum_I M_I^2/L_I,
    Z(P)=C(P)+sum_I V_*(a_I,L_I).

**Theorem M2 (PROVED).** For every interval refinement P' of P,

    C(P) <= C(P') <= S_X:=sum_[X,2X) J(n)^2 <= Z(P') <= Z(P).          (2)

Both bounds equal S_X at singleton leaves. All formulas use the same
current accumulated from source 1; an octave reset is forbidden.

**Proof.** The lower inequalities are the exact ANOVA identity. For a split
I=A|B with lengths l,r, its innovation is
`lr/(l+r)*(mean_A J-mean_B J)^2`. Averaging the cross-pair increment bound
from M1 gives

    |mean_B J-mean_A J| <= mean_B c_a-mean_A c_a.

The latter is nonnegative. ANOVA applied to the deterministic counter c_a
therefore says that the actual innovation is at most
`V_*(a,l+r)-V_*(a,l)-V_*(a+l,r)`. This is precisely Z(P')<=Z(P).
Induct over splits. ∎

This is a computable upper-certificate mechanism, not merely convergence
of lower projections. Adaptive splits and different terminal controllers
may be used, provided they produce the same actual source sums and retain
their phases and integer lengths. Unauthenticated replacement labels do
not satisfy that premise. No complexity claim for computing M_I is made.

## 3. All fine scales below a square-root mesh are paid

If b divides X and P_b is the uniform mesh, M1 gives

    0 <= S_X-C(P_b) <= X*b^2/27+X/3.                              (3)

For X=2^k and b=2^(floor(k/2)-j), j>=0 with nonnegative exponent,

    (S_X-C(P_b))/X^2 <= 1/(27*4^j)+1/(3X).                     (4)

For odd k the first term can be halved. If b=1 the remainder is exactly
zero. Thus only O(2^j sqrt(X)) authenticated block means are needed to
approximate normalized energy to this deterministic tolerance; computing
those means may still require many source evaluations.

**A sufficient remaining target (OPEN).** For every dyadic X>=2^34, prove

    C(P_(2^floor(k/2))) <= (4/9) X^2.                          (5)

The convenient rational 4/9 is a proposed budget, not a derived asymptotic
constant. Combining (5) and (3) gives
`S_X <= (13/27)X^2+X/3`. Since

    E_X=sum H(n)^2 <= (2 sqrt(S_X)+sqrt(X))^2,

this implies E_X<=2X^2 for X>=2^11. To check the explicit threshold,
divide by X^2: both boundary terms decrease with X, and at X=2^11
the rational square test below succeeds. Consequently (5) on its stated
range would discharge the infinite energy premise of B3 in
[the amplitude note](collatz_amplitude_boundary_20261013c.md).
It would not itself prove universal Collatz entry.

The method also works with any budget below 1/2, choosing a finer j to
leave room for (4). The proof has closed the fine-scale allowance, not
the coarse mean/current obligation. The full root mean is a summand of
C(P_b), so recentering cannot remove it.

## 4. Independent finite evidence

The native checker grounds every source 1<=n<=2^29 by an exact 128-bit
descent to a smaller already grounded source. Overflow and descent caps
fail closed. It then checks summation by parts at every n<2^29 and the
Haar decomposition on all 29 octaves. Five meshes per octave are computed
both by direct block sums and by summed Haar innovations.

The Python audit independently follows literal trajectories to 1 through
2^14; checks all 191,916 splits of all {-1,0,1} label words of lengths
1 through 8 in all three phases; checks the periodic extremizers; and
compares all 27 old break-clock octaves and all 29 overlapping octaves
of the earlier 2^34 discrepancy census. The latter uses (X,2X], so the
audit explicitly adds H(X)^2-H(2X)^2. This endpoint change is essential.

For an exact rational upper bound Z, the energy certificate is checked
without floating point: with A=2X^2-4Z-X require A>=0 and A^2>=16XZ.
All square-root-mesh certificates pass for the computed k=14 through28.
Selected values (decimals display exact-rational computations):

| k | b | actual S_X/X^2 | upper Z/X^2 |
|---|---:|---:|---:|
| 24 | 4096 | 0.442190256 | 0.479222033 |
| 24 | 2048 | 0.442190256 | 0.451446864 |
| 26 | 8192 | 0.270996435 | 0.308031002 |
| 27 | 8192 | 0.208782824 | 0.227300096 |
| 28 | 16384 | 0.043142579 | 0.080178391 |

Reproduce from the repository root:

```sh
c++ -O3 -std=c++17 04-computation/experiments/collatz_clock_mesh_20261015a.cpp -o /tmp/clock-mesh-20261015
/tmp/clock-mesh-20261015 29 05-knowledge/results/collatz_clock_mesh_20261015a.json
python3 04-computation/experiments/collatz_clock_mesh_audit_20261015a.py 05-knowledge/results/collatz_clock_mesh_20261015a.json 05-knowledge/results/collatz_clock_mesh_audit_20261015a.json
```

Data: [native census](collatz_clock_mesh_20261015a.json),
[independent audit](collatz_clock_mesh_audit_20261015a.json).

## 5. Connection contract and the next decisive test

Source: an actual bounded label sequence and its integer-clock current.
Target: a hierarchy of interval means with a deterministic variance budget.
Map: conditional averaging on the actual interval partition. Preserved:
the source address, length, current offset and exact orthogonal energy.
Lost by averaging: within-block fluctuations. Sidecar restoring a certified
upper bound: V_*(a,L), including the phase a mod3. Cheapest hostile test:
the monotone periodic labels saturate the envelope and have cubic energy.

The [Q--Pell model](golden_nested_trees_20261013c.md#6-incoming-q-covariance-work-a-genuine-two-step-pell-mechanism)
and [golden-tree controls](gold_heavy_child_20261014j.md) guide the search
for a state closed under refinement. They do not identify the actual
Collatz coarse means with that state. Recent golden comparison work shows
why one-step marginal quality can fail after conditioning; recent
[THM-4655](../../01-canon/theorems/THM-4655-strong-tournament-balance-lemmas.md)
shows a positive model where marked ports and a complete transfer matrix
do propagate through arbitrary context. The next test is a deterministic
transfer inequality for the vector of source block sums, with its full
boundary forcing retained. A projective fixed point alone cannot supply it.
