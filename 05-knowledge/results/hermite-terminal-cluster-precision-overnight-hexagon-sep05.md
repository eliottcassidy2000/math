# The full Smith list needs unit data; the sharp precision does not

**PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.** The
terminal-cluster argument gives an exact metric-only formula for the largest
two-jet Smith exponent at every finite integer node set. Independent full
proof and literal-Smith reviews passed. The canon route is
[THM-4439](../../01-canon/theorems/THM-4439-all-node-twojet-metric-precision-by-terminal-clusters.md).

## 1. Inheritance and the question recovered from a failed extrapolation

[THM-4435, metric blindness and universal Hermite precision](../../01-canon/theorems/THM-4435-four-node-metric-blindness-and-universal-hermite-precision.md)
refutes distance-only full four-node partitions. Its positive replacement
gives the largest exponent from F'(x_i),F''(x_i), retaining a unit-sensitive
coordinate rather than claiming that coordinate is always necessary.
[THM-4429, arbitrary three-node Smith forms](../../01-canon/theorems/THM-4429-arbitrary-three-node-two-jet-smith-form-and-metric-precision.md)
already has a metric-only largest factor, with a ternary full-residue wall.

The canonical hostile is the all-scale isometric pair
2^e(0,1,2,5),2^e(0,1,3,4): full partitions differ for e>=3, but the
largest exponent remains7e+5. The corrected near miss is to transfer the
intermediate-minor failure to the largest inverse denominator. The least-used
operation is now **maximize over a complete terminal residue cluster before
discarding unit data**. Pointwise inverse columns can cancel; the complete
cluster has a low-degree polynomial obstruction to simultaneous cancellation.

The live board is: full Smith partition; exact inverse columns; p-adic
distance tree; terminal residue clusters; local reciprocal cancellation;
maximization before quotient. The source-to-target map is the node set to
its unlabelled distance tree. It preserves the largest precision exponent
by the proof below, but demonstrably destroys intermediate
Smith factors. Keeping the target invariant distinct is load-bearing.

## 2. Precise terminal-cluster statement

Let p be prime and let X be a finite set of n>=2 distinct integers. Write
v=v_p and put S_x=sum_(y in X,y!=x)v(x-y). A **terminal cluster** is a
set C with a depth f>=0 such that, for any x in C,

```text
C={y in X : v(y-x)>=f},
|C|=m>=2,
v(y-z)=f for every distinct y,z in C.
```

These are exactly the nonsingleton vertices of the distance tree whose
children are all leaves. Their normalized residues are distinct, so m<=p.
Every finite tree with at least two leaves has such a vertex. For x in C,
the quantity S_x is constant; denote it S_C. Distances to any point outside
C are equal for all x in C, by the ultrametric inequality.

The exact sharp precision loss is

```text
L_p = max_(terminal C) [2S_C+max(0,f_C-[|C|=p])].       (T1)
```

The indicator is one precisely when the terminal cluster fills every
residue at its branching scale. For n=1 the loss is zero. As usual, for
every N>=1, observations modulo p^(N+L_p) recover all coefficients of a
degree<2n polynomial modulo p^N, and one less digit fails uniformly.
The data are values and first Hasse derivatives, with every node included.

In particular the sharp loss, and hence the largest integer Smith factor
prime by prime, depends only on the pairwise valuation tree at **all node
counts**. This does not restore a metric-only full partition.

## 3. Nonterminal leaves cannot control the worst denominator

For x in X define its nearest depth f_x=max_(y!=x)v(x-y) and
q_x=sum_(y!=x)1/(x-y). The exact inverse-column formula from THM-4435 is

```text
L_p=max_x max(2S_x, B_x),
B_x=2S_x-v(2q_x), with v(0)=infinity.                  (T2)
```

The valuation inequality gives v(q_x)>=-f_x, and therefore
B_x<=2S_x+f_x-v(2). Put T_x=2S_x+f_x.

Suppose the ball C_x={y:v(y-x)>=f_x} is not terminal. The child containing
x is a singleton, since no point is closer to x. Some other child contains
two or more points. Choose y in that child. Distances to points outside
C_x are unchanged; every inside distance at y is at least f_x and at
least one is greater. Consequently

```text
S_y>=S_x+1,   f_y>=f_x+1,   T_y>=T_x+3.               (T3)
```

Repeat along a deeper nonsingleton child until a terminal cluster is
reached. Thus every nonterminal leaf has both smaller S and a T value at
least three below some terminal cluster's T. In particular max_x S_x
is attained in a terminal cluster. Once Section4 proves that a terminal
cluster attains B at least T-1, no nonterminal leaf can maximize B.
This step concerns the worst inverse denominator, not every local column.

## 4. Exact cancellation budget inside a terminal cluster

Fix C of size m and depth f. Choose an integer center a in C, and write
x=a+p^f*u_x for x in C, with u_x integral and pairwise distinct modulo p.
Let

```text
G(U)=product_(x in C)(U-u_x),
Q_x=p^f*q_x=sum_(y in C,y!=x)1/(u_x-u_y)+E(u_x),
E(u)=sum_(z outside C) p^f/(a+p^f*u-z).               (T4)
```

Every outside depth is at most f-1. Hence E(u) belongs to p Z_p and
E(u)/p modulo p is constant as u varies: a term with outside depth f-1
has constant reduction after division by p, while all shallower terms
vanish. More explicitly, the difference between its values at u and v
is divisible by p^2, since it has numerator p^(2f)(v-u) and denominators
whose product has valuation at most2f-2. If f=0 there are no outside
points, so the same statement holds with E=0.

The internal reciprocal sum is G''(u_x)/(2G'(u_x)), where G'' is the
ordinary second derivative. All G'(u_x) are p-adic units.

### Odd prime, fewer than p children

If p is odd and2<=m<p, the reduction of G'' is nonzero of degree m-2,
with leading coefficient m(m-1). It cannot vanish at all m distinct
residues. At some x the internal reciprocal sum is a unit, and E(u_x)
cannot cancel it. At every x it is integral. Therefore

```text
min_(x in C) v(Q_x)=0,
max_(x in C) B_x=2S_C+f.                            (T5)
```

### Odd prime, all p children

If p is odd and m=p, reduction gives G(U)=U^p-U modulo p. Thus every
coefficient of G'' is divisible by p, and G''/p has degree p-2 with
unit leading coefficient p-1. Also G'(u_x)=-1 modulo p for every residue.
The reduction of Q_x/p, as a function of the residue u_x, is consequently
a polynomial of degree p-2 with nonzero leading coefficient, plus the
constant contributed by E/p. Since p>=3 that degree is positive, so the
constant cannot cancel the polynomial identically. A degree p-2 polynomial
cannot vanish at all p residues. It follows that

```text
min_(x in C) v(Q_x)=1,
max_(x in C) B_x=2S_C+f-1.                          (T6)
```

This proves an exact one-digit simultaneous cancellation budget. It does
not assert that every local column has exactly that loss.

### The dyadic boundary

At p=2 every terminal cluster has exactly two children. The internal sum
at each point is a single unit reciprocal, while E is even, so v(Q_x)=0.
The factor two in(T2) supplies precisely the one-digit subtraction:
max_(x in C) B_x=2S_C+f-1. This agrees with(T6)'s indicator, although the
mechanism is the ordinary derivative's factor two, not full-residue
reciprocal cancellation. Confusing the mechanisms would break this case.

Combining(T5),(T6) and the dyadic calculation gives
max_(x in C)B_x=2S_C+f-[m=p]. Section3 excludes all other leaves from the
global maxima of S and B. Substituting in(T2) proves(T1).

## 5. Controls, limits and decisive tests

The companion probe independently compares complete F'/F'' denominator
constants with an unlabelled distance-tree compiler. It compares the two
affine dilation costs separately, so a shallow maximum cannot hide a
different eventual precision. Its complete head uses all translated-to-zero
sets of4..6 nodes and diameter<=24 at p=2,3,5,7. A seeded lift mode uses
4..30 nodes and60,000 high-unit isometries at p=2,3,5,7,11,13. These
finite tests do not prove(T1); they challenge its sharp cancellation boundary.

```bash
python3 04-computation/hermite_metric_loss_probe_overnight_hexagon_sep05.py --height 24 --max-nodes 6
python3 04-computation/hermite_metric_loss_probe_overnight_hexagon_sep05.py --random 10000 --max-nodes 30 --primes 2 3 5 7 11 13
```

Normal and optimized versions of both primary commands have byte-identical
outputs. An independent referee separately reconstructed the reciprocal
constants and terminal tree on7,944 rows and compiled176 literal integer
Smith comparisons. Its complete small head uses n=1..6 and diameter<=12;
1,600 additional seeded irregular/deep rows and named n=7/9 controls retain
singleton/two-node cases, full terminal p-clusters with nearby outsiders,
both one-digit mechanisms, unequal branch sizes, signed nodes and metric
twins at depths0..4. The full proof independently passed, including the
degree-versus-full-residue cancellation obstruction and strict descent.

```bash
python3 04-computation/hermite_terminal_precision_referee_overnight_hexagon_sep05.py
python3 -O 04-computation/hermite_terminal_precision_referee_overnight_hexagon_sep05.py
```

The independent normal/optimized outputs match. The root also replayed
the referee successfully. Frozen raw LF-byte manifests:

```text
primary source 741a2cf27783d2320404b18e9e14e988ee61bb0929c86320137a5b06149457c4
head output d5c18a1849ed2aa525e42d00460d4f9e08c8d8fa9d7c79c738000c9530d85554
lifts output 8d58361d558f787f4368116611e894954c39bd71fb41b454acacd734d8491005
referee source 2a74b09d1a6f57f9f8e4bb42e2ac09138d8cd794ce7299a4303442861ffeb110
referee output 1ce7d0017bca294abe1af61703517542ee3f3609d85c19a00b4a2f3883d6eb85
```

Higher Hasse jets have different local coefficient sums and are not
asserted to obey this terminal formula. The first future test must preserve
jet multiplicity as well as the complete terminal residue cluster.
