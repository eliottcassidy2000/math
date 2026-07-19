---
id: THM-1134
title: The r=6 width-16 extrapolation and THM-1102's scalar sequential threshold are false — an exact covering family has T=1043/3>338, and an infinite covering progression has T/k5 tending to 28/27>1 while remaining lonely at t=5/16
status: PROVED as a correction/refutation.  The exact Fraction verifier reproduces THM-1102's window witness, proves the finite covering countercertificate, proves the symbolic unbounded-T lemma, and exactly integrates the rational two-dimensional limit atlas for the infinite progression.  THM-1121's weighted atlas on 92<=k<=332 remains fully valid.  This theorem does not refute LRC: every member of the infinite obstruction progression is explicitly lonely at t=5/16.  What is refuted is completion based only on THM-1102's scalar residual threshold T; stronger sequential geometry and THM-1135's order-free harmonic tail are not refuted
source: codex-2026-07-18-S67 (r6 unbounded-tail referee)
depends_on:
  - THM-1102   # finite width-16 scan whose tail extrapolation is corrected
  - THM-1121   # valid finite atlas; its explicit tail caveat is vindicated
related:
  - HYP-7602   # sequential residual tail is false; use harmonic/arithmetic hybrid dispatch
  - THM-1051   # earlier empirical all-scale R claim
  - THM-1061   # status correctly says the scaling law is measured, not proved
  - THM-1081   # later bottom-window/coarse-tail scan
  - THM-1101   # one-ray decay check, not a uniform tail theorem
  - THM-1135   # different order-free harmonic tail; not contradicted by this refutation
  - THM-1127   # exact fixed-ray torus polygon and eventual quasiperiodic formulas
script: 04-computation/r6_measure_tail_referee_codex_S67.py
output: 05-knowledge/results/r6_measure_tail_referee_codex_S67.out
---

# THM-1134 — the r=6 width-16 tail extrapolation is refuted

## Exact setup

For a seven-speed core `P`, let

\[
 E(P;k_1,\ldots,k_5)
 =\{t\in[0,1]:\|pt\|\ge 1/14\ (p\in P),\
                    \ \|k_i t\|\ge1/14\ (1\le i\le5)\}.
\]

Write `N(E)` for its number of positive-length interval components, `mu(E)` for its
measure, and `L(E)` for its largest component length.  The sufficient threshold used by
THM-1102 is

\[
 T(E)=\min\left\{\frac{N(E)}{6\mu(E)},\frac1{3L(E)}\right\}.
\]

If the sixth ordered killer `k_6` is larger than `T(E)`, the measure horn finds a point of
`E` outside the `k_6`-bad arcs.  This implication is valid.  The problem is the claimed
global bound on `T`, not the measure lemma.

## The finite-window computation remains valid

For the THM-1102 core and displayed window maximizer

\[
 P=\{1,2,4,7,9,11,12\},\qquad
 (k_1,\ldots,k_5)=(158,160,162,164,166),
\]

the exact verifier obtains

\[
 N=206,
 \quad \mu=\frac{191259974959}{1877944733280},
 \quad L=\frac{67}{61992},
 \quad T=\frac{20664}{67}=308.4179\ldots,
 \quad R=\frac{T}{166}=\frac{10332}{5561}.
\]

Thus the floating-point value printed in THM-1102 was accurate.  The scan exhaustively
establishes the maximum over the finite set it actually enumerated: for each core, all
five killers belonged to
`[13 max(P)+1, 13 max(P)+17)`.  Nothing here challenges that finite statement.

## Exact countercertificate outside the window

Keep the same core and take

\[
 (k_1,\ldots,k_5)=(290,292,294,296,298).
\]

Exact interval subtraction gives

\[
 N=370,
 \quad \mu=\frac{214450238441}{1981564300485},
 \quad L=\frac1{1043},
\]

and hence

\[
 \frac{N}{6\mu}
 =\frac{122196465196575}{214450238441},
 \qquad
 \frac1{3L}=\frac{1043}{3},
 \qquad
 \boxed{T=\frac{1043}{3}}.
\]

In particular,

\[
 T-333=\frac{44}{3}>0,
 \qquad
 R=\frac{T}{298}=\frac76>1,
 \qquad
 T-\frac{20664}{67}=\frac{7889}{201}>0.
\]

Now append `k_6=338`.  This speed is just outside THM-1121's certified interval
`[92,332]`, while the THM-1102 sufficient condition also fails because
`338<T=1043/3`.  Moreover, the full family is **Covering** in the exact sense used by the
sieve dispatch: carriers for `q=2,...,14` are respectively

```text
2, 9, 4, 290, 12, 7, 296, 9, 290, 11, 12, 338, 294.
```

Thus this is a literal point in the gap between the two proof pieces, and the non-covering
dispatch does not remove it.

It is **not** an LRC counterexample.  At `t=5/16`, an integer-residue check gives

\[
 \min_{v\in P\cup\{290,292,294,296,298,338\}}\|5v/16\|=\frac18>\frac1{14}.
\]

For comparison, appending `333` instead gives a non-covering family, because no speed is
divisible by 13.  That easier variant is dispatched at `t=1/13`, where

\[
 \min_{v\in P\cup\{290,292,294,296,298,333\}}\|v/13\|=\frac1{13}>\frac1{14}.
\]

The counterexample isolates a failure of the proposed certificate, not of the conjecture.

## Why an interior maximizer did not control the tail

THM-1102 says that its worst point lies at offset 9 inside a width-16 window and concludes
that the window is not truncating the answer.  Interior location only rules out truncation
**within that one finite window**.  It supplies no monotonicity, translation reduction, or
scale normal form relating that window to disjoint windows.  The translated exact witness
above has larger `T`, has `R>1`, and has largest removed killer 298 rather than the claimed
tail cutoff 172.

Consequently, the following phrases in THM-1102 must be read with finite-window scope:

- “max `T=308.4` over the `R>=1` region” means the scanned region only;
- “largest killer among them is 172” means the scanned region only;
- “inside (trustworthy)” verifies the local enumeration, not an unbounded maximum;
- `KB=333` does not follow globally from that scan.

THM-1121 already states this caveat correctly, so its theorem needs no weakening: it proves
the complete finite branch in which all six killers lie in `[92,332]`.

## A theorem: no constant-`T` tail bound can exist

The failure is not repairable by scanning a somewhat wider window and asserting another
constant maximum.  Let `m` run through positive integers with `m == 4 (mod 13)` and put

\[
 K_m=(m,m+2,m+4,m+6,m+8).
\]

For **any** core `P` contained in `{1,...,12}`, the point `t=1/13` is safe for the core and for
`K_m`: the five killer residues are `4,6,8,10,12`, all nonzero.  Hence the residual
`E_m=E(P;K_m)` is nonempty and has a positive-length component.

Because `E_m` is contained in the safe set for speed `m`, and every component of that safe
set has length `6/(7m)`,

\[
 L(E_m)\le\frac6{7m}.
\]

Also `mu(E_m) <= N(E_m)L(E_m)`.  Therefore both entries in the minimum defining `T`
satisfy

\[
 \frac{N}{6\mu}\ge\frac1{6L},
 \qquad
 \frac1{3L}\ge\frac1{6L},
\]

so

\[
 \boxed{T(E_m)\ge\frac1{6L(E_m)}\ge\frac{7m}{36}\longrightarrow\infty.}
\]

This is an infinite, rigorous obstruction to **every fixed numerical bound on `T`**.
This step-2 family alone is not an obstruction to a ratio theorem: the required comparison
is `T(E)<k_6`, and the exact probes at `m=1005` and `m=1733` have `R<1`.  The next normal
form does refute the uniform ratio statement.

## An infinite Covering obstruction to the sequential ratio theorem itself

The preceding step-2 family only kills a constant-`T` repair.  A different translated
normal form kills the proposed sequential max-`T` closure even after imposing Covering.
Keep the same core and put

\[
 A=\{0,10,12,14,18\},\qquad
 K_m=m+A,\qquad k_6=m+19.
\]

Introduce the fast phase `x={mt}`.  In the large-`m` limit, the five killer constraints
become the rational polygonal atlas

\[
 \Omega=\{(t,x):t\in S(P),\ \|x+at\|\ge1/14\text{ for every }a\in A\}.
\]

For a vertical fiber `Omega_t`, write `f(t)` for its measure, `c(t)` for its number of
components, and `g(t)` for its largest component.  The exact arrangement calculation
enumerates every crossing of the ten affine boundaries `x=-at +/- 1/14` and gives

\[
 \int f(t)\,dt=\frac{76004881}{627525360},\qquad
 \int c(t)\,dt=\frac{14813}{13860},\qquad
 \max_t g(t)=\frac9{28}.
\]

The one-dimensional residual is exactly the pullback

\[
 E_m=\{t\in S(P):(t,\{mt\})\in\Omega\}.                 \tag{9}
\]

Separate the finitely many vertical core boundaries, refine the remainder into
positive-width affine cells, and note that for `m>18` the graph `x={mt}` is transverse
to every nonvertical boundary of `Omega`.  Cell-by-cell slicing in (9) gives, as `m`
tends to infinity,

\[
 \mu(E_m)\longrightarrow\int f,
 \qquad
 \frac{N(E_m)}m\longrightarrow\int c,
 \qquad
 mL(E_m)\longrightarrow\max g.
\]

Between consecutive rational boundary crossings the
fiber endpoints are affine; the graph makes `m` turns, so its intersection lengths are
Riemann sums for `f`, its safe entries are Riemann sums for `c`, and its longest crossing
converges to the essential supremum of the horizontal fiber gap over positive-width
cells.  Core-boundary and endpoint-collision terms are finite and therefore contribute
`O(1/m)` after normalization.  The maximum `9/28` is not confined to a vertical boundary:
it occurs on the positive-width cells

```text
[29/126,13/56] and [43/56,97/126].
```

Thus the essential supremum equals the displayed maximum.  Equivalently, the first limit
is elementary equidistribution of `(t,{mt})` on the two-torus, while the other two are the
same finite affine-cell count with endpoints retained.  THM-1127's quasiperiodic machinery
provides the natural route from this asymptotic statement to an exact eventual formula.

It follows exactly that

\[
 \frac1m\frac{N(E_m)}{6\mu(E_m)}
 \longrightarrow
 \frac{111778898}{76004881}=1.47068\ldots,
\]

whereas

\[
 \frac1m\frac1{3L(E_m)}\longrightarrow\frac{28}{27}.
\]

Therefore

\[
 \boxed{\frac{T(E_m)}{m+18}\longrightarrow\frac{28}{27}>1.}
\]

In particular, `T(E_m)>m+19=k_6` for all sufficiently large `m`: the current measure horn
genuinely fails, rather than merely lacking a proved estimate.

This obstruction survives the Covering and clustered gap-family hypotheses and is
explicitly harmless to LRC.  First restrict to the cover-preserving progression

\[
 m=1360+3640n.
\]

The core together with `K_m` is Covering for every `n`: the missing core divisors
`5,8,10,13,14` are supplied periodically by `m,m,m,m+18,m+12`, and `3640` is divisible
by all five moduli.  The ratio limit remains `28/27` on this subsequence, so the measure
horn fails for every sufficiently large member.

Now take the odd-`n` subprogression `m=5000+7280r`.  Since `7280` is also divisible by 16,
a direct residue check at the base point persists along this subprogression:

\[
 \min_{v\in P\cup K_m\cup\{m+19\}}\|5v/16\|=\frac18.
\]

Thus there are infinitely many **Covering, explicitly lonely** families on which the
measure certificate `k_6>T(E_m)` fails.  Exact finite members already show the trend:

| `m` | `k_6` | exact `T(E_m)` | exact `T-k_6` |
|---:|---:|---:|---:|
| 1360 | 1379 | `13157480/9243` | `411383/9243` |
| 5000 | 5019 | `58590280/11271` | `2021131/11271` |
| 8640 | 8659 | `174626200/19461` | `6113401/19461` |

## Audit of the preceding R-ladder

The same logical gap appears earlier in this research arc.

| source | what was actually established | missing universal step |
|---|---|---|
| THM-1051 | exact bounded finite horn; later all-scale addendum | imports the empirical THM-1061 scaling assertion |
| THM-1061 | samples across five decades and three adversarial shapes | its own status correctly says “measured, not proved in closed form” |
| THM-1081 | bottom windows plus coarse tails | explicitly declines a general r=4 claim in its status |
| THM-1101 | a wider bottom window and one decay ray with two killers pinned | no reduction sends arbitrary triples to that ray |
| THM-1102 | every quintuple in one width-16 bottom window | no reduction sends arbitrary quintuples to that window |

Thus no theorem in this chain proves uniform `R` control at arbitrary height.  Finite
experiments are strong evidence and the bounded horns remain exact, but “the worst case
provably sits at the bottom” is not supplied by these scripts.  This audit matters even
when a later independent route may close the same runner family.

## Correct frontier

For ordered killers `k_1<...<k_6`, neither a constant bound on `T` nor a universal
inequality `T<k_6` is true.  The exact remaining statement must instead be a disjunction:

> either the analytic horn has `T(E(P;k_1,...,k_5))<k_6`, or a lift-compatible rational
> chart certifies the same family directly.

The natural variables are the ratio `T/k_6`, scale gaps `k_{i+1}/k_i`, and affine
residue/carry data for a bounded cluster.  Both explicit obstructions are informative:
their residual-component bound fails, yet the finite witness uses `q=16` instantly.  The
missing object is therefore a **hybrid dispatch**: analytic scale separation for low
`T/k_6`, and arithmetic charts for high-ratio clustered resonance.  THM-1102's scalar
max-`T` horn based on five exact removals and one final bound cannot close r=6 by itself;
this does not rule out symmetric discrepancy estimates such as THM-1135.

## Compatibility with THM-1135's harmonic tail

This theorem does **not** say that every measure argument fails.  THM-1135 keeps the
seven-speed core-safe set intact and charges all six danger sets symmetrically through the
oscillation of their centered primitives.  Its sufficient inequality

\[
 \sum_{i=1}^6\frac1{k_i}<\frac{7\mu(S(P))}{6C(S(P))}
\]

is order-free and does not pass through the fragmented five-killer residual `E` studied
here.  In fact, the infinite translated progression above has harmonic sum asymptotic to
`6/m`, so THM-1135 dispatches its sufficiently large members even though the sequential
threshold has `T/k_6 -> 28/27>1`.

The two theorems therefore sharpen the frontier together: THM-1134 rules out the old
max-`T` extrapolation and its uniform sequential-ratio repair; THM-1135 salvages the
unbounded direction by a different harmonic estimate and reduces it to a finite box.  It
does not yet verify every tuple inside that box.

## Assumption challenge and Tournament Analysis

Runner vertices conceal the relevant information.  Better candidate vertices are:

- residual interval components and their boundary events;
- rational-time obligations `(q,a)` as in THM-1121;
- affine translated cluster offsets `k_i-k_1`;
- proof obligations “measure horn succeeds” versus “modulus chart succeeds”.

The verifier nevertheless reports the required tournament diagnostic.  On the five
killer vertices, orient the one-killer-threshold difference and break ties by integer
label.  The result is transitive: score histogram `0,1,2,3,4`, zero directed cycles, five
singleton SCCs, and one Hamiltonian path.  It preserves only a scalar individual-removal
cost and destroys residual-component ownership and every higher bad-set intersection.
The exact interval complex — or its lift-compatible obligation hypergraph — is the
proof-facing object.

## Reproducibility

The dependency-free verifier uses `fractions.Fraction` throughout and freezes the exact
countercertificate, the rational limit-atlas integration, and broader deterministic probes
(translated, scaled, spaced, highly composite, near-consecutive, and geometric).  Live
output is byte-identical to the frozen artifact under both `python3` and `python3 -O`;
proof-bearing checks use explicit exceptions rather than optimization-sensitive `assert`.
SHA-256:

```text
script  e9a1efbfd5ccddf5006b8dd246b5aabbc329f3650b3475a6494b12a424e85599
output  db63edc0639e297670fdb2fe7bdedf2878a0e77407c4e0c6ccecacc6c2334fbe
```

## Named next

1. Prove a scale-separation lemma that dispatches the low-`T/k_6` regime and explicitly
   hands the complementary regime to arithmetic charts.
2. Classify the high-ratio bounded-cluster normal forms by affine offsets and residues;
   the infinite family here shows modulus 16 must be one chart.
3. Search for a finite lift-compatible atlas whose chart is allowed to depend on the
   cluster scale/residue, avoiding THM-1098's fixed-atlas obstruction.
4. Formalize the elementary infinite-family inequality `T(E_m)>=7m/36`; it cleanly
   prevents future finite-window extrapolations from being mistaken for tail theorems.
