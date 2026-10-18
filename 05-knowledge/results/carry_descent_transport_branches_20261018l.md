# Carry companions, transported orientation, and the complete local branch test

**Current scope, 2026-10-18.** The all-prime detector proofs are in
[THM-4711, Legendre/Zagier Wronskians and cube/Apéry carry divisibility](../../01-canon/theorems/THM-4711-legendre-wronskian-and-binomial-cube-lucas-carry-divisibility.md),
with an independent audit. Collatz concentration and six-colour forcing
remain **OPEN**. The computations below strengthen a bounded Collatz test
and prove that no one of the current geometry's 375 local branches can be
eliminated without additional constraints.

## 1. Inheritance and the live concept board

The anchor is the actual Lucas carry, with Collatz transported concentration
and six-colour boundary forcing as the two continuing consumers. The niche
is the older divided-reflection Deuring companion; the wildcard is finite
extraction through ideals, poles, and compactness.

Closest mechanisms: the [Deuring Wronskian](overnight9_20260906_jets_deuring.md),
HYP-9311's factorial-block shift, and the
[source-orientation/375-branch note](source_orientation_boundary_contexts_20261018k.md).
Canonical hostiles: a characteristic-p derivative-zero polynomial x^p;
the source tensor's zero trace despite nonzero direction; and a Petersen
minor that forgets the extendible port colourings. The corrected near miss
is MISTAKE-638's square-factor claim without its carry tangent. The least-used
sidecar is a reflected companion that descends through a quadratic map.

| Live object | Preserved coordinate | Decisive test |
|---|---|---|
| Lucas block (F,V) | Actual integer divided difference | Whole-polynomial Wronskian / common-root ideal |
| Clausen deformation | Parameter tangent, including denominator slopes | Exact coefficient derivative and p-integral tails |
| Apéry pullback | Both sheets, pole orders, analytic constant | Reflection invariance and degree below p |
| ROOT transport | Original source row, terminal and integer clock | Every stopped depth; no trace replacement |
| Six-way contacts | Original factor vertices and eight ports | Full-constraint branch witnesses |
| Compact refinement tree | Compatible prefixes and stopping predicate | Infinite path versus finite acceptance |

The method cards used are **Use redundant paths as detectors**, **Compose
witness fibres before taking their scalar shadow**, and **The same
representation is not the same carrier** in META-PATTERNS. No new generic
method card is needed: the new results supply concrete instances of those
existing rules.

## 2. The completed detector proof and its reusable mechanism

**PROVED, independently audited.** For every prime p>=7:

- HYP-9330(B) holds for central-binomial squares and all six Zagier
  sporadic families, including Franel and Apéry zeta(2).
- For central-binomial cubes, monic gcd(F_p,V_p)=rad(F_p), stronger than (C).
- For Apéry zeta(3), rad(F_p) divides V_p, including singular-fibre roots.

HYP-9330 remains OPEN for the remaining families and its other claims.
The classical Clausen, Euler, and CFVZ identities are credited in THM-4711;
these results do not claim those square factorizations as new.

The finite-degree lemma is a useful small reusable object. If a=-d modulo p
and 2d<p, take the equal-parameter hypergeometric period B and its harmonic
companion C through degree d. Then

```
B^2+z(BC'-B'C)=(1-z)^(2d).
```

The right side is nonzero at B's roots. Thus B is squarefree and C is a
unit at each of those roots. Legendre uses this pair directly; the cube
uses S*B^2 and S*B*C. The different common-root outcomes follow from the
operation on the pair, not merely from the word “weight”.

The cube's deformation uses lower parameter slopes 2 and 1. In dual-number
notation its first jet has the honest product law

```
(H+epsilon*J/2)^2 = H^2+epsilon*H*J,       epsilon^2=0.
```

This is a literal information-preserving lift of the period to its first
response. An arbitrary auxiliary J would not represent the integer carry.
The p-shift calculation establishes that identification before the square
is differentiated. The Euler singular factor is kept in both coordinates.

For Apéry, a direct global parameter deformation is unnecessary. With
Franel H,G and the quadratic deck involution sigma, the combination

```
B = G + eps*K*G(sigma) - eps*G(1/8)*H,
G_A(t(x))=H*B/K,              K=(1+x)^(p-1)
```

retains both sheets. Invariance makes the right side a rational function
of t; its pole divisor makes it a polynomial of degree<p. Only then does
the coefficient recurrence identify it with the actual normalized carry.
This is a concrete route around the truncation obstruction: symmetry alone
would give rational descent, and degree alone would not identify the sheet.
Both are needed.

Independent checks use the unchanged existing detector `mechanism()`:
77 odd primes through 397 for binomial identities, 450 Zagier prime-family
cases, plus the supplemental Apéry checks. The proof is all-prime algebra;
these finite banks are checks, not its justification. The parent rational
coefficient and exact carry checks pass with Python optimization enabled.
No new mathematical MISTAKE was found by the independent audit; the proof
header was synchronized with the expanded family scope before promotion.

## 3. Collatz: every source scale below 2^22, every stopped depth

**FINITE-EXACT.** The unchanged source-preserving C engine now tests every
integer 4096<=N<2^22, a total of 4,190,208 scales. It uses the same fixed
forward-closed 627-seed frontier and the same original odd-source weights,
ports, sibling pairing and integer clocks as the inherited note. Its route
pool contains 8,388,608 odd slots; all eligible routes enter within the
explicit cap. The maximum native depth is 237 and the maximum clock 391.
Overflow and route-cap checks stay enabled.

At every such N and every stopped native depth d:

```
D(N,d)=3E(N,d)-T(N,d) <= D(N,0),
E(N,d) <= P_N[0,0].
```

The prefix trace also never exceeds the original source trace in this
range. The largest prefix-trace ratio remains 0.991828045598 at N=5997, d=66.
The least packet-bound slack is 5115 at N=4096, d=52. The orientation bound
has equality at depth 0 by definition. These are finite statements, including
all d beyond the last arrival because the stopped state is then constant.

The proved fixed-boundary range -72<=D(N,0)<=1428 for all N remains inherited.
It does not make the transported inequality uniform. In particular, positive
increments in D at some later stages, already recorded in the previous note,
prohibit turning this into a stage-by-stage monotonicity claim.

The arithmetic orientation supplied in the question is retained exactly:

```
7=(2-omega)(2-omega^2),     (2-omega)/7=1/(2-omega^2),
omega^2+omega+1=0.
```

The three sibling-row leading cancellations are 10N/7, -8N/7, -2N/7.
Their sum is zero. The Hermitian trace discards the directional quadratic
quantity Re(sum W_u^2) that controls the actual row. Its loss is real even
though the trace is exact. The audit target therefore remains the actual
source row, not a rotated average.

The next sufficient question is unchanged: can the cumulative directional
inequality D(N,d)<=D(N,0), or the packet inequality E(N,d)<=P_N[0,0], be
proved for all N,d with their original source weights? A residue-shell
formula for the initial tensor alone gives no bound on how many original
sources later reach the same port with the same clock parity.

## 4. Six colours: all 375 branches survive the whole selected context

**FINITE-EXACT, with solver-free witness replay.** The universe is precisely
the earlier dense six-coordinate construction with factor sizes
1220,250,246,250,250,250, both mixed quadrilaterals, every axis unit edge,
and all eight named contact families. This is not a census of every ambient
unit-distance edge. The six independent coordinate directions and the
original eight boundary ports are unchanged.

For each of the eleven six-way contact groups, test all 375 choices:
choose two of the six factors, and on each factor require all directed
contact differences to lie in one of the five cosets

```
(p,a) in {(2,0),(2,1),(3,0),(3,1),(3,2)}.
```

Each branch assumption is imposed on the **entire** additive-colouring
constraint system, not merely its own voltage equations. The system has
16,422 original SAT variables and 364,517 clauses before branch selectors.
It retains 2,466 original factor vertices and 6,076 axis edges.

Result: **375 SAT in every group, hence 4,125 of 4,125 branches survive.**
The receipt contains 3,755 distinct full factor-colouring witnesses.
Every named contact relation and selected branch equation is checked on
the colours. An additional replay invokes no SAT solver: it reconstructs
the least support/sumset Horn closure from each colour assignment and
checks every original clause. That is 1,368,761,335 original clause checks,
plus all 4,125 branch assumptions.

This supplies a written obstruction to eliminating the current 375 local
branches: for each branch there is a colouring of the whole present context.
The eight port colours are retained in each witness. This does **not** say
that each branch extends every one of the 570 boundary patterns, or that all
choices of branches in different groups extend jointly. It says that none
can be ruled out individually by the current full set of constraints.

The next useful construction must introduce new, authenticated cross-interior
contacts, or glue contexts on shared original vertices/ports so that their
extension relations conflict. A fresh minor, a voltage quotient, or a local
coset count cannot erase the exhibited witnesses. For a candidate new edge,
the cheapest filter is to inspect which saved witnesses colour its endpoints
the same; the decisive step remains solving the full boundary extension
relation after the exact distance has been authenticated.

## 5. Ideals and compactness: useful certificates, with their missing hypotheses

The ring observation is correct for **commutative rings with identity**:
in a nonzero ring that is not a field, choose a nonzero nonunit a. The
principal ideal (a) is nonzero and proper, since 1 in (a) would make a a unit.
Without identity the statement fails: the additive group F_p with identically
zero multiplication has only the two ideals 0 and R, yet is not a field.

The direct connection to the detector is an ideal, not a shared metaphor.
In F_p[x], (F,V) is the unit ideal exactly when gcd(F,V)=1; it then has no
common zero over the algebraic closure. In the cube case the common ideal
is generated by rad(F). The Wronskian identifies which alternative occurs
and why the multiplicity drops. This turns an infinite family of possible
root evaluations into a single polynomial certificate.

For a countable graph with finite colour set, enumerate its original vertices.
The extendible assignments on the first n vertices form a finitely branching
tree. If every level is nonempty, König's lemma supplies a compatible infinite
colouring. Thus a six-colour obstruction in that countable graph is already
finite. The local branches must still be compatible across levels; a tally
of individually extendible branches is not such a tree path.

Compactness of [0,1] and the infinite-path principle justify this kind of
compatible extraction. They do not imply that a path hits a finite stopping
certificate. In a binary tree, an infinite path is precisely the alternative
to finite exhaustion. In a compactification, a limit can also leave the
intended source class: the positive integers 2^k-1 tend 2-adically to -1,
whose accelerated odd Collatz step (3n+1)/2 fixes it. This limit changes
the source; it cannot certify any fixed positive integer's arrival at 1.

A finite covering argument for Collatz would need independently proved
coverage by finite-entry neighbourhoods, or a uniform tail bound excluding
nonaccepting paths while preserving the marked positive source. Assuming
that coverage in order to invoke compactness would restate the open target.
For Apéry, the analogous finite extraction was available because the pole
orders independently bounded the polynomial degree; no convergence premise
was needed. That exact missing bound is the transferable research question.

## 6. Reproduction and freeze

The manifest records exact sources, receipts and inherited dependencies.
The compressed CNF and all colouring witnesses are retained, so replay does
not require regenerating the algebraic geometry.

```bash
# Frozen geometric certificate; no SAT solver used by this replay.
gzip -dc 04-computation/experiments/charge_contexts_20261018l/six_dense.cnf.gz > /tmp/six_dense_20261018l.cnf
python3 04-computation/experiments/replay_global_branches_20261018l.py --cnf /tmp/six_dense_20261018l.cnf --geometry 04-computation/experiments/charge_contexts_20261018l/six_dense_cnf.json --witnesses 04-computation/experiments/charge_contexts_20261018l/global_branches.json.gz --out /tmp/branch_replay.json

# Regenerate the exact-field geometry and complete branch search.
PYTHONPATH=/tmp/math-numba-20261017i python3 04-computation/experiments/export_six_coordinate_context_20261018l.py --six --dense-six --mixed-quads 2 --patterns 0 --cnf-out /tmp/six_dense_20261018l.cnf --out /tmp/six_dense_cnf.json
python3 04-computation/experiments/global_periodic_branches_20261018l.py --cnf /tmp/six_dense_20261018l.cnf --geometry /tmp/six_dense_cnf.json --out /tmp/global_branches.json.gz

# Same source/clock engine as the earlier independently checked scan.
cc -O3 -Wall -Wextra 04-computation/experiments/source_concentration_20261018k.c -o /tmp/source-concentration-20261018l
/tmp/source-concentration-20261018l 04-computation/experiments/first_arrival_20261018/V_closed_Y1000.txt 22 prefixes > /tmp/concentration_everyN_22.json

# Detector proofs: exact coefficients, actual integer carry, existing mechanism.
python3 04-computation/cm_detector_sporadic/proof_bc_20261018l.py --bound 199 --out /tmp/proof_bc.json
python3 04-computation/cm_detector_sporadic/apery_companion_descent_20261018l.py --bound 397 --out /tmp/apery_descent.json
python3 04-computation/cm_detector_sporadic/apery_pullback_operators_20261018l.py
```

The geometry exporter is a frozen copy of the inherited six-coordinate
producer with an optional CNF output added; the earlier producer and its
manifest remain unchanged. The audit directory contains the independent
binomial, Zagier and Apéry scripts and reports. These are source-preserving
certificates for the stated universes, not universal Collatz or plane-colouring
proofs.
