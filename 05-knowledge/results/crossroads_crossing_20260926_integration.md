# Crossing integration: corrections and surviving statements

**CORRECTED / PROVED elementary repairs / FINITE-EXACT witnesses, 2026-09-26.**
This audit repairs the parallel [crossing and seeds note](collatz_crossings_20260926_potential_and_seeds.md)
and its routes. It preserves the useful phase clock and source-law question.
It does not infer a flaw in every result of the attached poset manuscripts.

## 1. A balanced comparison is measure-dependent

For the positive-prefix multiplier cells, exact enumeration gives

| length | odd letters | extensions | best balance |
|---|---:|---:|---:|
| 5 | 4 | 3 | 1/3 |
| 7 | 5 | 7 | 3/7 |
| 10 | 7 | 30 | 7/15 |

Thus the displayed finite cells do not all have a perfectly balanced pair.
These witnesses do not settle an eventual asymptotic claim. The cells are
a sufficient multiplier-positive subclass, not every actual no-descent word:
affine carry can rescue a multiplier below one. Total word entropy also
does not equal the entropy of each fixed-weight cell.

A 1/3--2/3 comparison allows a worst-case identification bound of
`ceil(log_(3/2) e(P))` when such a comparison remains available after each
conditioning. It does not establish a `log_2 e(P)` upper bound. A source
residue can select a singleton, making a formerly balanced comparison
deterministic; [THM-4503](../../01-canon/theorems/THM-4503-collatz-poset-bridge-height-selection.md)
and the [cyclic carry-polynomial calculation](crossroads_crossing_20260926_dyadic.md)
retain that lost coordinate.

## 2. The carry clock needs its sign and scope

For successive positive odd states m_i, let
`C_l=product_(i<l)(1+1/(3m_i))`. Natural logarithms satisfy, for l>0,

    3 log C_l < sum_(i<l) 1/m_i <= 4 log C_l.

The lower strict inequality follows from `log(1+x)<x`; the upper follows
from `log(1+x)>=x/(1+x)>=3x/4` for 0<x<=1/3. Already27->41 has
`C_1=82/81`, refuting the former upper bound with coefficient3.
Bounded total carry requires the nonperiodic-orbit summability input from
[THM-4476](../../01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md).
It is not unconditional on the trivial cycle.

The rational address clock starts at Q_0=n. Its 3-adic norm is
`3^(l-v_3(m_l))`, simplifying to `3^l` for l>=1, or at l=0 if3 does not
divide n. Its increments are `Q_l/(3m_l)`. No uniqueness among all possible
clocks or potentials follows from this identity.

Same-band depth includes carry: D is a floor or ceiling of
`M log_2(3)+C_segment`. The exact segment4067->2051 has29 odd steps and
47 halvings, whereas `45<29 log_2(3)<46`. It visits2429 internally in the
same band, so this is not a counterexample specifically about consecutive
band visits. Band-exit counting additionally requires distinct states and
completed stays, with at most one unfinished stay.

## 3. Run resource and geometry

During a valuation-one rise the invariant is
`log_2(n+1)+log_2(3/2)v_2(n+1)`. The shift and positive sign are essential;
the former negative-sign height expression increases. A reset changes
both height and resource. After removing the self-loop at1, the reduced
odd inverse graph has the exceptional root edge5->1. It is not a regular
binary tree at its root. The stronger obstruction now has a separate
proof in [THM-4507](../../01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md).

The simple doubled-angle identity is for x^2-2; it does not identify the
dynamics of all three quadratic maps. A shared golden-ratio growth rate
is an analogy until a predicate-preserving map is supplied. Small measured
mutual information neither proves independence nor excludes every possible
structural relation between phase and valuation.

## 4. Data fidelity and a modularly invisible draft error

The user's35 colours are authoritative. The Wythoff/Fibonacci candidate
matches34 then predicts red where the supplied symbol is blue. Calling
the whole supplied word Fibonacci would discard evidence.

During preparation of this session's carry table, four unwrapped affine
constants were initially transcribed with an added320. Modulo64 this error
vanishes, so the source-residue check alone could not detect it. Direct
word composition gives287,251,227,211. The draft was corrected before its
checkpoint. The survivor is the residue projection; the missing coordinate
is the integer carry, which needs a separate unwrapped check.

## 5. Reproduction

Run `python -B 04-computation/experiments/crossroads_crossing_20260926_integration.py`.
The [script](../../04-computation/experiments/crossroads_crossing_20260926_integration.py)
enumerates precisely the three finite cells and the29-step segment above;
its [retained output](crossroads_crossing_20260926_integration.out) contains
the exact controls. The logarithmic inequalities and graph repairs are
proof arguments, not promoted from a finite numerical test. The session
[audit](crossroads_crossing_20260926_audit.md) replays this independently
from the incoming note's computations.
