# Exact interfaces: recovered principles and three proof extensions

**2026-11-03. PROVED scoped results, independently reviewed; FINITE-EXACT
implementation controls.** This packet reads recent Sessions176/177 and
the incoming gate-flow work against older proved mechanisms, then extends
them. It claims no literature priority or resolution of LRC(14), the
six-colour unit-distance problem, or the unrestricted tournament spectrum.

The final integration also recovered concurrent Session178: its fixed-support
finiteness proof and five-prime test overlap ours and are explicitly credited
as independent derivations. [Exact comparison and scope](INCOMING.md).

## 1. The uniform tournament profile has an exact information dimension

For tournaments with a fixed number n of vertices, let q(k) count ordered
spanning covers by k nonempty directed paths. The rational span of their
profiles has dimension **ceil(n/2)**, and their affine hull has dimension
**ceil(n/2)-1**. These are exact all-order statements, not extrapolated
census ranks. [R1--R4 and proof](profiles/PROOF.md).

The decisive recovery was August's
[THM-3181, half-grid reciprocity](../../../01-canon/theorems/THM-3181-tournament-half-grid-reciprocity-and-repeated-join-recurrence.md):
the path-colour polynomial Q(t)=sum q(k)binom(t,k) has pure parity.
The newest Session176 gave a valid but weaker triangle-based lower bound.
Reversing just the extreme edge of a transitive tournament supplies a
sharper family, with

    Q_(B_m)(t)=t^m+t[(t+1)^(m-2)-(t-1)^(m-2)], m>=3.

Odd-sized members, padded by ordered singleton sums, attain every allowed
least degree. This closes the gap between formal parity and actual
realizable rank. The known leading coefficient is what removes one
coordinate in the fixed-order affine problem.

There is a concrete consequence for the
[recent L/U gates](../circuit_flow_tournaments_20261103/PROOFS.md): readouts
at all prefixes through depth r have fixed-order linear rank
min(r+1,ceil(n/2)), and affine variation rank min(r+1,ceil(n/2)-1).
Successive L gates attain these ranks. For ceil(n/2)>=2, depth
ceil(n/2)-2 suffices and is sharp for affine reconstruction over the hull,
using all prefix readouts and known n. This is not a lower bound on
nonlinear encodings or a claim about one terminal readout.

Independent review: [rank audit](rank_audit/AUDIT.md).

## 2. One-step information can hide a whole module's orientation

If every run-count-specific path-cover start profile of a module equals
its end profile, reversing that module preserves **every one-vertex
insertion count in every ambient tournament quotient**. Vertex-transitive
modules satisfy the condition. Multiple disjoint such modules can be
reversed independently. [C1--C4 and proof](continuation/PROOF.md).

This extends the recovered global converse and ordered-triangle examples
to arbitrary contexts, including strong ambient tournaments. The proof
identifies the lost information: a skew endpoint matrix whose divergence
vanishes, so a single cut cannot read it. Regular degree alone is too weak;
the seven-vertex hostile gives responses403/407.

The quotient still fails under continuation. Two strong five-vertex
tournaments have H15 and identical entire one-step responses, yet the
same two successive attachments give **81 and87** (both intermediate
counts are31). More generally, any one-step response collision changes
arcs in a positive balanced circulation, entirely within the same ordered
SCCs. This is necessary, not sufficient: the capacities are endogenous.

Session27 already supplied full labelled reconstruction using n-1 external
vertices. It is recovered and credited rather than announced again.
Universal two-port reconstruction remains open. Independent review:
[continuation audit](continuation_audit/AUDIT.md).

## 3. Arithmetic feasibility needs valuations, then a complete support

Session177's LCM closure determines possible outgoing **supports** of
sigma transfers. The new exact-row criterion determines their prescribed
**valuations**, including zeros, role restrictions and the order-one
denominator exception. [A1--A4 and proof](arithmetic/PROOF.md).

A small decisive implication is

    11 and23 divide sigma(3^a), a positive even
        => 11^3 divides sigma(3^a).

The order of3 modulo23 is11, forcing an extra11 into a+1 on top of the
11-adic baseline from order5. Thus in a hypothetical odd perfect number
with both actual transfers, **11^4 must divide N**. Merely listing3,11,23
in its support does not imply this clause.

Exact incoming valuation budgets also settle the unbounded-exponent issue
for **each specified finite prime support S**. With A the maximum exponent,
s=|S| and an explicit computable constant C from the order baselines,

    A <= C+(s-1)log_3(A+1).

An integer inequality gives a certified finite exponent box. Complete
enumeration in that box, followed by the full sigma(N)=2N equality, decides
that specified support. The example S={3,5,7,11,13} passes maximal-graph
strong connectivity and the abundancy ceiling, but the bound gives short
valuation contradictions for both Euler roles. It is a mechanism control,
not an improvement to known general odd-perfect-number restrictions.

The same finiteness mechanism and test support were proved concurrently in
[Session178](../session178_prime_flow_honeycomb_20261104/primeflow/PROOF.md).
Our exact valuation-row iff and forced-overflow clause provide complementary
local constraints; its capacity iteration and strict cut intervals provide
further sound pruning. Finiteness is not counted as two separate advances.

The prime support itself still varies without a bound. Exterior factors
remain essential; internal budget equations alone are not asserted
sufficient for perfection. Independent review:
[arithmetic audit](rank_audit/ARITHMETIC_AUDIT.md).

## 4. What unifies the advances, and what does not

The common move is to identify the next permitted operation and the
information it can read, then distinguish formal states from realized
objects. The useful refinements differ by lane:

| Source and map | Preserved predicate | Lost information / required sidecar |
|---|---|---|
| Tournament -> full q profile | Every uniform-context Hamilton count | Named orientation; nonuniform continuation needs endpoint data |
| Module -> one-step response under converse | Every single extra-vertex count under the stated balance | Skew circulation; a second attachment can read it |
| Sigma exponent -> LCM-closed neighbor set | Named target support | Exact valuations, including joint forced overflow |
| Exact sigma rows -> internal budgets on fixed S | Necessary prime-power balance and a finite exponent bound | Exterior prime factors; full sigma equality finishes the fixed-support test |

The rank proof gets its force from **realized basis objects**. The arithmetic
bound gets its force from **the same exponents appearing as global demand
and logarithmically bounded local supply**. The module proof gets its force
from **a symmetry that annihilates exactly the permitted cut observer**.
These are exact mechanisms, not an identification of the three invariants.

This also preserves the earlier Hadwiger--Nelson lesson: a formally
incompatible boundary relation is not a planar unit-distance realization.
No map from the present tournament or valuation invariants to a genuine
six-colour obstruction has been established. The same discipline avoids
mistaking a faithful tournament response for an efficient A000568 census.
No new enumeration claim is made in this packet.

The next proof targets are now narrower: characterize two-port faithfulness
without erasing endpoints; control families of possible prime supports,
rather than only exponents on one support; and realize a forcing colour
relation geometrically before transferring a formal obstruction.

## 5. Evidence, replay and provenance

- Producer rank controls exhaust all1,099 labelled tournaments through
  order5, include2,120 actual uniform probe counts, basis objects through
  order11, and gate-rank saturation through order12/depth6. An independent
  permutation-based rank checker supplies21,624 checks.
- The continuation producer checks5,576 ambient comparisons,656 direct
  permutation counts, and all response collisions through order5. The
  parent independently checks840 context comparisons and all named hostiles
  by literal permutations, without importing producer code.
- Arithmetic controls include17,920 requested valuation vectors and176,400
  periodic lengths. The independent audit checks1,996,500 periodic lengths
  and all9,072 exponent vectors in the two proved support boxes by direct
  geometric sums and final sigma equality.

Run from the repository root, standard-library Python only:

```sh
python3 05-knowledge/results/proof_principles_20261103/verify.py
```

The verifier freezes the expected receipts before replaying all six engines
normally and under `-O`; it never regenerates expected results. It checks
the artifact manifest and exact inherited source pins. New proofs have
independent mathematical review, not just agreeing computations. These
are not Lean proofs or hosted-CI claims. See [board](BOARD.md),
[source pins](SOURCE_PINS.json), and [manifest](MANIFEST.json).
