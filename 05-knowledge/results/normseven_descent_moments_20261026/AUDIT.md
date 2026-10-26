# Scope, independent computational paths, and incoming-work audit

2026-10-26. Internal proof review plus independent implementations; no
independent human/agent proof audit and no formal proof-assistant build.
Every load-bearing gate raises an exception, so Python optimization cannot
remove it. Finite checks corroborate the proofs, not their infinite quantifiers.

## Inheritance and overlap

Startup was clean at `983768efde84`; the primary checkout's extensive unrelated
changes were preserved. The requested norm-seven, six/seven, twist/conductor
packets, navigation, guardrails, research protocol, core-paper index and bounded
messages were read. Searches covered determinant/congruence lifts, second
moments, normalization, inertia and relevant recent corrections. Session65/66
and the corrected MISTAKE-673 cycle claim were inspected for transfer limits.

The first local checkpoint `9e9b0c5e01` contained the completed moment and lift
proofs/checks, with the packet README explicitly still under construction.
Its push lost a race to `2daf5e00af`. The clean checkpoint was rebased onto that
incoming arithmetic/F1 packet. Its results fully overlap the main M1-M3 and
L1-L3 mechanisms, including the explicit integer matrix and both energy
polynomials. This is recorded as shared confirmation, not novelty.

The incoming correction to THM-4773 outranks the old open-central dichotomy:
eight rational nodes can exist without open rational points. Incoming boundary
C2-C3 then distinguishes rational nodes, divisor normalization and surface
resolution through the tangent norm conic. RESOLVED_FIBRES extends this to
nonsplit node fields and the two-versus-four minimum-degree comparison.

## Mathematical review of the new result

1. **Objects and maps:** the source is one rational point of the finite coarse
   quotient product, not an arbitrary point of either surface. Its preimage
   consists of eight nodes, not eight components of the whole surface.
2. **Arithmetic descent:** the eight nodes form four copies of Spec(Q(sqrt(D)))
   when D is nonsquare. Their exceptional conics are K-schemes regarded over Q;
   they are not Weil restrictions and have relative dimension one.
3. **Local equation:** in s3=1 the first three root equations give exactly
   s1^2+s2^2=delta*s0^2. The other roots have nonzero derivatives at each node.
   Blowing up inserts its projective tangent cone; no point on the conic is
   produced merely by normalizing the divisor.
4. **Necessity and sufficiency:** an F-point must first embed K in F, then solve
   the conic there. Conversely these data give a point over the original y.
   K=Q(sqrt(13)) embeds in Q3 by a simple-root lift, so the mod3 infinite-descent
   obstruction to A^2+B^2=3C^2 applies. No rational Gaussian-norm criterion is
   extrapolated to K and no Hasse principle is imported.
5. **Minimal degrees:** delta2 has a rational conic direction and still needs
   K, so degree2 is exact. Delta3 fails over K; every containing extension has
   even degree, and K(sqrt2) of degree4 supplies a direction. This rules out
   degrees1,2,3 and constructs degree4, not just an upper bound.
6. **Good-prime counts:** exclude 2,3,13 explicitly. The residue-node gate is
   unchanged and every nonsingular conic has q+1 finite-field points. The
   complete fibre zeta functions agree at every good prime, including every
   extension; the global surface zeta functions are not being identified.
7. **F1 scope:** torification is relative to K. Q2 splits over K, Q3 has no
   K-point and cannot be a union of nonempty split K-tori. Both complete
   fibres have no Q-point; neither has such a decomposition over Q.

## Exact finite universes and independent paths

| Checker | Declared universe / independent comparison |
|---|---|
| check_moments.py | Every inertia subspace and Frobenius translation in dimensions1..4, coordinate complements in every dimension split: 5,940 cases. Literal cosets/products compare with the stabilizer formulas. |
| check_moments.py | Both norm-seven pairs, all64 twists, every P2 point in F3,F5,F7,F11,F13,F9,F25. Direct root counts compare with full surviving-character algebras and independent coset membership. Fixed sqrt5 hostile at all primes7..300. |
| check_lifts.py | Every invertible 2x2 residue matrix at m=2..12, with constructive same-rank lifts or compensating-coordinate lifts. Also the displayed six-to-seven example and 903 marked root fibres in cyclic groups of orders1..42. |
| audit_independent.py | Imports no producer. BFS of reductions of elementary integer matrices plus a swap recovers the entire predicted determinant-sign subgroup for m=2..12. Independent Leibniz determinant at ranks2/7; 541 unit-root fibres in Z/N for N=2..42, including noncyclic unit groups. Independent 35 integer line determinants. |
| check_resolved.py | Direct projective conics for every nonzero D,delta over F5,F7,F11,F17,F19,F25,F49: 3,612 pairs. Separate tuple arithmetic; original seven equations at32 nodes; exact quadratic-field witnesses; all729 triples mod9 for the primitive obstruction; two roots of13 modulo3^e for e=1..9. |

No input filters beyond the stated fields, ranks, nonzero parameters and good
prime exclusions are hidden. The determinant checker's cases include composite
moduli and pivots that are nonunits until row addition; the second path does
not reuse its elimination. The conic checker includes both split and nonsplit
D and all delta classes; no search failure is used as a number-field proof.

## Rejected guesses and implementation repairs

| Guess / issue | Verdict and first failed step | Strongest survivor / missing datum |
|---|---|---|
| Every triple loses three bits | REFUTED by characteristic-three136: one kernel meets inertia | Use e=r-h-k; incoming reduction R1-R3 independently closes it |
| Equal twist means identify normalization fibres | REFUTED: (0,1) and (8,1) at an original triple | Retain nested stabilizers and labelled Frobenius, not just mass |
| Splitting the residue torsor completes a resolution lift | REFUTED by D13,delta3 | Retain the tangent conic over K; adjoining sqrt2 finishes this example |
| All good local zeta factors determine this arithmetic lift | REFUTED by delta2/3 over K=Q(sqrt13) | Bad-place solubility remains indispensable |
| Same-rank marked root automorphisms always lift | REFUTED by exponent2 mod7 (already exponent2 mod5) | Determinant +/-1; one nonidentity compensating character |
| Modular elimination finished after clearing lower columns | Implementation failed its own identity gate before any push | Added the final upper-column elimination; full independent BFS and determinant replay pass |

The temporary +349 characteristic-three coarse second-moment guess was
rejected before publication: two new triples each have d=2, not1, giving +351.
No inherited characteristic-zero freeness claim was false. No new correction
to the shared MISTAKES ledger is needed; the actual incoming THM-4773 correction
is retained in the inspected mainline history.

## Reproduction and integrity

Run `python3 -B verify_packet.py --replay` from this directory. The wrapper
checks the exact closed file list, SHA-256 source pins, four normal/optimized
receipt pairs and all eight replays. Saved text uses LF; subprocess universal
newline decoding makes comparison portable across Windows and Unix. It changes
no mathematical output. Run `python3 agents/check_docs.py` from repository root
for the maintained navigation surface. `--freeze` is an explicit author action
after reading incoming changes; ordinary verification never rewrites pins.

The manifest pins all proof files, programs, receipts, audit and source-pin
file. No new sources are assigned authority from packet titles alone; the
underlying proofs and their limits were read. External checks used Stacks
Tag04UV for the quotient-torsor convention and Deitmar sections1-3 for monoid
base extension. All new degree and count arguments are included in full.
