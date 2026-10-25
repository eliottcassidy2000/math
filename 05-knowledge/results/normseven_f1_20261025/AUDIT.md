# Scope audit and next targets

2026-10-25. Same-session adversarial proof review, plus separately implemented
integer-matrix and literal finite-field-root checks. No independent human or
subagent review is claimed. Source and receipt hashes accompany the packet.

## Proof review

**Root algebras.** Checked direction and transposes: w_j=u^(C_j), so recovery
uses BC=I+mH, not CB with the same indices. The basis proof works over an
arbitrary nonzero ring because all coefficient multipliers are units and
every residue monomial is a free basis vector. It proves necessity as well
as sufficiency. Characteristic two retains the algebra identity but loses
the etale sign-cover interpretation. Distinct branch valuations establish
generic degree64 for every constant quadratic twist in characteristic!=2.

**Fibre products.** Both factors have degree8, their characters span all six
bits, and tensor recovery is over the same K or same affine base ring.
The common point requirement is retained in point counts. Branch-free
algebras are already isomorphic; normalization belongs only to the finite
branched model. Resolutions are not substituted. The geometric free-action
classification is only imported in characteristic zero.

**Monomial lifts.** Direct rational inverse and an independent SymPy inverse
check the displayed determinant-one repair. The old determinant3 describes
an ambient torus isogeny, not a failure of marked quadratic recovery.
The base-dependent coordinate repair is explicitly distinguished from an
automorphism of the old unmarked output torus. The general binary lift is
proved by integer elementary matrices, with parity and determinant checked
by the elimination-word implementation. Negative exponents are units only
on the open torus. Coordinate axes give a sharp affine hostile.

**Integral incidence.** The reduction W4 is bijective over a ring without
dividing by its coefficient characteristic. Forest propagation has only
unit coefficients. Odd free paths eliminate a sum variable with coefficient
one. If none exists, the remaining scalar is annihilated by both N and
beta; Bezout gives gcd(N,beta). The prime hypothesis is used only to reduce
that gcd to one. All prime-power kernel sizes, rather than only field ranks,
determine the Smith exponent. Composite lengths in the theorem are cyclic
partition models; no projective-plane-over-Z/N terminology is used.

**Twist averages and F1.** The twist set is all64 vectors of square classes;
the field has odd order. Each source contributes64 after summing twists.
Triple minors identify all bad characteristics; characteristic3 has its own
formula. The negative torus coefficient excludes a finite split-torus
decomposition over the stated field, using counts over every extension.
It does not exclude all F1 frameworks, and positivity would only have been
a necessary test. Deitmar's unpointed monoid convention and the pointed
zero-adjoined convention are explicitly distinguished.

**Incoming port refinement.** The arithmetic row replacement changes a named
character, so its parity shear is a change of the measurement interface.
Both directions and all proper counting marginals were checked on all eight
states. Restoring x3=y1+y3 restores the original ports exactly; treating y3
as the old third port fails. The fixed-character Laurent lift provides the
alternative that preserves character names. This makes no unit-distance
claim and uses neither the incoming palette count nor a common-origin claim
as an arithmetic existence theorem.

## Candidate ledger

| Candidate | Outcome | First failure / strongest replacement |
|---|---|---|
| Determinant3 is intrinsic to the norm-seven pair | REFUTED | Replacing D mask4 by5 gives determinant1; keep base correction K4 |
| All Boolean bases have the same odd index | REFUTED | Absolute determinants1,3,5,7,9 all occur; 784-case census per pair |
| Every seven-branch complementary pair admits a Boolean unimodular basis | FINITE-EXACT, all120 pairs | An all-rank Boolean version is not claimed; all-rank Laurent quadratic repair is PROVED |
| Prime wrap torsion is exactly p, never p^2 | PROVED | The universal-ring kernel proves the full Smith form |
| The prime bo criterion extends unchanged to composite length | REFUTED | N8 has D2 with o=0; the gcd formula W1 replaces it |
| A single field-rank jump determines the torsion exponent | REFUTED | Full N49 has D49; retain tests over Z/7^e or W2 |
| Formal twist-average value7 supplies seven torus pieces | REFUTED | Coefficient -4 in q-1 prevents a positive split-torus decomposition |
| High clone cancellation supplies equally many legal toggles | Inherited REFUTATION | Incoming genuine tournaments have r1 and unbounded d; keep action and cycle profiles |

These were research guesses, not earlier proved claims. No inherited theorem
is retracted by this packet. The near-miss mechanism is also recorded in the
maintained mistakes ledger to prevent its later promotion.

## Next useful obligations (OPEN)

1. Extend the universal ring calculation to cyclic lengths divisible by3.
   The shift then has several cycles, so a vector of constants and their
   coupled sum constraints must replace the scalar c. First test N9 and N18
   in the *cyclic partition* model, with physical cross-product labels kept
   separate when multiplication by3 is noninvertible.
2. Retain boundary valuations in a compatible monomial atlas for the norm-seven
   covers. A Laurent inverse does not extend over all zero-root strata;
   identify needed charts and normalization maps rather than assert an affine
   inverse. This is distinct from the already closed torification obstruction.
3. Use the complementary arithmetic equations for a joint local-point search
   at a fixed rational base, retaining twist and denominator information.
   No local-global principle is supplied by the fibre-product identity.
4. Classify when higher-rank complementary binary codes admit **Boolean**
   unimodular character bases, and separate that bounded-exponent question
   from the already proved unrestricted Laurent lift. Seven-branch success
   alone does not predict all ranks.

The reusable research move is to specify which ambient changes are allowed
before interpreting a scalar invariant as an obstruction or a decomposition.
The kernel lift, integer Smith module, and actual toggle action are three
different forms of retained data, not interchangeable names for one object.
