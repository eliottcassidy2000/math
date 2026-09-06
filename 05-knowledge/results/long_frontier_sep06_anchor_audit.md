# Independent audit: residue envelope and smallest original Laurent phase

**Status: VERIFIED — full analytic and source audit PASS.** Root audited
[the C-only first-phase theorem](long_frontier_sep06_anchor.md), its exact
source and frozen output. No mathematical correction is needed.

The precise scope is five nonnegative roots with e1=13,e2=55 and the
prescribed monic quartic C weakly interlacing B. The prescribed D is used
in the response coefficient law but is not required to interlace. The
theorem therefore applies to the original two-interlacer class and a
larger weak class. The smallest positive zero of the original P(-s)
is unique and simple in (1/110,1/90), and sQ(-s)<-160 there. Other
original phases and actual general Laurent separation remain OPEN.

## Analytic checks

Weak interlacing of two monic real-rooted polynomials of consecutive
degrees gives nonnegative residues in C/B after cancellation. At any
repeated node, the multiplicity difference is at most one; hence no
higher pole remains. The remaining poles are nonnegative, and the
coefficient at infinity makes the residue sum one. The comparison
of rational series then gives exactly the four stated moments.

I independently checked the nonnegative-axis Cauchy inequality and the
three-by-three moment determinant. These imply 75<=e3<=135 and the
stated quadratic bound on e4, including e4<=175/2. The AM-GM e5 bound
holds for zero roots as well. The repeated boundary (0,0,3,5,5) has
the stated positive residue measure; it satisfies C interlacing while
failing the deliberately unnecessary D hypothesis.

The P endpoint and derivative bounds hold on complete intervals using
only these scalar envelopes. They imply no zero before1/110 and exactly
one simple zero before1/90. At an original zero, the e3 elimination is
an identity. The eliminated full response has precisely the six stated
monomials in e4,e5. The discarded square terms and the linear e5 term
are nonpositive throughout the interval. Both positive term bounds are
uniform, and H''<0 plus H'(1/110)<0 prove the full-interval quartic
minimum at1/90. The final rational margin is strictly less than-160.

The alternative e5-elimination Hessian is strictly indefinite for every
positive phase. This refutes a whole-plane quadratic-form certificate,
not the sign on the constrained residue locus. The note keeps that
distinction and does not promote its scalar envelope to an equivalence
with interlacing.

## Independent coefficient construction

The [standalone audit source](../../04-computation/long_frontier_sep06_anchor_audit.py)
imports no producer. It builds the ordinary polynomials
(1+u)^14 s^degree B(u^2/s), and analogously C,D. Coefficient9 yields
the original first equation; coefficient18 of the B square and
coefficient16 of the CD product yield the entire carried response.
This reconstructs the sign and Laurent normalization independently
of the producer's sparse four-variable Hadamard calculation.

Using SymPy exact symbolic algebra, the audit verifies all six
eliminated response coefficients, the original-root identity, the
residue series and Hankel determinant, the whole-interval rational
sign bounds and the alternate indefinite Hessian. It also checks the
weak repeated boundary as a positive control. These28 always-active
gates certify identities and interval bounds, not a sampled theorem.

Normal and optimized audit output are byte-identical. The producer's
106-gate normal replay also matches its frozen output; its independent
ordinary-carrier and two-shape controls were read and checked.

```bash
python3 -B 04-computation/long_frontier_sep06_anchor_audit.py
python3 -B -O 04-computation/long_frontier_sep06_anchor_audit.py
```

Raw LF audit hashes:

```text
source 6934a68c347b21af11a0f448ecb5d91e63786e1c795c300fde129584775309c2
output 2ad110db23aa09a06dfe657138211511391ae748c73e3fbe13b56056358d3b8d
semantic d89bae0937356cae66d8b13b483ae9fe83229b279d20066a8ecfe8b986dc4ddb
```

All scopes, repeated-root cases, inequality directions and original
coefficient consequences passed. No separate-root replacement or
sample-to-universal inference occurs.
