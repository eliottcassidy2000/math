# Independent audit: conductor-ten synthesis and publication copies

20 October 2026. **PASS in the stated primitive monic source scope.**

## Verdict

The theorem in `CONDUCTOR_TEN_SYNTHESIS.md` follows from its exact audited
dependencies: if a complex planar Keller pair has one positive-degree
coordinate monic in y, satisfies `C(x)(P,Q)=C(x)(y)`, and is not a
polynomial automorphism, its generic source-algebra codimension is at
least five and its generic normalization-conductor degree is at least ten.

The same statement with constant nonzero leading coefficient follows by
constant rescaling. This is not a bound on the field degree of the planar
map, not a uniform bound on every special fiber, and not a resolution of
the planar Jacobian conjecture.

Reviewed synthesis candidate SHA-256:

```
5ef901e577c30ee4c4baab7b7e0cc83c2fc557e27de761b8380ad0129eaaddc2
```

The publication copy and the reviewed source are byte-identical.
The final status-only promotion was independently reconciled: restoring
the candidate's first status line recovers exactly the hash above.
The promoted source and publication copy both have SHA-256:

```
4366ee86a86d4d8d435b18b4aa6ac886d89576b757fb76104fd07edea8270ee4
```

`SYNTHESIS_PROMOTION_RECONCILIATION.json` records that exact comparison.

## 1. Algebraic setup and the zero-defect argument

Let K=C(x), A=K[P,Q], B=K[y]. The positive y-degree of P makes B finite
over K[P]. Since `K(P,Q)=K(y)`, the two K[P]-modules A and B have equal
rank. Their quotient is a finite torsion K[P]-module and hence has finite
K-dimension. This justifies finite delta and the existence of a nonzero
conductor. The trace-dual proof in the low-delta dependency then supplies
`deg_y c=2 delta` and the symmetric numerical degree semigroup.

If delta=0, A=B. Membership of y in A means one finite polynomial
expression in P,Q with coefficients in C(x). A single nonzero common
denominator D clears every coefficient, giving a polynomial identity

```
D(x)y = Σ r_ij(x)P^iQ^j.
```

Choose a complex number a outside the finite zero set of D. Specializing
this identity at x=a shows
`C[P(a,y),Q(a,y)]=C[y]`. Thus the restriction to the genuine affine source
line x=a is injective. There is no illicit specialization of a denominator
at a pole and no assumption of saturation or constant fiber delta.

The original primary PDF of Gwoździewicz, *Injectivity on one line*,
Theorem 1.1, was freshly reopened for this audit:
https://arxiv.org/pdf/alg-geom/9305008v1 .
Its hypotheses are an algebraically closed characteristic-zero field,
a polynomial plane map with nonzero constant Jacobian, and injectivity
on an affine line. All match this application. Therefore delta=0 implies
an automorphism.

## 2. Complete coverage of positive delta at most four

The complete symmetric-semigroup list in the synthesis has seven cases.
The dependency applications are exact:

- `<2,3>`, `<2,5>`, `<2,7>`, and `<2,9>` are excluded by the
  independently audited all-quadratic-semigroup theorem. Its statement
  permits every odd m≥3, including composite m, and has exactly the
  required constant-leading actual-Keller hypotheses. It does not require
  an actual quadratic coordinate or a unique plane embedding.
- `<3,4>` is excluded by the independently audited constant-leading
  conductor-six packet, using its direct terminal-degree exclusion in
  Section 7. That endpoint does not rely on a general moving-embedding
  globalization or on the present synthesis.
- `<3,5>` is excluded by the characteristic classification and the
  separately audited direct cubic/quintic coefficient proof.
- `<4,5,6>` is excluded by the characteristic classification and the
  separately audited full-module normal form and quartic/sextic Jacobian
  argument.

The initial hypothesis that one coordinate is monic is retained at every
application. Constant-leading terminal reductions and polynomial source
shears are certified by the actual pair, not chosen from a generic
abstract presentation.

The all-quadratic-semigroup and conductor-six theorem statements and
their independent PASS audits were inspected for this dependency audit.
Their entire proofs are not claimed as newly re-audited here. The
classification, quartic/sextic, and cubic/quintic arguments have separate
full audits in this packet. There is no circular dependency.

No positive delta≤4 survives; together with the zero-defect case this
proves delta≥5 for a nonautomorphism. Multiplying by the exact
conductor identity gives the asserted bound ten.

## 3. Exact dependency and publication-copy checks

All five entries in `SYNTHESIS_DEPENDENCY_HASHES.json` were independently
matched to the corresponding actual proof files:

1. low-delta trace/conductor and symmetric-semigroup proof;
2. constant-leading conductor-six exclusion;
3. all-quadratic-semigroup exclusion;
4. publication copy of the genus-four/quartic-sextic proof;
5. publication copy of the cubic/quintic proof.

All **28** entries of `PUBLICATION_COPY_LINEAGE.json` were checked against
both their source and publication bytes. Exactly six copies differ:

- `CONDUCTOR_EIGHT_CLASSIFICATION.md` omits only administrative staging
  and startup paragraphs, and changes “staging directory” to “packet's
  directory” in the reproduction instruction.
- `CUBIC_QUINTIC_ADDENDUM.md` replaces internal provenance wording with
  the accurate independent-review statement and audit link.
- `review/AUDIT.md` omits its administrative-only final section.
- `review/CUBIC_QUINTIC_AUDIT.md` omits only administrative write-location
  wording, retaining its scope and promotion reconciliation.
- `review/REVIEW_SUMMARY.json` and
  `review/CUBIC_QUINTIC_REVIEW_SUMMARY.json` each omit only the administrative
  `external_mutations` boolean.

The complete diffs were read. No theorem, hypothesis, mathematical
argument, formula, source boundary, exact control, or substantive audit
finding was deleted or changed. All other listed copies are byte-identical.
Original reviewed hashes and the status-promotion receipts remain
available, so the publication hashes have explicit provenance rather than
silently replacing the reviewed evidence.

After the final wording cleanup, the publication classification proof
has hash `493163f9017a20ff467b043773c5b52369793e4696108c0d089d1411ae14f764`.
The synthesis dependency manifest was updated to that exact hash and
independently rechecked. No mathematical statement changed.

## 4. Reproduction and layout

The classification producer/independent checkers retain **57/42** controls;
the cubic/quintic producer/independent checkers retain **31/32** controls.
All four had already been rerun normally and with `-O` in this audit,
with byte-identical outputs matching their recorded receipts. Their
publication copies are byte-identical to those tested programs.

The sextic/nonic realization checker now uses explicit `RuntimeError`
enforcement rather than a bare Python assertion. Its current copied
program was independently rerun in both modes: all **11** controls pass,
the outputs agree, and they match the frozen receipt. The mathematical
realization and conductor formulas are unchanged.

`PUBLICATION_LAYOUT_CHECK.json` records all 18 relative-link resolutions
against existing repository files and the exact approved predecessor
files, as well as relocated normal/optimized reproductions. The target
mapping and evidence were inspected. These checks establish the intended
layout; they do not themselves assert that publication to a remote
repository has occurred.

The evidence is recorded in `PUBLICATION_DEPENDENCY_CHECK.json`.
The lower bound ten is a proof synthesis, not an inference from a bounded
coefficient search or a claim of historical priority.
