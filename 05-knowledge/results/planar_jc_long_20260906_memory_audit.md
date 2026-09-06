# Independent audit of the complete row-fifteen response

**Status: ACCEPTED; PROVED relative to the inherited row chart, with an
independent exact reconstruction.** The proof in
[the response report](planar_jc_long_20260906_memory.md) passes the type,
quantifier, full-universe, affine-action, and consequence checks below.

The source starts with all **128** polynomial coefficients in rows twelve
through fifteen of `(dA,dC)`, using the explicit first four background rows
from THM-4308 with `Phi=0,Delta=896/15`. It imports no repository mathematical
source, row solver, tangent parameterization, or frozen output. It forms
the literal differential of `P=C^2-A^3+3A/4+1/4` and the actual Jacobian.
The depth equations come from the independently proved complete diagonal
bank in [the depth report](planar_jc_long_20260906_depth.md), not the
producer's generator-matrix nullspace. The two paths consequently differ
both in their row variables and in their depth representation.

The original matrix has 78 coefficient rows from `P`, 70 from the
Jacobian, and all 91 depth constraints. Exact rational row reduction gives
tangent rank118, augmented rank121, and exactly the producer's three
source relations. All 239 original equations vanish after a newly computed
lift. The weight20/21 ranks are1 and2 after adjoining the high column;
the weight22/23 ranks are already3 and remain3. Thus this is a full
compatibility and minimality audit, not a check of a selected scalar.

The report's filtered linearization is valid: with perturbation order12,
quadratic `P` terms start at24 and quadratic Jacobian terms at23. The
finite equations through `P` row15 and Jacobian row14 are exactly affine
in the perturbation and depend only on background rows0..3. The tangent
description in the producer spans all capped solutions because the
boundary derivative entries are coprime. No nonlinear equation is dropped
inside the asserted range. The nonzero quadratic term at24 is retained
as a stopping boundary; it is not an all-row obstruction to new repairs.

There are exactly eight low-weight source monomials in the declared
valuation>=13, weight<=23 universe. Source monomials at the same valuation
have distinct leading x powers, so omitted earlier valuations cannot cancel
each other into this universe. Five odd source directions are neutral;
at weight<=22 only two remain. For a nonzero high response, the forced
`p^5*y^4` coefficient proves the scoped threshold22. A high response that
vanishes over an extension field is correctly exempted from minimality.

The `G_m x A^2` consumer uses the actual THM-4438 positive family and the
fixed first-four-row hypothesis. Its two odd coefficients are independent
of the inherited source parameters, and the terminal rank is constant10.
The weight22 transport is credited to incoming
`d0208173b2 / planar_jc48_sep06_weight14.md`; its independent recovery is
not counted as a second new transport theorem. Arbitrary lower valuations,
termination and chart entry remain open.

The script uses exact `QQ` row reduction through SymPy's DomainMatrix.
An initial run of the generic symbolic row reducer was stopped for excessive
intermediate rational growth; switching to the exact-domain implementation
changed the arithmetic engine, not the matrix or the tested predicate.

```text
python3 -B 04-computation/planar_jc_long_20260906_memory_audit.py
python3 -B -O 04-computation/planar_jc_long_20260906_memory_audit.py
```

Both modes byte-match the frozen
[output](planar_jc_long_20260906_memory_audit.out), with **248** live gates.
Source/output bytes are pinned by the session manifest. No claim of a Lean
formalization, arbitrary source minimality, or planar JC resolution is made.
