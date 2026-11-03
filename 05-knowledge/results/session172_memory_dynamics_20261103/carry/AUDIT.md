# Independent audit of carry read heads and positive contextual states

**Verdict: PASS. No mathematical correction required.**

Audited NOTE.md and check.py independently, then reran the author's checker normally and with `-O`. Both outputs match results.json byte-for-byte. The author receipt reports 105,275 checks.

## Proof checks

- Upper-triangular diagonal projection is multiplicative and preserves existing inverses. Hence all declared triangular word traces are carry-blind, even with repeated unknown matrices, inverses, and fixed two-sided contexts.
- Direct product computation gives tr(PG_b)=p11*a+p21*b+p22*d. Over a field, injectivity on a source containing at least two carry values is equivalent to p21 nonzero. J reads b, including in characteristic two. The ring caveat is necessary and stated.
- The external J experiment is explicitly not claimed to be a legal Collatz step, source-certificate operation, or physical tournament realization.
- For positive channels, an extreme deterministic output forces every positively weighted column to agree. Thus supports of sources from distinct contextual classes are disjoint. This proves h>=q without linearity or injectivity of initialization. Class point masses prove abstract attainment. The separate requirement of congruence for deterministic dynamical attainment is correctly stated.
- The lower bound counts simplex states, not bits or arbitrary signed dimensions. The 32/2048 versus 5/11 distinction is correctly preserved.

## Independent computational controls

`audit_independent.py` imports no author routines. It checks 134,067 matrix/probe configurations over F2, F3, F5 and F7, including all diagonal pairs and all 2x2 probes. It verifies the lower-left injectivity criterion and J readout. It independently solves exact odd two-step valuation guards modulo 32, obtaining only 19 for (1,3) and only 29 for (3,1).

Receipt: `audit_independent.json`. Author replay receipts: `audit_observed.json` and `audit_optimized.json`.

## Source scope

Read the exact supplied carry source excerpt; it supports inherited affine matrices, chronology, carry cocycle, trace blindness, and scoped faithful representation attribution. Read the published S33 source excerpt at the permanently linked S33 source, especially section 2a: it explicitly establishes the inherited 38-positive-state versus 19-signed-linear-dimension order-five example and pure-source condition. Read the guard source at the permanently linked information-observers source, especially equations around lines 95–117: it explicitly supplies the two exact guards and source-verification distinction.

The S171 framework/information/synthesis packets had already been read independently in this audit session. Their RM dimension, 17/19 readout, and contextual-class facts support the stated application. No independent re-enumeration of S33 or S171 tournaments is claimed. The note correctly marks these as inherited and makes no external priority or convergence claim.

For final publication, retain permanent repository source paths/links and exclude bulky fetched historical source copies unless needed for the archive. No forbidden history file was read or scanned.
