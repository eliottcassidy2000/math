# Earlier valuation-twelve source response: exact scout

**Status: FINITE-EXACT SCOUT; independent proof audit pending.** This extends
the declared linear universe of the audited
[valuation-thirteen classification](planar_jc_long_20260906_memory.md).
The full statement and global boundary consumer will be promoted only after
the independent reconstruction. JC(2) remains open.

Allow perturbations in unknown rows11..15 and source valuation>=12.
The first five background rows are fixed to the actual THM-4438 boundary;
the newly used row is

```text
A4=224x^4/9-256x^2/75-37376/405,
C4=-184x^5/15-2048x^3/25+98944x/675.
```

The script retains all 155 unknown polynomial coefficients and all literal
`P`/Jacobian equations, then all91 depth relations. The complete source
bank contains14 monomials of valuation12..15 and weight<=23, and the
designated weight24 high channel. Exact rational elimination gives tangent
rank145 and source rank5. Weights<=21 still cannot replace the high channel;
weight22 can. The replacement solution at weight<=22 has dimension5,
with four odd directions and one even direction. This is an exact matrix
consequence, not a proof of global source entry or termination.

The actual source coordinates, five raw relations, and replacement fibres
are printed in the [output](planar_jc_long_20260906_memory_valuation12.out).
The exploratory [source](../../04-computation/planar_jc_long_20260906_memory_valuation12.py)
imports no repository mathematics and builds the literal differential from
the fixed background. The separate audit uses the inherited tangent row
operator and generator depth matrices instead.

```text
python3 -B 04-computation/planar_jc_long_20260906_memory_valuation12.py
```

Quadratic terms begin at `P` row22 and Jacobian row21, beyond the tested
rows. This identifies the exact range where an affine interpretation could
be promoted. Valuation11 has a new background coefficient and is a separate
scout; its hypotheses must be retained, not specialized silently.
