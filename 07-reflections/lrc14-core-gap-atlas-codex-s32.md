# LRC14 Core-Gap Atlas Reflection

**Source:** codex-2026-06-19-S32, HYP-2651 / T898.

The main gain is a correction in proof hygiene.  HYP-2638-style Freiman
normalization is the right language for sector functionals where translation
has already been quotient-preserved.  THM-523's `G_C` crux is fixed-observer
gap measure, and translating a positive core changes that predicate.  So the
proof has to keep exact positive cores first, then route to Freiman or
state-word templates only after preserving the predicate.

The exact atlas through `B=19` is encouraging for the uniform fattening lemma:
`50,388` primitive positive 12-cores all satisfy
`meas(G_C) >= 7/858`, with equality only at
`(1,2,3,4,5,7,8,9,10,11,12,13)`.  The next distinct row is the endpoint drop
`(1,2,3,4,5,6,7,8,9,10,11,13)` at `426/35035`, separated from the collar by
`841/210210`.

The single-hole AP-window ledger is the first proof target.  Prove directly
that `[1,13]\{6}` has the smallest safe measure among `[1,13]\{e}`.  Then use
the state-word machinery to prove that anything below `426/35035` must have
the same drop-6 template.  That would turn the bounded near-collar problem into
one explicit wall calculation plus a template-exclusion lemma.

The tail rows that enter near the top through `B=19` are not random strangers.
They remain AP-cluster rows with a few holes, often still carrying the central
`6` hole.  This is the right place to combine HYP-2648 state words with
HYP-2644 far-element plateau: AP-cluster tails need a template compression
lemma, while genuinely dissociated tails should be handled by plateau
decorrelation.

Tournament-analysis lesson: the vertices should be proof obligations, not
runners.  The winning path is exact core gap, AP drop profile, state-word
template, sumset-excess classifier, far-element plateau, raw runner bound.
