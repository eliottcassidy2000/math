# Minimal historical correction text for parent integration

No source repository was edited. These are proposed insertions into the parent's
owned publication worktree. Preserve the historical body and original output;
link the correction to the Session67 final integration path selected by the parent.

## HYP-3602

Target:
`05-knowledge/hypotheses/HYP-3602-tournaments-are-intransitivity-among-n-things.md`.

At source pin `b4234f728c0516b5fc521c001d69e26f3d2508d5`:
- blob `c3a740953e8a4196e28b9af5c3844089d4a78694`
- SHA-256 `e0be287c02d7f93b93f9c8d10b4592576094ca7fed88066f5e9dae9794d8ebc9`

Proposed frontmatter status replacement:

```yaml
status: >
  HISTORICAL REFRAME / PARTIALLY CORRECTED (Session67, 2026-10-26).
  Relation-first interpretation and the triangle obstruction to transitivity survive.
  Literal zero-curl/cohomology and unlabelled order-times-shape identifications do not.
  See the current-scope correction immediately below the frontmatter.
  No LRC theorem or reduction is established by this historical reframe.
```

Insert immediately AFTER the closing frontmatter delimiter, before the existing
title; do not move or delete any later historical prose:

```markdown
> **CURRENT-SCOPE CORRECTION — Session67, 2026-10-26.** This is historical
> idea provenance, not a complete invariant or an LRC proof. Binary
> transitivity means `S_ij=sign(p_i-p_j)` for an injective ranking, not
> `S_ij=p_i-p_j`: a transitive triple already has signed curl `+1` or `-1`.
> The exact cyclic-triple indicator is `(c_ijk^2-1)/8`; zero real curl and
> binary transitivity must not be identified. The full clique complex has
> zero first cohomology, and GF(2) switching/cycle coordinates are a different
> quotient (switching can turn a transitive triple into a cyclic one).
> The labelled count identity `2^binom(n,2)=2^(n-1)2^binom(n-1,2)` does not
> factor A000568 into unlabelled order and shape counts; HYP-3595 already
> warned that stabilizers couple the factors. The repaired lossless family
> datum is a pair-free tournament core with positive transitive-module sizes,
> modulo its automorphism group. Alternatively, a cyclic-triple hypergraph
> needs its full realization orbit under its own automorphism group.
> “Only quadratic invariant” is not “only invariant of any degree.”
> The danger-relation and odd-cycle language below establishes no LRC
> implication; LRC(14) remains open. See Session67's information recovery
> note, §§2, 4–6, and Session=36's exact curl-realization results.
```

Dependencies and checks:
- HYP-3595 exact source path and hashes are S06 in `SOURCE_MANIFEST.json`.
- Hodge identities are directly proved in the recovery note and do not depend
  on the unaudited portions of THM-4739.
- Theorem 4.1 (maximal transitive-module partition and weighted-core uniqueness)
  and formula 5.2 (hypergraph realization orbits) have complete proofs there.
- All labelled n<=5 cases, all possible C3 inputs through n=5, and switching
  and source/sink hostiles independently replay with normal and optimized Python.

## March-header information script and stored output

Targets:
- `04-computation/information_theory_deep.py`: blob
  `49d118b366dd4fbc0da709f6ed6697b36b0513d1`, SHA-256
  `144513d4bc10f62ab8203e4666805d0100533e7f886696dc6cf497c731ef7f7d`.
- `05-knowledge/results/information_theory_deep.out`: blob
  `c94de2106516d62edb1defa9c8789e25bdd83166`, SHA-256
  `e7db503deb30b3ce537dbc312276087576d16286fe7ccb4bf80cae02de341a96`.

Both hashes are at the same fixed source pin above. The script's own header
says March 13, but the accessible file archive commit is June 12; do not convert
that internal header into an independently authenticated creation date.

This is not a speculative diagnosis: lines 128–130 explicitly claim
`H(pi|T)=E_T log2 H(T)` because “each consistent ranking [is] equally likely”;
line 141 sets `MI=H_pi-expected_log_H`, and line 145 prints it as `I(pi;T)`.
Hamiltonian paths only constrain consecutive comparisons, and the noisy
posterior on all rankings is neither this set nor its uniform distribution.

Proposed script-docstring warning (directly after opening triple quote):

```text
CURRENT-SCOPE CORRECTION (Session67, 2026-10-26): historical exploratory
script. Section 2's quantity labelled I(pi;T) is only the surrogate
log2(n!)-E log2 H(T), NOT channel mutual information. H(T) counts directed
Hamiltonian paths (consecutive-edge consistency), not rankings satisfying
all comparisons or posterior support. At comparison accuracy 1/2, true
I(pi;T)=0, while the n=3 surrogate is 2.188721875540867 bits. Section 5's
claim that Redei's theorem makes all tournaments strongly connected is
also false (the transitive tournament is a counterexample). The original
computations are retained as historical diagnostics; see the Session67
information-recovery note for repaired statements and independent tests.
```

A runtime warning before the first computation is recommended so rerunning
this historical script cannot silently revive the obsolete labels.

Proposed stored-output banner (byte zero, before the old output):

```text
HISTORICAL OUTPUT — CURRENT-SCOPE CORRECTION (Session67, 2026-10-26)
Section 2's I(pi;T) values are path-entropy surrogates, not mutual information.
Hamiltonian paths do not enforce every pairwise comparison; at 50% comparison
noise the true ranking/output mutual information is 0. The Section 5 universal
strong-connectivity claim is false. Original output follows unchanged for
provenance; corrected identities and hostile checks are in Session67's
information-recovery note. This banner is not a rerun receipt.

```

The exact valid posterior is

```
P(pi | T) = p^agreements(pi,T)(1-p)^(m-agreements(pi,T))
            / sum_sigma p^agreements(sigma,T)(1-p)^(m-agreements(sigma,T)),
```

for a uniform ranking prior, interpreted by its limiting values at endpoints.
True mutual information requires its conditional entropy. The n=3 p=1/2
hostile is an exact argument, not a small floating-point discrepancy.
