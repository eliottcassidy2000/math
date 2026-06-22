# LRC14 Bonferroni needs the missing-depth parity guard

*codex-2026-06-22-S115*

The corrected Legendre/Venn picture made the far-runner Newton packets look
like a natural Bonferroni-3 upper bound:

```text
p0(B union F) <= p0(B) + T_1 + T_2 + T_3.
```

That statement is too coarse.  The right object is pointwise and depth-labelled.
At slow time `x`, let `M(x)` be the set of inner sectors not hit by the base
`B`, and put `d(x)=|M(x)|`.  If far runner `i` has sector color `c_i(x)`, then
the far runners cover the base-missing sectors exactly when

```text
prod_{a in M(x)} OR_{i : c_i(x)=a} z_i = 1.
```

Expanding this product gives the r-far Newton packet

```text
T_r(x) = (-1)^(r+d(x)) C_{d,r}(x),
```

where `C_{d,r}(x)` counts r-subsets of far runners whose colors all lie in the
missing sectors and cover all of them.  After integrating over `x`,

```text
T_r = sum_d (-1)^(r+d) A_{d,r}.
```

So the `r>=4` tail is not controlled by Venn containment alone.  It is a
ledger of positive even-depth high packets against negative odd-depth high
packets.

Exact checks make this unavoidable.  The KPS binding-style row

```text
B={0,1,2,3,4}, F={20,21,22,23,24}
```

has negative tail `-13141/212520`; its depth-3 negative contribution dominates.
But

```text
B={0,1,2,3}, F={16,28,29,32}
```

has positive Bonferroni-3 tail `19/68208`, because the positive depth-4 packet
beats the negative depth-3 packet.  More importantly, activation depth alone is
not enough: the edge-active row

```text
B={0,1,7,10,13}, F={15,23,24,31}
```

has `T_2>0` and still has positive tail `307/598920`.

The proof target is therefore:

```text
positive even-depth high packets
  <= negative odd-depth high packets + cap slack.
```

If a class violates that internal domination, it is not a failure of the LRC14
cap route; it means that class must be routed through the existing
doublet/triple/decorrelation atlases rather than through a universal third-order
truncation.

Tournament Analysis for this session uses proof carriers as vertices:
`pointwise_depth_formula`, `depth_parity_tail_inequality`,
`doublet_plus_triple_cap`, `bonferroni3_unconditional`,
`raw_venn_containment`, and `raw_runner_vertices`.  The transitive Hamiltonian
path ranks the pointwise formula first and raw runner vertices last.  The
challenged assumption is that runners, or even first live far-packet depth, are
the right vertices.  They are not; the preserved predicate is the
missing-sector depth-labelled Newton packet.

Artifacts:

- `04-computation/lrc_bonferroni_depth_guard_codex_s115.py`
- `05-knowledge/results/lrc_bonferroni_depth_guard_codex_s115.out`
- `05-knowledge/hypotheses/HYP-2903-lrc-depth-corrected-bonferroni-obligation.md`
