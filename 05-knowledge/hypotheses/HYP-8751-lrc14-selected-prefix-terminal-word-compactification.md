---
id: HYP-8751
title: SELECTED-PREFIX TERMINAL-WORD COMPACTIFICATION
status: DERIVED / COMPUTER-EXACT / FORMALIZATION-PENDING (the selected-prefix component argument and endpoint-alternation caps have been independently derived, and the exact recurrence referee passes normally and under Python optimization; a full paper proof and Lean arithmetic consumer remain to be assembled)
source: codex-2026-07-19-S82 endpoint-word continuation
depends_on: [THM-1233, THM-1253, THM-1275, THM-1283]
related: [THM-1274, THM-1277]
script: 04-computation/lrc14_selected_prefix_terminal_word_compactification_hyp8751.py
output: 05-knowledge/results/lrc14_selected_prefix_terminal_word_compactification_hyp8751.out
script_sha256: 756946cb036761f53cb43d90e7d2c24e46aadc8076e55cf2f48e19fcfda7b387
output_sha256: 0ba1fcf3f060dcbda6b96a44df8ffae2fc0da9da60a5e249922ca043e49b9c83
---

# HYP-8751 -- selected-prefix terminal-word compactification

The following statement is derived.  If a globally
single-occurrence terminal endpoint owner is `d_r`, rerun THM-1233's
prefix-survivor/component-span proof after deleting only the **selected**
prefix teeth.  The same survivor mass remains, while the component count is
`1+sum n_i` in the selected occurrence counts.  With `n_r=1`, this yields
uniform selected-word caps.  Since the terminal occurrence is at one end of
the chronological word and adjacent selected teeth have different owners,
the fastest count is at most the non-fastest count.

The currently derived target table is

```text
r     selected caps (n_1,...,n_5,K)   word cap   d_6/c cap
2     (2,1,15,54,165,237)                 474       1659
3     (2,7,1,31,95,136)                  272        952
4     (2,7,7,1,41,58)                    116        406
5     (2,7,7,7,1,24)                      48        168
6     (2,7,7,7,7,1)                       31          7
```

The exact referee freezes the table and all strict/weak arithmetic ceilings;
the formal arithmetic consumer is still pending.  The result bounds
the selected **word**, not the carrier, the normalized speed packet, or the
phase/address stalk; it is not yet a finite enumeration and does not prove
LRC(14).
