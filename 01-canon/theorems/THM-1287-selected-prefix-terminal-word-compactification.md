---
id: THM-1287
title: SELECTED-PREFIX TERMINAL-WORD COMPACTIFICATION
status: RESERVED / PROOF IN PROGRESS (the selected-prefix component argument and endpoint-alternation caps have been independently derived; exact referee, paper proof, and Lean arithmetic consumer are still being assembled).  This namespace is claimed before the concurrent proof packet is completed
source: codex-2026-07-19-S82 endpoint-word continuation
depends_on: [THM-1233, THM-1253, THM-1275, THM-1283]
related: [THM-1274, THM-1277]
---

# THM-1287 -- selected-prefix terminal-word compactification

Reserved for the following derived statement.  If a globally
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

The table and every strict/weak endpoint convention remain to be frozen by
the exact referee and formal arithmetic consumer.  This theorem will bound
the selected **word**, not the carrier, the normalized speed packet, or the
phase/address stalk; it is not yet a finite enumeration and does not prove
LRC(14).
