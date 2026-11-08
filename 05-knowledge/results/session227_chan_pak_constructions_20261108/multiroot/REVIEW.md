# Independent audit and status promotion

The independent audit (source workspace: `session227-construction-audit`) returned **PASS** for all claims, including the exact minimal rational linear gadget summary (§3.1) and rooted-graft-indistinguishable strong collision family (§5).

The auditor ran 170,515 explicit finite checks, including 12,078 literal grafts, 10 full-state collision graphs through order 15, 25 literal response entries, and response ranks 2, 3, 4, 5. Isolated normal and optimized outputs byte-match. Independent output SHA256: `931818bee93e8198838c6861eea026c2793a1c3fbaed9723aed5999c9f072626`. Independent script SHA256: `e8fcef992f937252be92f89e4321bdc59378df144f571c97594abc45754863fa`.

After that signoff, the only changes to the audited proof and README were replacing their pending-audit status lines with independently audited PASS. No theorem, proof, scope, implementation, or receipt changed. This review note and refreshed producer hashes record the promotion. The audit reported no unresolved mathematical blocker.

The producer normal and optimized receipts remain byte-identical, SHA256 `fdae4b708949c956cf0e27c509f8b4c28f9e077d62715f93658f6523988cf981`.

The initial unequal-shore coefficient overstatement and its prepublication correction remain documented in README.md; the audit reviewed the corrected lower-bound statement and exact equal-shore formula.
