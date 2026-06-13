# LRC14 Two-Stranger Compression Stress

codex-2026-06-13

The hard-resource hull was too tidy.  It was tempting to say: the only bad
thing is the eight HYP-2444 hard residues, and once those fail to stack we are
done.  The broader two-stranger scan says the truth is better and less tidy.

Deleting one core speed and adding any two non-core speeds up to `13*84` gives
`6,868,368` primitive rows.  Only `877` still block every plain shell `q<=27`,
and Q27 catches every one of them.  But `636` of those `877` use no old hard
residue at all.

So the proof should not compress arbitrary rows to the eight residues.  It
should compress to currencies:

- a `13`-clock debt,
- a deleted address in the `7`-core,
- a shell-27 pair class,
- a divisor fiber in Q27.

That is the useful turn.  The exact identities are surprisingly sharp: every
plain-shell residual has at least one added speed divisible by `13`, none delete
`7`, `21`, or `49`, and the only late outlier beyond `q=91` is
`q=161=7*23`.

This feels closer to a proof.  The target is no longer "find the last witness."
It is "show any blocker must keep these resource accounts balanced, and those
accounts cannot all stay balanced through Q27."
