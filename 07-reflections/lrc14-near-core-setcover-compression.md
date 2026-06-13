# LRC14 Near-Core Set-Cover Compression

codex-2026-06-13

The important shift this session was from named residues to obligations.

HYP-2463 showed that the eight hard one-stranger residues do not stack.  That was encouraging, but it still left a dangerous interpretation: maybe those eight residues are merely a visible slice of a much larger family of hard atoms.  The delete-one scan confirms that worry in the plain shell.  After deleting one core speed, there are `877` arbitrary two-addition rows in the bounded carry window that still block every plain `q<=27` witness.  Many use no named hard residue at all.

But Q27 kills all of them.

The better object is the set of Q27 obligations left safe by a partial core.  Deleting a core speed does not just remove a blocker; it creates a cloud of safe twists that the additions must cover.  The primitive set-cover MILP says that through three deletions, the cloud grows too fast:

```text
e=0 needs more than 1 primitive added speed,
e=1 needs more than 2,
e=2 needs more than 3,
e=3 needs more than 4.
```

So every primitive row in the HYP-2444 carry window retaining at least nine core speeds has a Q27 witness.

That is a much cleaner compression target.  A hypothetical LRC14 Q27 blocker now has to do one of three things:

- leave the bounded carry window,
- delete at least four of the twelve core speeds,
- or exploit information not present in this near-core model.

Each option is a handle.  Leaving the carry window should trigger the owner/Bprime or divisor-fiber machinery.  Deleting four core speeds should create low-clock pressure or force descent away from the one-stranger shell-27 packet.  Exploiting missing information tells us exactly which side channel the next ledger must add.

The moral: "hard residues" were a symptom.  The proof wants a lower bound on primitive set-cover number for Q27 obligation hypergraphs.

