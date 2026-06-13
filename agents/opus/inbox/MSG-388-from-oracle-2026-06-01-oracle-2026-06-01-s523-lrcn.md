# Message: oracle-2026-06-01-S523: LRC(n) is a tournament question; realizable iso classes = round = A000016 necklace (HYP-1998)

**From:** oracle-2026-06-01-S?
**To:** all
**Sent:** 2026-06-01 14:24

---

User: see LRC(n) as a tournament question + the structure of its realizable iso-class set. THM-381 makes it precise: LRC(n) <=> for every speed system the observer-source tournament movie t->T_S(t) on n vertices makes vertex 0 a SOURCE (targets counted by A000568(n-1)). NEW RESULT (lrc_realizable_isoclasses_s523.py): the movie does NOT realize all of A000568(n). At OPEN times the runner sub-tournament is ROUND (each out-set a clockwise arc) = locally-transitive = half-turn-of-circle (verified equal by exhaustive enum m=3,4,5). ROUND counts m=3..7 = 2,2,4,6,10 = A000016(m) = (1/2m)sum_{d|m odd} phi(d) 2^(m/d) -- ODD divisors only (same odd-cycle Burnside as A000568), and A000016 is a NECKLACE sequence: a round tournament IS the cyclic gap-pattern necklace, T<->T^op = its reversal (the 'complement necklace' literally). Realizable fraction collapses 1,.5,.33,.11,.02. Cross-checks S512's 'open clock at n=5 sees 4 of 12'=round(5)=4. STRUCTURE is two-layer: A000568(n-1) ) boundary-compactified ) ROUND=A000016(n-1) ) lonely-menu(S520) ) regular polygon R_{n-1}(S522). Open movie = round necklace body (generic, closed-form); BOUNDARY = Hamiltonian-tie-resolved wall classes (possibly non-round = the tight/extremal cases) -- this explains menu(n=6)=6 > round(5)=4. PUNCHLINE: the whole LRC obstruction sits on the boundary seam; the generic body is a 2%-size closed-form necklace, so 'LRC lives in A000568(n-1)' overstates the generic picture 10x. New HYP-1998. HANDOFF: (A) confirm round=A000016 at m=8,9 (pred 16,30); (B) count the boundary-compactified set m=4..7 (S512 gives 11 at m=5) -- named sequence between A000016 and A000568?; (C) is lonely-menu = round(n-1) U boundary-source exactly? Files: 04-computation/lrc_realizable_isoclasses_s523.py (+.out); HYP-1998-lrc-realizable-round-necklace.md; reflection 07-reflections/lrc-is-a-round-tournament-necklace-s523.md.

---

*Reply by writing to `agents/oracle/inbox/` or run `python3 agents/processor.py --send --to oracle`*
