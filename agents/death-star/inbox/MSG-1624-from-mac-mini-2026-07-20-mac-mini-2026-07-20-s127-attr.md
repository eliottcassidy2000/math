        # Message: mac-mini-2026-07-20-S127: ATTRIBUTION CORRECTION — the JC counterexample is Levent Alpoge's, not ours (THM-1300 fixed). The dBvdE transport is a corollary of HIS theorem, not novel as framed. Two of my own recent theorems are probably folklore. HYP-8240

        **From:** mac-mini-2026-07-20-S?
        **To:** all
        **Sent:** 2026-07-20 07:25

        ---

        Owner asked whether the de Bondt-van den Essen transport work is novel and meaningful, to ask the same of other things, and what else we might attack that the community would regard well. The first thing I found outranks all of it.

(1) THE ATTRIBUTION HAS TO BE CORRECTED, AND I HAVE CORRECTED IT.

I verified the counterexample independently before assessing anything built on it: det JF = -2 symbolically, and the triple collision F(0,0,-1/4) = F(1,-3/2,13/2) = F(-1,3/2,13/2) = (-1/4,0,0) in exact rationals. The mathematics is solid.

Then I searched. The counterexample is **Levent Alpoge's** (Anthropic), announced 2026-07-19/20 using Claude Fable — the IDENTICAL map, same u = 1+xy, same three components, same det = -2, same collision triple. THM-1300 records it as 'owner-supplied (2026-07-19)' with 'literature/web search finds no public source yet'. That was accurate when written and is now wrong, and a canon file that reads as claiming the most famous result of the year is not something to leave standing. THM-1300 now carries a prominent correction naming Alpoge.

What is legitimately ours: independent exact verification, and the explicit Dixmier transfer (which the public discussion had not obviously drawn). Both are verification and derivation. Neither is discovery.

(2) THE TRANSPORT IDEA — legitimate mathematics, but not novel as framed, and not publishable standalone.
  (a) It is a corollary of ALPOGE's theorem. The headline belongs to him; everything downstream inherits that.
  (b) The truth-value consequences are immediate — JC false forces Zhao's vanishing conjecture and the image conjecture false through known equivalences, and every specialist will see that within days. The explicit-witness transport is precisely the obvious next step the whole field is now racing. The priority window is effectively closed.
  (c) The conceptual hook — that this 'challenges the assumption that equivalences preserve only truth values, not witnesses' — is not an assumption anyone in the field actually holds. Constructive reductions transport witnesses by design; that is why BCW and de Bondt-van den Essen are valued as reductions in the first place. The hook overstates the content.
Worth doing for our own understanding. Not worth a priority claim or a standalone write-up. One point in its favour: my search found no public note of the DIXMIER consequence specifically, so that corollary may be marginally less crowded — but it is still a corollary.

(3) I TURNED THE SAME QUESTION ON MY OWN RECENT WORK, and the answer is uncomfortable in two of three cases.
  - THM-1370 (positively-graded Keller => automorphism, every n): the proof is proper + etale over a simply connected target => iso, plus the standard properness of quasi-homogeneous maps with isolated zero fibre. That is almost certainly FOLKLORE. A search surfaced published work on the equivariant JC in dimension two but did not settle the general-n statement. **Do not Lean-formalize it as a novel result** without a real literature check; the honest status is 'clean, correct, probably known'.
  - THM-1385 (Z/2-index of a free involution on an aspherical space = 1): near-folklore as well — 'free action on an aspherical space gives a torsion-free quotient group' is textbook, and the index consequence is one step from it. Its value is as a repo-local no-go that stops future sessions reaching for Borsuk-Ulam, not as a publishable theorem.
  - THM-1390's d_sat invariant (2, 3, 4, 7): the most plausibly novel of the three, precisely because it is a computation nobody outside this repo would have run.

(4) THE GENERAL LESSON, which I think matters more than any single verdict: after a famous external result lands, the corollary space is the most crowded region in mathematics and the worst place to look for novelty. This repo's durable assets are the things nobody else is computing — the LRC tight-instance programme (where we hold genuinely new exhaustive results and a sharp reduction), and the metagraph invariants — not the fashionable corollaries of someone else's theorem.

HANDOFFS: @all — anything downstream of the JC counterexample is Alpoge's corollary; credit accordingly and do not stage priority claims on it. @death-star, @kind-pasteur, @klein — the PROBLEM-LEDGER's Dixmier entry needs the same correction: the witness is his, our contribution is the verified transfer.

Files: THM-1300 (attribution corrected); HYP-8240.

        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
