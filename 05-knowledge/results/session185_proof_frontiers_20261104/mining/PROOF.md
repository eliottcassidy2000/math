# A negative-constant barrier to ordered-context continuation

2026-11-04. **PROVED:** a general obstruction to continuing a marked suffix by adjoining arbitrary ordered contexts. Applied to the new Session185 suffixes at addresses24 and44, neither next digit can be obtained by this kind of extension, even with arbitrary fixed clock shifts. This does not obstruct the replacement witnesses supplied by the closure packet or internal edits and module inflation.

## Inherited setting

Use Session182's actual seven-coordinate signature z=(h,p,F,G,J,e,s), c=h−J, and ordered-sum law. The full root attachment has H_t=c+2^t F. A uniform digit d requires c_target=3c+d, independently of the clock shift. The seven-coordinate law and this necessary condition are inherited, not claimed as new.

For two marked tournaments X and Y, their ordered sum satisfies

    c_(X→Y)=c_X h_Y − h_X J_Y − e_X s_Y − G_X F_Y
             =h_X c_Y − J_X h_Y − e_X s_Y − G_X F_Y.

All signature coordinates are nonnegative integers. Every tournament has a positive odd Hamilton count. Moreover H=1 exactly for transitive tournaments: a directed triangle gives at least three paths on that triangle, and the inherited insertion monotonicity extends this lower bound to the whole tournament (equivalently the odd-cycle formula contributes1+2 for that cycle), while triangle-free tournaments are transitive.

## General theorem

Let B be a fixed marked tournament with c≤−3 and

    min(F_B,G_B)>−2c−2.

Then no marked ordered extension X→B→Y has constant3c+2 or3c+4. Either context may be empty; the marking on B is retained and the new context markings are arbitrary.

**Proof.** Applying the displayed identity twice and discarding nonnegative subtractions gives

    c_(X→B→Y) ≤ h_X h_Y c.

If either context is nontransitive, h_X h_Y≥3. As c<0, the new constant is at most3c, strictly below either requested target. Consequently any successful contexts must both be transitive.

Now h_X=h_Y=1. If X has a port, G_X≥1, so the left-attachment subtraction is at least F_B. Appending Y cannot increase the negative constant. Thus c_target≤c−F_B<3c+2≤3c+d.

If X has no port but Y has one, F_Y≥1 and G_(X→B)=G_B. The right-attachment subtraction is at least G_B, again giving c_target≤c−G_B<3c+2≤3c+d.

If neither context has a port, both are unmarked transitive tournaments, their J,e,s,F,G coordinates vanish, and c_target=c. Since c≤−3, c>3c+d for d=2,4. All cases are excluded. ∎

No assumption about the size, SCC order, or marked structure of X,Y is used. Strongness of the resulting full attachment need not be assumed because the necessary constant equation already fails. Fixed positive or negative clock shifts cannot change that equation.

## The two escaping suffixes

The new literal witnesses supplied by the companion closure packet are:

- Address24: rows(2,20,9,3,45,79,31), port mask42, signature(35,676,252,267,70,13,14). Here c=−35; both F252 and G267 exceed68. The target constants are−103 and−101.
- Address44: rows(2,4,33,23,71,155,47,95), port mask52, signature(89,2144,504,627,118,17,24). Here c=−29; both F504 and G627 exceed56. The target constants are−85 and−83.

Therefore both literal carriers require something other than adding ordered context if their next requested digits are to be realized. Their existence and the next-digit replacement witnesses are credited to the closure packet. This lemma is only the new all-context obstruction.

The first carrier has the unique proper nontrivial homogeneous module{0,3}, with3→0, port3 and nonport0. In the order(5,6,4,2,3,0,1), its three back arcs are4→5,1→4,1→2; exhaustive permutations establish minimum feedback arc number3. Inflating this mixed-port module or changing its internal arc arrangement is a concrete next research direction, not a construction proved here.

## Source and novelty check

Current-main S184.md was read directly via GitHub on2026-11-04; it still states only the necessary SCC-order≥5 boundary. Targeted GitHub code searches indexed at ad09a24c8fd25eca9ea320d1b240c04990ada96b found no exact occurrence of the address24 row tuple or signature substring252,267,70 (both with and without spaces). This is not a proof of historical novelty. Earlier S180 phases already supplies other clock-uniform paired constructions, so no first-ever uniform two-digit claim is made.

Prior related structural tools: THM-4748's increasing-run count for a bounded number of back arcs and its three-back-arc Fibonacci fan; S182's ordered-signature law; S184's SCC≤4 return-tree theorem. The suggested module-inflation route is distinct from adding ordered context.
