# Complete carry tests and a source-selected exact continuation

Status: **FINITE-EXACT complete classification in the stated grammar**, plus
**PROVED** recurrence and carry-transport identities. Global closure is OPEN.

## C1. Universe and the all-order cutoff

The parent must be a nonempty word in the original65 alphabet augmented by
s78. The return is a nonempty word in the786-block catalog monoid. For each
of the seven supplied input states v, seek both j=1,2 with

    u=3v+j t+7(0,cx,cy)=v(parent) P_return.           (C1)

There is NO bound on cx,cy or on word length. H is fixed at3h+2j, and all
F49 residue/charge data are retained. We minimize |cx|+|cy|, not graph order.
The nonempty-return requirement is essential: the inherited exact31->95/97
seed replacements remain valid, but do not supply such a return.

For a block P with first column(p0,p1,p2), the physical cone implies

    H_after >= [p0+min(p2,(p1+p2)/2)] H_before.

Minimize p1*x+p2*y over x>=0, y>=x, x+y>=h to obtain this bound. Every
catalog block has multiplier at least5. Absorb all but the final catalog
block into the parent; it remains a66-bank word. The largest target is889,
so the last parent must have H<=177. Monotone forward closure has exactly939
states including the formal empty root, which is excluded as a parent.

The proof of completeness is inherited norm descent plus forward closure:
H cannot leave this bound and reenter; the physical C16 bounds are invariant
also for s78. An independent extended inverse reconstructs every state.
All423 distinct return matrices are tested with the bound respected. There
are1444 feasible state/matrix evaluations. Repeated matrices are merged for
this counting-state existence test, with an actual word retained for each.
They are not declared identical graphs or minimum-cost presentations.

## C2. Exactly two of fourteen branch requests survive

| Input state | j | Output state | Minimum carry (0,cx,cy) | Parent word; return |
|---|---:|---|---|---|
|(97,40,75)|1|(293,231,295)|(0,16,10), norm26|a r41; b|
|(293,210,281)|1|(881,692,885)|(0,9,6), norm15|b r25; b|

The other twelve branches are absent in this grammar for EVERY carry, not
just for small carries. Neither minimum is a minimum vertex-order assertion.
This is absence of a return-compatible source presentation, not absence of
the scalar H values from the unrestricted tournament spectrum.
The inherited H97->293 repair with carry(0,13,8), norm21, remains valid in
its larger operation language; its output has no permitted terminal return
in the grammar used here. Our norm26 is not a minimum over arbitrary repairs.
The explicit outputs have9 and11 vertices and designated ports8 and10.
Their retained descending spines prove positive port-start feasibility.

The native compiler transforms the actual supplied inputs into these outputs
using14 and22 legal path-complement moves, with peak orders10 and12.
An independent graph-level replay checks every path, birth, temporary
deletion, and endpoint. These are source replacements changing old arcs,
not induced extensions. All old vertex names survive these two programs;
the input and output ports and the padding/temporary clocks are recorded.
Intermediate H and strongness are not asserted to be preserved.

## C3. Selecting the right H293 representative removes the next carry

The new H293 state differs from the previously supplied one by7(0,3,2).
They share H and the full residue data, but have different literal sources
and vertex clocks. For the new representative,

    (881,692,885)=3(293,231,295)+t.                  (C2)

Thus the source-selected chain is a genuine H97->293 carry replacement
followed by an exact H293->881 replacement. The latter has its own replayed
22-move program, peak order12. No representative is silently substituted.

The elementary re-encoding rule explains the cancellation. If an edge has
carry kappa, choosing actual representatives v+7eta and u+7xi changes it to

    kappa' = kappa + xi - 3eta.                      (C3)

Here eta=(0,3,2), xi=0, so the old carry(0,9,6) disappears. A formal choice
of eta,xi need not have a tournament representative. On a tree one can
formally eliminate carries recursively; the earlier unrealizable exact
affine states show why this algebra is not a realization theorem.

For unrestricted H-preserving multiples-of7 carries, acceptance of a target
in this fixed grammar depends only on H and x,y modulo7. The integer carry
cost still depends on the actual representative. Consequently the second
H293 branch to883 remains absent after this reselection, so the chain does
not close the ternary family.

## C4. An infinite actual offset family recovered from the positive branch

Let A=`a r41` and B=`b r25` in processing order. Their states satisfy

    v(A)=(37,29,35), v(B)=(109,88,105)=3v(A)-t.

Since tP_q2=t and tP_b=-t, every k>=0 gives

    H(B q2^k b)=3H(A q2^k b)+2.                       (C4)

These are actual words with vertex clocks9+4k and11+4k. Put M_k=H(A q2^k b).
The exact operator identity

    (P_q2^2-18P_q2+5I) P_b e_h=12(1,0,-1)^T

and the conserved charge h-y=2 give

    M_(k+2)=18M_(k+1)-5M_k+24,
    M0=293, M1=5223,
    sum M_k z^k=(293-344z+75z^2)/((1-z)(1-18z+5z^2)).

The first pairs are(293,881), (5223,15671), (92573,277721),
(1640223,4920671), independently counted from adjacency rows. The new family
is a transported one-offset identity, not an all-branch induction or a
prime-producing theorem.
