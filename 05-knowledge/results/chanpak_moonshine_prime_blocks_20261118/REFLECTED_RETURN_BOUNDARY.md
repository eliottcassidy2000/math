# Equal-gain recontact does not yet commute with the free reflected seam

## Source and definitions

The earlier packet is `05-knowledge/results/chanpak_wheel_square_selector_20261117`,
commit `43a77c33ac335958ce32b234b737cc0e9ef0a6a0`, files FREE_SEAM.md,
CONGRUENCE.md, REFLECTION.md and check_free.py (read directly).
The current physical pair is from `niche_bridge/EQUAL_GAIN_RECONTACT.md`.
No new priority for the inherited flank law or arc-reversal identity is claimed.

For a central tournament with a marked apex, let D be its deletion tournament,
mask mark the deletion vertices beaten by the apex, h=H(D), P=P_2(D), and
kappa the full central count. If E_i(D[S]), S_i(D[S]) are endpoint/start
counts and H(empty)=1, define L_even/L_odd by summing
H(D\S) sum_{i in S, mask_i=1} S_i(D[S]) over even/odd |S|.
Define R_even/R_odd similarly using E_i and mask_i=0.
L=L_even+L_odd, Lminus=L_even-L_odd, and likewise R,Rminus.
The seven-coordinate profile is (h,P,L,R,Lminus,Rminus,kappa).
The inherited exact read is

 H(u,w)=kappa+h(F(u)+F(w)-2)+R*u+L*w+P*u*w
        -Rminus*[u odd]-Lminus*[w odd].

F(1)=1, F(2x+1)=F(x), F(2x)=F(x)+x+(x mod2).
Here the right physical word is the reverse-complement of its address.

## Exact commutation criterion (all orders)

For u odd and w even,

 H(w,u)-H(u,w)=(R-L)(w-u)+Rminus-Lminus.

Thus reflection is always +2 if and only if L=R and Rminus-Lminus=2.
Necessity follows by varying w-u; sufficiency is direct substitution.
For a source edit with profile change delta, its count change is invariant
under the reflection for all such addresses iff delta R=delta L and
 delta Rminus=delta Lminus. This is numerical commutation of the edit
with the reflection. Disjoint graph operations may commute while failing
this stronger count requirement.

Moreover, a reflection-preserving edit is a constant baseline translation
for every odd u/even w iff
 delta h=delta P=delta L=delta R=0.
Its constant is delta kappa-delta Rminus. To prove necessity, compare
u=3 and7, which both have F(u)=1, for arbitrary even w. This forces
 delta P=delta R=0, hence delta L=0. Then varying F(w) forces delta h=0.
The common signed increment need not be zero: the two parity indicators
sum to1. This identifies a precise weaker invariant than freezing all
seven coordinates.

## Today's literal pair fails every reroot test

The current equal-gain construction supplies 131/133 on nine vertices.
At its natural reader (vertex8), the marked profiles are

 (49,1504,1175,92,-47,-10,131),
 (39,1512,1252,94,-38,0,133).

The central +2 becomes +144 at addresses(1,2) and +70 after reflection
at(2,1). Literal complete-tournament counts are 5689/5833 and4643/4713.
This is not a rounding error: the source-change residual is exactly

 2-10(F(u)+F(w)-2)+2u+77w+8uw-10[u odd]-9[w odd].

All nine choices of apex in each of the two tournaments have L!=R.
So simply rerooting either of these actual sources never gives the
free reflected switch. This finite statement does not exclude a larger
host, a different attachment, or a different equal-gain pair.

## Exhaustive small edits of the ACTUAL free core

The inherited deletion masks are (54,32,2,119,6,20,55), apex mask26,
profile (5,256,142,142,-2,0,91). The full center has eight vertices.
Exhaustively reverse each subset of 1,2,or3 of its28 arcs, allowing
both deletion arcs and apex contacts to change, retaining the same
marked apex and the same external flank construction. This is3682
edited tournaments, no isomorphism quotient or random sample.

Exactly1,1,4 edits at distances1,2,3 respectively retain L=R and
Rminus-Lminus=2. Every one has the IDENTICAL full profile as the
starting core. Hence every such edit is invisible to ALL addresses,
not merely to the test words. None supplies a nonzero baseline change.
The replay prints all six edit sets. In particular one-arc recontact
and any coordinated two/three-arc variant within this fixed marked
center do not provide the desired second knob.

The conclusion is deliberately bounded: a fourth-or-larger edit, more
vertices, a changed marked-interface definition, or edits involving
the actual flank vertices remain untested here. The all-order result
is the algebraic commutation/translation criterion, not a general
nonexistence theorem.

## Replay

Run `python3 reflected_return_replay.py`.
The standalone code uses no inherited experiment imports. It independently
computes subset endpoint DPs, all seven coordinates, all3682 edits,
all18 reroots, and eight literal flank tournaments. The symbolic law is
checked against their actual H counts. Output: reflected_return.out.
