# A literal handoff: chirality, module relocation, and the primitive return

**PROVED IN NOTES:** a common named intermediate for the two constructions,
an exact neutral module relocation, and their composite transition.
**FINITE-EXACT:** tiny complete quotient coefficients and literal controls.
This closes the handoff for the stated arbitrary-module family; it is not
a selector for arbitrary prime or7p targets.

Let M be any nonempty tournament, h=H(M),p=Q_M(2). The two previous
notes independently produced the scalar9p−h. Scalar equality would not
justify composition. There is, however, an exact labelled connection.

## 1. The common graph, with different module roles

Start with [CHIRALITY](CHIRALITY.md)'s positive graph
C3[TT2,M,1], with names r,a in TT2, root r, and b in the final singleton.
Apply N then R: add x with r→x→all other old vertices, select x, then add
y with x→y→all other current vertices and retain x. The resulting count
is9p−h. Retain the named pair(r,x) after the first birth and keep that
pair after the second, discarding y from the port set. This is not the
usual sealed balanced pair(x,y): the original r is explicitly retained
through both births and is never recovered from a discarded name.

The [selection repair](selection/PROOF.md) has intermediate skeleton

    I=(46,8,58,48,35,2), root4.

There is a literal graph equality under the roles

    (r,a,M,b,x,y) ↦ (4,5,1,3,0,2).                  (1)

Thus the NR graph is I with M in slot1, not slot0. Equation(1) is an
equality of all block incidences for arbitrary M, retaining every internal
arc of M. The scalar agreement has an actual source witness.

At root4, the two placements have states

    I[M at1]: (9p−h,p+h,p+h,11p−h),
    I[M at0]: (9p−h,3h,3h,9p+3h).                 (2)

The second line is inherited from the new selection lane. For the first,
M has exactly one outgoing exterior neighbor, role3. In a Hamilton path,
only one M-run can leave the module and at most one can be terminal.
For F the first path ends at root4 outside M, so at most the second path
can have a terminal M-run. The same two-run bound holds after deleting
root4 and for paths starting there. Hence every observation has at most
two module runs. Complete quotient-word enumeration gives

| observation | I with M at1: coefficients of h,c2 | I with M at0 |
|---|---|---|
| H | (17,9) | (17,9) |
| s | (3,1) | (3,0) |
| q | (3,1) | (3,0) |
| F | (21,11) | (21,9) |

Here c2=p−2h. This proves (2) at all orders; the finite coefficient
calculation is complete because the run bound is proved first.

## 2. Move the module and retain the original root star

In I[M at1], let u be the singleton in role0. Keep the same vertices of M,
the same u, and all internal M arcs. Exchange their exterior roles, so
M occupies role0 and u occupies role1. This is the following actual edit:

- reverse all edges between M and each of u, role2, and role5;
- reverse the edges u–2 and u–5.

There are exactly3|M|+2 reversed arcs. No arc at root4 changes: root4
beats both roles0 and1, before and after. The complete root star survives,
and H is unchanged. The response change is exactly

    (ΔH,Δs,Δq,ΔF)=(0,−c2,−c2,−2c2).              (3)

This is a neutral marked preparation, not a full rooted congruence. It
uses a genuinely changing module boundary, with a precise response cost.
For |M|≥2, c2>0. The root-deletion count then changes, so the move
necessarily leaves the labelled cyclic-triangle fibre even though H and
the root star are retained. Its next R-birth target decreases by exactly
c2: from19p−3h to18p−h. The singleton is the equality boundary.
Both endpoints are strong by their expanded skeleton cycles. We do not
require each individual edge reversal in the batch to preserve H or
strongness; (3) concerns the completed simultaneous edit.

Now role1 is again the singleton u. Add a new helper z which beats only
u and loses to every other vertex. The graph is precisely selection's
T[M]. Its rooted state becomes

    (15p−h,5h,5h,15p+5h).                          (4)

The resulting actual counts therefore follow

    3p−h → 9p−h → 9p−h → 15p−h.                  (5)

The first arrow is NR, the middle arrow is the root-star-preserving
module relocation, and the last arrow is the named sole-exit helper.
If X=3p−h, the first and final counts are3X+2h and5X+4h. These are
actual transitions with the stated module and port data, not a free
arithmetic action on every odd X.

## 3. What the composition retains and what it does not

All old vertices survive this handoff. The internal graph M and the
original named root r survive. The relocation retains the whole star of r,
but changes other old arcs and the rooted response by (3). The two births
and final helper have their explicitly specified new incidences. The
old active NR vertex x changes exterior role; its prior response is not
claimed to survive.

At the endpoint, the selection lane's neutral S[M]↔T[M] repair and its
commuting R digit apply. A common fresh seal gives the stated complete
future two-port congruence. It remains necessary to supply the module
responses(h,p) for a prescribed scalar target. The new result closes a
physical inter-frame handoff, not this joint arithmetic realizability gate.
It is not an automatic return to the initial C3 frame. The relocation
uses the retained full frame labels, including roles2 and5; it is an
explicit operation beyond the six-command two-port language. Only after
choosing the endpoint and its fresh seal do we invoke that language’s
future-response congruence.

For M=TT_m, m≥1, examples of (5) are

    5→17→17→29, 11→35→35→59, 23→71→71→119.

For the H11,P2=100 module D in the selection lane, the counts are
289→889→889→1489. The same mechanism works with arbitrary nontransitive
M and with the unbounded ordinal module families already constructed.

Run `python3 -B 05-knowledge/results/zeckendorf_target_selection_20261113/handoff.py`.
The checker preserves the exact names in (1), every changed edge in (3),
the two complete quotient tables, and all four states. It tests all
labelled modules of orders1..3 plus two order4 controls. Larger counts
follow the proved response identity, not a large tournament census.
