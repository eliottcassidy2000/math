# Two independent arc choices: a persistent +2 and a stored doubling digit

2026-11-28. The old +2 and ×2 sessions do not all double the same quantity. This note joins several of their mechanisms at one named root and builds four literal tournaments with two **independent** arc switches. One switch always adds 2; the other can be invisible in the Hamiltonian count, stored, doubled for arbitrarily many steps, and then exposed.

## 1. Put the old operations in one set of coordinates

Use the [moving-root registers](../chanpak_moving_root_and_module_interfaces_20261128) \((h,e,u,g)\): total Hamiltonian paths, paths ending at the root, ordered root-starting split sum, and the root-deletion count. The old root stays fixed under P (append a sink) and K (append a last vertex beating only the root). M uses the same arcs as K but moves the root to the newcomer. All three actions are homogeneous linear:
\[
P(h,e,u,g)=(h,0,2u,g),\quad K(h,e,u,g)=(h-e+u,g,2u-g,g),\quad
M(h,e,u,g)=(h-e+u,h-e,h+u,h).
\]

The [restarted rows](../ready_restarts_rounding_20261116) track the signed path sums \(r,s\) from the root, split by even versus odd number of vertices in the rooted path, and prove that K changes h by \(2r\). Since \(u=r+s\) by the definitions and the same change is \(u-e\), the exact translation is
\[
r=(u-e)/2,\qquad s=(u+e)/2,\qquad s-r=e.
\]
In particular the earlier renewal guard \(s-r=g\) is the K chart \(e=g\). There, with \(x=u-g,\ a=h-u+g\), the height after \(k\) K appends is
\[
h_k=a+2^k x .
\]
For two hosts on that chart under the same continuations, the marked height gap is \(\Delta a+2^k\Delta x\). Thus the gap is constant iff \(\Delta x=0\), or doubles at each step iff \(\Delta a=0\). An offset-doubling law for one host, a literal doubling of a gap across two hosts, and doubling the *number of vertices* remain distinct, as in [Session 222](../session222_tiling_recursions_20261107). Also, the two directions of the [path-splice array](../path_splice_rows_continuation_20261116) \(h+2^kB\) become two physical response channels below.

A pair keeps the same height gap under **every** common P/K continuation iff its difference has the form \((\Delta h,0,0,0)\); call this a flat pair. Sufficiency is immediate from the linear actions. Necessity follows by asking the chronological suffixes K, PK and KK: preserving the gap under K gives \(\Delta u=\Delta e\), PK gives \(\Delta u=0\), and KK then gives \(\Delta g=0\). This is the exact common-port condition for this interface, analogous in role to the [renewable-bank checks](../session282_chanpak_20261112), not a claim about all contexts.

There is a useful obstruction: starting from a flat pair, a common P/K/M word can never turn a +2 gap into a +4 gap. Before the first M the difference is fixed. That M sends \((d,0,0,0)\) to \(d(1,1,1,1)\). For any subsequent word W the gap is exactly \(dH(W)\), with \(H(W)\) read from the singleton seed. This is a physical, same-old-arc extension of the readout \(d\mapsto dH(W)\), akin to the [external one-arc module reader](../chanpak_mixed_defects_20261115/ARC_PORT.md), but uses a single shared root chronology. The earlier [two-block bridge law](../chanpak_moving_root_and_module_interfaces_20261128#3-a-two-block-construction-and-a-finite-negative-test) extends the same response to an arbitrary tail tournament B: make all old sites beat B except for one arc from B's marked root to the old root. In its formula the e and u terms cancel between the flat hosts, leaving exactly dH(B). Thus even a one-chord edit wholly inside the tail changes that coefficient without touching the original marked arc. Since [Rédei parity](https://arxiv.org/html/2307.05569) makes \(H(W)\) odd, a *separate* digit switch is needed for pure doubling.

## 2. Cancel the high switch once, then expose it at any dyadic scale

**General gate.** Start from a flat pair with common \(e,u,g\) satisfying \(u-e=2^r g\) with \(r\ge1\). In each host, append one new vertex b, choosing either P or K on the root–b arc. This is one literal switch, not a change of root. Taking K minus P its differential state is
\[
(2^r g,\ g,\ -g,\ 0).
\]
Now append the same \(P^r K\) to both choices. The P's erase the end-count difference and double the split difference r times. The final K cancels the whole height difference, leaving
\[
(0,\ 0,\ -2^{r+1}g,\ 0).
\]
The original flat switch remains exactly \((d,0,0,0)\) throughout. The two responses are independent because all common maps are linear and the first P/K choice changes the same linear read of e,u,g for both original hosts. Equivalently the four register vectors have zero mixed finite difference.

Append a further common \(P^k\). The second switch is still invisible in total H, and its stored difference has doubled k times. One final K exposes the gap \(2^{r+k+1}g\) while leaving the first gap d fixed. This is not the impossible operation \(H\mapsto2H\) on an odd Hamiltonian count; it really doubles a stored and subsequently observable difference.

**Literal two-arc family.** Take the word PKKKMP on the singleton, root 5 (vertices labelled by birth beginning with 0). The outgoing masks are
\[
(66,124,121,113,97,65,0).
\]
Reverse the off-tree pair 1–3 for the lower member. This is the [one-chord carrier](../chanpak_trijection_one_chord_20261128) interface applied to a unit gap rather than a single height:
\[
(15,0,36,9),\qquad (17,0,36,9).
\]
Here \(36=4\cdot9\), so use \(r=2\). The second edit is the pair 5–7 (old root and first newcomer), disjoint from 1–3. After the gate followed by PPK, the four literal states are

| First switch | Gate letter | h | e | u | g |
|---|---|---:|---:|---:|---:|
| lower | K | 303 | 9 | 495 | 9 |
| upper | K | 305 | 9 | 495 | 9 |
| lower | P | 303 | 9 | 567 | 9 |
| upper | P | 305 | 9 | 567 | 9 |

For every \(k\ge1\), the common \(P^kK\) gives four tournaments on \(12+k\) vertices:
\[
B_k,\quad B_k+2,\quad B_k+72\cdot2^k,\quad B_k+72\cdot2^k+2,
\qquad B_k=303+495\cdot2^k .
\]
Any switch changes **only** its original arc. In particular the first rectangle is \(1293,1295,1437,1439\). The verifier directly counts both seed hosts, all stored states, and all eight revealed hosts for k=1,2. The formula, not the tests, establishes the infinite family.

This supplies two independent edits but not two freely adjustable absolute coordinates. In particular the [free reflected seam](../chanpak_wheel_square_selector_20261117/FREE_SEAM.md) still offers a much wider placement mechanism, and neither method gives arbitrary integer targeting.

**Why this gate does not automatically rearm.** In the literal example \(g=9\) is fixed under P/K and each branch after storage has integer \(u/g\), odd and greater than 3, with \(e=g\). If P appends occur, \(u/g\) doubles and e becomes 0, so \(u/g\) keeps an odd part greater than 1. On a K chart, successive K moves just double \((u-g)/g\). If a run \(P^tK\) with \(t\ge1\) follows a K chart, the new value of \((u-g)/g\) is \(2^{t+1}(u_{\rm old}/g)-2\), with odd part \(2^t(u_{\rm old}/g)-1>1\). Hence if the initial \((u-g)/g\) already has odd part greater than 1, the guard \((u-e)/g=2^r\) never returns by P/K alone. The stored example has these odd parts 27 and 31. A restart must move or redesign the root, rather than just make the same fixed-root word longer.

## 3. An entire marked suffix block still closes on the four registers

Here is a separate reusable transfer that avoids guessing which internal data an inserted gadget needs. Glue a nonempty suffix B to a rooted host T, making every old nonroot vertex beat every site in B. The old root a has arbitrary contacts to B, specified by a pointed tournament \(C=B\cup\{a\}\). The root of the enlarged tournament stays at a. Write \(b=H(B)\), \(j=H(C)\), \(s_C,e_C\) for the numbers of Hamiltonian paths of C that start or end at a. Let \(u_C\) be the outgoing split register of C, and \(v_C=u(C^{\rm op},a)\) its reversed split register. Put
\[
Q(B)=\sum_{J\subseteq B}H(B[J]) H(B\setminus J).
\]
Then the exact transformed registers are
\[
\begin{aligned}
h'&=(h-e)b+e s_C+(u-g)(v_C-b)+g(j-s_C),\\
e'&=g e_C,\qquad
u'=(u-g)Q(B)+g u_C,\qquad
g'=gb .
\end{aligned}
\]

For h', split by whether any site of B occurs before a. If not, all of T is walked first, and only paths ending at a need the suffix to start along an a-contact. If some of B occurs before a, split the old nonroot set into its prefix and suffix around a. A nonempty old suffix supplies the register \(u-g\), after which the remaining B suffix can start freely. A singleton a-suffix in the old set supplies g, and the B pieces must form a C-path with some site before a, counted by \(j-s_C\). The endpoint and deletion formulas follow from old nonroots beating B. For u', if the root-starting part contains an old nonroot it can append its B-sites freely and sum over both B parts, giving \(Q(B)\); otherwise the gadget contributes \(u_C\). This is a marked-block law, not the unmarked substitution law.

In particular every common block with transitive B leaves a flat differential unchanged. This gives a precise enlargement of the moves available when looking for a gate restart, while the coefficient table specifies exactly which additional boundary data to compute. The verifier replays the whole map on 180 independently generated rooted host/block pairs.

## Replay

Run python3 verify.py. No search output is needed for the witness: the literal seed, outgoing masks, edit labels, and arithmetic are given above.
