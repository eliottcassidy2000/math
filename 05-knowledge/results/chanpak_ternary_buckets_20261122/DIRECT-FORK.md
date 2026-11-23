# A physical tripling edge with two single-arc exits

**PROVED.** For every integer \(s\ge1\), on the same \(s+8\) labelled vertices there are tournaments of heights
\[
                      h,\quad 3h,\quad 3h-2,\quad 3h+2,
                 \qquad h=21\cdot2^s-13 .
\]
Each of the three outside tournaments differs by one specified arc from the tournament of height \(3h\). The first forks have sources 29, 71 and 155. The arc from 29 takes it to 87, whose two exits lead to 85 and 89. There is an arbitrary-address extension of the positive exit and a different arbitrary-address family for the negative exit.

## Four seven-vertex readings

Let the ordinary vertices be labelled 0 through 6. The rows of \(P\) are \((0,1,11,19,7,31,63)\), with bit \(j\) in row \(i\) denoting \(i\to j\). The apex beats the ordinary vertices in mask 11. Obtain \(C\) by reversing \((4,6)\) in \(P\), and obtain \(L,U\) by reversing \((2,3),(4,5)\) respectively in \(C\). Each of these three choices is independent of flank length and address.

The [two-flank formula](../chanpak_gap_star_ports_20261122/README.md) gives the height of each tagged core at binary addresses \(x,y\) as
\[
 h(f(x)+f(y))+qxy+tx+cy+r[x\text{ odd}]+e[y\text{ odd}]+k.
\]
Here \(f(0)=0,\ f(2v+1)=f(v),\ f(2v)=f(v)+v+[v\text{ odd}]\). The entries below are finite counts of paths and ordered two-path covers. The [replay](direct_fork.py) reconstructs all four profiles, builds the actual tournament rows at sample addresses, and counts paths by subsets.

|Core|\(h\)|\(q\)|\(c\)|\(t\)|\(e\)|\(r\)|\(k\)|
|---|---:|---:|---:|---:|---:|---:|---:|
|\(P\)|3|192|21|33|3|3|5|
|\(C\)|9|336|63|54|9|8|15|
|\(L\)|5|256|63|50|5|4|17|
|\(U\)|9|336|63|60|9|8|17|

Write \(P,C,L,U\) also for their full heights. Subtracting rows gives
\[
 C-3P=-240xy-45x-[x\text{ odd}],\qquad U-C=2+6x .
\]
Hence exact tripling holds if and only if \(x=0\). On this entire face the two successive arcs give
\[
       P=3 f(y)+21y+3[y\text{ odd}]+5\quad\longrightarrow\quad 3P
       \quad\longrightarrow\quad 3P+2 .
\]
For the other exit, \(L-C=2-4[y\text{ odd}]-4f(y)\). The recurrence makes \(f\) even everywhere, and \(f(y)=0\) precisely at zero and at positive integers whose binary digits are all ones. Indeed an even positive argument has a positive increment, while an odd argument loses its final digit. Thus \(L=C-2\) holds precisely for \(y=2^s-1,\ s\ge1\). This proves the common-host fork and its exact address domain. The [finite-core classification](DIRECT-CENSUS.md) allows other fixed arcs and arbitrary odd prefixes in both moving-suffix modes.

## Substitution and the obstruction at a perfect number

**PROVED.** Attach any nonempty tagged ordinary tournament \(B\) to the right of the seven-vertex core by directing all cross arcs toward \(B\); retain the same apex and any further right address \(y\). Write \(H_B\) for the ordinary height of \(B\). Then for every \(B\), its tags and \(y\),
\[
 C_B=3P_B,\qquad U_B=C_B+2H_B,\qquad P_B\equiv2H_B\pmod3 .
\]
In fact the right-face profile of a join is \((hH,Qc+hC,hE,Hk+h(K-H)+cT+eR)\) in coordinates \((h,c,e,k)\), where capitals denote the right block. The four coordinates are linear in the left data \((h,c,e,k)\). Those data triple from \(P\) to \(C\), and from \(C\) to \(U\) only the final coordinate changes, by two. In \(P\) every coordinate except \(k\) is divisible by three. This proves all three claims. [Shared-apex ordinal law](../chanpak_e31_physical_20261122/README.md)

The designated upper exit stays at distance two exactly when \(B\) has ordinary height one, equivalently when it is transitive. For the equivalence, any subtournament Hamiltonian path extends by insertion of each missing vertex, retaining its old order. A nontransitive tournament contains a directed triangle with three path orders.

For odd \(p\ge3\), put \(N_p=2^{p-1}(2^p-1)\). Its residue modulo nine is one: \(A=2^{p-1}\) equals 1, 4 or 7 modulo nine and \(N_p=2A^2-A\). Thus a center of the form \(3P\) on the two sides of \(N_p\) must be \(N_p-1\), and its parent must be divisible by three. On the naked right face \(P\equiv2\pmod3\). Under arbitrary right substitution, a compatible center requires \(3\mid H_B\), which makes the designated outgoing step divisible by six. This applies to exponent 61 even if the search allows unbounded binary right words.

**PROVED.** An arbitrary tagged *left* ordinal prefix \(A\) cannot compensate this defect. In the join \(A\to C\to B\), relative to three times the parallel parent construction, the difference at any right address is
\[
   -240c_A(Q_B y+T_B)-H_B(45c_A+e_A).
\]
Indeed after adjoining only \(A\), the difference in right-face coordinates \((h,c,e,k)\) is \((0,-240c_A,0,-45c_A-e_A)\); now apply the join law to \(B\). An empty \(B\) uses \(Q_B=H_B=1,T_B=0\). All quantities in the brackets are nonnegative. Exact tripling at even a single address thus forces \(c_A=0\), equivalent to an empty tag set since every vertex starts a nonempty two-slot cover with its singleton first. Conversely such a prefix points toward the apex and every later vertex. It multiplies both heights and the outgoing step by \(H_A\). The same argument permits an arbitrary further left binary flank by including it inside \(A\).

**PROVED with a finite seven-vertex check.** If \(H_B>1\), *no* fixed reversal within the ordinary join \(C\to B\) or from its apex to the join gives an exit of \(+2\) at every right address \(y\) when \(3\mid H_B\). There are three cases.

1. A change inside \(B\) or in its apex tags preserves divisibility of the entire central height by three, by the substitution formula for \(C\).
2. A change inside \(C\) or in its apex tags must preserve the joined ordinary height and right slope. Indeed the four functions \(f(y),y,[y\text{ odd}],1\) are linearly independent, by evaluation at \(0,1,3,2\). Since \(H_B,Q_B>0\), the join formula forces preservation of \(h_C,c_C\). Among the 28 such single reversals, these two coordinates survive only at \((0,1)\), preserving the entire right-face profile, and \((4,5)\), giving its designated increment \(2H_B\). This is the seven-vertex check in the replay.
3. For a cross arc \(u\to v\), \(u\in C,v\in B\), ordinary height survives only if \(u\) is a sink in \(C\) and \(v\) a source in \(B\). To see this, losses equal the number of full paths in \(C\) ending at \(u\) times the number of full paths in \(B\) starting at \(v\). Gains use \(v\to u\). They are counted by splitting both blocks into their ordered before/after paths, with the second \(C\) path starting at \(u\) and the first \(B\) path ending at \(v\). These choices are independent. The first factor is at least \(H(C-u)\), by putting \(u\) alone in the second slot, and \(H(C-u)\) bounds the first loss factor from above. If \(u\to w\) inside \(C\), a two-vertex second slot \((u,w)\) makes the first inequality strict. The dual argument holds in \(B\), and both gains factors are positive. For a sink and source the reversal exchanges adjacent ordinary vertices with identical external incidences. Identical tags therefore give an isomorphism. Different tags swap the two tag weights in the original join. The contribution to right slope \(c\) of a tag at the sink \(u\) is \(H_C Q_B\), where \(Q_B\) counts ordered two-slot covers. A tag at the source \(v\) contributes \(H_C Q_{B-v}=H_C Q_B/2\). The latter identity assigns the source to the beginning of either slot. For the former, starting at the sink uses that single \(C\) vertex in the first slot; all the other \(C\) vertices occupy the second slot, and the two \(B\) slots are free. Thus a swap of distinct tags changes the slope, ruling out a persistent offset.

This obstruction allows a correction at a *specified* right address, and leaves reversals touching the exterior flank eligible.

## An unrestricted lower exit in a different chart

**PROVED.** Take center rows \((0,5,17,7,75,31,47)\) with mask 5. The parent differs at \((4,6)\), and the lower child at \((3,4)\). Their right-face profiles \((h,c,e,k)\) are respectively \((5,15,5,9),(15,45,15,27),(15,45,15,25)\) in parent, center, child order. Therefore for every \(y\ge0\) two physical edges realize \(P\to3P\to3P-2\) with \(P=5f(y)+15y+5[y\text{ odd}]+9\). Across this chart and the first one, both ternary signs have unrestricted binary address families; requiring both signs from the *same* center gives the stricter domain above.

Run from this directory with `python3 -B direct_fork.py`.
