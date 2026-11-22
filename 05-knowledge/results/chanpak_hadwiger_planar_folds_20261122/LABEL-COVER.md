# Relabelling holonomy: the flip quotient's two-sheet cover

**PROVED.** For order at least two, fix the orientation bit for each unordered pair of labels in their numerical order. Let \(\epsilon(T)\in\mathbb F_2\) be the parity of the sum of all bits. A single physical flip toggles \(\epsilon\). Relabelling by \(\pi\) permutes bit positions and complements precisely the coordinates for label-order inversions. Consequently
\[
          \epsilon(\pi T)-\epsilon(T)
                =\operatorname{inv}(\pi)\pmod2
\]
for every tournament, with no dependence on its orientation. In particular every automorphism is even (alternatively its group has odd size), so each ordinary tournament isomorphism class has precisely two classes under even relabellings, distinguished by \(\epsilon\).

Pass from the labelled flip cube to quotient graphs by the alternating group and by the full symmetric group. On either graph keep one edge when some physical arc flip joins its endpoints; on the full quotient this may include a loop if the flip preserves the isomorphism class. The alternating quotient is exactly the two-sheet bipartite cover of the full quotient: every edge lifts between opposite parities in both ways, obtained by applying an odd relabelling to a literal lift, and there is only one vertex per full isomorphism class on each parity sheet. This description concerns the simple adjacency support, not edge multiplicities. If loops are omitted in the base, omit their matching lift edges.

Thus an odd cyclic walk in the ordinary flip quotient lifts to a path that returns to the old ordinary class on the *other* sheet. Equivalently any relabelling that closes its lifted orientation is odd. This is the precise distinction between the labelled bipartite cube and a geometric odd cycle in the unlabelled quotient. It specializes to the marked Moser spindle comparison in this packet.

There is a second, independent parity. Rédei makes all ordinary heights odd, and along any closed quotient walk the signed increments telescope, so
\[
   \#\{\text{steps with } \Delta H=2\bmod4\}\equiv0\bmod2.
\]
In particular every odd closed walk contains an **odd** number of steps with \(\Delta H=0\bmod4\). A graph consisting solely of literal \(\pm2\) edges stays bipartite even after full relabelling, now colored by \(H\bmod4\). Geometric unit length on a quotient edge is therefore not a reliable label for unit height. The two parity observables play different roles at the boundary, much like the distinct sheet data in the [square-sum mirror](../chanpak_middle_out_square_sum_20261122/SQUARE-SUM.md).
