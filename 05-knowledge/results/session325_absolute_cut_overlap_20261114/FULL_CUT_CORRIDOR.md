# The ready triple survives every strong binary cut, not just narrow right rails

For a rooted strong left block keep (H,h,L,E), where E is the Hamilton end count at the root and the first three coordinates are as in [Ready overlaps](READY_OVERLAPS.md). With an arbitrary rooted strong right block B, write h_B for its root-deletion count, S_B and E_B for its rooted start/end counts, D_B for its incoming defect, and e_B=h_B−S_B. The [ordered seam excess identity](../session324_busch_inverse_and_gap_structure_20261113/BUSCH_INVERSE.md), collecting left coordinates, is

    H' = h_B H − e_B E + (H(B)+D_B−h_B)L − h(D_B−e_B).

Also h'=h h_B, E'=h E_B, and L'=Q_(B0)(2)(L−h)+h L_B, where L_B is the outgoing relaxation at B's root. Indeed a path ending at the root first follows the two ordered old blocks; the optional-start formula splits its marked start between the two, giving the stated L transfer. The actual start count at the shared root likewise obeys S'=S h_B, since the first nonroot vertex of any spanning path must lie in the left deleted block. Thus the triple also preserves both named endpoint counts under common attachments.

Thus, for left sources with matched h,L, differences transform by ΔH'=h_B ΔH−e_B ΔE. After a common step, E also matches automatically. The order-five pair of Ready overlaps has E=3 in both sources. Its strong order-seven triple has E=15 (and start count 30) in all three. Every common pointed right factor therefore multiplies their gaps by its h_B; gap two persists whenever h_B=1.

## All strong one-apex-over-chain right factors

When h_B=1, B0 is transitive, so B is precisely a one-apex tournament over a chain. Strongness requires its root beat the first chain vertex and lose to the last; both actual root endpoint counts equal one, so e_B=0. Read the incoming cut in reverse chain order as a binary integer b, of bit length n. Then b is positive even; conversely every positive even b is a physical strong right factor. Equivalently it is the outgoing cut of the converse in its increasing transitive order.

Write K(0)=0, K(2b)=K(b), K(2b+1)=K(b)+b+(b mod 2), and rev_n for n-bit reversal. The literal one-apex formula gives H(B)=1+b−K(b), while its incoming relaxation is 1+b. Its outgoing relaxation is 1+(2^n−1)−rev_n(b). Hence the full transfer is

    (H,h,L,E) → (H+L b−h K(b), h, 2^n L−h rev_n(b), h).

In particular every even b>0 yields the three strong, same-order, consecutive odd targets

    135/137/139 + 303 b − 45 K(b),

with one common active interface after construction. A further common right factor of this entire class preserves the same gap.

There is no hidden history parameter here: the ordered cross arcs join successive root-deleted chains into one chain, and the common root sees their concatenated cuts. If the second chosen cut is c (read in reverse order as above), the combined cut has integer 2^n c+b. This also follows algebraically from the binary seam identity; its possible parity correction vanishes because c is even. Different factorizations of one cut are assembly histories of one literal tournament, not distinct covers of the spectrum.

The three coordinates H,h,L are a complete and minimal scalar response for this right grammar. For example the cuts b=2,4,6 give H+2L, H+4L, H+6L−2h, which recover all three. E is not read at all here (it is refreshed to h). In particular scalar H alone cannot substitute for this marked interface even when its difference has been synchronized.

## Absolute selection and one remaining wall

The increments above are multiples of six, and the three starting counts represent the three odd classes modulo six. Thus the target itself determines which of the three initial tracks to use. Write b=2z. Exact attainment in this corridor asks for a positive integer z solving

    (N−H_base)/6 = 101z−15·K(z)/2.

In particular z is fixed modulo 15 by this equation, and K(z)≤z bounds z between N-shifted positive linear endpoints. These are exact inverse filters, not interchangeable with realization. The [two-sided square-root inverse](../session317_exact_anchor_atlas_20261113/EXACT_INVERSE.md) applies directly to the binary scalar formula, keeping the parity guard from the split.

Even this enlarged corridor remains absolutely incomplete: for z≥1 the increment is at least (606−45)z=561z, so its three tracks give at most X/187+O(1) distinct values through X. This is stronger than the density bound obtained by counting all positive b. The ready triple solves the continuation problem locally and preserves all three consecutive labels at every legal cut; it does not make those charts overlap at arbitrary absolute targets. The missing ingredient is a target-driven change of the active (h,L) state or an independent class with much larger capacity.
