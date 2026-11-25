# The complete unit-gap spectrum of flips inside a one-apex spine

**PROVED.** Start with any transitive tournament, add an apex with arbitrary incident orientations, and reverse one arc of the original transitive spine. If the two Hamilton path counts differ by two, the possible midpoints form the set
\[
 \boxed{\{2,4,10,18\}\ \cup \{6\cdot2^j,\;9\cdot2^j-2:j\ge1\}.}
\]
Every displayed midpoint occurs. A unit flip crosses at most **two interior vertices** of the original spine. With two interior vertices the only possible midpoints are 10 and 18. With one interior vertex the two infinite families occur and both hosts can be strong. In particular the family \(9\cdot2^j-2\) supplies physical midpoints equal to 1 modulo 3. The old transitive spine may become cyclic after the edit.

Consequently **no** such reversal crosses \(E_p=2^{p-1}(2^p-1)\) at any exponent \(p\ge3\). The odd part of \(6\cdot2^j\) is three. For \(j\ge2\), \(9\cdot2^j-2\) is two modulo four, whereas \(E_p\) is zero modulo four. Its \(j=1\) value and the four isolated values are all below \(E_3=28\). In fact **each endpoint of any exact straddle must need at least two vertex deletions to become transitive**. If an endpoint has a transitive one-vertex deletion, every incident flip is either a spine flip covered here or an apex-contact flip. For an apex flip strip any common ordinal sources or sinks. Two positive flank addresses give two strong hosts, already excluded by the cited apex-contact theorem. If one flank address is zero, reflecting and taking the converse exchanges the two sides, so take the left address zero. The [one-contact identity](../chanpak_apex_flip_local_global_20261122/README.md) then gives gap \(2\lceil w/2\rceil\) and original count \(1+f(w)\). A unit flip forces \(w=1,2\), with midpoint 2 or 4. This covers nonstrong endpoints as well. Earlier finite sampling of spine flips found small changes infrequently. [Apex contact and transitivity-preserving flips](../chanpak_support_birth_grammar_20261122/ODD-PERFECT-PORT.md) · [Earlier numerical spine-flip survey](../chanpak_waves_20261113/README.md)

The scalar triple \(9\cdot2^j-5,9\cdot2^j-3,9\cdot2^j-1\) was previously attained with several hosts. The same midpoint sequence appears as a specialization of another module formula. Here one marked arc of a one-apex host gives the last two counts. [Scalar arm triple](../chanpak_chronological_core_changes_20261114/ONE_OUTLET.md) · [Previous module family II](../chanpak_absolute_liouville_leech_20261122/MODULE-OFFSETS.md)

## Why only one endpoint can have a deficit

Write the core from the first edited endpoint to the second as \(D=\mathrm{TT}_{m+2}\), with \(m\ge1\) interior vertices. Mark each site that the apex beats. Denote the edited core by \(D'\). The exterior contacts are the integer \(x\) on the left and the reverse-complement reading \(y\) on the right. Removing leading zero contacts or terminal one contacts deletes only common ordinal sources or sinks and preserves both counts.

Use the two-arm formula with \(h=H(D)\), \(Q\) the number of ordered two-slot path covers, \(k=H(D+\text{apex})\). The second slot can be empty. Write \(c,o\) for the total and odd-first-slot subcount of ordered covers whose nonempty first slot starts at a marked vertex. Write \(t,o'\) for the analogous end-unmarked reads. Denote a change from \(D\) to \(D'\) by \(\delta\). With \(f(0)=0, f(2z+1)=f(z), f(2z)=f(z)+z+[z\text{ odd}]\), the signed change in the full count is
\[
\begin{aligned}
\Delta={}&\delta k+\delta h(f(x)+f(y))+\delta Qxy\\
 &+\delta t(x-[x\text{ odd}])+2\delta o'[x\text{ odd}]\\
 &+\delta c(y-[y\text{ odd}])+2\delta o[y\text{ odd}].
\end{aligned} \tag{1}
\]
Every \(f(z)\) is a nonnegative even integer. [Two-arm count and definitions](../chanpak_gap_star_ports_20261122/README.md)

Let \(h_*=2^m\). Then \(\delta h=h_*\), \(\delta Q=2(3^m-1)\), and \(\delta k\ge h_*-2\). In the edited spine each odd cycle uses the reversed arc and an odd-sized subset of its interior, in increasing order. Such cycles meet one another. The odd-cycle packing formula gives the first identity. In counting the change in two-slot covers, both endpoints must occupy the same slot, which can be chosen in two ways. If \(s\ge1\) interior sites accompany them, that slot gains \(2^s\) paths, while its complement stays transitive. Summing over all subsets gives the second identity. In the pointed core every previous odd cycle passes through the apex. The only previous cycle that can use the edited arc is the triangle apex–first–last. All other previous cycles survive. Counting the new spine cycles alone, and allowing for that single lost triangle, proves the lower bound. [Odd-cycle packing, THM-002](../../../01-canon/theorems/THM-002-ocf.md)

Consider the individual rooted-start contribution from a vertex \(v\) to the cover count. For the first core vertex, the change is precisely **minus one even first slot and zero odd first slots**. Starting there precludes later traversal of the reversed arc. The only lost path has first slot equal to the pair of endpoints. For every other start vertex both parity counts are nondecreasing. If both endpoints are in its slot, the old transitive path could never start at \(v\). If both endpoints belong to the other slot, that slot still has at least one Hamilton path. If the endpoints split, both slots are unchanged. Moreover every other start has at least one strictly positive gain: the last vertex can now start on the endpoint pair, and an interior vertex can now start on its triangle with the two endpoints. No restriction on the rest of the tags enters this argument.

It follows that \(\delta c\ge-1\) and \(\delta o\ge0\), with \(\delta c<0\) precisely when **only the first vertex is tagged**. Reflection and reversing all orientations gives the same statement for the end reads, with last and untagged in place of first and tagged. For \(x,y\ge1\), (1) gives at least
\[
 4xy-(x-[x\text{ odd}])-(y-[y\text{ odd}])\ge4.
\]
Indeed \(x-[x\text{ odd}]\le2(x-1)\) for every positive \(x\), and similarly for \(y\). Thus every unit gap lies on one of the two axes. Reflection exchanges them, so take \(x=0\).

## Solving the exceptional mask at every span

If only the first vertex is tagged, it belongs to every possible apex cycle. No new such cycle can use the reversed arc, since it would revisit that first vertex. The lost triangle exists, and no spine cycle is disjoint from an apex cycle. Thus \(\delta k=h_*-2\) exactly. In the old core all \(2^{m+1}\) subsets containing the first vertex give one starting cover, half of each parity. Its standalone apex count is \(2^{m+1}+1\), by the transitive one-apex formula. Put \(F=1+f(y)\), \(\eta=[y\text{ odd}]\). Then both the signed gap and the original count become
\[
 \Delta=h_*F-y+\eta-2,\qquad H_{\rm old}=2h_*(1+y)+F. \tag{2}
\]
In particular \(F\) is positive and odd. Also \(f(w)=0\) precisely when the canonical binary word is all ones or the address is zero, by its recursion.

* For \(\Delta=-2\), we need \(y=h_*F+\eta\). If \(m\ge2\), appending the \(m\) low bits contributes at least \((h_*/2-1)F\ge F\) to \(f(y)\), impossible. This follows by iterating the even recursion on the first \(m-1\) appended zeros, then taking the cheaper last bit one. If \(m=1\), the even case also gives \(f(y)\ge F\). The odd case gives \(F=1+f(F)\le1+(F-1)/2\), so \(F=1,y=3\). Its midpoint is 16.
* For \(\Delta=2\), we need \(y=h_*F-4+\eta\). If \(m\ge2\), this appends the \(m\)-bit suffix \(1^{m-2}0\eta\) to \(F-1\). As \(F-1\) is even, the coefficient of \(F-1\) in the appended read is \(3h_*/4\) for \(\eta=0\) and \(h_*/4\) for \(\eta=1\). The remaining terms \(f(F-1)\) and the suffix read are nonnegative. These coefficients follow from the zero-bit append recurrence, or equivalently from summing the prefix weights at the zero positions. For \(m\ge3\) this forces \(F\ge1+2(F-1)\), hence \(F=1\). The remaining suffix contains 10, so its read is positive, a contradiction. For \(m=2\), the even case likewise forces \(F=1\). The odd case reduces to \(f(F-1)=0\), again forcing \(F=1\) since \(F-1\) is even. The resulting \(y=0,1\) give midpoints 10 and 18.
* For \(\Delta=2,m=1\), put \(z=F-2\), which must be positive and odd. The odd choice \(y=2z+1\) yields \(F=1+f(z)\le1+z<F\). For \(y=2z\), the even read yields \(f(z)=0\), so \(z=2^j-1\), \(j\ge1\). Formula (2) gives exactly the midpoint \(9\cdot2^j-2\). Its original contact word is \(10010^j\); reverse the arc between the first and third spine sites.

## The remaining masks and sharpness

Here \(\delta c,\delta o\ge0\). Equation (1) gives \(\Delta\ge2^m-2\), leaving only \(m=1,2\) for a unit edge. For \(\delta c>0\), it also gives \(y-[y\text{ odd}]\le2\), hence \(y\le3\). The case \(\delta c=0\) has just three possibilities. An empty mask at either span gives \(\Delta=2^m(1+f(y))\) and \(H_{\rm old}=1+f(y)\), giving only midpoint 2. Otherwise the mask is 110 at \(m=1\). Its exact data are \(\delta k=0,\ \delta o=1,\ H_{\rm old}=5+6y+f(y)\), so \(\Delta=2(f(y)+[y\text{ odd}])\). Its unit solutions are precisely \(y=2^j-1, j\ge1\). They give the midpoints \(6\cdot2^j\) by flipping the first-to-third arc of contact word \(1100^j\).

The short table below enumerates *all* masks and \(0\le y\le3\), including exceptional masks for clarity. The last column lists every attained midpoint. These finite values also show which isolated midpoints really occur.

| Interior vertices | Low-word and mask cases | Unit midpoints |
|---:|---:|---|
| 1 | \(4\cdot 2^3\) | 2, 4, 12, 16, 24 |
| 2 | \(4\cdot 2^4\) | 10, 18 |

When \(m=0\) the reversed arc swaps adjacent vertices in the underlying transitive order. The inherited adjacent-swap rule requires both active flank addresses to be 0 or 1. The four choices give midpoints 2, 4, 4 and 10. This completes the classification. [Earlier adjacent-swap rule](../chanpak_exponent31_seams_20261122/README.md)

For both infinite families the original word starts at a tag and ends at a nontag, so it is strong. After the edit the initial triangle is strong and the apex has contacts in both directions there. Every later vertex is reached from the triangle, and it reaches the apex directly or via the final nontag. This verifies strong connectivity of the actual single-arc pairs.

Run `python3 -B spine_arc.py` in this directory for the small mask table, every one-apex word and every spine flip through spine length seven, and larger literal examples from both infinite families.
