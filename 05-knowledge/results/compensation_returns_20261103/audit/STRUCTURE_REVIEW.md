# Independent logical review of deletion-cap structure

2026-11-03. **VERIFIED:** the repaired all-order entry-height theorem,
its thirteen-name consequence, and the transported-auxiliary endpoint law
in [`structure/PROOF.md`](../structure/PROOF.md). This is a logical review,
not a new graph census. The existing finite carrier audits remain useful
cross-checks; the final obstruction does not depend on extrapolating them.

The final reviewed proof has SHA256
`8f3e7df51ac58adbde8952b50b9c95194752c5a4d6682a6a6d3eb94ce7d19a7d`.
The inherited exact low-degree classification was checked directly in
[`two_entry_funnels_20261103/READY.md`](../../two_entry_funnels_20261103/READY.md),
Sections1--3, with converse used for start defects. Ordered path-cover
coordinates retain their empty-secondary-slot convention throughout.

## The deleted-component classification and positive bounds

If H(T-w)=3, Hamilton counts multiply over the linearly ordered strong
components of T-w. A nontrivial strong component has a Hamilton cycle;
its cuts give at least as many distinct paths as its vertices. Therefore
a strong component of count3 must be the directed triangle. A count1
component must be a singleton. There is exactly one triangle component,
and the full deleted graph is

```
TT_a -> C3 -> TT_b.
```

This uses the inherited Camion consequence also used in READY; no claim
about a general prime Hamilton count determining a component is made.

For a,b>0, strongness forces w to beat the first prefix vertex and lose
to the last suffix vertex. Joining ordered covers at w gives the positive
contribution `12*2^(a+b-2)=3*2^(n-4)`. The12 is P2(C3), including empty
cap pieces in either slot. The extreme prefix/suffix vertices make the
whole slots nonempty, so there is no empty-slot overcount.

For exactly one nonempty tail, choose a cap port required by strongness.
The corresponding rooted C3 two-slot count is3. Free assignment of the
other tail vertices gives `3*2^(n-5)`. These are actual path injections;
neither is a signed count or a statistical estimate. In particular every
strong q3 graph of order at least14 has H at least1536, irrespective of
entry.

## Terminal cap: all vertices are accounted for

For `TT_m -> C3`, m>=3, the transitive vertices have outdegree at least3,
excluding start defect1/2 by the converse rooted-defect bound.

At the apex, strongness forces the first transitive vertex to be an
outport. One outport gives defect0. With two outports, a second transitive
outport leaves a nontransitive C3 residual; a cap outport leaves at least
two identical non11 letters in the converse READY word. Both violate the
complete exceptional language. At least three outports are excluded by
the degree bound. This treats the apex directly.

For a cap vertex beating w, the two-outport residual word similarly has
at least two identical non11 letters unless a forbidden00 already occurs.
For a cap vertex beaten by w, its sole-outport law gives exactly
`[z->w]*(C+[w->u])`, with `C>=2^(m-1)>=4`. Thus it has defect0 or at least4.
All vertex types are exhausted, and the m>=3 boundary is essential.

## Initial cap: the exception repairs the false implication

The first failed implication in the draft was extending the valid
non-apex argument to the apex. It is false that positive entry forces an
interior cap. In `C3 -> TT_b`, let w beat precisely two cap vertices and
lose to the third and every tail vertex. This is strong for every b>=1,
has apex start defect2, and has

```
H=6*2^b+1.
```

It is exactly the converse READY word `11^b 10`, so it also has an inherited
all-order consistency check. This is not a new exceptional family.

Here is an independent direct derivation, retaining the two ports. Let
the two cap outports be p,q with p->q. In a bad cover ending at p and
starting at w, the second path continues at q. The third cap vertex can
lie in either slot, giving two covers. All tail vertices must be in the
second slot. Reversing the port roles gives none. Thus sigma_w=2. There
are2 full paths starting w,3 ending w,2 internal joins at the third cap
vertex, and `6*2^i` joins at tail position i. Their sum is `6*2^b+1`.

The mixed case also admits a direct count independent of READY's bit
orientation. If the apex outports are cap vertex c and tail vertex t_i,
indexed0 through b-1, then

```
sigma_w = 3*2^i + [i=0] + 3[i=b-1].
```

For covers whose first path ends at t_i, the second starts `w,c` and the
deleted-graph join is `3*2^i`; the singleton-w second path contributes3
only when i is last. For covers ending at c, all cap vertices must be in
that first path and all tail vertices in the second starting `w,t_i`,
which is possible exactly when i=0 and contributes1. These exhaust the
two outports. The result is at least4 for b>=2. One cap outport gives0;
at least three outports are excluded by the rooted-defect bound.

For non-apex vertices, cap vertices have degree at least3. The last tail
vertex has defect0 or3. At the penultimate vertex, either the two-outport
residual retains C3, or the sole-outport formula is
`H(U+w)-end_w(U+w)>=3`. The latter follows by inserting w nonterminally
into each of U's three paths. Earlier low-outdegree tail cases retain C3
in the residual, while all others have degree at least3. Thus the stated
initial-cap classification covers every vertex for b>=2.

## The surviving all-order consequence

For strong T of order n>=7, q=H(T-w)=3, and any vertex v with start defect1
or2, the repaired conclusion is

```
H(T) >= 3*2^(n-4).
```

Interior caps satisfy the positive cover bound. Terminal caps have no
entry. Initial caps have the precise apex exception above, whose count
is larger than the bound. Therefore a strong target preserving all13
original physical vertices cannot retain H1409 or1427 and have both q3
and a positive entry, even with arbitrary extra births or changed arcs.
The final graph alone gives the contradiction; no native-guard assumption
or outer-root condition is needed. Intermediate operations are not being
classified. Together with the separate all-order q1 theorem, the result
closes these two deletion resources for that source-preserving target,
not arbitrary deletion counts or unrestricted ternary closure.

## Fixed endpoint anchor versus transported auxiliary

For the repeated C2 cap D_j, the induction is literal:

```
H(D_j)=3^j,
D_j-a0 = y1 -> (a smaller copy of D_(j-1) rooted at x1),
end_a0(D_j)=1.
```

Only x1 beats a0. Every path ending a0 comes from a path ending x1 in
that smaller copy, with y1 forced first. This proves the endpoint law;
the count law uses C2 at the transported first newborn x_j, whose start
defect is0. In a wrapper `T_j-w=U->D_j`, only a0 in D_j beats w, hence
`H(T_j-w)=3^j H(U)` while `end_w(T_j)=H(U)`. Endpoint mass stays bounded
as the deletion count grows.

The lone0 birth formula also checks by disjoint path types: paths ending
at the new vertex contribute `H-end_w`; all others contribute S_w(2).
It creates start defect end_w and changes the new endpoint at w to
H(T-w). Thus the positive entry is real but has an explicit count cost.
Applying later C2 to the fixed anchor a0 instead of transported x_j is
a different operation; at depth2 its endpoint count becomes3 instead
of1. The proof does not silently identify these physical vertices.

The initial orientation lemma additionally requires the stated r!=w.
Without it, `TT3+w` with w beating the first two vertices and losing to
the last has q1 and an end-zero root at w itself; w is not in the deleted
order. The final proof now states the distinction explicitly, consistent
with physical root99 and deleted mark0.

## Sparse-backedge sanity check

The parent's proposed cap-width observation follows from the same literal
construction. In order `x0=a0,y1,x1,y2,x2,...,yj,xj`, every edge is forward
except `x_i->x_(i-1)`. The odd-index triangles
`{x_(i-1),y_i,x_i}` are pairwise vertex-disjoint, giving feedback vertex
number at least ceil(j/2). Deleting every odd x_i removes all backward
edges and leaves a transitive tournament, attaining that bound.

Thus the cap has feedback vertex number exactly ceil(j/2), and the full
wrapper contains that cap as an induced graph. The bounded-endpoint family
cannot stay in a representation that becomes transitive after deleting
two vertices as j grows. This verifies the proposed connection; it makes
no priority claim about older sparse-defect grammars.
