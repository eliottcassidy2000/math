# The701 return changes ready lane

**PROVED all-size ready-image consequence of the inherited classification
and a complete finite head; FINITE-EXACT actual native return.** Local
labels R1--R4. Same-agent proof review and independent counting controls.

## R1 No one-ready source and exactly one strong two-ready class

Use h=H(T), s=#Hamilton paths starting at r, q=H(T-r), and F=E_r(2),
the number of ordered spanning two-bin path covers whose first path is
nonempty and ends at r. A b-ready root has F-h=b. The inherited NR reset
then produces h'=3h+2b.

At h=701 the complete classification is:

* No finite tournament has a1-ready root and H701.
* There is exactly one rooted strong2-ready isomorphism class with H701.
  It has order9 and rooted profile `(701,81,83,703)`.
* Every weak2-ready example is that strong core followed by a transitive
  tail. There is one rooted class at each order n>=9, including the strong
  example at9. This does not classify arbitrary second-port markings.

**Proof and complete finite reduction.** [READY](../two_entry_funnels_20261103/READY.md)
classifies every1-ready root as the one-in-neighbour newborn of a strong
one-start parent. A2-ready root is either the analogous two-start newborn
or one of five explicit exceptional families. The newborn count is
`h_child=h_parent+F_parent-b >=2h_parent-b`.
For h_child=701 and b=1 or2, h_parent<=351. The sharp strong-order bound
in [Busch2006, Corollary1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf)
gives minimum H375 at order12 and increasing minima thereafter. Parent
order is therefore at most11. This is a cited published bound.

The inherited unique-start rank/unrank grammar has counts
`1,3,8,21,55,144,377,987,2584` at orders3--11:4180 parents total.
None produces701. The independent two-start grammar has counts
`0,1,5,18,60,191,589,1775,5257`:7896 parents total.
Exactly one produces701. Its parent, in the grammar's canonical frame, is

```
(2,4,25,147,227,207,143,7), root0,
(h,s,q,F)=(83,2,49,620).
```

The newborn graph is

```
(258,4,25,147,227,207,143,7,254), root8.
```

Its count is83+620-2=701 and its defect is2. The five exceptional values
`3*2^k-1,6*2^k+1,5*2^k-1,13*2^k+3,14*2^k+3`, k>=1, omit701.
For each family it suffices to check until its first value exceeds701,
since it increases strictly. These computations use the constructive
grammars, not a sampled tournament census; [ready701.py](ready701.py) and
[two_ready701.py](two_ready701.py) cover every parent allowed by the bound.

The SCC argument in [O1](../ready_return_pentagon_20261103/OBSTRUCTION.md)
shows that a positive ready defect1 or2, smaller than H, puts the root in
the first strong component and forces every later component to be a
singleton. This proves the weak-source statements and completes R1.

Consequently an H-preserving return into the1-ready lane is impossible
at701, regardless of the proposed operation. The2-ready lane remains
available and is used below. This is an interface classification, not a
missing Hamilton-path value.

## R2 The actual retained-source return

Start with the saved ready233 source and its NR child in
[the preceding packet](../joint_port_geometry_20261103/RETURN_233.md).
The701 child's physical names are

`[0,1,2,4,5,6,7,8,13,14]`, with root13 and spare6.

The complete program in [native_return.json](native_return.json) performs:

1. One native complement off a supplied legal Hamilton path: H701->1373.
2. Delete original name0: H1373->193.
3. Six isolated old-arc reversals, each implemented by two guarded native
   complements with a temporary vertex, giving counts
   `193->625->723->747->803->905->701`.

The final root is still name13 and the spare is still name6. Original
names4,5,7 are also retained. The returned graph is strong and has profile
`(701,81,83,703)` at that actual root. The program uses13 native moves,
six temporary helper births/deletions, one old-vertex deletion, and peak
order10. H is not constant during the program. The initial native graph
and every post-repair endpoint are recounted in the receipt; strongness
is required for the deletion endpoint and each completed arc repair.

**Arc-reversal mechanism.** Choose a Hamilton path P containing the arc
a->b, with b not its last vertex. Add a helper h which beats only a.
The first native path is P followed by h; after complement, the second
inserts h between a and b. Both are legal. Their old consecutive-edge
sets differ only at a--b, so composing the complements and deleting h
reverses exactly that old edge. Every path is checked before use. This
extends the earlier first-edge implementation to any nonterminal arc.

The search compares all701 supplied native paths, five deletions outside
the protected names, and8! rooted identifications with the nine-vertex
target:141,321,600 pairs. The minimum old-arc discrepancy is six.
That is a minimum in this fixed target/one-complement/one-deletion bank,
not a global lower bound on native moves. A guarded strong repair order
is then found and replayed literally. Target identifications specify
which old names carry which incidences; they are not free graph edits.

## R3 Actual branches survive the obstruction

The returned2-ready root has an ordinary NR continuation to2107=3*701+4,
with the first new vertex as the next root and both new names recorded.

There is also a different one-complement/one-birth escape from the original
701 source to **each** of2103,2105,2107, with the same old root13 and spare6
retained. Complete enumeration gives38,43,20 strong labelled programs,
respectively. These count programs, not unlabelled classes. The selected
programs, all masks, full endpoint rows and charts are saved in
[ternary_escapes.json](ternary_escapes.json). Their large rooted defects
do not constitute a recursively ready bank.

The earlier complete one-complement/one-birth bank returning H701 had29
programs and no ready root. Its pruning H(native)>701 is justified by
Hamilton-path insertion monotonicity. A one-deletion search likewise did
not find the desired retained-root ready interface without further repairs.
These finite failures motivated R2; R1 explains the failure of the first
ready lane at every size.

## R4 Source phases and clocks

The retained root's golden coordinates change from
`(701,3634,3648)` to `(701,1042,1071)`; its F7^3 chart changes from
`(1,1,0)` to `(3,0,1)`. No phase neutrality is claimed. All five Q4
presentations, the full eighteen-polynomial two-port endpoint observer,
the exact native paths, source names and helper/deletion clocks are saved.
An independent endpoint DP recounts the original, returned and next source.
The81/27 hostile pair is replayed separately, without identifying its
automorphism order81 with the returned source's starting-path count81.
