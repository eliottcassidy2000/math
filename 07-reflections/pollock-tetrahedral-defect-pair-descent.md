# Pollock tetrahedral defect-pair descent

The productive turn was to stop looking directly for five tetrahedra.

Let `D_4` be the set of integers that are not sums of at most four
tetrahedral numbers.  Pollock says every number has one more tetrahedron to
spend.  If `n` sits in the shell `[Te_k, Te_{k+1})`, write `n=Te_k+r`.
Then either `r` is already four-tetrahedral, or we can step one anchor back:

```text
n = Te_{k-1} + (r + tri(k)).
```

So the obstruction is not a random failure of five summands.  It is a pair of
four-sum defects separated by a triangular number.  That feels like a real
object.

The computation found the known frontier again: through `10^6`, exactly `241`
numbers need five tetrahedra, with largest `343867`.  Among those defects, the
last triangular separation is

```text
3142 -> 343867 = 3142 + tri(825).
```

The shell stencil makes that statement geometric.  For every shell through
`k=1200`, one of the anchors `Te_k, Te_{k-1}, Te_{k-2}, Te_{k-3}` plus four
smaller tetrahedra covers the shell.  After `k=825`, only `Te_k` and
`Te_{k-1}` are needed.  This is exactly the kind of split one wants in a proof:
a finite certificate plus a tail lemma.

The offset tournament was also telling.  Offsets `3 -> 2 -> 1 -> 0` form a
transitive dominance order by private shell coverage.  Older anchors win
because they shift the remainder deeper into the four-sum set, where coverage
is denser.  But the proof does not need the strongest old anchor generically;
it needs the one-back anchor to stop seeing defect pairs after the last
triangular correlation.  Dominance is a guide, not the invariant.

This connects back to the repo's irreducibility and tiling language.  A
four-tetrahedral representation is a convolution lift from four cubic atom
indices to a boundary total.  A defect is a boundary total with no lift.  The
one-back descent asks whether two no-lift boundary totals can be separated by
the special triangular carrier `Te_k-Te_{k-1}`.  That is much closer to a
tiling obstruction than to a naked Waring statement.

The strongest route is probably still analytic: prove every sufficiently large
integer is a sum of four tetrahedra, then verify the finite tail.  But the
weaker route is more Pollock-shaped and may need less: prove only that the
four-defect set has no long triangular self-correlation.  It is plausible that
modular information useless for explaining single four-defects becomes useful
for explaining defect pairs.  In the first residue spot checks, four
tetrahedra already covered every residue class for many small moduli, so the
single-defect obstruction is not local.  The pair obstruction might be.

The next scout should enumerate defect-pair residues, not just defect residues:
for each modulus `m`, look at pairs `(d, d+tri(k))` with both endpoints in
`D_4`, and rank moduli by how fast they destroy the observed pair classes.  If
there is a small family of moduli that kills all `k>825` triangular pair
templates, Pollock becomes a finite certificate plus a modular no-pair lemma.
