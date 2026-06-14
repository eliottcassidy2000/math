# Pollock, Sierpinski, and the carry-pair lift

The Sierpinski instinct was worth following, but it turned inside out.

For Waring-style problems, the classic lower-bound picture is a carry trap:
numbers just below the next large atom can force many smaller atoms. So the
natural Pollock question was whether the tetrahedral four-defects are hiding a
dyadic/Sierpinski obstruction.

The scout says: not at the single-residue level. The tetrahedral atom
`Te_k=C(k+2,3)` hits every residue class modulo `2^e` for every `e<=12` in the
scan. Lucas parity is still there, cleanly:

```text
Te_k is odd iff k == 1 mod 4.
```

But that Sierpinski-looking parity fractal is not a missing local class. One
tetrahedral atom already fills the dyadic residue circle.

The lift from HYP-2491 is where the structure returns. Pollock is not asking
whether a single four-defect exists. In shell form, it asks whether two
four-defects can appear at the special separation

```text
tri(k)=Te_k-Te_{k-1}.
```

That is a two-endpoint carry condition. Once the scout measures pair residues
instead of single residues, dyadic levels start behaving like a real proof
coordinate: for pairs with `k>=100`, observed pair classes saturate at `168`
by `2^8`, while the full pair universe keeps growing as `4^e`. The resulting
dyadic-level tournament is completely transitive:

```text
12 > 11 > 10 > 9 > 8 > 7 > 6 > 5 > 4 > 3.
```

So the moral is not "dyadic methods fail." It is more precise: scalar dyadic
methods fail, lifted pair/carry dyadic methods may still be exactly right.

The carry-window data fits the same split. Many four-defects are close below
the next tetrahedral number: `85/241` are within `100`, and `240/241` are
within `5000`. That is the Waring/Sierpinski smell. But the largest defect
`343867` is `5637` below `Te_127`, so "near the next atom" is a powerful
diagnostic rather than a complete theorem.

The proof route I now trust most is:

1. Prove the 2-adic surjectivity lemma for single tetrahedral atoms. This
   intentionally shuts down a tempting but false local-obstruction route.
2. Keep HYP-2491's one-back descent as the main theorem:

```text
r and r+tri(k) cannot both lie in D_4 for k>825.
```

3. Attack that no-long-pair lemma with a lifted ledger: dyadic pair address,
   triangular gap address, and four-atom convolution feasibility.
4. Pair the tail lemma with the finite width-3 stencil through `k<=825`.

This feels close to the irreducibility work in the repo. A coefficient row is
irreducible not because one coefficient has a forbidden residue, but because
no hidden convolution lift exists. Here a Pollock obstruction is not a single
bad integer; it is a pair of boundary totals with no four-tetrahedron lift,
separated by the triangular shell carrier.

That is the shape I would push next: build a SAT/ILP or modular convolution
oracle for pairs, not for individual defects. The endpoint pair is the object.
