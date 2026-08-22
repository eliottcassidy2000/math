# The first residue packet closes by a birational root plus a quadratic

**Status: proved synthesis from THM-2473, THM-2570, THM-2576, THM-3530,
THM-3538, THM-3539, and THM-3540.**  The depth-two newest-prime
decomposition image is point-and-pair saturated.  No full decomposition-
centralizer equality, depth-three saturation, or all-level index result is
claimed.

## 1. Inheritance pass

The closest mechanism is THM-3539's local sandwich

```text
I <= D <= C_W(I),                                      (1)
```

whose missing coordinate is the residue action `D/I`.  The canonical hostile
is the trivial predecessor action, which leaves all six depth-two THM-3538
factors separate despite full global wreath monodromy.  The corrected near
miss is that strict henselization cannot answer the question: it splits the
predecessor cover by erasing exactly this residue action.  The least-used
sidecar is THM-2570's normalization axis of the old `L` surface.

At depth two those ingredients happen to fit perfectly.  The image map
`V(L)->V(H)` is birational, so it supplies one rational predecessor over the
generic point of `H`.  The other two predecessors are controlled by one
quadratic discriminant.

## 2. The `1+2` anatomy is now a residue theorem

Pull the inverse cubic over `H` back to the normalized old component.  Its
birational ancestry point has `x`-coordinate `X`, and exact division gives

```text
E(W)=(W-X)Q(W).                                        (2)
```

The product-discriminant identity and THM-2473's cubic discriminant give

```text
[Disc(Q)]=[-L(target)].                                (3)
```

On the normalization divisor `lambda=0`,

```text
(0,0,tau) -> (tau,0,0),
L(target)=16tau,
Disc(Q)=-256tau.                                       (4)
```

The odd `tau` valuation in `(4)` proves that `(3)` is nonsquare, so `Q` is
irreducible.  The residue action therefore fixes the ancestry block and swaps
the other two.  This is the full `S2` stabilizer of a point in `S3`.

The familiar `1+2` motif has acquired a precise new role: one is the
birational component section, and two is an irreducible residue orbit.  It is
not merely a collision picture or Chebyshev factorization.

## 3. Four packets, but not the Boolean `V4`

The marked-point `S2` has two block orbits and two unordered-pair orbits:

```text
blocks:          {marked}, {off_1,off_2};
pairs:           {marked-off_1,marked-off_2}, {off_1-off_2}.        (5)
```

Hence the six raw THM-3538 factors descend to the four valuation types

```text
i_(theta,2)
 =v(A_marked)+2v(A_off)
  +2v(R_marked,off)+v(R_off,off).                      (6)
```

The number four should not be merged with the XOR carrier.  In the Boolean
atlas, four states form the group `(Z/2)^2` and have a bidirected Cayley
square.  Here four is only the count of two point-orbit types plus two
pair-orbit types.  There is no packet addition law, no canonical `C4`, and no
tournament orientation.

This is a useful anti-numerology control: equal cardinality does not identify
operations or preserved predicates.

## 4. Relation to the three cubic discriminants

The three literal coordinate cubics at the first Keller level share the
square class `[-L]`.  THM-3540 shows a second appearance of the same class:
after one birational root is removed from the predecessor cubic over `H`, the
remaining quadratic also has class `[-L(target)]`.

The common mechanism is sign parity under removal of a rational sheet.  A
linear factor contributes only a square resultant to the discriminant, so
the quadratic retains the cubic's sign class.  What changes is the base:
the first statement lives over the ambient target field, while the residue
statement lives over `kappa(H)` and uses the normalized ancestry section.
They are compatible views, not the same extension.

## 5. What changed and what remains hard

THM-3539's conditional `n^2` compression is now unconditional at `n=2`.
THM-3538 already proved unitness there, so the whole depth-two local passport
can be read on four orbit representatives.

At `n=3`, the analogous predecessor algebra has degree nine.  One birational
ancestry section still exists, but “the other eight roots are irreducible” is
neither necessary nor obviously true.  The economical target is only its
orbits on the nine blocks and 36 pairs.  Full marked-LCA saturation would
give

```text
3 block types + 6 pair types = 9 packets.              (7)
```

A factorization that splits an LCA shell would be equally valuable: it would
identify the extra residue coordinate that the strict-henselian atlas loses.

## Connection contract

| field | exact answer |
|---|---|
| source | normalization of `V(L)` and birational map to `V(H)` |
| target | predecessor-block residue action at the newest depth-two prime |
| map | remove the rational ancestry root from the inverse cubic |
| preserved | discriminant square class, marked block, point/pair orbits, valuations |
| destroyed | last-step root labels, full residue kernel, unit residues |
| decisive sidecar | `lambda=0` divisor with odd `tau` valuation |
| proved gate | THM-3539 point-and-pair saturation at `n=2` |
| still open | full `D=C_W(I)`, `n>=3` saturation, `n>=5` unitness |

This result changes neither `JC(2)` nor the arbitrary-map classification.
It supplies no current, tournament equivalence, or LRC(14) implication.
