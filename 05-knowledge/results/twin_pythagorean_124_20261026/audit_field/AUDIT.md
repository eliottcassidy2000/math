# Independent audit of the paired finite-field design

**Verdict: PASS.** The proof, reciprocal slope conventions, named block
systems, exact groups, branch formulas, and finite continuation boundary
were independently checked. No mathematical defect remains. The producer
was not edited by the reviewer.

## Frozen files and replay

Reviewed producer files in `../finite_field/`:

| File | SHA-256 |
|---|---|
| `README.md` | `38ea6edc4f92c8d132d16efe1cbbaa4f006e7499baa0a8ec9d2bad91a2e33413` |
| `check_finite_field.py` | `f3082927900a979924590aec0806f9d38516690a29ec192025eb12e205335730` |
| `receipt.json` and optimized receipt | `b84c315371893ffd6910b8a479adf37c30da0bb19b440e63339484ef91b68f8c` |

The producer passes **3,616 gates** under ordinary and optimized Python,
both byte-identical to its final receipt. The independently written
`check_independent.py` passes **1,902 gates**, also byte-identical in both
modes, with receipt SHA-256
`630bc2911cdf526083ed094c67ebd417ac4cf8d8ec41f3f1b4572707e17761ab`.

One nonmathematical correction was requested and made before final freeze:
the producer's receipt now refers to the actual proof file `README.md`
instead of `RESULT.md`. The initial checker hash was
`26e01e9b3e5c0d13299d161fbafcccaf210af2c6fd363bfa7e3954364686438b`.
No arithmetic or proof changed in that amendment. The final proof's THM-1370
link was reread and has the correct `h-spectrum-omits-7-21-all-n` slug.

## Proof and type review

### Coordinates, point set, and marking

The common slope `r=n/m` and the auxiliary slope `t=m/n` are reciprocal,
with zero and infinity interchanged. The root gives `r=4,t=2` in F7;
the endpoint hostile gives `r=t=6`. Both homogeneous parametrizations in
the note are correct. Since `-1` is nonsquare modulo seven, normalization
by the last coordinate is valid at every conic point, including the
separate infinity parameter. The removed pole is not the integer tree root.

The exponential marking `E(t)=alpha^t`, `E(infinity)=0` has the exact listed
eight values. In the r chart it is `E(1/r)`. It is correctly typed as a
chosen set bijection with a proved incidence compatibility, rather than a
field map or circle-group homomorphism. The residue companion's witness
uses `r=2,t=4`; this is consistent with the same chart conventions.

### Two designs and the determinant character

The unique completion of three distinct F8 points by their sum proves the
Steiner property for D+. Cubing is bijective on F8, proving it for D-.
The characteristic-two identity

    x^3+y^3+z^3+(x+y+z)^3=(x+y)(x+z)(y+z)

is correct. Distinctness makes its right side nonzero and proves disjointness
of the two block systems. Thus there are 14 blocks in each, not two names
for one design.

For the generators, translation becomes multiplication by alpha and
`S:t -> -1/t` becomes `x -> L(x)+1` with the stated binary-linear L.
Its values at both zero/infinity parameters are included. L is not a field
automorphism; the note explicitly retains that boundary. The standard
translation and inversion generators generate PSL2(F7), either by elementary
matrix elimination or the complete finite closure checked here.

For `J:t -> 3t`, the transported map is cubing. Its square is Frobenius,
which preserves D+. Together with D- being the inverse image of D+ under J,
this proves that J exchanges the two systems. Normality of PSL then gives
preservation of D- by PSL. The determinant square class is well-defined
projectively and its two cosets give precisely preservation versus exchange.
No implication reverses the role of a point marking and a field operation.

### Image groups versus full automorphism groups

Every nonsingular 2-by-2 matrix over F7 was independently normalized by its
first nonzero coefficient. This gives 336 distinct projective permutations,
of which 168 have square determinant. Their action on all 70 four-subsets
has orbit sizes `14,14,42`. The two small orbits agree exactly with the
producer's named D+ and D- after evaluating its exponential marking.

A separate exhaustive oracle over all `8!` permutations obtains full
automorphism counts `1344,168,336` for one design, the ordered pair, and the
unordered pair respectively. The producer additionally compares the full
sets of permutations; its membership tests are correct. The proof that
the one-design group is AGL3(F2) is valid: after translation to fix zero,
preservation of every plane `{0,x,y,x+y}` forces additivity, including the
trivial equal/zero cases separately. An affine marking independently found
by this audit makes every square-determinant projectivity affine.

The note therefore does not mistake the 168-element generated subgroup
for the 1344-element full one-design group. Its full paired-group equality
has an explicitly complete finite universe, rather than an order-only test.

### Branches, lost endpoint bit, and pole restriction

The U,A,D matrices agree with the integer tree convention. Conjugating
by coordinate inversion gives the stated r-chart formulas. U and D have
square determinant and generate PSL; A has nonsquare determinant and
enlarges it to PGL. The polynomial formulas over F8 agree, including the
zero cases of the sixth power. Counting A letters gives the swap character.

The integer words UU and A give the stated different primitive triples,
the same projective endpoint, and opposite swap bits. Their affine triple
residues differ by scalar two, so equality is projective, not an equality
of all three residues. The four states at depths zero and one are distinct;
the depth-two hostile is therefore genuinely the first such collision.
The claimed 16 reachable point/bit states are checked exhaustively. A state
does not recover the full branch word or its integer height; the text does
not claim that it does.

The pole stabilizer and Paley statements agree with inherited THM-4553,
which was reread. Square affine slopes preserve square differences and
nonsquare slopes reverse them. The derived Fano blocks, full Paley group,
and Hamilton count are all tested on complete stated finite universes.
Two-transitivity prevents an invariant tournament on all eight points.
The signed circle determinant is a different observable: its formula and
the surviving antipodal zeroes are correct, so no cosmetic tournament
orientation is imposed on those ties.

### Form and larger-field hostiles

The F7 norm circle has an element of order eight, whereas additive F8 has
exponent two and multiplicative F8 has order seven. The stated group and
homomorphism obstructions follow. Removing one point does not make the
remaining seven points a subgroup of the order-eight circle.

The naive characteristic-two quadratic is indeed a double line with zero
gradient. The other conic has gradient `(z,0,x)`, which cannot vanish at a
projective solution. Equal counts of nine rational points do not identify
the schemes; the note explicitly keeps this distinction.

For fields of sizes 4,8,32, the audit independently implemented polynomial
coefficient convolution and long division, constructed the complete power
tables, and tested every primitive-generator exponential marking. The
results are respectively `2/2`, `6/6`, and `0/30` affine transports of S.
For the displayed generator at 31, `{0,1,2,3}` maps to `{1,0,18,14}`,
whose XOR is 29. This matches the producer's literal failed plane.
The note carefully limits the result to those three fields and those
markings. It does not classify all Mersenne primes or arbitrary bijections.

## Independent path and closure

The independent script imports no producer code. Its main universes are all
normalized projective matrices, all four-subsets, all `8!` point permutations,
the entire projective F7 quadratic locus, and all primitive markings in the
three specified binary fields. It derives the block systems first as PSL
orbits, then identifies the author's signs by the explicit marking. This
checks labels as well as cardinalities.

Run from this audit directory:

```bash
python3 -B check_independent.py
python3 -B -O check_independent.py
```

The paired-design proof is self-contained elementary finite-field algebra.
The all-parameter integer twin-prime predicate is not inferred from it;
the separate CRT obstruction explains precisely why fixed residue data
loses that predicate. No twin-prime, abc, RH, or universal Hamilton-spectrum
conclusion is introduced by this audit. A later status/link-only audit
promotion is covered as lineage; a mathematical change requires review.
