# Independent audit: abstract F20 embedding and restriction

**Disposition: PASS, conditional on the explicitly cited Baby Monster subgroup,
class and character inputs.** No Baby Monster matrices, standard-generator
words, or representation intertwiner were constructed by this audit.

## Product notation and class selection

The [ATLAS B entry](https://brauer.maths.qmul.ac.uk/Atlas/v3/spor/B) displays
`5:4 x HS:2`, of order 1,774,080,000. The
[Breuer–Magaard–Wilson protocol](https://www.math.rwth-aachen.de/~Thomas.Breuer/ctbllib/htm/ctblbm.htm)
uses table identifier `5:4xHS.2`; its nearby prose omits the final 2. These
were independently re-opened. The [ATLAS HS entry](https://brauer.maths.qmul.ac.uk/Atlas/spor/HS/)
independently confirms order 44,352,000.

The final producer proof correctly avoids deciding whether an outer involution
commutes with the first factor. It uses only the cited normal index-two base
`L=F20 x HS` in the listed subgroup M. Since HS is nonabelian simple,
`O5(L)=C5`; therefore the distinguished C5 is characteristic in L and normal
in M. A nonidentity element x of that C5 centralizes `C5 x HS`, of order
221,760,000. This exceeds the cited 5B centralizer order 6,000,000. With only
two order-five classes, x must be 5A.

Now `C_B(C5)=C_B(x)` has cited order 443,520,000. Conjugation embeds the
normalizer/centralizer quotient in Aut(C5), of order four. M already lies
in that normalizer and has order exactly four times this centralizer order,
so equality of the normalizer with M follows. This is a valid saturation
argument; neither maximality nor a direct-product description of the entire
centralizer is needed. These are consequences of cited group facts, not facts
proved merely by the arithmetic checker.

## Abstract embedding

The earlier explicit puzzle group has order twenty and generators g,h
satisfying the standard C5:C4 presentation with exponent two. The cited F20
factor admits generators a,b with the same relations. The common presentation
has at most twenty normal forms and has these faithful twenty-element
realizations. Thus g->a and h->b define an abstract isomorphism onto that
factor, followed by its inclusion in B. The image of g is 5A.

This proves abstract subgroup existence, not a computed embedding into a
particular Baby Monster matrix realization. It does not identify the full
thirteen-point puzzle action or arbitrary M12/M13 elements with a fixed
five-point action of B.

## Character and multiplicity audit

An independent implementation enumerates F20 as actual five-point permutations
by generator closure. It verifies the quotient C4 and constructs the four
linear characters with values in the Gaussian integers. The fifth character
is the actual deleted permutation module, with values 4 at identity, -1 on
nonidentity translations, and zero for nontrivial multipliers.

Exact character inner products form the identity matrix. The squared degrees
sum to twenty, proving the complete complex irreducible degree list
`1,1,1,1,4`. In particular there is no ten-dimensional complex irreducible of
this specified group.

The cited Baby Monster degree 4371 and 5A value 21 then force 870 copies of
the four-dimensional module and 891 one-dimensional constituents in total.
The separate four linear multiplicities are not determined. The complement
values `1,i,-1,-i` show that rationality forces precisely the additional
constraint `a1=a3` used in the report.

Crucially, a five-dimensional affine permutation module is `1_H + W`.
Having 870 copies of W does not ensure 870 copies of this whole module:
it would require at least 870 trivial H-constituents. The formal rational
H-module `870W + 891lambda_2` has the same known C5 data and no trivial
constituent. This is a counterexample to the proposed inference from those
data, **not** a claim that this formal module is the actual B restriction.

## Verification boundary

The independent checker executes 447 explicit checks with exact integer-pair
arithmetic; ordinary and optimized outputs match. Three mutations are rejected
in both modes: invalid 5B-centralizer input, corrupted deleted character, and
an unjustified whole-permutation-module conclusion. The producer's frozen
manifest and its normal/optimized baseline and hostile replays are also
verified. No current checker uses optimization-removable assertions.

`reviewed_source_pins.json` binds the exact producer files and earlier puzzle
and sporadic manifests. The prior scientific packets are untouched. This
addendum closes only this specified existence/restriction corollary; no new
geometric, Hecke, FRACTRAN, or universal sporadic correspondence is asserted.
