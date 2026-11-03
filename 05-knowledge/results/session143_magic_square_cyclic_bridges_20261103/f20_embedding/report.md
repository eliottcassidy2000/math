# A scoped abstract embedding into the Baby Monster's 5A normalizer

**Status:** PROVED from the cited ATLAS subgroup/centralizer facts and the previously verified puzzle `5:4`; FINITE-EXACT arithmetic and ordinary-character checks. No explicit Baby Monster generator certificate, intertwiner, or originality claim. Earlier frozen packets are unchanged.

## 1. The cited factor necessarily contains 5A elements

The already checked [ATLAS Baby Monster page](https://brauer.maths.qmul.ac.uk/Atlas/v3/spor/B) lists a subgroup of order `1,774,080,000` with displayed label `5:4 x HS:2`, and exactly two order-five classes, with centralizer orders

`|C_B(5A)|=443,520,000`, `|C_B(5B)|=6,000,000`.

**Parentheses discipline:** the already downloaded [Breuer–Magaard–Wilson verification protocol](https://www.math.rwth-aachen.de/~Thomas.Breuer/ctbllib/htm/ctblbm.htm) uses table identifier `5:4xHS.2`, while nearby prose abbreviates the normalizer as `5:4 x HS`. We do not infer that the final outer involution commutes with the `5:4`, and do not assert a direct-product structure for the whole centralizer. The argument needs only the index-two base

`L=(5:4) x HS`, with `|HS|=1,774,080,000/40=44,352,000`,

inside the listed subgroup `M`. Both possible outer-involution readings have this base. Let `C` be the normal order-five subgroup of its `5:4` factor. It is the largest normal 5-subgroup of `L`, hence characteristic in `L`; as `L` is normal in `M`, so is `C`.

For every `1!=x in C`, its centralizer contains `C x HS`, of order

`5*44,352,000=221,760,000 > 6,000,000`.

Thus `x` cannot be in `5B`, and must be in `5A`. We now use the cited `5A` centralizer order, without identifying its group structure.

Conjugation gives an injection

`N_B(C)/C_B(C) -> Aut(C)`, with `|Aut(C)|=4`.

Since `x` generates `C`, `C_B(C)=C_B(x)`. Thus `|N_B(C)|<=4*443,520,000=1,774,080,000`. But `C` is normal in `M`, so `M<=N_B(C)` and its order attains that bound. Consequently **`N_B(C)=M`**. Maximality is not required for this argument.

## 2. The explicit puzzle group embeds abstractly

The separately checked projective-puzzle permutations satisfy

`g^5=h^4=1`, `hgh^-1=g^2`, and `|<g,h>|=20`.

The cited first factor of the index-two base `L` is the same faithful semidirect product `C5:C4`. Choose its generators `a,b` with `bab^-1=a^2`; the order-four generator can be inverted if necessary. The assignment `g->a`, `h->b` therefore gives an **abstract injective homomorphism** `<g,h> -> B`, and the image of `g` lies in class `5A` by section 1.

This is existence relative to the cited subgroup, not a computed word in the standard generators of `B`. It does not identify the puzzle's thirteen-point action with a Baby Monster action, nor make the five-tournament an orbital diagram of `B`. The order-four puzzle element acts on the four tournament orientations, not as an automorphism of one selected orientation.

## 3. The 4371-dimensional restriction to this factor

Let `H=5:4`, `C=<a>`, and let `rho` be the ordinary complex Baby Monster representation of degree `4371`. The already checked [Breuer–Magaard–Wilson character verification](https://www.math.rwth-aachen.de/~Thomas.Breuer/ctbllib/htm/ctblbm.htm), sections 6–7, supplies `chi(5A)=21` and the rationality data. Its restriction to `C` has trivial-character multiplicity `891` and multiplicity `870` for each of the four nontrivial linear characters.

The complex irreducibles of `H` have degrees **`1,1,1,1,4`**. Indeed the four linear characters inflate from `H/C=C4`. Inducing any nontrivial character of `C` gives a degree-four character `W`; the complement acts freely transitively on the four such characters, so Frobenius reciprocity gives `<W,W>=1`. The squared degrees sum to `4+16=20`, exhausting all irreducibles. Equivalently, `W` is the deleted five-point affine permutation representation: it has values `4` at the identity, `-1` on nonidentity translations, and `0` on elements with nontrivial multiplier.

Each linear irreducible restricts trivially to `C`; `W` restricts as the sum of its four nontrivial linear characters. Hence

`Res_H(rho) = W^870 + lambda_0^a0 + lambda_1^a1 + lambda_2^a2 + lambda_3^a3`,

where `a0+a1+a2+a3=891` and each `ai>=0`. Here `lambda_t(a)=1`, `lambda_t(b)=i^t`. Rationality additionally forces **`a1=a3`**, since the nonreal linear characters are complex conjugates. The separate `ai` are not determined by the imported order-five data; further information about the chosen order-four element is needed. No new character-table lookup is made.

**Whole permutation-module boundary:** the five-dimensional affine permutation module is `1_H + W`. Having `870 W` does **not** imply `870` copies of that whole module: this would additionally require `a0>=870`, which is not determined. For a concrete non-inference witness, the rational `H`-module `870 W + 891 lambda_2` has the same imported `C5` restriction but `a0=0`, so it contains no whole affine permutation module. This is a valid abstract `H`-module, **not** a claim that it is the actual Baby Monster restriction. On `C5` alone, `870` regular modules plus `21` trivial modules is a valid decomposition.

Thus a proposed **10-dimensional ordinary complex irreducible of this `5:4` subgroup is impossible**. This does not exclude reducible ten-dimensional representations, or make a claim about an unspecified different stabilizer. Over `Q`, the two nonreal linear characters combine into one degree-two irreducible; the four-linear-character statement is specifically over `C`.

## Verification and scope

`python3 check_corollary.py` executes 77 explicit checks of the arithmetic, builds the exact affine twenty-element group, and verifies all five character functions and their orthogonality using integer-pair Gaussian arithmetic. It confirms the restriction dimensions and runs unchanged under `-O`. The optional mutation `--mutate bad-centralizer` must exit nonzero in either mode. `python3 test_corollary.py` passes all four normal/optimized positive/mutation cases; the harness itself also produces identical output under `-O`.

External facts remain CITED; the checker does not construct `B` or independently derive its character table. The new conclusion is a typed chain of **abstract subgroup existence plus representation restriction**, not a canonical geometric or spectral identification.
