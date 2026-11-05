# A small physical context whose positive matrix is not totally nonnegative

Session 197. **PROVED counterexample.** Universal entrywise positivity in central-factorial coordinates does not imply total nonnegativity, even for a parity-preserving actual tournament substitution context. This separates the general positive context calculus from the reviewed periodic D/X Stieltjes construction.

## 1. Physical context and three literal clone profiles

Let

    C[T]=C3[T,Tr_3,Tr_2],

where the quotient is the directed triangle T->Tr_3->Tr_2->T. For a singleton input, in block order the actual six-vertex rows are

    [14,60,56,48,33,1], root0.

The rooted context has a unique full bin-profile linear operator A. This is the reviewed whole-module cover-linearity theorem; it is not inferred from the following three values alone. Since B_(Tr_m)=x^m, the following actual clone profiles give its monomial columns:

    A(x)   = x^6+12x^4+10x^2,
    A(x^3) = x^8+36x^6+240x^4+232x^2,
    A(x^5) = x^10+60x^8+830x^6+2862x^4+1580x^2.           (1)

These profiles belong to graphs with respectively6,8,10 vertices. Their Hamilton counts are23,509,5333. The accompanying verifier constructs the graphs and derives their full profiles by exact subset path covers; a second independent cyclic-skeleton computation gives the same columns.

## 2. Exact negative minor

Use the usual parity bases

    q_k=x product_(j=1)^k(x^2-j^2),
    r_k=x^2 product_(j=1)^k(x^2-j^2).

The context adds five vertices and maps odd q columns to even r columns. Linearity and (1) give exactly

    A q_0 = 23r_0+17r_1+r_2,
    A q_1 = 486r_0+424r_1+49r_2+r_3,
    A q_2 = 2880r_0+6220r_1+1571r_2+85r_3+r_4.

Thus the leading3-by-3 block of its column-action matrix is

    M = [ 23   486   2880
          17   424   6220
           1    49   1571 ],

and

    det(M) = -468310.                                    (2)

Every displayed coefficient is positive, consistently with universal entrywise positivity. Nevertheless (2) disproves total nonnegativity.

## 3. A parity-preserving physical counterexample

Append one final transitive vertex to the output: C'[T]=C[T]->K1. This is again a single actual rooted context. Its singleton-input rows are

    [78,124,120,112,97,65,0], root0.

The operator is XA, where X multiplies by x, and it now maps the odd basis to itself. Since

    Xr_k=q_(k+1)+(k+1)^2q_k,

the first3 rows of XA depend only on the first3 rows of A, with the exact left factor

    P_3 = [1 0 0
           1 4 0
           0 1 9].

The leading block is therefore

    P_3 M = [23   486   2880
             91  2182  27760
             26   865 20359],

whose determinant is36 det(M)=-16859160. No truncation term from a higher basis index is omitted, because X maps index k only to k or k+1.

## 4. What the counterexample establishes

The reviewed periodic D/X theorem factors every parity-returning word in triangle wrapping D and transitive one-vertex padding X into TN matrices. Therefore C' cannot have the same full-profile operator as any such word. This is an exact physical operator-expressibility obstruction, even though C' is an allowed general rooted substitution.

It does not follow that the repeated-vacuum Hamilton sequence of C' has a negative Hankel minor. Non-TN production matrices can still have TN output moment arrays. Nor does this minor obstruct a single equal scalar count or an unrelated physical construction. The counterexample is precisely to extending the TN matrix claim from the specified D/X monoid to all actual contexts.

No minimality or classification claim is made. In particular the proof does not depend on a census of smaller contexts.

## 5. Reproduction and dependencies

Run `python3 verify_physical_non_tn_context.py`, and again with `-O`. The checker validates both tournament row lists, constructs the three actual clone graphs, computes all their bin coefficients, checks the independent literal skeleton operator on the three q columns, and evaluates both exact determinants. The zero and positive clone boundaries are also compared through m=6.

Named frozen dependencies:

- `UNIVERSAL_CONTEXTS.md`, SHA256 `8615845e7bf06a6cdb1925fdbd8198a77b1a8864aa9c22704b0459678ff2ed51`: full-profile context linearity and positivity.
- `PERIODIC_STIELTJES_WORDS.md`, SHA256 `82b941c9ce53b41ca73ee0224b240fe59c61df701df5fde6f9d29283f0fad166`: the TN word theorem used only for the non-expressibility consequence; published at `8c181a199d3f699347a3c76a885277d6324e2bd8` under `session197/path-polynomials/`.
- `verify_extension.py`, SHA256 `3b537a9cbfea66aabe8b1d4cf59c84de45cc41d2ea6f80d053df8c083b1cad76`.
- `verify.py`, SHA256 `10c2a35ce3f506aefb031645e81cffc1a2fb4d73bf9bbb38a163ccbf8e35dfcb`.
