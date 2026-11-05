# Stieltjes positivity for repeated physical wrapper words

Session 197. **PROVED all-order Hankel total nonnegativity for a physical construction monoid.** Fix any finite word W in the following two operations:

- P: take the ordered join of the current tournament with Tr_2.
- D: replace the large vertex of a directed triangle by the current tournament, leaving the other two triangle vertices singletons.

Apply the letters in their written order, and repeat exactly this same word m times. Start either at K1 or at Tr_2. Let mu_m be the resulting Hamilton-path count, including mu_0=1. Then every finite minor of the infinite Hankel matrix

    H=(mu_(i+j))_(i,j>=0)

is nonnegative. In particular mu is a Stieltjes moment sequence. The empty word and padding-only words are allowed and give the constant sequence 1; strict positivity is therefore not asserted for the entire class.

This extends the single-triangle-wrapper vacuum theorem. It does not permit an arbitrary extra initial tournament or an arbitrary unabsorbed prefix before the repeated word.

## 1. Exact physical profile matrices

For the labelled-bin polynomial f=B_T, the two physical operations act as

    P f(x)=x^2 f(x),
    D f(x)=x^2 f(x)+x[f(x+1)-f(x-1)].                    (1)

These are the previously proved ordered-join identity and triangle-wrapper identity. Define

    q_k(x)=x product_(j=1)^k (x^2-j^2),
    r_k(x)=x^2 product_(j=1)^k (x^2-j^2).

Use the q basis for odd input order, the r basis for even input order. Both operations preserve parity. The initial polynomials are respectively q_0 and r_0, and evaluation at x=1 picks exactly the coefficient of basis vector 0.

In either basis, the column matrix P is lower bidiagonal:

    P_(k,k)=(k+1)^2, P_(k+1,k)=1.

The column matrix J for D is tridiagonal:

    J_(k,k)=a_k, J_(k+1,k)=1, J_(k-1,k)=b_k,

where

    odd:  a_k=k^2+6k+3, b_k=2k^2(2k+1),
    even: a_k=k^2+6k+5, b_k=2k(k+1)(2k+1).

The b_0 term is absent. The scalar count after repeating W is therefore

    mu_m=(T^m)_(0,0),                                  (2)

where T is the corresponding product of P and J matrices. The product order is the reverse of the written operation word because these are column-action matrices. Every finite word has finite bandwidth; all matrix products and the entries used below are finite sums.

## 2. Why the production matrices are totally nonnegative

A matrix is totally nonnegative (TN) when every finite minor is nonnegative. This assertion is stronger than entrywise positivity.

The reviewed proof `STIELTJES_VACUA.md` establishes the strict positivity, in both parity cases, of the pivots

    ell_0=a_0,
    ell_k=a_k-b_k/ell_(k-1), k>=1.

Consequently J has the exact infinite bidiagonal factorization

    J=L diag(ell_0,ell_1,...) U,
    L_(k,k)=U_(k,k)=1,
    L_(k,k-1)=1/ell_(k-1),
    U_(k-1,k)=b_k/ell_(k-1).

All displayed entries are nonnegative. A lower or upper bidiagonal matrix with nonnegative entries is TN, as is a nonnegative diagonal matrix: each minor is zero or a sum of nonnegative order-preserving path products. Cauchy-Binet therefore proves J is TN. The same elementary bidiagonal argument proves P is TN. Every finite product T in these matrices is TN, as is its transpose.

The infinite-matrix statements are unambiguous. For any specified minor, choose a finite prefix large enough to include every intermediate index in the banded product. The usual finite Cauchy-Binet formula then applies exactly. No claim of arbitrary algebra preserving positivity is used.

## 3. From a TN matrix to all Hankel minors

For a row-finite matrix A define its output array

    O(A)_(m,k)=(A^m)_(0,k), m,k>=0.

**Lemma.** If A is row-finite and TN, then O(A) is TN.

Take a minor with distinct increasing row indices and column indices. Induct lexicographically on the sum of its row indices and its size. If row0 is selected, that row is the unit row: the minor is zero if column0 is absent, and otherwise expands with positive sign to the smaller minor obtained by deleting row0 and column0. If row0 is absent, subtract1 from every selected row index and use O(A)_(m+1,*)=O(A)_(m,*) A. Cauchy-Binet writes the minor as a finite sum of minors of the earlier output rows times minors of A. Every summand is nonnegative by induction. The initial empty minor equals1. Row-finiteness makes the sum finite.

Our T and T^T are both row-finite and TN, so both output arrays are TN. Matrix multiplication gives the exact factorization

    H=O(T) O(T^T)^T,

because its (i,j) entry is sum_k (T^i)_(0,k)(T^j)_(k,0)=(T^(i+j))_(0,0). Each required minor uses only finitely many states. Cauchy-Binet now proves that every Hankel minor is nonnegative, establishing the theorem.

In particular every finite principal minor of both (mu_(i+j)) and its shifted matrix (mu_(i+j+1)) is nonnegative; the shifted matrix is itself a submatrix of H. By the classical Stieltjes moment criterion this is equivalent to a positive measure on [0,infinity) with moments mu_m. Equivalently, the standard Stieltjes continued-fraction correspondence supplies nonnegative S-fraction coefficients, with termination allowed in degenerate cases. The combinatorial theorem proved here is the all-minor Hankel assertion; no explicit universal formula for those continued-fraction coefficients is claimed.

## 4. Single-vertex padding between wrappers

Let X be ordered join with a single vertex, so Xf(x)=xf(x). Its two cross-parity matrices are particularly simple:

    X q_k=r_k,
    X r_k=q_(k+1)+(k+1)^2 q_k.                         (3)

Thus the odd-to-even matrix is the identity, and the even-to-odd matrix is exactly the positive lower bidiagonal matrix P from Section1. Both are TN. A triangle wrapper preserves parity and uses the corresponding odd or even J matrix at that stage.

Consequently the theorem extends to any fixed word in D and X containing an even number of X letters. Such a word returns to the starting parity. Its column-action matrix is a finite product of the appropriate TN matrices, so Sections2-3 apply without modification. The operation is physical at every intermediate stage, even when the parity changes between wrappers. The initial tournament remains K1 or Tr_2.

For a word containing an odd number of X letters, this argument proves the assertion for the even subsequence obtained by repeating W twice per period. No assertion about all moments of an odd-X word is needed here: its one-period matrices act between different parity spaces, so it cannot be treated as one of the fixed matrices T above without a further argument.

## 5. A physical transitive-mirror application

Take a transitive side S=Tr_t with marked vertex having j predecessors, 0<=j<t. Form the rank-one mirror context: a copy O of S, a converse copy I, a substituted input T between I and O, and the sole reversed cross edge O_j->I_j.

The vertices with side labels greater than j form an initial transitive block in I and a final transitive block in O, each of size t-1-j. The remaining central piece is the directed triangle whose large module is

    Tr_j -> T -> Tr_j,

and whose other two vertices are O_j and I_j. Hence its full profile operator is exactly

    P^(t-1-j) D P^j.                                  (4)

This follows from a physical block decomposition and the labelled-bin product identity, not a merely formal polynomial factorization. Repeating any fixed composition of these mirror contexts and even transitive paddings therefore falls under the theorem. The starting input must still be K1 or Tr_2, as stated.

## 6. Examples and essential scope boundaries

Starting at K1, the first five moments are

    word D:  1,3,15,123,1479;
    word PD: 1,9,345,35241,7808889;
    word DP: 1,3,51,3171,483603.

Here PD means pad first and then triangle-wrap; DP means the opposite order. Their second Hankel determinants are respectively6,264,42.

An arbitrary initial prefix does not inherit the assertion. Starting at Tr_4 and repeating only D gives moments1,17,241, whose second Hankel determinant is -48. By contrast, starting at Tr_2 and repeatedly applying the entire word PD gives1,17,913 and determinant624. The difference is that padding is repeated as part of the period in the second case. One cannot silently move a one-time prefix into the period.

The theorem gives no real-rootedness claim for the input bin polynomial, no Hankel assertion for arbitrary positive initial profiles, and no guarantee of prime or7p counts in these sequences. It is an exact positivity bridge for the specified physical periodic constructions.

## 7. Reproduction and dependencies

Run `python3 verify_periodic_stieltjes_words.py` and `python3 -O verify_periodic_stieltjes_words.py`. The checker computes moments directly through the exact profile operators for every P/D word through length4 in both parity vacua, and every D/X word through length4 with an even number of X letters. It checks finite production-matrix minors, Hankel minors and initial rational S-fraction coefficients. It also builds actual tournaments for short words and directly counts their Hamilton paths, checks the transitive-mirror block identity on physical inputs, and reproduces the one-time-prefix hostile. These finite controls supplement Sections2-4 rather than proving the unbounded result by extrapolation.

Named frozen dependencies:

- `PROOF.md`, SHA256 `b88c86ac0afa64d3a1bf399c56fa7f8aad046fb93424c70a9ace8d164708183b`.
- `STIELTJES_VACUA.md`, SHA256 `dbe7c0e68e77e5e14e389c57221578191a22381b086cb4c5864c7b2d1abc6fe9`.
- `verify.py`, SHA256 `10c2a35ce3f506aefb031645e81cffc1a2fb4d73bf9bbb38a163ccbf8e35dfcb`.

The standard TN output-array and Stieltjes correspondences are used explicitly above; no literature-priority claim is made.
