# A positive split-state repair for two prescribed endpoints

Session 197 continuation of `MARKED_CONTEXTS.md`. **Exact positive two-state closure for one linked endpoint pair under whole-module contexts.** Nonuniform port births and arbitrary multiport connectivity interfaces remain outside this statement. No literature-priority claim is made.

## 1. Two different connectivity states are necessary

Fix distinct named vertices u,v in T. Define

- j_T(x): covers with one reserved path starting at u and ending at v, and x ordinary labelled bins;
- k_T(x): covers with two distinct reserved paths, one starting at u and the other ending at v, and x ordinary labelled bins.

Thus j_T(x)=J_T(x+1) and k_T(x)=K_T(x+2); their reservation shifts differ. Every reserved path is nonempty. The definitions distinguish one connected path from two components, rather than recording only the endpoint degrees.

The j-state alone does not close even under W(T)=C3[T,K1,K1]. In the transitive input v->u, j=0 and k=1. In v->w->u, j=0 and k=x. The two j-polynomials are identical, yet their W outputs have j=1 and j=x respectively. In particular the prescribed Hamilton endpoint count at x=0 is one in the first output and zero in the second. Both outputs are actual strong tournaments. A connectivity-sensitive split coordinate is genuinely required.

## 2. Universal triangular closure

For a fixed rooted tournament (Q,r), substitute the entire input T at r and retain u,v. Let A_Q be the unmarked bin-profile operator from `UNIVERSAL_CONTEXTS.md`, and C_Q the one-start/one-end reserved-bin operator from `MARKED_CONTEXTS.md`.

**Theorem 1.** There is an operator B_Q such that the entire pair of shifted profiles obeys

    j_(Q[r:T]) = A_Q j_T + B_Q k_T,
    k_(Q[r:T]) = C_Q k_T.                                    (1)

All three operator matrices are entrywise nonnegative integers in the extended central-factorial basis {1,q_0,r_0,q_1,r_1,...}.

**Proof of the decomposition.** Cut the output's reserved u-to-v path into maximal input-module runs. If u and v lie in the same run, that run must be the entire reserved path: nothing may precede u or follow v. Choosing that internal path and then covering the remaining context gives A_Q j_T by unmarked full-profile linearity.

Otherwise, its first module run starts at u and its last module run ends at v. These are distinct marked fragments; all other module fragments are unmarked. Cutting the entire output cover therefore gives exactly the local k-state. The recombination kernel depends only on the number of unmarked fragments and Q, not their internal vertices, since all exterior incidences are uniform. This defines B_Q linearly on the falling-factorial coordinates of k. The two cases are disjoint and exhaustive. The second equation is precisely the already proved marked-bin closure. QED.

### Positivity of the connecting block

Use the actual transitive test input

    v -> Tr_m -> u.

It has j=0 and k=x^m. Thus B_Q(x^m) is exactly the output j-polynomial on this physical family. These monomials span the whole polynomial space.

In the output, choose the reserved u-to-v path, then use the odd-cycle bin formula on its induced complement. Delete all free-chain runs. The resulting distinguished open-path skeleton starts at u, ends at v, and necessarily traverses exterior vertices. Both anchor fragments are singletons in the test input. Every intervening free-chain run contributes a nonempty-bin factor nabla=S-I. Crossing odd-cycle skeletons contribute the same L_(t,epsilon) kernels as before. A skeleton with t_0 free-chain run positions in the open path, c_0 pure exterior cycles, c_1 crossing cycles and a unused exterior vertices therefore contributes

    2^c_0 X^(a+c_0+c_1) nabla^t_0 product_i L_(t_i,epsilon_i).   (2)

Every factor is positive on the extended central-factorial cone. Canonical cycle rotations, the fixed endpoints of the open path, and distinct run positions give the same cutting/reconstruction bijection as in the marked proof. The skeleton identity holds on every monomial and hence, by linearity, on all polynomials. Its positive factors then prove positivity on every formal basis vector, not only physical k-states. Together with the prior A_Q and C_Q theorems, this proves Theorem 1.

The boundary u=v is separate: a path starting and ending at the same vertex is its singleton path. It cannot split into two fragments without repeating a physical vertex. Its shifted polynomial is B_(T-u), and its transport is simply A_Q B_(T-u). No off-diagonal split term is introduced.

## 3. An exact physical 2-by-2 transfer

For W(T)=C3[T,K1,K1], with a->T->b->a,

    j_new = D j + k,
    k_new = (D+2X+2S+I)k.                                    (3)

The first diagonal and the second equation are already established. In the split case the only possible connection is

    start-fragment, b, a, end-fragment.

It consumes both added vertices and leaves no intermediate input-module run, so B_W=I. This proves (3) as a full-profile identity.

At x=0, the prescribed u-to-v Hamilton count in W(T) is exactly k_T(0), the original spanning mixed two-path-cover count. Formula (3) is a usable positive coupled state; replacing it by j alone loses necessary connectivity information.

## 4. The connecting block consumes one unit of matching capacity

Let nu be the root cross-matching number from the companion matching theorem. If nu=0, B_Q=0. If nu>=1, B_Q has exact largest forward shift nu-1, with a nonzero coefficient polynomial having nonnegative integer coefficients.

Every open-path skeleton from u to v has one more exterior O-to-I gap than its number t_0 of free-chain runs. Every crossing cycle has as many such gaps as runs. Selecting one O-to-I transition from each gap gives a matching, so

    1+t_0+sum_i t_i <= nu.

The total forward shift in (2) is at most t_0+sum_i t_i, proving the upper bound and the nu=0 case. Conversely order the edges of a maximum matching into one open path from u to v, with a free-chain run between consecutive matched pairs. This gives nu-1 run positions and a positive highest shift. All highest shift contributions in (2) are nonnegative, so no cancellation removes it.

In the balanced case N=2nu, with mu perfect cross matchings, the highest forward coefficient is exactly

    [S^(nu-1)] B_Q = mu nu (x+1)(x+2)...(x+nu-1).              (4)

Indeed every exterior gap must be one matched O-to-I pair. Choose j>=1 matched pairs, order them in the distinguished open path, and arrange the remaining pairs into permutation cycles. The sum for one matching is

    sum_(j=1)^nu binom(nu,j) j! x^(overline(nu-j))
      = nu (x+1)^(overline(nu-1)).

The rising-factorial identity follows from the usual binomial hockey-stick identity, or by inserting the next labelled pair into the distinguished list or a permutation cycle. This is the open-path counterpart of the balanced unmarked coefficient mu x^(overline nu).

These exact matching formulas still leave lower coefficients and physical realizability constraints to be controlled. They do not imply a universal digit-return construction.

## 5. Verification

Run `python3 verify_split_endpoints.py` and the optimized version. Independent rooted subset-path DP computes j directly, while the previously audited marked subset-cover routine computes k. The checker compares complete actual context profiles against the triangular decomposition, tests positive connecting-block images on formal basis vectors, verifies (3) on actual rooted inputs, and checks the j-only hostile. Selected higher-matching templates are checked against (4) by extracting shift coefficients from actual monomial columns using exact polynomial interpolation, not by reading them off the skeleton formula.
