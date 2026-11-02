# Independent audit: parity-coset tensor bound and matching uncertainty

Status: **accepted within the matching-only domain stated below**. This is a
proof audit plus separate finite exact controls, not an inference from the
controls and not a novelty claim. No conclusion is asserted for arbitrary
uncertainty graphs or the outstanding finite-order stability gap.

## 1. Tensor lemma: induction is valid

For `r >= 0`, let every coefficient `b_U` be integral, with `b_empty` odd
and every other `b_U` even. For `|x_i| <= 1`, define

    f_S = sum_{J subset S} b_(S\J) product_{j in J} x_j.

Then

    sum_S f_S^2 >= product_i (1+x_i^2).

Split the output vector into `(a, b+x_r a)`, with `a` the transform of
coefficients not containing `r` and `b` the transform of coefficients
containing `r`. The coefficient vectors defining `a`, `a+b`, and `a-b`
all have odd constant and even other coefficients. Consequently, by the
induction hypothesis, with `L=product_{i<r}(1+x_i^2)`,

    A=||a||^2 >= L,
    A+C+2D >= L,   A+C-2D >= L,

where `C=||b||^2 >= 0`, `D=<a,b>`. Hence `2|D| <= A+C-L`.
Writing `t=|x_r| <= 1`, the required norm is at least

    (1-t+t^2) A + (1-t) C + t L >= (1+t^2)L.

Both discarded coefficients are nonnegative. The base case is the square
of an odd integer. The proof does not require `b` itself to be in the odd
coset, nor does it require any sign restrictions on the coefficients.

The domain and parity assumptions are essential: the mutation controls
below give exact counterexamples when any one is dropped.

## 2. Pfaffian grouping: no missing sign or parity assumption

Let a nonempty weighted tournament have `n=2r+t` vertices, with all
nonbinary edges contained in `r` disjoint designated pairs. After a
simultaneous vertex permutation, put each pair in adjacent positions,
write its upper skew entry as `x_i`, and append the binary observer last
when needed. Here `q_i=(1-x_i^2)/4`.

Partition principal subsets according to:

1. which paired blocks contribute exactly one endpoint, and which endpoint;
2. which of the `t` unmatched vertices are selected;
3. for every other paired block, allow either zero or two endpoints.

For the even principal-Pfaffian sum, retain the groups whose fixed base
has even size. For the observer-containing sum, retain groups whose base
of old vertices has odd size, and add the observer to that base. Thus the
fixed base `R` is always even, and every edge inside it is binary.

In each such group, let `F` be the set of free pairs. For `S subset F`,
expand the Pfaffian on `R` plus the pairs of `S` in their internal skew
variables. Full paired endpoints are adjacent even after taking a
principal restriction; the derivative sign is

    (-1)^(p+(p+1)+1) = +1.

Repeated pair differentiation therefore gives exactly the tensor form
in Section 1, not a version with uncontrolled signs. Its coefficient
`b_U` is the Pfaffian on `R` plus pairs `U`, with the `|U|` internal paired
entries set to zero. `b_empty=Pf(R)` is odd (including `Pf(empty)=1`).
For nonempty `U`, reduction modulo 2 turns every other edge into 1.
Inclusion-exclusion over the missing disjoint pair edges gives

    b_U ≡ sum_{J subset U} (-1)^|J| (|R|+2|U|-2|J|-1)!! ≡ 0 (mod 2),

because each double factorial is odd, with `(-1)!!=1`, and there are
`2^|U|` summands. Therefore every group satisfies the tensor lemma.

This is a whole-group estimate. An individual full Pfaffian can vanish
while the corresponding uncertain-edge monomial is nonzero; one exact
six-vertex witness has polynomial `3x-2` at `x=2/3`.

## 3. Signed margins and normalizations

Write `u=det(I+A)` and `v=(-1)^n det(A)`. The inherited skew identities are

    u+v = 2^(1-n) sum_(even old subsets) Pf^2,
    u-v = 2^(1-n) sum_(odd old subsets) Pf_(subset+observer)^2.

The weighted count of group lower bounds has parity polynomial

    Q(z) = product_i [(1+x_i^2)+2z] (1+z)^t.

Its even and odd sums are `(Q(1)+Q(-1))/2` and
`(Q(1)-Q(-1))/2`, respectively. Multiplying by `2^(1-n)` yields:

- If `t>=1`, both `u+v` and `u-v` are at least `product_i(1-q_i)`.
- If `t=0`, then `r>=1` and

      u+v >= product_i(1-q_i) + (-1)^r product_i q_i,
      u-v >= product_i(1-q_i) - (-1)^r product_i q_i.

Thus

    u-|v| >= product_i(1-q_i) - 1_(t=0) product_i q_i.

For `r=0`, the nonempty-domain assumption implies `t>=1`, and the bound
is the binary floor 1. The empty tournament must not silently be placed
in that binary case: its margin is zero. For `n=2,r=1`, the no-isolate
formula is exactly `1-2q`.

Ordinal sums of weighted pairs and binary singletons attain these
absolute floors; in the no-isolate case they attain both signed floors
with the displayed parity signs. No general equality classification is
asserted.

## 4. Rounding reserve and edge cases

Use the accepted one-edge inequality

    Delta(A) >= E_edge Delta(A_endpoint) + q_e (u_R-|v_R|).

Round the matching edges in any fixed order `1,...,r`. At step `i`, the
pair-deleted remainder has `r-i` still uncertain pairs and exactly
`t+2(i-1)` binary isolated vertices in its uncertainty graph. All cross
edges remain binary, regardless of previous rounding outcomes.

If `t>=1`, every correction is at least

    q_i product_{j>i}(1-q_j).

If `t=0`, only the first step lacks a binary isolate, and it subtracts
`q_1 product_{j>1} q_j = product_j q_j`. Subsequent rounded pairs supply
isolates. Therefore telescoping gives

    Delta(A) >= E Delta(T)
                + 1-product_i(1-q_i) - 1_(n=2r) product_i q_i.

For the rounding proof using nonempty remainders, take `n>=3`. Order two
is handled directly: `Delta(A)=Delta(T)=0`, and the displayed reserve is
`q-q=0`. Order one or `r=0` reduces to the tautology `Delta(A)=E Delta(T)`.
The empty tournament is outside this statement.

The binary expectation is over the same independent product-law
rounding for every fixed order. It is retained once, not multiplied by
the number of edges or orders. Determinant signs may change between
branches; the absolute-margin lower bound already handles that.

## 5. Exact controls and adversarial controls

`verify_exact.py` is standalone, uses Python `Fraction` throughout,
implements Gaussian determinants and independent recursive Pfaffians,
and raises explicit errors rather than relying on optimization-sensitive
assertions. Its completed output is `exact_results.json`.

- 12,604 exhaustive tensor controls for dimensions 0, 1, 2;
- 300 additional rational tensor controls for dimensions 3, 4, 5;
- 876 exact rational tournament controls for orders 1 through 9, including
  every binary cross orientation at order four and 16 nonbinary pair
  probability combinations for each orientation;
- both signed Pfaffian identities and signed/absolute floors verified in
  every tournament control;
- nonlinear baseline-containing reserve checked by exact independent
  enumeration of matching-edge binary completions; using `E H(T)=H(A)`,
  this evaluates `E[u(T)+|v(T)|]-(u(A)+|v(A)|)` directly;
- 7,394 complete grouped Pfaffian transforms, coefficient parity checks,
  and group lower bounds;
- 20 ordinal-block exact floor equalities.

`verify_mutations.py` rejects eight false strengthenings exactly:

1. dropping the coordinate restriction `|x|<=1`;
2. dropping evenness of the nonconstant coefficients;
3. dropping oddness of the constant coefficient;
4. omitting the no-isolate product-q correction;
5. reversing the signed parity correction in order two;
6. replacing the group estimate with an individual-Pfaffian estimate;
7. replacing an absolute matching floor by a conditional `(1-q)` endpoint
   interpolation factor (pair plus two directed three-cycles);
8. permitting unaccounted nonbinary cross edges.

Both scripts are also run with `python -O`; their JSON outputs must match
byte for byte. Controls are finite verification evidence only; Sections
1 through 4 supply the all-order proof audit.
