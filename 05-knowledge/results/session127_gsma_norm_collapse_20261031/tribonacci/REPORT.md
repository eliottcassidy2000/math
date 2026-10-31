# Tribonacci 305: exact recovery, norm mechanism, and scope

## Result first

**PROVED, with exact computational certificates.** Let `tau` be the positive
root of `f_3(X)=X^3-X^2-X-1`. The historical label **305** is the continued
fraction digit `a_5` (zero-based indexing):

\[
\tau=[1;1,5,4,2,305,1,8,2,1,4,6,14,3,1,13,5,1,7,23,\ldots].
\]

The mechanism is not a numerical coincidence:

\[
\tau^{-16}=56\tau-103,\qquad
N_{\mathbf Q(\tau)/\mathbf Q}(103-56\tau)=-1,
\]
\[
\tau-[1;1,5,4,2]=\tau-\frac{103}{56}
 =\frac{\tau^{-16}}{56}.
\]

The previous convergent is `46/25`; the complete quotient after `103/56` is

\[
t=\frac{1}{56^2(\tau-103/56)}-\frac{25}{56}
 =\frac{\tau^{16}-25}{56},\qquad 305<t<306.
\]

Thus the next convergent is `31461/17105`, not `31451/17105`.
The exact tail is about `305.8931204978` (decimal for orientation only).

**CITED + PROVED REDUCTION.** For every integer `m>0`,

\[
\tau^{-m}\in\mathbf Z+\mathbf Z\tau
 \quad\Longleftrightarrow\quad m\in\{3,16\}.
\]

This global cubic conclusion uses a published cubic Thue classification;
it is not inferred from a bounded recurrence search.

## 1. Recovery and a historical correction

Repository sources at commit `522535a7b6f519fd0c1e4304a2c5b1c909fce4c9`:

- `05-knowledge/results/rational_irrational_symbiosis.out`, Chapter 2;
  script `04-computation/rational_irrational_symbiosis.py`.
- `05-knowledge/results/free_exploration_s90g.out`, Connection 6;
  script `04-computation/free_exploration_s90g.py`.

Both contain the prefix through `305` and the convergent `103/56`. Both
generate a twenty-term continued fraction from the machine number
`1.8392867552141612`. The first wrong digit is **a_12: 17 should be 14**.
Everything from that position onward must be recomputed, not inherited.
The first twelve digits, including `305`, survive exact audit.
`sources.json` records source paths, byte counts, and SHA-256 digests.

The certified first forty digits (starting with `a_0`) are:

```text
1, 1, 5, 4, 2, 305, 1, 8, 2, 1, 4, 6, 14, 3, 1, 13, 5, 1, 7, 23,
1, 16, 4, 1, 1, 1, 1, 1, 2, 17, 1, 3, 1, 1, 1, 29, 1, 6, 1, 3
```

The known tournament references describe a different object with the same
characteristic polynomial. In particular, THM-4731,
`01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md`,
section 6, proves the tournament Hamiltonian count `h(0)=h(1)=h(2)=1` and
the blue multiplicity `h(floor(n/2)+1)` for `n>=4`. Those counts are not the
meaning of `305`, and their normalization must not replace the algebraic
root or standard Tribonacci sequence silently.

The exceptional unit identity is already in the literature: Adegoke,
Olatinwo, and Oyekanmi, *New Tribonacci recurrence relations and addition
formulas*, 2020, give the corresponding relation for every bilateral
Tribonacci recurrence. This packet independently verifies the polynomial
identity and makes **no priority claim**.
[Primary paper, pp. 165–167](https://nntdm.net/papers/nntdm-26/NNTDM-26-4-164-172.pdf).

## 2. Exact 305 certificate

Set

\[
F_3(p,q)=p^3-p^2q-pq^2-q^3.
\]

Integer evaluation gives

| `p/q` | `F_3(p,q)` |
|---|---:|
| `1/1` | -2 |
| `2/1` | 1 |
| `11/6` | -7 |
| `46/25` | 61 |
| `103/56` | -1 |
| `31461/17105` | 83326 |
| `31564/17161` | -10037 |

There is exactly one positive root of `f_3` by Descartes' rule. Its sign is
negative below that root and positive above it. Hence

\[
\frac{31564}{17161}<\tau<\frac{31461}{17105}.
\]

Once the exact prefix through `103/56` is established, the Möbius map
`t -> (103t+46)/(56t+25)` is decreasing. Its values at `305` and `306` are
exactly the two fractions displayed above, proving `305<t<306` without
any floating calculation.

The sparse identity is certified in either of two ways:

1. reduce `X^{-16}` in the integer quotient ring
   `Z[X]/(X^3-X^2-X-1)`, using `X^{-1}=X^2-X-1`;
2. reduce `(56X-103)X^16` and obtain the constant polynomial `1`.

Since `Norm(tau)=1`, taking norms gives the signed norm `-1` for
`103-56tau`. The sign is load-bearing: `56tau-103` is positive and has
norm `+1`.

## 3. Global classification of the cubic inverse-power collapse

**Cited input.** Knapp, Lemos, and Neumann, *Integral values of generating
functions of recurrence sequences*, Lemma 2.4 and Example 4.1, apply the
Delone–Nagell theorem to the discriminant `-44` cubic. In the present
coordinate order, the complete integer solution set of `F_3(P,Q)=1` is

\[
(P,Q)=(0,-1),(1,0),(2,1),(-103,-56).
\]

Indeed the negative-discriminant cubic bound allows at most four solutions
when the discriminant is `-44`, and these four are directly verified.
The solutions of `F_3(P,Q)=-1` are their negatives.
[Primary application, Lemma 2.4 and Example 4.1](https://math.loyola.edu/~mknapp/papers/General%20Quadratic%20-%202-10-23.pdf).

**Reduction.** If `tau^{-m}=a+b*tau`, its norm is `1`, so apply the above
classification to `(P,Q)=(a,-b)`. The four possible real values are
`tau`, `1`, `2-tau`, and `56tau-103`. For `m>0`, only the last two lie in
`(0,1)`. Their values are `tau^{-3}` and `tau^{-16}`. Strict monotonicity
of `m -> tau^{-m}` completes the iff proof.

This is a classification of the two-coordinate slice of the **single
unit-power orbit**, and in fact of all affine norm-one elements in this
cubic order. It does not classify arbitrary `S`-units, all Thue–Mahler
right sides, or any quartic/quintic slice.

## 4. Uniform next-quotient identity for k=3,4,5

Write

\[
f_k(X)=X^k-\sum_{j=0}^{k-1}X^j,
\quad F_k(p,q)=q^kf_k(p/q),
\]

and let `alpha_k` be the unique positive root. The polynomials are
irreducible: one may reduce modulo `3,2,5` for `k=3,4,5`, respectively.
The script exhausts every possible monic factor of degree at most
`floor(k/2)` over the corresponding finite field. Thus `F_k(p,q)` is the
field norm of `p-q*alpha_k`, not just a formally named quantity.

Suppose `p/q` is a convergent, `r/s` its predecessor, and `t` the next
complete quotient. Put `x=p/q`, `N=F_k(p,q)`, and

\[
D_k(x,\alpha_k)=\frac{f_k(x)-f_k(\alpha_k)}{x-\alpha_k}.
\]

The convergent determinant has absolute value one, giving

\[
|\alpha_k-p/q|=\frac1{q(qt+s)}.
\]

On the positive branch `D_k>0`. Substituting
`|alpha_k-x|=|N|/(q^k D_k)` proves the **exact identity**

\[
\boxed{t=\frac{q^{k-2}D_k(x,\alpha_k)}{|N|}-\frac{s}{q},
\qquad a_{n+1}=\lfloor t\rfloor.}
\]

This holds wherever the stated convergent and positive-secant conditions
hold; `k=3,4,5` is the audited common carrier, not an accidental cubic
formula. It is stronger than the approximation
`t ~ f_k'(alpha_k) q^{k-2}/|N|`, which must not be substituted as equality.

For a certified rational enclosure `L<alpha_k<U`, expand

\[
D_k(x,y)=\sum_{h=0}^{k-1}x^{k-1-h}y^h
 -\sum_{j=1}^{k-1}\sum_{h=0}^{j-1}x^{j-1-h}y^h.
\]

At `x>0`, a rigorous lower bound uses `L` in the positive sum and `U` in
the subtracted sum; the upper bound reverses them. Substitution in the
boxed identity supplies exact rational next-digit bounds. The script
pins every next digit of the first forty terms for each audited degree.

The common seed identity

\[
\alpha_k^{-k}=2-\alpha_k
\]

follows by multiplying `f_k(alpha_k)=0` by `alpha_k-1`. Its rational
shadow is `2/1`, with norm `1`, but it does not give an infinite family of
affine returns. The next digits after `2/1` are respectively `5,12,28`.
The cubic exponent `16` is the exceptional additional affine return.

## 5. Hostiles and nonclaims

- **A large digit does not imply norm one.** For `k=4`, `1091/566` is a
  convergent with norm `-179133` and next digit `21`. For `k=5`,
  `15941531077/8108825441` has norm
  `212788744741425601030399524296` and next digit `69`. The factor
  `q^{k-2}` and the secant are indispensable.
- A norm is not a Euclidean length or a unit-ball volume. No root-size,
  positivity, or metric fact may be recovered from the norm alone.
- The cubic global conclusion above relies on its cited classification.
  Quartic/quintic bounded searches are not global classifications.
- The unit case `|F_k|=1` is the `S=empty` corner of a possible
  Thue–Mahler program; no nonempty prime set is silently selected.
- The shared characteristic polynomial bridges tournament asymptotics
  and algebraic root arithmetic, but it does not identify their counts,
  initial conditions, or proof obligations.

## 6. Verification and reproduction

`check.py` uses only Python's standard library and active exception-based
checks, never optimizable `assert` statements. Two independent exact
continued-fraction methods are compared:

1. dyadic rational root isolation followed by interval inversion;
2. integer evaluations of `F_k` at Möbius images of integer tail bounds.

Norms are separately computed by an explicit permutation determinant of
the multiplication matrix. The exceptional unit identity has independent
positive- and negative-power checks. Historical digit `17` is a hostile
control. The externally cited completeness theorem is explicitly marked
as cited, not represented as a machine-verified result.

```sh
python3 check.py > receipt.json
python3 -O check.py > /tmp/tribonacci-receipt-optimized.json
cmp receipt.json /tmp/tribonacci-receipt-optimized.json
```

Both modes and an isolated-directory replay pass with byte-identical
output; `RUNS.json` records their hashes. Only one canonical `receipt.json`
is retained. No repository files were modified by this lane.
