# Independent audit of the complete degree-three counting series

October 22, 2026. **PASS.** The all-six closure supplies exactly

$$R_6(z)=z^5+z^6+z^7.$$

The coefficients count ordinary isomorphism types of **strong self-converse
deficit-six centers**. The three inherited small centers are retained in the
older exact receipt; [the incidence proof](README.md), the
[order-15 closure](../nontransitive_fifteen_20261022/README.md), and the
[VT21 companion](../nontransitive_fifteen_20261022/VT21.md) remove all six
remaining orders. This audit does not repeat those classifications.

## Exact object counted

We use precisely the inherited formula in
[the three-presentation classification, §5](../metagraph_classification_20261021a/README.md#5-the-next-trivalent-layer-has-one-remaining-finite-center-polynomial).
Its graph $Q$ is the **fundamental converse-merged paired-black graph**.
A non-self-converse pair has degree $h$; a self-converse class has degree
$(h-b)/2$. This is not the literal version with doubled paired-black edges.

Set $A(z)=(1-z-z^3)^{-1}$. The number $g_n$ of degree-three vertices of
$Q$ at order $n$ has the inherited generating function

$$G(z)=\frac{z^5}{2}\bigl(A(z)^2-A(z^2)\bigr)
       +A(z^2)R_6(z)+z^{10}A(z^2)^2(1+z+z^3).$$

The three disjoint contributions count, respectively, non-self-converse
$h=3$ converse pairs; rigid palindromic padding around a strong deficit-six
center; and nonrigid $h=3$ halves with middle $\varnothing,K_1,C_3$.
Self-converse pure-blue $h=3$ classes have degree zero, and are not included
in the first term. These are neither the tournament counts A000568 nor
labelled tournament counts. The formula also does not reconstruct neighbours,
port matching, edge multiplicities or black connectivity.

## Independent algebra and recurrence

Put $p=1-z-z^3$, $u=p(z^2)=1-z^2-z^6$ and $D=p^2u^2$. Direct exact
polynomial arithmetic over $\mathbb Q$, independent of the parent's SymPy
calculation, gives

$$
\begin{aligned}
N(z)={}&-z^{18}+z^{17}-z^{16}+2z^{15}-z^{14}
        +3z^{10}+z^9-2z^8-2z^7+z^5,\\
G(z)={}&N(z)/D(z).
\end{aligned}
$$

Indeed the numerator before simplification is

$$N=\tfrac12z^5(u^2-p^2u)+p^2u(z^5+z^6+z^7)
                +z^{10}p^2(1+z+z^3).$$

An exact extended Euclidean computation gives $\gcd(N,D)=1$, and the
receipt stores the rational-coefficient Bezout multipliers proving it.
Thus the minimal eventual constant-coefficient homogeneous recurrence has
order **18**. In ascending coefficient order,

$$
(D_0,\ldots,D_{18})=
(1,-2,-1,2,1,2,-4,2,0,0,-1,-4,3,-2,3,-2,2,0,1).
$$

With $g_n=0$ for negative $n$, the exact inhomogeneous identity is

$$\sum_{j=0}^{18}D_jg_{n-j}=[z^n]N(z).$$

Consequently the homogeneous recurrence holds for every **$n\ge19$**.
At $n=18$ its left side is $-1$, not zero. This is the exact onset for the
normalized minimal order-18 relation, not merely the earlier safe onset
deduced from an upper bound on the unknown center degree.

The requested indexed prefix is

$$
(g_{15},g_{16},g_{17},g_{18},g_{19},g_{20},g_{21})
=(114,175,276,425,665,1019,1583).
$$

## Reproduction

[gf_check.py](gf_check.py) performs two independent coefficient calculations:
exact rational division by $D$, and direct convolution of the three
combinatorial series using $a_n=a_{n-1}+a_{n-3}$. They agree through degree
256. The all-order statement follows from the exact polynomial identity;
the finite prefix comparison is a control, not the proof of rationality.
The script also checks the Bezout identity, the failure at 18, the required
prefix, the three inherited small-center deficits and the fixed source
receipt hashes. It does not invoke SymPy or enumerate tournaments.

```text
python -B 05-knowledge/results/center_incidence_20261022/gf_check.py --write
python -B -O 05-knowledge/results/center_incidence_20261022/gf_check.py --verify
```

Normal and optimized saved-output comparison both pass **1,017 exact checks**.
LF [gf_receipt.json](gf_receipt.json) SHA256:
`8627a0abafb835fecd1c89ddc4804eb73952a8b6020720e97183dacbcc259257`.
This is an independent audit of the parent's final consequence; no shared
formula source or navigation file was edited.
