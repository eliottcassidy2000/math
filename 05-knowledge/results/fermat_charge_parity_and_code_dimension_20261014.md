# Charge parity fails; the Golay code mechanism is excluded uniformly at Fermat primes

Date: 2026-10-14. Status: exact finite parity counterexamples and a proved uniform code-dimension lemma. This is a bounded successor to the [Golay obstruction](path_charge_golay_obstruction_20261014.md) and [all-class p=17 charge-coverage certificate](fermat_p17_path_charge_surjectivity_20261014.md).

## 1. The stronger parity carrier fails

For a cyclic tournament T_S, write C(P)=sum i*zeta^(v_i) in R=Z[zeta]/(eta), where eta=1+sum_(s in S) zeta^s and P=(0,v_1,...,v_(p-1)) is a normalized Hamiltonian path. Surjectivity is proved by the predecessor certificates for every prime-cyclic tournament at 3,5,17. A stronger candidate would require every charge fibre to have odd size, connecting the carrier more directly to the corpus's parity mechanisms. That stronger statement is false.

| p | S | Exact quotient | Odd fibres | Even fibres | Total normalized paths |
|---:|---|---|---:|---:|---:|
| 17 | {2,3,4,5,6,7,8,16} | F_103, zeta=13 | 51 | 52 | 793,904,749 |
| 11 | {1,6,7,8,9} | F_23, zeta=3 | 11 | 12 | 8,401 |

The quotient identifications are proved and certified in the predecessor packet. At p=17, charge 1 is an even fibre, although it contains the explicit path

\[
 (0,16,2,4,3,7,6,11,1,5,8,10,9,13,12,15,14).
\]

Thus this is a nonempty even fibre, rather than a failure of charge coverage. At p=11, direct enumeration gives exactly 344 paths of charge 1. All 23 fibres in that example are nonempty. Complete parity supports and the p=11 integer histogram appear in the result JSON.

### Exact subset/end recurrence

Let q be the displayed quotient prime, r the image of zeta, and M a subset of the nonzero vertices. In F_2[X]/(X^q-1), let F_(M,v) encode the charges of directed paths starting at 0, using exactly M after 0, and ending at v. Then

\[
 F_{\{v\},v}=\begin{cases}X^{r^v}&0\to v,\\0&\text{otherwise},\end{cases}
\]

and, for |M|>1,

\[
 F_{M,v}=X^{|M|r^v}
 \sum_{\substack{u\in M\setminus\{v\}\\u\to v}}
 F_{M\setminus\{v\},u}. \tag{1}
\]

The final fibre-parity polynomial is the sum over endpoints for M=F_p*. Representing its q coefficients by a q-bit integer makes addition XOR and multiplication by a monomial a cyclic bit rotation. This computes exact parity without enumerating all Hamiltonian paths.

The verifier simultaneously maintains ordinary integer path counts and checks parity agreement in every state. It repeats the bitset calculation under the relabelling v -> -v, with S -> -S and r -> r^-1, obtaining the identical charge polynomial. For the p=11 case it additionally enumerates all 8,401 paths directly (52,167 visited prefixes), checking the entire integer histogram against the DP parity. At p=17 it uses subset DP, not a full path census.

## 2. Half-sized support is not a translated-square formula

The odd-fibre supports have sizes (103-1)/2 and (23-1)/2. This numerical coincidence does not yield an affine quadratic-residue description.

For each field, the verifier checks both square classes of scalars and every translate. None of the 206 candidates at q=103 or 46 candidates at q=23 equals the odd-charge support. More strongly, it tests every nonidentity invertible affine map x -> ax+b: 10,505 maps at q=103 and 505 at q=23. Neither support has a nonidentity affine stabilizer.

Consequently neither support can be an affine image of a union of cosets of any **nontrivial** multiplicative subgroup. Such a union is preserved by that subgroup, and an affine image would retain its conjugated nontrivial affine symmetry. This rules out that family of simple cyclotomic support formulas; it does not classify every possible arithmetic description of the parity polynomial.

## 3. Uniform Fermat code-dimension lemma

The Golay counterexample at p=11 arises from a small-characteristic cyclic code whose one-coordinate extension is self-orthogonal. The following lemma rules out that precise global mechanism at every Fermat prime, independently of a finite census.

Let p be an odd prime, m=(p-1)/2, and retain eta as above. Let ell be a norm prime, and put

\[
 f_\ell=\gcd(\eta(x),\Phi_p(x))\in\mathbb F_\ell[x],
 \qquad d_\ell=\deg f_\ell.
\]

Here eta(x)=1+sum_(s in S)x^s. Since ell differs from p, Phi_p modulo ell is squarefree. Factoring the integral ideal (eta) over ell gives

\[
 \ell^{d_\ell}\mid N(\eta). \tag{2}
\]

Indeed the degree of each selected residue factor is its prime ideal's residue degree, and its valuation in eta is at least one. Multiplicities can increase the norm exponent but cannot reduce it.

**Norm bound.** For every prime cyclic tournament,

\[
 N(\eta)\le\left(\frac{p+1}{4}\right)^m. \tag{3}
\]

To see this, the real skew matrix K has nonconstant eigenvalues +/- i*mu_j, j=1,...,m, and a constant eigenvector with eigenvalue 0. Since every off-diagonal entry has square 1, its squared Frobenius norm gives sum mu_j^2=pm. The norm is product_j (1+mu_j^2)/4. Applying the arithmetic-geometric mean inequality proves (3).

Combining (2)-(3), whenever ell>(p+1)/4 one has d_ell<m. In particular this holds at every Fermat prime: the predecessor's norm-prime theorem gives ell>=2p+1. Hence

\[
 d_\ell\le m-1,\qquad
 \dim_{\mathbb F_\ell}\mathcal C_\ell=p-d_\ell\ge m+2,
 \quad
 \mathcal C_\ell=(f_\ell)\subset\mathbb F_\ell[x]/(x^p-1). \tag{4}
\]

**Conclusion.** At a Fermat prime, no linear one-coordinate extension of this length-p cyclic code can be self-orthogonal for the standard dot product in length p+1.

**Proof.** Any such extension c -> (c,L(c)) retains dimension at least m+2. A self-orthogonal subspace for a nondegenerate bilinear form in dimension p+1=2m+2 has dimension at most m+1. This contradicts (4). All relevant residue characteristics here are odd. QED.

This excludes the global identity

\[
 \sum_v c_v^2+\left(\sum_v c_v\right)^2=0
 \quad\text{for every }c\in\mathcal C_\ell,
\]

which produced the ternary Golay charge exclusions. For this particular identity, there is also a short diagnostic: Phi_p belongs to C_ell, so the all-ones word would require p(p+1)=0 modulo ell. This is impossible when ell>=2p+1. The dimension argument is stronger, excluding every linear one-coordinate extension, not only the coefficient-sum extension.

The lemma concerns a global identity on the full cyclic code. It does not rule out other composition restrictions, identities on smaller subsets, or failures of Hamiltonian admissibility. It establishes no new uniform surjectivity theorem.

## 4. Reproduction and remaining target

The standard-library script `fermat_charge_parity_check.py` emits the frozen `fermat_charge_parity_check.json`. Both ordinary and `python3 -O` runs produce the same output, with 1,058,833 active checks and approximately 1.4 seconds per run in the staging environment.

```sh
python3 04-computation/experiments/fermat_charge_parity_check.py > /tmp/fermat_charge_parity_check.json
cmp 05-knowledge/results/fermat_charge_parity_check.json /tmp/fermat_charge_parity_check.json
python3 -O 04-computation/experiments/fermat_charge_parity_check.py > /tmp/fermat_charge_parity_check_O.json
cmp 05-knowledge/results/fermat_charge_parity_check.json /tmp/fermat_charge_parity_check_O.json
```

The live question remains a uniform composition-and-order argument for charge coverage at Fermat primes. Universal oddness of the fibres and affine translated-square support are both disproved routes, while the Golay self-orthogonality mechanism is uniformly excluded. These boundaries do not replace the missing coverage argument.
