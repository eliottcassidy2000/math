# Actual endpoint responses connect Hamming parity to integer port memory

**PROVED conditional all-prime module theorem / FINITE-EXACT endpoint and
Smith certificates / CITED inherited endpoint spectrum.** Same-agent
review. This is a weighted two-endpoint result; general multi-run strong
attachments still need their complete ordered component data.

## A1 Recover the prior operator before adding its integer module

The niche source is Session54's [Paley endpoint operator](../session54_reframed_numerology_20261025/paley-endpoint-observer/RESULT.md)
and its [independent audit](../session54_reframed_numerology_20261025/audit-paley-endpoint/AUDIT.md).
They already prove the two-channel law, spectrum, quadratic order and a
congruence describing the full generated integral algebra. Their endpoint
counts and the appearance of 1729 are inherited, not new findings here.
The added invariant is the **cokernel of the endpoint map**, which is a
different object from its generated algebra.

For a prime q=3 mod4 let A be the Paley adjacency matrix, r=(q-1)/2, and
E_ij the number of directed Hamiltonian paths from i to j. The inherited
law is

\[
E=aA+bA^{\mathsf T},\quad s=a+b,\quad\delta=a-b,\quad
N={s^2+q\delta^2\over4}.
\]

The constant line has eigenvalue rs; on the augmentation subspace,
E^2+sE+NI=0 and the characteristic polynomial is (X^2+sX+N)^r.
Consequently det E=rs N^r. Redei parity makes s odd, and hence delta
odd and nonzero. The earlier audit's warning remains essential:
Z[E] need not be the direct product of its rational factors.

## A2 A complete integer cokernel under two explicit guards

**Theorem.** If gcd(a,b)=1 and q does not divide s, then

\[
\boxed{\operatorname{coker}E\cong
\mathbb Z/(rs)\ \oplus\ (\mathbb Z/N)^r.} \tag{10}
\]

With d=gcd(rs,N), its q Smith factors are

\[
\boxed{1\ (r\text{ copies}),\ d,\ N\ (r-1\text{ copies}),\ rsN/d.} \tag{11}
\]

The gcd and q-divisibility conditions are hypotheses, not asserted
all-prime facts about Hamiltonian path counts.

**Proof.** Fix a prime l different from q. Over Z_l, the idempotent J/q
splits the lattice into its constant line and its augmentation lattice.
The constant contribution to the cokernel is Z_l/(rs). If l does not
divide N, the augmentation map is invertible.

If l divides N, then s is an l-adic unit. For l=2 this follows from odd s.
For odd l different from q, simultaneous divisibility of s and N would
force l|delta, hence l|a,b, contrary to the gcd hypothesis. Thus the
quadratic X^2+sX+N has two distinct roots modulo l, namely 0 and -s.
Hensel lifting gives roots alpha,beta in Z_l with beta a unit and
v_l(alpha)=v_l(N). Their difference is a unit, so the polynomial
projectors split the augmentation lattice into two free summands. Its
characteristic polynomial gives rank r to each summand. E acts on them
as alpha and beta, contributing (Z_l/(N))^r.

At l=q, both rs and N are units because q does not divide s and
N=s^2/4 modulo q. There is no q-primary cokernel. Combining all primary
parts proves (10). Combining the cyclic factors rs and N by gcd/lcm
gives (11), with d|N|rsN/d and the remaining unit factors as displayed.

This argument classifies a finite abelian group. It does not contradict
the inherited nontrivial integral fiber-product relation for Z[E].

## A3 Exact strong atoms and precision layers

[atom_precision.py](atom_precision.py) counts the actual endpoint matrices
by a subset recurrence, checks all paths by permutations through order7,
and checks selected entries through order12 by independent inclusion-exclusion
over walks. Every reported Smith form has exact integer certification.

| Object | Hamilton paths | Nonunit Smith factors |
|---|---:|---|
| C3 | 3 | none |
| T5 with steps1,2 | 15 | 3 |
| Paley7 | 189 | 106,106,2862 |
| Paley11 | 95095 | 5,1362615,1362615,1362615,1362615,2355961335 |

The endpoint pairs for q3,7,11 are (0,1),(1,8),(628,1101). They satisfy
the hypotheses of A2, and the computed Smith forms equal (11). The
inherited q19 pair (3040229428,3817640437) also satisfies the hypotheses;
an independent integer Smith computation on its reconstructed matrix
agrees with the formula. The q19 Hamiltonian enumeration itself is inherited
from Session54, not rerun here.

For Paley7, N=106=2*53 and rs=27. The exact response-image cost modulo
2^k is therefore

\[
\log_2|\operatorname{im}(E\bmod2^k)|=4k+3\max(k-1,0).
\]

Its parity response has four independent bits; higher precision introduces
three more linear directions. The same valuation pattern occurs at prime53.
At prime3 the cost exponent is 6k+max(k-3,0): one direction is invisible
until the fourth ternary digit.

The prime53 has a concrete inherited spectral source. In Q(sqrt(-7)),

\[
{ -9+7\sqrt{-7}\over2}
={1+\sqrt{-7}\over2}(5+2\sqrt{-7}),
\]

whose factors have norms 2 and 53. This explains the determinant factor,
while the Smith certificate determines how that norm is distributed among
individual integral response directions. It is not a conjectural
prime-distribution pattern.

## A4 A literal map to the seven-prime square-total bank

For Paley7, E=A+8A^T, so E mod2=A. Its row span is the Hamming [7,4,3]
code, with weight enumerator 1+7z^3+7z^4+z^7. Label its vertices 0,...,6
by the primes

\[
(7,2,19,5,59,29,89).
\]

The square classes of p+1, in the basis (2,3,5), are respectively the
nonzero binary vectors

\[
(1,2,4,3,6,7,5),
\]

with integer codes for their three parity bits. These are the powers of
a root of x^3+x+1 in F8. The three out-neighbors at differences1,2,4
sum to zero in this field; hence every Paley adjacency row lies in the
kernel of the three-by-seven square-class matrix. Both spaces have
dimension four, so they are equal.

For a subset I of the seven vertices, the associated squarefree integer
n_I is the product of its prime labels. Multiplicativity gives
sigma(n_I)=product_(i in I)(p_i+1). Therefore

\[
\boxed{I\text{ is a parity endpoint-response codeword}
\iff \sigma(n_I)\text{ is a square}.} \tag{12}
\]

The checker exhausts all128 subsets and confirms the same sixteen words.
For example, the out-neighbors of vertex0 have labels2,19,59, and their
sigma product is 3*20*60=3600. This is a specified coordinate map preserving
the square predicate and XOR of subsets. It does not transport integer
Hamiltonian multiplicities or arbitrary later complement operations to
aliquot dynamics. The extra integer directions in A3 show exactly why
retaining only this parity code is insufficient for weighted responses.

## A5 The 81 and 27 pair separates at the second ternary digit

The exact inherited twelve-vertex pair is reconstructed with its original
labels. It has identical labelled cyclic triples, identical total
Hamiltonian path count462915, and complete automorphism orders81 and27.
Every entry of both endpoint matrices is divisible by3. Every entry of
the 81 matrix is divisible by9, but the 27 matrix has

\[
E_{9,1}=3171\equiv3\pmod9,
\]

where the corresponding 81 entry is2781, divisible by9. Labels are
zero-based. Vanishing of the entire matrix modulo9 is invariant under
any simultaneous relabelling, so this is not an artifact of choosing
different vertex names.

The full 3-adic Smith exponent lists are

```text
81 control: 2,2,2,2,2,2,2,2,2,4,4,6
27 control: 1,2,2,2,2,2,2,2,3,4,4,8
```

Thus their weighted endpoint images modulo3 each contain one response;
modulo9 the images have respectively one and three responses. This is a
sharp precision separator in the declared 3-power channel. Both full
transporter groups are independently enumerated and their actions on
every endpoint entry are checked.

These controls are strong composed tournaments, not prime atoms. The
endpoint matrix suffices for the stated separator; it does not repair the
known six-vertex counterexample to composing arbitrary revisits from only
two endpoints. Joint ordered run ports and realizability remain necessary.

## Reproduction and next classification problem

```text
python -X utf8 -B 05-knowledge/results/overnight_rank_return_20261029/atom_precision.py --write
python -X utf8 -B -O 05-knowledge/results/overnight_rank_return_20261029/atom_precision.py
```

Python's standard library handles endpoint counts and verification.
SymPy supplies compact Smith witnesses for larger matrices, checked by
plain integer multiplication and an independent Bareiss determinant for
unimodularity; a separate elementary Smith producer agrees. The complete
receipt stores matrices, transforms, source rows and the prime-label map.

The next precise target is the analogous integer response module for
ordered multiple-run attachments, with shared endpoints retained. The
all-prime formula A2 suggests testing the exceptional cases gcd(a,b)>1
or q|s separately instead of extrapolating a formula beyond its hypotheses.
