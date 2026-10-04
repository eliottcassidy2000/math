# The 63-clock: tensor coordinates, subfield norms, and a retained ternary digit

Status: **PROVED** elementary finite-field statements below; **FINITE-EXACT** exhaustive controls over all 64 field elements. The sextic field, its residue field, and the generator (b) are inherited from [sixth_clock_branches_20261004.md](sixth_clock_branches_20261004.md). This note does not assert a conjugacy with Collatz iteration or transfer a finite-field clock to real sixth-root branches.

Reproduction:

```text
python 04-computation/experiments/sextic_subfield_clock_20261004.py
python -O 04-computation/experiments/sextic_subfield_clock_20261004.py
```

The two commands have identical output, saved in [sextic_subfield_clock_20261004.out](sextic_subfield_clock_20261004.out). Checks use explicit exceptions and remain active under `-O`.

## 1. Inheritance and the exact object

The closest inherited mechanism is the primitive multiplicative residue clock in the incoming sixth-root work. If (t^3=1+3\phi), then the integer minimal polynomial is (t^6-5t^3-5). Its reduction modulo two is

\[
 f(t)=t^6+t^3+1=\Phi_9(t),\qquad
 F=\mathbf F_2[t]/(f)=\mathbf F_{64}.
\]

Indeed the roots of (f) have order nine, and the order of two modulo nine is six, proving irreducibility. The inherited element (b=1+t) has order 63: (b^{63}=1), while (b^{21}\ne1) and (b^9\ne1). The finite script independently represents elements by six binary polynomial coefficients and checks all 4096 products against the incoming sextic power-basis implementation.

Our live objects are the six-bit additive vector space, the multiplicative 63-clock, its two proper subfields, the subfield norm quotient, and the ternary branch coordinate. Their operations must remain distinct. The closest hostile is (t): it uses the same six-bit field representation but has order nine, not 63. The corrected near miss is that an element of (F) need not be a single product of one element from each proper subfield. The previously unused sidecar is the second ternary digit of a phase modulo nine.

Set

\[
 \varphi=b^{21}=1+t^3,\qquad
 \gamma=b^9=t^5+t^2+t.
\]

Here \(\varphi\) denotes the residue of the golden integer, not a positive real number. The subfields are

\[
 F_4=\{x:x^4=x\}=\{0,1,\varphi,1+\varphi\},\quad
 F_8=\{x:x^8=x\}=\{0,1,\gamma,\ldots,\gamma^6\},
 \quad F_4\cap F_8=F_2.
\]

The multiplicative orders in (F^\times) have the complete census

| order | 1 | 3 | 7 | 9 | 21 | 63 |
|---|---:|---:|---:|---:|---:|---:|
| number of elements | 1 | 2 | 6 | 6 | 12 | 36 |

## 2. A tensor product is not a single product

Multiplication gives an isomorphism of (F_2)-algebras

\[
 F_4\otimes_{F_2}F_8\longrightarrow F,\qquad u\otimes v\longmapsto uv.
\]

To prove this directly, (X^2+X+1) has no root in (F_8), whose multiplicative group has order seven. Thus adjoining \(\varphi\) to (F_8) gives all 64 elements, with basis (1,\varphi) over (F_8). Since \(\gamma\notin F_2\), it has degree three and (1,\gamma,\gamma^2) is an (F_2)-basis of (F_8). Consequently

\[
 (1,\gamma,\gamma^2,\varphi,\varphi\gamma,\varphi\gamma^2)
\]

is a six-element binary basis. The natural algebra map is surjective between six-dimensional vector spaces and hence bijective.

Each field element therefore has a unique binary (2\times3) coefficient matrix. The matrix ranks have counts (1,21,42) for ranks zero, one, two. In particular:

* The **tensor product** has dimension (2\cdot3=6) and 64 elements.
* The **additive Cartesian product** (F_4\times F_8) has dimension (2+3=5) and 32 elements.
* The set of **single products** (uv) has 22 elements, including zero.

For the last claim, nonzero factorizations are unique: if (uv=u'v'), then (u/u'=v'/v\in F_4^\times\cap F_8^\times=\{1\}). The 21 nonzero products form

\[
 H=F_4^\times F_8^\times=\langle b^3\rangle.
\]

Hence, for (x\ne0),

\[
 \operatorname{rank}_{\otimes}(x)=1
 \iff x\in H
 \iff x\text{ is a cube}
 \iff x^{21}=1.
\]

The three cosets (H,bH,b^2H) partition the 63 nonzero states into three sets of 21. Multiplication by (b) is binary-linear, so these are three linear images of the same rank-one tensor set. Rank two consists of the six elements of order nine and the 36 primitive elements. Therefore rank two alone does not certify primitivity; (t) is an explicit hostile.

## 3. The two norms lose exactly one ternary phase digit

The relative norms are

\[
 N_4(x)=x^{1+4+16}=x^{21},\qquad N_8(x)=x^{1+8}=x^9.
\]

For (x=b^k), their joint value is \((\varphi^k,\gamma^k)\). Thus

\[
 \pi:F^\times\longrightarrow F_4^\times\times F_8^\times
\]

is a surjective homomorphism with kernel (F_4^\times=\langle b^{21}\rangle), and each of its 21 fibres contains three elements. The norm pair retains the cube exactly:

\[
 \boxed{x^3=N_4(x)N_8(x)^{-2}.}
\]

It even retains the predicate of being primitive:

\[
 \operatorname{ord}(x)=63
 \iff N_4(x)\ne1\ \text{and}\ N_8(x)\ne1,
\]

because these conditions are (3\nmid k) and (7\nmid k). It does not retain every multiplicative order: the fibre over ((1,1)) contains (1,\varphi,\varphi^2), of orders (1,3,3).

There is no multiplicative section of \(\pi\). The target element \((\varphi,1)\) has order three, but all three of its lifts have order nine: in the script they are (t^2,t^5,t^5+t^2). A homomorphism section would have to lift it to an element whose order divides three. An arbitrary choice of representatives is possible, but it is not multiplicative; the carry below is unavoidable for that choice.

Nor is this an additive projection. For example (N_4(1+\varphi)=1), whereas (N_4(1)+N_4(\varphi)=0). This is a quotient of multiplicative phase, not deletion of one binary vector coordinate.

Cubing and taking sixth powers have the same three-element kernel and image (H), since squaring is a field automorphism. Thus (X^6-h), for (h\in H), has three distinct field roots, each with multiplicity two; outside (H) it has none in (F). In particular (t^6=\varphi), with exactly the roots

\[
 t,\quad t\varphi,\quad t\varphi^2
 \quad (=t,\ t^4+t,\ t^4\text{ in some order}).
\]

This finite-characteristic root fibre is a different object from six complex sixth-root branches or the unique positive real sixth root.

## 4. A complete 3-by-3-by-7 clock with a carry

Choose

\[
 \rho=b^{28}=t^5,\qquad
 \delta=b^{36}=(1+t)t^4.
\]

Then \(\operatorname{ord}(\rho)=9\), \(\operatorname{ord}(\delta)=7\), \(\rho\delta=b\), \(\rho^3=\varphi\), and \(\delta^2=\gamma\). Each nonzero element has a unique code

\[
 (i,j,\ell)\in\{0,1,2\}^2\times\{0,\ldots,6\},
 \qquad x=\rho^{i+3j}\delta^\ell.
\]

Its norms are \((\varphi^i,\gamma^\ell)\); (j) is precisely the lost branch. This decoder needs only the three- and seven-element subfield logarithms: obtain (i,\ell) from the norms, then the residual (x/(\rho^i\delta^\ell)=\varphi^j) determines (j).

For (x=b^k), the code says

\[
 k\bmod9=i+3j,\qquad k\bmod7=\ell.
\]

One multiplication by (b) is therefore

\[
 (i,j,\ell)\longmapsto
 ((i+1)\bmod3,\ (j+[i=2])\bmod3,\ (\ell+1)\bmod7).
\]

Multiplication of two codes is

\[
 (i,j,\ell)(i',j',\ell')=
 \left((i+i')\bmod3,
 \left(j+j'+\left\lfloor\frac{i+i'}3\right\rfloor\right)\bmod3,
 (\ell+\ell')\bmod7\right).
\]

All 3969 nonzero products are checked. After 21 steps the norms return and (j) advances by one; after 63 the complete state returns. The code is a Cartesian set of size (3\cdot3\cdot7), but its group operation is not the componentwise operation of (C_3\times C_3\times C_7). It implements (C_9\times C_7\cong C_{63}).

This is the exact missing depth in (63=3^2\cdot7). The prime three already occurs in (2^2-1), and seven in (2^3-1), while \(\operatorname{ord}_9(2)=6\). The statement is elementary here; no general primitive-divisor theorem is needed. The (\delta\) coordinate also agrees with the (a=2) normalization (c=(1+\zeta)\zeta^{(3^a-1)/2}) of the companion [cyclotomic_depth_towers_20261004.md](cyclotomic_depth_towers_20261004.md).

## 5. The faithful six-bit Fourier interface

Let \(\operatorname{Tr}(x)=x+x^2+x^4+x^8+x^{16}+x^{32}\). The 64 additive characters are

\[
 \chi_y(x)=(-1)^{\operatorname{Tr}(yx)},\qquad y\in F.
\]

The trace pairing is nondegenerate: the nonzero trace polynomial has degree 32, so cannot vanish on all 64 elements; choose (w) of trace one, and for (y\ne0), use (x=w/y). A nontrivial character has zero sum by translation through an element on which its sign is negative. Consequently all 64 character rows are orthogonal, and a full additive Fourier transform retains the complete function on six-bit states.

The multiplicative clock acts faithfully on this character indexing:

\[
 \chi_y(bx)=\chi_{by}(x).
\]

It fixes the trivial character and cycles all 63 nontrivial ones. This is an exact interface between a linear six-bit clock and its additive Fourier representation, not an identification of addition and multiplication.

For (u,u'\in F_4) and (v,v'\in F_8), the tensor pairing factors:

\[
 \operatorname{Tr}_{64/2}(uu'vv')
 =\operatorname{Tr}_{4/2}(uu')\operatorname{Tr}_{8/2}(vv').
\]

Indeed the six Frobenius exponents enumerate all pairs of exponents modulo two and three. In the basis from section 2, its binary Gram matrix is (G_4\otimes G_8), where

\[
 G_4=\begin{pmatrix}0&1\\1&1\end{pmatrix},\qquad
 G_8=\begin{pmatrix}1&0&0\\0&0&1\\0&1&0\end{pmatrix}.
\]

This is a six-by-six matrix of a bilinear form. It is **not** the 32-by-32 additive Fourier matrix (H_4\otimes H_8) for the Cartesian group (F_4\times F_8). The additive group of (F_4) is the same abstract Klein four group used in the tetrahedral cut coordinates, after a choice of binary basis; its order-three field multiplication is an additional operation, not the Klein four group operation.

A single trace bit loses information. Six successive trace bits recover it. The observation map

\[
 x\longmapsto(\operatorname{Tr}(x),\operatorname{Tr}(bx),\ldots,\operatorname{Tr}(b^5x))
\]

is an invertible binary-linear map: (b) has degree six, so (1,b,\ldots,b^5) is a basis, and the trace pairing is nondegenerate. The 63 cyclic length-six windows in (s_k=\operatorname{Tr}(b^k)) are exactly all nonzero six-bit words. The single-bit balance is 31 zeroes and 32 ones. The inherited primitive polynomial gives

\[
 s_{k+6}=s_{k+4}+s_{k+3}+s_{k+1}+s_k\quad\text{in }F_2.
\]

## 6. Connection ledger and verification boundary

| source and map | target and preserved predicate | lost data and required sidecar |
|---|---|---|
| Six-bit polynomial state; tensor-basis change | Binary (2\times3) coefficient matrix; all field operations transported exactly | None, with the stated bases retained |
| (F^\times\); pair of relative norms | 21 norm pairs; cube and primitivity predicate | One (C_3) branch; retain (j) and generator anchors |
| (F^\times\); three-coordinate code | 63 phase states; exact multiplication and clock action | None; preserve the ternary carry rather than componentwise addition |
| A function on the additive six-bit space; full trace-character transform | 64 Fourier coefficients; exact invertibility and clock permutation | None for the full transform; magnitudes alone are not asserted sufficient |
| A field state; six consecutive trace outputs | Six bits; exact state | None with clock and observation basis; one trace bit alone is insufficient |

The finite universe is the entire 64-element field, all 4096 products, all 3969 nonzero code products, every norm fibre, every sixth-power fibre, all 4096 character inner products, and all 63 six-bit trace windows. Positive controls include the inherited primitive element, the two explicit subfields, and the exact decoder. Hostile controls include (t) of tensor rank two and order nine, failure of additive norms, failure of a multiplicative norm section, and loss of the (j) digit after a 21-step visible return. The general explanations above prove why these checks succeed or fail; they are not inferences from a sample of integers.

The map is ready to carry finite phase memory. It supplies neither an integer Collatz route nor a real/integer lift from residue data. Those would require a separately specified source map and guards. All field arguments here are elementary derivations; no literature-priority claim is made.
