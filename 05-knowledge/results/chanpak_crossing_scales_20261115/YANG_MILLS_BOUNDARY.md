# A reflection-positive tournament seam and its two scaling channels

**Proved here:** an exact rank-two reflection-positive *centered seam*, its square determinant, a uniform all-size boundary factorization, and a cross-scale exclusion criterion on the physical `(H,Q)` profile. These are one-apex tournament statements, not Yang–Mills results or a completeness result for Chan–Pak.

## Source read and precise analogy

The [October 28 OpenAI Lean construction](https://github.com/yang-mills-group/OpenAI-YM-20261028#the-public-challenge) (public November release) selects a Wilson family, takes appropriate fixed-volume mesh limits, and permits a subsequence in the subsequent increasing-volume limit. Its constructed hierarchy provides Wightman reconstruction and separated OPE estimates; a positive mass gap and uniqueness of the general volume limit are not asserted. The [November 10 independent validation](https://alexisdpc.github.io/yang-mills/references/assessment-1-2026-11-10.pdf), §§4–8 and 10, distinguishes the reconstructed correlation hierarchy from a continuum probability measure on connections, documents the OS positivity checks, and does not infer a mass gap. The analogy here is to separate an observable, retain the actual finite models, and control the error uniformly before passing to any limit. Our elementary finite product law does not use their theorem. The repo's [October energy-ladder comparison](../unit_path_ring_kuratowski_pythagorean_yang_mills_20261018.md#5-a-yangmills-reading-sources-the-repository-readme-only) compared observed ladder slopes and a reflection-positivity idea. Here we use the [exact polynomial transfer](../chanpak_heisenberg_target_20261115/GAP_HEISENBERG.md) in place of a fitted law.

## 1. Isolate the seam before taking a positivity test

A spine word (w) of length (n) describes the relation to a single apex: A means the apex beats that vertex, B means the vertex beats the apex; the spine is transitive in word order. For $n\ge1$, set
```
c_w(t) = sum over i<j and (w_i,w_j)=(A,B) of t^(j-i-1).
a_w(t) = sum over A positions i of t^(n-1-i).
b_w(t) = sum over B positions j of t^j.
H(w) = 1 + c_w(2) + c_w(0).
Q(w) = 2^(n+1) + 2^(n-1) [c_w(3/2)+c_w(1/2)].
```
Indices start at zero and `0^0=1` in polynomial evaluation. Indeed, a directed odd cycle uses the apex, an eligible A-to-B pair, and an even subset of the vertices in between. Its contribution to (H) is 2 and to (Q) is $2^{n-k}$ when (k) middle vertices are used: assign its entire odd cycle to one slot and assign all unused vertices individually to either slot. The empty packing supplies the constant terms. This uses the [tournament odd-cycle expansion](https://arxiv.org/html/2412.10572v3#S3.SS4.Thmcorollary20).

Write (w^*) for reverse-complement, and $R(w)=(Q(w)-2^{n+1})/2^{n-1}$. A reverse-complement preserves both (c_w) and the count observables; also $b_{w^*}=a_w$. The seam identity is $c_{uv^*}=c_u+c_v+a_u a_v$, by splitting pairs at the seam. Thus
```
K_H(u,v) = H(uv*) - H(u) - H(v) + 1
         = a_u(2)a_v(2) + a_u(0)a_v(0).
K_R(u,v) = R(uv*) - R(u) - R(v)
         = a_u(3/2)a_v(3/2) + a_u(1/2)a_v(1/2).
```
Both matrices are positive semidefinite on every finite list of halves, with rank at most two. Every two-half principal minor is exactly the square of the two-dimensional oriented area of the corresponding feature vectors. In particular for (K_H) it is an **integer square**. This tests actual splicings, not arbitrary formal operations. The statement concerns the centered seam; no positivity assertion is made about the uncentered count matrix. If halves start A, the glued tournament is strong, since its word ends B as well.

**Exact diagnostic.** Writing $z=a_u(2)$, one has $K_H(u,u)=z^2+(z\bmod2)$; at fixed width the diagonal selects the literal half-word uniquely. The function is strictly increasing on nonnegative integers and the binary digits of (z) encode all letters. Thus fixing reflection-diagonal data itself cannot leave a spare half-word to adjust an absolute checksum. The two (Q) channels are separately visible in the Gram representation, but the joint feature vector still comes from one physical word.

## 2. Two boundary scales, full-sequence limit

Now restrict to strong words starting A and ending B, with $n\ge3$, and use a cutoff $1\le r\le\lfloor n/2\rfloor$. For (t>1), define
```
L_r(t) = sum over 0<=i<r of [w_i=A] t^(-i);
R_r(t) = sum over 0<=k<r of [w_(n-1-k)=B] t^(-k).
```
An exact uniform bound is
```
0 <= c_w(t)/t^(n-2) - L_r(t)R_r(t)
  <= 2 t^(-r)/(1-1/t)^2.
```
After division by (t^{n-2}), every allowed pair has geometric weight $t^{-i-k}$ with (k=n-1-j). The rectangle (i,k<r) has every pair in the correct order. Every pair outside the rectangle has $i\ge r$ or $k\ge r$; bounding by the two full geometric strips proves the result. This controls all words, not merely average words.

For uniform random strong words, couple left A-bits and right B-bits to two independent iid fair bit sequences with fixed first digit 1, taking $r=\lfloor n/2\rfloor$; if necessary fill the middle vertex independently. Define their infinite geometric sums (L(t),R(t)). Then jointly along the **full** sequence of orders,
```
(H/2^(n-2), Q/(2*3^(n-2))) → (L(2)R(2), L(3/2)R(3/2)).
```
The per-coordinate coupling errors are respectively (O(2^{-r})) and (O((2/3)^r)), deterministically. For the corrections use (0\le c_w(0)\le n-1) and (0\le c_w(1/2)\le 2n) by summing the allowed gap weights for each start, and use the geometric tails of (L,R). The limits share the same underlying bits at both scales; they are not claimed independent of each other. In the first coordinate the two endpoint sums are independent uniforms on ([1,2]); its limiting density is (log x) for $1\le x\le2$, (log(4/x)) for $2\le x\le4$, zero otherwise, by integrating (1/u) over the allowed $u\in[1,2]$ with $x/u\in[1,2]$. A smooth macroscopic marginal law does not imply unit-resolution coverage.

## 3. Turn the scaling relation into an exact rejection test

Set $d_t=c_w(t)/t^{n-2}-1$, the defect after removing the compulsory pair from the first to the last spine vertex. Its first active boundary shell (s) is the least positive depth of an A from the left or a B from the right, excluding the two compulsory endpoint letters. There is such a shell for $n\ge3$, and $s\le\lfloor n/2\rfloor$ by inspecting a middle position. Pairing the shell letter with the opposite endpoint gives $d_t\ge t^{-s}$. Each non-endpoint tail is at most $t^{-s}/(1-1/t)$. Bounding by the unconstrained product of full boundary sums and using $s\ge1$ gives
```
2^(-s) <= d_2 <= 6*2^(-s),
(2/3)^s <= d_(3/2) <= 12*(2/3)^s.
```
In particular the weighted channels obey
```
(d_2/6)^beta <= d_(3/2) <= 12*d_2^beta,
beta = log_2(3/2).
```
The exponent is sharp: words (A B^r A^r B) have $s=r+1\to\infty$, and the same bounds force both defects to be of their respective shell scales.

For an alleged physical pair (h=H,q=Q) at spine length (n), compute
```
X = h - 1 - 2^(n-2),
Y = q - 2^(n+1) - 2*3^(n-2).
```
Search $1\le s\le\lfloor n/2\rfloor$. At least one shell has to pass **both** integer tests
```
2^(n-2-s) <= X <= 6*2^(n-2-s) + (n-1),
2^(s+1)*3^(n-2-s) <= Y <= 12*2^(s+1)*3^(n-2-s) + n*2^n.
```
The additive widths account for (c_w(0)) and $2^{n-1}c_w(1/2)$; all exponents are nonnegative. This works unchanged after inferring $Q=(N-(M-L)h)/L$ on the bounded progression of the [anchor construction](../chanpak_quadratic_anchor_atlas_20261115/README.md). It is a warranted reject test before attempting a physical core, not a characterization. For example, take the actual (h) of (AB^rA^rB), but choose (q=6\cdot3^{n-2}+O(24)) even and satisfying ordinary order-mod-3 congruence as well as divisibility by 8. Such pairs pass the elementary size lower bounds for large (r), yet cannot be profiles: the first test forces $s\to\infty$ while the normalized (Q) defect tends to 2. Thus macroscopic accessibility of each coordinate separately loses a real constraint on **joint** target placement.

Run `python replay_ym_boundary.py` for independent direct path/cover checks at small order, exact-rational boundary/Gram checks, and the shell-exclusion family.
