# A next-order cusp phase with a controlled matching error

2026-10-19. **PROVED; independent audit completed.** This extension resolves hypothesis (H) of the [initial phase checkpoint](../parabolic_escape_precision_20261019/parabolic_escape_precision_20261019.md). The argument has received mathematical review and a [separate independent audit](review/independent_phase_review.md). No worldwide novelty claim is made.

## 1. Statement

Use the initial note's conventions: x₀=0, xₖ₊₁=xₖ²+1/4+ε, R>1/2 fixed, and εn is the unique positive parameter with xₙ(εn)=R. Let δn=√εn and let C_R be the convergent endpoint phase in equation (10) of that note.

**Theorem.** As n→∞,

\[
n=\frac\pi{\delta_n}+C_R+\frac\pi3\delta_n
 +O_R\!\left(\delta_n^2\log\frac1{\delta_n}\right).
\tag{1}
\]

In particular, putting M=n−C_R,

\[
\epsilon_n=\frac{\pi^2}{M^2}+\frac{2\pi^4}{3M^4}
 +O_R\!\left(\frac{\log n}{n^5}\right),
\tag{2}
\]

and the actual adjacent escape-band widths satisfy

\[
\epsilon_n-\epsilon_{n+1}
=\frac{2\pi^2}{n^3}
 +\frac{(6C_R-3)\pi^2}{n^4}
 +O_R\!\left(\frac{\log n}{n^5}\right).
\tag{3}
\]

The logarithmic remainder here is an upper bound; no claim is made that a logarithmic term is actually present.

## 2. An improved exact clock and its defect

Write e=δ², a=−1/2, b=R−1/2, q=y²+e, and fδ(y)=y+q. Set

\[
\widetilde\Psi_\delta(y)
=\frac{1+e/3}{\delta}\arctan\frac y\delta
 +\frac12\log(y^2+e)-\frac y2.
\tag{4}
\]

By the first note's exact arctangent identity,

\[
\widetilde\Psi_\delta(f_\delta(y))
 -\widetilde\Psi_\delta(y)=1+E_e(y),
\tag{5}
\]

where, for w=1+y,

\[
A_e(y)=\frac{\arctan(\sqrt e/w)}{\sqrt e}
       =\int_0^1\frac{w}{w^2+e t^2}\,dt,
\]

\[
E_e(y)=(1+e/3)A_e(y)+\tfrac12\log(w^2+e)
                  -1-\tfrac12(y^2+e).
\tag{6}
\]

The integral defines A smoothly also at e=0. All estimates below are uniform on the fixed compact interval y∈[−1/2,b] and 0≤e≤e₀, for e₀ small enough.

There is a constant K_R such that

\[
|\partial_yE_e(y)|\le K_R(y^2+e),\qquad
|\partial_eE_e(y)|\le K_R(|y|+e),
\tag{7}
\]

\[
|E_e(y)|\le K_R(|y|^3+e|y|+e^2)
          \le K_R(|y|+e)(y^2+e).
\tag{8}
\]

**Verification.** The exact derivative is

\[
\partial_yE_e(y)
=\frac{y-e/3}{(1+y)^2+e}-y
=\frac{-2y^2-y^3-ey-e/3}{(1+y)^2+e},
\]

which proves the first estimate in (7). Smoothness and ∂eE₀(0)=0 prove the second by the mean-value theorem on a fixed compact rectangle. For completeness,

\[
\left.\partial_eE_e(y)\right|_{e=0}
=\frac1{3w}-\frac1{3w^3}+\frac1{2w^2}-\frac12,
\]

which vanishes at w=1. Also

\[
E_0(y)=D_0(y)-\tfrac12y^2=O_R(|y|^3),
\]

since E₀(0)=E₀′(0)=E₀″(0)=0. Integrating the second estimate in (7) from 0 to e proves (8).

The key improvement over the first clock is cubic, rather than quadratic, vanishing of the unperturbed defect, together with cancellation of the order-e defect at y=0.

## 3. Uniform endpoint comparison out to order 1/δ

Let uₖ=f₀ᵏ(a)<0 and vⱼ=g₀ʲ(b)>0, with g₀ the positive inverse of f₀. Their reciprocal identities give

\[
\frac1{2k+2}\le |u_k|\le\frac1{k+2},
\tag{9}
\]

\[
\frac1{b^{-1}+j}\le v_j
\le\frac1{b^{-1}+j/(1+b)}.
\tag{10}
\]

For (9), tₖ=−1/uₖ starts at 2 and obeys tₖ₊₁−tₖ=1+1/(tₖ−1)∈[1,2]. For (10), 1/vⱼ−1/vⱼ₋₁=1/(1+vⱼ)∈[1/(1+b),1].

At a threshold parameter e=δ²=εn, write y₀=a,…,yₙ=b for its full orbit. Choose

\[
m=\left\lfloor\frac1{10\delta}\right\rfloor.
\tag{11}
\]

For all sufficiently small δ, depending only on b, m≥1, 2m<n, and the following estimates hold for 0≤k≤m and 1≤j≤m:

\[
u_k\le y_k\le u_k+k e\le u_k/2<0,
\tag{12}
\]

\[
0<v_j/2\le v_j-j e\le y_{n-j}\le v_j.
\tag{13}
\]

The factor 1/2 in (12) is only asserted for k≥1; at k=0 the displayed weak inequalities still hold because u₀<0. The condition 2m<n follows from the already proved nδ→π.

**Proof of (12).** Until crossing zero, f₀ is increasing and 1-Lipschitz on [−1/2,0]. Consequently the nonnegative forward perturbation difference is at most ke. Inductively it remains on that interval because (9) gives

\[
ke\le |u_k|/2
\quad\hbox{whenever}\quad k(k+1)e\le1/4.
\]

For k≤m the latter is true for sufficiently small δ: m(m+1)δ²≤1/100+δ/10.

**Proof of (13).** The last points are positive inverse iterates for as long as the estimates show positivity. The upper comparison follows from fδ≥f₀ and monotonicity. For the lower comparison use L=vⱼ−je>0. A direct expansion gives

\[
f_\delta(L)
=v_{j-1}-(j-1)e-2jev_j+j^2e^2
\le v_{j-1}-(j-1)e
\]

provided je≤2vⱼ. By induction the last expression is at most yₙ₋ⱼ₊₁. Since L>0 and fδ is increasing there, its positive inverse exists and gives yₙ₋ⱼ≥L. Finally je≤vⱼ/2 follows from (10) as soon as j(j+b⁻¹)e≤1/2, which holds uniformly for j≤m at sufficiently small δ. This completes the induction, including existence of the positive inverse branches.

In particular both cut points have magnitude O_R(1/m)=O_R(δ), and all paired endpoint coordinates have magnitude O_R(1/(k+1)) or O_R(1/(j+1)).

## 4. Summing the endpoint, middle, and tail errors

From (7), (9), and (12),

\[
|E_e(y_k)-E_0(u_k)|
\le K_R\left((u_k^2+e)ke+e(|u_k|+e)\right)
\le K_R\left(\frac e{k+1}+(k+1)e^2\right).
\]

Thus the first m defects differ from their unperturbed counterparts by

\[
O_R(e\log(m+1)+e^2m^2)
=O_R\!\left(\delta^2\log\frac1\delta\right).
\tag{14}
\]

Exactly the same estimate, using (10) and (13), holds for the last m defects, paired with vⱼ for j=1,…,m.

On the middle block m≤k<n−m, monotonicity and the cut-point bounds give |yₖ|=O_R(δ). From (8) and the **exact** telescoping of q,

\[
\sum_{k=m}^{n-m-1}|E_e(y_k)|
\le O_R(\delta)\sum_{k=m}^{n-m-1}(y_k^2+e)
=O_R(\delta)(y_{n-m}-y_m)=O_R(\delta^2).
\tag{15}
\]

Finally the omitted unperturbed tails are O_R(δ²). For example,

\[
\sum_{k\ge m}|E_0(u_k)|
\le K_R|u_m|\sum_{k\ge m}u_k^2
=K_Ru_m^2=O_R(\delta^2),
\]

and the reverse tail follows from Σⱼ>mvⱼ²=vₘ.

Combining (14)–(15) proves the quantitative matching statement

\[
\sum_{k=0}^{n-1}E_e(y_k)
=\sum_{k\ge0}E_0(u_k)+\sum_{j\ge1}E_0(v_j)
 +O_R\!\left(\delta^2\log\frac1\delta\right).
\tag{16}
\]

This is the substantive step absent from a merely formal modified-vector-field expansion.

## 5. Matching the constant and extracting the correction

At the fixed endpoints a,b, (4) yields

\[
\widetilde\Psi_\delta(b)-\widetilde\Psi_\delta(a)
=\frac\pi\delta+h(b)-h(a)-\frac{b-a}{2}
 +\frac\pi3\delta+O_R(\delta^2).
\tag{17}
\]

Since E₀(y)=D₀(y)−y²/2 and the squared increments telescope,

\[
\sum_{k\ge0}E_0(u_k)+\sum_{j\ge1}E_0(v_j)
=\sum_{k\ge0}D_0(u_k)+\sum_{j\ge1}D_0(v_j)-\frac{b-a}{2}.
\tag{18}
\]

Sum the exact identity (5), then substitute (16)–(18). The constants combine to the **same** C_R as in the initial phase theorem, and (1) follows.

## 6. Inversion and adjacent differences

Put t=π/δn and M=n−C_R. Equation (1) reads

\[
M=t+\frac{\pi^2}{3t}+O_R(\log t/t^2).
\]

Thus

\[
t=M-\frac{\pi^2}{3M}+O_R(\log M/M^2),
\]

and expanding π²/t² proves (2). Taking adjacent differences in (2) gives

\[
\epsilon_n-\epsilon_{n+1}
=\frac{2\pi^2}{M^3}-\frac{3\pi^2}{M^4}
 +O_R(\log n/n^5).
\]

Here no unproved derivative estimate for the remainder is used: the difference of two O(log n/n⁵) terms is still O(log n/n⁵). Expanding M=n−C_R proves (3).

## 7. Relation to the local complex period

At the split fixed point iδ, the multiplier is 1+2iδ. Its local linearizing time period has expansion

\[
\frac{2\pi i}{\log(1+2i\delta)}
=\frac\pi\delta+\pi i+\frac\pi3\delta+O(\delta^2).
\]

The real coefficient π/3 agrees with (1). The elementary proof above does not assume a global embedding of the quadratic map into a vector field, nor does it infer endpoint matching from this local period. Instead, the improved clock captures the local correction and estimates the global matching defect directly.

The initial checkpoint's numerical controls are compatible with (1), but none is used in its proof. A useful next question is whether the logarithm can be removed and the order-δ² coefficient computed from the two endpoint Fatou coordinates. That stronger assertion is not established here.
