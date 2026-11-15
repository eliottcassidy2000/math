# Bounded CRT amplification: when exact tournament windows compose

This is an exact target selector for a specified physical deck, not an existence theorem for all targets. The missing resource is now quantitative: matched quotient windows for new coprime spine counts, and enough movement of the *absolute* baseline. Reusing finitely many cores at greater heights cannot supply it.

## The tournament interface

Take a nonempty tournament spine $R$, mark the vertices beaten by a common apex, and write $n_R$ for the Hamilton-path count after adding that apex. Keep the full exterior profile $(h,Q,s,e,t,f)$: ordinary spine count $h$, ordered two-path cover count $Q$, total and signed marked starts $(s,e)$ and unmarked ends $(t,f)$. The [two-sided recontact](../chanpak_tail_return_20261115/TWO_SIDED_RECONTACT.md) describes the profiles and proves their composition rule. In particular, equal exterior profiles permit changing $n_R$ without changing any subsequent or preceding block response.

Join $r$ such spines in order, using **one** common apex. For block $i$, all its palette variants have the same exterior profile and provide values $b_i+2D_i$ for a finite integer set $D_i$ containing 0. Let $H=\prod_i h_i$, $w_i=H/h_i$, and let $N_0$ be the literal apex count with all base variants $b_i$. Allowing painted transitive tails outside the blocks changes $N_0$ but not the weights. Then the possible full counts are **exactly**
$$
N_0+2\sum_i w_i D_i.
$$
Indeed replacing $n_i$ adds the change times the old spine count to the left, and every later spine multiplies it by its own count. All cross-seam terms are unchanged at every step. This is the weighted form used by the [bounded Bézout pair](../chanpak_tail_return_20261115/RECONTACT_ENERGY.md); the result below supplies a whole interval criterion and an exact decoder.

If a prime divides at least two $h_i$, it divides every $w_i$. More exactly $\gcd_i(w_i)=H/\operatorname{lcm}_i(h_i)$. This immediately rules out adjacent odds within such an independent product fiber. It does not rule out several fibers having different baselines that interlace. A block with $h=1$ is transitive by the positive odd-cycle packing expansion. For a transitive order, all individual weighted starts are distinct powers of two, so $s$ uniquely identifies the marker set. It has no nontrivial matched palette. Thus a repeated finite stock of nontransitive blocks cannot give an indefinitely deep sequence of coprime amplifiers of this kind.

## Exact prescribed-target test

Assume the $h_i$ are pairwise coprime, omitting trivial $h_i=1$ modules. For a target $N$, take $y=(N-N_0)/2$ and put
$$
r_i \equiv w_i^{-1}y \pmod{h_i},\quad 0\le r_i<h_i,\qquad
L=\frac{y-\sum_i w_i r_i}{H}.
$$
Then $L$ is an integer. Set $U_i(r_i)=\{q:r_i+h_i q\in D_i\}$. The target is realizable in this fiber if and only if $L\in\sum_i U_i(r_i)$. Reduction modulo $h_i$ forces each digit into this single fiber, and substituting $r_i+h_i q_i$ proves sufficiency. The output gives every local **physical** variant to use; it is not a modular witness.

It is enough to prove each selected quotient fiber contains a certified interval and that their summed endpoints bracket $L$. No other residue fibers need be filled for this particular target. Thus the strategy can seek *target-selected* local palettes rather than an unnecessarily ambitious complete palette at every module.

## Sharp common window for box palettes

If $D_i=[0,m_i]\cap\mathbb Z$ and $m_i\ge h_i-1$, write
$$
R=\sum_i w_i(h_i-1),\qquad
M=\sum_i w_i m_i,\qquad C=R-H+1.
$$
Every half-offset in the interval $[C,M-C]$ is realized. For each residue vector its least representative $B=\sum_i w_i r_i$ lies in $[0,R]$, and its largest representative $E=\sum_i w_i(r_i+h_i\lfloor(m_i-r_i)/h_i\rfloor)$ lies in $[M-R,M]$. Every point of the same class mod $H$ between $B$ and $E$ occurs, by summing intervals of allowed quotient digits. If $y\ge R-H+1$, then $y>B-H$ forces $y\ge B$; the upper endpoint is identical by reflection. Moreover $C-1=R-H$ would force all $r_i=h_i-1$, whose least representative is $R$. This lower adjacent point is therefore impossible, as is $M-C+1$ by reflection. The assertion about the interval is vacuous when its endpoints are reversed. For two blocks $C=(h_1-1)(h_2-1)$, the usual two-generator Frobenius conductor. If some $m_i<h_i-1$, the missing residues mod $h_i$ prohibit any interval of $h_i$ half-offsets.

A scale-stable sufficient threshold falls out:
$$
m_i\ge 2h_i-2\ \ \forall i
\quad\Longrightarrow\quad
M-2C=2H-2+\sum_i w_i\,\big(m_i-(2h_i-2)\big).
$$
The same two-sided exterior profile is inherited by every variant of the whole glued spine. Thus a window that meets this threshold can be used as one new block with the same property; the genuine surplus scales by the complementary spine product. This is a *local-to-absolute* statement: full CRT residues are necessary, but the bounded integer carries are what fill the interval.

This theorem does not position $N_0$. For a proposed family $\alpha$, check both (i) the decoder yields a filled interval, (ii) the actual baseline plus that interval contains the requested $N$. In an ordered list of odd windows it suffices to bound successive lower endpoints by two plus the preceding upper endpoints. Neither widening at an uncontrolled baseline nor dense residue sets prove that inequality. The published two-sided witnesses at $h=15$ and 45 have half-widths 1 and 2, far below even the first residue threshold $h-1$. The theorem reduces an induction architecture to a new physical supply question, rather than asserting such a supply exists.

## A proof-strategy obstruction for fixed fan windows

A separate distinction concerns the [parallel-fan windows](../chanpak_leaf_fans_20261115/PARALLEL_FANS.md). For any finite bank of $b$ rays $L2^j+c$, all its values occupy at most $bk$ residue classes modulo $M_k=2^k-1$. Given any fixed reduced prime progression $a\bmod m$, choose a sufficiently large **prime** $k$ so that $\gcd(M_k,7m)=1$; for any odd prime $p$ dividing $7m$, its order of 2 is greater than 1 and cannot divide arbitrarily large prime $k$. The unit group modulo $M_k$ has size at least $\sqrt{M_k/2}>2bk$, by the product formula for the totient (and in fact at prime $k$ it has asymptotic density 1). Since multiplication by 7 permutes that group, some unit residue $u$ has both $u$ and $7u$ outside the ray bank. CRT and [Dirichlet's theorem](https://arxiv.org/abs/0808.1408) put infinitely many primes $p$ in the combined progression $p\equiv a\bmod m,\ p\equiv u\bmod M_k$. For these, both $p$ and $7p$ miss the bank *even after ordinal products of arbitrarily many bank members*: a product realizing $p$ or $7p$ must supply the full target as a bank member, since 7 has no physical witness. Finite fixed-count extra seeds can simply be added to the banned residues.

This identifies an explicit prime/7p lane needing genuinely different sources; it is not a theorem that any of those primes is missing from the unrestricted spectrum. It complements the [earlier least-missing analysis](../chanpak_least_missing_jacobian_20261114/BLACK_BOX_AND_SYMMETRY.md) and the fixed-grammar modular saturation counterexample: one proof architecture has visible congruence obstructions, while another has none but still cannot place absolute targets. The CRT lift avoids both diagnoses **once** suitable physical quotient windows are supplied.
