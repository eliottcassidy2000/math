# What survives in the primitive lanes

A least extra hole in the tournament path-count spectrum is prime or seven
times a prime. This reduction is inherited: see the
[earlier least-missing analysis](../chanpak_least_missing_jacobian_20261114/BLACK_BOX_AND_SYMMETRY.md).
It only uses multiplicative closure under ordinal join and the attained
values 49, 63, and 343. The factorization explains why 7 is central but does
not put 7 or 21 back into the spectrum: products requiring a factor 7 are
not available, so $7p$ is an extra primitive lane.

Here are two different obstructions to extrapolating exact windows in
these lanes. Neither is a missing-value result for unrestricted tournaments.

## Fixed fan windows leave a common prime/7p progression

Use a finite physical bank of $b$ rays $L2^j+c$, such as the
[parallel-fan windows](../chanpak_leaf_fans_20261115/PARALLEL_FANS.md),
and any finite collection of exceptional single counts. Modulo
$M_k=2^k-1$, a ray occupies at most $k$ residue classes. Given a fixed reduced
prime progression $a\bmod m$, choose a large prime $k$ so that
$\gcd(M_k,7m)=1$. For each odd prime $q\mid7m$, the multiplicative order of 2
modulo $q$ is greater than 1, and cannot divide arbitrarily large prime $k$.
The unit group modulo $M_k$ has size at least $\sqrt{M_k/2}$:
$\varphi(q^e)^2/q^e=q^{e-2}(q-1)^2\ge1$ for odd primes, with the only
smaller possible factor being $1/2$ at $q^e=2$.
Choose $k$ so large that two copies of the banned set still do not exhaust
the units. Then some unit $u$ has both $u$ and $7u$ outside the banned set.

CRT and [Dirichlet's theorem](https://arxiv.org/abs/0808.1408) give
infinitely many primes $p\equiv a\bmod m,\ p\equiv u\bmod M_k$.
Both $p$ and $7p$ miss the bank, even after arbitrary ordinal products
of its members: to factor either without using it as a bank member would
require a nonunit factor 7. Neutral transitive padding is free but changes
no count. This proves actual common missing progressions for this proof
architecture, not merely a zero-density bound.

## Repeating cores can destroy a primitive lane even when the baseline varies

The independent-switch obstruction uses gcd of *differences* within a fixed
palette product. In contrast, the following congruence is about the full
Hamilton count and holds with **any** apex masks and outer ordinal contexts.

Let one apex be marked over an ordinal succession of nonempty blocks.
For a modulus $M>1$, suppose blocks $i,j$ satisfy
$H(B_i),H(B_j)\equiv0\pmod M$, and some block $k$ strictly between them
satisfies $Q_2(B_k)\equiv0\pmod M$, where $Q_2$ counts ordered two-path
covers allowing empty slots. Then the full Hamilton count is
$0\pmod M$, for every choice of marks and exterior ordinal blocks.

Indeed, write $(h,N,S,d)$ for the state of the prefix, with ordinary count,
full count, total and signed starts. The established [exact cut
transfer](../chanpak_selective_target_20261114/CUT_STATE.md) for a next
block with ordinary count $a$ and two-cover count $q$ is
$$
(h',N',S',d')=(ah,\ aN+(n_B-a)h+tS+fd,\ qS+hs,\ he).
$$
After $i$, $h=0$ modulo $M$ and remains so. At the next block $d$ is
zero and remains so. At $k$, $S$ is zero and remains so. At $j$,
$N$ is therefore zero, and every subsequent row preserves the all-zero
state. No conditions on the marked coefficients are needed.

An immediate corollary is an **absolute** obstruction for repeated
$p$-divisible cores with $p\mid Q_2$: three copies at arbitrary positions
suffice, even surrounded by arbitrary other painted blocks. Three copies
of the regular pentagon spine $(H,Q_2)=(15,120)$ force divisibility by
15. Three copies of the eight-vertex $H=9$ spine of the [literal
example](COPRIME_SWITCH.md), with $Q_2=672$, force divisibility by 3.
Thus repeating either particular reservoir cannot supply large prime or
prime-times-seven counts by widening a one-apex product; changing marks
or changing the surrounding baseline does not repair this specific
obstruction.

The contrast is useful: a modularly complete grammar can still miss
absolute targets, as in the earlier [boundary
example](../chanpak_least_missing_jacobian_20261114/README.md), but a
scaled repetition may fail earlier by killing the entire state. A
minimal-counterexample proof needs an operation that generates a new
absolute witness in the appropriate primitive lane, not more copies of
an ordinal factor.

## Fixed count bounds the hidden ambiguity of the start readback

There is a separate all-order bound on the supply of independent digits.
At any fixed ordinary spine count $h$, the number of different apex
counts in a complete exact profile fiber with two-cover count $Q$ is
$$
O_h\big((1+\log_2 Q)^{\lfloor\log_3 h\rfloor}\big).
$$
For prime $h$ the bound is $O_h(1+\log_2 Q)$. This counts possible
literal variants, not primes or numerical congruence classes.

To prove it, use the ordinal factorization of strongly connected
components. Every nonsingleton component has count $\ge3$ by the positive
odd-packing expansion; for a fixed count its size is bounded by
[Busch's theorem](https://doi.org/10.37236/1141). For fixed $h$ there are
therefore only finitely many ordered lists of nonsingleton components
and choices of their masks. Their length $r$ is at most $\lfloor\log_3
h\rfloor$ (at most one for prime $h$). Once this marked list is chosen,
the two-cover count factorizes ordinally, so the given $Q$ determines
the total number of singleton components: each contributes a factor
two. There are $O((1+\log_2Q)^r)$ distributions of those singleton
positions over the $r+1$ gaps.

Fix a distribution and the marked nontrivial components. Their
contributions to total weighted start $S$ are known. Marking a singleton
contributes the product of the ordinary Hamilton counts of preceding
blocks and the two-cover counts of following blocks. Between two
successive singletons the earlier such weight is at least twice the
later: the ratio is two times the product of $Q_2(C)/H(C)$ over the
intervening nonsingleton components, and each ratio is at least one.
Consequently every earlier weight exceeds the sum of all later singleton
weights, even in the absence of intervening components. The fixed value
of $S$ uniquely determines these singleton marks as a superincreasing
subset sum. This proves the bound even if the end response and both
signed responses are ignored. Changing the count-nine core location in
the literal example is precisely one of the residual ambiguities left
by this decoder; freely moving singleton marks at fixed positions is not.

This provides a finite parametrization of the possible digits for a
given observed $(h,Q,S)$, not a uniform bound in $Q$, and says nothing
about the effect of a nonordinal recontact.
