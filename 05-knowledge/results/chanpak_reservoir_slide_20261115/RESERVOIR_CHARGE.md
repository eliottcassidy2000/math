# Expose a single internal apex: moving charges and shifted product decoding

The [previous reservoir admission test](../chanpak_window_descent_20261115/GROWING_WINDOWS.md) asks for a left ordinary count $g$, an apex count $n$, and the two signed start observables $s,d$. One way to choose all four physically, not independently, is to take a full member of an earlier one-apex palette as the new left spine, and expose only its old apex to the new apex. This yields a shifted product rather than a guessed endpoint vector.

## Exact wrapping law

For an ordinary painted spine $A$, write $a=H(A)$, $S,D$ for its total and signed weighted starts at marked vertices, and $g=H(A+u)$ for its count after attaching its original apex $u$. For a variable one-apex palette with fixed interface, $g$ varies while $(a,S,D)$ does not. Use $L=A+u$ as a **new** left spine, and mark only $u$ for the exterior apex. Then its actual state is

$$
(g,n,s,d)=\bigl(g,g+(a+S)-(a-D),a+S,a-D\bigr).
$$

A nonempty weighted path starting at $u$ is $u$ alone, of weight $a$, or $u$ prepended to a marked start in $A$, which reverses its parity. Hence $s=a+S,d=a-D$. The new apex beat set is the singleton $u$; the one-leaf identity $n=g+s-d$ follows by separating odd cycle packings involving the new apex. Every such cycle is the new apex followed by an even-length path (in vertex count) starting at $u$, and contributes twice its complement count. This is the same one-leaf law used in [the fan analysis](../chanpak_leaf_fans_20261115/FAN_CHANNELS.md).

If the right painted spine has ordinary count $b$, local apex count $B$ and ending data $t,f$, the [ordinal cut transfer](../chanpak_selective_target_20261114/CUT_STATE.md) yields

$$
N=gB+(b+t)(a+S)+(f-b)(a-D).
qquad (1)
$$

In particular independently variable left and right count palettes can be joined with both exposed responses held fixed. If the right spine is the converse of another copy of the left painted bank, with its marks complemented, reversal gives $(b,t,f)=(a,S,D)$. Consequently

$$
N=gB+(a+S)^2-(a-D)^2.
qquad (2)
$$

Here $g,B$ can be selected independently from two copies of the same fixed-interface palette. For a fixed charge $C$, target inversion is now a literal shifted factorization: find a divisor $g$ of $N-C$ in the left palette with $(N-C)/g$ in the right palette. Congruence alone is not sufficient. For example, if both banks contain only counts $>1$, a shifted prime $N-C$ defeats that charge altogether. This is a statement about one fixed bank, not a claim about missing primes from the tournament spectrum.

## Neutral tapes change the charge without changing the factors

Append a nonempty all-marked ordinal tail $T$ inside the left painted bank, with $c=H(T)$ and $Q=Q_2(T)$ its two-path-cover count. The original apex and earlier spine all dominate $T$. The ordinary painted spine has count $ac$, the full inner apex count becomes $cg$, the total start is $QS+a(Q-c)$, and the signed start becomes $ac$ (all endpoints in the last ordinary strong component are marked). Thus after exposing the original apex, its two weighted starts are $s'=Q(a+S)$, $d'=0$. This uses the generic fixed-profile transfer, not a property of a special tail shape.

Likewise one may append a marked tail inside the bank whose mirrored, complement-marked spine supplies the right block. If that original right bank has charge $b+S_B$, the two appended tails have $(c_L,Q_L)$ and $(c_R,Q_R)$, and at least one is nonempty, then

$$
\boxed{N=c_Lc_R\,gB+Q_LQ_R(a+S)(b+S_B).}
\qquad (3)
$$b+S_B$.}
qquad (3)
$$

For an empty tail take $c=Q=1$. With a nonempty left tail the leftover signed term vanishes via $d'=0$; a nonempty right tail instead makes $f'=b'$ and does the same job. An all-marked transitive tape of length $m$ has $c=1,Q=2^m$: the entire factor window survives unchanged, while its shift acquires a factor $2^m$. Equal symmetric taping gives $N=gB+[2^m(a+S)]^2$, $m\ge1$. Thus each prime dividing either factor must see $N$ as a square (allowing zero), a fast necessary filter, not a sufficiency test. Asymmetric tapes preserve the same factors but give $2^{m+n}(a+S)(b+S_B)$ for $m+n\ge1$, which is not generally a square.

Cyclic tags break the pure powers-of-two charge restriction. For example the two literal six-vertex ordinal tails with rows $(30,44,56,48,34,1)$ and $(22,44,56,49,34,1)$ both have $c=37$, while $Q$ is respectively $352$ and $376$. Choosing one replaces the charge without changing the factor scale or order; unlike a count-three tag, their two-cover counts need not be divisible by their path counts. For an exact $c$-modulus debt, $Q$ must be accounted for: a tag with $gcd(c,Q)>1$ forces that gcd into (3), independent of how the inner banks are chosen.

## A target-selected physical modulus 37

Use as the left inner bank just a singleton marked 1 (a second vertex, the internal apex, dominates it), followed by the $c=37,Q=352$ tag and $m\ge0$ marked sink singletons. This is a left reservoir of ordinary count $37$, with the exposed source carrying $s=704\,2^m,d=0$. It can vary $m$ from the prescribed target. Take the right mirror of the [unit slide of ordinary count 45](UNIT_SLIDE.md), with height $k$, exact apex-count interval $[U_k-2k,U_k]$ by twos, and $\eta_k=45+S_k$, where $S_k$ is the marked start response of that original slide. Formula (3) gives
$$
N=37(U_k-2j)+704\eta_k\,2^m,
quad 0\le j\le k,\quad m\ge0.
qquad (4)
$$
The combined tournament is strong: the exterior apex reaches the exposed left source, the left source reaches everything else on the spine, and the mirrored right bank has an unmarked final sink that returns to the exterior apex. A count-five cyclic tag $Q=32$ is another physical modulus in exactly the same way.

The explicit unit-slide start recurrence has
$\eta_kmod37=0,4,34$ for $k=0,1,2mod3$. Thus if $37\nmid N$ and $3\nmid k$, the residue equation in (4) has a unique $m_0mod36$: 2 generates the units modulo 37, because $2^{18}\equiv-1$ and $2^{12}\not\equiv1$. This services **every unit residue** at the chosen pump height, including both primitive lanes when the target is coprime to 37. It does not recenter the interval. For a prescribed $N$, restrict $m\equiv m_0mod36$ further by the actual integer inequalities

$$
\left\lceil\frac{N-37U_k}{704\eta_k}\right\rceil
 \le 2^m \le
\left\lfloor\frac{N-37(U_k-2k)}{704\eta_k}\right\rfloor .
$$

For an odd target, each surviving $m$ produces the literal slide position
$j=[37U_k+704\eta_k2^m-N]/74$; the residue, parity and bounds make it integral and admissible. At height 1, the prescribed target $11721300267$ returns $(m,j)=(0,1)$. The test is necessary and sufficient at fixed $k$ in this exact bank (for target units and $3\nmid k$). A fixed modulus 37 is a concrete prototype; allowing target-selected cyclic tags and two asymmetric charges learns which placement constraints are actually required but does not prove universal supply.

## Certified small controls in the primitive lanes

- With the $(37,352)$ tail at $m=0$, and a mirrored count-three spine (directed triangle with one mark before mirroring: $B=5,\eta=6$), the literal order-12 tournament has $N=37\cdot5+704\cdot6=4409$, prime.
- With the $5,32$ tail and a mirrored marked singleton ($B=1,\eta=2$), the literal order-eight tournament has $N=5+64\cdot2=133=7\cdot19$.
- The replay checks the unpadded mirrored charge with its subtractive square, symmetric and asymmetric transitive tapes, the two distinct cyclic tags of count 37, and a two-sided cyclic-tail assembly by full subset DP. These fixtures are concrete counts, not a primality-density claim.
