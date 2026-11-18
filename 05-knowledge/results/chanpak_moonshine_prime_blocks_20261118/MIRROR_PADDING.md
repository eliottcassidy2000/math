# Growing a neutral cut requires growing contact complexity

Two restrictions on the [prime-support route](PRIME_SUPPORT.md): (i) fixed pad shapes cannot be made to work at indefinitely many heights by padding each independently with ordinal singletons, even when the mirror-marked exterior is allowed to change topology; (ii) if the exterior is entirely transitive, any unbounded sequence of nonmirror neutral cuts must have an unbounded number of mismatched mirror pairs. These are all-order results for the indicated architectures, not an obstruction to general tournament construction.

## 1. A dual obstruction: fixing the pad shapes and varying the exterior

Use the [two-pad transfer](PADDING_OBSTRUCTION.md) for a nonempty marked core $U,I$ with equal marked end and start totals at readings 1 and 2. This includes every self-converse core under an anti-isomorphism preserving $I$, not just a transitive order with mirror contacts. Denote the equal totals by $\mathcal E_1,\mathcal E_2$ and put $K=\mathcal E_2/\mathcal E_1$.

For fixed nonempty pad shapes $L,R$, let $x=L_1/L_2$, $y=R_1/R_2$, $z=L_3/L_2$, $w=R_3/R_2$. Add $a,b\ge0$ ordinal singletons to the left and right pad shapes respectively, writing $d=b-a$. Ordinal insertion of singletons multiplies $f(q)$ by $q$ per singleton, independently of their position in the ordinal list. Dividing the neutrality equation by the two-slot pad readings gives, for a nonzero neutral move,
$$
K=3^{-a}\frac{x-y2^{-d}}{z-w(3/2)^d}.
$$
The numerator is nonzero by nonzero shift, and the denominator is therefore nonzero in a solution. Positivity forces $d$ strictly between $\log_2(y/x)$ and $\log_{3/2}(z/w)$, regardless of their order. For each of these finitely many integer $d$, the factor $3^{-a}$ and $K\ge1$ bound $a$; then $b$ is fixed by $d$.

In fact a **fixed ratio** $K$ bounds the core order even for nontransitive cores. For order $m$, the marked slot totals obey
$$
\frac{\mathcal E_{q+1}}{\mathcal E_q}\ge\frac{m+q}{q+1}.
$$
An input is a marked nonempty path ending in $I$, with $q$ ordinary possibly empty slots. Split any ordinary slot, or split off a proper prefix of the marked path and place that prefix in a new first ordinary slot. An input supplies $m+q$ choices: each ordinary slot contributes its length plus one, and the marked path contributes its length, allowing an empty prefix but not an empty marked suffix. An output has at most $q+1$ preimages: at most $q$ merges of adjacent ordinary slots, plus one way to prepend the first ordinary path to the marked path. Thus $K\ge(m+1)/2$. We obtain finite orders, hence finitely many marked cores and full switch realizations for these fixed shapes. This generalizes the transitive-exterior observation from the background pass.

For the published 347/349 **unpadded** order-eight pad shapes, the exact inverse on transitive mirror-marked cores gives $(a,b,m,I)=(0,0,7,20),(1,1,5,10),(2,2,3,5)$ with $I$ a contact bit mask. All three have order 24 including $v$; they are the same labeled figure when surplus singletons are placed on the inner sides of the pads, by passing the unused ends of the middle transitive block to the pads. **Placement qualifier:** if the right pad's extra singletons are instead placed outside its nontransitive part, as in the first publication, the three displayed pre-switch tournaments have different outdegree multisets. They are not all one isomorphism class. The three readings stay the same because ordinal placement within a pad does not change its slot counts. There is no new height under either placement convention.

## 2. A wholly transitive exterior: an exact signed-digit equation

Let $m$ be the exterior size, $n=m-1$, and label its vertices $0,\dots,n$ in transitive order. Put $b_i=1$ if $v\to i$, and let $Q(z)=\sum_{i=0}^n (b_i-b_{n-i})z^i$. The coefficient difference is anti-reciprocal and vanishes at a midpoint. For a nonempty marked path ending at $i$, any chosen predecessors can be in that path. Thus the unsigned end-minus-start read at general $q$ is $D_q=q^n Q((q+1)/q)$ for positive integers $q$. The adjacent-read identity gives the **exact** neutrality criterion
$$
 2^n Q(3/2)+Q(2)=0.
$$
The full ordinary change is $Q(2)+b_n-b_0$ (the last difference is the full-support rooted read $D_0$). A nonmirror contact mask automatically has nonzero ordinary change: the top nonzero term of the signed binary sum $Q(2)$ dominates the lower terms, and $b_n-b_0$ is either zero or has that same sign. So the displayed equation for a nonmirror mask is exactly the missing neutral move.

Take one index $i<n/2$ per mismatched pair and put $j=n-i$. The neutral equation is equivalently a signed sum of positive weights
$$
 V_i=2^i(3^j-1)-2^j(3^i-1)
$$
with digits $-1,0,1$. There is no solution with at most three mismatched pairs. Indeed, $V_i>0$ and $V_{i+1}<(2/3)V_i$ wherever both indices are before the midpoint; directly $2V_i-3V_{i+1}=4\cdot2^i+(5\cdot3^i+1)2^{j-1}>0$. Two weights cannot cancel. With three, only three consecutive indices can have the two smaller weights add to the largest: otherwise their ratio sum is at most $2/3+8/27<1$. But for $j-i\ge5$ the exact value of $(V_i-V_{i+1}-V_{i+2})/2^i$ is
$$
5-3^{j-2}+2^{j-i-2}(11\cdot3^i+1).
$$
This cannot vanish modulo 8, since $3^{j-2}$ is 1 or 3 modulo 8. More generally the first two mismatch indices must be at distance at most 2 (a gap of 3 makes even the whole remaining geometric tail less than the first weight), and if all mismatch indices are at pairwise distance at least 2 the sum cannot vanish (the tail is at most $4/5$ of the first weight).

### All-order consequence: the number of mismatches must grow

**Theorem.** For every fixed $s$, only finitely many exterior orders admit a nonmirror neutral mask with at most $s$ mismatched pairs.

Suppose otherwise, and choose a sequence with $n$ tending to infinity and a fixed positive count of mismatches. Write each used weight as $A_i-B_i+C_i-D_i$, where
$$
 A_i=2^i3^{j},\quad B_i=2^{j}3^i,\quad C_i=2^{j},\quad D_i=2^i,\qquad j=n-i.
$$
Retain its signed coefficient in the neutral equation. Partition the resulting signed terms into minimal vanishing blocks. There are finitely many partitions of the fixed list of term positions (number the mismatches in increasing order), and finitely many sign patterns. Pass to a subsequence with those fixed. Within each block, divide by one term. There are only finitely many normalized tuples by the [Evertse–Schlickewei–Schmidt theorem for finite-rank groups](https://annals.math.princeton.edu/2002/155-3/p04), applied to the positive $\{2,3\}$-units with fixed signs in the coefficients; a two-term block has ratio 1 directly. Pass to a further subsequence on which all term ratios within blocks are fixed.

A block containing any $A_i$ cannot contain a pure power of two, since $v_3(A_i)=j\ge n/2$ and a fixed ratio would require that exponent bounded. It also cannot consist solely of $A$ terms: their 3-adic valuations are all distinct, so the term of least valuation cannot cancel. Thus each $A_i$ must share its block with some $B_k$. Their fixed ratio makes $j_i-i_k$ bounded. But
$$
j_i-i_k=\tfrac12\big[(n-2i)+(n-2k)\big]
$$
is a sum of nonnegative deficits; each participating $i$ and $k$ lies a bounded distance from the midpoint. In particular **every** mismatch index $i$ lies a bounded distance from the midpoint, as every one contributes an $A_i$. Now both $A_i$ and $B_i$ have 3-exponents tending to infinity. The pure powers can therefore be grouped only with pure powers in the fixed vanishing blocks. This is impossible: their exponents are all distinct (one on each side of the midpoint per pair), and the largest signed binary term strictly dominates the sum of all smaller powers. The contradiction proves the theorem. No effective upper bound on the order for fixed $s$ is asserted.

## 3. Replay, finite scope, and the remaining route

Run **python3 mirror_padding_replay.py** alongside the earlier replays. It checks the exact three decompositions and their inner and outer placement distinctions. It checks the literal transitive contact equations by subset DP through exterior order 7 and uses exact meet-in-the-middle on the signed weights to exclude nonmirror neutral masks through exterior order 48. The finite check does not rule out larger complicated masks.

This removes fixed pad shapes and bounded mismatch complexity as routes to indefinite neutral switches in their respective models. The previous 347/349 witnesses remain valid. New pad shapes, changing marked asymmetric exteriors, or a genuinely growing and highly interleaved contact mask remain possible routes. None of these reductions places an arbitrary prime at the required height.
