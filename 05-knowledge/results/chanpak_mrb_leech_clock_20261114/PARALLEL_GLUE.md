# Crossing matched ports: one exact inverse, and a second saturation obstruction

**Scope.** The [exact corridor decoder](https://github.com/eliottcassidy2000/math/blob/edf094e6bc983f5e8c8599aecdd2e0c00783773c/05-knowledge/results/chanpak_frame_decoder_pressure_20261114/RESETS_AND_DECODER.md) compresses two neighboring counting states only at seams of length 0 or 13. A question prompted by the Leech-lattice glue analogy is whether overlapping several nontrivial connections, rather than factoring at a reset, supplies the missing counts. Here there is an actual tournament construction and an exact obstruction: not even arbitrarily many parallel matched ports, nested arbitrarily and interspersed with ordinal sums from the terrace, yield enough scalar counts. This does not apply to ports with a freely chosen matching, nonuniform connections to the middle, or general overlapping gadgets.

## Exchange independent factors for interacting ports

For a nonempty tournament $B$, put $h$ left terminals $a_1,\ldots,a_h$ before $B$ and $h$ right terminals $b_1,\ldots,b_h$ after it. The terminals on either side are transitive in increasing index; all arcs between regions are forward except $b_i\to a_i$, $1\le i\le h$. Call the result $P_h(B)$. Thus $P_1$ is the ordinary wrap. For each $h\ge2$, $P_h(B)$ is neither a base terrace nor an ordinal sum nor a single outer wrap:

* It is strong. From $a_1$ one reaches $b_1$ through any vertex of $B$ and all other vertices directly. Every left terminal can reach $a_1$ via $B,b_1$; the remaining vertices can reach those left terminals.
* Among $n+2h$ vertices, where $n=|B|$, the unique outdegree $n+2h-2$ is at $a_1$, and the unique outdegree 1 is at $b_h$. A wrap would require $b_h\to a_1$, but $a_1\to b_h$.
* Its triangle count is $C_3(B)+h(n+h-1)$. A triangle spanning several regions must contain one backward arc; the reversed arcs have disjoint ends, so it cannot contain two. For each backward arc its intermediate vertices, and only those, complete the triangle. This gives excess $(h-1)(n+h-2)+C_3(B)>0$ over the terrace minimum.

For an integer $m\ge0$, use the $m$-color count
$$
F_m(B)=\sum_{c:V(B)\to[m]}\ \prod_{i=1}^m H(B[c^{-1}(i)]) .
$$
The empty-set Hamilton count is one. $F_m$ is multiplicative for ordinal sums, increasing with $m\ge0$, and at least $m^{|B|}$ for $m\ge1$ by tournament path existence. Write $H=F_1$. Then the exact two-port formula is
$$
\boxed{H(P_2(B))=2F_3(B)+4F_2(B)-F_1(B).}
$$
It uses a negative lower channel but every final path is counted positively.

### Why crossing ports are an honest signed-color computation

The [Grinberg–Stanley odd-cycle packing identity](https://arxiv.org/abs/2307.05569) says that the path count is the sum of $2^{\#\mathcal C}$ over disjoint directed odd-cycle packings $\mathcal C$. Applying it in each color gives the extended polynomial identity, for slot weights $t_i\in\{0,1,-1\}$,
$$
\sum_{c:V(B)\to[m]}\prod_{i=1}^m
        t_i^{|c^{-1}(i)|} H(B[c^{-1}(i)])
=
F_{\sum_i t_i}(B).
$$
On the right, $F_x(B)$ is the polynomial
$$
F_x(B)=\sum_{\mathcal C} 2^{\#\mathcal C}
                 x^{|B|-|V(\mathcal C)|+\#\mathcal C}.
$$
Indeed a cycle must be in a single color and has odd length, so its slot contribution is $\sum_i t_i$; each unused vertex gives the same sum separately. Use $0^0=1$. The exponents have parity $|B|$, giving $F_{-j}=(-1)^{|B|}F_j$. In particular $F_0(B)=0$ for $B$ nonempty.

To evaluate $H(P_h(B))$, fix a packing skeleton outside $B$. Every cycle using both a terminal and $B$ runs through one or more contiguous $B$-segments of the cycle. Record each such nonempty segment as a labeled slot with its adjacent terminals and the directed order of the terminals; one can anchor each cycle at its least terminal to eliminate rotations. Distinct slots consume disjoint subsets of $B$. There are at most $h$ slots: between leaving one $B$-run and reentering $B$ the cyclic terminal itinerary must use some $b_i\to a_i$, since no arc goes from the right terminal side to $B$, or from $B$ to the left terminal side. All runs in all cycles require distinct such arcs.

As the arcs connecting a terminal to $B$ are uniform, a slot on a given subset has exactly its induced Hamiltonian count as weight. The leftover color, unconsumed by slots, carries the entire odd-cycle packing inside $B$. For each fixed external cycle, oddness is one congruence on the sum of its slot sizes, shifted by its terminal count. Expand these parity tests into half the sum or difference of sign characters; impose each nonempty-slot guard by inclusion-exclusion (weight zero forces that slot empty). Each term now has the signed-color form above with the leftover color assigned sign +1. Thus at any fixed parity of nonzero $|B|$, the whole ordinary count is a linear combination of $F_1,\ldots,F_{h+1}$, with coefficients depending only on $h$ and that parity.

The coefficient of $F_{h+1}$ is exactly $h!$ at either parity. To reach a color sum of absolute value $h+1$ requires using all $h$ slots with every sign positive. Each left and right terminal is then used exactly once as the endpoint of a $B$-run, and all the arcs $b_i\to a_i$ are used to close the cycles. Thus the only terminal itineraries contributing at the top are the permutations $a_i\to B_i\to b_{\pi(i)}\to a_{\pi(i)}$. A permutation gives coefficient one there: each cycle's factor 2 cancels the denominator of its parity test. This is a precise parallel-glue interaction, not the product factorization at the 13-zero corridor.

For $h=2$, the finite linear coefficient computation can either enumerate the port skeletons or use three independent probes per parity now that the span is proved. Take $B=TT_n$, so $F_j=j^n$. Literal subset end-DP for $n=2,4,6$ gives $33,225,1713$ and for $n=3,5,7$ gives $85,613,4885$, respectively. Each 3-by-3 matrix of $(1,2^n,3^n)$ rows is invertible by Vandermonde in the bases 1,4,9. Both produce coefficients $(-1,4,2)$. The empty case, if allowed, is a four-terminal count of 5 and agrees separately. Those probes are not a guess by interpolation without a span theorem.

## Two-port arithmetic and an exact prescribed-target check

For $B$ nonempty, signed-color parity yields
$$
H(P_2(B))\equiv
\begin{cases}
0 \pmod3 & |B|\text{ even},\\
H(B)\pmod3 & |B|\text{ odd}.
\end{cases}
$$
Indeed $F_3\equiv F_0=0$ and $F_2\equiv F_{-1}=(-1)^{|B|}F_1$ modulo 3. Hence prime or $7p$ targets above 21 that are not divisible by 3 require odd reservoir size. The smallest transitive middle gives 13; the transitive middle of size 5 gives 613.

More generally fix $B_0$ and append a transitive block of $r$ vertices within the middle, ordinally. The target must satisfy
$$
N=2F_3(B_0)3^r+4F_2(B_0)2^r-F_1(B_0).
$$
Every monomial of the color polynomial for nonempty $B_0$ has positive degree and nonnegative coefficient, so $F_3(B_0)\ge\frac32F_2(B_0)$. It follows that
$$
3^r\le\frac{N+F_1(B_0)}{2F_3(B_0)}<3^{r+1}.
$$
Thus at most one padding length needs an exact test: get $r$ by the floor of the base-three logarithm of that rational (use exact comparisons), reject $r<0$, and check the displayed count.

For a literal marked reservoir on five vertices take out-row masks $(18,8,3,5,14)$, bit $j$ meaning an arc to vertex $j$ with indices starting at zero. Its $(F_1,F_2,F_3)$ is $(13,104,471)$. The rail becomes
$$
942\cdot3^r+416\cdot2^r-13.
$$
For prime-cofactor $7p$ with $p\ne3,7$, $r$ must be even and modulo 7 it must be $r\equiv2\bmod 6$. The $r=2$ result is $10129=7\cdot1447$, with 1447 prime. This is an exact family of candidates with a target inverse and one prime cofactor, **not** a theorem that the rail contains infinitely many prime cofactors. It also does not claim the count itself previously lacked a witness.

## Arbitrary width and arbitrary nesting of this glue is still sparse

**Theorem.** Start with singletons and the strong minimum-triangle terrace. Close under arbitrary finite ordinal sums and under $P_h$ for every $h\ge1$, nesting arbitrarily. For all distinct ordinary counts $N>1$ in this closure,
$$
\sum_N N^{-119/125}<\infty.
$$
So the number at most $X$ is $o(X^{119/125})$. In particular this closure cannot provide either all but finitely many eligible primes or all but finitely many $7p$: their reciprocal sums already diverge.

Here is why the unbounded width and depth do not break the estimate. Write $s=119/125$. An atom description is a base seed or $P_h$ of an ordered nonempty sequence of atoms. Descriptions may coincide as tournaments; they overcount, which is safe. Let $U(z)$ be their ordinary order generating function and let $D(z)=U(z)/(1-U(z))$ count nonempty ordered reservoir descriptions. The bases have series $A(z)=z+z^3+z^4/(1-2z)$ from the terrace bit-word classification, so the positive grammar satisfies
$$
U(z)=A(z)+\frac{z^2}{1-z^2}\frac{U(z)}{1-U(z)}.
$$
At $z=7/25$ the right side, evaluated with $U=2/5$, is strictly less than $2/5$. At $z=11/50$ the right side at $U=7/25$ is less than $7/25$. Monotone iteration by maximum nesting depth proves $D(7/25)\le2/3$ and $D(11/50)\le7/18$.

The weighted $F_3$ count is appreciably smaller than this order-only grammar at base 3 because a new glue invokes larger color channels:
* $F_m(P_1(B))\ge mF_{m+1}(B)$ for $m\ge1$, from the [single-wrap color recurrence](https://github.com/eliottcassidy2000/math/blob/edf094e6bc983f5e8c8599aecdd2e0c00783773c/05-knowledge/results/chanpak_frame_decoder_pressure_20261114/WRAP_CHANNELS.md).
* $F_m(P_2(B))\ge2mF_{m+2}(B)$. In the $m$-color sum restrict all four terminals to the same color $i$. Split the part of $B$ in color $i$ into its three, two, or one distinguished subcolors according to the proven two-port formula. For each $i$ the result is $2F_{m+2}+4F_{m+1}-F_m\ge2F_{m+2}$. This is valid even if the original color has no middle vertices, as the two-port formula is valid there too.
* For $h\ge2$, $F_m(P_h(B))\ge m^{2h-4} F_m(P_2(B))$: delete all but two matched pairs, then use color-by-color induced path count monotonicity for each of the $m^{2h-4}$ extensions of a coloring. Monotonicity follows even by taking the same odd packings in the larger color class.

Consequently the $F_3^{-s}$ sum over **all atom descriptions** is at most
$$
A(3^{-s})+3^{-s}D(4^{-s})+
       \frac{6^{-s}}{1-3^{-2s}}D(5^{-s}).
$$
Exact rational upper estimates replace $3^{-s},4^{-s},5^{-s},6^{-s}$ by $9/25,7/25,11/50,19/100$, respectively. The result is bounded by
$$
A(9/25)+(9/25)(2/3)+
       \frac{19/100}{1-(9/25)^2}(7/18)
       =\frac{3390927427}{4284000000}<1.
$$
Thus $F_3^{-s}$ has a finite sum even over **ordered** ordinal reservoir descriptions, by multiplicativity and the geometric series on atom length.

Next the sum of $F_2^{-s}$ over atom descriptions is finite. For base seeds, use the earlier $0.93$ [single-seed $F_2$ estimate](https://github.com/eliottcassidy2000/math/blob/edf094e6bc983f5e8c8599aecdd2e0c00783773c/05-knowledge/results/chanpak_frame_decoder_pressure_20261114/WRAP_CHANNELS.md), which is explicitly over terrace words and includes finite smaller seeds separately. For $P_1$ atoms use $F_2(P_1(B))\ge2F_3(B)$ and the ordered reservoir sum just proved. For $h\ge2$ use $F_2(P_h(B))\ge2^{2h-2}F_4(B)\ge2^{2h-2}4^{|B|}$, with the convergent geometric sum over $h$ and $D(4^{-s})<\infty$.

Now control the **ordinary** outputs. For $P_1(B)$, the count $F_1(B)+F_2(B)$ is determined by the multiset of its constituent atoms, since each coordinate multiplies. Its inverse-power sum over distinct outputs is no more than the Euler product over atom descriptions with weights $F_2^{-s}$; this converges since each $F_2\ge2$ and their sum is finite. For $h\ge2$, keep the first two matching pairs together with $B$, and group the remaining pairs into disjoint doublets. Taking disjoint packs within these groups gives
$$
H(P_h(B))\ge H(P_2(B))\,5^{\lfloor(h-2)/2\rfloor}
           \ge 2F_3(B)\,5^{\lfloor(h-2)/2\rfloor}.
$$
Each four-terminal doublet has count 5 and the packings need not cross groups. So the $H^{-s}$ sum over all $P_h$ descriptions with $h\ge2$ converges using the ordered reservoir $F_3$ sum. The base distinct ordinary counts $>1$ have a finite sum at $s$ by the [terrace pressure theorem](https://github.com/eliottcassidy2000/math/blob/edf094e6bc983f5e8c8599aecdd2e0c00783773c/05-knowledge/results/chanpak_frame_decoder_pressure_20261114/PRESSURE.md). A final Euler product over the distinct atom counts $>1$ handles unrestricted ordinal products. To get little-o, split the summable series at a fixed count $M$, bound the count of outputs in $(M,X]$ by $X^s$ times the tail, and then send $M$ to infinity.

The conclusion is stronger than a fixed-width or one-wrap barrier: the same obstruction survives infinitely many possible widths and unbounded joint channel dimension. To make selection truly stronger, one needs a different kind of glue, for example nonuniform interaction with the middle or overlapping reservoirs whose intervening edges cannot be described by any parallel matching. Neither the numerical appearance of 13 nor the Leech analogy establishes such a construction.

Run `python cross_replay.py` for the rational grammar bounds, literal DP, congruences, and the exact width test.
