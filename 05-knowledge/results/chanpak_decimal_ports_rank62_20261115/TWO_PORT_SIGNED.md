# A physical two-deletion-port unipotent counter

There is one explicit strong tournament $T_k$ of order $16+10k$ with two distinguished vertices $a,b$ such that
$$H(T_k-b)-H(T_k-a)=2k\qquad(k\ge0).$$
For $k\ge1$, $T_{k+1}$ is obtained from $T_k$ by a fixed ten-vertex insertion retaining every old arc. Thus the signed integer readout
$$(u,1)=\bigl((H(T_k-b)-H(T_k-a))/2,1\bigr)$$
is transported by the genuine unipotent matrix $\begin{pmatrix}1&1\\0&1\end{pmatrix}$. This does not claim unipotent growth of either individual path count, or of the full two-apex count. It is a simultaneous measurement on two physical deletion ports of the **same** tournament.

## Literal construction and nesting

Start with a transitive spine. A 1 on an apex's mark word means that it beats the corresponding spine vertex. Give the two apices the mark words
```
A_k = 1101010010 (0001010101)^k 0110
B_k = 1010       (0000000001)^k 0110011100
```
and give the mutual arc either orientation, say $a\to b$. Both words have length $14+10k$ and begin with 1 and end with 0. Thus both apices reach the first spine vertex, the last spine vertex reaches both, and $T_k$ is strong. The same argument proves that each deletion is strong.

For each $k\ge1$ insert ten new transitive spine sites just after site 10. Their marks in $A$ are $0001010101$, the exact period of $A$. In $B$, site 11 is seven positions into a period (zero-based phase six), so their marks are $0001000000$, that period rotated by six slots. The cut is internal to the periodic region in $B$ for $k\ge1$. Inserting a full period in each stream changes neither subsequent phase nor any old mark or spine arc. This is a stationary action on ordered tournaments with two fixed apex labels and a fixed insertion cut. The formulas also give $u=0$ for the stand-alone $k=0$ tournament; nesting is claimed only for $k\ge1$.

## Proof that the exponential parts cancel exactly

For a painted transitive spine the [inherited singleton cut state](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/chanpak_selective_target_20261114/CUT_STATE.md) has $h=1$. On appending one bit $\mu$,
$$(N,S,d)\mapsto(N+(1-\mu)(S+d),\,2S+\mu,\,\mu).$$
Start with $(1,0,0)$. For any nonempty bit word $W$ let $q=2^{|W|}$, let $s$ be its usual binary value, let $t$ be the reversed complemented binary value, and let $n$ be its $N$ output starting from $(0,0,0)$. Iterating the displayed singleton rule gives
$$(N,S,d)\mapsto(N+tS+(1-W_1)d+n,\ qS+s,\ W_{|W|}).$$

Take a prefix of binary value $P$ ending in 0, a repeated tile $W$ beginning in 0 and ending in 1, and a suffix beginning in 0 with reversed complemented value $U$. After $k$ tiles, $S_k=(P+s/(q-1))q^k-s/(q-1)$. The input $d$ for the first tile is zero and every later $d$ is one; the single missing term is returned by reading the suffix. Summing the resulting geometric progression, including $k=0$, gives
$$H_k=H_0+(P+s/(q-1))(U+t/(q-1))(q^k-1)+k(n+1-ts/(q-1)).$$

In the two rows the data $(P,s,t,n,U)$ are $(850,85,343,30,9)$ and $(10,1,511,0,793)$, with $q=1024$ in both. In fact
$$850\cdot 1023+85=85(10\cdot1023+1),\qquad 793\cdot1023+511=85(9\cdot1023+343),$$
so the exponential coefficients agree. The suffix and prefix choices also make $H_0=8259$ in both. The linear drifts are $2558/1023$ and $512/1023$, whose difference is exactly two.

This is not a restatement of the matched-core slide: each single-apex deletion has an entirely transitive spine and the two marked probes are different; the construction does not claim identical full continuation profiles. Its virtue is a **fixed physical insertion producing the exact signed signal** inside one object; its limitation is that extracting that signed signal as a single positive Hamiltonian count would require an additional gadget outside append-only growth, by the [sharp no-go bound](UNIPOTENT_OBSTRUCTION.md).

The separately assembled partition replay checks the literal two deletion probes through $k=1$, as well as the exact state identities and adjacency nesting for longer controls.

