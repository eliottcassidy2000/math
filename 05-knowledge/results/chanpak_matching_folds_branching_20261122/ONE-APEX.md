# One-apex support bound: an endpoint-product law and a symmetry tax

**Status:** PROVED for the class of tournaments made transitive by deleting at most one vertex. Chan–Pak in general remains OPEN.

**New parts:** the joint real-scale and modular upper bound in Section 4, and the quantitative cutoff and residue law for the previously classified self-converse sub-spectrum in Section 5. Sections 1–3 record the inherited word formula and unconditioned density bound; see the repository's [binary-product density packet](../feedback_word_hamilton_spectrum_20261020/BINARY_PRODUCT_DENSITY.md). The unconditioned fixed-modulus cancellation is also known, from Session 193 and THM-4886; the argument below isolates the new simultaneous consequence.

**Claim.** Let $\mathcal A$ be the set of $f=(h(T)-1)/2$ for tournaments $T$ with a vertex whose deletion is transitive. Then
$$
\#(\mathcal A\cap[2^k,2^{k+1}))\le (\log 2)2^k+O(2^{k/2}).
$$
In particular $\mathcal A$ has upper natural density at most $\log2$. The same relative upper density holds inside every fixed arithmetic progression of $f$, so imposing the one-apex condition excludes infinitely many odd targets in every such progression. The upper bound counts symmetry orbits, not distinct outputs, so it does not claim equality.

## 1–3. Inherited word model and real-scale law

The [binary-product density packet](../feedback_word_hamilton_spectrum_20261020/BINARY_PRODUCT_DENSITY.md) supplies the following inputs. We state them in the orientation needed for the modular argument.

1. Add an apex \(v\) to a forward transitive order and encode \(w_i=1\) for \(v\to i\). Removing leading zeros and trailing ones preserves \(h\). Every nontransitive count has a normalized word \(w_1=1,w_L=0\), and
   \[
   F(w)=(h(T_w)-1)/2=\sum_{i<j}w_i(1-w_j)W(j-i),
   \quad W(1)=1,\quad W(d)=2^{d-2}\;(d\geq2).
   \]
   In particular \(F(1^a0^b)=2^{a+b-1}-2^{a-1}-2^{b-1}+1\). A normalized word representing a value in \([2^k,2^{k+1})\), \(k\geq2\), has length \(L=k+2\) or \(k+3\).
2. Put \(A_L=\sum_i w_i2^{-(i-1)}\) and \(B_L=\sum_j(1-w_j)2^{-(L-j)}\). Then \(F/2^{L-3}=A_LB_L+O(L^22^{-L})\). Under uniform normalized words, the first \(\lfloor L/2\rfloor\) bits at each end are disjoint sets of bits. Thus by dyadic rounding the product is uniformly within \(O(2^{-L/2})\) of a coupled product \(Z=UV\) of independent uniform variables on \([1,2]\). The density of \(Z\) is \(\log z\) on \((1,2)\) and \(\log(4/z)\) on \((2,4)\).
3. Converse acts on words by \(w_i\mapsto 1-w_{L+1-i}\), preserves \(F\), and has \(O(2^{L/2})\) fixed normalized words. Counting its orbits at lengths \(k+3\) and \(k+2\), the leading number of candidates in \([2^k,t2^k)\), \(1\leq t\leq2\), is
   \[
   2^k\!\left[\Pr(Z<t)+\tfrac12\Pr(2\leq Z<2t)\right]
    =2^k(t-1)\log2.
   \]
   The equality follows by differentiating in \(t\), using \(\log t+\log(2/t)=\log2\). Bounded density and the fixed words give error \(O(2^{k/2})\), uniformly in \(t\). These are word-orbit counts, hence only upper bounds for distinct target values.

## 4. Local equidistribution, global holes

For every fixed integer $q\ge2$ and residue $c$, the relative upper density of $\mathcal A$ inside $c\pmod q$ is at most $\log 2$. More precisely, uniformly for $1\le t\le2$, the number of occupied values in the residue and the truncated shell $[2^k,t2^k)$ is at most
$$
\frac{\log 2}{q}(t-1)2^k+O_q(2^{(1-\eta_q)k})
$$
for some $\eta_q>0$.

**State machine.** When the existing binary prefix has integer reading $R$, appending bit $e$ changes $R$ to $2R+e$ and adds $(1-e)\lceil R/2\rceil$ to the score $F$. This follows by summing the weights from all earlier ones to the new position. For the score modulo $q$ the state needed is only $R\bmod 2q$. In particular, modulo $q$, the formula for the increment is independent of the integer representative chosen for $R$.

**Fourier cancellation.** Fix $\ell$ with $2^\ell\ge2q$. An $\ell$-bit word can reset any prescribed input state to $R=0\bmod 2q$: choose its integer reading congruent to $-2^\ell R\pmod {2q}$. From zero compare two continuations of the same length: $0^{a+b}$ and $1^a0^b$. Choose $a\ge2$ with $2^a\equiv 1\pmod {q_{\rm odd}}$, where $q_{\rm odd}$ is the odd part of $q$, and choose $b\ge\max(2,1+v_2(q))$. Both continuations end at state zero, since $(2^a-1)2^b$ is divisible by $2q$. Their score difference, by the sorted-word formula, is
$$
D=2^{a+b-1}-2^{a-1}-2^{b-1}+1.
$$
It is odd, and $D\equiv 1/2\pmod {q_{\rm odd}}$ when $q_{\rm odd}>1$, so $D$ is a unit modulo $q$.

To make the resulting mixing precise, let $J=\ell+a+b$ and fix a nontrivial character $\chi$ of $\mathbb{Z}/q\mathbb{Z}$. The transfer matrix for $J$ random bits has row $R$, column $R'$ equal to the expectation of $\chi$(score increment) times the indicator of the transition from $R$ to $R'$. Every row has two particular equally likely words (the same reset and the two continuations), ending at the same state with phase ratio $\chi(D)\ne1$. Thus its absolute row sum is at most $1-2^{-J}(2-|1+\chi(D)|)<1$. The remaining bits contribute transfer operators of absolute row sum at most one. This proves exponential decay of the conditional Fourier coefficient whenever a contiguous block of free fair bits grows linearly with $L$, even if all bits before and after that block are fixed. This assertion follows by starting at the prefix state, multiplying by the block operators, and taking the final vector of unit phases associated with the fixed suffix. Fourier inversion gives uniform conditional residue probabilities $1/q+O_q(\rho_q^L)$ for some $\rho_q<1$, when the free block has length at least $L/2$.

**Separate local residues from scale.** Condition on the outer $\lfloor L/4\rfloor$ bits at each end, leaving the middle free. The corresponding truncated endpoint product is measurable under this conditioning. The uniform error from the complete normalized score is $O(2^{-L/4})$, and the independent truncated endpoint product lies in a boundary band of that width with probability $O(2^{-L/4})$ by the same dyadic coupling. Therefore each shell-membership event, including every truncated shell, has joint probability with $F\equiv c\pmod q$ equal to $1/q$ times its unconditioned probability plus $O_q(2^{-\eta_q L})$. Apply the same converse orbit count in both lengths. The fixed-word error is unchanged and smaller than the asserted error after reducing $\eta_q$ if necessary.

This proves exponential mixing at fixed local moduli but still a uniform fractional shortage of possible realizations within each progression. It does not estimate prime outputs: the constants and mixing rates depend on $q$, and fixed-modulus equidistribution does not amount to a sieve at the necessary scale.

### Growing-modulus version

The [short-voltage-loop bound](../chanpak_cyclic_mixing_20261115/README.md) applies to exactly the same address/checksum recurrence. After $r$ random bits the distribution of $(R,F)\pmod{(2q,q)}$ has total-variation distance from uniform at most $Cq\exp(-cr/\Lambda(q))$, where $\Lambda(q)=q^3\log^2(2q)$ (or $q\log^2(2q)$ if $q$ is a power of two). Under a fixed suffix, a conditionally uniform middle checksum stays uniform: the suffix only adds a deterministic function of the middle address to that checksum, and total variation cannot increase. The outer-quarter argument thus bounds the orbit-census error by
$$
O\!\left(2^k\left[2^{-k/4}+q\exp(-c'k/\Lambda(q))\right]\right).
$$
The same error works for a union of residues, using the union as the total-variation event, not summing errors. In particular, for $q\le k^\theta$ with fixed $\theta<1/3$ (or $\theta<1$ for powers of two), the main term $(\log 2)(t-1)2^k/q$ for an individual class has absolute error $O_A(2^k k^{-A})$ for each fixed $A$, uniformly in $t$. It is an asymptotic for presentation orbits and only an upper bound for distinct values.

For example, take $q$ to be the product of distinct odd primes with $q\le k^\theta$, and restrict to odd targets $h=2f+1$ coprime to $q$. Exactly $\varphi(q)$ residues of $f$ qualify; the full shell contains $2^k\varphi(q)/q+O(q)$ qualifying integers. By the same argument at least
$$
(1-\log 2+o(1))\,2^k\frac{\varphi(q)}q
$$
of these targets are absent from the one-apex class, uniformly in this growing-modulus range. This still does not isolate primes.

## 5. The blue diagonal is thin but locally uniform

The [existing blue one-apex classification](../chanpak_self_converse_core_20261118/encoder/ONE_APEX.md) gives exactly
$$
\{1\}\cup\{2^j+1:j\ge1\}\cup\{K(P):P\ge1\},\qquad
K(P)=P^2+2h(P)-1+(P\bmod2).
$$
Here $h(P)$ is the apex count of the binary word for $P$. The classification is on strong cores; source/sink padding cannot change a count and any converse isomorphism preserves the unique nontrivial strong component. Its proof already shows $P^2<K(P)\le(P+1)^2$, so $K$ is injective and replacing $K(P)\le X$ by $P\le\sqrt X$ changes at most one address. The dyadic exceptions add at most $O(\log X)$ other values. In particular the exact distinct-value cutoff is $\sqrt X+O(\log X)$.

**Local diagonal law.** For any fixed integer $q\ge1$ and odd residue $a\bmod 4q$, the number of distinct blue one-apex values up to $X$ in that residue is
$$
\frac{\sqrt X}{2q}+O_q\!\left(X^{(1-\eta_q)/2}\right)
$$
for some $\eta_q>0$. This is a law for distinct values, not merely presentations, by injectivity of $K$ and the sparse dyadic exceptions.

In fact $K(P)=P^2+4F(P)+1+(P\bmod2)$, where $F(P)$ is the score read by the same checksum machine. Its residue modulo $4q$ is determined by $P\bmod 2q$ and $F(P)\bmod q$. Half the address residues have the parity giving any requested odd value modulo 4, and for each such address exactly one checksum residue gives $a\bmod4q$. Thus the uniform product law on $(R,F)$ gives probability $1/(2q)$. Partition the initial address interval into aligned dyadic blocks (at most one at each power of two); on a block of size $2^r$ the high prefix is fixed and the last $r$ bits are fair. The short-voltage-loop bound in Section 4 therefore gives the result. More explicitly, the blocks of size below $Y^{1/2}$, where $Y=\lfloor\sqrt X\rfloor$, cost $O(\sqrt Y)$ addresses, and the others give total error $O(Y q\exp(-c''\log Y/\Lambda(q)))$. Including dyadic exceptions gives an additional $O(\log X)$. Hence the same uniform polylogarithmic range for $q$ is valid here as in Section 4.

Converse swaps the two boundary clocks. The two-dimensional bulk has the logarithmic product-volume law, while the fixed diagonal forces a square term and only square-root many targets, despite still reaching every fixed compatible residue. In both sectors the congruence laws can look completely unrestricted while the support is globally thin, and in the blue sector this is even true for distinct-value asymptotics. In particular no symmetric one-apex construction alone can cover the prime lane. These statements do not classify blue values outside one-apex tournaments.

## Interpretation and next obligation

A single feedback vertex might look like an unrestricted binary compiler. Its overlapping dyadic bands supply $\log 2$ presentation orbits per integer. Even when selected by congruence or small-prime sieving, there are not enough candidates, so a proof of Chan–Pak must use other families. Exact support density and nonsymmetric prime coverage remain open.

Source for the general problem: [Chan–Pak, Computational complexity of counting coincidences](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf), Conjecture 6.2, credited there to bof and Gordon Royle.
