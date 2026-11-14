# A fixed-size return with exact ternary growth

A ready marked tournament has \((H,h,L,E)=(H,e,H,e)\): \(h\) counts paths after deleting the root, \(L\) counts the relaxed rooted split covers and \(E\) counts paths ending at the root. Call the ready pair **primitive** if \(\gcd(H,e)=1\). A terminal reroot returns to this form. The moves below add thirteen vertices, return a primitive pair, and can be reversed from that pair alone. A finite two-state rule gives exact ternary growth of the number of ready profiles.

For a near-transitive attachment \(R_j\), use \(O_j\) for the old root, \(I_{j,k}\) for the \(k\)-th interior chain vertex, and \(T_j\) for the new terminal. Each step adds \(j+1\) vertices. Write \(P_j=I_{j,j}\) for the penultimate vertex. Use the corrected [interior-port transfer](FERMAT_INTERIOR_PORTS.md): in the first-interior case the response subtracts \(E\), including at \(I_{2,1}\). For \(P_j\) with \(j\ge2\), on the pair \((H,L)\) only,
\[
P_j:(H,L)\mapsto(H+dL,2H+dL),\qquad d=2^j.
\]
The terminal step acts into a ready pair by \(T_j:(H,L)\mapsto(H+dL,H)\). Both \(2\times2\) matrices have determinant \(-d\). The first coordinate stays odd; over every odd-prime field the matrices are invertible, so they preserve primitivity.

## One small prefix is always safe

From a ready input \((H,e)\) choose between two ways to use the same three-vertex attachment:

| prefix | \((H_{\rm pre},L_{\rm pre})\) | when it loses primitivity |
|---|---|---|
| \(O_2\) | \((5H,8H-e)\) | \(e\equiv3H\pmod5\) |
| \(I_{2,1}\) | \((5H,10H-e)\) | \(e\equiv0\pmod5\) |

Each determinant is \(-5\). The two bad lines intersect only in the zero vector modulo 5, impossible for a primitive input. Allow either prefix whenever it is safe. A deterministic subset may choose \(O_2\) unless it is bad. These prefixes' output ranges of \(L_{\rm pre}/H_{\rm pre}\) are respectively \([7/5,8/5)\) and \([9/5,2)\), because \(0<e\le H\). Thus the prefix choice is detectable from its two outputs.

Now append **either** \(P_2,P_4\) **or** \(P_4,P_2\) in that order, followed by \(T_1\). The two penultimate moves cost eight vertices in either order, and the terminal move costs two, for a total of thirteen. The only exponents used are the true Fermat exponents \(1,2,4\).

## The choice survives the return

Write \(r=L/H\ge1\) at a physical root. A penultimate step sends this ratio to
\[
1+\frac{1}{1+d r}.
\]
Last letter \(P_4\) makes \(r\le18/17\). Last letter \(P_2\) following \(P_4\) gives \(r\ge1+1/(1+4\cdot18/17)=106/89>18/17\). The terminal reading \(e'/H'=1/(1+2r)\) therefore identifies which order was used. Invert the terminal and penultimate matrices; the prefix ratio distinguishes the two prefixes, each itself invertible. Thus the output determines the input and the chosen prefix and suffix. In particular different words of fixed length give distinct ready profiles.

## Exact growth by three, with only two core states

Call a ready profile *free* if both prefixes are safe, and *forced* if exactly one is safe. After either allowed prefix \((H_{\rm pre},L_{\rm pre})\equiv(0,u)\pmod5\) with \(u\ne0\). The two penultimate orders give \((3u,2u)\) and \((0,u)\) modulo 5, respectively. After \(T_1\) these become \((2u,3u)\), whose \(e/H\) ratio is \(4\) mod 5 (both prefixes safe), or \((2u,0)\) (only the old prefix safe). The transition count matrix, with rows and columns free/forced, is therefore
\[
\begin{pmatrix}2&2\\1&1\end{pmatrix}.
\]
If the starting pair has \(p\in\{1,2\}\) safe prefixes, the exact number of reached ready profiles at depth \(m\ge1\) is \(2p\,3^{m-1}\), all at order \(n+13m\). The graphs are strong; at least \(2p\,3^{m-1}/(n+13m)\) unrooted isomorphism classes occur since one graph supplies at most its order many rooted profiles. This count is exact for the named grammar, not for tournaments at that order.

### Keep the count prime to three at a small size cost

Instead of fixed \(T_1\), choose \(T_1\) or \(T_4\) so \(H+dL\not\equiv0\pmod3\). Because \(d=2,16\equiv-1,1\pmod3\), at least one works: both failing would imply \(H\equiv L\equiv0\), contrary to primitivity. Allow one deterministic choice, preferring \(T_1\) when both work. The mod-5 growth rule does not change: \(T_4\) sends the two penultimate outputs \((3u,2u),(0,u)\) to \((0,3u),(u,0)\), which are still free, forced, respectively. The decoding remains injective. Before the terminal \(1\le r\le6/5\) in either suffix, and the intervals \(1/(1+2r)\) and \(1/(1+16r)\) are disjoint. The indicated terminal can be read back, then the suffix and prefix as before. Thus the \(2p\,3^{m-1}\) distinct ready profiles have 3-unit counts and lie in only \(m+1\) possible orders \(n+13m+3i\) for \(0\le i\le m\). One such order contains at least \(2p\,3^{m-1}/(m+1)\) profiles. This is a physical phase guard, not primality.

## Scalar outputs, and the price of this particular grammar

With the fixed \(T_1\) the two scalar child readouts are \(2201H-212e\) and \(2321H-272e\) under the old prefix, or \(2625H-212e\) and \(2865H-272e\) under the interior prefix. Each is strictly changing in \(e\) for fixed \(H\). Let \(D_m\) count distinct scalar values and \(S_m=2p\,3^{m-1}\) ready profiles at depth \(m\ge1\). Some \(H\)-fiber has at least \(S_m/D_m\) states. Assign a deterministic safe prefix per state again; at least half use the same prefix, and a fixed suffix for that half gives distinct counts at the next depth. Hence \(D_mD_{m+1}\ge S_m/2\). In particular there are exponentially many distinct strong counts across two consecutive size levels of this dictionary.

This is an actual return, not just a compatible unrooted count or a congruence representative. It still branches slowly compared with the count growth: every step in it costs at least a factor three, and there are four physical seam steps per return. Therefore a single seed in this exact return language cannot cover an interval of all large odds by itself. The target-selection problem remains. The [literal replay](replay_rooted_returns.py) checks all corrected ports on small graphs, the final readouts from a seed, and controls for both possible safe prefixes.

## A complete primitive 7-adic address space in the fixed-size family

Start at the literal order-four ready pair \((H,e)=(5,1)\), supplied by \(R_2\). Keep only the old prefix and \(T_1\), but retain both penultimate orders. This is allowed indefinitely: modulo 5 the returned ratios are \(4\) and \(0\), both old-prefix safe. Modulo 3 the two words have the same map \((H,e)\mapsto(2(H-e),H+e)\), so \(e/H\) alternates \(2\leftrightarrow0\) from the given seed. Every count reached by a nonempty return word is prime to 15, and every pair stays primitive.

The exact matrices (acting on column pairs) for the two return words are
\[
M_0=\begin{pmatrix}2201&-212\\709&-68\end{pmatrix},\qquad
M_1=\begin{pmatrix}2321&-272\\685&-80\end{pmatrix};
\]
both have determinant 640, prime to 7. The replay checks these three finite field data:
* the \(M_0\) orbit of \(x_0=(5,1)\) is all 48 nonzero vectors modulo 7;
* \(M_0^{48}\equiv I+7A\pmod{49}\) with \(A=\begin{pmatrix}6&2\\5&0\end{pmatrix}\pmod7\);
* putting \(B=M_1 A M_1^{-1}\pmod7\), the vectors \(Ax_0=(4,4), Bx_0=(0,1)\) are independent.

**Proved for every precision.** Positive words realize every pair modulo \(7^k\) that is not simultaneously zero modulo 7. For the inductive digit step, the power \(M_0^{48\cdot 7^{t-1}}\) acts by \(I+7^t A\) modulo \(7^{t+1}\), by binomial expansion. Its \(M_1\)-conjugate acts by \(I+7^t B\). At the residue \(M_0^i x_0\), conjugate both by \(M_0^i\); their two displacement vectors remain independent. Both corrections fix all lower digits and drive any desired next pair of digits, using powers 0 through 6. Every inverse used in a conjugation can be replaced by a positive power modulo the current modulus, because the two matrices are invertible over that finite ring. Hence the correcting words are literal positive physical returns, not inverse graph operations. This proves the lift.

One exact example is \(M_1(5,1)=(11333,3345)\) with \(11333=7\cdot1619\) and 1619 prime. Repeating \(M_0^{336}\) first fixes the pair modulo 49, then \(M_1\) gives \(H\equiv14\pmod{49}\) at arbitrarily large orders. This proves an unbounded exact valuation-one lane and complete local pair control, not infinitely many prime cofactors or an absolute integer selection theorem.

## Fixed-level branching with both 99 phase guards

The same seed and fixed old prefix also admit a branch language keeping \(H\) prime to \(3,5,11\) at every return, with no variable-size adjustment. When \(H\) is an 11-unit put \(r=e/H\pmod{11}\). The two matrices act on this ratio by
\[
0:\ r\mapsto (5+9r)/(1+8r),\qquad
1:\ r\mapsto (3+8r)/(3r).
\]
A step is legal for the guard when its denominator is nonzero. The initial ratio is 9; two steps \(11\) take it through 4 to 2. From 2, the word \(0000\) stays at 2, while \(1011\) goes \(2,5,9,4,2\). Every denominator in these words is a unit. Thus after the initial two steps, either four-step route can be repeated in any order, retaining both 99 guards, the existing modulo-5 guard, and a primitive ready state. By the decoding proof there are \(2^m\) distinct such profiles at the single order \(30+52m\). This transfers a usable phase rather than merely returning to a residue once.

### All addresses at one common order

The all-precision statement strengthens to **all sufficiently large orders of the form \(4+13m\)** for each fixed precision, in the same fixed-old-prefix sublanguage. The length-13 positive bit word \(0001011011111\) returns \(x_0\) modulo 49. Both matrices have order 336 modulo 49, by the finite replay.

Lift the endpoint of that short loop back to \(x_0\) at precision \(7^k\) using the two independent digit corrections above, now only at \(t\ge2\). The \(M_0\) pulse has length \(48\cdot7^{t-1}\), a multiple of 336. To implement the conjugated pulse, use a positive \(M_1\) step and its inverse represented by one less than its finite order at precision \(7^{t+1}\); together they add a multiple of 336. Thus the corrected loop at \(x_0\) has length \(a\equiv13\pmod{336}\). A second positive loop is the \(M_0\) identity word of length \(b=48\cdot7^{k-1}\). For \(k\ge2\) these lengths are coprime. Prefix any word reaching a prescribed pair with nonnegative combinations of the two loops. Coprime lengths generate every sufficiently large integer: choose the coefficient of \(a\) in \(0,\dots,b-1\) to achieve the desired residue modulo \(b\), then fill with \(b\). There are finitely many target residue pairs, so a single depth bound works for all of them. Project from precision 49 for \(k=1\).

One can take the bound \(m=O(49^k)\). The first lifted loop length is \(O(7^k)\) since a digit step uses bounded repetitions of pulses of length \(O(7^t)\); the finite-order conjugation overhead has the same bound by binomial lifting from modulo 49. A shortest reaching word has length smaller than the \(48\cdot49^{k-1}\) possible residue states, by cycle erasure. Add the coprime-loop padding cost \(O(ab)\). This is a common-size theorem for local addresses, **not** common-size exact coverage of the integer line.
