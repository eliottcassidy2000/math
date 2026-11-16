# Every residue from a single apex; exact children at 61 and 63

**New strengthening of the all-moduli theorem.** For every positive
modulus \(M\), every residue represented by an odd integer is the
number of Hamiltonian paths of a strong one-apex tournament at **every
sufficiently large order**. The constructed tournaments at an order
share a transitive deletion and the same number of reversed arcs in
their natural orders. The same construction can pin a leading
real-scale window and the residue simultaneously.

**Priority.** The underlying congruence completeness was already
proved by [Session 25's reset blocks](../session25_spectrum_wiring_20261021/session25-spectrum/PROOF.md),
not first proved here. The alternative marked-position inversion
proof below gives an explicit fixed-order statement and a construction
whose order bound grows like \(\sqrt m\), rather than \(m\), apart from
the dyadic and period parameters (\(m\) is the odd part of \(M\)).
The preceding [prime-palette result](../prime_palette_residue_tree_20261116/README.md)
remains useful for exact counts; its modular bound was not intrinsic.

## 1. A strong source and its arithmetic

Take the transitive spine \(1,\dots,n\). Apex \(0\) beats every spine vertex except at marked positions \(t\) in a set \(A\subseteq\{2,\dots,n\}\), and require \(n\in A\). The cycle \(0\to1\to\cdots\to n\to0\) makes it strong. The known [one-apex formula](../../../01-canon/theorems/THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md) is
\[
H=1+\sum_{t\in A}2^{t-1}
      -\sum_{s<t\in A}2^{\max(t-s-1,1)}. \tag{1}
\]
For a short derivation, the [Grinberg–Stanley odd-cycle expansion](https://arxiv.org/abs/2307.05569) gives \(H=1+2C\) here: all directed cycles pass through the apex, so no two can be disjoint. For a losing position \(t\), all smaller positions are allowed first contacts provisionally; the even subsets between that contact and \(t\) give \(2^{t-2}\) cycles in total before exclusions. Exclude each marked first contact \(s<t\), whose doubled contribution is \(2^{\max(t-s-1,1)}\). This proves (1). The companion code checks (1) by Hamiltonian-path subset DP.

Write \(M=2^am\) with \(m\) odd, put \(L=\operatorname{ord}_{m}(2)\) (set \(L=1\) when \(m=1\)), and call gaps at least \(D=\max(2,a+1)\) long. Start at or beyond \(D\) and separate gadgets by long gaps. Only the short internal pairs read at the dyadic coordinate. At the odd coordinate the contribution for every nonadjacent pair is \(2^{t-s-1}\) and depends only on the ordered phases \(s,t\bmod L\). A cross term between disjoint ordered gadgets depends only on each phase multiset: it is half the product of their sums of \(2^{-s}\) and \(2^t\). Adjacency \(t=s+1\) contributes one extra unit to the subtracted pair sum.

## 2. Three controls, all on the same spine

1. **Low binary bits.** For every bit \(i=1,\dots,a-1\), allocate a pair of marks. Give it gap \(i+1\) (activated) or a long gap congruent to \(i+1\bmod L\) (off). Both options have the same contribution modulo \(m\) and the same phase multiset. The activated option subtracts \(2^i\bmod 2^a\), and the off option subtracts zero. Every even dyadic correction can be read independently of the odd coordinate.
2. **One ternary carry, if \(3\mid m\).** Use two three-mark gadgets, each of phase multiset \(\{0,1,2\}\bmod L\). A gadget is either phases \((0,1,2)\) with gaps \(1,\text{long}\), or phases \((0,2,1)\) with gaps \(2,\text{long}\). Both subtract 2 at the dyadic coordinate. Modulo 3, in the formal no-adjacency sum every pair term is orientation-independent, since \(2^{t-s}=2^{s-t}\); the gap-1 option incurs one extra subtraction. The two independent switches therefore give offsets \(0,-1,-2\) modulo 3. Nothing here assumes \(3^a\) is the entire odd modulus.
3. **The remaining odd residue.** Take \(x\) tail marks of phase 0 and \(x\) of phase 1, all gaps long. Swapping an adjacent \(01\) to \(10\) changes the subtracted pair from \(1\) to \(1/4\bmod m\); other contributions are unchanged modulo \(m\). Thus \(H=H_0+(3/4)I\bmod m\), where \(I\) counts instances of a 1 before a 0. Every integer between 0 and \(x^2\) is an inversion count: start with all 0s before all 1s and bubble-sort to the reverse block, adding one each time. Take \(x=\lceil\sqrt{m/\gcd(m,3)}\rceil\) when \(m>1\) and none when \(m=1\). The step generates the multiples of \(\gcd(m,3)\); control 2 chooses the right coset if needed. The tail has no dyadic contribution.

Choose bit activations using the fixed dyadic baseline \(1\) or \(1-4\), depending on whether the ternary gadgets are omitted or present. Fix the terminal mark at the desired last spine position \(n\), far from all other marks. Include it when reading \(H_0\). It never changes any swap increment. Each long step can set its requested phase in a gap at most \(D+L-1\); each short step is smaller than \(D\). Consequently if \(B\) is the number of marks including the terminal one, every \(n\ge B(D+L)\) works, with the same \(B\) and common deleted spine at that order. This proves the theorem. The constructor takes \((M,r,n)\), and places marks without enumerating tournaments or searching over residues.

The verifier uses the pair sum (1) as a control for the linear-time-in-marks recurrence in the constructor. It replays all allowable residues for \(M\le150\), an extra fixed larger order in every case, and mixed high-prime-power shadows of both forbidden counts.

**Absolute boundary.** For a final losing mark \(n\), already the cycles starting \(0\to1\) and ending \(n\to0\) give \(H\ge1+2^{n-2}\) for \(n\ge3\). The long gaps genuinely make big counts, not small remainders hidden in another coordinate. No bound of one modulus on the absolute error follows. In particular this proves density in the congruence topology even in the strong one-apex family, not density in the ordinary integer sense. In fact [Session 26](../session26_random_bridges_20261022/session26-tournament/APEX_CLASSIFICATION_AND_SPECTRAL_LAW.md)
proves that the missing one-apex counts have lower relative density
at least \(1-\log 2\) among odd integers.

**Combined large-scale and residue placement.** For fixed \(M\)
the same controls also give \(H\equiv r\bmod M\) with
\(|H-N|=O_M(\sqrt N)\) for arbitrary sufficiently large real-scale
integer targets \(N\) and allowable residues \(r\). Take
\(n=1+\lfloor\log_2N\rfloor\) and \(t_0=\lfloor n/2\rfloor\).
Reserve a band of constant length \(C(M)\) above \(t_0\) for all
controls, followed by a long gap; set \(u=t_0+C(M)\).
At positions \(t\ge u+1\) mark precisely the 1-bits \(t-1\) of \(N\).
This suffix includes the terminal position \(n\); its singleton sum is
the truncation of \(N\) at bit \(u\). All control choices fit below the
gap (enlarge \(C(M)\) using the preceding bound). Read the known
internal suffix correction into the dyadic baseline, and include
the suffix in the odd baseline. The three control increments are
unchanged: all cross terms see only phase multisets on their control
side and vanish dyadically. Control singleton sums and the tail of
\(N\) below bit \(u\) are \(O(2^u)\). The entire subtraction in (1)
is \(O(2^{n-t_0}+n)\): bound the formal weights by
\(\frac12(\sum_{s\ge t_0}2^{-s})(\sum_{t\le n}2^t)\), plus at most
one adjacency correction per mark. Both errors are \(O_M(\sqrt N)\).
Here one can fix \(r=N\bmod M\), but the hidden constant depends on
\(M\); this is not a uniform absolute repair or a variable-modulus
proof of equality.

## 3. Literal dyadic children: two paths from 31, but only one induced

A three-loser fan has marks \(a,b,L_0=a+b\), with \(2\le a<b\). From (1), cancellation gives
\[
H=2^{L_0-1}+1-2^{\max(b-a-1,1)}.
\]
Insert an apex-beaten spine vertex between \(a,b\). The latter mark and \(L_0\) move right by one. A gap \(b-a=1\) gives \(2H+1\); gaps at least 2 give \(2H-1\). The latter renews indefinitely. In particular the \(\{2,4,6\}\) fan at 31 extends to marks \(\{2,5,7\}\) at 61, then to 121, 241, and so on, without changing any old arcs.

The balanced Mersenne fan \(F_k\) uses \(a=\lfloor k/2\rfloor\), \(b=k+1-a\), \(L_0=k+1\), and has \(H=2^k-1\). At even \(k\), the same induced insertion gives \(F_{k+1}\). At odd \(k\ge5\), put the new apex-beaten vertex **after** \(b\) and before the final mark, then reverse one old spine arc, between the loser at \(a\) and its immediate winning successor. This swaps their locations in the transitive spine, producing exactly \(F_{k+1}\). The same source thus gives either sign at every odd \(k\), with zero old arc changes for the minus sign and one for the plus sign. There is also an entirely induced two-vertex move \(F_k\subset F_{k+2}\): insert a winner before \(a\) and another between \(b\) and \(L_0\). It realizes \(H\mapsto4H+3\).

**Zero arc changes are impossible for \(31\to63\), even with arbitrarily many new vertices.** An \(H=63\) tournament has only one nontrivial strong component: any nontrivial factorization of \(63\) needs 7 or 21, which are [globally excluded](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md). The [Busch bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3) makes this component order at most eight (\(5^{8/3}>63\)). Pad it to order eight with ordinal singletons. The exhaustive [THM-344 census](../../../01-canon/theorems/THM-344-n8-h63-complete-omega.md) gives exactly the above \(F_6\) with marks \(\{3,4,7\}\) and its distinct converse, identified by their score profiles. Direct checking of every proper induced subgraph gives, at orders six and seven respectively,
\[
\{1,5,11,15,17,19,23\},\quad \{1,11,25,33,37\}.
\]
Smaller subsets are also checked in the script and give no 31. Converse has the same induced spectrum. Exterior singletons do not affect the count of an induced subset, so the obstruction holds at every order. One old arc change at this dyadic edge is therefore sharp; choosing a different parent works too: \(33\to63\) is a literal \(2H-3\) extension (delete spine position 3 from \(F_6\) to get 33). The count-only binary parent was too restrictive.

## 4. Both signs, no old edits, for infinitely many parents

Write a fan as a binary spine word, with 1 for a vertex beaten by the apex. For \(p=2^d\), \(d\ge1\), set
\[
A=1\,0^d\,1,\qquad B=1^d0.
\]
The word \(A\,000\,B\) has \(S=(4p+3)^2\) paths. Insert a loss in the central triple to get \(A\,0000\,B\) with \(2S-1\), or insert a win to get \(A\,0100\,B\) with \(2S+1\). The first instance is \(121\to241,243\). It is on the very same induced chain \(31\to61\to121\) above.

One quick computation: for a word \(W\), put \(X(W)=\sum W_j2^{|W|-j}\) and \(Z(W)=\sum(1-W_j)2^{j-1}\). Counting odd cycles by first and last spine contact gives \(H(W)=1+\#10_{\rm adjacent}+\sum_{i<j,W_i=1,W_j=0}2^{j-i-1}\); hence concatenation has cross term \(X(U)Z(V)+[U_{\rm last}=1,V_{\rm first}=0]\) in addition to \(H(U)+H(V)-1\). Use \(H(A)=H(B)=p+1\), \(X(A)=2p+1\), \(Z(B)=p\). For a middle word \(I\) beginning in 0 this yields
\[
H(AIB)=2p+H(I)+1+(2p+1)Z(I)+p\big(2^{|I|}(2p+1)+X(I)\big).
\]
The triples \((H,Z,X)\) for \(000,0000,0100\) are \((1,7,0),(1,15,0),(5,13,4)\). More generally the left middle \(0^{3+j}\) gives \((S-1)2^j+1\); the right middle \(0\,1^j\,00\) gives \(S2^j+1\) for \(j\ge1\). Both branches can continue by pure insertions, but this does not prove they can branch again at every generation.

## 5. The \(3+7=10\) microcosm and genuinely different meanings of 7

Specialize the archived [divisor braid](../../../05-knowledge/results/arithmetic_braids_20260917_divisors.md) to \(3^2\cdot7\cdot11=693\). Its 10 proper nontrivial divisors have seven squarefree entries, of which three are primes. The equality \(F=U+S=3+7\) is not a disjoint union: the three repeated prime entries compensate for the three nonsquarefree entries \(\{9,63,99\}\). In the column \(7,21,63\), the first two entries are the only absent exact counts in the whole divisor box. The second 3-step takes that column out of the squarefree face and into the compensating triple. The corresponding \(11,33,99\) column is complete. So the formal swap \(7\leftrightarrow11\) of the divisor box is **not** a symmetry of exact realization. The exact split is 8 realized and 2 absent; the archived 3+7 bookkeeping has a different meaning.

The same \(3+7=10\) also has a literal tournament interpretation.
In the coordinate \(t=(H-1)/2\), multiplication by 3 sends
\(t\mapsto3t+1=t+H\): the three addresses for \(7,21,63\) are
\(3,10,31\), giving \(3+7=10\) and \(10+21=31\).
For a one-apex tournament \(t\) is the number of vertices in the
complete odd-cycle conflict graph. Applying the modular theorem to
\(H\equiv2t+1\bmod2M\) shows that these graph orders also cover every
congruence class, even though complete orders 3 and 10 are absent
and 30 and 31 occur. This ties the two charts together at the same
integers without identifying the divisor-box counts with graph orders.

There is a supporting but family-specific explanation on the Mersenne ray. The orders of 2 modulo 3 and 7 are 2 and 3, respectively. Thus both divide \(2^k-1\) first at \(k=6\), where the count is 63, and whenever they both divide it, so does 9. LTE gives \(v_3=1+v_3(k)\) at even \(k\) and \(v_7=1+v_7(k)\) when \(3\mid k\). The 7-ray reads \(7(1+8+\cdots+8^{k/3-1})\) for \(3\mid k\): 7 is missing, but 63 and \(511=7\cdot73\) are realized. This does **not** explain 21's all-order exclusion by itself: \(105\), which has \(v_3=v_7=1\), does occur in other families.

In fact 511 is the **only** count on this ray of the form \(7p\) with
\(p\) prime. Write \(k=3m\). If \(3\nmid m\) and \(m>1\), then
\((2^k-1)/7=(2^m-1)(1+2^m+2^{2m})/7\) splits into two factors greater than
one, with the 7 in the second factor. If \(3\mid m\) and \(m>3\), the
quotient is a proper multiple of 73 because \(2^9-1\mid2^k-1\).
The other two cases are \(m=1\), with quotient 1, and \(m=3\), with
quotient 73. This is a concrete limit to extending a prime palette using
one renewable family.

Separately, literal witnesses at 497, 511, 553 and 581, together with
469 and the primes 67, 71, 73, 79, 83 in the preceding packet, enlarge
the exact palette to 83. The multiplicative proof for smooth targets is
unchanged. This is a finite seed extension, not an all-prime theorem.

The same micro-to-macro argument turns **any** verified prefix
through \(B\ge343\) into unbounded exact coverage for all odd integers
whose prime factors are at most \(\lfloor B/7\rfloor\), except
7 and 21. The factor of 7 guarantees the seed \(7p\) as well as \(p\).
There is no division by 7 needed away from \(v_7=1\): if \(v_7=0\) or
\(v_7\ge2\), every odd count with other primes at most \(B\) is covered
using the ordinary prime seeds and (when needed) 49 and 343.
For instance, the [corrected independent 32-bit prefix
certificate](../../../session203/coverage-recertification/README.md)
gives an exact smooth palette through 613,566,756, without certifying
any new primes past that original finite prefix.

Ordinal decomposition into strong components is canonical for a tournament; the resulting product on counts forgets the components. On the proved smooth palette the multiplicative atoms are the ordinary odd primes other than 7, plus \(49,63,343\) and \(7p\) for palette primes \(p\ne3,7\). Think of the ordinary and exceptional atoms as an analogy to infinite families and sporadic cases, not a classification of simple groups: there is no Jordan–Hölder uniqueness of arithmetic atoms (e.g. \(35^2=49\cdot5^2\)). Also arithmetic 7 is not sevenfold symmetry. A 7-order tournament automorphism has a 7-cycle orbit inducing a circulant. The eight direction choices have counts \(175\) (six) and \(189\) (two). Insertion monotonicity gives \(H\ge175\) in any tournament admitting such symmetry, so none of 35, 49 or 63 can have it. The Paley incidence symmetry of the Fano plane and tournament arc-preserving symmetry must not be conflated.

There is a stronger distinction for prime counts. Automorphisms act
freely on ordered Hamilton paths, so their order divides \(H\).
For a prime count \(p\ge5\), the only nontrivial strong component has
order at most \(1+3\log_5 p<p\), by Busch. It cannot admit an
element of order \(p\), and the canonical order of components fixes
singleton padding. Thus every prime-count tournament at \(p\ge5\)
is rigid. The count 3 is different: its triangle has automorphism
group \(C_3\).

**Remaining:** select exact sources at arbitrary targets or prove a physical \(\pm3\) or editable-child rule on a sufficient class of parents. The universal modular result eliminates any purely congruential obstruction to that selection, but does not make the selection for us.

Replay with `python verify.py`.
