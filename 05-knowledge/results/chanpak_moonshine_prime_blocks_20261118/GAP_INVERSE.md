# Inverting an unmarked gap between two marked flanks

**Scope:** exact inverse for indefinitely growing gaps, an effective order cap when the inverse fails, and an all-order exclusion for the two nearest-singleton contacts used by the 347/349 witness. An integer-weighted collision is possible, so there is a genuine distinction between the signed lattice and binary physical contacts. No binary collision or all-prime selector is claimed.

## 1. Exact three-read transfer

Take nonempty flanks $A,B$ and gap $Q$ in ordinal order $A\to Q\to B$. Add a vertex $v$ whose out-neighbors form fixed sets in $A,B$ and do not meet $Q$. Switch just $v$. Write $f_T(q)$ for the optional-slot polynomial, $\mathcal E_q(T),\mathcal S_q(T)$ for the marked end/start totals on a flank, and $F_j=f_Q(j)$. For each $q=0,1,2$ define
$$
 a_q=\mathcal E_q(A) f_B(q)-\mathcal S_q(B)f_A(q),\qquad
 b_q=\mathcal E_q(B)f_A(q+1)-\mathcal S_q(A)f_B(q+1).
$$
The unsigned marked end-minus-start reading on the entire exterior is $a_q F_q+b_q F_{q+1}$. A path ending in $A$ must leave $Q,B$ wholly outside that path; a path ending in $B$ can include an extra slot from each earlier block. The starts are reversed. This proves the reading by the ordinal product law. Since $Q$ is nonempty, $F_0=0$. By the [adjacent-read identity](PADDING_OBSTRUCTION.md),
$$
\begin{aligned}
 \Delta H &= (a_1+b_0)F_1+b_1 F_2,\\
 \Delta B &= 2\big(a_1 F_1+(a_2+b_1)F_2+b_2 F_3\big).
\end{aligned}
$$
The replay checks these formulas with whole-host subset DP on 100 literal marked triples.

Let $c=(c_1,c_2,c_3)=(a_1,a_2+b_1,b_2)$. If $c_3\ne0$, a neutral gap of order $n$ must satisfy
$$
 |c_3|(n+2)/2\le |c_2|+|c_1|/(n+1),
$$
by splitting slots to get $F_3/F_2\ge(n+2)/2$ and $F_2/F_1\ge n+1$. This is an effective order cap; the left side increases with $n$ and the right decreases. If $c_3=0$ and $c_2\ne0$, then $n+1\le|c_1/c_2|$. If both vanish and $c_1\ne0$, neutrality is impossible. Separately, a nonzero $c$ allows at most two distinct transitive gap widths, because $c_1+c_2 2^n+c_3 3^n$ has at most two real zeros (its derivative has at most one unless the function is constant). Thus unbounded neutral gap order is equivalent to $c=0$, in which case **every gap topology is neutral**. The problem of sustaining a switch reduces to this exact equality, rather than searching through widening gaps.

## 2. A three-number inverse and two hidden reads

Set $B=C^{\mathrm{op}}$ with the same marking as $C$ under reversal. Then $c=0$ is exactly equality between the marked flanks $A,C$ on the three-number key
$$
\left(
 \frac{\mathcal E_1}{H},\
 \frac{\mathcal S_2}{f(3)},\
 \frac{\mathcal E_2-\mathcal S_1}{f(2)}
\right).
$$
The normalized coefficients of $F_1,F_2$ in the remaining $\Delta H$ are differences on the two coordinates
$$
\left(\frac{\mathcal S_0}{H},\ \frac{\mathcal S_1}{f(2)}\right).
$$
More precisely, since $a_1=0$, they are respectively $b_0/(H(A)H(C))$ and $b_1/(f_A(2)f_C(2))$, with right-minus-left differences. To get a genuine switch at unbounded order, some hidden coordinate must differ. If the second differs, the magnitude grows with $F_2$ along arbitrary sufficiently large gaps. If only the first differs, it grows along gaps with $H(Q)\to\infty$, for instance ordinal words of directed triangles; it does not grow on a transitive gap.

Appending any ordinal unmarked tail to either normalized flank leaves all five fractions invariant: the end total acquires the tail factor at $q$, the start total at $q+1$, and the denominator has exactly the same factor. Thus one can trim past the last marked strong component before looking for collisions. For a wholly transitive flank the first fraction is just the binary contact address $\sum_{i\in I}2^i$. It determines all five readings after trimming the trailing zeros, so two transitive flanks cannot give a nontrivial collision.

The background search reported no qualifying collision through flank order eight, over all its pointed classes. Two additional discovery scans here found none among equal-slot-profile order-nine classes, or with two different markings on the same randomly sampled order 9–12 orientation. These runs are partial beyond order eight, and are not needed for any all-order result.

**Signed-lattice versus binary contacts.** On the literal five-vertex flank with rows $(20,5,24,19,2)$, take signed contact weights $(-9,-5,8,2,4)$. Direct evaluation gives difference zero on all three unnormalized first reads, but changes the two unnormalized hidden reads $(\mathcal S_0,\mathcal S_1)$ by $(12,15)$. Thus no linear identity forces the hidden coordinates to agree. This is **not** a physical collision: a difference of two actual binary contact masks has entries only in $\{-1,0,1\}$. The constructive question is to lower these multiplicities to binary ones or to implement them with a genuine graph gadget without creating unwanted slot terms.

## 3. The published prime-bank escape has exactly one gap

Let the left flank be $A=L\to a$ and the right flank $B=b\to R$, with $a,b$ singleton vertices and the sole out-contacts of $v$ at $a,b$. Keep $Q$ between them. Write
$$
 u=f_L(2)H(R)-f_R(2)H(L),\qquad
 w=f_L(3)f_R(2)-f_R(3)f_L(2).
$$
Then $c=(u,2u+2w,3w)$ and $\Delta H=u(F_1+2F_2)$. Indeed the marked end and start on $L\to a$ are $f_L(q+1)$ and $f_L(q)$, respectively; on $b\to R$ they are $f_R(q)$ and $f_R(q+1)$. Hence infinitely many neutral gaps force $u=w=0$, which also makes $\Delta H=0$. Putting the two contacts just at the inside singleton endpoints cannot sustain a nonzero switch, even with an arbitrary growing gap.

For the [347/349 matched pad witness](PRIME_SUPPORT.md), the replay gives $u=-134784$, $w=51840$, and $-u/w=13/5$. A neutral gap must satisfy
$$
\frac{2F_2+3F_3}{F_1+2F_2}=\frac{13}{5}.
$$
For order $n\ge1$ the ratio on the left is at least
$$
\frac{2+3(n+2)/2}{\,2+1/(n+1)\,}.
$$
This bound is strictly increasing and equals $13/5$ at $n=1$. The one-vertex gap attains equality; **no other nonempty gap of any topology** works with these flanks. The earlier two prime witnesses remain valid. This specifically rules out varying that unmarked gap while keeping the marked flanks fixed.

Replay: **python3 gap_replay.py** in this packet. The remaining positive obligation is a physical collision of the three-number key that changes a hidden read, rather than another attempt to grow the singleton-contact seed.
