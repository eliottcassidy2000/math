# Prime slopes require a different construction

**Proved.** In a tournament obtained by adding one apex to a transitive
spine, fix a binary contact word and let up to two of its spine positions
grow into constant-letter blocks. Their widths may be arbitrary fixed
positive affine functions of a common integer parameter. Such a family
cannot give an exact ray $1+p2^{ak+b}$ with $p$ an odd prime and $a>0$.
In particular this method cannot uniformly close the outstanding
defect-14 ray $38\cdot2^j+1$, even on a progression of heights. This does
not prohibit individual one-apex witnesses or changing the word itself.

Also certified: the full $46\cdot2^j+1$, $58\cdot2^j+1$, and
$74\cdot2^j+1$ rays, $j\ge0$. The first two close the entire defect-20
layer; the third closes one branch of defect 22. These witnesses come from marked quotients that
need not have the one-apex structure.

## One growing block: an exact divisor criterion

Write spine vertices in forward order. A bit 1 denotes an arc from
the apex to that spine vertex. For a nonempty word $w$,
$$
H(w)=[w_1=1]+[w_\ell=0]+\#01_{\rm adjacent}
+\sum_{i<j,\ w_i=1,w_j=0}2^{j-i-1}.
$$
Indeed a path either starts or ends at the apex, or it enters from
$i$ with bit 0 and exits to $j$ with bit 1. With $i<j$ only consecutive
positions work; when $i>j$, each intermediate spine vertex can go on
either side. Equivalently define $F(u)=H(w)$ for the binary integer
$u$ represented by $w$, padding by leading zeros if desired. Then
$F(0)=1$, $F(2u+1)=F(u)$, and $F(2u)=F(u)+u+[u\text{ odd}]$ for $u>0$.

Consider $L0^mR$ with $m\ge1$. Let $A$ be the binary integer of $L$,
and $B$ the binary integer of the reverse-complement of $R$
(empty words give zero). If $A=0$ the count is independent of $m$.
Otherwise split the inversion sum at this block to obtain
$$
H(L0^mR) = A(B+1)2^m+F(A)+F(B)-A-[A\text{ even}].
$$
In detail, the pairs from $L$ into the zero block contribute
$A(2^m-1)$ and the pairs across it contribute $AB\,2^m$; the other
terms give $F(A)+F(B)-[A\text{ even}]$. This holds if $R$ is empty
as well, using $F(0)=1$. Converse and spine reversal handle a
growing block of ones.

Write $A=e2^a$, $B+1=f2^b$ with $e,f$ odd. The recurrence gives
$$
H=ef\,2^{m+a+b}+F(e)-e+F(f-1).
$$
For example, positive zero padding of odd $e$ adds
$e(2^a-1)+1$ before the parity correction, and appending ones
to $f-1$ leaves $F$ unchanged. Thus the odd parts $q$ of
positive slopes that occur with intercept 1 are exactly those with
an ordered factorization $q=ef$ satisfying $F(e)-e+F(f-1)=1$.
For sufficiency take $A=e$, $B=f-1$ and their associated finite words.

No odd prime $p$ passes: for $(e,f)=(1,p)$ the intercept is
$F(p-1)>1$; for $(p,1)$ it is $1+F(p)-p<1$.
For the latter, $F(u)\le u+1$ by induction (use the stronger
$F(2v+1)=F(v)$ when the half input is odd), hence
$F(p)\le(p+1)/2$. The inequality $F(p-1)>1$ follows from
the last-bit recurrence at the positive even integer $p-1$.
There are infinite positive cases: with $u=2^t$, $t\ge1$,
take $e=2u+1,f=u+1$. Here $F(u)=u+1$ and $F(e)=F(u)$;
the criterion gives the family $1+(2u+1)(u+1)2^m$.

## Two growing blocks: bit reversal blocks the apparent escape

For $Lx^sy^mR$ with a fixed intervening word $M$ between the two
blocks, put $X=2^s,Y=2^m$, $K=2^{|M|}$,
$d=\operatorname{bin}(M)$ and
$z=\operatorname{bin}(\operatorname{revcomp} M)$.
Let $A,B$ be as in the one-block formula, now from the outer words.
The count has form $C_0+UX+VY+WXY$. Its coefficients are:

| $(x,y)$ | $W$ | $U$ | $V$ |
|---|---|---|---|
| 00 | $AK(B+1)$ | $A(1+z-K)$ | $(B+1)d$ |
| 01 | $AKB$ | $A(1+z)$ | $B(d+1)$ |
| 10 | $(A+1)K(B+1)$ | $(A+1)(z-K)$ | $(B+1)(d-K)$ |
| 11 | converse of 00 | | |

Here the notation means $L\,x^s\,M\,y^m\,R$.
For derivation, scan left to right and keep the current binary
prefix value $v$. A zero block adds $(X-1)v$ to the inversion
sum and updates $v\mapsto Xv$. A one block updates
$v\mapsto X(v+1)-1$ without adding to the inversion sum.
Scanning $M$ adds $zv$ plus a fixed amount and updates
$v\mapsto Kv+d$. The endpoint and adjacent terms never depend
on the positive widths. This gives the table.

If $W=0$ the table reduces to one active block. Otherwise let widths
be $s=ak+b$, $m=ck+d_0$ with $a,c>0$. The leading base of the
result is $2^{a+c}$. With unequal speeds $a\ne c$, both linear
terms must vanish individually for a pure exponential plus a
constant. With equal speeds they may cancel, but 01 has strictly
positive linear terms and 10 strictly negative terms.
For 00, if $M$ contains no 1, merge the zero blocks; otherwise put
$e_M=K-1-z>0$. The alleged cancellation says
$$ A e_M\,2^b=(B+1)d\,2^{d_0}. $$
The binary words of $d$ and $e_M$ are reversals. After removing
powers of 2, both have the same bit length, that of $M$ after its
initial and terminal zero strings are removed. Their ratio is between
$1/2$ and $2$. If $W$ has odd part the prime $p$, the odd parts of
$A,B+1$ must be $1,p$ in some order, so that ratio would have to
be $p$ or $1/p$. Contradiction. The 11 case follows by converse.
The merged and inactive cases are covered by the one-block criterion.

This is an obstruction to the **form of the growth rule**. It is
not a count exclusion for arbitrary one-apex words: for example
$39$ already occurs at $w=101000$. Any all-height solution of
the $38$ ray via this kind of stabilization must change the core,
couple more than two growth blocks, or use a non-spine structure.

## Constructive comparison: two higher-defect rays

Use the following three quotients and replace the specified mark
with a transitive block of width $m$.

| Slope | Width $m$ | Outgoing row masks | Mark | Visits $(c_1,c_2)$ |
|---:|---:|---|---:|---|
| 46 | $j-3$ | 64,5,49,39,75,403,174,31,223 | 8 | (737,368) |
| 58 | $j+1$ | 254,124,120,48,224,192,8,78 | 6 | (59,29) |
| 74 | $j+1$ | 254,204,248,112,98,194,128,24 | 6 | (75,37) |

Each mark has in- or outdegree one, so at most two runs are
possible. A two-run word allows $2^m-2$ allocations into nonempty
increasing runs. Exact Hamilton subset DP at widths 1 and 2
certifies the profiles. The resulting counts are
$46\cdot2^j+1$ for $j\ge4$, and $58\cdot2^j+1$, $74\cdot2^j+1$ for
$j\ge0$. The four smaller values $47,93,185,369$ have direct arc-patch
witnesses in the replay.

The three words are $1000\,1^{j+1}$, $0010\,1^{j+1}$ and
$11010\,1^{j+1}$ respectively. The last has defect 22. For the
defect-20 classification, use the run costs from
[signed cancellation](SIGNED_CANCELLATION_DESCENT.md).
An initial two-zero run costs 6 and its least later cost is 14,
giving exactly $0010\,1^{j+1}$ in this case. An initial longer run
either costs 14, after which a further run exceeds the budget,
or already costs at least 30. A single initial zero followed by
one later run costs $2+3\cdot2^s$ for a single later zero, or
at least 22 for a double. Three runs there exceed 20. With no
initial zeros, a lone later run costs $2^s$ or $3\cdot2^s+2$ for
length one or two; neither is 20. The length-three minimum is
exactly 20, uniquely after one initial 1. A longer run exceeds 20.
Two later runs either have the minimal costs $2+10=12$ or cost
at least 22, and a third run costs still more. Thus precisely the
46 and 58 families constitute defect 20.

**Scope.** Defect 14 still has the unresolved $38$ family in this
argument; there is no claim of a full Chan–Pak proof. The point of
the prime-slope theorem is to change the representation used against
that ray, rather than trying longer fixed one-apex templates
of the same form. Run python3 prime_slope_replay.py.
