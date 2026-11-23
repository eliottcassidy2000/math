# Exact minimized residue automata for three strong apex tiles

Fix an odd modulus $q\ge3$. In a singleton apex word, the old vertices form a transitive order, and the letter at a vertex is $O$ if the apex beats it or $I$ if it beats the apex. A guarded word begins $O$ and ends $I$, so its tournament is strong. Use exactly three allowed macro tiles $C=OI,\ P=OOI,\ S=OII$. Every nonempty macro word is an actual strong tournament. The language for residue $r$ accepts precisely when its Hamilton path count is congruent to $r$ modulo $q$.

## Literal accepting Hamilton histories: a five-state sidecar verifier

On each old vertex, attach a sidecar bit $L$ or $R$ indicating which side of the apex it occupies in the path. Since the old tournament is transitive, each side necessarily occurs in increasing order. Thus the candidate associated with a bit string is uniquely “increasing $L$, apex, increasing $R$.” It is a Hamilton path iff the last $L$ is $I$ (or there is no $L$), and the first $R$ is $O$ (or there is no $R$). In particular accepted bit strings and actual Hamilton paths are in bijection.

This is a deterministic finite verifier on the paired alphabet consisting of an $O/I$ tile and one sidecar bit. It needs two Boolean registers: whether an $R$ has yet occurred, and whether the last $L$ was legal (initially true). If an $R$ with letter $I$ arrives before any other $R$, move to an absorbing failure state. An $L$ updates the last-$L$ legality. The final test is last-$L$ legality and no failure. The five states are all necessary for the unrestricted paired alphabet: all four live states and the failure state are reachable, good and bad live states differ on empty continuation, and two live states of the same legality having different $R$ registers are distinguished by $IR$ followed, if necessary, by $IL$. An $IL$ always repairs any live state but never the failure state.

It also gives an independent unsigned integer count for guarded words. The all-$L$ and all-$R$ histories both accept. Otherwise put $j=\max L$, $i=\min R$, so $W_j=I$ and $W_i=O$. If $j<i$, they must be adjacent, giving one history per adjacent $IO$. If $i<j$, all vertices before $i$ are forced left, all after $j$ are forced right, and the $j-i-1$ interior vertices are free. The total is
$$
2+\#\operatorname{adj}(IO)+\sum_{i<j:\,W_i=O,W_j=I}2^{j-i-1}.
$$
Since a guarded word begins with $O$ and ends with $I$, it has one more adjacent $OI$ than adjacent $IO$, giving the equivalent transfer form below.

## The tile representation

For guarded words let $d=2^{\text{number of singleton vertices}}$, $a$ be the positive start path read, $b$ the positive end path read, and $g=H-1$. The following representation is multiplicative on guarded concatenations:
$$
G(W)=\begin{pmatrix}1&a&g\\0&d&b\\0&0&1\end{pmatrix}.
$$
Indeed the unsigned two-ended transfer has middle coordinate $d$, and its $13$-entry counts each positive pair $O_i,I_j$ with weight $2^{j-i-1}$; the signed term is exactly the number of adjacent $OI$. This follows from projecting the selected subset between the pair onto even total path order (for an adjacent pair the empty internal subset occurs twice). The signed term adds across guarded concatenation because its junction is $IO$, not $OI$; the unsigned transfer multiplies by independent choice of which vertices between a pair are selected. Thus one may put the signed term in the $13$-entry without changing the multiplication on guarded words.

Work over $R=\mathbb Z/q$. Write
$$
O=\begin{pmatrix}1&1&0\\0&2&0\\0&0&1\end{pmatrix},\quad
I=\begin{pmatrix}1&0&0\\0&2&1\\0&0&1\end{pmatrix},\quad
Z(t)=\begin{pmatrix}1&0&t\\0&1&0\\0&0&1\end{pmatrix}.
$$
Here the single-letter matrices are unsigned. The three physical generators are $G(C)=Z(1)OI,\ G(P)=O G(C),\ G(S)=G(C)I$. Their monoid in the finite invertible matrix group is a group $\Gamma_q$, since positive powers produce every inverse. It therefore contains $O,I$, and $Z(1)=G(C)(OI)^{-1}$; conversely these three matrices generate it. Every group element, including the identity, has a nonempty macro-word realization.

## Exact syntactic group

Put $m=\operatorname{ord}_q(2)$. The diagonal coordinate $d$ has image $\langle2\rangle$ of order $m$. The unipotent quotient modulo its $g$-center is the additive pair $(a,b)$. In that quotient $OI^{-1}$ has pair $v=(1/2,-1)$, and its conjugate by $O$ has pair $w=(1/4,-2)$: conjugation scales the two coordinates respectively by $1/2$ and $2$. Their determinant is $-3/4$, and elementary elimination leaves a unit first pivot and a second pivot $-3/2$. Therefore their span is all $R^2$ if $3\nmid q$; when $3\mid q$ it is exactly the pairs with $a\equiv b\pmod 3$ (both generators lie in that congruence).

This last relation is also necessary for every matrix with $d=1$ when $3\mid q$. Reduce modulo three. There $O,I$ are involutions, $OI$ has order three, and its unipotent pair is $(2,2)$. Thus the group of $O,I$ modulo three is the nonabelian dihedral group of order six, whose diagonal-one subgroup is generated by $OI$ and has $a=b$. All central translations are independent for this condition. Since $Z(1)$ generates the entire $q$-element center, it follows that
$$
|\Gamma_q|=\operatorname{ord}_q(2)\,\frac{q^3}{\gcd(q,3)}.
$$
More explicitly, if $3\mid q$, diagonal-two matrices reduce to $a+b=1\pmod3$, by multiplying a diagonal-one matrix on the right by $O$; there are no further restrictions. Thus for every $(a,g)\in R^2$ some matrix with $d=1$ exists.

## Exact deterministic and two-sided minimization

Only $(a,g)$ is needed for rightward processing. A tile $B$ updates these coordinates by
$$
(a,g)\longmapsto(d_B a+a_B,\ g+g_B+a b_B),
$$
and the readout is $1+g$. The $q^2$ pairs are all reachable by nonempty physical words. Distinct states give distinct outputs already without continuation if their $g$'s differ, or after appending $C$ if only $a$ differs, since $b_C=2$ is a unit. For the Boolean language for any fixed $r$, append in addition a pure central shift $Z(t)$, realized by a physical word, sending exactly one of the two distinct outputs to residue $r$. Hence the minimal deterministic automaton has exactly $q^2$ states. This includes the empty word as the one-apex tournament; if one explicitly excludes the empty word in the Boolean language, the exceptional entrance adds one state exactly for $r=1$, and no states for other residues: nonempty suffixes already distinguish the original state from every other state.

Two-sided contexts distinguish every group matrix. The bare word measures $g$; right concatenation by $C$ gives $a$ and left concatenation gives $b$. Subtracting the two one-sided values from the two-sided $CWC$ value and adding the bare value isolates $a_C d b_C=4d$, up to the known copy constants which cancel in the mixed difference. Thus two distinct matrices have different $H$ in some physical left/right context, and the central shift turns this into Boolean separation exactly as above. Conversely equal matrices have identical contexts. Therefore the syntactic monoid (with the one-apex empty-word convention) is exactly the group $\Gamma_q$. Positive words suffice, and any fixed physical prefix and suffix may be imposed while realizing any reachable matrix, by solving their group equation for the middle word.

## A length obstruction invisible to the unrestricted language

For $3\mid q$, reduce modulo three. The unsigned $O,I$ group has order six and meets the central translations trivially (its center is trivial). Each physical tile is one central $Z(1)$ times an element of that unsigned group. In exactly $k$ macro letters, the reachable matrices modulo three therefore lie entirely in $Z(k)\langle O,I\rangle$. The unsigned right-state orbit of $(a,g)=(0,0)$ under $O,I$ is
$$
\{(0,0),(1,0),(2,1)\},
$$
by the elementary right update. Consequently
$$
g\equiv k+\mathbf1_{a\equiv2}\pmod3,\qquad H\not\equiv k\pmod3.
$$
This forbids one of the three residues at every fixed macro length despite full residue reachability with unbounded length.

The irreducible random walk on $\Gamma_q$ that appends any of the three tiles with positive probability has period precisely $\gcd(q,3)$. For an irreducible Cayley walk the period is the largest cyclic grading of the group sending every tile to $1$: define the grade by any word length to an element, using a word back to the identity to see independence modulo the period; multiplication then adds grades. Any such grading kills $O=G(P)G(C)^{-1}$ and $I=G(C)^{-1}G(S)$, so sends $Z(1)$ to $1$. The commutator of the two chosen unipotents is central $Z(-3/4)$, up to choice of commutator order; together with $Z(q)=1$, this forces the period to divide $\gcd(q,3)$. For $3\mid q$, the central grading in $\Gamma_3$ proves the lower bound. In particular the walk for $3\nmid q$ converges to uniform on the group, yielding uniform $H\bmod q$; if $3\mid q$ the three cyclic classes must be kept distinct. This is a statement about the positive guarded macros, not about all tournament counts or Hamilton counts as exact integers.

## This same complete modular family has a sparse, decidable integer spectrum

A nonempty macro word expands to $n\ge2$ singleton vertices. Its first $O$ and last $I$ alone contribute $2^{n-2}$ to the positive pair sum. Thus if $H=M>1$, necessarily $n\le\lfloor\log_2(M-1)\rfloor+2$. The number of macro words of expanded length $n$ has ordinary generating function $1/(1-x^2-2x^3)$. In particular it is $O(\lambda^n)$, where $1<\lambda<2$ is the positive root of $\lambda^3-\lambda-2=0$; alternatively this bound follows by induction from the recurrence by parts of weights 2 and 3. Enumerating just these lengths and updating the integer representation therefore decides exact membership in $O(M^{\log_2\lambda}\operatorname{polylog}M)$ bit operations. The number of distinct realized integers below $X$ is $O(X^{\log_2\lambda})$, with $\log_2\lambda<1$.

In particular, this one strong family realizes every residue modulo every odd modulus but omits infinitely many integers in each odd arithmetic progression having infinitely many positive terms. It gives a literal tournament example where complete modular automata and a sparse exact spectrum coexist, not a claim that the full tournament spectrum is sparse.
