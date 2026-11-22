# Two independently branching rooted flanks with a six-count switch

**PROVED.** Two independently grown and decoded rooted tournaments can be inserted into opposite ports. A seven-vertex middle switch changes the absolute count by exactly six, independent of both branch histories. The two images occupy disjoint residues modulo 32, and have a complete inverse by positive divisor tests followed by the two individual branch inverses. They are sparse, not a solution of Chan–Pak.

## Four-channel graft

For a tournament $T$, let $H(T)$ count its Hamilton paths, with $H(\varnothing)=1$. Write $R_x(T)$ for the count of choices of a nonempty path ending at $x$ and a complementary Hamilton path, allowing the latter to be empty. Write $L_x(T)$ for its start version, and $R_x^-,L_x^-$ for the corresponding signed versions, with weight $(-1)^{|P|}$ on the selected path $P$.

The signed end sum is minus the number of full paths **starting** at $x$. Indeed, deleting the final $x$ from the selected path gives the recursion
$$
R_x^-(T)=-H(T-x)-\sum_{u\to x}R_u^-(T-x).
$$
Inductively the right side is minus the number of paths in $T-x$ starting at a vertex beaten by $x$, which is the asserted count. The singleton case gives $-1$ on both sides. Converse gives the start version.

Take the ordinal order $I\to B\to O$, with also $I\to O$, and reverse only the two arcs to $x\to i$ and $o\to y$, where $i\in I,o\in O,x\ne y\in B$. Assume $L_i^-(I)=R_o^-(O)=0$. Put
$$
h=H(I),\quad k=H(O),\quad U_I=L_i(I),\quad U_O=R_o(O).
$$
Let $P_2(T)=\sum_{S\subseteq V(T)}H(T[S])H(T-S)$, and $D_{y,x}(T)$ count full paths starting at $y$ and ending at $x$ (zero unless both marks are present). Define the middle signature
$$
\begin{aligned}
g&=H(B),& r&=R_x(B),& l&=L_y(B),\\
M_1&=\sum_{S\subseteq B}D_{y,x}(B[S]) P_2(B-S),\\
M_2&=\sum_{S,T\subseteq B,\ S\cap T=\varnothing}
 H_{\to x}(B[S])\,H_{y\to}(B[T])\,H(B-(S\cup T)),& m&=M_1+M_2.
\end{aligned}
$$
Here $H_{\to x}$ and $H_{y\to}$ are the full marked-end and marked-start counts, zero on sets without the mark. The graft count is
$$
H=g hk+r U_I k+l h U_O+m U_I U_O. \tag{1}
$$

**Proof.** Use the odd-cycle packing expansion of $H$ (the $q=1$ case of THM-4940): each disjoint nontrivial directed odd cycle contributes a factor of two; unused vertices have unit weight. If no mixed cycle uses a reversed arc, the contribution is $g hk$. If only $x\to i$ participates, the mixed cycle concatenates a path in $I$ starting at $i$ and a path in $B$ ending at $x$; the untouched $O$ contributes $k$. The oddness filter $2[|P|+|Q|\text{ odd}]=1-(-1)^{|P|+|Q|}$ leaves the unsigned product $r U_I k$ because $L_i^-=0$. The would-be closed walk on the two endpoints alone is automatically excluded by the same parity filter; every other formal pair uses a valid forward joining edge. The case of only $o\to y$ is analogous.

If the reversed arcs lie in the same odd cycle, start at $y$, follow a path in $B$ ending at $x$, follow the reverse arc to $i$, go through one path in $I$, then through an optional disjoint path in $B$, then through one path in $O$ ending at $o$, and close via $o\to y$. There can be no other runs in either port because each reverse edge is unique. Weighting the optional middle path and its complement is $P_2$ on the unused middle vertices. The oddness filter again leaves only the unsigned product, giving $M_1 U_I U_O$. If the two reversed arcs lie on distinct odd cycles, the selected paths in $B$ ending at $x$ and starting at $y$ must be disjoint, and each oddness filter separately eliminates its signed term. This gives $M_2 U_I U_O$, proving (1).

The graft is strong in our applications: $i$ reaches all of $I$, and all of $O$ reaches $o$, via the unit-end paths of the rooted flank constructions. The flanks have at least two vertices. Then $i$ can reach all of $B$, including $x$ via an additional vertex of $I$; all of $B$ enters $O$, using another vertex if needed to avoid the missing $y\to o$. If $y$ reaches $x$ in $B$, the two reversed arcs close these reachability chains, giving strength.

## Two literal sources, both absolute counts

Start each flank from the transitive pair rooted at its last vertex, and apply any word in the repository's [three rooted migrations](../chanpak_pubsub_branching_20261117/PROOF.md). Independently on each side the input Hamilton count starts at one, the marked weighted end is $3h-1$, precisely one full path ends at the root, and none starts there. The letters act in chronological order by
$$
A(h)=82h-21,\quad B(h)=5976h-1579,\quad C(h)=96264h-24407.
$$
Take the converse of the left flank and leave the right one unchanged. Then the two signed sums required in (1) are zero, and $U_I=3h-1,\ U_O=3k-1$.

Use the following outgoing masks, with $x=0,y=1$; the bit at column $j$ records an arc to vertex $j$. In the first source $y\to x$ directly; in the second $y\to2\to x$.
 
| Source | Outgoing masks | $g$ | $r$ | $l$ | $M_1$ | $M_2$ | $m$ |
|---|---|---:|---:|---:|---:|---:|---:|
| lower | $64,93,17,101,41,7,52$ | 87 | 283 | 283 | 310 | 321 | 631 |
| upper | $18,116,97,39,12,81,25$ | 141 | 265 | 265 | 254 | 383 | 637 |

For replay, compute full end counts by the subset recurrence $E(S,v)=\sum_{u\in S-v:u\to v}E(S-v,u)$, starting at 1 on singleton subsets. A $y$-starting copy instead initializes only $\{y\}$ to 1. Apply full end counts to the converse for start counts, and use the displayed finite sums. The companion script independently checks the literal graft counts for four transitive flank pairs, and a nonempty $A$-migrated flank with the other flank the seed.

The signature difference is $6(9,-3,-3,1)$ in $(g,r,l,m)$, i.e. the tensor of the preserved one-port row $(3,-1)$ with itself, multiplied by six. In fact a constant switch of size $d$ on the whole independent grammar is possible if and only if its signature change is $d(9,-3,-3,1)$: after substitution the $hk,h,k$ coefficients must each vanish on an infinite Cartesian product, and solving those three linear equations gives exactly this direction. Before imposing either invariant, the switch changes the output by
$$
6(3h-U_I)(3k-U_O).
$$
Therefore the actual counts, for all independent legal histories including arbitrary repetitions, are
$$
N_s=7464hk-2176(h+k)+631+6s,\qquad s=0,1. \tag{2}
$$
For example $h=k=1$ gives $3743,3749$; $h=A(1)=61,k=B(1)=4397$ gives $1992271711,1992271717$. These are count witnesses, not an assertion about primes.

## Complete inverse and a persistent local hole

Every reachable one-port count is $1\bmod4$, since all three branch maps preserve that residue. Consequently (2) gives $N_0\equiv31\pmod{32}$, $N_1\equiv5\pmod{32}$. In particular the images are disjoint already modulo 8 and together omit every odd residue except these two out of the sixteen modulo 32. This holds across all depths on both flanks, not merely for a fixed pair. Each port affects the output even at the smallest point: the coefficient in $h$ with $k$ fixed is $7464k-2176>0$, and vice versa.

In general (1) on the two invariant lines gives $N=a hk-bh-ck+d$, where $a=g+3r+3l+9m,\ b=l+3m,\ c=r+3m,\ d=m$. Its exact completion is
$$
(ah-c)(ak-b)=a(N-d)+bc. \tag{3}
$$
Both sides of the product are positive at physical $h,k\ge1$. Here $a=7464,\ b=c=2176,\ d_s=631+6s$. To decode a target, first identify $s$ from the residue, if possible. Factor the positive right side of (3), take each ordered positive divisor pair congruent to $(-c,-b)\bmod a$, and recover $h,k$. Test each independently by the [exact rooted-history inverse](../chanpak_pubsub_branching_20261117/PROOF.md). This lists every ordered pair of physical histories and no others; it is a finite inverse, not a polynomial-time factoring claim. If one prefers a literal binary alphabet, restrict both ports independently to $B,C$. Then the final letter is $B$ at $h\equiv2\bmod3$ and $C$ at $h\equiv1\bmod3$ after testing the seed $h=1$; use exact affine division and reject invalid or nonpositive predecessors. The unrestricted inverse uses the known valuation tag $v_2(27h-7)$ for the trailing $A$-run.

As another minimal middle, $B$ can consist of the two marks with $y\to x$. Its $(g,r,l,m)=(1,2,2,2)$ gives $N=31hk-8(h+k)+2$ and
$$
(31h-8)(31k-8)=31N+2.
$$
Thus any positive odd $N$ with $31N+2$ prime is excluded from *this* simple architecture independently of the branch-language test. Its omitted progression of transformed primes is infinite by Dirichlet; this does not say $N$ itself is prime.

## Capacity and the two outside prompts

The one-port grammar has $\Theta(X^\sigma)$ addresses through $X$, where $\sigma>0$ satisfies $82^{-\sigma}+5976^{-\sigma}+96264^{-\sigma}=1$ (about $0.14509$). Formula (2) is bounded above and below by positive constants times $hk$. Dyadic grouping gives $O(X^\sigma\log X)$ ordered pairs through $X$; to get the matching lower bound, use shells with a sufficiently large fixed ratio so the difference of the one-port upper and lower bounds leaves $\gg Y^\sigma$ addresses per shell, and let the other port range through $cX/Y$. This counts representations, not distinct outputs. Fixing one flank at the seed gives $\Omega(X^\sigma)$ distinct values. A second independent branch changes the number of presentations by a logarithm, not a new power of $X$; this explains why even an exact two-port repair and inverse are not close to complete coverage on their own.

[Pasten–Vidaux, §§4–6](https://arxiv.org/pdf/2611.24748) separates local approximation from the propagation of an actual valuation witness. Here the rooted suffix valuation is one exact witness, and (3) makes $v_\ell(a(N-d)+bc)=v_\ell(ah-c)+v_\ell(ak-b)$ a second, necessary local test for all primes $\ell$ before inversion. The positive integer divisor check is indispensable: local admissibility alone cannot turn a candidate into a physical history. No arithmetic conjecture is used.

The [RH/GRH announcement remark](https://alpo.ge/remark.pdf) suggests following information that survives repeated lifts and completing a determinant before taking a scalar read. Here repeated migrations preserve a signed functional at each physical port; their tensor supplies the exact six-count switch, while the bilinear determinant completion gives an inverse rather than a positivity claim. The count of presentations explicitly distinguishes multilevel address multiplicity from distinct target coverage. No result asserted in the announcement is a dependency.

**Next target:** a wider middle signature bank changing $a,b,c$ and the low residue without destroying independently decodable ports; merely growing the same two flanks preserves the obstruction above. A switch of gap two in the specific four-channel model remains open here.

Reproduce the middle signatures and new literal graft checks with `python3 check_two_flanks.py`.
