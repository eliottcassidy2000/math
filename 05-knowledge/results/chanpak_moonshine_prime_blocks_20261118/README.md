# Prime selection: orbit quotients and physically realizable lattice points

**Proved:** a sevenfold-symmetry quotient that is in the signed profile lattice and the nonnegative coefficient cone, but is not a tournament profile; the first entry of genuine sevenfold symmetry on the prime-partner axis; a signed window with an exact physical decoder; and a polynomial lower bound on the offset scale needed by any finite signed-reader bank for almost all primes. **Open:** the arbitrary-prime selector and Chan–Pak.

Let $h(T)=H(T)$ and $B(T)=f_T(2)$, where $f_T(q)$ counts ordered covers by $q$ directed paths, allowing empty slots. Set $f_\varnothing(q)=1$. The [Rédei–Berge identity of Grinberg–Stanley](https://arxiv.org/abs/2307.05569) expresses $f_T(q)$ as $G_T(2q,q)$, marking disjoint directed odd cycles by $2q$ and uncovered singletons by $q$. All degrees have the parity of $n=|T|$, the leading coefficient is one, and lower coefficients are nonnegative and even. Consequently
$$
B(T)\equiv(-1)^n h(T)\pmod3,\qquad 2h(T)\le B(T)\le 2^n h(T),\qquad B(T)\ge2^n
$$
for nonempty $T$. The last bound counts assignments of vertices to slots, using Rédei on both induced sets. The upper bound uses coefficient positivity. Ordinal joins multiply both counts: concatenate the paths in each slot in component order.

## 1. What the Monster–Baby Monster analogy actually transfers

A Monster involution of class $2A$ has centralizer $2.B$, the double cover of the Baby Monster. The [old Sylow fixed-point calculation](../session32_square_monster_20261023/SYLOW_FIXED_POINT_PARTITION.md) distinguishes unchanged, enlarged, and missing prime support by stabilizer valuations, not by matching the lengths of prime lists. The [level-two follow-up](../session32_square_monster_20261023/monster/MONSTER_BABY_PRIME_PARTITION.md), using [Aricheta's theorem](https://doi.org/10.3842/SIGMA.2019.007), finds a Frobenius-conjugate pair of coarse 2-isogeny edges over each of $29,41,59,71$ whose underlying endpoints are rational. Rational endpoints do not imply rational marked data. Similarly, [Griess and Lam's rootless $EE_8$ pairs](https://arxiv.org/abs/0806.2375) use distinguished sublattices and their involutions, not just numerical lattice invariants. The transferable idea is to **retain the mark when passing to a quotient**. There is no claimed moonshine identity for tournament counts.

A tournament automorphism fixing an ordered Hamiltonian path or ordered spanning cover fixes each vertex. Thus $|\operatorname{Aut}T|$ divides the count of either object. By [Busch's lower bound](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf), a strong tournament satisfies $n\le1+3\log_5 h$. For $h=p$ prime and $p\ge5$ this is smaller than $p$; a nonidentity automorphism would have order $p$ on fewer than $p$ vertices. Hence the prime core is rigid. For $h=7p$, with $p\ge11$ prime, again $n<p$. So the group is trivial or $C_7$. Both kinds of targets have exactly one nontrivial strong component: SCC counts multiply and 7 is excluded. Neither sporadic simple group is realized as an automorphism group here.

If $n<14$ and $C_7$ acts nontrivially, it moves a single orbit $O$ of size seven. The rest of the vertices are fixed and see $O$ uniformly, so it is a module. Its cyclic step set takes one of each pair $\pm1,\pm2,\pm3$. There are two types, the six unit multiples of $\{1,2,3\}$ (carousel) and the two unit multiples of $\{1,2,4\}$ (Paley). The exact $(h,B)$ readings are $(175,1988)$ and $(189,2016)$. Paley's translations and dilation by 2 give an internal group of order 21, which extends by fixing everything outside. Thus Paley is excluded at $7p$ for $p\ge11$. This is why the old [deletion-compatible Paley sheets](../chanpak_fano_square_compression_20261118/SEVEN_SHEETS.md) cannot produce prime quotients.

Through order 9, the only strong quotients after contracting $O$ are the singleton and the directed triangle, giving 175 and $3\cdot175+(1988-350)=2163$. Neither has prime quotient by seven. At order 10, substituting the carousel at the four possible ports of the unique strong four-vertex quotient gives exactly $4151=7\cdot593$ and $5789=7\cdot827$, by literal enumeration of the four ports. These are the only symmetric strong prime-partner cores at this order. The 4151 carousel instance already occurred in the [prime rail](../chanpak_prime_frontier_selector_20261118/PRIME_RAIL.md); the result here is the classification of the symmetry sector, not a priority claim for its counts.

In any $C_7$ action, even with further moving orbits, choose a seven-vertex orbit $O$, choose $a\in O$, and put $F=V\setminus O$. Each orbit of Hamiltonian paths has exactly one representative whose first vertex of $O$ is $a$. Splitting at that marked vertex gives
$$
 h(T)/7=\sum_{S\subseteq F} \operatorname{end}(T[S\cup\{a\}],a)\,
                \operatorname{start}(T[V\setminus S],a).
$$
The two paths meet only at $a$. This is an exact positive count with first-contact marking. No ordinary tournament with exactly this marked family of paths has been constructed in general.

A physical supply on the upper axis: replace two vertices of a directed triangle by carousel and a nonempty $A$, leaving the third a singleton. The [cyclic two-slot formula](../seven_gap_atlas_20261117/README.md) gives $h_{\rm out}/7=259B(A)-209h(A)$. The replay also tests this via the first-contact formula. In particular $A=\mathrm{TT}_{2t}$, $t\ge1$, gives quotients $259\,2^{2t}-209$, prime at $t=1,5$ (827 and 265007). For $3\nmid h(A)$, odd order of $A$ forces divisibility by 3 of the quotient, while even order avoids it. This family realizes the upper axis with actual symmetry; no all-order lower-axis realization of the quotient is inferred.

## 2. An integral lattice point with no physical lift

**The divided carousel profile $(h,B)=(25,284)$ is not the profile of any tournament.** Its scalar quotient 25, by contrast, is a count: join two count-five tournaments.

Indeed $B\ge2^n$ forces $n\le8$, and the mod-3 rule forces odd $n$. Orders 1 and 3 have too few paths. At order 5, writing $f=q^5+u q^3+v q$ gives $u+v=24$ and $80+6u=284$, so $u>24$. At order 7, writing $f=q^7+u q^5+v q^3+wq$ gives $u+v+w=24$ and $5u+v=18$. Nonnegative even coefficients allow only $u=0,2$. But $u$ is twice the number of directed triangles. With no triangle the tournament is transitive. With a single triangle, every outside vertex sees it uniformly, since mixed contacts would create another triangle. Every vertex dominating the triangle dominates every vertex it beats, by the same triangle test, and both outside parts are transitive. This tournament has exactly three paths, not 25.

The signed lattice nevertheless accepts the profile. At fixed order $n\ge5$, the integer span of polynomial differences has basis $2q^d$ for $d=n-2,n-4,\dots,d_0$, with $d_0=1$ at odd order and 2 at even order. For the reverse inclusion, reverse the extreme arc of a transitive block of odd length $\ell$, pad with $n-\ell$ ordinal singletons, and subtract the transitive polynomial. Its unique spanning cycle contributes 2 at the lowest degree $n-\ell+1$. All shorter odd cycles use the same reversed arc and only contribute higher degrees. Induction on odd $\ell$ isolates the monomials. This is the two-reading specialization of the earlier transverse packing-lattice argument.

At readings 1 and 2 the resulting lattice has basis $(2,2^{d_0+1})$ and $(0,3\cdot2^{d_0+1})$. Generate the second by subtracting the two lowest-degree vectors; higher ones differ by multiples of it because $4^j-1$ is divisible by 3. At order 7,
$$
(25,284)-(1,128)=(24,156)=12(2,4)+9(0,12).
$$
Even the nonnegative formal polynomial $q^7+2q^5+8q^3+14q$ has the readings. What fails is compatibility between the triangle and all the other cycles. This gives an all-order counterexample to “divide legitimate readings by the symmetry factor, then choose a block from the lattice.”

## 3. Real lifts: use the residue to choose a block, then check height

Use the [signed wrapper](../chanpak_six_signed_height_pairs_20261118/README.md). For $b\to A\to D$ and $I\to b\to D\setminus I$, let $G$ count ordered two-part covers of $D$ with nonempty first path ending in $I$ and second optionally empty. Let $Q$ count those with second path empty or starting outside $I$. The literal read is $G B(A)+(h(D)+Q-G)h(A)$: there are $h(D)+G+Q$ templates with one $A$-run and $G$ with two. Every $A$-run except a possible initial one must start immediately after $b$.

Four checked marked $D$'s give $(G,\beta)=(6,1),(42,7),(96,-1),(672,-7)$ with $\beta=h(D)+Q-G$. The positive pair uses four additional ordinal singletons on the input, which multiplies $B$ by 16 and leaves $h$ fixed. Thus the **same arbitrary $A$** gives $96B\pm h$ and seven times each, by changing readers rather than uniformly inflating the lower output. The positive readers were already recorded; the negative pair comes from the signed-height packet's discovery table and is checked on literal arcs here.

Eight nine-vertex sources share $B=32852$ with $h=1945+6j$, $0\le j<8$. Their masks are independently replayed; we do not need completeness of the discovery census. Pad by $t$ ordinal singletons. With $C=96\cdot32852=3153792$ we realize
$$
p=C\,2^t\pm(1945+6j),\qquad t\ge0,\quad 0\le j<8,
$$
and every corresponding $7p$. The 16 signed offsets run over every unit modulo 48 once. A prescribed $p$ therefore chooses the source and sign by its residue. The remaining test is whether the residual divided by $C$ is a positive exact power of two. For example, the prime $51,671,730,073$ chooses $j=0$, plus, and $t=14$; the literal sevenfold partner has $361,702,110,511$ paths. The window is an exact constructor, not a congruence witness. Its fixed-center support is still only logarithmic.

## 4. Height tax for any finite bank of fixed readers

**Theorem.** Fix finitely many readers $p=\alpha B(A)\pm h(A)$ with fixed $\alpha\ge1$. For every
$$
0<\theta<\frac1{2+3\log_5 2}\approx0.3038,
$$
the number of distinct targets $p\le X$ produced with source $h(A)\le X^\theta$ is $O(X^{1-\delta})$ for some $\delta>0$. In particular, if the readers cover all primes, almost every prime must use a source with $h(A)>p^\theta$. This rules out fixed or polylogarithmic offsets *even with completely free choice of source topology* in this architecture.

Proof: split the SCCs of $A$ into singletons, finitely many nontrivial small types, and large components. Given $\eta>0$, choose the threshold so that Busch gives $m\le(3+\eta)\log_5 h$ for every large component of order $m$. Since $h$ and $B$ multiply, the ordinal order does not affect them. If $h(A)\le Y$, the finitely many small types have multiplicities at most $\log_3 Y$, giving $(1+\log Y)^{O_\eta(1)}$ choices. The join $L$ of large components has total order at most $(3+\eta)\log_5 h(L)$, hence $B(L)\le h(L)^{1+(3+\eta)\log_5 2}$. There are at most $O(Y^{2+(3+\eta)\log_5 2})$ distinct pairs of large-component readings, including the empty join. Singleton multiplicity contributes at most $O(\log X)$ choices, since $p\ge(\alpha-1/2)B(A)$ and $B(A)\ge2^{|A|}$. Substitute $Y=X^\theta$ and choose $\eta$ small; the exponent remains below one after absorbing the polylogarithms. The prime number theorem gives the prime corollary. Also for integer $\alpha>1$, prime outputs require $\gcd(h,B)=1$, since a common divisor divides $p$ and is at most $h<p$. In particular the source components must be rigid. This extra condition is not needed for the power saving.

**Next proof target:** choose a growing rigid source from the target, or change the reader itself. On the group route one must preserve the first-contact mark or alter the reader to compensate. On the lattice route the point on the target height line has to be realizable, not merely in the signed span. These two restrictions explain why extracting the right prime in a quotient and finding the right residue in a window do not yet constitute a Chan–Pak proof.

Replay with **python3 replay.py**. To build the example above, run **python3 replay.py 361702110511 --scale 7**. The proof and replay do not need an external census.
