# All fixed-depth symmetry corrections to A000568

The [double-triangle asymptotic](MODULES.md) is one layer of a general hierarchy. Define the defect of a permutation $g$ on $n$ vertices to be $d(g)=n-\#\{\text{cycles of }g\}$, counting fixed vertices as cycles. Every automorphism of a tournament has odd order, hence its defect is even. For fixed positive even $D$, let $B_D(n)$ count ordinary tournament classes with a nontrivial automorphism, but **none** of defect less than $D$.

For each odd prime $p$ with $p-1\mid D$, put
$$
r=\frac{D}{p-1},\qquad s=pr,\qquad
\tau_{p,r}=\Pr(K=C_p),
$$
where $K$ is the block-preserving automorphism group of a uniform fixed table on the $s$ moved vertices, and $C_p$ is the designated simultaneous rotation. Then for fixed $D$ there is $c_D>0$ with
$$
B_D(n)=
\left(
\sum_{\substack{p\text{ odd prime}\\p-1\mid D}}
\frac{p\,\tau_{p,r}}{(p-1)p^r r!}\,
\frac{2^{\binom n2-Dn+D(s+1)/2}}{(n-s)!}
\right)\bigl(1+O_D(e^{-c_D n})\bigr).
$$
The same formula counts just the strong classes in this layer. The finite sum retains terms of the same exponential scale with different polynomial factors. Its largest-support term, always $p=3$, dominates polynomially.

## Why exactly these coupled prime actions appear

An automorphism $g$ of minimum nonzero defect has prime order. For composite order, a proper nontrivial subgroup of prime order refines the vertex-orbit partition strictly: otherwise its generator would generate the same orbits on every cycle of $g$, which requires its exponent as a power of $g$ to be coprime to every cycle length and hence to the full group order. Thus an action at defect $D$ has $r$ cycles of length $p$, with $s=pr$ moved vertices.

The automorphisms of the finite moved table that preserve each block separately form $K$. They all extend to the full tournament, since each fixed vertex sees an entire block uniformly. Each such element has defect at most $r(p-1)=D$. If $K>C_p$, the action of $K$ on any one of the $p$-vertex blocks cannot be free, since $|K|>p$. Some nonidentity element therefore fixes a vertex of that block. That block contributes at most $p-2$ to its defect, so the element has defect strictly below $D$. Thus **exactly** the finite carriers with $K=C_p$ can survive in this stratum.

The mixed-mask interaction graph of [TARGETS](TARGETS.md) has independent edge probability $1-2^{1-p}$ under a uniformly chosen fixed table. Its connectivity is necessary: disconnected components allow independent rotations of smaller defect. For $p=3,5$ it is also sufficient. In a regular $p$-vertex block any point stabilizer preserves its in- and out-neighborhoods of size at most two; an odd-order permutation fixing both sets acts trivially. So each internal group really is just the translations, and connectivity rules out other blockwise shifts. For higher $p$ one must actually test $K$; for example a Paley-7 block has connected one-vertex interaction graph but carries an order-three automorphism of defect four.

## Entropy bound and normalization

For a designated defect-$D$ action on $s$ vertices the number of pair orbits is
$$
\binom n2-Dn+D(s+1)/2.
$$
Indeed each fixed vertex loses $r(p-1)=D$ bits of incidence to the moving set; the contributions within a moving block and between moving blocks are $(p-1)/2$ and $p$, respectively. The finite fixed table defining $K$ is independent of choices on all other pair orbits.

An extra automorphism $h$ outside $K$ cannot preserve every $g$-orbit separately. If its support has size $t$, then $\langle h,g\rangle$ merges at least $\max(1,\lceil2t/3\rceil-D)$ of the $g$-vertex-orbits (counted as the decrease in their number). Each vertex fixed by both permutations loses at least that many independent incidence bits. At $t\le n/2$, summing over at most $n^t$ maps yields exponential decay: split bounded $t$ from $t$ large enough in terms of $D$, when the merge bound is at least $t/3$. At $t>n/2$, the $h$-fixed condition alone costs $\Omega(n^2)$ bits among all orientations, hence still costs $\Omega(n^2)$ among $g$-fixed orientations. This is uniform in the finite moved carrier. Therefore almost all tables have group exactly $K$.

There are $n!/[p^r r!(n-s)!]$ choices of this designated action. Every admissible typical class has exactly $p-1$ chosen nonidentity actions and $n!/p$ labelings. Multiplying the table count by the action count and by $p/((p-1)n!)$ gives the displayed term. Classes with no action of defect $D$, or atypical tables contributing a wrong generic group, are exponentially smaller: at bounded support, defect $>D$ costs at least $(D+2)(n-t)$ pair bits from fixed-vertex incidences alone; at larger support $t$ the general bound $t(2n-t-1)/3$ dominates $Dn+cn$ after splitting off a sufficiently large constant support depending on $D$. Normalizing nonidentity Burnside mass controls the residual class count since each nonrigid class contributes at least $2/3$. The asymptotic has a positive leading term from the connected $p=3$ tables.

The strong result follows from the independent fair bits between the $n-s$ singleton orbits and on their incidences to each moving block. Directed cuts at the orbit level have exponentially small total probability (for any cut the number of independent bits across it is at least $k(m-k)-\binom r2$, $m=n-D$). Each block is strong, so orbit-level strong connectivity lifts.

## Consequences

* After excluding single rotating prime modules and the coupled double-triangle actions, the leading residual term is the coupled **triple triangle**, $D=6,p=3$. Other defect-six actions are single $7$-cycles and are excluded by the module filter. Here $\tau_{3,3}=27/32$, giving $2^{\binom n2-6n+30}/(128(n-9)!)$.
* The full hierarchy differentiates isolated rotation modules from genuinely coupled finite cores. For each fixed layer one can compute the exact finite carrier weight without enumerating the growing tournament.
* A prime factor supplied by symmetry still has to divide $H$, so none of this establishes realizability at every allowed odd target. In particular the dominant $p=3$ correction is invisible at prime-to-three targets. The finite carrier classification explains why very large improvements to scalar A000568 do not, by themselves, describe its exact-$H$ fibers.
