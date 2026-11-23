# Prime-power control of the two physical fork readings by ordinal reagents

For a nonempty tournament $M$, let $D=H(M)$ and $Q=Q_M(2)$, the number of ordered two-path covers with empty slots allowed. Both multiply under ordinal join $X\to Y$: a Hamilton path cuts uniquely at the ordinal boundary, and each path slot does so independently. Thus each ordinal reagent $(a,b)$ acts as the diagonal map $(D,Q)\mapsto(aD,bQ)$. In the seven-fold fork the parent read is $7(40Q+D)$, the center is $21(40Q+D)$, and its two reversals move by $\pm2D$.

## Full unit-pair reachability

Take the finite bank of all distinct $(D,Q)$ among tournaments on one through six vertices (39 profiles). The attached checker constructs it without a database: at each order fix a spanning Hamilton path by relabeling, then enumerate each orientation of the remaining arcs; a standard subset recursion gives $H(M[S])$, and $Q=\sum_S H(M[S])H(M[\bar S])$. Every bank entry retains a literal witness tournament.

**Theorem (finite-certified prime range, all powers).** For every odd prime $p<1000$ and every $e\ge1$, all pairs in $((\mathbb Z/p^e)^\times)^2$ occur as readings of nonempty finite ordinal joins of bank tournaments. The certificate is finite at each $p$, independent of $e$.

Keep only $(a_i,b_i)$ with $p\nmid a_i b_i$. Let $m=p-1$, fix a primitive root modulo $p$, and write $v_i=(\log a_i,\log b_i)$ mod $m$. The mod-$p$ span is full exactly when the gcd of the two-by-two minors of the integer matrix with columns $v_i,(m,0),(0,m)$ is one. This is the index of the generated sublattice (or use Bézout on the minors and each two-column adjugate). The second test is that the Fermat quotient pairs
$$
F_i=\left(\frac{a_i^{p-1}-1}{p},\frac{b_i^{p-1}-1}{p}\right)\pmod p
$$
span $\mathbb F_p^2$. The checker computes the first gcd and a nonzero determinant of the second pairs at every one of the 167 indicated primes.

These tests suffice for all powers. A pair in the mod-$p$ image has a word from the first test. At level $j$, powers $(p-1)p^{j-1}$ of the bank pairs are trivial modulo $p^j$ and their next digits modulo $p^{j+1}$ are precisely $F_i$, by binomial lifting. Two suitable generators therefore correct each next pair of digits independently. Negative exponents used in the mod-$p$ solve cause no difficulty: at modulus $p^e$ each unit generator has order dividing $(p-1)p^{e-1}$; reduce each integer exponent modulo that number to get a nonnegative-word realization. If all exponents vanish one may concatenate one bank tournament that many times to realize the identity nonemptily.

### Deterministic compressed witness construction

For a requested pair modulo $p^e$, use the target log vector modulo $p-1$. For any two columns $v_i,v_j$, their determinant times the target is $\det(\mathrm{target},v_j)v_i+\det(v_i,\mathrm{target})v_j$. Extended gcd on the minors expresses the target as an integer combination; the two artificial lattice columns vanish modulo $p-1$. Then correct each higher digit with the inverse of one fixed two-by-two Fermat quotient matrix modulo $p$, and finally reduce all exponents modulo the unit exponent. The output is a compact list of adjacency rows and nonnegative repetitions; it need not build or search a group graph of size roughly $p^{2e}$. Run:
$$
\texttt{python chanpak-ordinal-reagents-check.py --target P E D Q}.
$$
The literal ordinal join may of course be much larger than its logarithmic-size description.

For example, for $p\ne5$ in the certified range, choose $Q\equiv1$ and $D\equiv-40+p^e\pmod{p^{e+1}}$ to make $v_p(40Q+D)=e$ while both coordinates remain units. This controls an exact valuation in the physical fork read; it neither makes the integer cofactor equal to $p$ nor requires persistence under every uniform context. At $p=5$, unit inputs cannot make the read divisible by five.

## CRT caveat and the unavoidable finite-bank obstruction

Full projections to individual prime powers do not prove full reachability at their product: the same word must satisfy all factors and orders can share divisors. Independently, the included finite group closure proves full reachability for each of the 127 odd moduli $3\le q\le255$, including composites; it does not assume a CRT lift from separate tests. For any other fixed composite modulus, discard nonunit profiles and use the group they generate, potentially a proper subgroup. Independent CRT combinations require an additional joint test.

In fact **no finite bank** of positive integer pairs can have full unit-pair reachability at every odd prime. Let $L$ be eight times the product of all odd prime radices occurring in its entries. For primes $\ell\equiv1\bmod L$, each such radix is a square modulo $\ell$, by quadratic reciprocity, and 2 is a square because $\ell\equiv1\bmod8$. Every generated pair therefore lands in the product of the square subgroups, which is proper. Dirichlet gives infinitely many such primes. The point of the explicit certificate is broad controlled reachability in a fixed range, not universal modular generation from six-vertex reagents.

## Real geometry and a nonreturn inequality

The determinant-normalized real ordinal map is diagonal with reciprocal positive eigenvalues $\sqrt{a/b},\sqrt{b/a}$. A nonempty reagent $V$ has $b/a=Q_V/H(V)\ge |V|+1>1$: cutting a full path at each interior position gives distinct two-nonempty-slot covers, and either slot can hold the full path. Hence every such normalized physical graft is hyperbolic (a Lorentz boost), never a nonidentity parabolic or elliptic motion. The only nondegenerate homogeneous quadratic form preserved conformally even by the singleton reagent $\mathrm{diag}(1,2)$ is, up to scaling, $DQ$: the three monomials $D^2,DQ,Q^2$ have distinct eigenweights $1,2,4$, and the other two eigenforms are degenerate. In the fork coordinates $(C,\Delta)=(840Q+21D,2D)$ this split form is proportional to $\Delta(C-21\Delta/2)$. This is a classification of actual ordinal transitions, not of other projective chart changes.

Operationally $R=C/\Delta-21/2=420Q/D$ is multiplied exactly by $b/a$. After $k$ nonempty reagents it grows by at least $2^k$, while $\log R$ records their additive positive displacements. Thus no nonempty physical ordinal sequence returns the fork to its original projective ratio, even though modulo a certified odd modulus the same diagonal group may contain every invertible ratio and a positive word realizing its identity.

### One finite CRT shell suffices for all exponent tuples on a fixed support

Let $S$ be a set of odd primes. Keep only bank profiles whose two entries are units at every prime of $S$. For each $p\in S$, put $K_p=\max_{r\in S,r\ne p}v_p(r-1)$ (empty maximum zero), and set
$$
B_0=\prod_{p\in S}p^{K_p+1}.
$$
It suffices to verify full unit-pair reachability modulo $B_0$, and rank two of the Fermat quotient vectors separately for each $p$, **using only those jointly eligible profiles**. These finite tests imply full joint reachability modulo $\prod_{p\in S}p^{e_p}$ for every exponent tuple.

Proof: extend any requested exponent below $K_p+1$ to that boundary and ultimately project back. Start by attaining the pair simultaneously modulo $B_0$. To correct the $p$ digit from level $j$ to $j+1$ with $j\ge K_p+1$, let $L_p=\operatorname{lcm}_{r\in S\setminus\{p\}}\varphi(r^{e_r})$. Its $p$-valuation is $K_p$. A power $L_p(p-1)p^{j-1-K_p}$ of any jointly eligible reagent is entirely trivial in every other component; in the $p$-component it acts on the next pair of digits by its Fermat quotient vector multiplied by the unit $L_p/p^{K_p}$. Solve the two-dimensional digit correction at $p$, and repeat at each component. All exponents are nonnegative, so all resulting words are actual ordinal joins.

The checker certifies this stronger joint all-power assertion for supports $\{3,5\},\{3,7\},\{5,7\},\{3,5,7\},\{3,11\},\{5,11\},\{7,11\},\{3,5,11\},\{3,7,11\}$. For example, support $\{3,5,7\}$ uses the finite base $315$, with $20736$ unit pairs; the same bank thereby controls every $3^a5^b7^c$ jointly, not just its three projections.

There is also an explicit CRT failure for this very bank. Take the squarefree product of all primes dividing its first coordinates:
$$
m=3\cdot5\cdot11\cdot13\cdot17\cdot19\cdot23\cdot29\cdot31\cdot37\cdot41\cdot43.
$$
Every prime factor separately has full unit-pair range by the certified prime theorem. But the only jointly unit bank generators have $D=1$, and products containing any nonunit $D$ cannot later regain unit status. Thus every jointly unit reachable pair for this $m$ has first coordinate 1. This is a concrete obstruction to replacing the joint test by independent local tests.
