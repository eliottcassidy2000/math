# An exact half-cut fold for the joint path and blue statistics

The [paired-site model](PAIRED_SITES.md) folds blue paths to a single
signed transversal. Full Hamilton counts include paths that are not
anti-palindromic. They too admit an exact fold, but **nontransversal
half-cuts must be kept**. This is an exact per-chart transfer; it does
not sum all black incidences by a two-state matrix.

Write $\rho$ for the anti-involution. For a subset $S$ and
$i\in S$, let $E_S(i)$ count Hamilton paths on $S$ ending at $i$.
Use the ordinary endpoint subset recurrence, stopping at half size.
Starting-path counts on $R$ are $S_R(j)=E_{\rho R}(\rho j)$
by reversing and reflecting.

The formulas below are for $m\ge1$. At orders 0 (empty identity) and 1,
use $H=F=1$ directly.

## Even order $2m$

Split each path after $m$ vertices, so its cut is $S\mid R$,
with $R$ the complement of $S$, and
$$
 C(S)=\sum_{i\in S,j\in R} E_S(i) [i\to j] E_{\rho R}(\rho j).
$$
The mirror sends $S$ to $\rho R$ with the same value.
Sum $2C(S)$ over nonfixed cut orbits. Fixed cuts are exactly the
$2^m$ transversals of the paired sites. At such a cut the matrix
$B_{ij}=[i\to\rho j]$, $i,j\in S$, is symmetric. Its total
contribution to $H$ is
$$
 \sum_i B_{ii} E_S(i)^2+
 2\sum_{i<j}B_{ij}E_S(i)E_S(j);
$$
its fixed-path contribution to $F$ is
$\sum_i B_{ii}E_S(i)$. The first diagonal squares count independent
choices of left and reflected right paths. Only the diagonal in the
actual path-choice product, not all choices with equal endpoint, is
fixed.

## Odd order $2m+1$

Split each path as $S\mid x\mid R$ around its *middle path vertex* $x$,
which need not be the anti-involution's fixed center $c$.
The cut contribution is
$$
 C(S,x,R)=\left(\sum_{i\in S}E_S(i)[i\to x]\right)
 \left(\sum_{j\in R}[x\to j]E_{\rho R}(\rho j)\right).
$$
Again mirror cuts in pairs under
$(S,x,R)\mapsto(\rho R,\rho x,\rho S)$. A fixed cut has $x=c$
and $S$ transversal; its contribution is $f_S^2$ with
$f_S=\sum_i E_S(i)[i\to c]$, and $F=\sum_S f_S$.
In particular any proposed small absolute target $h$ must pass
$\sum_S f_S^2\le h$ before the nonfixed cuts are evaluated.
This retains exact port multiplicities; knowing only that a cut is
live or knowing only $F$ does not supply the missing doubled terms.

The proof is a bijection: cut an actual Hamilton path at the stated
position(s), and glue any valid endpoint-matched subpaths. The companion
replay compares both $H$ and $F$ with literal endpoint DP and signed
half-word enumeration. This adds a full Hamilton coordinate to the
paired chart and preserves its blue coordinate; it is not a scalar
recurrence for the spectrum.
