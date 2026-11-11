# Ready fibres, rigidity, and actual constructor transfer

Status: **PROVED** structural consequences below, using the explicitly
inherited readiness classification and source families. **FINITE-EXACT /
INHERITED** the complete strong-ready fibre at 221 and its absence at 1997.
No tournament census is added.

The main distinction is between the ordinary fibre $H(T)=h$ and the subfibre
with the actual root needed for a $3h+2$ constructor. The latter has a precise
form and a different growth law. The inherited strong unique-start grammar is
the closest mechanism; 221 at order eight is the hostile scalar fibre; the
failed implication is that ordinary multiplicity forces readiness; the useful
sidecars are the sole in-neighbour and the unique starting Hamilton word.

Sources:

- [Complete readiness classification, READY section 1](../../two_entry_funnels_20261103/READY.md).
- [Unique-start grammar and Fibonacci count, FUNNELS F1–F2](../../polarized_path_calculus_20261103/FUNNELS.md).
- [Actual $S_k,R_k$ families and all-order ready counts, TARGET_RETURN sections 1 and 4](../../berggren_target_return_20261111/reset/TARGET_RETURN.md).
- [S258 arithmetic interpretation](../../session258_endpoint_seed_20261111/arithmetic/RESEARCH.md).
- [Ordinary/labelled multiplicity compiler in this packet](PROOF.md).

## 1. Definitions and inherited inverse

Write $F_T(r)$ for the number of ordered spanning covers by two paths in
which the first path is nonempty and ends at $r$; the second may be empty.
Call $(T,r)$ **1-ready** if $F_T(r)-H(T)=1$. Write $Q_T(2)$ for the same
two-bin path-cover count with both bins allowed to be empty.

The inherited classification says that a strong 1-ready tournament is
exactly the following newborn construction:

1. Take a strong parent $P$ and a vertex $a$ with exactly one Hamilton path
   starting at $a$.
2. Append $r$, with $a\to r$ and $r\to v$ for every $v\ne a$ in $P$.

In particular $r$ has indegree one, $P=T-r$ is strong, and the construction
has a unique rooted inverse. The inherited two successive births, each
beaten only by the current root, send a 1-ready source of count $h$ to an
actual tournament of count $3h+2$.

## 2. Every strong 1-ready class is rigid and has a unique ready root

**PROVED.** In the newborn construction above, every old vertex other than
$a$ has child indegree at least two. The indegree-one vertices are therefore
either just $r$, or the pair $\{a,r\}$ when $a$ had parent indegree one.
In the second case their oriented edge $a\to r$ fixes both vertices under
every automorphism. In either case every automorphism fixes $r$ and its sole
in-neighbour $a$. Its restriction to $P$ fixes the unique Hamilton word
starting at $a$, hence fixes all vertices. Thus $T$ is rigid.

Any other ready root would have indegree one, so the only possible other
root is $a$. But deleting $a$ makes $r$ a global source, contrary to the
strong-parent requirement in the inverse theorem. Therefore $r$ is the
unique 1-ready root.

Consequently forgetting the mark gives a bijection from rooted strong
1-ready isomorphism classes to ordinary strong classes admitting a ready
root. Every such ordinary class contributes exactly $n!$ labelled tournaments.
The inherited rooted unique-start count at parent order $n-1$ is
$F_{2n-6}$, with Fibonacci convention $F_0=0,F_1=1$. Hence

\[
\#\{\text{ordinary strong 1-ready classes of order }n\}=F_{2n-6},
\quad n\ge4.
\tag{1}
\]

This is an unmarked-class consequence of the marked inverse, with no
division by an unknown root multiplicity.

## 3. Complete weak-padding classification and fixed-count stabilization

**PROVED.** If $H(T)>1$, then $(T,r)$ is 1-ready if and only if

\[
T=R\ \mathbin{\triangleright}\ TT_b,\qquad b\ge0,
\tag{2}
\]

where $R$ is strong and 1-ready at $r$, and every arc points from $R$ toward
the transitive tail. In particular every such $T$ is rigid and has exactly
one ready root.

To prove necessity, write the ordered SCC decomposition around the root's
component as $T=A\triangleright C\triangleright B$, allowing empty $A,B$.
An ordered cover whose first path ends at $r$ assigns all vertices of $B$ to
its second path. Restrictions to the three ordinal blocks give the exact
factorization

\[
F_T(r)=Q_A(2)F_C(r)H(B),\qquad
H(T)=H(A)H(C)H(B).
\tag{3}
\]

Since their difference is one, $H(B)=1$, so $B$ is transitive. If $A$ were
nonempty, $Q_A(2)\ge2H(A)$: put the whole path in either bin. Also
$F_C(r)\ge H(C)$: split each Hamilton path just after $r$. It follows that

\[
F_T(r)-H(T)\ge H(A)H(C)=H(T)>1,
\]

a contradiction. Thus $A$ is empty, $C$ is strong of count $H(T)>1$, and
$F_C(r)-H(C)=1$. Conversely (3) shows that every transitive suffix preserves
readiness. Unique SCC order and Section 2 prove rigidity and uniqueness.

Let $M_m(h)$ count ordinary strong 1-ready classes of order $m$ and count
$h>1$, and let $J_n(h)$ count all ordinary 1-ready classes of order $n$ at
that count. Then the exact transfer is

\[
\boxed{J_n(h)=\sum_{m\le n}M_m(h)},\qquad
\boxed{L_n^{\rm ready}(h)=n!J_n(h)}.
\tag{4}
\]

The marked labelled count is also $n!J_n(h)$; with a specified root label
it is $(n-1)!J_n(h)$. By the inherited strong-order bound, $M_m(h)=0$ for
$m>1+3\log_5 h$. Thus $J_n(h)$ eventually becomes constant, while ordinary
$N_n(h)$ has positive polynomial degree at every attained $h>1$. Its ready
proportion tends to zero. If no strong ready class exists, no weak ready
class exists either, at any order.

Summing over $h>1$ in (4) gives the all-order unmarked count

\[
\sum_{h>1}J_n(h)=\sum_{m=4}^nF_{2m-6}=F_{2n-5}-1,
\quad n\ge4.
\tag{5}
\]

The excluded $H=1$ case consists of the transitive tournament; at order
$n\ge2$ its second vertex is its unique 1-ready root. This follows directly
from $F_{TT_n}(r)=2^{j-1}$ when $r$ is the $j$th vertex. It is separate from
the strong-parent classification.

## 4. A genuine injective constructor on ready fibres

For a ready tournament with $h>1$, append $x$ beaten only by $r$, then $y$
beaten only by $x$. The inherited count identity gives $3h+2$. The two
newborns are precisely the two indegree-one vertices in the final graph:
every old vertex other than $r$ acquires two in-neighbours; $r$ already had
indegree one and acquires $y$. The arc $x\to y$ canonically distinguishes
the two newborns. Deleting them recovers $T$, and the sole old in-neighbour
of $x$ recovers $r$.

This proves that the operation is injective even after forgetting all marks.
Its outputs are rigid, since an automorphism fixes $x,y,r$ and restricts to
an automorphism of the rigid source. Each birth also preserves strongness
when the input is strong: the old graph reaches the newborn through its sole
in-neighbour, and the newborn reaches the old strong graph through its other
vertices. Therefore

\[
\boxed{N_{n+2}(3h+2)\ge J_n(h)},\qquad
\boxed{L_{n+2}(3h+2)\ge(n+2)!J_n(h)}.
\tag{6}
\]

This is a literal transfer from the appropriate multiplicity to constructor
existence, with a graph-level inverse and exact vertex cost. Replacing
$J_n(h)$ by the ordinary $N_n(h)$ is unjustified. The output is not asserted
to be ready again.

## 5. Universal arithmetic rigidity criterion

**PROVED, using the inherited strong-order floor.** For an odd $h>1$, if

\[
\operatorname{lpf}(h)>1+3\log_5 h,
\tag{7}
\]

then every tournament with $H=h$, at every vertex order, is rigid. Hence
$L_n(h)=n!N_n(h)$ for all $n$.

The automorphism group acts freely on Hamilton words, so
$|\operatorname{Aut}(T)|$ divides $h$. Automorphisms preserve each SCC because
the condensation order is rigid. A nontrivial automorphism group would
therefore have a prime divisor $p\mid h$ in some strong component $C$.
As $\operatorname{Aut}(C)$ embeds in the symmetric group on $|C|$ vertices,
$p\le|C|$. But the inherited bound gives

\[
|C|\le1+3\log_5H(C)\le1+3\log_5h<\operatorname{lpf}(h)\le p,
\]

a contradiction. Singleton SCCs have no automorphisms. The criterion can
be sharpened by replacing the analytic bound with the largest order allowed
by the exact Busch minimum. It is sufficient, not necessary.

Every prime path total $p\ge5$ satisfies (7): the function
$x-1-3\log_5x$ is positive at 5 and strictly increasing for $x\ge5$.
The $p=3$ triangle is the sharp small exception. This general result does
not assert that any particular prime has a realizing tournament.

## 6. The two actual sources at $h_k=99\cdot2^k+23$

Use the inherited $k\ge1$ families $S_k$ of order $k+7$ and $R_k$ of order
$k+8$, both strong and both of count $h_k$.

### $R_k$: a rigid, uniquely marked usable source

The parent of $R_k$ is obtained from

    (34,244,249,243,225,192,129,1)

by replacing vertex 5 with $TT_k$. Its unique Hamilton path starting at
vertex 0 is the increasing path; the appended root is beaten only by 0.
Section 2 proves rigidity without any factorization of $h_k$. The parent
indegree of 0 is exactly 5, so the newborn is the unique indegree-one vertex.
The inherited marked coordinates are

\[
H=h_k,\quad H(R_k-r)=27\cdot2^k+3,\quad
s_{R_k}(r)=27\cdot2^k+2,\quad F_{R_k}(r)=h_k+1.
\]

The actual map in Section 4 produces count $297\cdot2^k+71$ on $k+10$
vertices. This source has the extra endpoint data required for the map.

### $S_k$: rigid, but with no 1-ready root

The other source is the fixed template

    S_1=(222,76,16,116,226,199,132,14)

with vertex 0 replaced by the increasing block $x_0,\ldots,x_{k-1}$.
Retain names 1 through 7 for the other template vertices. The block
vertices have indegrees $1,2,\ldots,k$; the remaining indegrees are

\[
(k+3,k+5,k+2,k+2,2,k+4,k+3)
\tag{8}
\]

at old names 1 through 7. Thus $x_0$ is the unique indegree-one vertex,
and its sole in-neighbour, old vertex 5, is fixed by every automorphism.
All block vertices are now fixed by indegree. Old vertices 2 and 6 have
unique remaining indegrees, and each equal-degree pair $\{3,4\}$ and
$\{1,7\}$ is fixed pointwise by its internal oriented edge. This proves
that $S_k$ is rigid for every $k\ge1$.

Any ready root must have indegree one, so only $x_0$ is a candidate. After
deleting it, the parent has the two distinct Hamilton paths

\[
5,x_1,\ldots,x_{k-1},1,2,4,6,7,3,
\]
\[
5,x_1,\ldots,x_{k-1},1,2,4,7,3,6.
\]

The block is empty when $k=1$. Every displayed arc follows directly from
the template and homogeneous block rule. Both paths start at 5. For an
indegree-one root, its readiness gap equals the parent's start count at its
sole in-neighbour. Hence this gap is at least two, and $S_k$ has no 1-ready
root.

### Exact size-dependent multiplicity lower bounds

Insert arbitrary transitive singleton prefixes and suffixes around either
strong source. Uniqueness of the SCC word prevents collisions, including
between the two sources, whose strong-component sizes differ. With
$t_+=\max(t,0)$ this gives, for all $n\ge1$,

\[
N_n(h_k)\ge(n-k-6)_++(n-k-7)_+,
\]
\[
L_n(h_k)\ge n!\big[(n-k-6)_++(n-k-7)_+\big].
\tag{9}
\]

For $n\ge k+8$ the bracket is $2n-2k-13$. Cumulatively, the ordinary
count through size $n$ is at least $(n-k-6)_+^2$. These bounds are
unconditional on the prime factors of $h_k$: both cores are rigid.

In this particular constructed subfamily, exactly **one** class at each
order $n\ge k+8$ is 1-ready: $R_k\triangleright TT_{n-k-8}$. All padded
$S_k$ classes and all padded $R_k$ classes with a nonempty prefix fail
readiness by Section 3. Thus the scalar fibre gains linearly many classes
while this usable marked subfamily contributes one per order.

## 7. An exact quantitative hostile and a complete stopping gate

At $k=1$, the complete inherited strong-ready inverse finds precisely one
strong 1-ready class of count 221, namely $R_1$ of order 9. Therefore (4)
gives

\[
J_n(221)=0\ (n\le8),\qquad J_n(221)=1\ (n\ge9).
\tag{10}
\]

The complete strong inventory in this packet has 29, 22, 7 ordinary strong
classes at count 221 in orders 8, 9, 10, respectively. Its four maximal
two-component words use the $13\times17$ factorization. Consequently

\[
N_8(221)=29,\qquad
N_n(221)=L_n(221)/n!=2n^2+20n-262\quad(n\ge9).
\tag{11}
\]

Every class is rigid, also directly covered by (7), since the least prime
factor of 221 is 13. Exactly one of the 80 order-nine classes is ready.
The proportion of ready classes (ordinary or labelled) is exactly
$1/(2n^2+20n-262)$ for every $n\ge9$.

The inherited all-order strong-ready exclusion at 1997 is even sharper:
Section 3 extends it to **every** tournament of count 1997, including weak
ones. Yet the inherited strong ready $Q_1$ of count 665 has order 11; the
injective two-birth map produces a rigid strong 1997 class of order 13.
Singleton padding therefore gives

\[
N_n(1997)\ge n-12,\quad L_n(1997)\ge n!(n-12)\quad(n\ge13),
\qquad J_n(1997)=0\quad\text{for all }n.
\tag{12}
\]

Thus no threshold of ordinary fibre growth can force this particular ready
interface. The exact positive bridge is (4) followed by (6), retaining the
actual strong core and its unique starting-word inverse. The universal
existence of suitable interfaces for all odd targets remains open; other
constructor interfaces are not excluded by (12).

The incoming checkpoint `b305e1d03f` supplies precisely such an escape:
a same-count source at 1997 with indegree one and defect three, followed by
actual 5993 and 5995 branches using source replacements and gap transports.
See the [incoming synthesis, section 2](../../berggren_target_return_20261111/SYNTHESIS.md)
and [defect-three construction](../../berggren_target_return_20261111/reset/DEFECT_THREE_ESCAPE.md).
These are compatible with (12): they change the constructor interface rather
than producing a missing 1-ready root at 1997. The triangle filtration in
that synthesis supplies another exact sidecar: equal path totals can have
different extension-gcd lifetimes depending on zero-gap adjacency.

## Reproduction and scope

Run `python3 -B 05-knowledge/results/spectrum_multiplicity_signed_tables_20261111/multiplicity/rigidity.py`
from the repository root. The structural receipt pins the four inherited
sources, checks the degree identities and the two explicit parent paths at
$k=1,2,7,16$, and reads the exact 221 core-word polynomial from the compiled
inventory. It performs no tournament enumeration, path-count computation,
or automorphism search. The all-order proofs use the template identities,
ordered SCC factorization, and the named inherited classification; the finite
controls are not a substitute for those arguments.
