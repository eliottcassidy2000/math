# Transfer from LRC(14): divisor clocks, exact witnesses and run complexity

**Proved.** An LRC-style divisor sieve removes 23,004 of the 27,720
exponent classes of the remaining defect-14 family from consideration
for a *least* Chan–Pak counterexample. Four explicit tournaments
and exact compositeness descents imply that any least counterexample
on this branch must have exponent $j\ge365$. Also, if either of two inherited +14 relays hits
this family at unbounded height, the number of runs in its binary
address must grow at least logarithmically in the exponent. Neither
statement closes the family.

## What transfers from the runner work

The [LRC(14) covering reduction](../../../07-reflections/the-disproof-search-builds-the-proof-lrc14-covering-reduction-mac-mini-S3.md)
uses a finite list of rational candidates $1/q$. Missing all speeds
divisible by $q$ immediately certifies that candidate. The
[filter/witness analysis](../../../07-reflections/the-filter-witness-asymmetry-and-why-no-fixed-template-closes-lrc14.md)
warns that reducing to a set of residues does not itself supply the
witness at each surviving height. The transfer here uses two separate
objects: an exponent clock with smaller-count descents, and an exact
Hamilton witness construction. No implication between the conjectures
is used.

## The exponent clock

Write $N_j=38\cdot 2^j+1$. For any odd modulus $q$ coprime to 38,
the test $q\mid N_j$ depends only on $j\bmod\operatorname{ord}_q(2)$.
Here are useful exact clocks:

| $q$ | Period | Residue for divisibility |
|---:|---:|---:|
| 3 | 2 | 0 |
| 5 | 4 | 3 |
| 11 | 10 | 1 |
| 23 | 11 | 8 |
| 67 | 66 | 22 |
| 199 | 99 | 76 |
| 463 | 231 | 147 |
| 3697 | 1848 | 1769 |
| 4621 | 4620 | 3467 |
| 49 | 21 | 13 |

For $j\ge12$ each listed factor and its cofactor, when present,
are strictly smaller than $N_j$ and neither is a hole: indeed
$N_j>21\max q$. If $N_j$ is the least additional counterexample,
both factors would therefore be attained. Ordinally joining
their tournaments multiplies their path counts, a contradiction.
In particular this is a descent under minimality, not a
standalone witness for every term in those residues.

The joint period is 27,720 and its complement has 4,716 classes
(17.013 percent, rounded). Among those classes 1,392 allow
$j\equiv1\bmod3$, exactly where $7\mid N_j$. Outside those classes
a least counterexample would have to be prime, since its proper
factors would be 7-free. In the 1,392 classes it would have to be
$7p$ with $p$ prime: 49 is sieved out, and a composite $N_j/7$
other than 9 admits a proper nontrivial factorization with one
factor different from 3, avoiding both forbidden cofactors.
The value 63 is not in this large-height range.

A useful exact transfer of the LRC clock idea: for $j>i$,
$N_i\mid N_j$ iff $j\equiv i\bmod\operatorname{ord}_{N_i}(2)$.
Moreover $N_j/N_i$ cannot be 7 or 21 since $N_j\equiv N_i\equiv
1\bmod38$. Thus a verified target on the ray can seed another
clock for the least-counterexample reduction, without needing
that target to be prime. The new clock alone still does not
realize its cofactors.

## Exact witnesses and descents below exponent 365

Use the inherited [F1 bridge](../chanpak_moonshine_prime_blocks_20261118/POLYGON_COMPLETION_AND_PRIME_COUNTS.md).
Its address pair $(b,e)$ builds two one-apex tournaments with a
single exceptional arc from the left apex to the right apex, and
reads $be+F(b-1)F(e-1)-[b,e\text{ odd}]$. Here $F(0)=1$ and for
positive arguments $F(2u+1)=F(u)$,
$F(2u)=F(u)+u+[u\text{ odd}]$. The following literal addresses
are on the clock's surviving side:

| $j$ | $N_j$ | $(b,e)$ |
|---:|---:|---|
| 29 | 20401094657 | (2,8283070098) |
| 33 | 326417514497 | (2,133767957234) |
| 37 | 5222680231937 | (2,2058347761263) |
| 45 | 1337006139375617 | (71,17552393816392) |

These are actual witnesses, not a use of minimality or a prime
claim. They also give a separate checkpoint beyond residues.
The earlier corrected certificate through $2^{32}-1$ covers
$j\le26$. The replay then checks every surviving $27\le j<365$ not covered by
the table. Put $Q=N_j$ when $7\nmid N_j$ and $Q=N_j/7$ otherwise.
A certificate $a^{Q-1}\not\equiv 1\pmod Q$ proves $Q$ composite,
which rules out this height for a least counterexample by the exact
factorization criterion above. For each checked term a small base
suffices, without factoring the term. Hence a least additional
counterexample on this branch must have $j\ge365$. This is **not** an
unconditional spectrum census through that height, and it does not
assert nonrealizability at $j=365$. The [old cap and literal bridge proof](https://github.com/eliottcassidy2000/math/blob/main/session203/coverage-recertification/README.md)
and the current replay have different scopes.

## The inherited +14 relays cannot use a bounded number of runs

The [fourteen-axis single-arc and original two-contact relays](../fourteen_axis_20261117/README.md)
read a normalized nonzero binary address $z$ as
$$
R(z)=C+A(z+1)+gF(z)+d,\qquad d\in\{0,14\},
$$
with $(A,g,C)=(13190,37,-5884)$ or $(10542,175,-3374)$.
The exact count at $z=0$ is irrelevant for an unbounded family.
Use the [endpoint/low-bit identity](../inverse_fences_landing_20261116/INVERSE.md)
in the form
$$
F(LU+r)=F(U)+ZU+F(r)-1+\epsilon[U\text{ odd}],
$$
where $L=2^\ell$, $0\le r<L$, $Z$ is the reverse-complement of
the $\ell$ low bits, and $\epsilon=[r<L/2]$. For $\ell=0$ this
means $Z=0$ and no correction. The identity follows by iterating
the displayed recursion for $F$.

Suppose the run lengths of $z$ in ascending order are
$a_1,\ldots,a_s$, with partial sums $S_i$. If there is an index
$i<s$ with $a_{i+1}\ge S_i+20$, call those at least $a_{i+1}$
long. Let $P$ be the length of the prefix before the first long run,
and $\ell$ the length of the suffix after the last long run. Both
consist only of short runs, so $P+\ell\le S_i$. Put $D$ equal to
the length from the start of the first long run to the end of the
last long run, and $U=\lfloor z/2^\ell\rfloor$.
If the first long run has bit $b$ and length $\sigma$, and
the prefix has binary value $p$, put $v=p+b$.
Then $1\le v\le 2^P$ and
$$
\left|U/2^D-v\right|\le 2^{-\sigma}.
$$
All normalized suffixes after that first long run contribute at
most $2^{-\sigma}$, also for a run of ones. The leading bit of
$z$ is 1, which guarantees $v>0$.

Let the last long run have length $\tau$, and put $\delta=1$ if
its bit is 0 and 0 otherwise. The identities
$F(x2^\tau)=F(x)+x(2^\tau-1)+[x\text{ odd}]$ for $x>0$
and $F((x+1)2^\tau-1)=F(x)$ imply
$$
|F(U)-\delta U|\le2U\,2^{-\tau}.
$$
Here $1\le F(x)\le x+1$ follows by elementary induction.
The exact suffix formula and $D\ge\max(\sigma,\tau)$ now give
$$
\left|R(z)/2^D-vB\right|
\le (2(A+g)L+2g(v+1)+|C|+14)\,2^{-(P+\ell+20)},
\qquad B=AL+g(Z+\delta).
$$
For example, bound the suffix constant by
$(A+g)L+|C|+14$, its linear part by $(A+g)L$, and
$U/2^D$ by $v+1$. For either relay the right-hand side plus
$2^{-D}$ is less than $1/2$.

If $R(z)=N_j$, we obtain $|38\cdot2^{j-D}-vB|<1/2$.
When $j<D$ this is impossible because $vB\ge A>19$.
Otherwise both entries are integers and must be equal. In that
case the odd part of $B$ divides 19. But
$$
B/L\in [A,A+g],
$$
and neither interval contains a power of 2 or 19 times a power
of 2 (including negative exponents). This gives the contradiction.

Therefore an exact hit has no gap of the form
$a_{i+1}\ge S_i+20$. Induction on the cumulative lengths gives
$\operatorname{bitlen}(z)<20(2^s-1)$. On the other hand
$N_j=R(z)<(A+g)2^{\operatorname{bitlen}(z)}$ forces
$\operatorname{bitlen}(z)\ge j-8$ for either pair of constants.
Thus any relay hit obeys
$$
j-8<20(2^s-1).
$$
This rules out every bounded-run template for those relays, even
if its individual lengths are allowed to depend arbitrarily on
$j$. It does not forbid a sufficiently complicated changing
address, a different relay donor, or an independent tournament
construction. This is the precise analogue of the earlier
runner warning: a clock or a bounded menu filters a problem,
but the certificate that survives may require growing complexity.

Run python3 lrc14_defect_replay.py for the divisor period, literal
addresses, and small profile checks.
