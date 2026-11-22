# Three apex contacts, Mersenne mantissas and perfect-number neighbors

**PROVED.** Every integer
\[
             1+2^j(2^{2k}-1),\qquad j\ge1,\quad k\ge2,
\]
is the Hamiltonian path count of an explicit strong tournament of order $j+2k+2$. Within the single growing-block construction, all prime-square mantissas are also classified below. A different pair of marked tournaments gives both odd neighbors of half of every even perfect number greater than 6. A Wieferich prime is distinguished by an exact repetition of the two endpoint residues in the same physical width phase.

## A two-contact face with opposite single changes

Use a transitive spine and one apex. Number the spine positions by descending binary exponents $x,x-1,\ldots,0$. The apex beats precisely the marked positions. If the positive marks form a set $S$ with pairwise distance at least two, its Hamilton count is
\[
 H(S)=1+\sum_{t\in S}2^t-\sum_{\substack{t,u\in S\\t>u}}2^{t-u-1}.
\]
Indeed a mark at $t$ would contribute $2^t$ if every later position were unmarked. A later mark at $u$ removes the indicated contribution since the difference is at least two. This is the isolated-position specialization of the inherited exact weighted-inversion count. Adjacent marks have a different weight. [Single-apex seam, §2](../chanpak_graceful_ternary_20261115/GRACEFUL.md)

Fix $j\ge1,k\ge2$ and choose
\[
                 x=j+2k,\qquad y=j+k,\qquad z=k-1.
\]
Keep the mark at $x$ and allow marks at $y,z$ independently. All four words start in 1 and end in 0, so a directed cycle through the apex and the spine proves strong connectivity. Their counts make one two-dimensional face of the labelled arc-flip metagraph. For $h=1+2^x$ and $a=2^y-2^z$, the four corners are

| Extra marks | Count |
|---|---:|
| neither | $h$ |
| $y$ | $h+a$ |
| $z$ | $h-a$ |
| both | $h-2^j$ |

Indeed $x-y-1=z$, $x-z-1=y$, and $y-z-1=j$. The two individual changes cancel in the aggregate. Their joint contribution measures the separation between the two new contacts and is exactly $-2^j$. This produces the claimed all-height count in one fixed prescription. For $2k=2^b$ the odd part $2^{2k}-1$ is the product of the first $b$ Fermat numbers, by repeated difference of squares. This statement includes composite Fermat numbers.

Among triples of isolated exponents $x>y>z>0$, opposite single changes happen **if and only if** $x=y+z+1$. The cancellation equation is $2^y+2^z=2^{x-z-1}+2^{x-y-1}$. Uniqueness of binary expansion identifies the two unequal exponents in decreasing order. Thus the all-height prescription recovers exactly the even gap $x-j=2(z+1)$ in this three-mark mechanism. For an odd gap this specific cancellation cannot occur.

An equivalent factor-level description chooses $e=2^k+1,\ f=2^k-1$ in the inherited single-block inverse. In the canonical retained-root encoding, the target $1+(2^{2k}-1)2^j$ starts with $2k-1$ zeros followed by $j$ ones and has defect $2^{2k}-2$. Thus these rays have unbounded defects, and the construction applies at every height of each ray. Other targets at those same defects can have different encodings. [One-block divisor criterion and prime obstruction](../chanpak_moonshine_prime_blocks_20261118/ONE_APEX_PRIME_SLOPE.md) · [Canonical defect recurrence](../chanpak_moonshine_prime_blocks_20261118/LOW_DEFECT_DESCENT.md)

## Exactly which prime squares pass a single growing block?

For the single growing zero block, let $F(0)=1$, $F(2v+1)=F(v)$ and $F(2v)=F(v)+v+[v\text{ odd}]$ for $v>0$. The inherited criterion requires odd factors $q=ef$ satisfying
\[
                   F(e)+F(f-1)=e+1.
\]
It gives a physical word at every height $j\ge1$. The previous source proves that $q>1$ cannot pass via a unit factor and proves a stronger prime obstruction for two affine growing blocks. In the single-block criterion, the factors also satisfy $(e+1)/2\le f\le2e-1$: $1\le F(e)\le(e+1)/2$ and $(f+1)/2\le F(f-1)\le f$ when $f>1$. [One-apex prime slope](../chanpak_moonshine_prime_blocks_20261118/ONE_APEX_PRIME_SLOPE.md)

**Theorem.** For an odd prime $p$, the single-block construction gives the full ray $1+p^2 2^j$ precisely when $p=3$ or $p=2^t+3$ for some $t\ge3$.

**Proof.** The only surviving ordered factorization is $e=f=p$. Writing $p=2r+1$ gives $2F(r)=r+2-[r\text{ odd}]$. If $r=2s>0$ this requires $F(s)+[s\text{ odd}]=1$, impossible for positive $s$ since a positive even address already has count greater than one. If $r=2s+1$ we need $F(s)=s+1$. The solutions are zero and the powers of two at least two. Indeed no positive odd $s$ satisfies equality by the odd bound. For $s=2u$ equality requires $F(u)=u+1$ if $u$ is even and $F(u)=u$ if $u$ is odd. In the latter case the odd bound forces $u=1$, which works. Downward induction proves the classification. Substituting $p=4s+3$ finishes. Conversely the divisor criterion supplies the required physical realizers by taking the prefix to have binary value $p$ and the reverse-complement suffix to have value $p-1$. ∎

In particular $1+49\cdot2^j$ never arises from this fixed-block ray construction. The admitted prime-square examples start at $9,121,361,4489$ (primes $3,11,19,67$). The exact classification concerns an entire growth rule. Individual counts can arise from other addresses or tournaments.

## Full and half neighbors of an even perfect number

Put $M_r=2^r-1$, $E_r=2^{r-1}M_r$. The even perfect numbers are exactly those $E_r$ for which $M_r$ is prime. The formulas below also make sense at composite exponents. [Euclid, IX.36](https://mathcs.clarku.edu/~djoyce/java/elements/bookIX/propIX36.html) · [GIMPS, mathematical background](https://www.mersenne.org/various/math.php)

1. Every even perfect $P$ has realizable $P-1$ and $7(P-1)$. Indeed $E_r-1=(2^{r-1}-1)(2^r+1)$ and the signed-height constructions give both factors and the sevenfold first factor at exponents $r\ge5$. Ordinal joining multiplies the counts. At $P=6,28$, the required counts are $5,35,27,189$: $189=3\cdot63$. The same proof gives both counts for every $r\ge5$. [Signed height selectors and small witnesses](../chanpak_six_signed_height_pairs_20261118/README.md)
2. For every $r\ge4$, the odd numbers $E_r/2-1,E_r/2+1$ are both physical counts on order $2r$. In Family I of the two-port construction take $m=r-3,\ k=m+1$. The center $8\cdot2^{m+k}+2^{m+1}-2^{k+1}$ is exactly $2^{r-2}(2^r-1)$, and replacing two marked vertices of the seven-vertex core yields order $5+m+k=2r$. At $r=3$, the claims at $P=28$ are 13 and 15, both realizable at order five with outgoing masks $(18,8,3,5,14)$ and $(16,9,3,5,14)$ respectively. At $P=6$ its half-neighbors are even. These are two specified physical sources; an arc reversal between them is not assumed. [Two-port neighbors, Family I](../chanpak_moonshine_prime_blocks_20261118/TWO_PORT_NEIGHBOR_FAMILIES.md)
3. For every even $r\ge4$, the full upper neighbor $E_r+1$ has the three-contact construction with $j=r-1,k=r/2$. A uniform construction for all odd prime exponents above three remains open here; $E_2+1=7$ is the first global hole. An arithmetic restriction on a hypothetical *least new odd hole* of this shape follows from
\[
               8(E_r+1)=(2^{r+1}-1)^2+7.
\]
An odd prime divisor other than 7 must be a quadratic residue modulo 7 by quadratic reciprocity for $-7$. The factor 7 occurs iff $r\equiv2\pmod3$ and then exactly once, with $(E_r+1)/7\equiv1\pmod7$. Thus for prime $r>3$ a least new hole of this shape is either prime $\equiv2\pmod7$ (when $r\equiv1\pmod3$), or seven times a prime $\equiv1\pmod7$ (when $r\equiv2\pmod3$). Indeed any other factorization can be grouped into two smaller allowed numbers: 3 is a nonresidue modulo 7 and the factor 7 appears at most once. Ordinal joining gives the contradiction. This restricts absolute targets without asserting that a full upper neighbor is absent.

**How quickly the two scales separate.** If $p<q$ are any two Mersenne-prime exponents and $d=q-p$, direct subtraction gives
\[
 E_q=4^d E_p+2^{q-1}(2^d-1).
\]
After $p=2$ the exponents are odd primes, so $d\ge2$ and $E_q-E_p>15E_p$. This lower bound applies in particular whenever a next even perfect number exists. Meanwhile the two-port realization at $r\ge4$ uses exactly $2r$ vertices, one more than the binary length $2r-1$ of $E_r$. Moving between two such exponents adds just $2d$ vertices to the construction. The arithmetic candidate sets have different size cutoffs. For $X\ge3$, $2^{2^b}+1\le X$ means $b\le\lfloor\log_2\log_2(X-1)\rfloor$, while $2^p-1\le X$ means $p\le\lfloor\log_2(X+1)\rfloor$, restricted to prime exponents for a Mersenne-prime search. These are exact cutoffs on candidates. The gaps between exponents actually giving primes remain a separate problem. The previous small-spectrum microcosm records how the first 32 odd values are cut into finite runs by 7 and 21, while the face above produces an all-height count by tracking the two individual arc changes. [THM-5271, flip potential and small-spectrum microcosm](../../../01-canon/theorems/THM-5271-the-count-is-a-potential-on-the-flip-metagraph-the-seam-ladder-implies-chan-pak-the-microcosm-3-1-6-1-21-the-square-sum-cycle-at-32-and-the-order-10-near-join-reach.md) · [GIMPS, Mersenne exponents and known perfect numbers](https://www.mersenne.org/primes/)

## Two explicit offsets at the full perfect-number scale

**Theorem.** For every integer $r\ge5$, the count $E_r-5$ has a strong tournament on $2r$ vertices. For every $r\ge7$, $E_r-9$ has one on $2r-1$ vertices.

Both constructions replace vertices 1 and 9 of a ten-vertex tournament by transitive blocks of sizes $m,k$. Rows give outgoing masks in zero-based vertex order. The columns $c_{ij}$ count compressed directed words through the eight unmarked vertices once each and the two marked vertices $i,j$ times.

| Offset | Outgoing rows | $(c_{11},c_{12},c_{21},c_{22})$ | $(m,k)$ |
|---|---|---|---|
| $-5$ | 990,1020,1016,976,928,137,688,768,608,32 | $(491,294,210,128)$ | $(r-4,r-4)$ |
| $-9$ | 998,1020,1016,977,929,456,656,768,576,32 | $(755,526,368,256)$ | $(r-6,r-3)$ |

**Proof.** The finite coefficients in the table, with all other coefficients zero, have the following exact enumeration. Start a word at any vertex. At its last letter extend by any outneighbor whose label is marked or still unused. Record one contribution with its two marked multiplicities whenever all labels have been used, then continue exploring legal extensions. This terminates because every segment that visits a mark twice must use at least one unmarked vertex between those visits. The executable calculation is in the accompanying check.

Each visit of multiplicity one admits one filling of its block and each of multiplicity two admits $2^m-2$ or $2^k-2$ ordered fillings. This gives the polynomial
\[
 c_{11}+c_{12}(2^k-2)+c_{21}(2^m-2)+c_{22}(2^m-2)(2^k-2).
\]
For the $-5$ row it reads $128\cdot2^{m+k}-46\cdot2^m+38\cdot2^k-5$. For the $-9$ row it reads $256\cdot2^{m+k}-144\cdot2^m+14\cdot2^k-9$. Substituting the stated widths gives the desired counts. The quotient cycles in the check visit every vertex and can be lifted by running through each block consecutively. The inflated tournaments are therefore strong. ∎

## The first repeated width phase detects Wieferich primes

Fix an odd prime $\ell$ and the **two actual counts** $h_\pm(r)=E_r/2\pm1$ of the two-port construction, $r\ge4$. Their simultaneous residues $h_-\equiv-1,\ h_+\equiv1\pmod{\ell^a}$ occur exactly when
\[
              r\equiv0\pmod{d\,\ell^{\max(0,a-s)}},\qquad
         d=\operatorname{ord}_{\ell}(2),\quad s=v_\ell(2^d-1).
\]
Indeed the common center is $2^{r-2}(2^r-1)$. If $d\mid r$, the binomial expansion of $(2^d)^{r/d}$ gives valuation $s+v_\ell(r/d)$. Raising a unit $1+\ell^s v$ to $\ell$ adds exactly one to the valuation of its difference from one, while raising it to an exponent coprime to $\ell$ leaves that valuation unchanged. The exceptional case $s\ge2$ is equivalent to $2^{\ell-1}\equiv1\pmod{\ell^2}$ because $(\ell-1)/d$ is prime to $\ell$. Thus the exact phase spacing in $r$ stays fixed on passing from $\ell$ to $\ell^2$ precisely for base-two Wieferich primes. Actual tournament orders are $2r$, so their spacings are doubled. The three-mark face independently gives the same integral test on its final count: for fixed $j$, that count is $1\pmod{\ell^a}$ exactly when the displayed period divides $2k$. [Crandall, Dilcher and Pomerance, *A search for Wieferich and Wilson primes*](https://math.dartmouth.edu/~carlp/PDF/paper111)

For a Mersenne prime $\ell=2^p-1$ one has $d=p,s=1$. For a Fermat prime $\ell=2^{2^b}+1$, one has $d=2^{b+1},s=1$ because $2^{2^{b+1}}-1=\ell(\ell-2)$. Both force multiplication of the physical period by $\ell$ at the first lift. The primes 1093 and 3511 instead have $(d,s)=(364,2)$ and $(1755,2)$ respectively. No frequency claim about perfect numbers, Fermat primes, or Wieferich primes enters these predictions. The mobile two-port pair and the single-apex face have different vertex operations. Only the stated counts and width phases have been transported.

Run the checks from the repository root with

    python3 05-knowledge/results/chanpak_perfect_wieferich_seams_20261122/check.py

The checker counts all four corners by subset dynamic programming at the indicated small orders, enumerates the divisor criterion for odd slopes through 2001 and prime squares through $p=500$, and computes the two-port and prime-power examples independently. All-height conclusions use the proofs above.
