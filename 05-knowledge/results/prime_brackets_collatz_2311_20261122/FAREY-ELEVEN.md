# At 11, the Farey slope lives on a one-digit-delayed graph

Use the closed form for the canonical Farey spine $n>0$, $4\mid n$, from the two-clock packet. Write $S=h/145=2^n+1$. The integer read
$$V=117t-85680S+986112=1024000\cdot3^n$$
removes the binary clock from the left slope. A pair $(S,t)$ of prescribed modular coefficients is equivalent to prescribing $(S-1,V/1024000)$ at 11.

**Exact 11-adic image.** For $r\ge1$ and a fixed residue $n_0\in\{4,8,12,16,20\}$ modulo 20, the binary coefficient $S$ ranges over all residues congruent to $1+2^{n_0}$ modulo 11 at precision $11^r$, uniquely in the parameter $j=(n-n_0)/20$ modulo $11^{r-1}$. The corresponding $V$ is forced. More precisely, for distinct exponents in the same class,
$$v_{11}(S(n)-S(n'))=1+v_{11}(j-j'),\qquad v_{11}(V(n)-V(n'))=2+v_{11}(j-j').$$
Thus $S$ modulo $11^r$ already fixes the ternary component $V$ modulo $11^{r+1}$. Within a class the ternary projection alone covers every residue congruent to $1024000\cdot3^{n_0}$ modulo 121 at each precision $11^r$ for $r\ge2$. The full profile modulo $11^r$ has exactly $5\cdot11^{r-1}$ realizable pairs $(S,t)$, not independent coordinates.

Proof: $v_{11}(2^{20}-1)=1$ and $v_{11}(3^{20}-1)=2$; the two valuation equalities are LTE, since both $1024000$ and the base powers are units. The induced maps attain every claimed value by injectivity and equality of finite cardinalities, and the five binary classes are disjoint since the order of $2^4$ modulo 11 is five. Nonnegative representatives $j$ suffice at every finite level. Also $85680$ is a unit modulo 11, so $v_{11}(t(n)-t(n'))=v_{11}(S(n)-S(n'))$ within each class. Consequently the left slope itself can serve as the advancing clock, but must not be confused with its ternary component.

The formal equation $1000\cdot3^n=963$ has $n=3\pmod 5$ at 11 and has **no** solution modulo 121: $963/1000=49\pmod{121}$, while $3^3=27\pmod{121}$ and $3^5=1\pmod{121}$. The binary obstruction $S=0\pmod{11}$ never occurs at $4\mid n$. This rules out 11 as a cause of a Farey axis-gcd failure, but still gives an actual blocked lift on its ternary coefficient. It does not exclude any unrestricted E61 flank pair.

As a single-prime analogue of the shared-exponent compatibility at spine 8200, this calculation distinguishes three readings: $S$ and $t$ move by one digit at each advance, but after their exact linear cancellation the surviving ternary reading moves a digit later. This distinction persists to all depths, not merely modulo 121.

## Coupling the actual Euclid target, where the two exponents advance together

If one deliberately varies the spine and Euclid exponent together by $n=p-1$, with $p=1+20m$ and $m>0$, then the prescribed lower target is itself a binary clock on the same branch:
$$N_p=E_p-1=2(2^n)^2-2^n-1=(S-2)(2S-1).$$
Since $S=2$ modulo 11 on this branch, the difference quotient
$2(S+S')-5$ is a unit modulo 11. Hence $S\mapsto N_p$ is an isometry at every finite precision on this branch, and exhaustive counting makes it a bijection from the class $S=2$ onto the class $N_p=0$ modulo 11 at each precision. In particular $v_{11}(N_p)=1+v_{11}(m)$ and two such targets congruent modulo $11^r$ force equal slope $t$ modulo $11^r$ and equal canceled ternary read $V$ modulo $11^{r+1}$. This ties the read to an actual moving Euclid target (including $p=61$), not to an arbitrary right-hand side. It only describes the diagonal choice $n=p-1$; it does not authorize identifying the two exponents in the unresolved spine-four or spine-eight inverse.

## Replay

Run `python3 replay.py`. The finite controls check the affine
coefficient identities on the existing small charts and enumerate
the complete finite clock images through \(11^4\). The proof above
establishes all levels.
