# The seven-fold fork on arbitrary modules, and its prime obstruction gates

Use the explicit five-edit route from [Physical forks](../chanpak_ternary_buckets_20261122/PHYSICAL-FORKS.md). It has one uniform growing transitive block of size $s\ge1$, and exterior hosts with heights
$$
280u+7,\quad 780u+45,\quad468u+45,\quad312u+67,\quad840u+19,\quad840u+21,\qquad u=2^s.
$$
The center also has an upper child at $840u+23$. The route's fifth host is the lower child. The indicated arcs are all in the fixed exterior.

## Every exterior works verbatim with an arbitrary replacement

Replace the uniform block by any nonempty tournament $M$, keeping all block-to-exterior incidences. Write $D=H(M)$ and let $Q$ count its two ordered path slots, either of which may be empty. Then the seven literal heights (six on the path and the upper child) are
$$
\begin{gathered}
280Q+7D,\quad780Q+45D,\quad468Q+45D,\quad312Q+67D,\\
840Q+19D,\quad840Q+21D,\quad840Q+23D.
\end{gathered}
$$
In particular the first is $h=7(40Q+D)$, the sixth is $3h$, and the two marked reversals give $3h\pm2D$.

Proof: the inherited exterior permits at most four occupied slots. For each fixed exterior, deleting and reinserting those slots expresses the count as $\sum_{j=1}^4 b_jP_j(M)$, where $P_j$ counts ordered nonempty path covers in $j$ slots, and the coefficients are fixed nonnegative integers. On a transitive block of size $s$, $P_j=j!{s\brace j}$. The values on sizes 1 through 4 give a lower triangular matrix with nonzero diagonal; hence an identity $A2^s+B$ there forces the functional $AQ+BD$ on every module, since $P_1=D$ and $P_2+2P_1=Q$. This argument also applies to the upper child, with the same four-slot bound since its edit stays in the exterior and does not affect any block incidence. The exterior arc edits never touch the replacement module. Finally $Q/D\ge2$ (the two placements of a full path). Divide all heights by $D$ and substitute $u=Q/D\ge2$ in the original linear inequalities; these hold on the whole real interval $u\ge2$, as can be checked coefficient-by-coefficient. Thus the intermediate heights stay between $h$ and $3h$, not just on transitive inputs. The edits give literal unit children only when $D=1$, equivalently when $M$ is transitive; this is not a claim of unit edits for general $M$.

## Exact arithmetic filters on the prime-cofactor lane

For a minimal missing odd height larger than 21, a possible single-seven obstruction has the form $7p$ with $p\ge5$ prime and $p\ne7$: if the part coprime to seven is composite and has a divisor $a\ge5$, split off $7a$ and the remaining factor, choosing $a$ proper; if the cofactor is a pure power of three of exponent at least two, use $63$ and powers of $3$. The usual inherited exclusion of two or more sevens uses the attained heights $49,343$. Conversely this multiplicative reasoning by itself does not imply $7p$ when $p$ is prime. The fork therefore supplies new literal $7p$ witnesses whenever $40Q+D=p$ is prime.

For the original transitive module, this cofactor is $5\cdot 2^{s+3}+1$. Odd $s$ gives a multiple of 3. For even $s\ge0$, write $p_s=5\cdot 2^{s+3}+1$; one has the exact base-3 certificate
$$
p_s\text{ is prime}\quad\Longleftrightarrow\quad3^{(p_s-1)/2}\equiv-1\pmod{p_s}.
$$
If prime, $p_s\equiv5\pmod{12}$, so 3 is a quadratic nonresidue by reciprocity. Conversely, the displayed congruence forces every prime divisor $r$ to have $2^{s+3}\mid\operatorname{ord}_r(3)$, hence $r\ge 2^{s+3}+1>\sqrt{p_s}$. No composite integer can have that property.

For general replacement of order $n$, the disjoint odd-cycle-packing formula gives
$$
 D=\sum_{\mathcal C}2^{|\mathcal C|},\qquad
 Q=\sum_{\mathcal C}2^{\,n-\sum_{C\in\mathcal C}(|C|-2)}.
$$
The second identity assigns each odd cycle or uncovered singleton to either ordered slot. Reducing the summands modulo 3 gives $Q\equiv(-1)^n D\pmod3$. Thus for $p=40Q+D$ to be prime it is necessary that $n$ be even and $3\nmid D$. Also $\gcd(D,40Q)=1$ is necessary because $p>D,40Q$. For even $n\ge4$, each summand of $Q$ is divisible by 8: the empty packing has exponent $n$, a singleton-cycle packing has an odd number of uncovered vertices and exponent at least 3, and any packing with at least two cycles has exponent at least 4. Consequently
$$
 p\equiv D\pmod{320}\quad(n\ge4\text{ even}).
$$
For a fixed filler $V$, adding $s$ ordinal singleton vertices leaves $D$ unchanged and multiplies $Q$ by $2^s$, giving entire affine lanes $p_s=40Q(V)2^s+D(V)$. Only $s+|V|$ even can yield a prime. No infinitude of primes of any such shape is asserted.

For an exact inverse on this family, $Q\ge2^n$ follows from the empty packing and $Q\ge(n+1)D$ by cutting each Hamilton path at each of its $n-1$ internal positions into two ordered nonempty paths, and adding the two empty-slot cases. Hence for a given prime $p$ search only orders $n\le \lfloor \log_2((p-1)/40)\rfloor$, and only $D\le p/(40(n+1)+1)$ in the indicated congruence. This is not an effective inverse for primes outside the family.

## Six prime witnesses off the transitive ray

In each row take the six outgoing adjacency masks numbered from zero. The indicated $D,Q$ are reproducible by ordinary Hamilton-ending subset DP and $Q=\sum_S H(M[S])H(M[V\setminus S])$; the prime status can be certified by trial division up to the square root.

| Outgoing rows | $D$ | $Q$ | prime $40Q+D$ |
|---|---:|---:|---:|
| (38,60,56,49,33,0) | 11 | 200 | 8011 |
| (22,60,56,49,32,1) | 19 | 256 | 10259 |
| (28,53,56,50,32,1) | 23 | 272 | 10903 |
| (14,60,56,48,33,1) | 23 | 296 | 11863 |
| (26,44,57,48,34,1) | 31 | 328 | 13151 |
| (22,44,56,49,34,1) | 37 | 376 | 15077 |

Thus all six products $7p$ are exact heights by this one exterior, with no need to infer them from congruences or divide an already realized count. The six instances are new to this chart, not a claim that these integers were absent from the repository's full prior certified spectrum.

The checker reconstructs the four literal endpoint exteriors (parent, lower, center, upper), verifies that each one-vertex-module base host is strong, and directly recounts the four inflated hosts for each row. Replacing a vertex of such a strong host by any nonempty tournament preserves strength: paths out of and back into the replaced vertex can use the common external incidences, including to reach a different vertex of the module. The padded five-edit route need not itself be strong. The checker also verifies the omitted transitive-ray parent at $s=0$: in a singleton-apex host take the contact word OOIOIOIOI, which has count $287=7\cdot41$. No fork edits at $s=0$ are asserted.


## The endpoint and the larger census

At \(s=0\), the original parent alignment is not used. The
nine-contact singleton apex word \(\mathrm{OOIOIOIOI}\) has
literal height 287. With the inherited seven-vertex center
at right address 1 the same two core flips give
\(859,861,863\), so the scalar fork extends to \(s=0\).
This does not claim the five-edit aligned route at that endpoint.

The companion anchored census tests every labelled filler
with a fixed spanning path (every tournament can be relabelled
in this way). It has 1727 distinct \((D,Q)\) profiles at order
eight, giving 276 distinct prime cofactors, all off the dyadic
rail. Order six gives seven more; the ranges of prime labels at
the two orders are disjoint. The six examples in the table
are representative, not a complete order-six list. The
maximum order-eight prime tag 354581 has filler masks
\((226,100,137,19,39,76,156,58)\) and profile \((661,8848)\).
Its parent and three children are physically recounted in
the replay, without assuming the larger chart is onto.
