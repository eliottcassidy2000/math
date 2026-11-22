# A fixed arc realizes unit gaps within cube-root distance of every target

**PROVED.** For every sufficiently large real target \(N\), there are strong tournament pairs with path counts \((L_-,L_-+2)\) and \((L_+,L_++2)\), where
\[
L_-<N<L_+,\qquad \max(N-L_-,L_+-N)=O(N^{1/3}).
\]
In each pair **the same specified internal arc of a six-vertex core is reversed**. One can construct both pairs using logarithmically many evaluations of binary-address formulas. In particular the number of distinct single-arc unit seams with lower endpoint at most \(X\) is \(\Omega(X^{2/3})\).

**Two independently selected targets are exact.** The first twin primes above \(2^{32}\), namely **4,294,967,387 and 4,294,967,389**, are joined by one arc reversal on 34 vertices. The first prime above \(10^{12}\), **1,000,000,000,039**, is joined to **1,000,000,000,041** on 38 vertices. Explicit constructions for both appear below. Chan–Pak and the [full seam ladder, THM-5271](../../../01-canon/theorems/THM-5271-the-count-is-a-potential-on-the-flip-metagraph-the-seam-ladder-implies-chan-pak-the-microcosm-3-1-6-1-21-the-square-sum-cycle-at-32-and-the-order-10-near-join-reach.md) remain **OPEN**. Chan and Pak record the spectrum conjecture and credit its formulation to bof and Gordon Royle. [Conjecture 6.2](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf)

## Inheritance and the change in mechanism

The [single-apex contact theorem](../chanpak_apex_flip_local_global_20261122/README.md) has square-root spacing and square-root capacity because its fixed gap forces the two flank addresses within three choices of one another. The [complementary-carry theorem](../chanpak_cone_collision_placement_20261115/STEERING.md) already gives cube-root spacing for scalar counts with two free addresses. The [star product](../chanpak_moonshine_prime_blocks_20261118/STAR_SEAM_PRODUCT.md) retains two tagged path masses, while the [rooted triple splice](../fresh_unit_prime_gluing_20261116/triple_splice/README.md) transports a unit difference with one directed continuation. A [near-twin Farey inverse](../least_hole_near_twins_20261117/README.md) shows why a small-gap equation on two free profiles can itself rule out prescribed inputs.

Here the differing arc is *inside* the tagged core. At an actual pair of cores the unsigned changes vanish on both sides and the signed changes depend only on the two end bits. Hence two independently chosen flank addresses can set the height while a fixed parity choice suffices for the gap. The preserved object is a literal edge in the single-arc flip graph. The resulting pairs need not be consecutive edges on one increasing walk.

## Two freely variable flanks and one audited six-vertex core

Take a tournament \(D\) in ordinal order between two transitive flanks. Add an apex beating the tagged vertices in all three parts. A 1 bit means the apex beats that vertex. The left word is the standard binary representation of \(x\); the right word is the reversed complemented standard representation of \(y\). A zero address denotes an empty flank. When \(x,y>0\), the first vertex is tagged and the last is untagged. The apex can reach every vertex via the left flank (and its own outgoing arcs), and every vertex can reach the apex via the right flank. The tournament is strong.

For \(D\) define \(h=H(D)\), \(k=H(D+\text{apex})\), and \(Q=\sum_{S\subseteq D}H(S)H(D-S)\), with \(H(\varnothing)=1\). To form \(e,o\), count paths on nonempty subsets that **start tagged**, weighted by the complementary Hamilton count, separating even and odd subset size. Let \(e',o'\) be the corresponding counts for paths that **end untagged**. Set
\[
p(0)=0,\qquad p(2z+1)=p(z),\qquad p(2z)=p(z)+z+[z\text{ odd}].
\]
The path count with the two flanks is
\[
k+h(p(x)+p(y))+Qxy+(e'+o')x+(e+o)y
 -(e'-o')[x\text{ odd}]-(e-o)[y\text{ odd}]. \tag{1}
\]

For the derivation, combine consecutive marked blocks using their common apex. An odd cycle using that apex and both blocks uses a nonempty path starting tagged in the first, followed by a nonempty path ending untagged in the second. Equal path-size parity is required. In [THM-002, odd-cycle packing formula](../../../01-canon/theorems/THM-002-ocf.md) its weight two turns this contribution into the product of unsigned sums plus the product of signed sums. Cycles avoiding the apex stay in their individual blocks. The unsigned start sum after appending a deletion block \(D\) is multiplied by \(Q\), with the new starts in \(D\) added with the old deletion count as weight. The signed contribution from an optional path in nonempty \(D\) vanishes, by the slot convolution evaluated at \((-1,1)\). These are precisely the [ordinal arm rules, including the signed term](../chanpak_moonshine_prime_blocks_20261118/ALIGNED_ARMS.md). A left transitive word gives apex count \(1+p(x)\), unsigned start sum \(x\), signed start sum \(-[x\text{ odd}]\); the right word gives the reversed end analogues. Applying the arm rules twice gives (1). When \(x\) is even and positive, its final untagged vertex is literally a width-one returning block in the star product, with the remaining left flank in the earlier arm and \(D\) and the right flank in the later arm.

The two actual cores have common tag mask 22, with vertex 0 the low bit. Outgoing masks and the complete read are

| Outgoing masks | \(h\) | \(Q\) | \((e,o)\) | \((e',o')\) | \(k\) |
|---|---:|---:|---:|---:|---:|
| \((4,17,34,7,13,27)\) | 29 | 320 | (59,70) | (73,89) | 113 |
| \((4,17,42,3,13,27)\) | 29 | 320 | (58,71) | (72,90) | 111 |

The sole changed arc joins vertices 2 and 3. The finite identities in this table use the integer recurrence for paths on \(S\) ending at \(v\): sum the counts ending at predecessors of \(v\) on \(S-\{v\}\), initializing the permitted singleton starts to one. The replay performs this entire certificate directly from the listed rows. Subtracting (1) for the second core from the first gives
\[
2-2[x\text{ odd}]-2[y\text{ odd}].
\]
Thus the unit switch works in either like-parity chamber, with opposite orientation on the two chambers. On mismatched parity the flip preserves the count exactly. Discarding the signed profile would predict the wrong edge on half the input pairs.

Take \(x=2a, y=2b\). On \(a+b\) odd, the second core has the smaller count
\[
L(a,b)=1280ab+353a+287b+29(p(a)+p(b))+140,\qquad a,b\ge1. \tag{2}
\]
The first core has exactly \(L(a,b)+2\). The flanks in (2) are genuinely independent.

## Explicit two-sided bound and decoder

The following constants apply to integer \(N\) with \(z=\lfloor\sqrt{N/1280}\rfloor+1\ge64\). Both lower endpoints can be chosen within \(11200z^{2/3}\) of \(N\). The asymptotic statement for real targets follows by bracketing neighboring integers.

Choose the largest power of two \(q\le z^{1/3}\). Let \(t=qu+(q-1)/2\) be the smallest such half-integer at least \(z+q\). Then \(z+q\le t<z+2q\). Move \(a=t+d, b=t-d\) with \(d\) a positive half-integer, so their sum is odd. Write \(a=q(u+j)+v\), \(b=q(u-j)+q-1-v\), \(0\le v<q\). Unrolling the recurrence on a width-\(\log_2 q\) suffix gives
\[
p(qw+v)=p(w)+(q-1-\operatorname{rev}(v))w+p(v)+[v<q/2][w\text{ odd}],
\]
where the reversal uses the full padded suffix width. The pair of low suffixes is complementary. Since \(0\le p(w)\le w\) and \(|j|q\le d+q/2\), this gives
\[
\left|L(a,b)-G(d)\right|\le E(d)=29(2u+2q+d+1),
\quad
G(d)=1280(t^2-d^2)+640t+66d+140+29u(q-1). \tag{3}
\]
Specifically, the high-word linear terms sum to \(u(q-1)+j(q-1-2\operatorname{rev}(v))\), the two low reads are at most \(q-1\) together, the two high reads at most \(2u\), and one parity bit remains. This uses the prior complementary-carry mechanism on the two-sided *edge* formula.

Let \(D\) be the first half-integer at or above \(4\sqrt{zq}\). We have \(q\le z/16\), \(q\ge2\), \(D<4\sqrt{zq}+1\le z+1\), and \(D\le t-1\). At \(d=1/2\), \(G(d)-N\ge1280(2zq-1/4)\), while \(E(d)<38z\). At \(d=D\), \(t^2-(z-1)^2\le4zq+4q^2+2z\le6zq\). Therefore \(1280(t^2-D^2)-N\le-12800zq\), and all remaining terms of \(G(D)+E(D)\) together are less than \(1000z\) (use \(t<9z/8\), \(u\le t/q\), \(D<z+1\)). Thus the exact lower count crosses the specified target strictly.

Put \(W=z^{2/3}\). In this half-slice \(u\le2W+2\), \(q\le W\), \(d\le4W+1\) and \(W\ge16\), whence \(E(d)\le310W\). Consecutive grid values of \(G\) differ by at most \(1280(2D+1)+66\le10485W\). Consecutive *actual* lower counts consequently differ by at most \(11200W\). Bisect the integer index of the grid with opposite signs at its ends, once against \(N-1/2\) and once against \(N+1/2\). The strict low and high witnesses at the two crossings have the asserted error. Bisection requires only a sign bracket, so the fluctuations in the actual sequence cause no monotonicity assumption. Packing intervals of this length in \([X/2,X]\) yields the \(\Omega(X^{2/3})\) count of distinct rungs.

For the previously fixed scale \(N=10^{24}\), the decoder gives lower counts \(N-32,746,460,615\) and \(N+2,951,942,547\). Their respective \((a,b)\) are \((27964641611,27937064628)\) and \((27964641610,27937064629)\). Each has its literal \(+2\) mate. This bound leaves an unbounded residual at a prescribed absolute integer.

### The cube-root exponent is sharp for the chosen configurations

**Theorem.** Let \(\mathcal L\) contain every \(L(t+d,t-d)\) on the full grids in the two-sided proof: for every integer \(z\ge64\) take its specified \(q,t\) and all half-integers \(1/2\le d\le D\). Then
\[
            |\mathcal L\cap[X,2X]|=\Theta(X^{2/3}).
\]
The longest subinterval of \([X,2X]\) containing no member of \(\mathcal L\) has length \(\Theta(X^{1/3})\). Both statements also hold after adjoining every upper endpoint \(L+2\). The full independent two-flank family (2) permits other choices of \(a,b\).

**Proof.** Uniformly on a grid, \(q\asymp z^{1/3}\), \(t=z+O(z^{1/3})\), \(d=O(z^{2/3})\). Since \(0\le p(v)\le v\), equation (2) gives
\[
                L=1280(t^2-d^2)+O(t)=1280z^2+O(z^{4/3}).
\]
Thus any grid contributing a point in \([X,2X]\) has \(z\asymp\sqrt X\). There are only boundedly many possible dyadic \(q\asymp X^{1/6}\) in this range. For fixed \(q\), the allowed centers \(t\) lie in one progression of spacing \(q\) over a range of length \(O(\sqrt X)\). This gives \(O(X^{1/3})\) centers. Each permits \(O(X^{1/3})\) values of \(d\). The product bounds the number of address pairs before any coincidences of counts. Adjoining each upper endpoint at most doubles this bound, applying the same estimate in \([X-2,2X]\).

Conversely, sample targets separated by more than twice the bracketing error throughout \([5X/4,7X/4]\). Their lower witnesses are distinct and lie in \([X,2X]\), giving the matching lower bound. Splitting \([X,2X]\) at all observed endpoints leaves at least one empty interval of length \(\Omega(X^{1/3})\) by the upper bound on their number. A longer interval than a sufficiently large constant times \(X^{1/3}\) would contain a witness by bracketing its midpoint. This proves the claimed scale in both directions.

The choice of suffix width explains the exponent. For variable dyadic \(q\) the dominant terms in the step estimate have orders \(z/q\) from the unresolved high bits and \(\sqrt{zq}\) from the quadratic change in height between successive displacements. Equating them gives \(q\asymp z^{1/3}\). The counting argument shows that even perfect selection among *all* configurations on these grids cannot uniformly beat the resulting error. Supplying further addresses at the same height is mathematically necessary for that improvement. A separate complementary-digit chart also bounds its near-balanced image by counting admissible centers and suffixes before projecting to counts. [Complementary-digit inverse](../chanpak_affine_sieve_metagraph_20261122/DIGIT-INVERSE.md) · [Independent branching switches](../chanpak_affine_sieve_metagraph_20261122/UNIT-SWITCHES.md)

**A necessary imbalance even off the grids.** For the whole chart (2), let \(Z\ge1\) be an integer. The number of its lower counts in \([X,2X]\) obtainable with \(|a-b|\le Z\) is at most
\[
              2\lfloor\sqrt{2X/1280}\rfloor\lceil Z/2\rceil.
\]
Indeed \(L\ge1280ab\), so the smaller address has at most the indicated square-root number of possibilities. The positive difference is odd, leaving at most \(\lceil Z/2\rceil\) choices and two orders. The same bound applies to upper endpoints separately, since \(L+2\le2X\) still implies \(1280ab\le2X\). Thus even with *every* binary prefix and displacement, addresses with \(Z=o(\sqrt X)\) supply density zero in a dyadic target interval. At \(Z=\lfloor\sqrt X\rfloor\), the proportion of odd integers in that interval obtainable from either endpoint is at most \(4/\sqrt{640}+o(1)<0.159\). To cover all prescribed counts in this chart one would have to use substantially unequal flanks. The successful two-ended address \((17,79709)\) at the perfect-number band illustrates this choice in a different marked core. [Exact band certificate](PERFECT-BAND.md)

## One-arm target inversion and two exact landings

Set \(x=0\) and take even \(w=y>0\). The complete read (1) reduces to \(k+c w+h p(w)\), where \(c=e+o\). Two cores with the same \(h,c\) and apex counts \(k,k+2\) therefore carry a physical unit flip on every even suffix, if their rows differ on one arc. This asymmetric profile condition finds paired witnesses beyond the two-sided six-vertex chart.

For a specified lower target, split \(w=Bz+r\), \(0\le r<B\) and \(B\) a power of two. For each even \(r\) put \(A=cB+h(B-1-\operatorname{rev}(r))\) and \(T=N-k-cr-hp(r)\). The exact inverse is
\[
T=Az+h\bigl(p(z)+[r<B/2][z\text{ odd}]\bigr).
\]
Discard negative \(T\), then test the integers between \(\max(0,\lceil(T-h)/(A+h)\rceil)\) and \(\lfloor T/A\rfloor\). Taking \(B\) on the scale \(\sqrt{N/c}\) gives \(O_{c,h}(\sqrt N)\) elementary candidate evaluations and finds every even address in the given chart. It is the binary split inverse specialized to a persistent unit edge.

| Prescribed lower target | Lower core rows | Flip | Tags | \((k,h,c)\) | \(w,\ p(w)\) |
|---:|---|---|---:|---|---|
| \(4,294,967,387\) | (4,5,48,7,11,27) | 2–3 | 1 | (37,23,24) | 113404326, 68402762 |
| \(1,000,000,000,039\) | (18,24,67,69,44,15,51) | 3–5 | 126 | (419,131,1311) | 712585682, 502291378 |

For the first row the flip changes the tagged even/odd reads from \((7,17)\) to \((8,16)\) and \(k\) from 37 to 39. The second row keeps \((605,706)\) fixed while sending 419 to 421. Thus substitution into the integer equation proves exact equality at the two targets. Appending the corresponding reversed complemented word constructs the full rows. Each marked core plus apex is strong and each nonempty right word ends at a vertex beating the apex, so both larger hosts are strong. Trial division certifies the prime claims and the indicated minimality above the two round cutoffs. The upper count in the trillion example is divisible by three.

## Positive-density dyadic slopes with an explicit intercept

In the odd chamber of the two-sided core pair put \(y=2^m-1\). For every odd \(x\ge1\) and \(m\ge1\), the lower count is
\[
(320x+129)2^m+11-158x+29p(x).
\]
Consequently a prescribed slope \(q\ge449\) with \(q\equiv449\pmod{640}\) selects \(x=(q-129)/320\) immediately and gives an unbounded ray of one-arc rungs. The intercept is at most \(11-129x\). The [defect-4096 star selector](../chanpak_moonshine_prime_blocks_20261118/STAR_INVERSE_4096.md) fixes intercept one and seeks its slope. Here every slope in the indicated progression is available, with a different specified intercept. Exact coverage of arbitrary target intercepts remains open.

## The remaining absolute correction

Suppose the lower pair supplied by the bracketing theorem at each sufficiently large odd target \(N\) can be continued at *the same absolute heights* into a family with lower counts
\[
                L_-+2\sum_{j=1}^r w_j b_j,\qquad b_j\in\{0,1\}.
\]
Here each bit choice must have a physical witness. A persistent edge to its upper mate is only needed if the conclusion sought is the entire single-arc ladder. For positive sorted integer weights, all half-offsets through \(S=\sum w_j\) occur exactly when \(w_1=1\) and \(w_j\le1+\sum_{i<j}w_i\) for every later index. Inductively the next translated interval then touches or overlaps the previous one. A first failure leaves its first skipped integer beyond the reach of every later weight. Thus \(S\ge\lfloor5600 z^{2/3}\rfloor\) suffices here, since \((N-L_-)/2\) is an integer. Supplying this for all sufficiently large targets proves eventual coverage; Chan–Pak would further require the remaining finite targets through the new threshold.

This specifies the size of a sufficient correction at the improved scale. A bank covering that whole interval with binary choices requires at least \(\log_2(1+\lfloor5600z^{2/3}\rfloor)=\tfrac13\log_2 N+O(1)\) bits. If each bit is one independently additive arc reversal, the arcs must form a matching. Two arcs sharing a vertex have a nonzero mixed difference in the one-apex cut form, whichever orientations label their starting choices. Such a bank therefore needs at least twice as many vertices as bits. A fixed six-vertex core can steer the two free flank sizes while retaining a unit edge. To use a bank of this particular form for every residual, the pool of switchable vertices must grow. Other target inverses can change the flank addresses. [THM-5272, one-head cut form](../../../01-canon/theorems/THM-5272-cut-form-and-bank-rigidity-for-one-head-over-a-base-the-commutator-carry-of-the-dyadic-and-triadic-steps-and-27-63-1473-1729-in-the-descent-tree.md) · [Existing carry-bank experiments](../chanpak_carry_banks_ports_20261115/README.md)

The signed star intercept itself can be read by a single five-vertex arc; a Farey determinant selects the opposite arm for a prescribed pair of intercepts. This differs from appending a uniformly tagged block, which forces multiples of six. Nine consecutive even centers starting at the perfect number with exponent 17 have literal straddling edges, as does the perfect number with exponent 19, and a centered symmetric attempt at the full upper neighbor forces a unique binary prefix of \(1/\sqrt2\). At order nine four marked arcs suffice for every even center phase modulo 24. [Star-arm formulas and exponent 31](STAR-ARMS.md) · [Perfect-number target witnesses](PERFECT-TARGETS.md) · [Mersenne and Wieferich scales](../chanpak_perfect_wieferich_seams_20261122/PROOF.md) · [Centered inverse](../chanpak_perfect_wieferich_seams_20261122/CENTERED-INVERSE.md) · [Marked switching and moonshine](../chanpak_mod24_moonshine_ports_20261122/README.md)

## Reproduction

For a general positive odd input, `python3 replay.py --exact 8589869057` searches either endpoint and prints both outgoing-mask lists on success. It returns stored certificates directly and tries both suffix parities in the persistent fixed-left charts supplied by the target and band witnesses. Failure excludes this search only. [Perfect-number band and templates](PERFECT-BAND.md) The integer inverse can sieve its high part before testing the recurrence: with \(g=\gcd(c,h)\), a parity-adjusted intercept \(k\) requires \(g\mid N-k\), and then \((c/g)Bz\equiv(N-k)/g-(c/g)r\pmod{h/g}\). Since \(h\) is odd this gives one residue of \(z\). At \(B\) on the scale of \(\max(1,\sqrt{(N-k)g}/c)\), the previously proved interval has boundedly many candidates in this residue for the fixed core. Every survivor is checked by the full integer equality.


Run `python3 replay.py` here. The replay verifies the complete finite core reads and 50 direct subset-DP controls at \(x,y\in\{0,1,2,3,4\}\) for both main cores, recovers the two prime-selected edges and the three perfect-number edges from the inverse, checks their strongness and literal arc difference, and exercises the two-sided decoder at six target heights. Each theorem over all widths follows from the displayed recurrences. Use `python3 replay.py --target 1000000000000000000000000` for a bracket and `python3 replay.py --witness 4294967387` for both full outgoing-mask lists.
