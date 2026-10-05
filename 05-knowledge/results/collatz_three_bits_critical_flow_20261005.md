# Three bits, three colors, and a critical Collatz flow

2026-10-05. **PROVED, scoped and author-audited:** the independent code
factorization, the climb/cofactor factorization, the critical-flow criterion,
the finite ascent resolvent, the refuel counterfamilies, the cycle-anchor
closures, and the finite-stage debt identities below. **FINITE-EXACT:** the
accompanying checks. **OPEN:** a
summable flow paying both kinds of edge, a uniform bound for the debt closures,
and universal Collatz. No claim of literature priority or canon promotion.

[Program](../../04-computation/experiments/collatz_three_bits_critical_flow_20261005.py)
and [exact output](collatz_three_bits_critical_flow_20261005.json).

## 1. Inheritance and the research board

The anchor is universal entry into rooted Collatz certificates. The niche is
the three-color Zeckendorf reader; the wildcard is the sequence 2,10,42 and
the information retained by redundant integer descriptions.

Closest mechanisms are [P6–P8 of the atomic-prefix note](collatz_effective_prefix_mass_20261005.md)
and the independently derived [gamma-code measure](algorithmic_collatz_measure_20261005.md).
The complementary [paid-guard route](paid_guard_budget_20261005.md) retains
its four-credit rule on 219 mod256 and its named 3-power-family census;
the present flow construction does not replace that certificate machinery.
The hostile examples are an arbitrarily long initial climb and the minus
cycles with minima 5 and 17. The corrected near miss is treating a finite
color or a complete residue tree as an ordinary-integer coverage proof.
The least-used sidecar here is **ternary ancestry depth v3(n+1)**.

The board has six objects: **source atom / description alias / binary climb
depth / ternary ancestry depth / golden carry / aggregate incoming flow**.
Incoming commit `4afc2fcaeb26` supplies the exact measure comparison. Older
[ternary Berggren work, B11–B13](ternary_berggren_20260925.md) supplies the
actual source of 0,2,10,42. The [guard automaton](zeckendorf_guard_automaton_20261003.md)
and [marked-unit construction](duck_zeckendorf_20260925.md) keep the different
three-color conventions separate. No additional remote commits appeared at
the mid-session fetch.

Use U(n)=oddpart(3n+1) on positive odd integers, ell(n)=bit_length(n), and
V={3,5,7,...}. An edge reaching 1 is killed in every flow operator below.
The measure names are standardized here; earlier documents used different
letters for the same measures.

## 2. Three appearances of exactly the same probability law

Let

\[
 \nu_n=\frac8{3\,4^{\ell(n)}},\qquad
 \mu_n=2^{-(2\lfloor\log_2((n+1)/2)\rfloor+1)},\qquad
 q(a)=\frac3{4^a}\quad(a\ge1).
\]

Both integer measures have total mass one. The incoming comparison is

\[
 \frac{\nu_n}{\mu_n}=
 \begin{cases}4/3,&(n+1)/2\text{ is a power of two},\\
               1/3,&\text{otherwise}.
 \end{cases}                                                   \tag{1}
\]

Consequently they have the same null sets and equivalent mass-one coverage
criteria. The exceptional locus in (1) is precisely the all-one family
n=2^h−1, which carries arbitrarily long initial climbs. Its total mass is
8/9 under nu and 2/3 under mu. Thus the comparison's special factor already
points to the canonical hostile family; its bounded factor cannot pay an
unbounded climb.

**C1 — independent description alias.** Choose a header height D>=1 with
probability 2^-D, then uniformly choose an odd N<2^D. This is the complete
prefix code of the inherited note, with length 2D−1. For a fixed n and
j=D−ell(n)>=0,

\[
 \Pr(N=n,J=j)=\frac2{4^{\ell(n)+j}}
             =\nu_n\frac3{4^{j+1}}.                           \tag{2}
\]

Hence N and A=J+1 are independent, and A has law q. The extra coordinate is
real information about a description of n, although it changes no integer.

This is the same q obtained by normalizing the square of the Haar valuation
law p(a)=2^-a:

\[
 \sum_a p(a)^2=1/3,\qquad q(a)=p(a)^2/(1/3),\qquad
 \frac{q(a)}{p(a)}=\frac3{2^a}.                               \tag{3}
\]

Thus q is also the common value of two independent p-valuations conditioned
on their equality. For a word of length r and total valuation A, (3) gives
3^r/2^A, the affine coefficient in the existing likelihood calculation.
It is a change of law on words, not a distributional assertion about the
successive valuations of one fixed integer.

**C2 — independent climb depth and cofactor.** Write uniquely

\[
 N=2^H T-1,\qquad H=v_2(N+1)\ge1,\quad T\text{ odd}.
\]

Under nu, H and T are independent, with

\[
 \Pr(H=h)=q(h),\qquad
 \Pr(T=t)=\begin{cases}8/9,&t=1,\\ \nu_t/3,&t>1.\end{cases}    \tag{4}
\]

Proof: if t=1 then ell(2^h t−1)=h; if t>1 is odd then its bit length is
h+ell(t). Substitution into nu gives (4), including its normalization.
The cofactor law can also be written (2/3)delta_1+(1/3)nu.

Equations (2) and (4) imply that the description alias A and the actual
climb depth H are independent copies of q. They do **not** identify A
with H. One is free coding redundancy; the other is fixed by the integer.

There is a literal 3+1 sampling model for q: at each level choose one of
four equally likely outcomes, continue on one outcome, and stop on the
other three. The first stopping level has probability 3/4^a. Remembering
which of the three stopping outcomes occurred adds a separate ternary
label. Dropping that label preserves the stopping level, not the whole word.

### What “three bits” measures

The gamma code is one-to-one on source indices, so H(mu)=3 bits. The
redundant code in (2) also has entropy and expected length 3, but

\[
 H(\nu)=\log_2 3+\tfrac13,\quad
 H(q)=\tfrac83-\log_2 3,\quad H(N,A)=H(\nu)+H(q)=3.             \tag{5}
\]

Indeed E_nu ell=5/3 and E_q A=4/3; inserting the log probabilities proves
(5). Three expected bits describe an infinite-support distribution. They
are not three fixed bits holding every integer. This distinguishes an
entropy identity from an eight-state representation.

## 3. The 2,10,42 ray is both certified and ternary-complete

For j>=1 set

\[
 a_j=\frac{2(4^j-1)}3,\quad b_j=a_j/2,\quad m_j=(b_j+1)/2.
\]

| j | a_j | b_j | m_j |
|---:|---:|---:|---:|
| 1 | 2 | 1 | 1 |
| 2 | 10 | 5 | 3 |
| 3 | 42 | 21 | 11 |
| 4 | 170 | 85 | 43 |

The binary expansion of a_j is `10` repeated j times. Exact identities are

\[
 a_{j+1}=4a_j+2,\quad 3b_j+1=4^j,\quad 3m_j-1=2^{2j-1}.        \tag{6}
\]

Thus b_j is a one-edge plus-root family and m_j a one-edge minus-root
family. The latter integers are the gamma source indices of the former.
The general seam is U_+(2m−1)=oddpart(3m−1); it is a one-edge identity,
not a conjugacy of the two full dynamical systems.

The plus-root ray has exact atomic masses

\[
 \nu\{b_j:j\ge1\}=32/45,\qquad
 \mu\{b_j:j\ge1\}=19/30.                                   \tag{7}
\]

For nu the summands are (32/3)16^-j. For mu the first is 1/2 and the
others are 32*16^-j. The remainder of all integers is still the problem.

Including a_0=0, the inherited identity
v3(a_j−a_i)=v3(j−i) makes this a bijective odometer modulo every 3^d.
It is simultaneously a sparse ordinary sequence and a complete ternary
residue tree. The missing sidecar when passing back to integers is the
ordinary nonnegative phase j. This is precisely why residue completeness
cannot supply universal Collatz coverage by itself.

## 4. The three-color reader supplies a typed finite coordinate

The owner-confirmed 35-symbol stripe is
`RKBRKRKBRKBRKRKBRKRKBRKBRKRKBRKBRKB`.
The marked-unit note shows why a tempting Wythoff continuation agrees only
through position 34. Its full infinite continuation is not supplied here.
The marked representation has two or three representatives of each value;
its marker is part of the data. It must not be silently replaced by the
following auxiliary XOR charge.

Put alpha=phi^-2 and b(n)=floor(alpha(n+1)). The separately defined charge
is Q(n)=(n−2b(n),b(n)), reduced modulo two when used as a color. Its high-to-low
Zeckendorf digit reader is

\[
 (X,Y)\longmapsto(Y+d,X+Y+2d),\qquad
 (X,Y)=(n,2n-b(n))\text{ at the end}.                         \tag{8}
\]

The previous digit rejects `11`. Modulo two, (8) uses two charge bits and
one guard bit: exactly eight reachable states. The linear part
M(x,y)=(y,x+y) fixes zero and cycles the other three charges, with M^3=I.
This is a precise 1+3 structure, or the action of a nontrivial element of
F_4^* after identifying the two-dimensional vector space with F_4.

For actual multiplication by three and addition of one, the missing carry is

\[
 \eta(n)=b(3n+1)-3b(n)
        =\lfloor3\{\alpha(n+1)\}-\alpha\rfloor\in\{-1,0,1,2\},
\]
\[
 Q(3n+1)=3Q(n)+(1-2\eta,\eta).                              \tag{9}
\]

All four cases occur on odd inputs: n=7,5,1,9 give eta=−1,0,1,2.
Two carry bits distinguish these cases; the charge modulo two retains
only eta modulo two. Computing the odd endpoint still requires the exact
dyadic valuation of 3n+1.

The [three-bit Fano calculation, F3–F4](creation_fano_20260925.md) provides a
useful adverse comparison. Its plus parity map modulo eight is
[0,5,2,3,4,1,6,7], with a nonlinear carry term x0*x1. The linear golden
reader and the nonlinear Collatz address map therefore cannot be identified
just because each uses three bits. Their exact integer and carry coordinates
give a legitimate combined observer; no ranking inequality follows merely
from the number of its labels.

## 5. Reroute: a uniform discount is unnecessary

For a nonnegative mass vector v on V define

\[
 (K v)(m)=\sum_{n\in V:U(n)=m}v(n).
\]

**C3 — critical-flow criterion.** Universal positive Collatz is equivalent
to the existence of a strictly positive, summable v on V satisfying

\[
                         Kv\le v.                           \tag{10}
\]

The older sufficient criterion required Kv<=rho*v for one rho<1. Equation
(10) needs neither a uniform discount nor a strict inequality at each vertex.

Proof of sufficiency: along any non-root orbit,
v(U(n))>=v(n). An infinite orbit with distinct vertices would give infinitely
many terms at least v(n)>0 in a summable series. Thus a non-root orbit must
reach a finite cycle C. Every vertex m in such a cycle is nonzero modulo
three. Its odd predecessors (2^a m−1)/3 exist for infinitely many positive
a of one parity. Choose one outside C. Summing (10) over C then gives

\[
 \sum_{m\in C}v(m)\ge\sum_{m\in C}(Kv)(m)
                    \ge\sum_{m\in C}v(m)+v(p),
\]

which contradicts v(p)>0. This use of an **external predecessor** is the
reason equality in the local payment rule is allowed. Without that property,
an isolated fixed point with positive stationary mass would be a counterexample
to the general graph statement.

Conversely, assume all sources root. Enumerate n_j=2j+1, j>=1, let tau_j
be their odd-step root times, and inject a_j=2^-j/tau_j along each of the
tau_j non-root states of its path. The resulting v has total mass one,
v(n_j)>=a_j>0, and Kv=v−a. Its tail over source indices j>J is at most
2^-J. This is an existence proof conditional on convergence, not an
independent construction of the desired flow.

The minus cycles {5,7} and {17,25,37,55,41,61,91} fail (10) when only
the root 1 is killed: the external incoming edges 27->5 and 363->17
already force the cycle-sum contradictions. Killing certified cycle sets
instead of just 1 gives the corresponding criterion for entry into them,
provided every hypothetical remaining cycle has an external predecessor.

### Why the third payment bit can be dropped

On a rising edge h=v2(n+1)>1, h decreases by one and ell increases by
at most one. Thus nu(n)*4^-h is nonincreasing along that edge in the
incoming-flow direction. The earlier 8^-h factor paid two bits for the
possible height increase and one extra bit for a uniform factor 1/2.
C3 removes the need for that extra discount bit. It does not remove the
need to pay the next refuel, or to sum all incoming edges.

## 6. An exact finite-mass payment for every climb

Split K=R+D: R contains n=3 mod4 rising edges, and D contains falling
edges from n=1 mod4, with the root still killed. Define B=I+R+R^2+... .
For a point m its backward R-chain has exactly s=v3(m+1) edges:

\[
 p_j=\frac{2^j(m+1)}{3^j}-1,\qquad 0\le j\le s.             \tag{11}
\]

Every p_j is a positive odd integer; for m>1 none is the killed root.
The other direction has h(n)=v2(n+1) vertices before its first fall.

**C4 — complete ascent closure.** With a=nu restricted to V,

\[
 W=B a,\qquad
 W(m)=\sum_{j=0}^{v_3(m+1)}\nu\left(\frac{2^j(m+1)}{3^j}-1\right),
\]
\[
 W-RW=a,\qquad \sum_{m\in V}W(m)
 =\sum_{n\in V}h(n)\nu_n=\frac43-\frac23=\frac23.             \tag{12}
\]

The last equality uses C2 and subtracts the root's mass. All sums are
nonnegative, so exchanging the order is justified. This construction is
unconditional and pays **all** rising edges, not just a finite census.

It achieves an unbounded correction to nu with finite total mass. Along
m=2*3^t−1, the term p_t=2^(t+1)−1 alone gives
W(m)/nu(m)>=4^(ell(m)−t−1), which is unbounded. This is exactly the kind
of correction forbidden to finite color reweightings but allowed by the
unbounded ternary ancestry coordinate.

### The first failure and an unbounded refuel family

The first falling edge failing W(n)<=W(U(n)), in increasing source order,
is 17->13. Here the ascent chain is 7->11->17 and

\[
 W(17)=\nu_7+\nu_{11}+\nu_{17}=7/128,
 \quad W(13)=1/96,\quad W(17)/W(13)=21/4.                    \tag{13}
\]

The program checks all preceding falling sources directly. The family
2*3^(h−1)−1 -> (3^h−1)/2 for odd h>=3 yields unbounded failures of the
same form: the target has v3(m+1)=0, while the source's ancestry contains
2^h−1, giving ratio at least 4^(ell(m)−h), which diverges. This is an
obstruction to this candidate W, not to critical flows.

A different family pinpoints the refuel's ternary cost. Put

\[
 h=2\,3^s-1,\quad n_h=(2^{h+2}-5)/3,\quad m_h=2^h-1,
 \qquad s\ge1.
\]

Then U(n_h)=m_h, v3(n_h+1)=s, and W(m_h)=nu(m_h). In (11) the deepest
ancestor of n_h satisfies

\[
 p_s+1=(2/3)^s(n_h+1)<\tfrac43(2/3)^s2^h.
\]

Since 2^ell(p_s)<=2(p_s+1),

\[
 \frac{W(n_h)}{W(m_h)}\ge\frac{\nu(p_s)}{\nu(m_h)}
              >\frac9{64}(9/4)^s.                           \tag{14}
\]

The ratios for s=1,...,7 are 5/4,9/4,25/4,41/4,105/4,361/4,617/4.
The first edge is 41->31. Thus the simple repair relocates the obstruction
from current binary fuel to an unbounded ternary ancestry bill. The two
depths are now exact coordinates of one calculation.

## 7. A computable debt ledger after any finite number of refuels

For n=2^h t−1 the last vertex of the initial climb is
P(n)=2*3^(h−1)t−1, and its falling endpoint is

\[
 F(n)=U(P(n))=\operatorname{oddpart}(3^h t-1).                \tag{15}
\]

This block uses exactly h odd steps. Let C=DB: it pushes each source mass
once to F(n), killing it if F(n)=1. Every non-root block has finite length,
without any Collatz convergence assumption.

**C5 — exact boundary debt.** Define

\[
 W_r=\sum_{j=0}^r B C^j a.
\]

Then

\[
                 (I-K)W_r=a-C^{r+1}a.                     \tag{16}
\]

Indeed (I−K)B=I−DB=I−C, and the sum telescopes. Equation (16) retains the
whole unpaid boundary distribution. Its label is the actual endpoint,
not just the number of successful local rules.

Each W_r is summable and computable with an explicit source-tail bound.
Since ell(F(n))<=2ell(n), its cost from source n is at most
(2^(r+1)−1)ell(n). Also

\[
 \sum_{n\ge2^L}\ell(n)\nu_n=\frac{2(L+2)}{3\,2^L}.
\]

Consequently truncating the source head at 2^L leaves at most

\[
 (2^{r+1}-1)\frac{2(L+2)}{3\,2^L}                           \tag{17}
\]

uncomputed total flow mass. This is an ordinary-integer error bound, not a
Haar null-set argument.

There are two different remaining targets:

- **Exact coverage target:** ||C^r a||_1 -> 0. This is equivalent to every
  source reaching 1, because every source has a positive atom and every
  block is finite. It requires no bound on mean stopping time.
- **Stronger sufficient flow target:** sup_r ||W_r||_1 < infinity. The
  monotone limit is then summable and satisfies KW=W−a, giving C3. This
  target asks for finite nu-weighted mean odd stopping time; universal
  termination by itself does not establish that moment bound.

For the finite head of odd n<2^14, the complete weighted odd-step cost is
8838787/6291456, with maximum 101 steps at 13255. The remaining input mass
is 1/24576. That small input mass does not bound its unknown stopping-time
cost. The output records exact partial closures through r=8 and the
rigorous error (17), without extrapolating a finite plateau to a limit.

## 8. A second reroute: resum the three inherited cycle anchors

**C6 — fixed-anchor closures.** The actual odd valuation words of the three
known minus cycles with minima 1,5,17 give these plus-sheet macros:

| Anchor | Valuation word | P | Q | B in (Pn+B)/Q |
|---:|---|---:|---:|---:|
| −1 | (1) | 3 | 2 | 1 |
| −5 | (1,2) | 9 | 8 | 5 |
| −17 | (1,1,1,2,1,1,4) | 2187 | 2048 | 2363 |

Let r be the word length and A its valuation sum, so P=3^r, Q=2^A.
For anchor −c the macro satisfies

\[
 G_c(n)+c=(P/Q)(n+c),\qquad B=c(P-Q).                         \tag{18}
\]

Its exact source guard is n=−c modulo 2^(A+1). The word is legal at −c,
and the exact-valuation cylinder has that modulus; the program also replays
every intermediate valuation on positive inputs. All three macros expand
positive inputs, so they cannot pass through the root.

The number of consecutive admissible repetitions from positive odd n is

\[
 J_c(n)=\left\lfloor\frac{v_2(n+c)-1}{A}\right\rfloor.        \tag{19}
\]

Every repetition spends A binary digits of n+c. Backward predecessors are

\[
 p_j=\frac{2^{Aj}(m+c)}{3^{rj}}-c,
 \qquad0\le j\le\left\lfloor\frac{v_3(m+c)}r\right\rfloor.   \tag{20}
\]

For these three anchors all such predecessors lie in V: for j>=1 the
even quotient (m+c)/3^(rj) is at least two, giving p_j>=2*2^(Aj)−c>1.
Their source guards follow from v2(p_j+c)>=Aj+1. Thus (20) is an exact
finite ancestry formula, not just an integrality test.

Let R_c be the macro pushforward and B_c=I+R_c+R_c^2+... . Then W_c=B_c a
satisfies W_c−R_c W_c=a. Its total mass is particularly simple. Since
c<2^A for each row, the cylinder J_c>=j has canonical representative
2^(Aj+1)−c of bit length Aj+1. The atomic cylinder formula gives

\[
 \nu\{J_c\ge j\}=4^{-Aj},\qquad
 \sum_{m\in V}W_c(m)=\frac13+\frac1{4^A-1}.                 \tag{21}
\]

The root contributes no repetitions. These are three unconditional finite
closures, of masses

\[
 2/3,\qquad 1/3+1/63,\qquad 1/3+1/4194303.
\]

This gives **63 a specific new role**: the −5 pattern has valuation cost
three, so its repeat probability under nu is 1/4^3=1/64, and the total
extra mass for all repeat depths is 1/(64−1)=1/63. The integer 11 appears
as the actual valuation sum of the −17 cycle word. These are derived
identities, not identifications with every other occurrence of 63 or 11
in the repository.

The −1 row is C4; the −5 row recovers the earlier 9/8 growth law anchored
at −5. The −17 row starts at positive source 4079 and reaches 4357 in
seven odd steps. A single ternary factor is insufficient for either
larger macro: at m=1,c=5 the alleged predecessor is 1/3, not an integer.

The broad −1 guard overlaps the other two. Here is an explicit disjoint
policy that supplies **universal legal block coverage**:

| Letter | Source guard, in priority order | Block | Exact nu mass on V |
|---|---|---|---:|
| H | n=−17 mod4096 | seven-step −17 word | 1/4194304 |
| G | n=−5 mod16 | two-step −5 word | 1/64 |
| R | remaining n=3 mod4 | one rising step | 1/4−1/64−1/4194304 |
| D | n=1 mod4, n>1 | one actual falling step | 1/12 |

The H and G guards are disjoint already (15 and 11 modulo16). Every
positive odd n>1 has exactly one legal next block. The first three blocks
increase n; D decreases it. Grouping an actual orbit this way neither
creates nor removes root entry, since block lengths lie between one and
seven and no expanding block can pass through 1. Thus this is a concrete
three-expansion/one-fall representation covering every input, with its
root-entry assertion still OPEN. Its four letters are newly defined here;
they are not identified with another four-letter algebra by label count.

There is an exact transfer back to the required flow. If the policy's
endpoint operator is P_* and Lw puts w(n) on every preterminal state of
its chosen block, then

\[
 K Lw=Lw-w+P_*w,\qquad \|Lw\|_1\le7\|w\|_1.                \tag{22}
\]

Internal edges cancel, leaving the source and endpoint. Hence a positive
summable w with P_*w<=w would supply the original critical flow v=Lw.
The retained word and its exact source guard make this a sound reduction;
the program independently verifies the boundary identity on a finite head.

Adding the separate flows does not automatically pay this policy, nor
does closing each constant-anchor run pay arbitrary switching. The R
domain restriction only reduces its internal-run resolvent, so (12)
remains an upper bound there. The unhandled target is the distribution of
**anchor switches and residual falls**, after internal repeats have been
summed exactly. No completeness of the list of negative cycles is assumed.

## 9. What the enriched integer should retain next

A useful research state is

\[
 (n;\ h=v_2(n+1),\ s=v_3(n+1),\ t=(n+1)/2^h;
       Q(n),\text{golden carry};\ \text{labelled unpaid frontier}).
\]

Projection to n is exact. The p-adic depths pay/describe ascent transport;
Q and its carry provide a small compatible observer; the frontier records
the still-unpaid mathematical obligation. The independent code alias may
schedule additional search, but it cannot manufacture arithmetic credit.

**Projection guard.** For any lifted deterministic system projecting to U,
with killing allowed only over the certified root, a summable subinvariant
lifted mass pushes forward to a mass satisfying (10). Therefore adding
labels while keeping the marginal uniformly comparable to nu cannot
evade the inherited long-climb obstruction. Finite positive reweightings
are such bounded corrections. W in (12) succeeds on climbs precisely
because its marginal correction is unbounded.

A false positive explains why the killing condition matters. Over an
unrooted identity map, attach j>=0 with mass 2^(-j−1), send j to j−1,
and kill at j=0. The incoming mass is half the stored mass everywhere,
yet the integer never reached a root. Budget exhaustion is not a terminal
certificate. Our block operator kills only on an actual U-path to 1.

Three concrete helper questions follow from the successful closure:

1. Can the ancestry corrections (12) and (20) be enlarged to a summable
   cone stable under a disjoint refuel/anchor-switch policy? Test it first
   on 17->13 and both infinite counterfamilies, not just typical inputs.
2. Can a prefix decomposition of the joint binary/ternary carry prove
   a bound on ||C^r a|| for all r, retaining each exceptional integer atom?
   This is weaker than demanding a uniformly bounded mean-time ledger.
3. Can golden carry states reduce the complexity of that exact boundary
   distribution? A useful answer must preserve the ordinary source guard
   and improve a stated bound in (16) or (17). An eight-state color count
   alone is not that answer.

These questions reroute a paid-controller search from uniform contraction
to critical flow, and from repeated checking of individual climbs to the
unpaid refuel distribution. The new exact bridge is q across coding,
climb statistics, and word likelihood; the new constructive object is the
finite-mass ascent closure. Universal coverage remains the open target.

## 10. Reproduction and scope

Run from the repository root:

```sh
python3 04-computation/experiments/collatz_three_bits_critical_flow_20261005.py
python3 -O 04-computation/experiments/collatz_three_bits_critical_flow_20261005.py
```

The program uses exact integers, fractions and integer square roots; checks
remain active under -O. It separately accumulates full ascent paths to
check (12), checks every measure comparison below 2^16, all code alias
offsets 0..6 there, all golden inputs below 20000, ternary odometers through
depth seven, forty root-ray terms with exact infinite tails, and the stated
finite root head. Negative controls include both minus cycles, the parity
carry, and the overdraft pseudo-proof. The JSON includes the source hash.
No finite check is used as the proof of an unbounded statement above.
The final run contains **899,089 exact checks**; normal and optimized runs
produce identical output.
