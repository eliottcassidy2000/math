# Arithmetic guards, squarefree growth, and the prime 101 reset obstruction

**Date:** 2026-09-21. **Status:** PROVED elementary and sieve statements;
FINITE-EXACT modular certificates. Positive Collatz convergence remains
open. This note audits the attachment headed “The Affine Arithmetic Guard
Framework and Global Collatz Descent.”

The map is `T(n)=(3n+1)/2^{v_2(3n+1)}` on positive odd integers. An exponent
word specifies the exact valuations, not merely divisions that happen to
produce rational values. The source residue modulo nine will be five.

## Inheritance and concept board

The closest mechanisms are the [finite-word completion theorem](arithmetic_braids2_20260917_inverse_completion.md),
the [squarefree linear-form sieve](arithmetic_braids2_20260917_squarefree_symmetry.md),
and the [short-interval energy obstruction](collatz_blueprint_20260921_energy.md),
§5a. The hostile is a long all-one growth word; the corrected near miss
is conjoining separate existence claims without a common sieve. The
least-used sidecar here is a prime square dividing a repeated affine
orbit, retained along with the **final** endpoint.

| Lane | Exact object | Preserved coordinate / decisive test |
|---|---|---|
| Anchor | Exponent word `(1,1,5)` and its repetitions | Ordered carry 19; source cylinder; final endpoint |
| Niche | Squarefree nodes in residue five modulo nine | Local roots modulo every prime square, then a controlled sieve tail |
| Wildcard | Prime sums and square indices | Missing odd-composite count; no asserted map to a Collatz path |
| Boundary | Real energy near a dyadic minimum | Ordinary height and a short interval, in addition to congruences |

The outcome is an exact distinction: arbitrary growth length remains
compatible with squarefreeness, whereas sufficiently many repetitions
of this particular contracting reset word force a square divisor.

## 1. The actual reset family and its source guard

The attachment's path is correct:

```text
23 --1--> 35 --1--> 53 --5--> 5.
```

The labels are halving exponents. Its three-step map is

```text
R(n)=(27n+19)/128.                                         (G1)
```

The exact source guard for this word is `n=23 mod256`. Indeed, every
`n=23+256t` gives the odd nodes

```text
23+256t, 35+384t, 53+576t, 5+54t,
```

and `3*(53+576t)+1=32*(5+54t)` has exact valuation five.
Conversely, the finite-word congruence gives this unique class modulo
`2^(1+1+5+1)=256`.

Intersecting with `n=5 mod9` forces `t=0 mod9`, producing

```text
23+2304s, 35+3456s, 53+5184s, 5+486s.                      (G2)
```

These inputs have relative density `1/128` among positive odd integers
congruent to five modulo nine. Every one satisfies `R(n)<n`, because
`128n-(27n+19)=101n-19>0` for positive integer n. The endpoint equals
five only when s=0. The row condition alone does not enforce this word.

The four forms in (G2) can all be squarefree: there is no local square
obstruction. At two and three they are units; at any `p>=5` their four
forbidden classes occupy at most four of `p^2` classes. The linear-form
sieve therefore gives a positive density of all-squarefree realizations.

## 2. Residue five modulo nine still allows unlimited squarefree growth

**PROVED.** For every fixed integer `L>=1`, there are infinitely many
positive odd `n=5 mod9` for which the first L exponents are all one and
all L+1 displayed nodes are squarefree. Such inputs occur with positive
relative density in the corresponding source progression.

Put `M=2^(L+1)` and write

```text
n_j=3^j*2^(L+1-j)*q-1,       0<=j<=L.
```

The required source row is precisely

```text
q = 6*M^(-1) mod9.                                         (G3)
```

It forces `3|q` and thus every node is `-1 mod3`. At two there is no
square obstruction either. For `p>=5`, the number of forbidden residues
of q modulo `p^2` is

```text
nu_L(p)=min(L+1,ord_(p^2)(3/2)).
```

Imposing (G3) does not change this count, because nine is invertible
modulo `p^2`. In fact q=0 modulo `p^2` avoids every forbidden root, so
there is no local obstruction for any L. The relative density within
(G3) is the positive product

```text
P_L=product_(p>=5)(1-nu_L(p)/p^2)>0.                        (G4)
```

The inherited proof counts a finite set of prime conditions by CRT, then
bounds the omitted-square tail by `(L+1)X/z+O_L(sqrt(X))`. Every L here
is fixed before the limit. No lower bound uniform in L is asserted.

The earlier energy obstruction also holds in this row **simultaneously**
with squarefreeness. Restrict q to (G3) and the interval determined by

```text
2^m <= Mq-1 <= 2^m*(1+1/m).
```

The interval has length `H_m=2^m/(Mm)`. The fixed-row short-interval
sieve gives `(P_L/9+o(1))*H_m` surviving q: its error is controlled by
`O_(L,z)(1)+(L+1)H_m/z+O_L(sqrt(2^m))`. Divide by H_m, take m to
infinity with z fixed, then take z to infinity.

For the earlier proposed envelope `V(n)=ln(n)*A({log_2 n})`, the initial
excess above ln(n) is `O(m^(-3))` in this interval, while every one of
the first L height ratios exceeds `(3/2)^j`. Thus for all sufficiently
large m the surviving sources have

```text
n=5 mod9, all n_0,...,n_L squarefree,
v_2(3n_j+1)=1 for j<L,
V(T^j(n))>V(n) for every 1<=j<=L.                           (G5)
```

This directly tests the proposed “five modulo nine forces a reset”
mechanism. It establishes arbitrarily long finite obstructions, not an
infinite divergent orbit or an all-squarefree infinite orbit.

## 3. A general affine prime-square resonance lemma

**PROVED.** Let p be an odd prime and `F(x)=alpha*x+beta` an affine map
over the p-adic integers, with `alpha=1 modp` and beta a unit. Then for
every x and every positive integer t,

```text
v_p(F^t(x)-x)=v_p(t).                                      (G6)
```

Indeed, `F(x)-x=(alpha-1)x+beta` is a unit, and

```text
F^t(x)-x=((alpha-1)x+beta)*(1+alpha+...+alpha^(t-1)).
```

For alpha=1 the sum is t. Otherwise the elementary odd-prime lifting
identity gives `v_p(alpha^t-1)=v_p(alpha-1)+v_p(t)`, so the geometric
sum again has valuation `v_p(t)`. The lifting identity follows by the
binomial theorem, first for a p-fold power and then by separating the
power of p in t.

Consequently F is a **single cycle modulo `p^a`** for every `a>=1`:
the least positive return time at any x is `p^a`, the entire size of
that ring. In particular p^2 consecutive orbit values cannot all be
squarefree integers, since one is zero modulo p^2.

For a rational affine branch `F(x)=(a*x+b)/d`, a sufficient and explicit
criterion is `p|d-a`, `p` dividing neither d nor b. This is a connection
from slope gap and carry to square divisibility. It concerns repetitions
of that guarded branch; it does not assert that every Collatz trajectory
must use the branch.

The same coordinate already appears in integer closure: a fixed point
must satisfy `(d-a)*x=b`, so every prime dividing the gap d-a must divide
the carry b. A prime with `p|d-a` and `p` not dividing b obstructs that
integrality. The resonance lemma strengthens this obstruction to a
finite-squarefree-run restriction. This is one mechanism with two
consequences, not an independent classification of Collatz cycles.
Conversely, absence of such a prime does not itself prove integrality;
prime-power multiplicities still matter. In particular, the words of
any actual integer cycle, including the known `3n-1` cycles, necessarily
fail this unit-carry resonance test at every gap prime.

## 4. The sharp 101-square obstruction for repeated resets

For (G1), `128-27=101` and 19 is a unit modulo 101. Apply (G6): R is
one cycle of length `N=101^2=10201` modulo N. Index this cycle by
`x_i=R^i(0)`, so `R(x_i)=x_(i+1)` with subscripts modulo N.

Within a reset block, the three pre-final phases are

```text
P_0(x)=x,        P_1(x)=(3x+1)/2,        P_2(x)=(9x+5)/4.
```

Their zero roots modulo N, and their indices along the R-cycle, are

| Phase | Root modulo 10201 | Index i with x_i equal to the root |
|---|---:|---:|
| P_0 | 0 | 0 |
| P_1 | 3400 = -1/3 | 3795 |
| P_2 | 9067 = -5/9 | 5214 |

These are exact finite modular identities, certified independently by
iterating the cycle and by binary affine exponentiation. Together with
(G6), they determine every repeated-word length, not just the lengths
sampled in a trajectory census.

For r reset blocks, include **all `3r+1` nodes**, including the final
endpoint `R^r(n)`. The excluded source indices are

```text
{-t:0<=t<=r}
 union {3795-t:0<=t<r}
 union {5214-t:0<=t<r},                      modulo 10201.   (G7)
```

The three cyclic gaps between the anchors are 3795, 1419 and 4987.
The backward interval ending at zero is one longer because of the final
endpoint. Hence the exact number of admissible source residues is

```text
h(r)=(3795-r)_+ +(1419-r)_+ +(4986-r)_+,                    (G8)
where (t)_+=max(t,0).
```

For `r=4985`, exactly one source residue survives: `n=7332 mod10201`.
For `r=4986`, none survives. Omitting the final endpoint would give the
wrong cutoff, leaving one residue at r=4986.

Therefore any positive integer realization of 4986 consecutive copies
of `(1,1,5)` has a displayed node divisible by `101^2`. This does not
make the word itself impossible. Every finite halving word still has
positive integer realizations; those realizations simply cannot all
remain squarefree at this length.

## 5. All other primes permit the repeated word

The obstruction at 101 is the **only** local obstruction to simultaneous
squarefreeness, for any number of repeated blocks in the specified row.

For `p>=5`, `p!=101`, use the affine fixed point

```text
f=19/101,       R(f)=f.
```

The three phases at f are respectively

```text
19/101,        79/101,        169/101.                      (G9)
```

Each is nonzero modulo p^2 except that the third vanishes for p=13.
Thus f is an avoiding residue modulo p^2 for every such prime except 13,
uniformly for every number of blocks. The primes 19 and 79 cause no
problem: those numerators have valuation one, not two.

At p=13 take instead `n=f+13 mod169` (the residue 50). With
`alpha=27/128`,

```text
R^t(n)=f+13*alpha^t modulo169.
```

Since alpha is a 13-adic unit, P_0 and P_1 remain units modulo 13,
and P_2 equals `(169/101)+(9/4)*13*alpha^t`, which has valuation
**exactly one**. This avoids all square zeros forever. The modular
control has period twelve. At two and three, the original word guard
and the row `5 mod9` make every displayed integer odd and prime to three.

## 6. Exact squarefree-existence threshold, with positive density

**PROVED.** For an integer `r>=1`, the following are equivalent:

1. There exists a positive odd integer `n=5 mod9` realizing r successive
   copies of `(1,1,5)` with all `3r+1` nodes squarefree.
2. Such sources have positive relative natural density inside the exact
   word-and-row progression.
3. `r<=4985`.

The implication 1 to 3 is (G8). To prove 3 to 2, the word specifies one
odd class modulo `2^(7r+1)`. CRT intersects it with `5 mod9` to give one
class `n=a_r+Q_r*s`, where `Q_r=9*2^(7r+1)`. Its jth iterate is the
integer linear form

```text
n_j(s)=c_j+A_j*s,
A_j=9*3^j*2^(7r+1-K_j),      0<=j<=3r.                    (G10)
```

Here `K_j` is the cumulative exponent along the repeated word; c_j is
the jth iterate of a positive representative a_r. For every `p>=5`,
the source coefficient Q_r is a unit modulo p^2, so s runs over all
source residues in that ring. Sections 4–5 supply at least one common
avoiding residue for every prime. If `nu_r(p)` counts the excluded
s-classes, then `nu_r(p)<p^2` and `nu_r(p)<=3r+1` for `p>=5`.
The positive density is

```text
delta_r=product_(p>=5)(1-nu_r(p)/p^2)>0,                  (G11)
with its 101-factor exactly h(r)/10201.
```

For completeness, finite-prime CRT gives the truncated product. Up to
parameter height X, an omitted square divisor has `p=O_r(sqrt(X))`;
the number excluded is at most `(3r+1)X/z+O_r(sqrt(X))`. Taking X
then z to infinity proves the product formula. The sum
`sum_p nu_r(p)/p^2` converges for fixed r, so positive local factors
give a positive product. This justifies existence, rather than assuming
that compatibility at every finite collection of primes is enough.
Finally 2 implies 1 immediately.

Relative to **all** positive odd sources in `5 mod9`, the density in 2
is `2^(-7r)*delta_r`. No usable lower bound uniform in r is claimed;
even the unfiltered word cylinder becomes extremely thin. These sources
perform r contracting blocks, but their eventual convergence is a separate
question. Their terminal values have not been fixed to a known cycle.

The genuine transfer is now explicit: a guarded word with slope/carry
maps to an affine permutation modulo a resonant prime square; positions
of intermediate zero roots determine the longest permitted squarefree
run. Passing to just the real contraction factor loses this information.

### 6a. Constant-exponent words have an elementary sharp threshold

**PROVED.** Fix `k>=3`, and let p be the least prime factor of `2^k-3`.
The word consisting of r copies of the exponent k can have all r+1
displayed nodes squarefree, with positive relative density of sources
also satisfying `n=5 mod9`, **if and only if**

```text
r<=p^2-2.                                                  (G13)
```

For the one-step branch `F_k(n)=(3n+1)/2^k`, every prime q dividing
`2^k-3` satisfies the resonance lemma and gives one cycle modulo q^2.
Its r+1 displayed values can avoid zero in that ring precisely when
`r+1<q^2`. The most restrictive of these conditions is that for the
least prime factor p. At any other odd prime q the fixed point
`1/(2^k-3)` is a unit modulo q^2 and avoids zero forever. The source
word/row CRT handles two and three. The same fixed-r linear-form sieve
then proves positive density whenever every local condition is possible.

In particular, **24 consecutive exponents equal to three force a square
factor of five somewhere among the 25 nodes**. Twenty-three such steps
still have all-squarefree realizations. For exponent four, the sharp
impossible length is 168, caused by `13^2`.

For `k=1,2`, the slope gaps are respectively `-1,1`, with no gap prime.
Their fixed points `-1,1` avoid zero modulo every odd prime square; all
finite repetition lengths admit all-squarefree positive realizations
in the prescribed source row. Here the fixed point -1 is only a modular
avoiding seed; the actual finite realizations remain positive integers.

Thus squarefree **starting points** forbid no finite halving word,
whereas squarefreeness at **every intermediate node** imposes genuine
forbidden finite words. These restrictions still allow arbitrarily long
growing words and do not imply global descent.

### 6b. A finite local criterion for every word; the shortest obstruction

**PROVED.** Fix any word `w=(k_1,...,k_L)` of positive halving exponents,
including the empty word if `L=0`. Among its positive odd realizations
with source `n=5 mod9`, all L+1 displayed nodes can be simultaneously
squarefree if and only if the following **finite** local test passes:

```text
For every prime 5<=p<=sqrt(L+1),
the prefix zero roots do not cover all of Z/(p^2)Z.         (G14)
```

Whenever the test passes, such sources have positive relative density
inside the exact word-and-row progression. To compute the roots, put
`K_0=B_0=0`, and update

```text
K_(j+1)=K_j+k_(j+1),
B_(j+1)=3*B_j+2^K_j.
```

The jth node is `(3^j*n+B_j)/2^K_j`. For every `p>=5`, its unique zero
root in the source coordinate is

```text
rho_j=-B_j*(3^j)^(-1) modp^2,        0<=j<=L.              (G15)
```

Let `nu_w(p)=#{rho_0,...,rho_L}`. The source progression has modulus
`9*2^(K_L+1)`, a unit modulo p^2, so it samples all these source residues.
At two, the exact word guard makes every node odd. At three, the source
row is a unit and every successor stays a unit. Consequently those two
primes introduce no square obstruction. For `p>sqrt(L+1)`, the bound
`nu_w(p)<=L+1<p^2` rules out a local obstruction automatically. This
proves that (G14) is necessary and lists every possible obstructing prime.

If (G14) passes, all local factors are positive. Every displayed node
is a linear form in the source-progression parameter, with coefficients
having prime factors only two and three. The same fixed-word sieve tail
as in §6 proves density

```text
delta_w=product_(p>=5)(1-nu_w(p)/p^2)>0.                   (G16)
```

This establishes sufficiency, not just compatibility for finitely many
primes. Computing (G14) therefore decides simultaneous-squarefree
realizability of a given **finite** word. It is not a decision procedure
for an entire Collatz orbit, and does not fix a terminating target.

**Sharp length corollary.** Every word of length `L<=23` passes, because
its L+1 zero roots cannot cover a prime-square ring of size at least 25.
At `L=24`, the word with every exponent three has all 25 roots modulo
25 by §6a and fails. Thus **24 is the shortest possible length of a
halving word that forbids squarefreeness at every displayed node** in
the source row `5 mod9`. This is a statement about 24 transitions and
25 nodes; counting only the pre-final nodes would change the boundary.

## 7. The prime-sum identity at 196: true arithmetic and its actual sidecar

The numerical equality is correct:

```text
1+3+5+7+11+13+17+19+23+29+31+37 = 196 = 14^2.
```

There are eleven prime additions. The cumulative sequence, counting
the initial one as its first term, is

```text
1,4,9,16,27,40,57,76,99,128,159,196.
```

Thus 196 is its twelfth term under that convention. The extra final 7
in the attachment's parenthesized prime list is not consecutive and
would instead give 203 as a thirteenth term.

There is a useful exact comparison with square growth. Let p_i be the
ith odd prime and c_i the number of odd composites from three through
p_i. Counting the odd integers gives `p_i=2i+1+2c_i`, so

```text
1+sum_(i=1)^k p_i = (k+1)^2 + 2*sum_(i=1)^k c_i.           (G12)
```

The first three additions match successive squares because c_i=0 there.
For k>=4 the correction is positive. At k=11, the c_i are
`0,0,0,1,1,2,2,3,5,5,7`, summing to 26, so `196=12^2+52`.
The square coincidence is equivalently
`(14-12)*(14+12)=2*26`. This identifies the missing additive coordinate:
the cumulative count of skipped odd composites. It supplies no map to
Collatz exponent words, carry 19, or resonant prime 101. No causal
prime-square interpretation of the isolated number 196 is claimed.

## Reproduction and scope

```text
python 04-computation/experiments/collatz_guards_20260921_squarefree.py
python -O 04-computation/experiments/collatz_guards_20260921_squarefree.py --output squarefree-optimized.json
```

The [script](../../04-computation/experiments/collatz_guards_20260921_squarefree.py)
and [JSON](../../04-computation/experiments/collatz_guards_20260921_squarefree.json)
contain the complete 101-square orbit certificate, independent affine
powers for the three anchors, every survivor count through r=5100,
the p=13 repair, finite controls at other primes, and exact integer
squarefree trajectories for small r. Constant-exponent controls check
the full cycles modulo 25 and 169 and their endpoint-count thresholds.
The general root-union routine is checked against direct source-residue
enumeration for 23 and 24 exponents equal to three, and 24 equal to one.
Separate CRT/short-interval controls
check growing squarefree sources in `5 mod9`, with rigorous rational
energy-excess bounds. There is no enumeration of huge actual sources
for r=4985: their existence follows from the fixed-r sieve proof.

The independent summand lane audited the cycle, root indices, final-endpoint
off-by-one boundary, p=13 repair, all-prime sieve threshold in §6,
constant-exponent threshold in §6a, the finite local criterion and sharp
length-24 boundary in §6b, and the gap/carry integrality bridge. All passed.
The focused §6b review also checked the generic root routine against its
independent source-enumeration control. Before that extension, the
independent replay matched that version's ordinary/optimized JSON. The
final extended version's ordinary and optimized replays also match;
the additional independent §6b audit was proof and code review.
Parent integration handles the maintained
truth surfaces and git; this lane owns only this note and its experiment
artifacts.
