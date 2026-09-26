# Uniform growing-prefix multiplicative rank from incident carry primes

**Status: PROVED application of the CITED Beukers--Schlickewei bound and
the canon Terras bijection; FINITE-EXACT controls. This is a statement
about almost all starting integers, uniformly over intervals. It does
not prove G2 on every injective orbit or prove Collatz.** 2026-09-26.

Reproduce the finite controls with
`python 04-computation/experiments/crossroads223_20260926_geometry_growth.py`.
The matching `.out` includes the script SHA256. The argument below is the
proof; the small controls check its carry and counting mechanisms.

## 1. Inheritance and the new question

The parent requested an actual growing-horizon consequence of the
fixed-word S-unit argument in
`05-knowledge/results/crossroads223_20260926_geometry.md`, section 6.
The older OPEN G2 carrier in
`05-knowledge/results/crossroads_20260926_geometry.md` keeps both slots
of the transition pair (-3m,3m+1); its generic eight-node hostile shows
that distinct sources and successors alone do not ensure independence.

The new operation is to replace the global prime bank by the bank of
carries incident to a single node, then use the factorial cost of many
distinct primes. The concept board is: exact valuation cylinder,
incident resultant, S-unit solution count, source-interval length, and
chronological versus statistical quantifiers. The original 223 return
is one structured source of such carry primes, but the theorem below
does not single out that prime.

Incoming commit `a99d62c0f`, read without checkout mutation, contains
THM-4487,
`01-canon/theorems/THM-4487-dip-spectrum-entropy-curve-and-sharpness-of-thin-divergence.md`,
and `05-knowledge/results/constants_atlas_20260926_recurrent_numbers.md`.
The relevant common mechanism is an exact parity-word count followed by
a binomial tail, not equality of recurrent decimal constants. This proof
uses the established Terras bijection (THM-4476), and does not depend on
THM-4487's new lower bound. A false gamma=1 lower-bound display in the
incoming file was independently identified and sent to the parent;
the parent supplied the separately audited prepend-one-odd-step repair.
That correction is owned by the parent and is not repeated as a live
dependency here.

## 2. Exact finite-horizon bound

For a positive odd integer n, put m_0=n and

```
m_(j+1)=(3m_j+1)/2^k_(j+1),
k_(j+1)=v_2(3m_j+1),    S_j=k_1+...+k_j.
```

Let B_N(I) count odd n in a set I for which the N positive integers

```
3m_0, 3m_1, ..., 3m_(N-1)
```

are multiplicatively dependent. For N>=2 define the entirely explicit
integers

```
P_N = 96^(N(N-1)/2),
a_N = max{a>=0 : a! < P_N},
E_N = binom(3N,N) binom(N,2) 2^(16a_N+40).
```

**Theorem.** For every N>=2 and every interval I consisting of H
consecutive positive integers,

```
B_N(I) <= E_N + (H/2)(27/32)^N + (27/4)^N.               (1)
```

The estimate is uniform in the location of I. It counts odd starting
integers, which number H/2+O(1). Full independence of these first-slot
values implies independence of the N actual transition pairs
`(-3m_j,2^k_(j+1)m_(j+1))`; hence (1) also bounds paired-rank failure.
No injectivity condition is imposed on the selected words or their
integer realizations. Repeated-source trajectories are included among
the possible exceptions.

### 2.1 Prime support local to a node

Fix an exact positive exponent word w=(k_1,...,k_N) with S_N<=3N. Its
positive integer sources form one residue class modulo 2^(S_N+1).
Parametrize them by n=a+2^(S_N+1)t, t>=0, with a the least positive
representative. Write m_i=L_i(t)=A_i t+B_i. Then

```
A_i=3^i 2^(S_N+1-S_i),
B_i=(3^i a+C_i)/2^S_i,
C_i=sum_(ell<i) 3^(i-1-ell) 2^S_ell.
```

The quantities C_i/3^i strictly increase, so the roots of the L_i are
pairwise distinct. Thus every determinant

```
D_ij=A_j B_i-A_i B_j,    0<=i<j<N,
```

is nonzero. If C_(i,j) denotes the carry of the subword from node i to
node j, then

```
D_ij=-3^i 2^(S_N+1-S_j) C_(i,j).                         (2)
```

For each source-node index i, let S_i^prime be {2,3} together with all
primes dividing the incident carries C_(min(i,j),max(i,j)), j!=i.
This superscript distinguishes the prime set from the cumulative
valuation S_i. Write s_i=|S_i^prime|.

Suppose there is an integral first-slot relation
`product_i (3L_i(t))^c_i=1` with some c_i nonzero. If p outside
S_i^prime divides L_i(t), then p cannot divide any other L_j(t): a
common divisor would divide D_ij. Also p!=3. Taking the p-valuation
therefore forces c_i=0. Consequently every supported node L_i(t) is
an S_i^prime-unit. There must be at least two supported indices, since
each integer 3L_i(t)>=3.

For a supported pair i<j, the exact identity

```
x+y=1,
x=A_j L_i(t)/D_ij,    y=-A_i L_j(t)/D_ij                 (3)
```

has x an S_i^prime-unit and y an S_j^prime-unit. This uses (2) and the
fact that every A_i has prime support contained in {2,3}. The product
of these two rational unit groups has torsion-free rank s_i+s_j.

**CITED input.** Beukers--Schlickewei, Theorem 1.1, bounds solutions of
x+y=1 in the Q-closure of a rank-r finitely generated subgroup of
(C*)^2 by 2^(8r+8). See the verified
[author-hosted paper, Theorem 1.1, first page](https://webspace.science.uu.nl/~beuke106/s-units.pdf)
and the original Acta Arithmetica 78 (1996), 189--199
[publisher record](https://www.impan.pl/en/publishing-house/journals-and-series/acta-arithmetica/all/78/2/109078/the-equation-x-y-1-in-finitely-generated-groups).
Applying that theorem to (3), and using injectivity of t -> x, bounds
the exceptional parameters for this pair by 2^(8(s_i+s_j)+8).

This is the decisive improvement over the global bank: only two
incident banks enter the equation, rather than every pairwise carry.

### 2.2 The factorial cost of incident prime support

A subword of length r=j-i and total valuation K=S_j-S_i has carry C
coprime to 6. Its first term is odd and its last term is nonzero modulo
3. Since every exponent is positive, its partial valuations obey
S_ell<=K-(r-ell), and hence

```
C <= 2^(K-r) sum_(ell=0)^(r-1) 3^(r-1-ell)2^ell
  = 2^(K-r)(3^r-2^r)
  < 2^K (3/2)^r.                                       (4)
```

Let Q_i be the product of all N-1 carries incident to node i. The sum
of their valuation spans is at most (N-1)S_N<=3N(N-1). The sum of
their lengths is at most N(N-1)/2 (attained at an endpoint). Therefore

```
Q_i < 2^(3N(N-1)) (3/2)^(N(N-1)/2)
    = 96^(N(N-1)/2) = P_N.                              (5)
```

If t_i distinct primes divide Q_i, their product divides Q_i and is at
least t_i!. Thus t_i!<P_N, so t_i<=a_N and s_i<=2+a_N. By (3), each
fixed word has at most

```
binom(N,2) 2^(16a_N+40)                                 (6)
```

exceptional positive integer realizations, independently of their
heights. This count includes any source-node collision: equal nodes
have their common value dividing the corresponding determinant.
A node equal to 1 is a unit in every bank; the factor 3 in the first
slot prevents a one-index relation. Final-node collisions need not
cause rank failure and are not assumed absent.

There are exactly

```
sum_(S=N)^(3N) binom(S-1,N-1) = binom(3N,N)
```

positive words of length N with total valuation at most 3N. Multiplying
(6) by this number yields E_N. This is a union over every such word,
not a selection based on experimental injectivity.

### 2.3 Exact tail normalization on an arbitrary interval

Put L=3N. For an odd source, S_N>L exactly when the half-Collatz
parities at times 1,...,L contain at most N-1 odd terms. The initial
time-zero parity is fixed odd. Terras's bijection on L+1 bits therefore
makes this event a union of exactly

```
R_N=sum_(j=0)^(N-1) binom(L,j)
```

residue classes modulo 2^(L+1). Each such class meets an interval of
H consecutive integers at most H/2^(L+1)+1 times. The elementary
binomial estimate obtained by weighting at z=1/2 gives

```
R_N <= 2^N(1+1/2)^(3N) = (27/4)^N.
```

Thus the number of large-valuation sources in I is at most

```
(H/2^(3N+1)+1) R_N
 <= (H/2)(27/32)^N + (27/4)^N.                           (7)
```

Combining (6) and (7) proves (1). This uses a residue-class count,
not an assumption that a deterministic orbit has independent digits.

## 3. Growing horizon and uniform source density

The elementary integral bounds for log(a!) give
`log_2(a!)=a log_2(a)-O(a)`. Since
`log_2(P_N)=(5+alpha)N(N-1)/2`, where alpha=log_2(3), inversion yields

```
a_N ~ [(5+alpha)/4] N^2/log_2 N,
log_2 E_N ~ 4(5+alpha) N^2/log_2 N.                      (8)
```

The binomial factor and the factor binom(N,2) contribute only O(N)
and O(log N), respectively, to the second logarithm.

For H tending to infinity let Q=log_2 H and choose

```
N(H)=floor( sqrt(Q log_2 Q) / 8 ).                       (9)
```

Then (8) gives

```
E_N = H^(beta+o(1)),
beta=(5+log_2 3)/8 = 0.82312031259... < 1.
```

The term (27/4)^N is H^o(1). Consequently, uniformly for every interval
I of H consecutive positive integers,

```
B_N(I) <= (H/2)(27/32)^N + H^(beta+o(1)) = o(H).         (10)
```

**Therefore a proportion tending to one of the odd sources in every
sufficiently long interval has full first-slot multiplicative rank,
and hence full paired rank, through the growing number (9) of odd
transitions.** The o(1) and the threshold depend on H, not on the
location of the interval.

In particular, full rank through N(n) holds for a set of odd starting
integers of natural density one relative to the odds. Indeed, N(n) is
nondecreasing for large n, and dependence of a prefix persists when
more columns are appended. Failure at N(n) for n<=H is therefore
contained in failure at N(H), which is bounded by (10).

More generally any fixed constant
`0<c<1/sqrt(8(5+log_2 3))=0.1377775...` can replace 1/8 in (9), with
power 8(5+alpha)c^2<1. This is not claimed to be an optimal constant or
an optimal horizon. Without the factorial refinement the same local
bank argument gives a square-root-logarithm horizon; the global bank
alone gives only a cube-root-logarithm horizon.

The exact finite bound (1), rather than the limiting decimal exponent,
must be used for a concrete numerical window. No practical small
threshold is claimed. For each fixed N the rank-failure source set
also has upper Banach density zero: first truncate at an arbitrary
L>=N, use finitely many fixed-word exceptions and the uniform Terras
tail, let H tend to infinity, then L tend to infinity.

## 4. Conditioning on the no-descent candidates

There is a sharper connection to the gamma=1 endpoint of THM-4487.
Let D(X) be the positive sources n<=X whose half-Collatz trajectory
never falls below n during times 0,...,floor(log_2 n). Define

```
P_N^(2)=24^(N(N-1)/2),
a_N^(2)=max{a : a!<P_N^(2)},
E_N^(2)=binom(2N,N) binom(N,2) 2^(16a_N^(2)+40).
```

**PROVED exact bound.** For every N>=2, the number of sources in D(X)
whose first-slot rank through N odd transitions is deficient is at most

```
4^N + E_N^(2).                                         (11)
```

Proof: discard n<4^N, costing at most 4^N sources. For the others,
floor(log_2 n)>=2N. If S_N>2N, at most N odd letters occur among the
first 2N half-steps (including the time-zero odd term). The standard
affine carry bound gives

```
T^(2N)(n) <= (3/4)^N n + (9/4)^N - 1
          <= [(3/4)^N+(9/16)^N]n - 1 < n.
```

The last coefficient is at most 225/256 for N>=2. This contradicts
membership in D(X). Thus S_N<=2N. Repeat the local-bank argument with
2N in place of 3N: the incident product bound becomes 24^(N(N-1)/2),
and the number of words becomes binom(2N,N), proving (11). In this
conditioned statement there is no large-valuation tail to pay for.

With the growing horizon (9) and H=X, (11) is

```
X^(beta_2+o(1)),
beta_2=(3+log_2 3)/8 = 0.57312031259... .                (12)
```

The proof of (11)--(12) uses only the standard affine carry bound and
the S-unit result, not THM-4487's new lower bound. Comparing it with
the gamma=1 count `#D(X)=X^(h(log_3 2)+o(1))`, once the incoming
THM-4487 endpoint repair is incorporated, gives a useful limitation:
**almost all of the entropy-sized no-descent candidates already have
full first-slot rank through this growing horizon.** The exponent
0.5731203 is smaller than h(log_3 2)=0.949956. Imposing this rank
condition therefore does not remove the leading no-descent entropy
count. This is a comparison of explicitly different predicates, not
a new descent theorem.

**Sharper asymptotic cutoff.** Put K_N=ceil(alpha N)+1. The same reasoning
gives an exact alternative to (11). Define

```
P_(N,K)=2^((N-1)K-N(N-1)/2) 3^(N(N-1)/2),
a_(N,K)=max{a : a!<P_(N,K)},
E_(N,K)=binom(K,N) binom(N,2) 2^(16a_(N,K)+40).
```

For every N>=2, the rank-deficient sources in D(X) number at most

```
2^(K_N+1)+E_(N,K_N).                                   (13)
```

Indeed discard n<2^(K_N+1). The no-descent window includes time K_N.
If S_N>K_N, the first K_N half-steps contain at most N odd letters,
so their multiplier is at most 3^N/2^K_N<=1/2. Their carry is less
than (3/2)^K_N<=2^K_N<=n/2. This forces strict descent, a contradiction.
Thus S_N<=K_N. The same incident-bank proof, with that exact integer
budget, gives the displayed P and E.

Since K_N=alpha N+O(1),

```
log_2 P_(N,K_N)=((3alpha-1)/2)N^2+O(N).
```

At the growing horizon (9), (13) is consequently

```
X^((3log_2 3-1)/8+o(1)) = X^(0.46936093770...+o(1)).     (14)
```

Both (13) and (14) are proved without the lower bound of THM-4487.
Only the interpretation as a vanishing fraction of its entropy-sized
no-descent set uses that comparison. The script computes K_N exactly
as `bit_length(3^N)+1`, rather than rounding a floating-point logarithm.

## 5. Checks, interpretation, and stopping boundary

**FINITE-EXACT controls.** The script enumerates all 3,597 positive
words with N=2,3,4,5 and total<=3N. For every incident carry bank it
checks the coprimality, integer form of (4), span and length budgets,
(5), factorial bound, and exact word count. The largest observed
incident bank has seven carry primes; this finite statistic is not
used as a uniform bound. An independent direct orbit census checks
all odd residues modulo 2^(L+1) for N=1,...,5 and L=2N,3N, agreeing
exactly with the binomial tail. No random samples or decimal threshold
tests certify these checks.

This is a genuine transfer beyond the fixed-word finite-exception
statement: the word length now grows, every bounded-total word is
included, the discarded words have an explicit uniformly small source
density, and the conclusion holds in every source interval regardless
of its location. The arithmetic carriers are subword resultants, so
the observation around the 223 return has acquired a precise general
use rather than being elevated as a privileged prime.

It does not improve THM-4476's thin-divergence exponent. The total
exception bound in (10) is H times a stretched-exponentially small
factor in log H, not a fixed-power bound H^delta with delta<1; its
leading exponent is still 1. Moreover rank failure is a different
predicate from no descent. A single hypothetical exceptional orbit
may visit a source set of zero density forever. The argument extends
in form to other fixed odd multipliers, including maps for which
divergence is expected, after adjusting the fixed prime banks and
constants. Thus statistical full rank cannot alone distinguish the
convergent Collatz dynamics from such hostile maps.

The remaining OPEN target is a uniform chronological obstruction to
an actual dependent trajectory, or an additional consequence of full
rank that controls integer height. Neither is supplied by (10).

An explicit hostile shows why the all-words height quantifier cannot
be silently strengthened. For every integer a>=2,

```
n_a=(4^a-1)/3 -> 1 -> 1
```

is an actual positive odd trajectory. Its first three transition pairs
include two copies of (-3,4), so their rank is deficient at arbitrarily
large starting heights. The first valuation 2a grows, and the first
three valuations sum to 2a+4; these sources eventually enter the tail
discarded in (7). This is not a counterexample to injective-orbit G2.
It is a sharp warning that fixed-word finiteness and statistical
independence cannot be replaced by independence for every sufficiently
large start across all words.

For literature context, fixed-translates multiplicative independence
already has a developed theory: Dubickas--Sha, *Multiplicative
dependence of the translations of algebraic numbers*, Theorem 3.1,
[primary arXiv manuscript](https://arxiv.org/pdf/1608.05458), proves
eventual independence of distinct fixed algebraic translates and
identifies that statement as a special case of earlier work. Its
equal-translation formulation is not directly substituted for our
unequal affine slopes. No priority claim is made for fixed-word
finiteness; the result developed here is the explicit word-uniform
count, growing horizon, and source-interval quantifier.
