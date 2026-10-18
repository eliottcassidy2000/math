# Positive ROOT transport, a uniform logarithmic-depth remainder, and the remaining clock law

2026-10-17, clock_duality_oct15. **PROVED, scoped, self-audited:** the positive
band transport, its Syracuse factorization, the one-time source-rounding
bound, and the explicit depth/scale remainder below. **FINITE-EXACT:** tariff
brackets through depth10 and marked-clock checks. **OPEN:** the required
cancellation rate, ROOT energy, small-uncertainty entry, and Collatz.
No independent-referee or novelty claim. The uniform remainder is a genuine
improvement over an unspecified fixed-depth constant; it is not the missing
uniform current estimate.

## Inheritance and concept board

Anchor: cancellation at growing depth with the original integer source
measure. Niche: a golden pentagon shared by two geometric coefficient axes
([companion](shared_golden_contact_20261017g.md)). Wildcard: change the sign
gauge before changing the norm. Board: **source row / affine carry / ternary
conductor / band exponent / physical clock / shared geometric port**.

Closest proved mechanism: R8 of [marked tariffs](root_braid_logistic_20261017f.md).
Its fixed-depth coefficient recurrence and last-digit balance are inherited.
Canonical hostile: the ROOT-normalized doubling-and-braid profile with
current16N/7 in R4. Corrected near miss: squared branch norms are not the
norm of their sum, as explicitly corrected in
[THM-4658, two-octave-memory-energy-descent-critical-lines](../../01-canon/theorems/THM-4658-collatz-two-octave-memory-energy-descent-critical-lines.md).
Least-used sidecar: the terminal braid index, which restores the integer
clock after a convenient redistribution of its even part.

The inheritance search also recovered [refinement debt R2--R4](collatz_refinement_debt_20261016a.md):
transporting an error once is already a proved method for the Gram remainder.
Section3 applies that method to the affine rounding of the initial source;
it does not rename the old Gram theorem. The earlier
[five mirrors](collatz_five_mirrors_20260929.md), especially its cost-split
Fourier computation, already distinguishes total-cost conditioning from an
unmarked marginal. Its Fourier/renewal results are stronger than a numerical
guess from our ten tariffs, but concern a different observable.

Incoming work through a79870665e was read. The audited
[full ROOT census](root_current_ground_truth_contact_rotations_20261017.md)
now has exact integer energy numerators. Its positive covariance examples
remain hostile to a fixed-negative-covariance argument. MISTAKE-626's
Paley-order correction reinforces the use of explicit witnesses rather
than an unchallenged search. Its Farey sweep retains the actual event
position; the analogous retained coordinates here are band and conductor.
No ranking theorem is used as a Collatz premise.

## 1. A sign gauge makes the ideal transport positive

Use all definitions from R8: q_d=3^(d+1), bands
M_d/2^(b+1)<t/N<=M_d/2^b, M_d=4(3/2)^d, and K_d(r,b), zero for b<0.
Write A=a+2j, where a=1 for A odd and a=2 for A even. The forward residue
map corresponding to that inverse fibre is

    F_A(r)=(3r+c_A)/2^A mod 3q,
    c_A=(3+2^a-2^A)/3,                r mod q.                 (1)

The numerator defining c_A is divisible by3. For each A, F_A is an
injection from residues modulo q to residues modulo3q. It is precisely
B^(-j)((3r+1)/2^a), B(r)=4r+1. All divisions in residue coordinates use
units; this is not an assertion that every branch is legal at one integer.

**P1 (PROVED).** Define a finite signed measure on (r,b) by

    nu_d(r,b)=(-1)^(b+d) K_d(r,b)/(4*2^(d+b)).                (2)

Then nu_d is obtained from nu_(d-1) by the positive probability kernel

    (r,b) -> (F_A(r), b+A-1),       probability 2^(-A), A>=1. (3)

Its total variation is exactly gamma_d, and hence

    2^(-d) <= gamma_d <= gamma_(d-1) <= gamma_0=16/7.         (4)

Proof: in the old recurrence the sign is (-1)^a. Since
b_new=b_old+a-1+2j and d_new=d_old+1, the sign in (2) makes every branch
positive. The denominator changes by2^A, giving (3). Branch probabilities
sum to1. At b=0 only A=1 can occur, so the b=0 mass alone is2^(-d): the
initial row (-1,-1,2) has absolute mass1 after division by4, and each such
branch is injective. This proves the lower bound as well as contraction.

The total signed mass is zero, even separately at each b, by R8's
three-lift balance. Thus gamma_d is the variation of a difference of two
positive measures of equal mass8/7. This is a representation with no
norm amplification at any depth. Strict decay requires overlap of these
particular measures; positivity alone does not supply it.

**Sharp one-step hostile.** Put delta_0-delta_1 in a single band at modulus3.
After one transition, the band determines A, and within each A the residue
map is injective. Its norm stays exactly2. It already has the required
three-lift sum zero. Therefore a uniform one-step factor below1 fails even
on balanced band measures. The inherited initial band distribution carries
additional information which a successful invariant class must retain.

## 2. Exact connection to Syracuse, with the clock boundary restored

Let U_a(r)=(3r+1)/2^a. Define positive kernels on residues and band exponent:

* A: U_1 with probability2/3 and band increment0; U_2 with probability1/3
  and band increment1.
* C: B^(-j) with probability(3/4)4^(-j) and band increment2j, j>=0.

Composition acts on measures from right to left. Then

    P=C A,       S=A C,       P^d=C S^(d-1) A  (d>=1).       (5)

Indeed U_a B^(-j)=U_(a+2j) exactly, and
(2/3 or1/3)*(3/4)4^(-j)=2^(-(a+2j)). Hence S is the ordinary random
Syracuse affine kernel with geometric valuation and retained exponent
increment A-1. This is an algebraic bridge, not a shared-vocabulary analogy.
Residues, probabilities, and exponent sums are preserved. Forgetting the
last coordinate can create cancellation unavailable at a fixed band.

For a genuine integer path write x_i=B^(j_i)(s_i), with s_i primitive and
j_0=0 for a primitive source. If a_i=v2(3s_(i-1)+1), then

    A_i=a_i+2j_i,
    physical shortcut clock through d odd transitions
      =sum_(i=1)^d A_i - 2j_d.                              (6)

For a nonprimitive initial source add2j_0. Thus A_i is the homogeneous
scale exponent, not by itself the physical clock of that edge. Its parity
is correct, and the terminal correction restores the full clock. We retain
the endpoint and j_d in every actual-path check. Continuing the even ROOT
cycle in these algebra checks is harmless to signs; it is not a claim about
first-hit time. A first-hit controller stops and records its actual prefix.

**CITED, precisely scoped.** [Tao, 1909.03562v5, Proposition1.14 and Remark1.16](https://arxiv.org/html/1909.03562v5#S1.SS4)
prove fine-scale mixing of the Syracuse **offset**, and distinguish a stronger
mixing statement for the whole affine map. The latter is not proved there.
Equation(5) identifies the relevant classical operator, but our band readout
retains the exponent. We therefore cannot substitute offset-only mixing
for the joint estimate. The published theorem concerns logarithmic density,
not universal ROOT entry. No theorem of that paper is used to assert decay
of gamma_d here.

## 3. Round the source once, then transport its error

Let S={positive odd s:s!=5 mod8}. In H_N(s) replace B^j(s) by4^j s only
in the two size tests, retaining the residue phase s+j mod3. Call the result
Hhat_N(s). Both are finite vectors on the same literal integer sources.
Let L delta_s=(-1)^(v2(3s+1)) delta_(P(s)) be exact primitive transport.
This has l1 norm at most1 and fixes delta_1.

**P2 (PROVED).** For all N>=1 and all d>=0,

    ||L^d(H_N-Hhat_N)||_1 <= ||H_N-Hhat_N||_1
       <=4 floor(log_4(4N)).                                (7)

Proof: write each summand as

    w(s+j) 1_(B^j(s)<=4N)
      +[-3chi(s+j)-w(s+j)] 1_(B^j(s)<=2N).

Each coefficient has absolute value at most2. After dividing by4^j, each
threshold moves by(1-4^(-j))/3<1/3, so at most one integer source is crossed
per threshold. There are floor(log_4(4N)) relevant positive j. Transport
contracts the entire error vector, rather than billing it again at each
step. This is uniform for arbitrary changing common authenticated signed
maps too, by the same column-norm argument.

For N=2^k the result sharpens to the exact formula

    H_(2^k)-Hhat_(2^k)=e_k delta_1,
    (e_0,...,e_5)=(-2,-1,1,-1,1,2), e_(k+6)=e_k.            (8)

Every integral threshold except1 is even and cannot be a primitive odd
source. The threshold1 occurs at exactly one of the two cutoffs; evaluating
its residue gives the table. Thus its error remains a known ROOT atom at
every depth and can be corrected exactly.

Over an octave X<=N<2X, (7) also bounds the l2 norm of any scalar error
evaluation by4 sqrt(X) floor(log_4(8X)), uniformly after transport. This
does **not** compare the transported discrete Hhat_N with its homogeneous
density limit. That separate arithmetic comparison is addressed next.

## 4. An explicit remainder valid at logarithmically growing depth

Let A^(d)(N)=||L^d H_N||_1, the same complete odd-transport envelope as R8.
It includes the ROOT coefficient. Put

    m=ceil(log_2(4N)),  K=m+2d,
    C_d=(4/3)(3/2)^d-1,
    W(K,d)=sum_(k=0)^floor((K-d)/2) binom(K-2k,d).            (9)

**P3 (PROVED).** For every N>=1, d>=0 and 1<u<2,

    |A^(d)(N)-gamma_d N|
      <=36(C_d+2) W(K,d)
        + N (5/2)/(1-u^2/4) [u/(2-u)]^d u^(-K).            (10)

In particular, for **every** d<=floor(m/16),

    |A^(d)(N)-gamma_d N| <=87(m+1) 3^(m/4)
      =O(log(N+2) N^(log_2(3)/4)),                          (11)

with one constant independent of depth. The exponent log_2(3)/4 is about
0.39625, strictly below1/2. This replaces the old O_d(log^(d+1)N) by a
usable uniform statement on a growing depth range.
There is a slightly longer alternative horizon: for d<=floor(m/12),

    |A^(d)(N)-gamma_d N| <=87(m+1)(231/4)^(m/12).           (11b)

Its exponent log_2(231/4)/12 is still below1/2, because231/4<64.

### Proof of the word budget

An inverse word has initial braid k and d positive exponents A_i=a_i+2j_i.
Its homogeneous slope from terminal t to original odd source is
alpha=2^(2k+sum A_i)/3^d. Each actual primitive edge satisfies

    2^A(t+1/3)=3(s+1/3)+2^a/3 <=4(s+1/3),

since s>=1 and a<=2. Also B^k(s)+1/3=4^k(s+1/3). Hence every actual word
with original source<=4N and terminal>=1 has

    2^(2k+sum A_i) <=4^d(3N+1/4) <=4^d*4N <=2^K.          (12)

Every homogeneous word contributing at an integer t>=1 obeys the same
cutoff, since3^d<=4^d. For fixed k, the number of positive d-tuples with
sum A_i<=K-2k is binom(K-2k,d). The d=0 convention gives one empty tuple.
This proves W without a hidden depth constant.

### Affine displacement and lattice counting

One inverse edge has slope alpha_i=2^A/3 and normalized offset

    beta_i/alpha_i=(1-(2^a+3)/2^A)/3,
    |beta_i/alpha_i|<=1/2,       alpha_i>=2/3.

The initial B^k has normalized offset at most1/3. Composition, in the
actual order, therefore gives

    |beta/alpha| <= (1/2)sum_(h=0)^(d-1)(3/2)^h
                    +(1/3)(3/2)^d = C_d.                  (13)

An inverse word's integrality selects one terminal class modulo3^d; its
original coefficient phase splits that into three classes modulo3^(d+1).
Each has three primitive classes modulo8. Thus at most **nine** progressions
are involved per word. On each, the two cutoff coefficients have modulus
at most2. Moving both thresholds by at most C_d costs at most4(C_d+1)
per progression. Summing gives36(C_d+1)W for actual versus homogeneous
coefficients; taking their absolute values cannot enlarge this difference.

For the homogeneous truncated sum, counting each progression differs from
its density integral by at most1 for each prefix. Apply summation by parts
to the absolute value of the summed step function. Its variation is bounded
by the sum of the word variations: at most4 per progression per word.
This costs another36W. This argument merges coefficients before taking
absolute values and counts jumps rather than assuming fixed signs.

### The remaining continuous tail

Words beyond K have no integer contribution. They can still contribute to
the density integral at terminal size below1. For one initial braid k,
the sum of the absolute prefix/rim integrals over the three residue phases,
with their primitive density1/8, is exactly(5/2)4^(-k). For one fixed
following word it is multiplied by2^(-sum A_i). This follows either from
the change in interval length and progression density, or directly from
the inverse slope above. The unsigned mass of all omitted words is at most

    sum_(2k+sum A_i>K) (5/2)4^(-k) 2^(-sum A_i)
    <= (5/2)/(1-u^2/4) [u/(2-u)]^d u^(-K).                (14)

The geometric sums converge for1<u<2. The reverse triangle inequality
compares the full and truncated density norms. Multiplying (14) by N and
adding the preceding two bounds proves (10).

### Uniform constants in (11)

For d<=m/16, K<=9m/8 and

    W<= (K/2+1) binom(K,d) <=(m+1)54^(m/16).

Use binom(K,d)<=(eK/d)^d and e<3; the bound is increasing in d up to m/16
after replacing K by9m/8. The d=0 case is immediate.
Also C_d+2<=(7/3)(3/2)^d. The first term of (10) is at most
84(m+1)81^(m/16)=84(m+1)3^(m/4).
Choose u=7/4. Since N<=2^(m-2), the second term is at most

    (8/3)[(8/7)(16/7)^(1/16)]^m <=(8/3)3^(m/4).

The last inequality follows by taking sixteenth powers:
(8/7)^16(16/7)<81. Adding proves (11). No numerical data enters this proof.
For (11b), replace9m/8 by7m/6 and use e<11/4, obtaining
binom(K,d)<=(77/2)^(m/12). The carry factor changes this to
(231/4)^(m/12). The same tail is smaller since
(8/7)^12(16/7)<231/4. The constant87 still suffices.

## 5. What is closed, and what must change next

The old unspecified depth dependence is closed on d<=m/12. Source rounding
can separately be carried forever at the price in (7). The band transport
is positive and norm-bounded at arbitrary depth. None of these facts proves
a decreasing uniform bound for the actual ROOT current.

**P4 (PROVED, actual-source obstruction).** At every N>=1 and d>=0,

    A^(d)(N) >= N 2^(-d)-12.                                (16)

This does not require the asymptotic comparison. Write R=(3/2)^d. Above
the terminal threshold R(2N+1)-1, the only possible contributing word is
k=0, A_1=...=A_d=1. Indeed every other word has homogeneous slope at most
3^d/2^(d+1), including the initial braid factor, and forward affine offset
at most R-1: each branch constant c_A/2^A is at most1/2. Its original
source is at most4N, so it cannot exceed that threshold.

The unique outer word sends s to R(s+1)-1 and has physical clock exactly d.
Put t=3^d h-1 and s=2^d h-1. Its full outer-band mass is exactly

    E(N,d)=sum_( (2N+1)/2^d < h <= (4N+1)/2^d,
                  3^d h-1 in S) |w(2^d h-1)|.             (17)

All inverse steps are legal: starting from a primitive terminal in this
progression, each inverse exponent1 gives an odd primitive predecessor.
There are nine admissible classes of h modulo24, with total weight12.
The interval has length2N/2^d. Counting those progressions gives
|E(N,d)-N2^(-d)|<=12, proving (16). No other word can cancel these columns.

Thus an absolute-envelope bound O(sqrt(N)) requires at least
d>=(log_2 N)/2-O(1). Our proved small-error horizon is about(log_2 N)/12.
This is a limitation of the envelope, not a lower bound for the signed ROOT
current. In the outer strip, different terminal signs could still cancel
when actually evaluated. Cross-band or cross-row signed information could
therefore beat the envelope. Exponentiating the first tariff ratio would
again be invalid. Equations(16)--(17) identify the precise slow packet that
a new norm or adaptive stopping policy should isolate.

The next precise questions are:

1. Can the low-clock, fast-growth words be carried as a separate positive
   packet, while merging the other words before paying lattice error? The
   present binomial W counts every word even when many have the same jump.
2. Can one bound the **joint** offset/exponent law in (5), with the input
   nu_0 and three-lift cancellation retained? A marginal estimate loses the
   coordinate used in the current evaluation.
3. Can the transported octave mean/detail Gram from refinement-debt R2--R4
   exploit overlap between different source rows before rowwise absolute
   values are taken? P3 now gives a uniform error budget for that comparison
   on a nontrivial depth range, but supplies no mean-current bound itself.

These are sufficient-target questions, not claims that the corresponding
inputs follow from positivity. Single-band data rule out a generic one-step
gap. Finite tariffs alone do not decide an asymptotic rate.

## 6. Evidence and reproducibility

The [manifest](../../04-computation/experiments/root_positive_transport_20261017g/manifest.json)
lists exact scripts, receipts and hashes. The positive-table computation
keeps b<120. If t=3/2, the omitted variation is bounded by

    M_0(t)(2-t)^(-d)t^(-120),
    M_0(t)=sum_(b=0)^5 [sum_r|K_0(r,b)|/4](t/2)^b
                              /[1-(t/2)^6].               (15)

This follows from (3): each increment A-1 has generating function1/(2-t).
Thus the reported intervals are rational enclosures, not fitted decimals.
The first four reproduce the inherited exact rational tariffs. Every
computed clock row and last-digit fibre sums to zero. Independent inverse
recurrence checks and literal shortcut-clock replay are recorded separately.

| d | gamma_d, rounded from a rational interval |
|---:|---:|
|4|1.079843060221|
|5|0.967085959924|
|6|0.877335224707|
|7|0.802929800548|
|8|0.738891220381|
|9|0.682907714873|
|10|0.633171196420|

The certified tail at depth10 is below10^(-16). Dropping the band before
taking variation at that depth instead gives about0.151341378850, a
strictly smaller and inapplicable substitute for gamma_10. This explicitly
measures one loss from forgetting the clock coordinate. The audit compares
43560 signed coefficients through depth4 by the inverse recurrence,
19296 inverse words across25 source/depth cases, and94838 literal shortcut
steps. Normal and optimized runs agree byte-for-byte. These finite checks
support the formulas; the uniform quantifiers in P1--P4 rest on their proofs.

Source/target/map ledger: H_N -> L^d H_N preserves evaluation against every
compatible bounded ROOT profile; source arrays and integer clocks restore
the forgotten paths. K_d -> nu_d preserves the envelope exactly but changes
the evaluation gauge; b restores it. CA -> C S^(d-1) A preserves the complete
joint law, not a deterministic integer's randomization. The geometric
companion fixes literal ports and tests complete induced edges. Its colour
extensions and this kernel's contraction are separate theorems; their shared
research move is to preserve the interface before simplifying the interior.
