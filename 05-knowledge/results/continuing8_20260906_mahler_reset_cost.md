# Mahler reset cost, a sharp terminal clock, and a failed Rule 30 quotient

**Status: PROVED + FINITE-EXACT; [independent audit accepted](continuing8_20260906_mahler_reset_cost_audit.md).** A fixed nonzero safe carry prefix followed by zero resets has an exact least-launch cost spectrum, determined by its reduced ternary denominator. The same coordinate gives a sharp bound on the zero-reset runway after the actual last native bit. An explicit infinite family attains that bound and then provably rejects. Uniform phase on the finite cost cycle has a limiting two-head fair-coin waiting law, with an explicit error; this is not a randomness assertion about one Mahler orbit.

The attempted Rule30 bridge fails on an exact native object: two ordinary terminal prefixes have the same saturated inverse monoid, follower state, depth and native bit length, but opposite next reset decisions. Mahler's Z-number problem and the named Rule30 prizes remain OPEN.

## 1. Inheritance and the proposed map

The current proved inputs, with exact identities and type boundaries, are:

- [THM-4072-mahler-safe-terminal-fibre-product-and-finite-state-obstruction.md](../../01-canon/theorems/THM-4072-mahler-safe-terminal-fibre-product-and-finite-state-obstruction.md): the strict-safe follower graph, rooted carry/native isometry, and affine cylinder law A=r_m+2^m k implies T^m(A)=u_m+3^m k.
- [THM-4074-mahler-denominator19-postterminal-arbitrary-delay.md](../../01-canon/theorems/THM-4074-mahler-denominator19-postterminal-arbitrary-delay.md): reachable postterminal programming and arbitrary finite delay, with a changing denominator19 launch prefix.
- [THM-4077-mahler-denominator19-2adic-tangent-full-shift-isometry.md](../../01-canon/theorems/THM-4077-mahler-denominator19-2adic-tangent-full-shift-isometry.md) and [THM-4082-mahler-renormalized-linear-chart-and-exact-bit-defect.md](../../01-canon/theorems/THM-4082-mahler-renormalized-linear-chart-and-exact-bit-defect.md): parameter/output termination differ; a limiting isometry and every finite program do not produce one fixed safe ordinary orbit.
- [THM-4204-rule30-debruijn-reset-and-dyadic-prefix-saturation.md](../../01-canon/theorems/THM-4204-rule30-debruijn-reset-and-dyadic-prefix-saturation.md): the41-element spatial inverse monoid and its right-absorbing rank-one resets.
- [THM-4210-rule30-lossless-dyadic-block-current-cartier-tree.md](../../01-canon/theorems/THM-4210-rule30-lossless-dyadic-block-current-cartier-tree.md): the forward current retains both physical children and a transverse channel; its physical quadratic admissibility is not supplied by a finite inverse monoid.

THM-2228, `THM-2228-mahler-three-halves-carry-tail-and-integral-stabilization.md`, supplies the underlying exact strict-real-tail/ordinary-state criterion; THM-3848 supplies its finite safe language. Targeted searches of those current theorems and the reset/height/cost routes found the affine cylinder and arbitrary-delay results, but no prior fixed-prefix cost spectrum or sharp clock law stated below. These are corollaries developed from that carrier, not a claim to have rediscovered its construction or the elementary binary expansion of a rational number.

The canonical hostile is THM-4072's A8/A13: even a terminal follower state with a native-one marker misses the next acceptance decision. The corrected near miss is treating a reset of one observer as a reset of the entire native state. The least-used sidecar here is the reduced denominator of the finite **oriented** real tail, coupled to its precise native address.

The Anchor remains the separate open LRC14 graph-consumer lane. The Niche is the fixed Mahler zero-reset ray. The Wildcard compares its finite-cycle distribution with a simple fair-coin waiting distribution. The live board is:

| Object | Operation | Preserved predicate | Lost coordinate / required sidecar |
|---|---|---|---|
| Mahler safe follower | Zero reset at state0 | Continued finite safety | Native affine exit state |
| Native carry cylinder | Least ordinary completion | Exact finite carry prefix | Height growth and last native bit |
| Finite oriented tail D/N | Append zero carries | Same rational tail | Phase b_n of its ternary unit cycle |
| Rule30 inverse word | Multiply its boundary matrices | Cyclic spatial predecessor count | Temporal/native continuation after a reset |
| Finite unit cycle | Choose its phase uniformly | Exact combinatorial distribution | No sampling law for one selected infinite orbit |

## 2. Exact cost theorem

Use T(A)=ceil(3A/2), whose carry at a state A is A mod2. Fix a nonzero finite carry word w=c_0...c_(m-1), assumed follower-safe and ending at follower state0. Thus every appended zero is a reset. Let j be one plus the position of its last1, so1<=j<=m, and write

    C_m=sum_(i<m)c_i 2^i 3^(m-1-i),
    C_m/3^m=D/N in lowest terms.

Then

    N=3^j,   D=C_j,   gcd(D,3)=1,   0<D<N/2.          (1)

Indeed the last term of C_j is2^(j-1), and all preceding terms are divisible by3. Trailing zero carries multiply the numerator by3 without changing its reduced phase. The last inequality is the strict truncated suffix inequality at the initial position. This is the real phase of the formal completion w0^infinity, not automatically the real tail of an ordinary integer orbit.

For n=m+L, L>=0, put

    b_n=[2^(-n)D]_N in {1,...,N-1},
    A_n=(2^n b_n-D)/N.                               (2)

**Cost theorem.** A_n is the least positive ordinary integer whose first n carries are w0^L. Every ordinary completion and its exit state are exactly

    A=A_n+2^n k, k>=0,
    T^n(A)=3^(n-j)b_n+3^n k.                          (3)

In particular, for any such completion A<=X,

    2^n<=NX+D,
    n<=floor(log_2(NX+D)).                            (4)

The lower height bound A_n>=(2^n-D)/N is attained for infinitely many n, and (4) is then an equality with X=A_n.

**Proof.** The inherited cylinder congruence for w0^L is

    A=-D/N modulo2^n.

The representative in[0,2^n) satisfies NA+D=2^n b with0<b<N. ModuloN this forces b=b_n, giving (2). The representative is positive because a zero state has only zero carries. Substitution into the exact affine exit formula gives (3). Since b_n>=1 and A>=A_n, it gives (4).

The number2 generates the units modulo every3^j. For completeness, the binomial identity gives v_3(4^(3^a)-1)=a+1 by induction. Raising to a power not divisible by3 leaves this valuation unchanged. Therefore an exponent r with2^r=1 mod3^j must be even, r=2q, and3^(j-1) divides q; the exact order is

    P=2*3^(j-1).                                      (5)

Consequently b_n traverses every unit of Z/NZ once in a periodP, including1. This proves both the full cost spectrum and its infinitely repeated sharp boundary.

More explicitly, the **shifted** normalized cost is exactly periodic:

    (A_n+D/N)/2^n=b_n/N.                              (6)

The unshifted ratios A_n/2^n have these finitely many values as their subsequential limits, rather than being exactly periodic themselves. Their minimum and maximum limit values are1/N and1-1/N. The canonical native update is

    A_(n+1)=A_n+2^n(b_n mod2).                        (7)

To see this, division by2 moduloN sends an even b to b/2 and an odd b to(b+N)/2. Thus a longer zero-reset prefix frequently forces a new high native1. Exactly half the unit residues are odd, so one period has exactlyP/2 such native increments.

## 3. A sharp clock law at actual native termination

Let H_n=bit_length(A_n), so the last native1 occurs at positionH_n-1. Write

    b_n=2^r c,   c odd,
    K_n=n+r-H_n.

Then, for every n>=m,

    K_n=floor(log_2(N/c)),
    0<=K_n<=R:=floor(log_2 N).                        (8)

This identity uses the actual native bit length, not a newly chosen observation clock. Since N is an odd power of3 and3 does not divide c, N/c is never a power of2. The defining equality

    2^(n+r)c=NA_n+D

and0<D<N, together with2^(H_n-1)<=A_n<=2^H_n-1, give

    N 2^(H_n-1)<2^(n+r)c<N 2^H_n.

Taking binary logarithms and using that K_n is an integer proves (8).

At the observation depthn, the exit state is3^(n-j)b_n, so its next r carries are exactly0 and its following carry is1. If H_n>=m, then all carries from the actual terminal depthH_n through n+r-1 are zero resets at follower0. Hence K_n is exactly their number. This interpretation holds for every sufficiently large n; for example it holds after advancing the phase by one full period, since P>=R+1. When H_n<m, (8) remains the clock difference to the displayed odd state, but it must not be renamed a consecutive postterminal reset count without checking the earlier prefix.

The same bound has an all-orbit consequence. Suppose any positive ordinary A has terminal depthH>=m and its first H carries are w0^(H-m). Then A is the least representative at that depth, so (8) with n=H proves that its next zero-reset runway has length at mostR. Thus fixing the prefix and appending only resets cannot yield an unbounded postterminal runway; at least the address j of the last carry1 must change. In particular K postterminal zero carries in this class require

    2^K<3^j.                                         (9)

This does not contradict THM-4074: its denominator19 construction changes the launch prefix and that last-carry-one address as the desired delay grows.

## 4. Exact clock spectrum and a limited probabilistic bridge

Choose n uniformly modulo the periodP, at depths large enough that H_n>=m. This is a finite cyclic-phase experiment on different ordinary launches. It is not a stochastic model for the carry sequence of a fixed launch.

For k>=0 define

    a_k(N)=#{c odd:3 does not divide c,
                       N/2^(k+1)<c<N/2^k}.

An odd part c in this interval allows exactly k+1 unit residues b=2^r c<N. All of them give K_n=k by (8). Therefore the exact law is

    Pr(K_n=k)=(k+1)a_k(N)/P,   0<=k<=R.              (10)

In particular the maximumR occurs at precisely R+1 phases, corresponding to the powers1,2,...,2^R. Intermediate values can be absent: atN9, the complete spectrum is K=0 with frequency2 and K=3 with frequency4. The maximum alone would miss this arithmetic gap.

The count of odd integers not divisible by3 through an integerM>=0 is

    F(M)=floor((M+1)/2)-floor((floor(M/3)+1)/2).

Thus

    a_k(N)=F(floor((N-1)/2^k))-F(floor(N/2^(k+1))).   (11)

Uniformly in k, |a_k(N)-N/(3*2^(k+1))|<=2. Since P=2N/3, (10) gives, as j tends to infinity,

    Pr(K_n=k) -> (k+1)/2^(k+2),
    total_variation <= (R+2)^2/N.                    (12)

For the error bound, the pointwise discrepancy throughR is at most3(k+1)/N. The remaining mass of the limiting law is(R+3)/2^(R+2), so half the sum of discrepancies is at most[3(R+1)(R+2)+R+3]/(4N), which is bounded by the expression in (12).

The limiting law is exactly the elementary distribution of the number of failures before the second head in independent fair-coin tosses: the last toss is a head and one of the preceding k+1 positions is the first head. This supplies a genuine distributional correspondence with a declared measure. It does not identify the Mahler dynamics with coin tossing, assert independence of its native digits, or prove any infinite safe orbit. The exact finite residue cycle, rather than a sampling heuristic, supplies every probability in (10).

## 5. An explicit sharp family that always rejects

For arbitrary j>=1 and t>=2, set

    N=3^j, P=2*3^(j-1), R=floor(log_2 N), E=P t,
    A=2^(j-1)(2^E-1)/3^j,
    H=E+j-1-R.                                       (13)

These are ordinary positive integers. The order calculation in (5) proves the divisibility. Also E>=R+3, and

    A/2^H=(2^R/N)(1-2^(-E)) in(1/2,1),

so H is exactly its native bit length. The lower inequality follows from2^(R+1)-N>=1 and E>R+1. The simple bounds R<2j and E>=4j justify E>=R+3 for every stated parameter, including j1.

The complete carry word through the first rejection is

    0^(j-1) 1 0^(E-1) 11.                            (14)

Indeed after the first j-1 even steps the state is(2^E-1)/3. Its next state is2^(E-1), which emits E-1 zeros before reaching3^(E-1). The numberE is even, so the latter state and its next state(3^E+1)/2 are both odd, giving the final11.

The first1 is followed by at least two zeros, and100 is a safe reset loop. All intervening zeros reset at state0. Thus (14) is safe until the second1 of its final11, which rejects because the greedy boundary begins101. The exact rejection index, with carry indices starting at0, is

    H+R+1.                                           (15)

At the actual terminal depthH the state has exactly R subsequent zero carries, all reset edges. Hence the bound in (8) is attained for every j by infinitely many actual ordinary launches. Their postterminal delays are unbounded as j grows, but every displayed orbit provably rejects. These are explicit excluded integer parts for the Mahler candidate criterion, not a finite computation extrapolated to all integers.

For instance j2,t2 gives A=910, H10 and R3; the first rejection is at carry index14. The j1 family specializes to(2^(2t)-1)/3, with a long safe zero string but exactly one zero reset after native termination. The family distinguishes long reset evidence, the moving native clock, and an actual eventual failure.

## 6. The strict real endpoint is a separate carrier

For a fixed prefix in Section2 and any n>=m, define the rational companion

    alpha_n=A_n+D/N=2^n b_n/N.

Its first strict lower-half failure is exactly

    time n+v_2(b_n)+1,   fractional part1/2.           (16)

Up to timen, the affine carry identity says that the fractional parts are the oriented finite tails of w0^infinity, hence belong to[0,1/2). From timem onward those tails are zero. After timen the starting integer is3^(n-j)b_n, so the next r=v_2(b_n) phases are integers and the following one is an odd half-integer. This proves both firstness and the exact boundary value in (16).

The real companion must not be used to reject the ordinary carry orbit in general. At the carry indexn+r its next bit is1, read from follower0 and therefore accepted; its prefix of lengthn+r+1 remains safe even though the companion's phase at timen+r+1 is already forbidden. The infinite real tail forced by the actual future carry word need not equal D/N. Only the separately proved word (14) supplies actual rejection for the sharp family. This is the exact endpoint version of the denominator19 pseudo-orbit warning in THM-4074.

## 7. Why the proposed Rule30 reset quotient fails

There is a literal map on finite binary words: take the Rule30 inverse-boundary product M_w of THM-4204. Its generators, as maps on00,01,10,11, are

    g_0=(0,2,3,3), g_1=(2,0,1,1).

The word0010 yields the constant E_3, and E_3 A_0=E_3 A_1=E_3. Thus

    M_00100=M_0010011=E_3.

The first word is Mahler follower-safe and ends at state0; the second rejects on its final11. Consequently even the finite safe-carry predicate does not factor through this inverse spatial monoid. This is a failure of this concrete carrier map, not a theorem that no relation between the two subjects can exist.

A stronger native control retains the clock and actual positive termination:

| Carry prefix | Ordinary A | T^8(A) | Reduced zero-tail denominator | Next carry |
|---|---:|---:|---:|---:|
|00100000|180|4617|27|1, a match|
|00100100|148|3798|729|0, a reset|

Both prefixes have length8, follower state0, a seen native1, remaining native clock0, actual bit length8, and the same constant M_w=E_3. Their next reset events differ. Their exact cost spectra also have periods18 and486. The reduced phase and the native exit state cannot be recovered from the combined discarded quotient.

The distinction is structural. A Mahler reset block maps its native residual parameter k to nu+3^ell k, which is injective; it resets only the safe follower. Rule30's rank-one right ideal erases every later letter of its **spatial inverse observer**. Its lossless forward current in THM-4210 retains an infinite current sequence, the transverse channel and physical quadratic admissibility. Neither the41-state monoid nor the present scalar Mahler cost theorem reconstructs that physical carrier or solves a named-seed prize. Shared binary or dyadic notation is not an intertwiner.

## 8. Exact evidence and stopping scope

The standalone source imports no repository implementation. It uses declared reset prefixes100,00100,00100100 and0^(j-1)100 for2<=j<=6. On every complete unit phase cycle it compares the cost formula with literal ceiling orbits, native modular addresses, all-completion affine controls and actual terminal bit lengths. Small controls verify minimality by direct ordinary-integer search. Rational multiplication checks the first endpoint in (16), separately from the integer carry test.

The exact frequency formula and total-variation bound are checked without floating point. The sharp family is directly followed through its first rejection for1<=j<=7 and t2,3,4; its all-parameter claim is the proof above. Literal transformations check the Rule30 quotient hostiles, and the inherited A8/A13 and A1 corrections are retained. The source and evidence are prepared outside the repository; parent owns filing and promotion.

Both `python continuing8_20260906_mahler_reset_cost.py` and `python -O continuing8_20260906_mahler_reset_cost.py` pass **18,649 always-active exact gates**, with byte-identical raw LF output and JSON. Filed source paths resolve the certificate to the sibling results directory; the outside version writes it adjacent.

- Source SHA256:`77757ed311aac5246cf3a420dea282e84cb7aaaf428abcc928c633e8865f5646`.
- Output SHA256:`dfa7f2b47c883336c27aef04d55e9ac057e0a083c5a6a3841a68bc657d5c3e0c`.
- Certificate SHA256:`ac3f7f74a5ed80351218669c628a4bf068384f7e8132fa41af0b0828b6646d6c`.

The result is an exact reset cost and terminal-clock theorem, an infinite family with known first rejection, and a precise quotient obstruction. General nonzero reset skeletons, arbitrary continuations after the first post-runway1, the intersection of strict real safety with ordinary stabilization, and physical Rule30 asymptotics remain outside these conclusions. A useful next question is whether an analogous retained arithmetic cost can be controlled when a fixed nontrivial reset block is repeated, while preserving both its exact native address and its infinite-tail predicate.
