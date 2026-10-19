# Realized first-crossing sources cannot contain macroscopic expansion repetition

19 October 2026. **PROVED, with a CITED explicit logarithmic-form input.**
This is an unbounded-episode, all-word-length theorem in the specified bank
language. Arbitrary bank entry, general convergence, and complete description
of the realized 2-adic boundary remain **OPEN**.

## Main theorem

Use the chronological Terras bank

    K=110110110110110101011011101010,
    F=10, H=110101000, V=K^4 F.

K,V are expansion letters; F,H are explicit reset letters. Throughout this
note, a first-crossing bank word starts with a nonempty expansion run over
{K,V}. Stop a fully lawful such word W at its first explicit reset boundary
at which the cumulative affine multiplier is below one. This expansion-first
restriction matters for the standalone carry and short-range lemmas: F at
x=1 is a return, and is not in the language considered here. Let N be its Terras length, and let L count
all Terras bits belonging to its expansion letters, including internal bits
of V. There is no bound on the number of expansion/reset episodes.

Suppose W contains a consecutive factor U^k, where U is a nonempty word over
{K,V}, k>=2. Define the *duplicated length* r=(k-1)|U|. If

    100r >= L,                                      (1)

then T^N(x)<x for every positive integer x following W.

In particular it suffices to have 100r>=N. The theorem tests a repeated
factor anywhere in the actual macroblock word, not necessarily at its start.
It proves descent at the endpoint itself, with no extra repair time and no
periodic alternative. A repeated contracting whole word is not assumed.

This genuinely extends beyond bounded-episode repair and beyond the
single-episode zero-head theorem. An explicit family with arbitrarily many
expansion episodes is given in Section 7. It supplies no arbitrary-integer
entry theorem.

## 1. Exact inherited arithmetic and the source floor

A parity word of n bits and h odd bits has branch (a x+b)/M, with a=3^h,
M=2^n and b>=0. Its unique lawful residue is -b a^(-1) modulo M. The
inherited binary-fuel identity says

    U^k lawful at y iff M^k divides (a-M)y+b.        (2)

For expansion words U over {K,V},

    a>M, b/(a-M)<=C,
    C=83829391903560839836419813590128166431
      /157489106280555890461927692448373459 <533.

These are the scoped identities in
[the bank proof](../collatz_bank_reset_20261019.md) and
[the pumping-memory note, Section 5](../collatz_guarded_pumping_memory_20261004.md).

The inherited forced-prefix argument is independent of episode count. A
non-descending first-crossing endpoint must have source in the K^4 cylinder,
so

    x>=r4=588943162011782273059357097228697595.       (3)

Indeed an initial expansion run not beginning K^4 is K, K^2 or K^3. Its
first reset already contracts and descends on every lawful positive source.
Any V contains K^4 at its start. The present checker reproduces the exact
short-word inequalities and K^4 residue. See
[the interleaved proof, Section 4](../collatz_interleaved_height_repair_20261019/INTERLEAVED_HEIGHT_REPAIR.md).

## 2. A step-aware carry bound at every first crossing

Write a macroblock as s_j z+e_j and let S_j be the cumulative multiplier
through that block, with S_0=1. A direct exact check on the four letters gives

    e_j/s_j = b_j/a_j <= n_j/6,                    (4)

where n_j is its Terras length. Equality occurs for F; the other three
inequalities are strict. Explicitly the quantities b_j/a_j are

    K: 4912087879/1162261467 <30/6;
    V: 83829391903560839836419813590128166431
       /5474401089420219382077155933569751763 <122/6;
    F: 1/3=2/6;
    H: 101/81<9/6.

If a prefix P has affine map S_P x+e_P, the composition identity is

    e_P/S_P = sum over j in P of (e_j/s_j)/S_(j-1). (5)

Every macroblock input before the first crossing has cumulative multiplier
at least one: expansion letters increase it, and an earlier explicit reset
below one would already have stopped W. Equations (4)--(5) therefore imply

    e_P/S_P <= |P|/6.                              (6)

For the full first-crossing word, writing its map as Sx+e with S<1 gives

    e/S < N/6, hence e<N/6.                        (7)

Strictness follows because W includes an expansion letter. This is a new
step-aware bound, compatible with the inherited e<m(C+1) episode bound:
one may use the smaller of the two. It is not a statement about arbitrary
contracting words whose earlier coefficient multiplier has already crossed.

## 3. Exact short-range closure with arbitrarily many episodes

### Lemma: every expansion-first first-crossing word with N<=2^52 descends

Let H0=2^52, alpha=log2/log3, and use the rational certificate

    p=5750934602875680, q=9115015689657667,
    gcd(p,q)=1, q>2H0, |alpha-p/q|<1/q^2.           (8)

The last inequality is rigorously certified using rational lower and upper
atanh-series bounds for log2 and log3, not a floating-point continued
fraction claim. For 1<=N<=H0 and integer h, the integer Np-hq is nonzero,
since otherwise q would divide N. Hence

    |N alpha-h| >=1/q-N|alpha-p/q| >1/(2q).         (9)

Suppose a first-crossing endpoint is non-descending. By (3), x>=r4. Put
u=1-3^h/2^N>0. Equations (7) and (3) give

    u < H0/(6r4) <1/2.

Thus A=N log2-h log3=-log(1-u)<2u<H0/(3r4). But A>1/(2q) by (9) and
log3>1. The exact integer inequality

    2H0 q < 3r4                                  (10)

makes these bounds incompatible. This proves the lemma.

The checker also verifies the much stronger margin
2136 H0 q<r4, which would suffice with the older coarse carry bound e<534N.
The present short-range result therefore does not reuse the single-episode
short-length theorem outside its scope. It is proved here directly for
unrestricted interleaving.

## 4. Explicit all-length polynomial height ceiling

Set B=10^10. **CITED:** the positive-rational specialization of Matveev's
theorem gives, for 1<=h<=N and 3^h<2^N,

    1-3^h/2^N > exp[-B(1+log N)].                  (11)

The version checked directly is Theorem 2.1 of Languasco, Luca, Moree and
Togbe, [Sequences of integers generated by two fixed primes](https://link.springer.com/article/10.1007/s12188-025-00293-9#Sec2),
which credits Matveev and Bugeaud--Mignotte--Siksek. Apply it to 3^h 2^(-N)-1;
the bases are positive rationals, exponent bound is N, and unique prime
factorization makes the expression nonzero. Its coefficient is less than
2*30^5*32*2<B, using log2<1 and log3<2. This is a cited theorem application,
not a numerical verification of a logarithmic-form theorem.

A non-descending first-crossing source satisfies (1-S)x<=e<N/6. Therefore

    x < (1/6) exp(B) N^(B+1).                     (12)

Because x>=r4>C, the augmented source height

    Z=x+N/6+C

obeys Z<2x+N/6<(1/2)exp(B)N^(B+1). We shall only need

    log Z < B+(B+1)log N.                         (13)

## 5. Transport a repeated factor back to the original source

Write W=P U^k Q_tail and n=|U|, r=(k-1)n. P may contain arbitrarily many
previous episodes. Let y be the actual value at U's start and S_P the
cumulative multiplier there. Since this is before the first crossing,
S_P>=1. Equation (6) gives

    y+C <= S_P(x+|P|/6+C) <= S_P Z.               (14)

Write U(z)=(a z+b)/M, M=2^n. The positive integer in (2) is at least M^k.
Together with b/(a-M)<=C this implies

    y+C > M^k/a = 2^r/s_U.                        (15)

Combining (14)--(15),

    Z > 2^r/(S_P s_U).                            (16)

This is the source-preserving context cancellation used here. Reset letters
cannot be deleted from a guarded word; their effect is instead retained in
S_P and the carry debt |P|/6.

Put lambda=3^19/2^30. K has slope lambda and length 30, while V has slope
lambda^4(3/4) and length 122. Consequently every expansion subword of total
length t has multiplier at most lambda^(t/30). Resets have multiplier below
one, and the product of all expansion multipliers dominates S_P s_U. Hence

    S_P s_U <= lambda^(L/30).                     (17)

Exact arithmetic gives 1<lambda<13/12, so log lambda<1/12. Using log2>1/2,
condition (1) and (16)--(17) give the key exponential lower bound

    log Z > r log2-(L/30)log lambda
          > L/200-L/360 = L/450.                 (18)

No ordering of a 2-adic number is invoked. Equations (14)--(18) concern
actual positive integer orbit values and their lawful finite words.

## 6. Join the two ranges: proof of the main theorem

The first crossing also bounds N in terms of expansion length L, despite
unrestricted interleaving. Let t be the number of explicit reset letters.
Just before the last reset the cumulative multiplier is at least one. The
product of all expansion multipliers is at most (3/2)^L, and each reset
multiplier is at most 3/4. Therefore

    (4/3)^(t-1) <= (3/2)^L,
    t-1<3L/2,
    N<=L+9t<29L/2+9<=15L,                        (19)

where (3/2)^2<(4/3)^3 and L>=30 are sufficient for the last comparison.
The product calculation does not reorder any orbit or change a parity guard.

Suppose the main theorem fails. Set L0=2^48.

- If L<L0, (19) gives N<15L0<2^52. Section 3 excludes non-descent.
- If L>=L0, (13), (18), and (19) would imply

      L/450 < B+(B+1)log(15L).                    (20)

  At L=L0, because 15L0<2^52 and log2<1, the right side is less than
  53B+52. But the exact inequality L0/450>53B+52 holds. The derivative of
  the left side minus the right side is 1/450-(B+1)/L>0 for all L>=L0,
  because L0>450(B+1). Thus (20) is impossible throughout this range.

The two contradictions prove strict endpoint descent for every allowed
length, episode count, and lawful source. The constants 100, 450 and the
cutoff are sufficient choices, not claims of optimality.

## 7. Nonvacuity with unbounded expansion episodes

For L>=2 choose j>=1 least such that

    lambda^(L+j-1) (3/4)^j <1.

Such j exists because 3lambda/4<1. Use the actual bank word

    W_L = K^L (F K)^(j-1) F.                      (21)

Its earlier explicit reset boundaries all have multiplier greater than one,
by minimality of j and prime-factor uniqueness. Its final boundary is the
first crossing. Also lambda^2(3/4)<1, so j<=L; and j tends to infinity with
L, since any fixed j eventually fails the displayed inequality.

The repeated factor K^L contributes duplicated length 30(L-1), whereas the
total expansion length is 30(L+j-1)<=30(2L-1). Thus 100r>=L_exp with ample
margin. The theorem proves endpoint descent on the full lawful cylinder of
every W_L. Every cylinder has infinitely many positive integer lifts.

The parse has j expansion episodes and j is unbounded. The exact rewrite
K^4F=V can remove the first explicit reset when available, but every later
reset is preceded by only one K, leaving unboundedly many episodes. Thus
this family is not explained by raising a fixed episode cap or by the
one-episode repair theorem.

## 8. Consequences for S_d and its 2-adic boundary

For d>0, retain only sources satisfying the exact first-crossing endpoint
equation

    T^N(x)=x+d,
    x=(b-2^N d)/(2^N-3^h) in the positive integers. (22)

The theorem says that every witness for membership in S_d must avoid every
factor U^k in expansion letters with 100(k-1)|U|>=L_exp. This is an actual
endpoint-realizability exclusion, not a condition on all formal bank paths.
It also excludes d=0 on the same repeated-factor language.

Before making the coarse 1% estimate, (13), (16)--(17) give the more precise
necessary repetition budget for any non-descending first-crossing witness:

    r log2 < (L_exp/30)log lambda + B+(B+1)log N.   (23)

With a specified context P, L_exp in the first term may be replaced by the
expansion-bit length of P plus |U|. This keeps the actual position of the
repeated block, rather than charging unrelated later expansion to it.

In particular, initial K^L approximations to the formal source -c still
obey the sharper anchor constraint from the companion note,

    15L < B+38+(B+1)log N.                         (24)

So an endpoint-realizing sequence approaching K-infinity must use
exponentially long completions relative to its initial K depth. The main
theorem additionally forbids those completions from containing any
macroscopic expansion repetition as in (1).

Neither fact alone proves -c absent from Lambda_infinity(S_d): extremely
long and repetition-poor interleavings are still possible formal candidates.
A finite complete symbolic forest cannot drop this arithmetic distinction.
The refined open target is to rule out endpoint solutions (22) along those
long, repetition-poor completions, or to prove translated coefficient
stopping on exactly their realized limit set. Arbitrary pointwise bank entry
remains another, separate problem.

## 9. Connections recovered, with the lost information named

- **Pumping memory -> actual source exclusion.** Repetition's divisibility
  produces a positive integer fuel quotient, hence (15). Its real consequence
  survives an arbitrary earlier interleaved context only after retaining the
  normalized carry debt (6). A scalar multiplier alone loses this fact.
- **Interval/carry certificates -> uniform all-length theorem.** One exact
  rational bracket covers the full short interval. A cited explicit
  logarithmic-form bound covers its complement. Their verified overlap,
  not a long list of examples, supplies the universal quantifier.
- **Rational recurrence and p-adic trees -> source-sensitive boundary.** The
  repeated anchor names a formal branch; (22), context length, and fuel are
  the necessary sidecars. The result filters actual endpoint witnesses
  without claiming a universal rational coefficient-stopping theorem.
- **Ordered cancellation and tournament insertion.** The prior
  [context-cancellation note](../collatz_context_cancellation_20261007d.md)
  requires a matched ordered middle word, and the
  [tournament insertion note](../tournament_insertion_slots_20260926.md)
  refutes a recursion that deletes a vertex without preserving the path.
  Here (14) transports through, rather than deletes, an earlier context.
  There is no tournament invariant being promoted to a Collatz theorem.

The companion [power-expansion note](POWER_EXPANSION_ZERO_HEAD.md) records
the discovery route and its simpler single-episode corollary. The present
note is the stronger result relevant to the unrestricted episode barrier.

## 10. Reproduction and audit boundary

Run `python3 verify_interleaved_repetition.py` and repeat with `python3 -O`.
The outputs must match the receipt. The exact controls separately authenticate
all block identities, the rational logarithm bracket and overlap constants,
complete bounded prefix trees, source/endpoint guards, the transported fuel
inequality, and the unbounded-episode family at a declared finite control set.
Positive, zero-return, wrong-guard and earlier-crossing hostiles are retained.
The written proof supplies all-length quantifiers; the external logarithmic
bound remains a **CITED** input. The full S_d boundary is **OPEN**.
