# Every finite automatic fixed-prefix selector is arithmetically sparse

Session197. **All-order obstruction to a specified finite-state inverse strategy.** A synchronous regular subrelation of any fixed physical final-gadget function has a polynomially sparse selected domain. Consequently a finite explicitly tagged bank of such automatic selectors covers a zero proportion of primes, and a zero proportion of seven-times-prime targets. This is not an omission theorem for tournaments or for nonautomatic inverses.

## 1. Exact statement and encoding

Fix a physical rooted prefix (h,e,R) of the reviewed one-seam grammar. Its final-gadget count is

    F(L)=h g(L)+R L−e(L mod2), L>=1,

where h,R are positive integers, 0<=e<=h and

    f(n):=g(n+1),
    f(0)=1, f(2x+1)=f(x), f(2x)=f(x)+x+(x mod2).       (1)

These physical laws are inherited from ARITHMETIC_INVERSE.md, SHA2564ac1727eb27c85f8a39f5598eb9f0e1b020bc4020cbca0a94e6288cb0676c3b7. We use ordinary canonical binary encodings, most-significant-bit first, padding the shorter track on the left by zeros. Extra common leading-zero columns and reversal of all tracks do not change regularity.

A selected subgraph is any set R0 contained in {(L,F(L)):L>=1}; it need not contain every exact pair. Assume its padded encodings form a regular language. Let D0 be its selected numerical domain.

**Theorem1.** There exist constants epsilon>0 and C, depending on the chosen regular subgraph, such that

    #{L in D0:L<=X} <= C X^(1−epsilon)               (2)

for all X>=1. The same type of bound holds for the number of its distinct selected target values at most X.

**Corollary2.** A finite explicitly tagged bank of regular selected subgraphs, one for each of finitely many fixed physical prefixes, cannot cover a cofinite set of odd targets. It cannot cover all sufficiently large prime targets or all sufficiently large seven-times-prime targets. In fact the fraction covered in either prime lane tends to zero.

This strengthens the separately reviewed complete-graph nonregularity theorem. It does not alter that published snapshot or assert sparsity for arbitrary nonregular selected domains.

## 2. A rectangular family in a complete live automaton component

Shift L=n+1 by regular binary addition. After projecting away an auxiliary track, first take the regular left quotient by common-zero columns and then intersect with canonical encodings. This normalization is necessary when an auxiliary number had the longest track. Thus R0 gives a regular selected subgraph of

    J(n)=F(n+1),

and its canonical binary domain D is regular. Removing n0 or other finitely many small inputs does not affect the desired bound.

Call a state live if it is reachable from the start and can reach an accepting state. A complete live strongly connected component is one whose two outgoing binary transitions from every state stay inside that component.

Suppose D has such a component. Choose a nonempty leading-one word p reaching it. This is possible in a canonical positive-integer language; if the chosen prefix is initially empty, append a1 within the component. A prefix beginning with0 could not have arbitrarily long accepted continuations in the canonical language.

Repeated1 transitions eventually cycle. Choose A>=1 at a cyclic state and its period t1. From that state, repeated0 transitions eventually cycle; choose B>=1 and period t0. Let T be a common positive multiple of t1,t0. A fixed suffix s then reaches an accepting state. Hence

    p 1^a 0^b s belongs to D
    whenever a=A+Ti, b=B+Tj, i,j>=0.                (3)

This uses only deterministic finite-state transitions, not a pumping assertion about arbitrary data values. The suffix may be empty.

## 3. Exact arithmetic on the rectangle

Let P be the integer encoded by p. Let t=|s| and S be the value of its t-bit word, allowing leading zeros; for an empty suffix set t=S=0. Let q_s be the t-bit reverse complement of s, interpreted as an integer, with q_s=0 when t=0.

Put x=value(p1^a0^b). Since a,b>=1,

    x=(P+1)2^(a+b)−2^b,
    f(x)=f(P)+(P+1)2^(a+b)−(P+1)2^a−2^b+2.         (4)

Indeed the ones leave f unchanged, and the zeros add the geometric prefix sum with a parity correction only at the first zero. In particular x is even.

Appending the fixed suffix gives

    f(2^t x+S)=f(x)+q_s x+f(S)−1.                  (5)

To verify (5) directly, every zero at position j from the left adds 2^j x from the old prefix. The only parity term depending on x could occur at the first suffix bit, and it vanishes because x is even. All remaining additions are exactly those obtained by reading the same suffix from zero, namely f(S)−1. This includes the empty suffix.

Define positive integers

    K=h q_s+R2^t,
    A0=(h+K)(P+1), B0=h(P+1), C0=h+K,

and the fixed integer

    D0=h(f(P)+f(S)+1)+R(S+1)−e((S+1) mod2).

For n=value(p1^a0^b s), equations(4)-(5) give exactly

    J(n)=A0 2^(a+b)−B0 2^a−C0 2^b+D0,             (6)

with B0,C0>0. All coefficients are fixed on the chosen rectangle.

The auxiliary number Y=A0 2^(a+b) is regular to recognize together with the input word. Its canonical binary representation is bin(A0) followed by a+b zeros, and

    |bin(Y)|−|input|=|bin(A0)|−|p|−|s|,

a fixed length offset. Fixed-prefix/zero-tail patterns and a fixed track-length offset are synchronous regular conditions.

Consequently, if the selected graph of J were regular, regular addition and existential projection would recover a regular relation between the input and

    Z=Y+D0−J(n)=B0 2^a+C0 2^b>0.                  (7)

All auxiliary-track projections use the common-zero left quotient before canonical intersection, just as in Section2. No unbounded numerical register or reversal of only one track is being smuggled into the automaton.

## 4. A bounded odd part forces a nonregular diagonal band

For a positive integer z, write odd(z)=z/2^v2(z). For any fixed bound H, the canonical language

    {bin(z): odd(z)<=H}

is regular: it is the finite union bin(d)0* over positive odd d<=H.

The odd part of B0 2^a+C0 2^b depends only on delta=a−b. Choose delta0=A−B, the difference on the infinite diagonal i=j in(3), and choose H to be the odd part at that difference. The set E of integer differences passing this odd-part bound is therefore nonempty.

It is also finite. Once delta+v2(B0)>v2(C0), the odd part of B0 2^delta+C0 grows exponentially with delta. The negative direction is symmetric. More explicitly, with ell=floor(log2 H), every allowed difference lies in

    −max(0,v2(B0)−v2(C0)+ell)
        <=delta<=max(0,v2(C0)−v2(B0)+ell).           (8)

If the shifted valuations are unequal, factor out the smaller power of two and compare its positive remaining term with H. If they are equal, their difference is already the one fixed alignment value, also within the displayed bound.

Intersecting the recovered regular relation with this regular odd-part language and projecting onto its input would make the following language regular:

    {p1^a0^b s: a=A+Ti,b=B+Tj,i,j>=0,a−b in E}.     (9)

It is not regular. Remove fixed p,s by regular left/right quotients. Choose M so that TM exceeds every |delta−delta0| for delta in E. Then the infinitely many prefixes1^(A+TMi) are pairwise distinguished by the suffixes0^(B+TMi): matching indices give difference delta0 in E, whereas distinct indices differ from delta0 by more than the allowed bound. This is an infinite Myhill–Nerode family.

The contradiction proves that no complete live strongly connected component can occur in D.

## 5. From automaton leakage to a strict power bound

Let the trimmed live state set have s0 states. If it is empty, D is empty. Otherwise every live state has a word escaping the live set. If not, the states reachable from it would form a closed finite directed subgraph and contain a complete live strongly connected component, contrary to Section4.

A shortest escaping word has length at most s0, since its visited live states do not repeat. Pad it to length s0. Once a reachable state is no longer coaccessible, no further suffix can restore acceptance. Thus, from every live state, at most2^s0−1 of the2^s0 binary blocks of length s0 remain live.

The number of accepted words of length n is therefore bounded by a constant times

    (2^s0−1)^(floor(n/s0)) 2^(n mod s0).

Its exponential base lambda=(2^s0−1)^(1/s0) is strictly below2. Summing lengths up to log2 X gives O(X^(1−epsilon)) for some epsilon>0; in the degenerate lambda1 case, absorb the resulting logarithm into any positive power below1. Translation n=L−1 preserves this counting conclusion. This proves the domain bound(2).

For the selected targets, physicality gives F(L)>=RL because h g(L)−e(L mod2)>=0. Hence an output at most X comes from L<=X/R. The number of distinct selected targets is at most the number of those selected inputs, proving the second assertion of Theorem1.

Fixing each finite prefix tag in a regular tagged relation gives one regular subgraph. A finite sum of the power bounds has a common exponent strictly below1. Thus Corollary2 applies equally to a finite explicitly tagged automatic witness relation.

## 6. Prime lanes and exact scope

An elementary estimate suffices for the prime conclusion. The central binomial coefficient satisfies

    binom(2n,n)>=4^n/(2n+1).

For each prime p, its valuation is at most floor(log_p(2n)): each term floor(2n/p^j)−2floor(n/p^j) is0 or1. Therefore

    binom(2n,n)<=(2n)^pi(2n),

so pi(X)>=cX/log X for a positive constant c at all sufficiently large X. A set of size O(X^(1−epsilon)) is consequently a zero fraction of the primes up to X. The same argument with X/7 treats seven-times-prime targets. Removing finitely many atom exceptions does not change this conclusion.

The theorem excludes finite synchronous regular selected-witness banks with an explicit recoverable prefix tag. It does not exclude:

- arithmetic inverses with unbounded registers or target-dependent splits;
- a bank of physical prefixes that grows with the target;
- a different nonlinear or nonregular witness-selection rule;
- an untagged union where an individual prefix cannot be recovered by a regular restriction;
- unrestricted tournament constructions or cofinite low-cycle coverage itself.

The modular-universality and complete-congruence results are consistent with this obstruction. Their witnesses are supplied by arithmetic constructions, not by a fixed finite tagged synchronous selector covering the target values exactly.

## 7. Reproduction and dependency boundary

The standalone verifier checks the exact shifted-run formula, empty and nonempty suffixes, the finite odd-part difference bound, concrete complete-component rectangle extraction and live-state leakage bounds. It also checks the elementary prime-count inequality. These finite controls corroborate the formulas and automaton lemmas; the all-order sparsity proof is the preceding argument.

Run `python3 verify.py` and `python3 -O verify.py`. Neither writes files. The only inherited tournament ingredient is the physical one-seam function and recurrence. Regular arithmetic, quotient normalization and the original balanced-language argument are also discussed in the separately reviewed fixed-prefix-nonregularity proof, published at ec65fdb76f9014f6397a588573adb2de03bb6716, SHA256fe0c822872555e54db26f313967fd97cdf0939737f39bd379c5a432dbaec8058. No finite coverage experiment is extrapolated here.
