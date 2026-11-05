# Arbitrary-modulus boundary banks and conditional probability laws

Session 193, 2026-11-05. **Elementary all-order proofs with standalone exact controls.** Independently reviewed. The global tournament Hamiltonian-path spectrum conjecture remains open. No literature-priority claim is made.

## 1. What this adds, and what is inherited

The main addition is an exterior-independent, equal-length residue bank for **every positive modulus**, followed by exponential mixing and size/residue independence in the original uniform strong one-apex ensemble. This combines an exact continuation identity, finite modular arithmetic, and the existing two-boundary probability law on the same physical words.

S192 already proves simultaneous square-root nets modulo every power of two, and its follow-up improves dyadic bank lengths to `k^2+2k-2`. We do not claim the dyadic net as new, nor an improvement to that length. The bank here handles arbitrary odd factors, where the dyadic finite-range insulation no longer applies. It has a deliberately elementary and usually much worse length bound. S25 already proves modular completeness for each modulus without fixed exteriors or a common length. We use its block recurrence and supply the two missing synchronization properties.

A second result makes the topological support of the S192 fixed-apex boundary law explicit: for at least two apexes even the strong-conditioned limiting law has the full interval `[0,r!]` as support. This yields a deterministic strong two-apex macroscopic net without assuming a density. Neither result proves integer intervals. The inherited one-apex positive missing-density theorem remains a decisive control.

Sources:

- [Session 25 spectrum proof](../../session25_spectrum_wiring_20261021/session25-spectrum/PROOF.md), sections 2–4: binary recurrence and modular block counter.
- [Session 191 one-apex proof](../../session191_tournament_growth_20261104/one_apex/PROOF.md): exact continuation, strongness and one-apex H=13 exclusion.
- [Session 192 dyadic proof](../../session192_surprising_approaches_20261105/dyadic/PROOF.md), sections 1–5: two-ended selector and inherited S26 product law.
- [Session 192 dyadic compression](../../session192_surprising_approaches_20261105/dyadic_compression/PROOF.md): sharper dyadic banks and their information/error barrier.
- [Session 192 fixed-apex limits](../../session192_surprising_approaches_20261105/boundary_limits/PROOF.md), sections 1–4: fixed-apex boundary statistic, uniform coupling and strong-conditioned product law.

## 2. Exact physical continuation coordinates

For a binary word w, T(w) consists of an increasing chain and one apex x, with `x -> i` iff the i-th bit is 1. Put H(empty)=1 and define

    A(w) = binary integer represented by w,
    C(w) = sum_(j=1)^|w| (1-w_j) 2^(j-1),
    p(w) = last bit, or 0 for empty,
    epsilon(w) = [w is nonempty and begins in 0].

The inherited recurrence, which follows by classifying paths around the apex, is

    (H,A,p) --1--> (H,2A+1,1),
    (H,A,p) --0--> (H+A+p,2A,0).

Its exact concatenation consequence is

    H(uv)=H(u)+H(v)-1+A(u)C(v)+p(u)epsilon(v).        (1)

For nonempty v of length L, apply (1) twice:

    H(uvw)=H(u)+H(v)+H(w)-2
           +A(u)C(v)+(2^L A(u)+A(v))C(w)
           +p(u)epsilon(v)+p(v)epsilon(w).           (2)

This is an equality of actual Hamiltonian-path counts, retaining both boundary integers. H is always odd: appending 1 preserves H, while A+p is even. A nonempty word beginning in 1 and ending in 0 gives a strong tournament.

## 3. An equal-length bank for an arbitrary modulus

Fix M>=1. Write `M=2^k m`, m odd, and put

    g=gcd(2,M),  T=M/g.

Choose an integer ell>=max(1,k) with `2^ell=1 mod m`. For m=1 this congruence is vacuous. Let

    Q=2^ell, D=(Q-1)^2+1, B=1^ell 0^ell,
    L=2 ell T,
    v_t=B^t 0^(2 ell(T-t)),  1<=t<=T.                (3)

All T words have exactly length L, begin in 1, end in 0, and are distinct.

**Theorem 1 (universal arbitrary-modulus bank).** For every t,

    A(v_t)=C(v_t)=0 mod M,
    H(v_t)=1+tD mod M,
    gcd(D,M)=g.                                     (4)

Consequently, for every pair of exterior words u,w, including empty ones,

    H(u v_t w)=K(u,w)+tD mod M,
    K(u,w)=H(u)+H(w)-1+2^L A(u)C(w).                 (5)

As t ranges from 1 to T, these residues are exactly all residues if M is odd, and exactly all odd residues if M is even, each once.

**Proof.** For a state ending in 0, appending B has the exact S25 action

    A'=Q^2 A+Q(Q-1),
    H'=H+Q(Q-1)A+D.                                 (6)

The integer `Q(Q-1)` is divisible by M: Q contains its power-of-two part, and Q-1 its odd part. Iteration from the empty state yields `A(B^t)=0 mod M`, `H(B^t)=1+tD mod M`. Also

    C(B^t)=Q(Q-1)(1+Q^2+...+Q^(2t-2))=0 mod M.

The zero padding has length d=2 ell(T-t). It multiplies A by `2^d`, changes H by `A(B^t)(2^d-1)` because the old last bit is 0, and changes C by `Q^(2t)(2^d-1)`. The last expression is divisible by `2^k` since t>=1, and by m since d is a multiple of ell. This proves the first two lines of (4), including t=T when d=0.

Modulo m, D=1. Modulo `2^k`, D=2 (when k>=1). Hence its gcd with M is g. In (2), v_t has epsilon=0, p=0, and both boundary integers vanish modulo M. Substituting (4) gives (5). Its progression has period exactly T. For even M, K is odd because H(u),H(w) are odd and `2^L` is even. QED.

The t>=1 condition is essential. The tempting t=0 all-zero padding has `C=2^L-1`, which is not zero modulo an even M. The bank includes t=T for the zero increment instead.

**Length and decoding.** If m>1, take the least positive multiple of `ord_m(2)` that is at least max(1,k); for m=1 take max(1,k). Thus

    ell < max(1,k)+ord_m(2)   (m>1),
    L=2 ell M/g = O(M(k+m)).                           (7)

For odd M this is O(M^2); for a power of two it is `M log_2 M` when M>=2. This is not competitive with the inherited quadratic-in-log dyadic bank. The exact quantities ell and L, rather than an unspecified modulus-dependent constant, enter every estimate below.

Given an admissible desired residue r, compute K modulo M. Solve

    t = ((r-K)/g) (D/g)^(-1) mod T,                    (8)

using integer representatives with r-K divisible by g, and choose t in 1..T. For T=1 choose t=1. No search over tournaments is needed.

## 4. Magnitude and all finite congruences on the same tournament

We recall the proved two-ended selector, not as a new result. For any middle length L, a,b>=1, n=a+b+L and real h in `[2^(n-2),2^n]`, set `z=sqrt(h/2^n)`. Choose prefix u to be the a-bit expansion of

    min(2^a-1,floor(2^a z)),

and suffix w to be the reverse complement of the b-bit expansion of

    min(2^b-1,floor(2^b z)).

For every middle word v of length L, the word uvw is strong and

    |H(uvw)-h| <= E(n;a,b),
    E(n;a,b)=2^(n-a)+2^(n-b)+floor(n/2)+1.             (9)

For completeness, define X as the left binary fraction and Y as the right reverse-complement fraction of a strong word. The exact slot formula gives

    H=2^n XY+2+N_01-R,
    0<=R<=N_01,  N_01<=floor((n-2)/2).

The selector makes `|X-z|<=2^-a`, `|Y-z|<=2^-b` independently of the middle. Since all quantities lie in [0,1], (9) follows. The nonnegative correction is at most `floor(n/2)+1`.

**Corollary 2 (arbitrary fixed-modulus square-root net).** Use L from (3), n>=L+2, balanced a,b with a+b=n-L, any admissible residue r modulo M, and h in the above interval. There is an explicitly decoded strong one-apex word of chain length n satisfying

    H=r mod M,
    |H-h| <= 2^ceil((n+L)/2)+2^floor((n+L)/2)
                 +floor(n/2)+1.                     (10)

Choose u,w by the selector, then t by (8). Both assertions hold on the **same** physical word u v_t w. In particular, every fixed finite collection of congruences jointly compatible with each other and oddness can be imposed by taking their combined modulus, with a square-root-scale error as the target grows. For fixed M the error is `O_M(2^(n/2))`; the constant grows as `2^(L/2)`. This cannot be used with a growing modulus while hiding L.

The familiar exactness criterion would be `H=h mod M` and `|H-h|<M`. Here L>=2T and the displayed error is at least `2^(L+2)>M`, even at n=L+2. Thus this unchanged separate error estimate cannot force equality. This does not prove a lower bound on actual errors or preclude a sharper, correlated correction.

## 5. Exact modular mixing in the original random-word ensemble

Let W_n be uniform among all strong binary presentations of chain length n>=2: its first bit is 1, last bit 0, and all interior bits are independent fair bits. Let R_M be uniform on the T parity-admissible residues.

**Theorem 3 (all fixed moduli, explicit exponential rate).** With p=T/2^L,

    d_TV(law(H(W_n) mod M),law(R_M))
       <= (1-p)^floor((n-2)/L).                       (11)

**Proof.** Partition as many interior positions as possible into q disjoint L-bit blocks. Each block belongs to the T-member bank with probability p, independently. On the event E that at least one does, condition on the index of the first such block, on all other bits, and on the fact that it belongs to the bank. Its member is then uniform among the T distinct bank words. The universal-exterior property (5) maps these members bijectively to the admissible residues. Hence the residue is exactly uniform on E. The exceptional event has probability `(1-p)^q`, yielding the mixture bound (11). QED.

This is a theorem about the original uniform word ensemble, not a reweighted constructed subfamily. No independence of all bank occurrences is asserted after conditioning on the final residue. The first-occurrence conditioning avoids that mistake. For large M, the explicit mixing rate is extremely weak; no uniform-in-M mixing assertion is hidden here.

## 6. Every exact residue fibre has the same macroscopic law

Keep disjoint boundary windows of lengths a,b>=1, a+b<=n, and put q=floor((n-a-b)/L). Use bank tests only in the remaining middle positions. Define

    delta=2^-a+2^-b+(floor(n/2)+1)/2^n,
    eta=(1-T/2^L)^q.

**Theorem 4 (joint size/residue coupling).** There are independent `U,V ~ Uniform[1/2,1]` and an independent uniform admissible residue R* such that

    |H(W_n)/2^n - UV| <= delta,
    P(H(W_n) mod M != R*) <= eta.                    (12)

**Proof.** Read the left a bits and the right b reverse-complement bits. They are independent dyadic uniform grids with initial bit 1. Add independent uniform tails in their cells to obtain U,V. The product identity used in (9) gives the first inequality for every middle word. Let E be the first-bank event in the middle. On E take R* to be the actual residue; on E^c draw it independently uniform. Conditional on the boundary windows, E has the same probability and the actual residue on E is uniform by the first-block argument. Therefore R* is uniform independently of the windows and their added tails. Its mismatch probability is at most P(E^c)=eta. QED.

For n tending to infinity with, for example, a=b=floor(n/4), delta and eta tend to zero. Consequently

    (H(W_n)/2^n, H(W_n) mod M) ==> (UV,R_M),           (13)

with independent components. The limit density of UV is the inherited one

    f(z)=4 log(4z), 1/4<z<1/2;
         -4 log z, 1/2<=z<1;
         0 otherwise.

It is bounded by D0=4 log 2. If F is its CDF, for every admissible r and every real z, (12) gives

    F(z-delta)/T-eta <= P(H/2^n<=z, H=r mod M)
                         <= F(z+delta)/T+eta.        (14)

In particular, if T eta<1,

    sup_z |P(H/2^n<=z | H=r mod M)-F(z)|
        <= (D0 delta+2T eta)/(1-T eta).              (15)

This proves a conditioned probability law, rather than assuming that an unconditioned continuous distribution remains valid on exact arithmetic fibres.

**Decisive scope control.** S26's distinct-support bound, restated in S192 arithmetic, is

    #{one-apex attained H<=X} <= (log 2)X/2+O(sqrt X).

Thus even the conjunction of (10), uniform residue mixing, and the same smooth limit law in every fixed arithmetic fibre coexists with a positive proportion of missing odd integers. H=13 is the smallest relevant explicit class hostile. These results solve the macroscopic and fixed-congruence portions of a possible proof; exact unit-scale local support needs another operation or estimate. They do not claim that two-apex support has positive missing density.

## 7. The fixed-apex boundary law has full support when r>=2

Use S192's boundary statistic U_r for a word over r-bit letters. Slots 0,...,r are assigned uniformly to each position. U_r is the probability that the first occurrence of each constrained slot j>=1 has bit j equal to 1. A prefix of length k determines U_r to error at most `r alpha^k`, where `alpha=r/(r+1)`.

Fix all bits in coordinates 1,...,r-1 identically equal to 1, and allow only the last bit b_i to vary. Then every other first-boundary condition succeeds automatically, so

    U_r = (1/(r+1)) sum_(i>=1) alpha^(i-1)b_i.        (16)

The two maps `x -> alpha x` and `x -> 1/(r+1)+alpha x` send [0,1] to two intervals whose union is [0,1]. Choosing an inverse branch recursively gives a binary sequence representing every x in [0,1]. If r>=2, the first letter of every such word is nonzero because at least one permanently-1 coordinate exists.

**Theorem 5.** For r>=2, both U_r and its version U_(r,+) conditioned on a nonzero first letter have topological support exactly [0,1]. Consequently S192's unconditional and strong-conditioned normalized path-count limits both have support exactly `[0,r!]`.

**Proof.** Equation (16) realizes every x deterministically. The uniform prefix bound makes U_r continuous in the infinite word. Every finite prefix has positive probability under the iid uniform-letter model, and the constructed prefixes remain admissible under the strong first-letter conditioning. Thus every neighborhood of every x has positive probability. U_r lies in [0,1], so this is its entire support. Products of two independent such variables have support [0,1]; multiply by r! and apply S192's proved limit law. QED.

For r=1 the strong-conditioned boundary support is instead [1/2,1]. Dropping this exception would contradict the one-apex size window. No absolute continuity of U_r or its product is used or claimed.

### A constructive strong two-apex net

For r=2, take an arbitrary target z in [0,2]. Generate a binary expansion (16) for x=z/2, and use its first k letters, with the first apex bit fixed to 1, as prefix u. Set the suffix to `0^k`, and put a letter 2 somewhere in the intervening nonempty middle. For any m>=2k+1 this gives a strong two-apex word: its first letter is nonzero, last letter is not 3, and it contains 2.

S192's uniform boundary theorem and its bound `0<=U-G(prefix_k)<=3(2/3)^k` give

    |H(T(w))/3^m-z| <= 12(2/3)^k.                    (17)

Indeed the nonnegative transfer error is at most `12(2/3)^k`, while the deficit of its main term `2G(u)G(3^k)` from z is also between zero and this bound. The difference of those two nonnegative quantities has absolute value at most their common upper bound, not twice it. Choosing k=floor((m-1)/2) gives an explicit strong macroscopic net at the `3^m` scale. It does not give +2 steps, residue control in this two-apex family, or exact integer coverage.

## 8. Reproduction

Run `python3 verify.py` and `python3 -O verify.py`; stdout and receipt bytes must agree. The standalone standard-library verifier uses explicit exception gates. It checks graph-level path counts, concatenation, arbitrary-modulus banks and their conditioning interfaces, simultaneous selectors, complete small-word residue histograms, and the two-apex net against an independently implemented exact run-slot transfer and selected graph-level subset DPs. All infinite conclusions above have proofs; finite receipts are controls, not extrapolations.
