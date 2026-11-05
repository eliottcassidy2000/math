# Exact two-sided dyadic returns and their lacunarity barrier

Session 194, 2026-11-05. **PROVED elementary all-order statements from explicit finite seeds; independent audit PASS.** The unrestricted Hamilton-path spectrum conjecture with exceptions `{7,21}` remains **OPEN**. This packet gives a genuine core-changing escape from the S193 one-sided response obstruction, an exact classification of that escape, and a proof that it still cannot amplify a finite seed bank into unbounded consecutive dyadic-offset windows.

## 1. Inheritance and scope

Write H(V) for directed Hamiltonian vertex words, with H(empty)=1. The closest mechanism is [Session 193 verified guide](https://github.com/eliottcassidy2000/math/blob/169792fb5eeb483e351112b5a09706e6d7bd89f0/00-navigation/S193.md), exact cut-response packet: a nonempty attachment S on a fixed core has response R(S), and sink padding changes it by `R -> 2R+bH(V)`. Pure-dyadic response therefore permits only zero sink digits. Its finite bank proves 135 consecutive offsets around every sufficiently large power of two.

The apparent escape tested here changes the core by adding **sources as well as sinks**. The source law is different: `R -> R+bB`, where B counts ordered two-path covers. New pure-dyadic responses genuinely become possible. The hostile example to the stronger, false assertion “all transitive padding preserves the old dyadic response exponent” is the explicit seven-vertex core in section 5.

The corrected near miss is to read infinitely many new slopes as growing exact windows. The underused sidecar is the *opposite* two-path response L(S), which measures the simultaneous loss in the transported intercept. Its positivity forces lacunarity.

Concept board: physical cut responses; sources versus sinks; exact modular return; opposite-end boundary counts; growing windows versus isolated new offsets. Anchor: exact coverage under a wider operation family. Niche: the physical source-response law. Wildcard: exponential-plus-linear tournament families. No literature-priority claim is made.

The results apply to a fixed finite bank of **nonempty marked attachments**. They do not exclude arbitrary new cores, old-arc changes, additional independently variable apexes, or constructions whose response need not be a power of two.

## 2. The two-sided physical interface

Fix a nonempty tournament V and nonempty S subset V. The distinguished apex x beats exactly S. Set

    h = H(V),                 F = H(V+x),
    B = sum_(A subset V) H(V[A]) H(V\A),
    R = sum_(A subset V) sum_(i in A intersect S) s_i(V[A]) H(V\A),
    L = sum_(A subset V) sum_(i in A intersect S) e_i(V[A]) H(V\A),
    s = sum_(i in S) s_i(V),   e = sum_(i in S) e_i(V),
    b = F-R+e.

Here s_i and e_i count Hamilton paths beginning and ending at i. Thus B counts two ordered paths with empty slots allowed. The first/second choice of slot does not change B. Splitting every Hamilton path just before, or just after, a named vertex proves

    R >= h,   L >= h,   0 <= s,e <= h,
    R <= B-h < B,      L <= B-h.                       (1)

The totals over all vertices are B-h. Rédei parity gives h odd. Complementary subsets pair off, so B is even.

### Source law

Adjoin a new base source y, with y beating every old base vertex. Let a in {0,1} say whether x beats y. Then

    h' = h,         B' = 2B,
    R' = R+aB,      L' = 2L+ah,
    s' = ah,        e' = e,
    F' = F+a(B-L-s),
    b' = b-a(L+s).                                    (2)

**Proof.** A spanning two-path cover puts the source at the beginning of exactly one of its two slots, giving B'=2B. If the second path begins at an old vertex i, y must lie in the first path: old R_i is unchanged. If it begins at y, the rest is an arbitrary two-slot cover of V: R_y'=B. Dually, y can be put at the beginning of either slot of a cover whose first path ends at old i, so L_i'=2L_i; the first path ending at y must be the singleton, giving L_y'=h. Every Hamilton path of the base starts at y, while its end is unchanged. This proves the boundary transformations.

If a=0, y is a source of the complete tournament too, so F'=F. If a=1, the directed-cut capacity between y and i is L_i+s_i, by the source dual of the S192 sink-capacity formula. Only old i outside S contribute newly; their total is `(B-h-L)+(h-s)=B-L-s`. This proves F', and subtraction proves b'. Alternatively the same identity follows by reversing all arcs in the sink law. QED.

For comparison, a new base sink with apex digit a has

    R' = 2R+ah,     e' = ah,
    F' = F+(1-a)(R+e),       B'=2B.                    (3)

All of these are actual vertex additions, not an abstract sum of cut supports.

### Arbitrary mixtures can be put in normal form

Allow any finite sequence of source and sink additions, always relative to the current base and with arbitrary incidences to x. Sources are in reverse creation order before V, and sinks are in creation order after V. Their mutual arcs do not depend on how the two types were interleaved. Thus the same final labelled tournament is obtained by doing all source additions first, then all sink additions, preserving each type's internal order.

Write the source digits in creation order as a_1,...,a_l, and put

    J = sum_(j=1)^l a_j 2^(j-1).

Let t be the number of sink additions and A the ordinary binary integer read from their digits in creation order. Then

    R_final = 2^t(R+BJ)+hA,       0 <= A < 2^t.        (4)

Zeros at the end of the source word add neutral source vertices without changing R, F, e or b. Consequently J, rather than the padded word length, determines the intercept.

## 3. Exact classification of every dyadic return

**Theorem 1 (two-sided return criterion).** Under the hypotheses above, `R_final=2^k` if and only if every sink digit is zero and

    R+BJ = 2^q,   J>=0 an integer,   k=q+t.            (5)

For each admissible q the source word is uniquely determined by the binary digits of J, apart from irrelevant trailing zeros.

**Proof.** Since R+BJ>=h>=1, (4) implies k>=t. Modulo 2^t, odd h is invertible, so hA=0 forces A=0. Equation (5) follows, and conversely suffices. Uniqueness is ordinary binary expansion. QED.

This is already an exact, target-independent classification: one checks `2^q = R mod B` and the nonnegative quotient J. The following clean form applies to every S193 seed.

**Corollary 2 (a dyadic seed).** Suppose R=2^r, write B=2^v m with m odd, and recall R<B.

* If r<v, or if m=1, the only return is J=0, q=r.
* Otherwise m>1, and with `d=ord_m(2)` all returns are exactly

      q=r+nd,   J=(2^(r+nd)-2^r)/B,   n=0,1,2,... .  (6)

  Here d>=2. Thus there are genuinely new dyadic slopes at every positive n.

**Proof.** Nonnegative J implies q>=r. If r<v and q>r, the 2-adic valuation of `2^q-2^r` is r<v, impossible. If m=1, R<B implies r<v. In the remaining case divide the congruence by 2^r modulo m; powers of two are units modulo odd m, and the condition is exactly d|(q-r). The power-of-two part divides the difference because r>=v. Finally d=1 would force m|1, whereas m>1. QED.

More generally, two different admissible exponents can never differ by one: if B divides both `2^q-R` and `2^(q+1)-R`, it divides R, contradicting 0<R<B.

## 4. Why the new returns still cannot produce growing windows

For the source word associated with J, let D(J) be the total decrease in b. Directly iterating (2) gives the exact nonnegative expression

    D(J) = LJ + s a_1
         + h sum_(1<=i<j<=l) a_i a_j
                         (2^(j-1-i)+[j=i+1]).         (7)

In particular,

    LJ <= D(J) <= (L+h)J+s.                           (8)

For the upper bound, at each j>=2 the sum of the earlier powers is at most `2^(j-1)-1`; the extra adjacent term is at most 1. Their sum is at most `2^(j-1)`. Sum over the selected j. The initial s term is bounded by s. This also proves (8) for J=0.

After the sources, the intercept is b_J=b-D(J). After t>=1 zero sink additions the S192 law gives

    H = 2^(q+t)+b_J.                                  (9)

**Theorem 3 (lacunarity).** List the admissible exponents q increasingly and their source intercepts b_q. Apart from a finite initial part, b_q strictly decreases, and the gap between successive b_q tends to infinity.

**Proof.** For successive admissible exponents q<q', put J=(2^q-R)/B and J'=(2^q'-R)/B. Their difference is at least two, so

    J' >= 4J+3R/B.

Using (8),

    D(J')-D(J)
      >= (3L-h)J+3LR/B-s
      >= 2LJ-s.                                      (10)

Here L>=h>0. The right side tends to infinity with q. A sequence with only finitely many admissible exponents is harmless. QED.

**Corollary 4 (finite-bank no-growth theorem).** Start with any fixed finite bank of nonempty marked attachments. Allow every source/sink addition word, retain every resulting pure-dyadic response, and collect every intercept b_J in (9). This set has a finite maximum length of a run of consecutive odd integers. In particular it cannot yield unbounded exact dyadic-offset windows, despite allowing arbitrarily many changes of the base core by sources and sinks.

**Proof.** Let M be the number of bank entries. By Theorem 3, each entry has a finite exceptional initial set after which successive intercepts have gaps greater than 2M. Beyond the union of these finite exceptional sets, any interval of length 2M contains at most one point from each entry, hence at most M points. It therefore cannot contain M+1 consecutive odd integers. A run touching the finite exceptional region has bounded length as well: an unbounded such run would contain M+1 consecutive odd integers beyond that region. QED.

The offset set consists of odd integers because each offset is the difference of an odd physical H and an even power of two. Its positive part is finite: D>=0 bounds every entry above by its initial b. The argument does not assume a common modulus or a bounded return period across different bank entries.

This is not a no-go theorem for arbitrary growing windows of actual tournament counts. It concerns this precise method: pure-dyadic continuation slopes and the offsets supplied by a fixed marked-core bank under arbitrary two-sided transitive padding. No source of new arbitrary cores or additional variable apexes is included.

## 5. An explicit infinite escape, with strong physical witnesses

Take the S193 b=71 seed

    V = (12,65,26,50,35,7,61),    S={3,6},
    h=129, F=551, B=1536,
    R=512, L=335, s=47, e=32.                         (11)

The standalone verifier recomputes all these integers from literal subset path counts. Here B=3R, so d=2 and the admissible source words are `(10)^n` in creation order. The final zero is a neutral source and may be omitted, giving `1(01)^(n-1)` of length 2n-1 for n>=1. The final source is selected, so it can be reached from x.

**Theorem 5 (exponential-plus-linear dyadic family).** For every n>=1 and k>=10+2n, there is an explicitly constructed strong tournament of order k-2 with

    H = 2^k + 24+86n - 421(4^n-1)/3.                 (12)

**Construction.** Begin with (11). Add 2n-1 base sources, in creation order with apex digits `1(01)^(n-1)`. Then add `t=k-9-2n>=1` base sinks, all with apex digit zero. All arcs of the old marked tournament are preserved.

**Proof of the count.** Each two-letter source block `10` sends L to 4L+2h and subtracts L+s from b. The first block has s=47; every subsequent block starts with s=0. Summing this elementary geometric recurrence gives

    R_n=512*4^n,
    b_n=71-47 - (335+86)(4^n-1)/3 +86n
       =24+86n-421(4^n-1)/3.                         (13)

Omitting the final zero changes neither R nor b. Apply (9). The number of vertices is `7+(2n-1)+t+1=2n+t+7=k-2`. The leftmost source is selected and is beaten by x; it reaches all other base vertices. The final sink is not selected and beats x, and every other base vertex reaches it. These paths prove strong connectivity. QED.

At the minimum k the first counts are

    n=1: H=3785,   k=12, order=10, offset=-311;
    n=2: H=14475,  k=14, order=12, offset=-1909;
    n=3: H=56977,  k=16, order=14, offset=-8559;
    n=4: H=226727, k=18, order=16, offset=-35417.

This family escapes the original fixed-right-padding slope restriction. It supplies new infinite families, not a consecutive interval, prime-lane coverage, or a proof of infinitely many primes among its values. The same exact return test applies to the whole 135-seed bank; a separate finite inventory records which entries pass. The theorem itself uses only the one explicit seed (11).

## 6. Reproduction and remaining route

Run `python3 verify.py > receipt.json` and `python3 -O verify.py > receipt.optimized.json`, then compare the receipts. The verifier uses explicit errors, not removable assertions. It checks source and sink laws against complete boundary profiles on actual tournaments, normal forms, the return criterion across all short two-sided words, the exact decrement and its bounds, literal permutation controls, and the strong family on actual vertex subsets. Universal quantifiers follow from the proofs, not an extrapolation of the tested n.

The useful new proof obligation is now narrower. A source-based amplification must change more than the transitive margins of a finite marked bank, or must use non-dyadic target-dependent response equations. Retaining L explains exactly why this physical escape does not fill the missing intervals. The unrestricted `{7,21}` spectrum, an all-target construction, and an unbounded consecutive-window construction remain **OPEN**.
