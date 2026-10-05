# A four-credit guard from the 5/8/9 corridor

2026-10-05 (UTC). **PROVED:** a smaller common-future dependency on the whole
guard 219 mod 256, a sharp four-credit funding rule, exact guard composition,
and first-hit-safe common-future compilation. **FINITE-EXACT:** the declared
experiments. **OPEN:** guard coverage and root completion for arbitrary
positive integers.

Artifacts: [program](../../04-computation/experiments/paid_guard_budget_20261005.py)
and [saved output](paid_guard_budget_20261005.out).

## 1. Recover the objects before transporting a count

The closest proved mechanism is the [adaptive credit controller](adaptive_credit_potential_20261004.md):
H earns two credits, the source-sensitive four-slot bank P earns one or two,
and the actual word G=(1,2) spends one. The invariant is the rational account
E(x,k)=(x+5)(9/8)^k with immutable original source N.

The canonical hostile is 27: G(27)=31, whose smaller inverse-G predecessor
is just 27 again. The corrected near miss is to call that round trip payment.
The least-used useful sidecar is the guard of an intermediate dependency
after affine cancellation. Our board is **source/exponent type; signed
anchor; binary guard; ternary inverse depth; original debt; retained proof
expression**. Anchor: enlarge the paid library. Niche: lossless guard
composition. Wildcard: exact 5/8/9 clocks from other domains.

The [cross-repo recovery](five_eight_nine_transfer_20261005.md) distinguishes
the eight-unit golden clock in F9, the eight fixed-path tournament masks and
their five-path strong class, and the arithmetic G=(9x+5)/8. Each has a
specified map and a loss boundary. The clock at the prime 3 is invertible;
G's shifted ternary action is contracting. It cannot inherit that clock's
orbit coverage merely because both involve 8 and 9.

We used META-PATTERNS' **Turn certificate failure into an address, then change
sidecars** and **Treat pulls and collisions as synthesis opportunities**.
Three independent probes found the same 219 cell; the contributions below
separate its residue partition, the general exit bank, and credit/certificate
transport rather than counting the rediscovery three times.

## 2. Two meanings of 27 mod 64, and the actual branch split

The inherited [four-slot result](collatz_four_slot_compression_20261004.md)
pays n=3^e when **the exponent** e is 27 mod 64, through the source cell
187 mod 256 and actual word (1,2,1,3).

For **the integer source** n=27+64t, the first four valuations instead are
(1,2,1,1). In fact v2(n+5)=5, so exactly one complete G block is available
before the pattern changes. The [exact residue tree](paid_guard_cover_tree_20261005.md)
splits the next valuation into

| Source cell | Fifth valuation |
|---|---|
| 27 mod 128 | 1 |
| 91 mod 256 | 2 |
| 219 mod 256 | at least 3 |

The first two rows still grow at this point. This note pays the entire third
row. Its exact-3 half, 219 mod 512, pulls back to e=83 mod 128 for powers
of three, a new 1/16 of the specified exponent domain e=3 mod 8 against the
frozen earlier banks. The at-least-4 half already belonged to an earlier
bank. These are paid dependencies, not supplied root proofs for every child.

## 3. Pay once, then recover an inverse-growth credit

For every integer t>=0 define

    n = 219+256t,
    K(n) = (81n+53)/128 = 139+162t,
    L(n) = (9n-3)/16 = 123+144t.

All three are positive odd and L(n)<K(n)<n. The source's actual word
v=(1,2,1,1) gives

    F_v(n) = (81n+85)/32 = 4K(n)+1.

Since 3(4K+1)+1=4(3K+1), U(F_v(n))=U(K(n)). Also L(n)=11 mod 16,
so its actual word (1,2) reaches K(n). Hence L(n), K(n), and n have
explicit common futures. In particular, a supplied root proof for either
smaller child transports to n.

This is where applying the inverse-G rule actually helps:

    K(n)=4 mod 9,
    (8K(n)-5)/9=L(n),
    K(n)+5=18(8+9t),
    L(n)+5=16(8+9t).

Exactly one integral inverse-G step exists: the first shifted value has
ternary valuation two, and the second has valuation zero. Unlike the
27->31->27 round trip, the quarter-child operation already paid the original
source before taking this inverse step.

For a=v2(3K(n)+1), the actual common-future words are

    source n: (1,2,1,1,a+2),
    child L:  (1,2,a).

They have the same endpoint U(K(n)). No prefix pads a root: the intermediate
source states exceed n, while L and its first two odd images exceed 1.
If their common endpoint is 1, it is the first hit.

The original n is 3 mod 4. Consequently L(n)<n also pays the inherited
[proper rank](collatz_branch_toll_rank_20261004.md), by the same energy bound
used in the earlier controller. It is not merely a decrease relative to an
intermediate enlarged state.

## 4. Four credits, sharp payment, and finite adaptive episodes

Keep E(x,k)=(x+5)(9/8)^k. The new dependency satisfies

    L(x)+5=(9/16)(x+5)+2.

Its native least source is 219. Thus

    E(L(x),k+4)/E(x,k)
       =(9/8)^4 [9/16+2/(x+5)]
       <= (9/8)^4*(4/7)
        = 6561/7168 < 243/256 < 1.

The bound is attained at 219. The intermediate K can instead earn three
credits with exactly the same sharp factor. Indeed the **guarded identity**
K=G composed after L gives

    E(K(x),k+3)=E(L(x),k+4).

Add L to the earlier controller, retaining every native guard and every
actual credit award. If h,l,p count H,L,P funding actions and g counts G,
then every nonempty funded episode obeys

    6 <= E(x,k) <= (N+5)(243/256)^(h+l+p) < N+5,
    g <= 2h+4l+2p,
    x<N and R(x)<R(N).

The first action must be a funding action and puts N=3 mod 4. These bounds
prove that no episode contains infinitely many executed legal actions, under
arbitrary choices and pattern switches. They concern completed actions;
actual Collatz paths inside a dependency may exceed N. An unbounded search
for an unavailable next pattern is a separate algorithmic issue.

The new entry guard 219 mod 256 is disjoint from the H guard 155 mod 2048
and the four-slot P bank: its first four valuations are 1211, containing
three ones. It therefore enlarges the **previous H/P/G controller's** entry
domain. It does not mean every point of that guard was unhandled by all
previous repository constructions.

**A fifth credit really fails.** The composition G^5 after L is

    (531441n+1925333)/524288 > n

for every positive n. It is fully legal on 932059 mod 1048576; the least
source gives 932059->944779, increasing both numerical value and inherited
rank. Thus four is the greatest universal integer award per L for this
G-denominated size-payment account. This is a legal overspending example,
stronger than a failure of the chosen potential inequality alone.

Uniform H/L receipts now have increments +2,+4,-1, giving ordered trees
with ternary and five-child funding nodes. This follows from pending child
slots, not from identifying the five-child tree with a five-path tournament.
The full typed receipt retains each native guard. A state-sensitive P node
also retains its entry condition and actual one/two-credit award.

## 5. Cancellation loses three guard bits unless they are stored

The affine formula L(n)=(9n-3)/16 has positive odd integer outputs already
on n=27 mod 32. The proved dependency is only on 219 mod 256. For example
the formula sends 27 to 15, but 15 starts with valuations (1,1), so the
required child prefix (1,2) fails. This does not assert that 27 and 15 lack
some eventual common future; it refutes that local certificate without its
guard. The quotient has discarded exactly three binary source bits.

Store an affine map together with its exact source residue. If u is
(P_u x+B_u)/Q_u on r_u mod m_u, and v has guard r_v mod m_v, its chronological
composite has source domain

    x=r_u mod m_u,
    P_u*x+B_u=Q_u*r_v mod Q_u*m_v.

For these dyadic guards P_u is odd, so the second congruence has one residue;
intersect the two by agreement modulo the smaller modulus. An inconsistent
pair is an empty guard, not a discarded condition. Keep the larger modulus
when they agree. This gives an exact associative composition rule, because
it is literal intersection with a preimage. The guarded equality LG=K
restores the three missing bits through G's native input condition.

No claim is made that the enlarged H/G/L affine alphabet losslessly decodes
its entire history. Store its typed expression until such injectivity is
proved. The earlier pure H/G decoder is not silently extended.

The program additionally compiles any legal H/G/L expression into actual
common-future words. Traverse it backward. G prepends 12; H uses the inherited
quarter-child rewrite; L replaces a child prefix (1,2,a) by (1,2,1,1,a+2).
If the retained common-future tail is too short, extend both routes together
until the required prefix is present. Native L outputs force their first
two steps to be 12, with nonroot intermediate states, so this extension
cannot require padding a root. Relative to the terminal route, the compiled
source route has 6h+2g+2l additional odd steps and 10h+3g+4l additional halvings.

## 6. What broadens coverage, and what remains

The [all-parameter boundary bank](paid_guard_27_boundary_20261005.md) retains
the exact height inequality at the first coefficient-contractive exit of
(12)^q1^r. This adds infinitely many paid rows, including the new 219 mod512
half. Its carefully scoped density computation raises the named coverage
of power-three exponents e=3 mod8 from 42.7206 percent to 53.7448 percent,
including the single new 83 mod128 class exactly once. Its additional
odd-source density is approximately 0.00689011741468514; after the declared
ternary-bank overlap, the named residual is approximately 0.105842463387755.
Use that note's exact bounds and named comparison bank for the total;
the present credit theorem does not independently add those densities.

The stronger L operation improves the dependency and the available credit
on an already declared source cell. General paid-guard coverage and supplied
root proofs for arbitrary exit states remain open. The unresolved next
objects include the two growing branches of the 27 mod64 split, and the
residual exponent cells outside the expanded bank.

Reproduction:

    python -X utf8 -B 04-computation/experiments/paid_guard_budget_20261005.py
    python -X utf8 -B -O 04-computation/experiments/paid_guard_budget_20261005.py

All 1,092 H/G/L words of lengths one through six are checked at lifts 0,1,7:
3,276 native/composition checks and actual first-hit-safe compiled joins,
including 2,142 funded instances with every prefix compared against its
immutable source. Seven additional local lifts verify both dependency
words and exact ternary inverse depth. The tests retain unfunded words,
the lost-bit witness, a native fifth-credit size/rank failure, and four
malformed/guard/budget controls. Arithmetic is exact; checks use exceptions
and run unchanged under optimization. Computed carriers are caches with
provenance, not self-authenticating arbitrary supplied dataclass fields.
