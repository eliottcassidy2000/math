# Independent audit of the entry reduction and bounded-exception interpreter

**VERIFIED proofs and FINITE-EXACT reproduction, 2026-09-27.** Scope:
positive odd integers, U(n)=oddpart(3n+1), nonnegative integer fuel, and
the fixed171-row certificate bank in
[reset_20260926_swaplift.out](reset_20260926_swaplift.out). The review
accepts the new entry and interpreter claims with the distinctions below.
No universal Collatz conclusion follows.

## 1. Entry into the fixed hard-swap interface is valid

For n>1 with n=1 mod4, U(n)<n. Repeating this branch therefore ends at1
or a value m=3 mod4. In the latter case put k=v2(m+1)>=2. For
0<=j<=k-1,

    U^j(m)=3^j*2^(k-j)*((m+1)/2^k)-1.

All steps before j=k-1 have division exponent1. At j=k-2 the value is
3 mod8, because it is4 times an odd integer minus1. This establishes
entry into{1} union{n:n=3 mod8}; it makes no assertion that the endpoint
is below the original starting source.

For every n=3 mod8, the retained split q=r=1, u=n-2 is positive and
odd and has v2(3u+1)=v2(3n-5)=2. Consequently it is exactly the hard
swap with recreated precision R=1. This is an all-source reduction to
that interface, not a proof that the interface has a terminating rank.

The separate assertion that every odd orbit reaches1 mod4 is also valid:
for k=v2(n+1), apply k-1 valuation-one steps. The endpoint's next step
decreases it unless it is1. Repeatedly entering a locally decreasing
region does not establish descent below the original source. The script's
countermodel, even n->n+3 and odd n>1->n-1, gives the exact obstruction:
starting at2, every pair of moves raises the even state by2 while every
visited odd state locally decreases.

## 2. The interpreter is total and characterizes convergence by its union

Its deterministic policy takes a certified numerical decrease without
spending fuel, or takes one actual U step and spends one unit of fuel.
An unmatched state at fuel0 is rejected. All uses of the inherited bank
are on its proved positive residue classes, including the separately
verified small sources. The parsed row order q,J,A,P,b,clip,R,residue,K
matches the frozen table.

For n>1, U(n)<2n. Therefore the positive integer rank

    rank(n,d)=n*2^d

strictly decreases both on a certified step (n decreases) and on an
exception (d decreases and the new value is below2n). Thus the interpreter
halts for every positive odd input and every d>=0. The domain assumptions
matter: the proof does not concern negative fuel or arbitrary rational n.

Let C_d be its accepted sources. Acceptance provides a finite concatenation
of genuine U segments ending at1, so C_d is contained in the convergent
basin. Conversely a source whose orbit reaches1 uses only finitely many
policy decisions and therefore finitely many exceptions; skipping several
actual steps cannot skip past the absorbing U fixed point1. Enough fuel
accepts it. Consequently

    C_d subset C_(d+1),
    union_(d>=0) C_d = {positive odd n whose orbit reaches1}.

These are unconditional statements. Equality of the union with *all*
positive odd integers is precisely the remaining Collatz assertion.

## 3. The original Mersenne hostile and its logical limit

The bank has maximum modulus depth Kmax=65, maximum selected horizon41,
and no row's residue equals-1 modulo its own modulus. If H=d+65 and
n=2^H-1, its first d+1 policy states are

    n_j=3^j*2^(H-j)-1, 0<=j<=d.

They are3 mod4 and-1 modulo2^65, hence match neither the quarter branch
nor any bank row. The interpreter rejects n at fuel d. This proves that
no finite C_d covers all positive odd integers.

The Mersenne argument alone does not prove that a rejected witness
converges. In particular it alone cannot justify strict containment in
the convergent basin. The following explicit completion repairs that
gap without assuming Collatz or relying on the new47-completion theorem.

## 4. A known-convergent witness outside every C_d

For H>=2 set

    M_H=3^(H-1),
    t_H=(2^M_H+1)/3^H,
    N_H=2^H*t_H-1.                                 (R1)

These are exact finite integer expressions, though their sizes grow
rapidly. The numerator of t_H has3-adic valuation exactly H. Here is an
elementary proof that needs no external valuation theorem. For x=2 the
value x+1 has valuation1. Replacing x by x^3 multiplies x+1 by
x^2-x+1. Whenever x=-1 mod3, writing x=-1+3s gives

    x^2-x+1=3*(1-3s+3s^2),

whose valuation is exactly1. Induction gives the claimed valuation.
Thus t_H is a positive odd integer.

For 0<=j<=H-1 the actual states are

    U^j(N_H)=3^j*2^(H-j)*t_H-1.                    (R2)

Before the last state these are3 mod4 and the division exponent is1.
The phase coordinate H-j is v2(U^j(N_H)+1); it decreases by1 while
the odd coefficient changes from3^j t_H to3^(j+1)t_H. This is the correct
rank and update, not a claim that a fixed t survives that phase. At the
last state,

    3*U^(H-1)(N_H)+1
       =2*(3^H*t_H-1)=2^(M_H+1).

Hence U^H(N_H)=1, and none of the earlier states is1. Every N_H is
therefore proved convergent, with exactly H odd steps to1.

Now take H=d+65. At every 0<=j<=d, (R2) is-1 modulo2^65 and3 mod4.
The same bank-exclusion proof as above applies. Thus

    N_(d+65) is convergent but is not in C_d          (R3)

for every d>=0. This proves, unconditionally, that *each fixed C_d is a
proper subset of the convergent basin*. It does not assert that every
adjacent inclusion C_d subset C_(d+1) is strict. The exception cost is
unbounded even on this explicit known-convergent family.

The construction illustrates why the all-ones nonuniformity argument
needed its arithmetic sidecar: one can complete a growing initial word
to a known terminal value, while retaining enough trailing ones to defeat
the entire fixed finite bank. The astronomical final division is selected
by the exact numerator in(R1), not predicted for an arbitrary source.

## 5. The completed47-family gives a second valid strictness proof

The new [family note](entry_20260927_families.md) constructs convergent
sources of the form a*8^k-5 with positive even a and arbitrarily large k.
At an even phase with t>=21 blocks remaining, every one of its next2t
actual U iterates is above that current source. Since every bank horizon
is at most41, no bank certificate can match the phase source. It is also
3 mod4, so the interpreter must spend one exception. The intermediate
state is1 mod8; the prioritized quarter rule then takes its one-step
decrease to the next even phase.

There are k-20 such forced exceptions for k>=21. Choosing k=d+21 gives
another convergent source outside C_d. This argument is sound and uses
the exact even-coefficient guard: a*8^t-5 is11 mod16, so its next division
exponents are exactly1 and2. It is not enough merely to know that a family
has long total stopping time; the source-relative growth at *each* remaining
phase is what excludes every bank rule.

## 6. Reproduction and review boundary

The original entry script was run independently both normally and under
Python-O. Both runs exactly reproduced its saved output:49999 sources
below100000 completed the finite policy census, maximum observed exception
cost36 at52527, and the bounded interpreter matched the debt calculation
on2047 sources at8 fuel levels. No finite cap is asserted universally.

The separate [review script](../../04-computation/experiments/entry_20260927_entry_review.py)
checks the inherited bank constants and residue exclusions; constructs
N_H explicitly only for2<=H<=10; checks exact3-adic divisibility modulo
3^(H+1) for2<=H<=129; and checks2145 modular policy states for65
known-convergent witness formulas. It never constructs enormous values
2^(3^(H-1)) in the large-H controls. An additional128 even-coefficient
controls verify the42-step growth mechanism of section5.

Commands:

    python -X utf8 -B 04-computation/experiments/entry_20260927_entry.py
    python -X utf8 -B -O 04-computation/experiments/entry_20260927_entry.py
    python -X utf8 -B 04-computation/experiments/entry_20260927_entry_review.py
    python -X utf8 -B -O 04-computation/experiments/entry_20260927_entry_review.py

The [review output](entry_20260927_entry_review.out) matches under normal
and optimized execution. All checks use explicit exceptions, not removable
assertions. This review audits the entry/interpreter composition and its
logical scope; it inherits the already audited all-height bank theorem
rather than claiming a fresh proof of every underlying bank inequality.

## 7. Independent audit of the three-type board

For an input n=3 mod8, L=v2(n+5)>=3. The definitions
k=floor((L-1)/3), r=L-3k and b=(n+5)/2^L give exactly one representation
n=2^r*b*8^k-5 with r in{1,2,3} and positive odd b. The inherited guarded
block proves, for0<=j<=k,

    x_j=U^(2j)(n)=2^r*b*9^j*8^(k-j)-5.

Consequently every block preserves2v2(x_j+5)+3v3(x_j+5): the two
valuations change by-3 and+2. There are two exact bookkeeping choices.
One keeps the original b fixed and retains the elapsed counter j. The
other recomputes the canonical odd coefficient b_j=9^j*b and updates it
to9b_j as the remaining k decreases. Claiming that the canonical odd
coefficient itself stays fixed would lose the transferred ternary factor.

Put c=b*9^k at the boundary. The three exact returns are:

| r | Boundary | Verified return |
|---|---|---|
|1|m=2c-5=1 mod4|U(m)<m; here k>=1, so m>=13>1 |
|2|m=4c-5=7 mod8|exactly1+v2(c-1) exponent-one steps lead to1 mod4 |
|3|m=8c-5=3 mod16|U(m)=12c-7, and U^2(m)=oddpart(9c-5) |

The type2 residue is **7 mod8**, not3 mod8. The minimal original27-family
example is27->41->31, with boundary31. Its stated run-length formula was
already correct: v2(m+1)=2+v2(c-1), so the exponent-one run has one fewer
steps. The valuation is well defined since type2 has k>=1 and c>=9.

The new all-k successful family is the *guarded subset* of type3 satisfying
the additional congruence in(B2), not every type3 source. For example
786427 has L=18,k=5,r=3,b=3, but fails that guard and its proposed selected
return. The unguarded two-step boundary identity remains exact. These are
the two literal repairs and one coordinate clarification reported to the
board author; no change to the underlying block formulas is needed.

The board's section4 distinctions are accepted: universal entry into a
local descent region need not repay the original source; the bounded-fuel
rank proves total checking with rejection allowed; the accepted union is
exactly the convergent basin; and omitting convergent sources at every
fixed fuel level does not preclude a finite parameterized macro grammar.
This section is an algebraic audit; the existing review output is unchanged
and no new finite census is attributed to its script.
