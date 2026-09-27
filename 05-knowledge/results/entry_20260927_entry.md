# What universal entry must mean, and a total recursive certificate interpreter

**PROVED elementary reductions and certificate soundness / FINITE-EXACT
policy census / OPEN universal successful certification, 2026-09-27.**
No novelty claim. Write U(n)=oddpart(3n+1), for positive odd n.

## 1. Inheritance and the distinction that changes the target

The closest proved mechanism is the all-height cylinder certificate in
[the reset lifting note](reset_20260926_swaplift.md). The canonical hostile
is its expanding family 4*8^k-5. The corrected near miss is treating entry
into a region with a smaller *next* iterate as payment of the original
source's accumulated growth. The least-used sidecar is a checker rejection:
it can carry a precise missing proof obligation without declaring divergence.

The live board is original source, local descent, certified basin, exception
fuel, repeated-word macros, and the exact unread tail. The anchor is recursive
certification; the niche is a total verifier with controlled exceptions; the
wildcard is reducing the remaining convergence question to one fixed reset
interface. These are different predicates and must not be interchanged.

**Local-entry lemma.** Every positive odd orbit enters either 1 or the
one-step descending region {n>1: n=1 mod4} in finite time.

For k=v2(n+1), direct induction gives

    U^j(n)=3^j*(n+1)/2^j-1,   0<=j<=k-1.             (E1)

The first k-1 steps have division exponent1, and the last displayed value
is 1 mod4. If it is greater than1, its next U-image is smaller. This is
already a universal local-entry theorem. The last displayed value can be
much larger than n: for n=2^H-1 it is 2*3^(H-1)-1.

Consequently proving entry into a union of local descent cylinders does
not, by itself, finish Collatz. One needs a return below the retained
original source, or entry into a region whose entire future is certified
to terminate, or another well-founded rank across the whole recursion.

A minimal logical model illustrates the missing implication. Define
F(1)=1, F(even n)=n+3, and F(odd n>=3)=n-1. Every nonterminal orbit reaches
a locally descending odd state within one step, yet
2 ->5 ->4 ->7 ->6 ->9 ->... diverges. This is a countermodel to an inference,
not a counterexample to Collatz or to a source-preserving descent certificate.

## 2. A fixed low-precision interface already contains the global problem

**Entry reduction.** Every positive odd orbit enters {1} union {n=3 mod8}.
First apply U while the current value is 1 mod4 and greater than1. These
are strict decreases of positive integers, so either 1 is reached or a
value m=3 mod4 appears. Put k=v2(m+1)>=2. After k-2 valuation-one steps,

    U^(k-2)(m)=4*3^(k-2)*((m+1)/2^k)-1 =3 mod8.     (E2)

Thus convergence for all positive odd integers is equivalent to convergence
for all positive odd integers congruent to3 modulo8. This is an equivalence
about convergence, not an equivalence with merely having a locally
decreasing step somewhere in that residue class.

Every n=3 mod8 admits the exact split

    (q,r,u)=(1,1,n-2),  v2(3u+1)=2,  R=2-1=1.       (E3)

It therefore lies at the inherited difficult swap interface (q,R)=(1,1).
This is a legal exact representation; it is not claimed to be the state
selected by every canonical decomposition history. Its unbounded register
u is essential. Fixed interface labels do not make the residual arithmetic
problem finite, but they do identify a precise target for a recursive rule.

## 3. Certified jumps plus a rigorously controlled exception budget

Freeze the 171 rows in [the reset bank](reset_20260926_swaplift.out).
Each row has a proved residue guard and a finite actual U segment ending
below its input, at arbitrary height. All its guarded sources are 3 mod4;
every segment has length at most41. Define a deterministic policy:

1. Accept at1.
2. At n=1 mod4, apply the actual one-step decrease.
3. Otherwise use the first matching bank row, in its saved order, to jump
   along its actual orbit to the proved smaller endpoint.
4. If no rule matches, mark the current state unresolved by this policy.

For each integer d>=0, define C_d as the set accepted by this interpreter:
at an unresolved state with d>0, take one exact U step and replace d by d-1;
at an unresolved state with d=0, reject. Certified jumps keep d unchanged.

**Totality inequality.** Every executed recursive call strictly lowers
the positive integer

    W(n,d)=n*2^d.                                    (E4)

For a certified jump n decreases. For an exception n>1 and
U(n)<2n, hence U(n)*2^(d-1)<n*2^d. The only remaining actions halt.
This is an explicit inequality controlling insufficient-precision moves
inside a total recursive verifier. Its input fuel is not a proved bound
on the exceptions needed by an arbitrary orbit. Rejection is permitted.

**Soundness and exact scope.** The C_d are nested, decidable sets of
positive odd integers, each contained in the basin of1. Indeed acceptance
concatenates finitely many actual U segments. Conversely every convergent
orbit is accepted at some finite d: following the deterministic policy
can skip actual segments but still reaches1, with finitely many exceptions.
Therefore

    union_(d>=0) C_d = the actual basin of1.           (E5)

Neither totality of each checker nor (E5) establishes that this union is
all positive odd integers. A uniform source-computable successful budget,
or a different total successful grammar, would supply the missing theorem.
This states the coverage obligation rather than presuming it.

## 4. No fixed exception budget suffices, even on a known convergent family

The maximum bank modulus exponent is65 and none of its rows contains -1
modulo its own modulus. For H=d+65 and n=2^H-1, the first d+1 policy states
are

    n_j=3^j*2^(H-j)-1 =-1 mod2^65,  0<=j<=d.         (E6)

All are unmatched and 3 mod4, so C_d rejects n. This proves that no fixed
C_d covers all positive odd inputs; it does not assume convergence of
these arbitrary Mersenne sources.

The independent [entry audit](entry_20260927_entry_review.md) strengthens
this to a known convergent family. For H>=2 set

    M=3^(H-1),  t_H=(2^M+1)/3^H,  N_H=2^H*t_H-1.    (E7)

Factoring successive cubes gives v3(2^(3^(H-1))+1)=H, so t_H is a positive
odd integer. Its first H-1 U steps have exponent1, with
U^j(N_H)=3^j*2^(H-j)*t_H-1. At the last such state,

    3*U^(H-1)(N_H)+1=2*(3^H*t_H-1)=2^(M+1),

so U^H(N_H)=1. With H=d+65 the same unmatched prefix as (E6) forces
rejection by C_d. Thus every C_d is strictly smaller than the convergent
basin, unconditionally. This is a specialization of inverse-word
completion, not an obstruction to all finite recursive grammars.

There is a second witness within the new
[27-containing family](entry_20260927_families.md). With k blocks remaining,
its first2k iterates exceed the current source. For k>=21 no bank rule
of length at most41 can apply at that even phase. Its next state is 1 mod4,
so the quarter rule returns to the next even phase. At least k-20 exceptions
are forced. Choosing k=d+21 gives another convergent source rejected by C_d.
The [macro calculus](entry_20260927_recursive.md) certifies this same family
with a bounded number of parameterized rule nodes: atomic fuel and recursive
certificate size are fundamentally different measures.

## 5. A finite census, not a universal extrapolation

The [script](../../04-computation/experiments/entry_20260927_entry.py) and
[saved output](entry_20260927_entry.out) use exactly the inherited bank,
not the new cylinders introduced in this session. On all49,999 odd sources
3..99,999 the exact deterministic policy completed within the declared
cap of1,024 decisions per newly traced path. Counts accepted at each fuel:

| Fuel d | Accepted sources |
|---:|---:|
|0|15,045|
|1|18,303|
|2|21,191|
|4|26,218|
|8|32,510|
|16|47,740|
|32|49,982|
|64|49,999|

The largest needed fuel in this finite universe is36, first attained by
52,527; source27 needs13. The total interpreter was independently compared
with the finite memoized census for all2,047 odd sources3..4,095 and all
eight displayed fuel levels. There are also16,383 entry-reduction controls,
65 Mersenne rejection controls, and100 grow/drop pairs of the abstract countermodel.
All checks remain active under optimized Python. These finite counts give
no universal fuel bound.

Reproduce from the repository root:

    python -X utf8 -B 04-computation/experiments/entry_20260927_entry.py
    python -X utf8 -B -O 04-computation/experiments/entry_20260927_entry.py

The sharper next question is to replace repeated exception charges by a
source-decodable macro whose guard both holds on the actual integer and
pays the retained root threshold. An unbounded valuation counter is allowed;
a branch that silently changes the source, discards its unread tail, or
resets its comparison threshold is not. The all-height guarded cylinders
and the closed family through27 provide positive cases of this target.
