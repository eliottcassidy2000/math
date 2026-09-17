# Completing every finite halving word to a chosen target

**Status: PROVED elementary construction and obstruction; FINITE-EXACT controls.**
Collatz convergence and completeness of any signed cycle list remain **OPEN**.
No novelty claim is made. The independent algebra audit was performed by a
second research agent; exact controls are reproduced below.

## Inheritance and scope

The closest mechanism is the full ternary inverse-fibre odometer in
[the first Collatz braid note](arithmetic_braids_20260917_collatz.md), together
with its ordered-word carry. The hostile example is the three disjoint known
positive basins for `3n-1`. The corrected near miss is treating full residue
coverage as coverage of all integers. The needed sidecar is **ordinary size**,
in addition to binary prefix and ternary residue.

For odd `b` with `3∤b`, define on nonzero odd integers

    T_b(n)=(3n+b)/2^v2(|3n+b|).

The numerator never vanishes at an integer under these assumptions. Targets
`u` below are nonzero odd integers with `3∤u`. For positive Collatz use `b=1`;
for its reflected negative system use `b=-1` on positive integers. Other b
need not preserve a half-line under unrestricted iteration. Our constructed
finite path has the sign of u once its free parameter is sufficiently large.

## 1. A prescribed word can end at any admissible target

Fix any finite word `w=(k1,...,kL)` of positive integers, including the empty
word. Set `K0=0`, `Ki=k1+...+ki`, `K=KL`, and

    B=sum_(i=0)^(L-1) 3^(L-1-i) 2^Ki.                         (1)

The affine iteration identity is

    2^K v = 3^L n + bB.                                      (2)

Choose `k0` in `{1,2}` so `2^k0 u ≡ b (mod3)`. Every integer `t>=0` gives
an odd immediate predecessor

    v_t=(2^k0 4^t u-b)/3,          T_b(v_t)=u.                 (3)

To complete the word before this predecessor, put

    n_t=[2^(K+k0)4^t u-b(2^K+3B)]/3^(L+1).                  (4)

**Theorem.** Exactly one residue class of t modulo `3^L` makes (4) integral.
Every sufficiently large nonnegative t in this class gives an odd integer
n_t whose first L exact halving exponents are w, followed by exponent
`k0+2t` and target u. All nodes on this finite path have u's sign.

**Proof.** Integrality is equivalent to

    4^t ≡ b(1+3B 2^(-K))/(2^k0 u)  (mod3^(L+1)).             (5)

The right side is a unit congruent to1 mod3. By lifting the exponent,
`v3(4^m-1)=1+v3(m)` for nonzero m. Thus 4 generates precisely the
`3^L` principal units modulo `3^(L+1)`, proving existence and uniqueness.

The numerator in (4) is odd. To see that total integrality gives every
intermediate integer, work backwards from v_t. Modulo3, (2) implies
`2^kL v_t ≡ b (mod3)`, since `B ≡ 2^K_(L-1) (mod3)`.
Therefore `(2^kL v_t-b)/3` is an odd integer. Repeating proves each earlier
node integral and odd. Consequently each prescribed power of2 is exact,
since the next node is odd. Each intermediate node is a positive multiple
of `4^t u` plus a fixed constant, so eventually all have u's sign. ∎

This is an explicit construction, not an application of the Collatz
conjecture. The script finds t by lifting one ternary digit at a time.

## 2. Every ternary residue occurs inside the completion family

Fix one admissible t0. For `t=t0+3^L j`, subtraction in (4) gives

    n_(t0+3^L j)-n_t0
      =2^(K+k0)4^t0 u (4^(3^L j)-1)/3^(L+1).

Thus, for distinct integers j1,j2,

    v3(n_(t0+3^L j1)-n_(t0+3^L j2))=v3(j1-j2).             (6)

The map from j to n is a bijection modulo every `3^s`. Each desired
residue is attained infinitely often, with j in one class modulo `3^s`.
All these n have the same word w, which fixes a single odd residue modulo
`2^(K+1)`. Arbitrary additional binary data inconsistent with w cannot be
prescribed; this is a necessary compatibility condition.

**Mixed-residue corollary.** For every `H>=1,s>=0` and every odd residue r
modulo `2^H 3^s`, there are infinitely many same-sign ancestors of u in r.

For proof, take an integer representative of the desired odd binary class,
and read a finite orbit word until its total exponent K is at least `H-1`.
No zero numerator occurs. Its word fixes that binary class. Equation (6)
then sets the desired ternary residue without changing the binary class.
This proves density for the topology defined by these finite quotients.
For the full shortcut graph, adjoining even doublings gives all classes,
including even ones. It does not assert density for arbitrary moduli.

For `b=-1`, use `u=1,5,17` from the known cycles

    (1),        (5,7),        (17,25,37,55,41,61,91).

Their three disjoint basins **each** meet every odd class modulo `2^H3^s`.
Moreover every finite halving word occurs in each basin. Therefore no
classifier using only a fixed amount of binary/ternary residue information,
even together with a fixed finite halving prefix, distinguishes these basins.

Small direct examples, all congruent to1 modulo72, are

    1 -> 1,
    73 -> 109 -> 163 -> 61 -> ... -> 17,
    361 -> 541 -> 811 -> 19 -> 7 -> 5.

These are independent forward checks, not instances inferred from a census.

## 3. Why this does not close the global conjecture

Three separate losses are visible rather than hypothetical:

1. **Changing integers.** Completing a prefix produces a new integer. It
   says nothing about the infinite tail of a fixed starting integer.
2. **Size.** For fixed w and a prescribed ternary class, t advances by
   `3^(L+s)`, while `|n_t|` grows like `4^t`. This particular construction
   supplies only `log X/(3^(L+s) log4)+O(1)` sources with `|n_t|<=X`.
   The family has natural density zero despite full residue support.
3. **Missing primes.** For `b=-5`, all ancestors of u=5 are divisible by5.
   They are dense in the stated binary/ternary quotients but omit all units
   modulo5. The gcd stratum is a real additional invariant.

The source family has an informative boundary in the 2-adic integers:

    lim_(t->infinity) n_t = -b(2^K+3B)/3^(L+1).              (7)

Applying its first L affine branches reaches `-b/3`, where `3v+b=0`.
The forced terminal exponent tends to infinity. Thus these ordinary-size
diverging completions accumulate 2-adically on a preimage of the singular
point. This explains how arbitrary finite symbolic compatibility can coexist
with poor height control. The limit is not a fixed positive integer orbit.
More precisely, `v2(n_t-limit)=K+k0+2t`, whereas
`v3(limit)=-(L+1)`: the limit is not even a 3-adic integer. The completion
family is unbounded in ordinary size, convergent toward a 2-adic singularity,
and a residue-covering isometry in its reindexed 3-adic coordinate. These
three behaviors belong to different topologies and must be kept distinct.

Krasikov and Lagarias supply a substantially stronger **CITED** counting
result: for every fixed positive target a not divisible by3, at least
`X^0.84` integers at most X reach a, for sufficiently large X depending on a.
Their inequalities retain a size coordinate as well as ternary residue.
See [the authors' paper](https://arxiv.org/pdf/math/0205002), abstract and
Sections2,6. This is historical attribution, not a claim of a current best
bound. Neither their lower bound nor this completion construction proves
that every positive starting integer reaches1.

## 4. Reproduction and decisive controls

    python 04-computation/experiments/arithmetic_braids2_20260917_inverse_completion.py
    python -O 04-computation/experiments/arithmetic_braids2_20260917_inverse_completion.py

The frozen [JSON](../../04-computation/experiments/arithmetic_braids2_20260917_inverse_completion.json)
contains 7,623 completions for every word of length0..4 with exponents1..3,
seven specified `(b,u)` pairs and every source residue modulo9; 1,116
additional complete mixed-cylinder checks for H1..5 and s2; direct mod72
basin witnesses; and the mod5 exclusion control. Integer valuations and
every forward step are checked independently of the congruence solver.
All checks remain active under optimization. The maximum constructed source
in the first universe has 1,464 bits. The infinite claims follow from the
proof, not from extending this finite universe.

An [independent referee](../../04-computation/experiments/arithmetic_braids2_20260917_inverse_audit.py)
imports no completion code and replaces digit lifting by brute exponent
search. Its 4,000 parameter cases check 108,000 full trajectories, 104,000
valuation comparisons, and 16,000 singular-limit identities using exact
rational arithmetic. Its complete universe and separate JSON are frozen;
normal and optimized replays agree.

## 5. Next precise problem

Replace existence of an ancestor in every cylinder by a **height-controlled
cover** that applies to a fixed integer. A useful candidate certificate must
survive both hostile controls: three dense disjoint minus basins, and long
exponent-one growth prefixes. A finite residue statistic alone has already
lost the coordinate needed for such a certificate.
