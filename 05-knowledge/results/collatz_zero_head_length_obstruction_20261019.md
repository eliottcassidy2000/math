# A rigorous length obstruction to a Collatz zero-head counterexample

**Status, 2026-10-19: PROVED necessary condition; unrestricted conjecture
OPEN.** For the expanding bank {K,V=K^4F} and reset bank {F,H}, every fully
guarded word with a nonempty expansion-bank segment followed by a nonempty
reset-bank segment, contracting total multiplier, and
length

    N < 918644030413102796640970006656455          (1)

has zero excluded head: every positive source in its cylinder descends.
Any counterexample must have at least this many Terras steps. This is a
symbolic obstruction, not a census of words up to that length. No universal
zero-head proof or genuine counterexample is established here.

## 1. Inherited setup

The [bank theorem](collatz_bank_reset_20261019.md) and its
[modular refinement](collatz_modular_bank_head_20261019.md) concern exact
lawful cylinders, not arbitrary source entry. Use the inherited adjacent
branches

    K(x)=(Px+C_0)/Q,
    P=3^19, Q=2^30, C_0=4912087879,
    F(x)=(3x+1)/4,
    H(x)=(81x+101)/512,
    V=K^4 F.

The expansion bank is {K,V}; the final reset bank is {F,H}. Both segments
are required to be nonempty. This hypothesis is necessary: with an empty
expansion segment, F alone fixes its least positive source 1 and would
contradict zero-head. Words are read
chronologically. A nonempty expansion-bank word followed by a nonempty reset word
has an exact complete parity-cylinder map

    T^N(x)=(Sx+B)/M,
    S=3^h, M=2^N, Delta=M-S>0,
    x=rho+tM, 0<rho<M, t>=0.

The inherited carry envelope is

    0<B/M<E, E=C+1,
    C=83829391903560839836419813590128166431
            /157489106280555890461927692448373459. (2)

The zero-head conjecture asks whether Delta rho>B always holds. Since the
map has contracting slope, descent at rho implies descent at every later
lift. Conversely, any non-descending positive lift implies a non-descending
least lift. Adjacent companions have the same exact cutoff, as proved in
the bank theorem. We retain all full-word source guards.

## 2. A forced K^4 prefix, or immediate descent

Any nonempty expansion-bank word either is K, K^2, or K^3, or its chronological
parity word begins with four complete K blocks. Indeed, if it contains no
V it is a power of K. If it contains V, its first V starts with K^4,
preceded only by zero or more K blocks. This is a structural property of
the actual bank alphabet; no probability model is used.

The three exceptional expansion words already descend after their first
reset, at every lawful positive source. To see this, put r_1=815464443,
the least positive K source. For all x>=r_1,

    F(K^l(x)) <= F(K^3(x)) < x,  l=1,2,3.        (3)

The first comparison uses that K is increasing and strictly above the
identity on positive reals. For the second, the affine slope of F K^3 is
3P^3/(4Q^3)<1, and the exact evaluation at r_1 is below r_1. Both integer
inequalities are verified by the checker. Also H(y)<F(y) for y>0, since
81/512<3/4 and 101/512<1/4. Thus either first reset descends. Every later
lawful F or H block is nonincreasing on positive integers: F(y)<=y for
y>=1 and H(y)<y for y>=1. Once below the initial source, the endpoint
therefore stays below it throughout subsequent complete reset blocks.
No first-stopping-time assertion about intermediate H steps is needed.

Every potential counterexample must consequently begin with K^4. The
exact K^4 carry and least positive source are

    C_4=C_0(P^4-Q^4)/(P-Q),
    r_4=-C_4(P^4)^(-1) modulo Q^4
       =588943162011782273059357097228697595.     (4)

Its source is therefore at least r_4. The checker verifies the exact
residue and its range; a complete word's cylinder is contained in this
prefix cylinder. This source bound is much stronger than merely r_1.

## 3. Non-descent would require an exceptionally small linear form

Suppose a guarded positive source x fails to descend. Set s=S/M<1 and
e=B/M. Then

    (1-s)x<=e<E,
    0<1-s<epsilon:=E/r_4<1.                     (5)

Put eta=epsilon/(1-epsilon). Since

    -log(1-u)=integral from 0 to u of dt/(1-t)
               < u/(1-u)  for 0<u<1,

(5) implies

    0<N log 2-h log 3=-log s<eta.

We have log 3>1. With alpha=log 2/log 3, this gives the necessary condition

    0<N alpha-h<eta.                            (6)

The exact certificate bounds eta between 905/10^36 and 906/10^36. There
is no floating-point decision in (5) or (6).

## 4. An elementary neighboring-fraction obstruction

Take the explicit integers

    a=253438635377865553122983825848491,
    b=401690733307859070263635245069941,
    c=326161216350563845551964374972875,
    d=516953297105243726377334761586514.

Exact rational logarithm bounds, described in section 5, establish

    a/b<alpha<c/d,
    cb-ad=1,
    b alpha-a>eta.                              (7)

In particular, b alpha-a lies between 1314/10^36 and 1315/10^36, safely
above eta's upper bound 906/10^36. These neighboring fractions can be found
by continued fractions, but no unverified decimal expansion or external
best-approximation theorem is needed for the argument below.

**Lemma.** If positive integers b,d and integers a,c satisfy (7), then any
integers N>0,h with 0<N alpha-h<eta satisfy N>=b+d.

**Proof.** The determinant condition makes (a,b),(c,d) an integer basis.
Thus uniquely

    (h,N)=u(a,b)+v(c,d), u,v integers.

Let e_-=b alpha-a>eta and e_+=d alpha-c<0. Then

    N alpha-h=u e_-+v e_+.

If v<=0, positivity of N forces u>=1, and the last expression is at least
e_->eta. If v>=1 but u<=0, it is negative. Both contradict the specified
small positive error. Therefore u,v>=1, whence N=ub+vd>=b+d. ∎

Now b+d is exactly the integer in (1). Combining sections 2-4 proves the
claimed universal zero-head result below that length. It applies to every
bank order and every guarded lift in the stated range, without enumerating
them. Equality at the length bound and longer words are not excluded by
this argument.

## 5. Rational logarithm certification

For 0<z<1 and integer m>=1, the atanh series gives the rigorous enclosure

    L_m(z)=2 sum over 0<=k<m of z^(2k+1)/(2k+1),
    L_m(z)<log((1+z)/(1-z))
        <L_m(z)+2z^(2m+1)/((2m+1)(1-z^2)).        (8)

The upper remainder follows by replacing every remaining denominator by
2m+1 and summing the geometric series. Set m=160. Rational z=1/3 bounds
log 2, and z=1/2 bounds log 3. If their intervals are [l_2,u_2] and
[l_3,u_3], respectively, then

    l_2/u_3<alpha<u_2/l_3.

The checker verifies the stronger rational inequalities

    a/b<l_2/u_3<u_2/l_3<c/d,
    b(l_2/u_3)-a>eta,
    l_3>1.

Only sums, products, divisions, integer powers of modest exponents, and
exact rational comparisons are used. No powers at the astronomical length
in (1) are computed. The proof is not dependent on the number of decimal
digits displayed, machine precision, or a guessed continued fraction.

## 6. The remaining obstacle

This resolves all shorter words symbolically and explains why modest
finite schedules were unable to challenge zero-head. It does **not** prove
the unrestricted conjecture. Near-resonant powers at larger lengths are
not ruled out, and no lower bound on their actual least source residues
has been proved here beyond the forced K^4 prefix.

A possible next route is to exploit the restricted block language to bound
its least positive 2-adic cylinder representatives in terms of length, then
compare that bound to the near-resonance threshold. The restricted language
has fewer parity words than unrestricted binary strings. That counting
observation alone does not imply large least representatives for each
word; neither random-residue heuristics nor a density argument supplies
the missing deterministic inequality.

The result remains a theorem about explicit sufficient cylinders. It gives
no arbitrary source entry, no forced reset capacity at arbitrary exits,
no global trajectory-frequency statement, and no proof of the Collatz
conjecture. Earlier frozen packets and their scopes are unchanged.

## 7. Reproduction

The [standalone read-only checker](../experiments/collatz_zero_head_length_obstruction_20261019.py)
performs 17 exact certificate checks. This small count reflects a symbolic
proof certificate, not a reduced word census. It verifies the K^4 source,
the small-word descent inequalities, the bank constant, the determinant,
the rational logarithm brackets, the decisive error comparison, the exact
length bound, and the reported rational enclosures. Normal and optimized
receipts agree byte for byte.

Run from the repository root:

    python3 05-knowledge/experiments/collatz_zero_head_length_obstruction_20261019.py > /tmp/collatz-zero-head-normal.json
    python3 -O 05-knowledge/experiments/collatz_zero_head_length_obstruction_20261019.py > /tmp/collatz-zero-head-optimized.json
    cmp /tmp/collatz-zero-head-normal.json /tmp/collatz-zero-head-optimized.json
    cmp /tmp/collatz-zero-head-normal.json 05-knowledge/results/collatz_zero_head_length_obstruction_20261019.json
    cmp /tmp/collatz-zero-head-optimized.json 05-knowledge/results/collatz_zero_head_length_obstruction_20261019.optimized.json

The checker prints JSON to stdout and performs no file writes or network
operations. It uses only Python's standard library.
