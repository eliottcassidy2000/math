# A first-crossing repair theorem for 24 interleaved bank episodes

Date: 19 October 2026.

## Result in plain language

The original-height repair can survive several returns to expansion. For a
lawful word built from the four banks below, stop at the first complete reset
block whose cumulative affine multiplier drops below one. If at most 24
expansion episodes have occurred, then either the actual source orbit falls
below its starting value within 92 further Terras steps, or the source lies on
a cycle and there is a different nontrivial positive cycle containing a smaller
integer. Consequently a least nonconverging positive integer cannot start such
a word.

There is no bound on the length of an individual episode or on the total word
length. The restriction is the number of expansion episodes before this first
crossing, after any exact bank rewrites described below. This is a restricted
language theorem; it supplies no pointwise theorem that an arbitrary Collatz
orbit reaches this language.

The argument has a symbolic carry bound and a small computer-assisted lemma
about 12,798 endpoint offsets. Those offsets are not starting integers. Each
offset certificate applies to every positive lift of a prescribed 120-bit
source cylinder.

## 1. Map, bank words, and lawfulness

Use the Terras map

    T(x) = x/2                       if x is even,
           (3x+1)/2                  if x is odd.

Parity words are chronological. A word of length n with h odd bits has affine
map (3^h x+b)/2^n on its unique lawful residue class modulo 2^n. “Fully guarded”
means the source follows every prescribed parity, not merely that a formal
endpoint is integral.

The lower and companion words are:

    K  = 110110110110110101011011101010
    K+ = 001111111111111101010110000100
    F  = 10                 F+ = 01
    H  = 110101000          H+ = 001101001
    V  = K K K K F          V+ = K+ K+ K+ K+ F+

K and V are expansion blocks; F and H are reset blocks. In this note a V is one
expansion block, of length 122 Terras steps. Its internal F is not a separate
reset boundary in that parse. A K has length 30, F length 2, and H length 9.

Write

    P = 1162261467 = 3^19,
    Q = 1073741824 = 2^30,
    C0 = 4912087879,
    D = P-Q = 88519643.

Then

    K(x) = (P x+C0)/Q,      x = 815464443 modulo Q;
    F(x) = (3x+1)/4,        x = 1 modulo 4;
    H(x) = (81x+101)/512,   x = 43 modulo 512.

For each block with branch (a x+b)/M, its companion has branch
(a x+b+M-a)/M, and its lawful source residue is exactly one higher. Thus a
lawful lower block carries an adjacent pair (x,x+1) to (y,y+1), using the
separately lawful companion word on x+1. Induction composes this identity
through arbitrary bank concatenations. It is an exact integer identity, not
an appeal to almost-everywhere coalescence.

Set

    C = 83829391903560839836419813590128166431
        /157489106280555890461927692448373459,
    E = C+1.

C is the expansion center of V, and exceeds the center C0/D of K. Numerically
C is about 532.286924; only the displayed rational is used. The elementary
bank checks, reproducible by the accompanying script, are:

1. Every nonempty prefix of K or V has affine multiplier strictly above one.
2. An expansion word has map A x+B with A>1 and B <= C(A-1).
3. A reset word has map r x+d with 0<r<1 and 0<=d<=1-r.
4. Each proper prefix of F or H has multiplier greater than its full reset
   multiplier and, at input 1, value at least the full reset value. Hence its
   value is strictly greater than the reset endpoint at every real input >1.

For item 2, compose maps whose individual centers are at most C. The carry
bound is preserved because
A2 C(A1-1)+C(A2-1)=C(A2 A1-1). Item 3 follows by composing reset maps that fix
or decrease 1. Item 4 includes the empty prefix of the reset block; the global
empty source prefix will be excluded whenever strict comparison with the
original source is asserted.

## 2. Precise first-crossing theorem

Fix a finite, fully guarded macroblock parse over K,V,F,H, beginning with K or
V. Let N be the first endpoint of an explicit F or H block at which the
cumulative affine multiplier is strictly less than one. Assume this endpoint
exists. Every earlier endpoint of an explicit F or H block has cumulative
multiplier at least one. This definition refers to each individual reset
block, rather than only to the endpoint of a whole consecutive reset run.

Truncate at N and group maximal consecutive runs as

    E1 R1 E2 R2 ... Em Rm,

where every Ei is a nonempty expansion word over K,V and every Ri is a
nonempty reset word over F,H. The last reset run ends exactly at the first
crossing just defined. The number m counts expansion episodes, not macroblocks
and not Terras steps.

**Theorem.** If m<=24 and x is a positive integer following this word, exactly
one of the following occurs:

(a) T^t(x)<x for some integer 1<=t<=N+92;

(b) T^N(x)=x; this lower orbit has least period N and minimum x; and a different
nontrivial positive cycle contains (x+1)/4, an integer strictly between 2 and x.

In particular, if a least positive integer that does not reach the cycle
1,2 exists, it cannot begin a word satisfying these hypotheses.

No statement here says that a later bank entry falls below the original
source before that entry. The source in the theorem is x at the start of this
first-crossing word. No statement replaces actual-orbit descent in (a) by a
descent on the adjacent companion orbit in (b).

## 3. An episode-dependent carry bound

Let Ei(x)=Ai x+Bi and Ri(x)=ri x+di. Put

    si = ri Ai,       ei = ri Bi+di,
    S0 = 1,           Si = s1 s2 ... si,       S = Sm < 1.

By definition of the first crossing, Si>=1 for i<m. The complete word has map
Sx+e, where

    e = sum from i=1 to m of ei S/Si.

The bank bounds give

    ei <= C si+1-(C+1)ri < C si+1.

Consequently

    e < C S sum from i=1 to m of 1/S(i-1)
        + S sum from i=1 to m of 1/Si
      <= C m S + (m-1)S + 1
      < m(C+1) = mE.

The last inequality also holds for m=1, since C>0 and S<1. This argument does
not require the individual episodes to contract; their multipliers may be
above, below, or equal to one, provided the stated cumulative first-crossing
condition holds.

If the endpoint does not descend, its integer displacement

    d = T^N(x)-x = (S-1)x+e

satisfies 0<=d<e<mE. Exact arithmetic gives

    12798 < 24E < 12799.

Hence for m<=24, every failed endpoint descent has 0<=d<=12798.

This bound is linear in the number of episodes. Section 8 explains why a
constant independent of m cannot replace it by a carry estimate alone.

## 4. Failure forces four initial K blocks

An expansion word either is K, K^2, or K^3, or its literal parity word starts
with K^4. Indeed, any V begins with K^4, and any word of at least four K blocks
does so as well.

For ell=1,2,3, the map of K^ell F has slope below one and sends every
x>=815464443 strictly below x. This follows from its exact affine identity:
check the slope and its value at the indicated lower endpoint. The reset H
is pointwise smaller than F for positive inputs, so K^ell H also descends.
The first reset after such a short initial expansion is already the first
multiplier crossing. Thus it is already the N of the theorem, and alternative
(a) holds.

In the remaining case the source must lie in the K^4 cylinder:

    x = r4 + u 2^120,    u an integer >=0,
    r4 = 588943162011782273059357097228697595.

The script reconstructs this residue from the full 120-bit word and also
checks it by direct parity replay. In particular x>466400.

## 5. Finite all-lift repair of every positive offset

**Offset lemma.** For every integer 1<=d<=12798 there is a prefix of the actual
orbit of r4+d, of length j(d)<=92, whose branch is (a z+b)/2^j with a<2^j and

    (a d+b)/(2^j-a) < 466400.

Its parity guard holds for every z=x+d with x=r4+u2^120. Therefore

    T^j(x+d) = (a(x+d)+b)/2^j < x

for every such lift.

This is a finite computer-assisted lemma with an exact certificate procedure.
For each offset, the verifier does the following:

    Start z=r4+d and (a,b,M)=(1,0,1).
    Read the parity of z and advance the genuine integer map T.
    For an odd bit replace (a,b) by (3a,3b+M); then double M.
    Stop at the first a<M, within the inherited 120-bit prefix.
    Reconstruct the entire branch and its lawful residue independently.
    Check j<=92, M divides 2^120, and r4+d matches that residue modulo M.
    Check a*d+b < 466400*(M-a), using integers only.

The full offset interval is checked without omissions. The maximum exact
threshold is

    810692053755261919/1738366812781 < 466400,

attained at d=11808. The largest repair time is 92, attained only at d=4.
The first-contracting parity cylinders compress to 718 disjoint offset
classes, and the verifier separately checks their exhaustive coverage.
The receipt records the time histogram and a deterministic digest of all
12,798 exact rows, avoiding a large duplicated table in this note.

The direct replay at one extra lift is a control, not the justification for
all lifts. That justification is the modulus divisibility and the displayed
uniform affine inequality.

If 1<=d<=12798 in Section 3, apply the lemma at time N. This gives actual
source descent by time N+j(d)<=N+92.

## 6. Zero displacement and the smaller cycle

Suppose d=0, so T^N(x)=x. No nonempty proper lower prefix can have value at or
below x. To see this, consider all macroblock endpoints before N. They have
cumulative multiplier at least one. The cumulative carry is strictly positive
because the word begins with an odd bit. Hence their values are strictly
above x. A prefix within K or V is larger than that block's input by item 1
of Section 1, including prefixes inside the internal reset of V. A proper
prefix within an explicit F or H is strictly above that reset's endpoint
by item 4; all such endpoints are above x except the final one, which is x.
Every reset input is greater than 1. Thus the comparison covers every time
1<=t<N, not the empty prefix at time zero.

It follows that the lower cycle has minimum x and least period N. Separately,
the lawful companion identity gives T^N(x+1)=x+1. Since the forced lower
prefix is K^4, the companion starts with two even steps. Its cycle therefore
contains

    z = (x+1)/4.

This is an integer with 2<z<x. The companion cycle is not the cycle 1,2 because
it contains the periodic point x+1>2. It is not the lower cycle because the
latter has minimum x. This proves alternative (b).

The two alternatives are mutually exclusive: in (b), the lower orbit cycles
with minimum x forever. If x were the least positive nonconverging integer,
(a) would provide a smaller nonconverging iterate on its own orbit, whereas
(b) would provide the smaller nonconverging periodic integer z. Either is a
contradiction. The least-source qualification is essential.

## 7. Exact normalization can reduce the episode count

The rewrite

    K K K K F  ->  V

changes neither the lower parity word nor its companion word. It therefore
changes neither the actual orbit, total Terras length, guard conditions, nor
any affine coefficient. It can be applied before counting episodes.

Before the first crossing, a replaced segment starts with cumulative slope
at least one and has multiplier equal to that of V, which exceeds one. Its
F endpoint cannot have been the first contracting reset boundary. Also,
every prefix of V expands, so no intervening hidden contraction is introduced.
Thus the first-crossing time N is preserved. The rewrite never increases the
number of expansion episodes; if the absorbed F separated two expansion
runs, it merges them and reduces the count.

The theorem therefore applies whenever this normalization produces at most
24 episodes, even if the original syntactic parse had more. It also applies
to any other explicitly supplied parse into these same exact blocks that
satisfies the first-crossing hypotheses. No arbitrary approximate block
replacement or unguarded branch substitution is allowed.

## 8. Why arbitrary interleaving needs another ingredient

There are two distinct limitations.

### 8.1 A concrete hole in this inherited-prefix certificate

For m=25 the carry budget admits offsets up to floor(25E)=13332. At the
specific offset d=12852, none of the first 120 prefixes of r4+d has affine
multiplier below one. The verifier checks this exactly. The same is therefore
true for the inherited 120-bit parity prefix of every x+d in that source
cylinder. Positive carries and d>0 prevent descent below x along those
noncontracting prefixes.

This does not say that d=12852 is realizable as an endpoint displacement of a
25-episode word. It does not say that it can never be repaired. It says only
that the present carry bound plus this inherited prefix cannot supply the
same universal repair lemma. Later parity bits require more source
information; they cannot be extrapolated from a 120-bit cylinder.

### 8.2 No episode-independent normalized-carry bound

There is also a structural obstruction to making Section 3 uniform in m.
Put

    lambda = P/Q > 1,       mu = 3/4,
    c = C0/(P-Q),
    theta = log(lambda)/log(4/3).

The number theta is irrational: if theta=u/v with positive integers u,v,
then lambda^v=(4/3)^u would imply
3^(19v+u)=2^(30v+2u), impossible by unique prime factorization.
For integers L tending to infinity along a sequence with fractional parts
{L theta} tending to zero from above, let

    jL = floor(L theta),
    UL = K^L F^jL,
    sL = lambda^L mu^jL > 1.

Such a sequence exists by density of irrational rotations; equivalently it
follows from the elementary one-sided rational approximation lemma. Then
sL tends to 1 from above and jL tends to infinity. The exact affine carry of
UL is

    eL = c sL + 1 - (c+1)mu^jL,

so eL tends to c+1. Given any positive episode count m, choose L sufficiently
far along this sequence that eL>c and sL^m<4/3. The bank word

    UL^m F

has its first contracting explicit reset endpoint at the very last F. Inside
each UL, reset multipliers decrease to sL>1; the preceding copies multiply
this by further factors sL. Its normalized carry is

    mu eL (1+sL+...+sL^(m-1)) + 1/4
      > (3c/4)m + 1/4.

Thus normalized carries are unbounded as the number of episodes grows, even
within this same guarded bank language and even when the entire word is
stopped at its first multiplier crossing. Every finite word has positive
lawful integer lifts, so these are legitimate cylinders, not formal
unrealizable branches. This proves no assertion of actual non-descent on
those lifts: a large carry alone cannot overcome an arbitrarily large source.

A modest exact control already shows why the one-episode bound cannot simply
be reused: take L=11, jL=3, m=24. The word

    (K^11 F^3)^24 F

has 8066 Terras steps, 5089 odd bits, total multiplier about 0.9164741983, and
normalized carry about 657.7024741, strictly greater than 534. All reset
boundaries before its last F have cumulative multiplier above one. The
checker verifies the inequalities with exact rationals and reconstructs the
full word independently. The decimals are only for readability.

## 9. Scope, lineage, and reproducibility

This extends the single-expansion/single-reset carry and offset argument in
[the original-height repair note](../collatz_original_height_repair_20261019/ORIGINAL_HEIGHT_REPAIR.md).
It does not amend that result or the separate
[height-controlled density estimate](../collatz_height_controlled_entry_20261019/HEIGHT_CONTROLLED_ENTRY.md).

The central unresolved step remains coverage: no proof is given that every
positive source, or a putative least nonconverging source, has an admissible
first-crossing bank prefix with a bounded episode count. No finite census or
density-one occurrence statement closes that pointwise gap.

From this directory run:

    python verify_interleaved_height_repair.py
    python -O verify_interleaved_height_repair.py

Both must reproduce verification-receipt.json. The exact controls support the
bank interfaces, forced cylinder, all 12,798 offset inequalities, hostile
12852 prefix, and constructive large-carry example. The symbolic proof is
what permits arbitrary bank lengths, arbitrary lawful lifts, and the
unbounded family in Section 8.2.
