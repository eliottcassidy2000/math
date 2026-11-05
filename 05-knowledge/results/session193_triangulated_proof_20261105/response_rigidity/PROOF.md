# Lossless continuation decoding, radix obstructions, and exact quartic surgery

**PROVED by elementary arguments and finite exact invariant-space certificates.**
Session 193, 5 November 2026. Independent review is separate from this packet.

This combines S191's physical two-apex automaton with S192's common tensor
flag, then changes representation to two binary integers and to discrete
surgery derivatives. It does **not** prove the all-odd Hamilton-path spectrum
conjecture with exceptions 7 and 21, construct every required prime, or supply
an interval-filling theorem. The positive outputs are a lossless decoder and
an exact low-degree surgery compiler. The obstruction identifies a stronger
failure boundary for fixed radix constructions than independent-end commutation.

## 1. Inheritance and physical conventions

The closest proved mechanisms are [Session 191 two-apex grammar](../../session191_tournament_growth_20261104/two_apex/TWO_APEX.md) (the exact
36-state transfer and rational continuation dimension 17), and
[Session 192 tensor flags](../../session192_surprising_approaches_20261105/tensor/PROOF.md) (a common tensor flag, mixed-block spectra, and
Jordan counters). The canonical hostile is the physical word 120 with H=13:
it has no one-apex realization, and successful first ternary children do not
make a closed grammar. The corrected near miss is S191's `END_COMMUTATION.md`:
fixed left/right growth commutes, but two fixed **suffix blocks** generally do
not, so that obstruction alone does not exclude them.

The concept board is: full continuation response; binary incidence integers;
positive run allocations; mixed-block eigenvalues; noncommuting radix maps;
bounded-degree surgery. The sidecar retained throughout is the actual attachment
word, not its scalar H or its leading spectral term.

Fix a->b. The chain is x_1->...->x_m. A letter c_i in {0,1,2,3} records a->x_i
in bit 0 and b->x_i in bit 1. Let T(w) be this tournament and H(w) its number
of directed Hamiltonian vertex words. Empty w gives the edge a->b, with H=1.

For each apex order (p,q), a Hamilton word has runs L_0,p,L_1,q,L_2. Retain
states (x,y,z), with x,y=0,1,2 meaning empty, good last edge to the following
apex, bad last edge; z=0,1 records whether L_2 is empty. A letter with bits
P=[p->x_i], Q=[q->x_i] acts by the three allowed assignments:

    L_0: (x,y,z) -> (2 if P else 1,y,z);
    L_1: allowed if y!=0 or P=1, target (x,2 if Q else 1,z);
    L_2: allowed if z!=0 or Q=1, target (x,y,1).

Accept x!=2 and either y=1 or y=0 with p->q. These checks are precisely the
first and last boundary edges of the three runs. Thus this nonnegative
36-state transfer is exact by a run-assignment bijection.

Index states as 18o+6x+2y+z, with o=0 for (a,b), o=1 for (b,a). Write v(w)
for its state vector, M_c for a letter matrix, and f for the accepting row.
Then H(ws)=f M_s v(w), where M_s applies its letters in their reading order.

## 2. Two integer coordinates recover the entire physical word

Define the ordinary binary integers

    A(w)=sum_i (c_i mod 2) 2^(m-i),
    B(w)=sum_i floor(c_i/2) 2^(m-i).

These are physical coordinates, not guessed linear statistics. In the (a,b)
order, count allocations with L_1 nonempty and L_2 empty. If j is the first
name in L_1, its bit a->x_j must be 1; earlier names all lie in L_0 and each
later name independently chooses L_0 or L_1. Consequently

    A(w)=v_2+v_4+v_8+v_10+v_14+v_16.                  (1)

The same argument in the (b,a) order gives

    B(w)=v_20+v_22+v_26+v_28+v_32+v_34.              (2)

The unique allocation of every chain name to L_0 gives

    v_0+v_6+v_12=1.                                 (3)

In particular A(wc)=2A(w)+(c mod 2) and B(wc)=2B(w)+floor(c/2).

The new point is that these coordinates, and the constant in (3), are
observable from short **Hamilton-count continuations**. Abbreviate H_s=H(ws).
For every finite w, including the empty word, the following exact identities
hold:

    36 A = -5H_0 -3H_1 +6H_2 +15H_00 -9H_10 -18H_20
           +6H_30 +9H_01 -4H_000 +6H_100 +6H_200 -9H_010;             (4)

    12 B = -H_0 +3H_1 +3H_00 -3H_10 -3H_01 -2H_000 +3H_010;        (5)

    72 = 18H_empty +13H_0 -21H_1 -12H_2 -39H_00 +27H_10
         +36H_20 -30H_30 +45H_01 -18H_11 -36H_21 +18H_31
         +14H_000 +6H_100 -12H_200 -27H_010 +18H_210.               (6)

**All-word certification.** Each right side is an explicitly specified row
of the 36-state transfer. Starting at the two empty states, `verify.py` builds
a rational independent spanning set of reachable vectors and checks every
letter image of every added basis vector. The resulting 28-dimensional span
is therefore invariant, not a fixed-depth sample. On all 28 basis vectors,
the rows in (4), (5), (6) are checked to equal respectively 36 times (1),
12 times (2), and 72 times (3). The script also checks their two binary update
identities, the invariant constant, and the initial value. Linear equality
on the invariant reachable span proves (4)-(6) for all words. All arithmetic
is exact and the construction needs no imported matrix or witness bank.

### Theorem 1: short-response decoding and its exact ambiguity

The 17 counts with suffixes

    empty, 0,1,2, 00,01,10,11,20,21,30,31, 000,010,100,200,210

determine w up to an arbitrary number of leading 0 letters. In particular
this same ambiguity is the complete ambiguity of the full suffix response.

**Proof.** Obtain A,B from (4)-(5), write their binary expansions with equal
length max(bitlength A, bitlength B), and combine the two bits at each place
into its letter. This recovers w after removing its leading 0 letters. If
A=B=0, the canonical word is empty. Conversely, an initial 0 is a global
source in T(ws), for every suffix s; every Hamilton path must begin there.
Adding or deleting such source padding preserves every suffix count. QED.

This is a nonlinear two-unbounded-integer encoding, not a two-dimensional
rational linear representation. It does not contradict S191's exact minimum
linear dimension 17. The decoding theorem is for the named two-apex
presentation, not for arbitrary unlabelled tournaments.

### Corollary 2: affine rigidity of physical continuation responses

Let u,v be physical words and alpha,beta real numbers. If

    H(vs)=alpha H(us)+beta                            (7)

for the 17 listed suffixes, then alpha=1, beta=0, and u,v differ only by
leading 0 padding. The same conclusion holds if (7) is postulated for all
suffixes.

**Proof.** Every coefficient sum in (4), (5), and (6) is zero. Applying (6)
to (7) gives 72=72 alpha. Equations (4)-(5) then give A(v)=A(u), B(v)=B(u).
Theorem 1 gives identical canonical words and thus identical H. The empty
suffix in (7) forces beta=0. QED.

In particular no nonzero constant H correction can be hidden in an otherwise
unchanged future response. Nor can a nontrivial affine rescaling of the entire
response be carried by another physical word. This is not an obstruction to
nonlinear or state-dependent updates, or to a rule compatible only with a
restricted future language.

## 3. No fixed-suffix binary or ternary affine grammar

### Theorem 3

Fix q in {2,3} and distinct real offsets b_0,b_1. There are no two fixed
nonempty attachment blocks d_0,d_1 and nonempty word family F such that F is
closed under w->wd_i and, for every w in F,

    H(wd_i)=q H(w)+b_i,       i=0,1.                  (8)

The theorem allows arbitrary hidden state, arbitrary word lengths, nonstrong
words, arbitrary source/sink padding in the seed, and genuinely noncommuting
suffix blocks. No assumption of a small continuation quotient is made.

**Step 1: the common tensor flag forces single-letter branches.**
S192's explicit basis E,B,D=G-B in each bounded run and E,O in the final run
makes each letter matrix triangular. Its diagonal is the number j in
{0,1,2,3} of occupied level-one tensor factors, independently of the letter.
Thus a length-L block has only eigenvalues 0,1,2^L,3^L. This follows directly
from multiplying triangular matrices, not from multiplying separate spectra
without a common basis. On both orders its characteristic polynomial is

    X^8 (X-1)^16 (X-2^L)^10 (X-3^L)^2.              (9)

For either branch, distinct offsets and closure imply that some w in F has
H(w) unequal to that branch's fixed point -b_i/(q-1). Otherwise all H in F
would equal that fixed point, and the other branch would force equal offsets.
Repeating that branch gives a nonzero q^k channel. The exact finite-matrix
recurrence and (9) force q to be 2^L or 3^L. For q=2 or 3, L=1. Hence both
branches are single letters. Identical letters would force equal offsets.

**Step 2: q=3.** A constant nonzero letter c has at least one apex bit 1.
In either apex order, the run immediately preceding that apex cannot contain
any of the new c names, because its last new name would fail acceptance.
There are at most two remaining run choices per new name. For any fixed
prefix u, therefore H(uc^k)<=2*3^|u|*2^k. This cannot support a nonzero 3^k
channel. The only possible letter is 0; both branches would be identical.
This contradicts distinct offsets.

**Step 3: q=2, discard letter 3.** Appending 3 adds a global sink, so
H(u3)=H(u). Equation (8) would make H constant on F, impossible with distinct
offsets. The remaining distinct pairs are {0,1}, {0,2}, {1,2}.

We use a concrete lower bound for mixed tails. The construction fixes first
births and final endpoints of all three runs, then assigns each intervening
chain name freely among them. In each case below there are exactly 2k-4 free
names, so the lower bound is 3^(2k-4), for k>=2.

- **Tail (12)^k, arbitrary prefix u.** Put all prefix names in L_0. Put the
  first 1 in L_1, the following 2 in L_2, the penultimate 1 in L_1, and the
  final 2 in L_0. All other tail names choose freely. The first boundaries
  are legal, and the last 1/2 names give good edges to b/a respectively.
- **Tail (10)^k, if u contains a b-outgoing bit.** Put one such prefix name
  in L_2 and the other prefix names in L_0. Put the first tail 1 in L_1.
  Of the last three tail letters 0,1,0, put them in L_0,L_2,L_1. The other
  tail names choose freely.
- **Tail (20)^k, if u contains an a-outgoing bit.** Put one such prefix name
  in L_1 and the other prefix names in L_0. Put the first tail 2 in L_2.
  Put the last three letters 0,2,0 in L_0,L_2,L_1. The other tail names
  choose freely.

These are injections into actual Hamilton paths, not unsigned or asymptotic
transfer estimates. In contrast, repeated use of (8) for a two-letter tail
gives H=O(4^k). Since 3^(2k-4) is not O(4^k), the first case rules out pair
{1,2} entirely.

For pair {0,1}, the second case rules out any seed containing a b-outgoing
bit. If no such bit occurs, b is a global sink and every subsequent word
uses only 0,1. Appending 1 then inserts another sink immediately before b,
so H(u1)=H(u). This again contradicts distinct offsets.

For pair {0,2}, the third case rules out any seed containing an a-outgoing
bit. Otherwise every word in the closure uses only 0,2. Delete b: the
remaining transitive order is the chain followed by a. Thus this is a
one-apex presentation with binary word p0, where p records the b bits of u.
Let A be the integer represented by p, let t be its last bit (0 for empty),
and let h=H(u). The elementary one-apex append recurrence gives

    H(u0)=h+2A,
    H(u2)=h+A+2-t.                                  (10)

For completeness, if K=H_one_apex(p), then h=K+A+t; appending 0 produces
K+3A+t, while appending 2 produces H_one_apex(p10)=K+2A+2. These are (10).
Combining (10) with (8) forces

    h=2A-b_0,
    A+t=b_0-b_2+2.                                  (11)

Here b_0,b_2 denote the offsets attached to letters 0 and 2. But repeated
letter 2 sends A->2A+1 and t->1, so A+t is unbounded. It cannot satisfy the
fixed right side of (11). This excludes the final pair and proves the theorem.

**What survives.** The result does not exclude finite-state rules choosing
blocks based on the current state, variable/re-marked apexes, changing old
arcs, deleting or rearranging names, nonlinear scalar updates, or larger-radix
block grammars. For a general integer multiplier q>1, the spectral argument
only gives the necessary condition q in {2^L,3^L} for every fixed branch
length L. No universal impossibility for all such q is asserted.

## 3.1 Exact +2 cancellation is possible, but its hidden response grows

There is a useful positive construction at precisely the boundary the decoder
warns about. For every k>=1, both words below are strong, and

    H(213 0^k) = -11 +3*2^k +12*3^k,
    H(132 0^k) =  -9 +3*2^k +12*3^k.                 (14)

Thus changing the three-letter head gives an **exact +2** correction while
cancelling both exponential channels, for an unbounded family. This is not
only a congruence or asymptotic correction.

But after the same extra letter 1 their difference is

    H(132 0^k 1)-H(213 0^k 1) = 10-3*2^k.             (15)

For example k=1,2,3 give differences 4,-2,-14. The small correction is not a
reusable constant-response increment. Its two decoded integer differences
are already 3*2^k and -2*2^k respectively, so the lost information is explicit.

**All-order proof.** On the complete invariant reachable span each constant
letter satisfies M^4-6M^3+11M^2-6M=0, checked as a full-vector identity by the
verifier. Its distinct positive roots are 1,2,3; therefore any fixed-prefix,
fixed-suffix count along 0^k has the form A+B*2^k+C*3^k for k>=1. The three
exact initial values for 213 0^k are 31,109,337 and those for 132 0^k are
33,111,339, giving (14). The three exact differences after 1 are 4,-2,-14,
giving (15). The full-vector recurrence certifies the formulas at all k,
not a scalar fit. Strongness follows because both first letters are nonzero,
both words contain 2, and their final letter is 0. The initial five pairs
are also checked by vertex-subset DP.

This proves that spectral cancellation can produce exact microscopic
corrections. What fails is preservation of the correction under unrestricted
future growth. The quartic compiler below keeps the full response needed to
analyze a more elaborate correction bank.

## 4. A positive replacement: exact bounded-degree surgery

The following statement is not restricted to a transitive-chain complement.

### Theorem 4: marked-vertex surgery degree

Take a fixed vertex set and independent Boolean choices e_1,...,e_m. Each
edge orientation is either fixed or depends on exactly one choice. Suppose
all edges whose orientations can change are incident with a fixed set of
r marked vertices. Then H(e), as its unique multilinear polynomial on the
Boolean cube, has degree at most 2r.

**Proof.** For each permutation of the vertices, its Hamilton-path indicator
is the product of its consecutive directed-edge indicators. A path uses at
most 2r edges incident with the r marked vertices. Each changing edge
indicator is affine in one Boolean choice. Hence this product has degree at
most 2r. Reducing e_i^a to e_i for a>=1 cannot increase degree. Summing the
path indicators proves the result. This argument includes cancellation and
does not rely on a scalar extrapolation or positivity of derivative signs.

For two-apex words this applies to any bank of m disjoint fixed-length blocks,
each with two prescribed alternatives of the same length. Every choice
changes only apex-chain incidences on the same named vertex set. The same
degree bound holds for H(w(e)s) for **every** fixed suffix s, because its
varying edges are still incident with the same two apexes.

### Corollary 5: an exact quartic compiler for the entire response

Let R(e) be the 17-component response from Theorem 1. For S subset of the
switches define the vector

    C_S = sum_(J subset S) (-1)^(|S|-|J|) R(1_J).

Then for every active set E,

    R(1_E) = sum_(S subset E, |S|<=4) C_S.            (12)

All coefficients with |S|>=5 vanish exactly. Therefore evaluating only the
sum_(j=0)^4 binom(m,j) choices with at most four active switches determines
all 2^m physical continuation responses. Each evaluation consists of the
same 17 short suffix counts, or one exact transfer state followed by those
17 fixed readout rows. This is O(m^4) response evaluations, not a claim of
polynomial-time exact-target inversion.

Because the same identity holds for every suffix, it really preserves full
physical response. Theorem 1 also reconstructs the canonical physical word
from every compiled response. This avoids treating a quartic scalar formula
as a sufficient state for future constructions.

### Sharpness inside an entirely strong bank

Take the six-letter word

    w(e)=2, e_1, e_2, 2e_3, 2e_4, 0.

Every one of its 16 choices is strong: the initial 2 connects the apexes
through b->x_1->a, and the last 0 returns the entire increasing chain to
both apexes. Its mixed fourth difference is exactly

    sum_e (-1)^(4-sum e_i) H(w(e)) = 4.              (13)

The verifier independently counts all 16 tournaments by subset DP. Thus the
quartic order cannot be reduced to quadratic or cubic, even with strongness
uniformly enforced. Independent switches do not generally add: the fourth
interaction is a physically visible carry.

The compiler makes a proposed correction bank testable at all choices with
bounded-order interactions. It does not show that its integer image contains
an interval or that any coefficient can be selected independently. This is
the precise remaining target-control obligation.

## 5. Connection contracts and limits

- **Tensor flag -> radix obstruction:** preserve the full transfer and its
  spectrum under products; combine it with positive run injections. Eigenvalues
  alone eliminate branch lengths, not physical positive seeds. The latter
  are excluded by the explicit 3-choice path allocations and frozen-bit digit
  recurrence. This is stronger than literal commutation but still architecture
  specific.
- **Continuation observability -> physical reconstruction:** three explicit
  readouts recover two incidence integers and a normalizing constant. Linear
  dimension stays 17; the nonlinear decoder retains every attachment digit
  except exactly neutral source padding.
- **Local surgery -> Boolean polynomial:** preserve actual path indicators,
  and therefore all future responses, while replacing exponential choice
  enumeration by fourth-order interaction data. The missing theorem is an
  arithmetic statement about the polynomial image, not another leading-scale
  distribution or a generic eigenvalue computation.

The prime and 7-prime atom reduction discussed in the parallel session does
not alter these conclusions: neither the decoder nor the quartic compiler
proves realization of those infinite atom families. It identifies a rigorous
way to analyze candidate constructions without losing their continuation data.

## 6. Reproduction

Run from this directory:

    python3 verify.py
    python3 -O verify.py

The standard-library verifier rebuilds the transfer, checks 1,391 physical
Hamilton counts independently by subset DP (including the sharpness bank),
proves the decoder identities on the invariant reachable span, checks 21,845
canonical reconstructions, checks the positive run injections, verifies the
frozen-bit formulas, and reconstructs three full-response surgery cubes with
blocks of lengths one, two, and three. It writes a deterministic receipt;
normal and optimized outputs must match. The all-order proofs are above;
the finite probes are controls, not replacements for their quantifiers.

Prior source packets: S191 `grammar/TWO_APEX.md`, `grammar/END_COMMUTATION.md`,
and `arithmetic/PROOF.md`; S192 `algebra/PROOF.md`. The transfer and common flag
are inherited with attribution. The response decoder, its affine rigidity,
the binary/ternary fixed-suffix obstruction, and the quartic-response compiler
are the contributions of this packet. No literature-priority claim is made.
