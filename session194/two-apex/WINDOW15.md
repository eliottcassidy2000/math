# Twelve consecutive odd values at every 15-dyadic scale

**PROVED all-order construction from twelve finite exact seeds.** For every
integer k>=8 and every odd a with -13<=a<=9 there is an explicit strong
two-apex tournament with

    H=15*2^k+a.

Its order is k+4 or k+5. This is a fixed-width band of twelve values at an
infinite sequence of centers. Its width does not grow, and the result has
density zero. It does not establish the `{7,21}` spectrum, prime-atom
coverage, or exact arbitrary-target programming.

## Physical construction and exact continuation

Use the same physical two-apex notation as `PROOF.md`: a->b, increasing chain,
and letter c recording a's and b's outgoing incidences in bits 0 and 1.
For the seed prefix u selected below, use the word

    w=u 1^(k-q).

The bit convention means each new chain vertex lies after all old chain
vertices, is beaten by apex a, and beats apex b. The old tournament and both
apex marks are retained.

There is a particularly simple positive-tail identity, valid for every
physical prefix u and every integer t>=1:

    H(u 1^t)=A(u)+B(u)*2^t,
    A=2H(u1)-H(u11),   B=(H(u11)-H(u1))/2.           (1)

Here B can be zero in degenerate cases; our chosen B is positive.

A basis-free counting lemma proves more. Let (U,b) be any rooted tournament,
not necessarily two-apex. Add a transitive chain Z of t>=1 new vertices with
U minus b -> Z -> b. Write

    h=H(U),  e=number of Hamilton paths of U ending at b,
    D=H(U minus b),
    R=number of ordered spanning two-path covers (P_1,P_2) of U
      with P_1 possibly empty and P_2 nonempty and starting at b.

Then the enlarged tournament has

    H_t = (h-e-R+2D) + (R-D)*2^t.                  (2)

**Proof.** Every nonterminal block of Z vertices in a Hamilton path must
be immediately before b, because a Z vertex beats no old vertex other than
b. Hence there are at most two Z blocks: one immediately before b and one
terminal block.

With one Z block, a terminal block contributes h-e paths: an old Hamilton
path must not end at b. A block immediately before b contributes R paths,
by deleting the block and retaining the ordered two-path cover. The two
cases are disjoint and give h-e+R.

With two Z blocks, the first lies immediately before b and the second is
terminal. Deleting them leaves a cover counted by R whose second path starts
at b and has at least two vertices. Covers with that second path equal to
the singleton b are exactly D, so there are R-D eligible covers. Every one
permits 2^t-2 assignments of the ordered Z vertices to two nonempty increasing
blocks. Thus

    H_t=(h-e+R)+(2^t-2)(R-D),

which is (2). This is a bijection and requires only t>=1. QED.

For the two-apex words here, the named apex b is precisely the root: each
new 1-name is beaten by every old vertex except b and beats b. Therefore
(2) proves (1), with A=h-e-R+2D and B=R-D. In particular B is a nonnegative
integer. The verifier checks these rooted parameters independently for all
twelve prefixes as well as the two initial tail counts.

The t>=1 condition is genuine. Extrapolating (2) to t=0 gives h+D-e rather
than h. At k=7, the displayed fixed shifts for offsets -5,5,9 would use no
padding and instead give H=1855,1795,1825 rather than 1915,1925,1929. Their
deficits are exactly D-e=60,130,104. Thus k>=8 is the smallest common
threshold for this particular fixed seed prescription, not a claim about
the smallest k at which those values can be realized by other tournaments.

For a separate exact algebraic certificate, if f is the accepting row and
M_1 the explicit 36-state letter matrix, the verifier checks the full-row identity

    f M_1^3 - 3 f M_1^2 + 2 f M_1 = 0.

Right multiplying by any initial prefix state and any further power of M_1
proves the second-order recurrence from t=1 onward. Its roots 1,2 are
simple, and its two initial counts yield (1). This does not fit an
unproved exponential law to two samples.

## Complete finite bank

For each listed seed, the exact arithmetic is A=a and B=15*2^q.

| a | prefix u | q | order minus k |
|---:|---|---:|---:|
| -13 | 2222000000 | 7 | 5 |
| -11 | 2032 | 1 | 5 |
| -9 | 21300 | 3 | 4 |
| -7 | 323 | 0 | 5 |
| -5 | 1001133110 | 7 | 5 |
| -3 | 2020000010 | 7 | 5 |
| -1 | 230 | 0 | 5 |
| 1 | 110000022 | 6 | 5 |
| 3 | 31120 | 3 | 4 |
| 5 | 3000000012 | 7 | 5 |
| 7 | 130210 | 4 | 4 |
| 9 | 1000103210 | 7 | 5 |

Because every q<=7 and k>=8, t=k-q is positive. Substitution in (1) gives
H=a+15*2^q*2^(k-q)=15*2^k+a. The order is |u|+(k-q)+2,
which gives the final column.

Every prefix starts with a nonzero letter. Each also contains a 2, or
contains a 3 before a later 0. For example, 323 is covered by its 2 and
requires no 0. Hence the two apexes
lie in one strong component, that
component reaches the first chain vertex, and the last chain vertex (whose
letter is 1) returns to it through b. The increasing chain then belongs
entirely to that component. This proves strongness for every k>=8, not only
for the small checked instances.

## Why this companion is distinct, and its stopping boundary

The inherited S193 theorem supplies 135 consecutive offsets at centers 2^k.
This smaller bank supplies consecutive offsets at centers 15*2^k, which
lie far from any dyadic center for large k. Multiplying the old theorem by
15 through an ordered tournament join gives offsets spaced by thirty, so
it does not directly give the present consecutive-odd band. No optimality
or broad novelty claim is made.

The experiment searched equal normalized tail slopes rather than equal full
continuation responses. This is deliberately weaker than S193's decoder:
we require only the future language consisting of a positive run of 1s.
The bank retains its actual prefixes and the exact powers q, so the phase
alignment k-q is physical, not just a congruence. The fixed bank does not
show how to enlarge its width or carry an arbitrary target digit.

## Reproduction

    python3 verify_window15.py
    python3 -O verify_window15.py

The script imports only the local standard-library `verify.py` transfer
implementation. It checks the full-state row annihilator, both coefficient
readouts on all 36 coordinates, all twelve rooted (h,e,R,D) profiles and
exact coefficient pairs, 24 one/two-tail actual tournaments by independent subset
DP, all twelve physical witnesses at k=8 by the same independent DP,
six small seeds by literal path DFS, and strongness/order. The deterministic
receipt is `window15_receipt.json`; normal and optimized outputs agree.
