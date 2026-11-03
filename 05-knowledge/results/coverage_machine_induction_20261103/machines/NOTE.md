# A fixed integer reader with named tournament outputs and native returns

**PROVED** elementary constructions and bounds below; **FINITE-EXACT** small
controls in [normal.json](normal.json) and [hankel.normal.json](hankel.normal.json).
This encodes every nonnegative integer in actual marked tournaments. It does
not assert that every integer, or every allowed odd integer, is their raw
Hamilton-path count. Exact prime/7p spectrum coverage remains **OPEN**.

## Inheritance and the ordered-list comparison

The exact metagraph carrier is inherited from
[THM-796, three-sorted recursive incidence](../../../../01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md).
A Hamilton path gives an address chart with
`m=binom(n-1,2)` free arc bits. All-tile complement `a(t)=t+1` gives a native
line; reflection `r` gives converse after relabeling. The line is blue iff
`r(t)=t`, and black otherwise. The endpoint nodes are converse-merged
isomorphism classes. The complete state is the tiling/path plus named ports:
a node alone generally does not determine its next complement partner.

There are **two different roles for order**. FRACTRAN instruction order is a
priority rule: the first fraction giving an integer executes. Reordering
`(3/2,1/2)` changes the output on input2 from3 to1. The fixed Hamilton path
instead supplies addresses and a legal native move; changing to another path
of the same tournament preserves its raw H. The useful dictionary is:

| Machine datum | Tournament construction datum | Information that must survive |
|---|---|---|
| Ordered guarded instructions | Ordered branch selector | Failed earlier guards |
| Tape/counter addresses | Hamilton-path positions and physical names | Chosen chart and relabeling |
| Current configuration | Actual marked source | Entry program and attachment data |
| Executed instruction | Legal native/birth program | Selected path, ports, and intermediate states |

An ordered instruction list is therefore analogous to the *controller of*
the tournament construction, rather than being identified with the Hamilton
path itself. Both carry ordering information, but their preserved predicates
differ. An irreversible arithmetic step also is not an involutive metagraph
line. For reduced fraction `a_i/b_i`, retaining the selected index i lets one
invert a legal step by `N=N' b_i/a_i`; a history stack yields a reversible
lift. Native path selection similarly cannot be discarded when composing
line moves.

The inherited Session126 packet already proves that coordinate-labelled
free-bit interventions need the full tiling carrier, and that its small
ring CAs do not descend to merged nodes. Session132/133/135 supply earlier
machine compilers. [Session158's one-apex compiler](../../session158_creative_20261103/binary_bank/FRACTRAN.md)
specializes a FRACTRAN list to each word; the present reader instead uses one
fixed list for all inputs. The underlying one-apex count identity is inherited
from its [binary bank](../../session158_creative_20261103/binary_bank/NOTE.md).
MISTAKE-701's PRIMEGAME one-hot error is avoided: every rule here has a
control-state prime, including all fallback rules.

## An all-integer two-port code

Let a base tournament be transitive in the order `1,...,m`. Add apex0 with
`0->j` iff word bit `w_j=1`. Begin with `w=10`, the directed triangle.
Read the binary digits of N from least to most significant; for N=0 read the
empty word. For each read bit b, append the two bits `b0` to w. Write T_N for
the resulting actual tournament. More generally T_w denotes a source after
any finite input word w, including a word with terminal zero digits.

Each bit uses two legal negative-polarity births: the newborn loses to every
old base vertex, while its orientation to apex0 is selected by its bit.
Retain `(apex0,newest base vertex)` after each birth. The old terminal name
is explicitly discarded; apex0 persists. Thus every forgotten vertex has
the same orientation to each newborn. The two-port observer allowing both
birth polarities applies. This is an enlargement of the original
positive-polarity seventeen-coordinate language.

Every completed pair returns to a strong tournament. The base begins with a
vertex beaten by the apex and ends with a vertex beating the apex; the base
order then proves mutual reachability. The intermediate b=1 birth is a sink
extension and can be nonstrong. That intermediate is not declared a strong
return. The path `0,1,...,m` is always a Hamilton path.

**Marked injectivity.** Deleting the named apex leaves a transitive tournament,
whose total order is unique. Its apex incidence word recovers `10`, followed
by pairs `b0`. Reading the b's least significant first recovers N. This is
injectivity up to apex-preserving isomorphism; unmarked class injectivity is
not asserted. The second named port records the continuation endpoint.
The order is exactly `3+2 bitlength(N)`, with bitlength(0)=0.

Let P be the ordinary binary value of the entire word and h the current
completed-pair count. Initially `(P,h)=(2,3)`. Then

    P' = 4P+2b,
    h' = h+(3-b)P+2b.                              (1)

Proof: the inherited one-apex formula is2 plus the sum of `2^(i-j-1)` over
every earlier1/later0 pair, plus the number of adjacent01 pairs. Appending0
adds P; appending1 adds1 precisely after0. Thus00 adds `P+2P`, whereas10
adds `1+(2P+1)`. This proof concerns the accumulator at completed strong
returns, without falsely reading the advance credit as an intermediate H.

The outputs for N=0,...,7 are `3,9,27,31,99,121,97,117`. Every h is odd.
This gives actual prime sources31 and97 without claiming prime infinitude.

**Supply bound.** If N has positive bitlength l, before its final (necessarily1)
digit, `P>=2*4^(l-1)` and `h>=3`. Hence `H(T_N)>=4^l+5`. Consequently, for
X>=9, at most `sqrt(X-5)` inputs can give `H(T_N)<=X`: sum the `2^(l-1)`
positive inputs at each allowed length, together with input0. This is also
an upper bound on distinct outputs. The code faithfully stores every integer
while its raw-count image is sparse. The bound does not include native exits.

## One fixed FRACTRAN program and an exact CA/Turing interface

Use register primes `R=2, P=3, H=5, Q=7, S=11`. Their exponents respectively
store unread integer N, the prefix P, current h, quotient scratch, and prefix
scratch. Distinct control primes are recorded in [normal.json](normal.json).
The following single ordered list has25 fractions:

    17/26, 79/13, 38/17, 161/76, 53/38, 29/19, 19/23,
    62/203, 29/31, 37/29, 75035125/111, 37/41, 43/37,
    141/473, 43/47, 13/43, 118/371, 53/59, 61/53,
    24523675/183, 61/67, 71/61, 219/781, 71/73, 2925/71.

It starts at `13*2^N*3^2*5^3` and halts at `79*3^P*5^h`, where P and h
are the exact tournament registers after all digits of N. It is not Conway's
prime game and is not claimed to be a universal interpreter.

Proof by phases: the main state tests whether R>0, consuming and restoring
one unit before division; its fallback halts. At division, the two-unit
transfer `R-=2,Q+=1` precedes the one-unit branch and then the zero branch.
The latter two choose b=1 or0 respectively after exhausting R. Q is restored
to R, so R now equals its previous floor-half. In branch b, exhaust P into
S at rate4 while adding `3-b` units to H per old P unit. Restore S into P,
then add2b to both P and H. This is exactly(1), with zero scratch at return.
Every loop uses a separate helper state; no fraction cancels its state guard.
All inner loops exhaust a nonnegative counter; positive R strictly decreases
at each completed digit. This proves termination on every declared input.

A genuine radius-one CA reads the same word. Its alphabet is
`{blank,0,1} x {unmarked,marked}`. Initially place the finite least-significant-
first binary word in cells0,...,l-1 and mark cell0. The uniform rule is

    (data_i,marker_i)'=(data_(i+1),marker_i).

The marked cell outputs one input digit per step and then blank forever.
This CA is a finite-word reader, not a claimed universal CA; the CA continues
shifting after its output has ended. A fixed two-tape Turing transducer
initially writes10; state Read copies the current input bit to output and
enters Zero, which writes0 and advances the input head. Read on blank halts.
Both heads move only right or stay. The output is exactly the graph word.

The preserved predicate is the entire ordered birth word, with different
clocks: one CA step, two writing phases, and a variable finite FRACTRAN block
per digit. Raw microsteps are not identified. The compiler retains actual
graph rows, ports and births alongside the arithmetic state.

For external context only, Endrullis--Grabmayer--Hendriks,
[*Complexity of Fractran and Productivity*, v3, Section2](https://arxiv.org/html/0903.4366v3#S2),
provides the first-applicable definition and an explicit Turing-machine
translation. Its Theorem2.1 concerns uniform halting over *all positive
integers*, not merely the initialized states above. That theorem is not used
to infer spectrum coverage or to justify an unverified program translation.

## Native complements give a second actual return family

Let X_w be the native complement of T_w along `0,1,...,m`. The retained
physical ports are unchanged. A bit extension has a literal conjugated
program

    X_w --native complement--> T_w --birth b, birth0--> T_wb
        --native complement--> X_wb.                       (2)

Thus the enlarged alphabet has an actual closed source-program family.
After any input word containing1, X_w is strong: its base has forward arcs
only at adjacent indices and backward arcs at all other pairs; the directed
triangles on consecutive triples overlap and span the base. The apex beats
vertex1, and an input1 supplies a base vertex beating the apex. In particular
every canonical positive-integer completed source X_N is strong. X_empty
is transitive, and all-zero input prefixes can also be nonstrong.

Independent path counts give

| Input N | H(T_N) | H(X_N) | Order |
|---:|---:|---:|---:|
| 2 | 27 | 29 | 7 |
| 3 | 31 | 53 | 7 |
| 4 | 99 | 97 | 9 |
| 8 | 387 | 329=7*47 | 11 |

These are concrete prime/7p exits and reusable literal sources; they do not
constitute target selection for an arbitrary requested prime. Equation(1)
does not evaluate H(X_w).

The exit29 is a genuine source-family escape: the inherited binary-bank
membership proof excludes29 from *every* strong transitive-deletion source.
Thus X_2 cannot simply be another one-apex source in a different chart.

The native family still has a proved supply limit. Its base contains
`floor((2+2l)/3)=ceil(2l/3)` pairwise disjoint consecutive directed triangles,
where l is the input length. By
[THM-002, odd-cycle formula](../../../../01-canon/theorems/THM-002-ocf.md),
their subcollections alone contribute at least `3^ceil(2l/3)` paths. Hence
`H(X_N)>=3^(2l/3)`. For X>=1 the number of canonical integer inputs with
native output at most X, including input0, is at most

    X^alpha,  alpha=3 log(2)/(2 log(3))=0.9463946303...<1.

Indeed sum the `2^(l-1)` positive inputs through the largest permitted
length L and add input0, obtaining `2^L<=X^alpha`. Allowing every finite
bit word, including terminal zero digits, changes the bound by a factor
less than2. Thus even the native exits miss infinitely many prime and7p
targets, using the standard elementary prime-count lower bound. This bound
concerns these declared sources and canonical-path native exits; other
native paths, additional gluing operations, or a larger controller language
are not covered.

There is also an exact colour observable. In the fixed-path chart, backward
free arcs of T_N are apex incidences at zero bits. All but the final apex arc
reflect to forward base arcs. Therefore

    weight(t+r(t))=2(2 bitlength(N)-popcount(N)).             (3)

Input0 is blue; every positive input lies on a black line. This is a defined
machine-to-metagraph map preserving an explicit reflection-defect statistic,
not an identification of FRACTRAN instructions with metagraph nodes.

## Seven coordinates are already necessary for a linear native observer

Anchor the code at the input word1 so all subsequent T and X sources are
strong. For words u,v in `(empty,0,1,00,01,10,11)`, set
`K[u,v]=H(X_(1uv))`. The exact matrix is

    9     49    53     285    265    371    325
    49    285   265    1817   1553   2011   1659
    53    371   325    2533   2055   2639   2073
    285   1817  1553   12137  9785   12207  9557
    265   2011  1659   14225  11153  13997  10627
    371   2533  2055   17575  13805  16993  12925
    325   2639  2073   19365  14875  18485  13697

Its determinant is `-66127356928`. If this complete word-response family had
a d-dimensional linear realization `alpha A_w beta` over the rationals, K
would factor through that d-dimensional space, so d>=7. Thus an inherited
two- or four-coordinate linear count transfer cannot silently be reused for
these native exits. This does not exclude nonlinear or state-dependent encodings, nor prove
that seven dimensions suffice. The required all-order native observer and
an exact target-selector remain **OPEN**.

The31 graphs behind the matrix have orders5 through13. [hankel.py](hankel.py)
counts them both by subset endpoint DP and by a separately implemented
enumeration/packing of directed odd cycles, using the inherited odd-cycle
formula. This is a small complete certificate for the stated matrix, not an
expanded scalar-spectrum census.

## Reproduction and limits

Run `python3 check.py` and `python3 -O check.py`; outputs match normal.json.
These verify register/CA/Turing/marked-source data for N=0,...,31 and literal
big-integer FRACTRAN plus independent path DP for N=0,...,15 (order at most11).
Run `python3 hankel.py` and `python3 -O hankel.py`; outputs match
hankel.normal.json. All verification gates remain active under optimization.
The simulators cap individual runs at2,000,000 microsteps; the all-input
termination proof is not a claim that the capped script completes every N.

The nearest hostile is scalar compression: the distinct marked words1000
and1010 both have H9. Only the latter is a canonical positive-input code;
the former is a legal prefix/return in the construction language. The
seven-dimensional certificate tests the stronger continuation loss after
native transport. The constructive missing ingredient is target-directed
selection among actual sources, rather than the ability to encode arbitrary
input integers or execute a machine.
