# The five-count suffix cannot close the 6015/6017/6019 triple

**Finite certificate with an all-order reduction, scoped to one physical
repair.** The open 6015 endpoint for the known (1185,8051) sink key
cannot be obtained by a nontrivial ordinal suffix. The single-strong-core
6015 strategy remains open, as does the different key (729,12217).

[The phase note](README.md) shows that the only possible nontrivial suffix
for this key has H=5 and ordered two-slot count 32. In particular, proving
that the reduced sink key (H,g,oddpart X)=(1203,237,8051) has no physical
witness eliminates the suffix alternative, not just H=5 at a chosen order.

## Why the inputs are finite

1. Here gcd(1203,237)=3. No nontrivial ordinal prefix can divide 8051,
   and a suffix of H=3 has two-slot count 12 up to a power of two, again
   incompatible with the odd mantissa. Thus the full host can be taken
   strong by stripping ordinary singleton padding.
2. A split of g=237=3·79 into nonsingleton deletion SCCs would include
   a triangle. Its contact-content divisor 3 would divide X, excluded by
   the odd mantissa. So the deletion has one nontrivial SCC of H=237.
   Busch's lower bound allows at most eleven vertices in that SCC.
3. The exact weighted cut is c=(1203−237)/2=483. The first deletion SCC
   must contain an outcontact and the last an incontact in a strong host.
   If j≥1 singleton SCCs lie strictly between these endpoint components,
   fixing two suitable contacts and choosing any even subset of the j
   singletons makes 2^(j−1) distinct even paths of weight at least one.
   Hence j≤9, with at most two endpoint singletons. Eleven singleton
   slots therefore suffice at all orders.

## Enumeration

The exact filtered bank is grown by one-vertex extensions. All
intermediate counts are at most 237 by induced monotonicity. The
canonical key puts degrees in sorted order and takes the minimum
upper-triangle word across all equal-degree permutations, so retains
every isomorphism type. At each extension the cut calculation is exact.
The code also searches the excluded three-by-seventy-nine split,
which merely overcovers and is not needed.

For each layout of the components and singleton slots the program computes
R, the maximum number of two-slot choices, and rejects R<8051. A
single-contact start weight adds to X; their sum on the not-yet-assigned
vertices bounds how much X can still grow. Only X in {8051·2^k:k≥0}
can land. The cut contributions between already assigned vertices are
all nonnegative, so a partial cut exceeding 483 can be rejected.
Unassigned vertices are tried on both sides. These filters preserve every
possible witness. The path weights are assembled by parity convolution
along the prescribed ordinal component order; unvisited components
supply their full Hamilton count, and components strictly between the
two endpoint components of a selected path supply their parity-split
two-slot count.

Replay:
```sh
c++ -O3 -std=c++20 reduced_compiler.cpp -o /tmp/reduced237
/tmp/reduced237 11 11
```

The full compiler reports 84 strong deletion types of count 237 and
11 of count 79. It checks 14,492 layouts surviving the R screen.
No complete contact choice of cut at most 483 has the requested X.
All hypotheses above concern the stated sink key; the experiment does
not settle arbitrary common signatures for the old triple.

## Sharp minimum: the same signature first occurs at H=2701

The obstruction is stronger than the exclusion of 1203. Over all rooted
tournaments with g=237 and oddpart X=8051, the minimum H is exactly 2701.

The triangle-content argument still gives a unique nonsingleton deletion
SCC of H=237. If it lies outside the root SCC, it is on the ordinal
left or right. The left would force 237|X; on the right the odd part
of its P₂ would have to divide 8051. None of the 84 certified cores has
that property, as checked directly in check_decoder.py. This also
excludes the H=g root-singleton case. Thus discard exterior singleton
padding and keep at most thirteen singleton slots when seeking
H≤2701: then c≤1232 and 3+floor(log₂ c)≤13.

There is an independent inversion of the contact choices. For each
strong deletion block let R=P₂, and let L(A) be the sum of weighted
path starts in the chosen core contacts. Write P and Y for the binary
contact words to its left and right, and q for the length of the right
word. The two-slot count is
```
X = 237(Y+1) + 2^q(L(A) + R P).
```
This follows by whether the selected nonroot segment is empty or
starts on the right, in the core, or on the left. For each q and
candidate X=8051·2^k, invert 237 modulo 2^q to recover Y. Then each
choice of A forces P by exact division. The final retained right
singleton, if present, is an incontact, so Y is even and X odd;
otherwise the singleton would be ordinary padding. Leading zeroes
on the left are also ordinary padding. This inversion generates
68 addresses within the complete 13-slot bound. Direct parity-DP
and even-cut evaluation give minimum c=1232. The minimizer has
outgoing rows with the root at index 11:
```
2046 3580 4088 2544 3040 1472 384 256 1024 3562 216 1505
```
Literal subset DP gives H=2701, g=237, X=16102, verifying attainment.

Run `python3 check_decoder.py` in this directory with BANK_237.txt.
The decoder is independent of the C++ contact enumeration; it consumes
the literal canonical SCC bank associated with that generator.
