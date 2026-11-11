# Translation by 30000: arithmetic predicates and actual rooted sources

**PROVED IN NOTE:** the translation-wheel identities, complete fixed-difference
square-prime classification, and a one-period H13 sidecar classification.
**FINITE-EXACT:** the displayed prime table, all 15 square divisor pairs, and
independent literal subset-DP checks on the source tournaments.
The new finite-certificate claims have internal independent-algorithm controls;
no literature-priority or independent-human-review claim is made.
**OPEN:** installation into a specified source, a target-preserving repeatable
return, Chan–Pak, and Landau's conjecture.

## 1. Inheritance and the live comparison

The closest inherited arithmetic mechanism is the prime-wheel transfer in
[Session51, prime-offsets/THEOREM.md](../../session51_structural_errors_relational_packets_20261024/prime-offsets/THEOREM.md).
Its hostile example `30030+19=151*199` separates wheel survival from primality.
The corrected near miss is the finite Fortunate/square-sum word correspondence:
it fails at its first extension, as recorded in
[Session52, fortunate-square-pattern/THEOREM.md](../../session52_saturation_residues_metagraph_20261024/fortunate-square-pattern/THEOREM.md).
The least-used relevant source sidecar is the full `(E,U,d)` response in
[PORTABLE_GAPS.md](../../h13_paired_descent_20261111/dynamics/PORTABLE_GAPS.md),
rather than only a scalar difference between two counts.

The board has five objects: a translated prime wheel; signed distance from
30030; square-difference factorization; actual pure-H paired sources; and
the modulo-13 phase of their right-word lengths. The first four are linked
below by explicit maps. Equality of the number 30000 is not itself a map
between prime certificates and tournament sources.

## 2. What addition by 30000 preserves

Put `M=30000=2^4*3*5^4`. For every integer x,

    gcd(x+M,M)=gcd(x,M).

For a prime ell, addition by M preserves every residue modulo ell exactly
when ell divides M. Hence the globally fixed prime channels are precisely
2,3,5. For any other prime, the translated number is rejected at
`x=-M mod ell`, which generally differs from the original rejection `x=0`.
This proves why the 2/3/5 wheel is preserved without asserting that primality
is preserved.

Let W be the product of the first 25 primes, through 97. Then

    W = 2305567963945518424753102147331756070,
    ord_additive(M mod W) = W/30
                         = 76852265464850614158436738244391869.

The exact number of residues x modulo W for which both x and x+M are units is

    product_(ell|W,ell|M)(ell-1)
      * product_(ell|W,ell not dividing M)(ell-2)
    = 117728574781173737675399914433625000.

The proof is the Chinese remainder theorem: a fixed channel forbids one
residue, and each other channel forbids two distinct residues. Relative to
the unit residues this fraction is `5132966455845/12094627905536`. This is an
exact finite residue fraction, not a density of prime pairs.

The first 25 primes and the extra value 107 give the complete table:

| x | prime rank of x | factorization of 30000+x | prime rank if prime |
|---:|---:|:---|---:|
|2|1|2*7*2143|—|
|3|2|3*73*137|—|
|5|3|5*17*353|—|
|7|4|37*811|—|
|11|5|30011|3246|
|13|6|30013|3247|
|17|7|13*2309|—|
|19|8|11*2729|—|
|23|9|7*4289|—|
|29|10|30029|3248|
|31|11|59*509|—|
|37|12|7^2*613|—|
|41|13|11*2731|—|
|43|14|13*2311|—|
|47|15|30047|3249|
|53|16|41*733|—|
|59|17|30059|3250|
|61|18|23*1307|—|
|67|19|107*281|—|
|71|20|30071|3251|
|73|21|17*29*61|—|
|79|22|7*4297|—|
|83|23|67*449|—|
|89|24|30089|3252|
|97|25|30097|3254|
|107|28|7*11*17*23|—|

Thus all three requested Berggren outputs become composite:

    17 -> 30017 = 13*2309,
    73 -> 30073 = 17*29*61,
    107 -> 30107 = 7*11*17*23.

Their prime ranks 7,21,28 do not implement the Berggren graph branch. The
inherited actual collision certificate is the three words `OIOIO`, `OOOI`,
`OIII`, all with H=9 but A-child counts 107,73,17; see
[Session256 spinor/PROOF.md](../../session256_chan_pak_triples_20261111/spinor/PROOF.md#5-the-scalar-h-obstruction-is-visible-at-h9).
It retains cut order, binary word and last-bit response. Forgetting that
information cannot be repaired by recording the rank of H among primes.

On the entire interval `0<=x<=107`, shifted primes occur precisely at

    x = 11,13,29,47,59,71,89,91,97,103.

The first-25-prime wheel also accepts x=49 and x=67, although
`30049=151*199` and `30067=107*281` are composite. The latter exposes the
requested number 107 as a missing sieve factor. Certification in this window
requires primes up to 173, since `sqrt(30107)<174`; stopping at 97 is not
enough. The composite original offset 91 becomes the prime 30091, explaining
the omitted prime rank between 30089 and 30097 in the displayed table.

## 3. The old 30013 example uses the nearby primorial center

The exact change of coordinates is

    30000+x = 30030+(x-30),
    30030 = 2*3*5*7*11*13.

Consequently `30013=30030-17` and `30047=30030+17`. These are the two first
prime endpoints at distance greater than one from 30030. The inherited
Session51 proof checks all smaller offsets and both prime certificates.
The primary treatment by
[Čejchan–Křížek–Somer, Section 3](https://cs.uwaterloo.ca/journals/JIS/VOL25/Krizek/krizek3.pdf)
distinguishes prime offsets near a primorial from a theorem generating prime
endpoints. No general Fortune-conjecture claim is imported here.

This recentering explains two of the new factorizations immediately:
`30017=30030-13` is divisible by 13, and `30107=30030+77` is divisible by
both 7 and 11. It does not predict whether the remaining cofactor is prime.
The signed offset is a required sidecar; 30000 and 30030 have different
prime support, despite being close as real numbers.

## 4. A complete, all-height A002496 intersection

Let S be the set of primes `n^2+1` with integer n. Then

    {x in S : x+30000 in S} = {547601}.

The unique pair is

    547601 = 740^2+1,
    577601 = 760^2+1,
    577601-547601 = 30000.

**Proof.** Choose nonnegative roots n<m. Subtraction gives

    (m-n)(m+n)=30000.

Conversely every positive divisor pair `(u,v)`, `uv=30000`, `u<=v`,
`u=v mod2` produces exactly one pair of nonnegative roots
`n=(v-u)/2`, `m=(v+u)/2`. There are exactly 15 such pairs. Because the
product is even, both u and v are even. Write u=2a,v=2b; then ab=7500.
For even roots, a and b have the same parity, hence both are even. Thus
`a=2c,b=2d` with `cd=1875=3*5^4`. There are precisely five factor pairs:

| (c,d) | (n,m) | n^2+1 | m^2+1 |
|:---|:---|:---|:---|
|(1,1875)|(3748,3752)|5*193*14557|5*13*216577|
|(3,625)|(1244,1256)|1547537, prime|13*121349|
|(5,375)|(740,760)|547601, prime|577601, prime|
|(15,125)|(220,280)|29*1669|78401, prime|
|(25,75)|(100,200)|73*137|13*17*181|

The other ten divisor pairs have odd n>1, so n^2+1 is an even integer
greater than 2. The exceptional prime n^2+1=2 would require n=1 and
`m^2=30001`, which is not a square. The prime assertions are certified by
trial division and independently by primes from an Eratosthenes sieve.
Every candidate has been included, so this is not an extrapolation from a
height cutoff. QED.

More generally, for any fixed positive M, simultaneous membership of
`x=n^2+1` and `x+M=m^2+1` reduces to the finite divisor pairs of M followed
by primality tests. This classifies fixed additive returns inside the
near-square-prime predicate. It gives no bound on the number of primes
`n^2+1` as M varies and does not prove Landau's conjecture. There are no
members of S in `[30000,30107]`: the adjacent candidate values are
`173^2+1=29930` and `174^2+1=30277`.

## 5. The same translation has actual tournament witnesses

Use the inherited H13 core `Q=(18,8,3,5,14)` with cuts 6 and 3, left
address A=3333, and right address B=0. Its two count formulas are

    H- = 13(g(A)+g(B)-1)+104AB+47A+9(A mod2)+32B+6(B mod2)+18,
    H+ = H-+3A-(A mod2)+2.

Here `g(0)=g(1)=1`, `g(2q+1)=g(q)`, and
`g(2q)=g(q)+q+(q mod2)`. We have g(3333)=1243, so H-=172837,
H+=182837. The common rooted sidecars are `(E,U,d)=(7,346677,13)`.

Prepend one directed triangle in ordinal order, retaining the old root.
Every full path visits this triangle before the old source, multiplying
H,E,d by 3. In U, the distinguished path starts at the old root and cannot
enter the triangle, so the triangle lies wholly in the other path and also
multiplies U by 3. Therefore the resulting actual 21-vertex tournaments have

    s- = (518511,21,1040031,39),
    s+ = (548511,21,1040031,39),
    s+-s- = (30000,0,0,0).

The adjacency rows are recorded in receipt.json. A newly implemented full
subset dynamic program counts every coordinate directly on both 21-vertex
graphs; it does not use the H13 formula or the ordinal multiplication rule.
This independently verifies these concrete source witnesses.

For the literal O/Z/J operations, the physical four-state maps are

    O(h,e,u,d)=(h,0,2u,d),
    Z(h,e,u,d)=(h+u-e,d,2u-d,d),
    J(h,e,u,d)=(h+u-e,h-e,h+u,h).

These linear maps send `(M,0,0,0)` to itself under O and Z, and to
`M(1,1,1,1)` under J. Thus the 30000 gap survives every common O/Z program;
after the first J its full response must be transported. For any later
common word W the H difference is `30000 H(W K1)`, with the rooted singleton
state `(1,1,1,1)`.

This is a literal source-compatible translation for the specified pair.
It is not a graph operation that adds 30000 to every source. Both baselines
here are divisible by 3 because of their ordinal triangle; neither is prime.
It also supplies no source with H=0: every tournament has a Hamiltonian
path. Arithmetic `0+x` is a coordinate origin, not a physical zero-count
tournament that can be substituted into this construction.

## 6. A complete one-period H13 sidecar classification

We now retain the sidecar rather than merely comparing formulas. For B>0,
put b=bitlength(B), r=rev_b(B). Both cuts have

    E=d=13,
    U=2^b(104A+45)-13r.

If two such sources have the same U, reduction modulo 13 gives
`2^b=2^b' mod13`. The order of 2 modulo 13 is exactly 12 (check the proper
divisors 1,2,3,4,6), so their right lengths differ by a multiple of 12.
For equal lengths, equality gives `8*2^b*(A-A')=r-r'`; since the latter
has absolute value less than 2^b, the addresses A,B are identical. Thus
same-length pure-H pairs in these two charts differ only by their cut.
Their nonzero width `3A-(A mod2)+2` is never divisible by 3.

For right lengths differing **exactly 12**, write the short-right source as
`(A,B)` and the long-right source as `(a,B')`. Equality of U gives

    A = 4096a+c,
    r' = r+2^b t,
    t = 14175-8c.

The bound `0<=r'<2^(b+12)` is equivalent to `0<=t<4096`, hence

    1260 <= c <= 1771.

Conversely every one of these 512 integers works. With `s=rev_12(t)`,
reversal of the concatenated bits gives exactly

    (A,B)=(4096a+c,B),
    (A',B')=(a,4096B+s),       a>=0, B>=1.

This proves that the parameterization is complete, including all actual
nonempty right words, with E,U,d preserved. The free parameters are actual
addresses, so no formal nonrealizable vector is being introduced.

For every one of these families and all four choices of the two cuts,

    H(long-right)-H(short-right) >= 49918.

This last number has a finite integer certificate. For each padded 12-bit
suffix z, the binary recurrence proves

    g(4096a+z)-g(a) = u_z a+v_z(a mod2)+w_z.

To compute these coefficients, read z from left to right with current
prefix p. A zero at position j (starting at zero) contributes `2^j` to
u, p to w, and either one to v at j=0 or `(p mod2)` to w at j>0. A one
contributes nothing. Then update p to twice p plus the bit.
Substitution into the displayed H13 formula yields

    H(long-right)-H(short-right)
      = K a+L B+P(a mod2)+Q(B mod2)+C.

All 2048 rows have K>=0 and L>=0. Thus each fixed-parity family is
nondecreasing and its minimum occurs among `a=0,1`, `B=1,2`. Enumerating
these explicit rows and four boundary corners gives minimum 49918.
The certificate is regenerated by `equal_sidecar_atlas()`; receipt.json
records its hash, coefficient ranges, all minimizers and 18432 direct
recurrence controls. No numerical search cutoff is used for a or B.

The minimum is attained at cut3/address `(1280,1)` and cut6/address
`(0,8111)`, whose full states are

    (213805,13,266317,13),
    (263723,13,266317,13).

Their literal 18- and 19-vertex adjacency graphs are independently checked
by the same subset DP. The small gap 30000 therefore cannot be realized
by same-length chart twins or by a single 12-bit sidecar-preserving length
change. **No exclusion is asserted for length differences 24,36,...,
other H13 charts, other cores, or other physical graph operations.**

This is a scoped source-fibre classification useful for the return search:
the loss of the modulo-13 right-length phase is now an explicit obstruction,
and the first possible length-changing phase has a known count cost.

## 7. Connection contract and reproduction

| Source | Target and map | Preserved predicate | Lost information / required sidecar |
|:---|:---|:---|:---|
|Integer x|x+30000|Residues modulo divisors of 30000|Other prime residues and prime certificates|
|Offset x from 30000|Signed offset x-30 from 30030|The actual endpoint integer|The center's factorization and offset sign must remain attached|
|Two square-plus-one values at distance30000|Factor pair `(m-n,m+n)`|Exact square condition, bijectively|Primality still needs its own certificates|
|Specified rooted H13 twins|Common O/Z continuation|Full difference `(30000,0,0,0)`|Absolute baseline and original source must remain attached|
|Equal-U H13 sources|12-bit address carry|E,U,d and literal source realizability|H changes by the certified affine response|

Run from this directory:

    python -B check.py --write
    python -B check.py
    python -O -B check.py

The first command regenerates receipt.json; later commands compare every
byte to recomputation. `g++` compiles the independent auditor in a temporary
directory. All validity gates remain active with optimization. Its arithmetic
has checked 64-bit additions and products; overflow aborts rather than
silently changing a count. The positive controls are K1, TT3 and C3, and a
two-vertex tie is rejected as an invalid tournament.

The arithmetic uses all first25 primes plus107, every x in0..107, and all
15 divisor pairs, with independent primality and inverse-square oracles.
The source audit uses six concrete large graphs and all 2048 one-period
address/cut families. The last two graphs independently check the base case
of [the growing-defect compensation obstruction](COMPENSATION.md): the known
actuator realizes the right gap but exceeds the target source's entire count,
so it cannot even occur there as an induced subtournament.
These scopes neither enumerate all tournaments nor
change the published values of A000568 or A135244. The A002496 gain is an
exact fixed-difference join, and the Chan–Pak gain is a realizable translation
pair plus a complete first-period sidecar obstruction.
