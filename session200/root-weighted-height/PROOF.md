# Root changes as additive carries and a linear-prime word budget

Session200 continuation. **PROVED from the physical O/Z/J grammar:** a
root-change-charged graph-order bound and an O(log n)-variable restricted
word circuit with a linear-in-n prime threshold. **FINITE-EXACT:** a hostile
to an over-restrictive two-length construction. **OPEN:** a constructor for
every odd target. This note does not replace the prior four-block lemma or
claim that its arithmetic range problem is solved.

## 1. The missing height ledger

For an actual source s=(D,E,U,d) let W=U-d. The physical cone gives W>=d>=1.
The exact source moves imply

    O: W'=2W+d, d'=d;
    Z: W'=2W,   d'=d;
    J: W'=W+d,  d'=D.

The remaining coordinates retain their full source update; they are not
functions of W alone. Thus O/Z are binary digit emissions with current
step-size d, while J is an additive carry injection that changes that
step-size to the previous path count D.

A block of k retained-root letters, with chronological bit a_i=1 for O
and0 for Z, has the exact physical address

    W'=2^k W+d A,
    A=sum_(i=0)^(k-1) a_i 2^(k-1-i), 0<=A<2^k.

Every A is realized by its literal word. This controls W, not the terminal
path count: D,E must still be updated by the same ordered word. Choosing W
and forgetting D is not an exact target construction.

Since W starts at1, each retained-root letter at least doubles W, and each
J increases W. A word with k retained-root letters and r J letters therefore
has W>=2^k. Under balanced OJ closure its count is

    n=D+2U=D+2W+2d >=2W+3.

For n>=5, writing j=floor(log2(n-3)), this proves

    k<=j-1, total word length<=j+r-1.             (1)

The closed graph has two seed vertices and two closure vertices, so its
order is at most j+r+3. In particular every represented four-block value
(r<=3) has an actual graph of order at most floor(log2(n-3))+6.

This is a bound on *every supplied word* with the stated root budget, not
only on a specially optimized representative. The singleton balanced value3
is a separate boundary. For k>=1, the root-zero words Z^k attain the graph-order bound:
their source is (2^k+1,1,2^k+1,1) and closed count3(2^k+1). The empty-word seed instead has count5
and also attains its bound separately.

## 2. A small modular constructor target with no hidden height condition

For n>=5 put j=floor(log2(n-3)) and L=j+2. Consider words of length at most L
with at most three J letters, padding by identity I. Their equations are the
ordinary Boolean one-hot source compiler: select I/O/Z/J, update the four
raw integer coordinates, and impose G=D+2U=n modulo a prime p.
The count of selected J letters is constrained to0,1,2,or3, for example
by the single degree-four equation c(c-1)(c-2)(c-3)=0 with
c=sum e_(i,J). Since p>L, this genuinely imposes the ordinary count bound.
The remaining compiler equations have degree at most two and O(log n)
variables. Alternatively keep a finite disjunction of four root budgets.

Every allowed word, regardless of the number of J letters, has

    U<=2^(L+1)=2^(j+3),
    0<G<=3*2^(j+3)<=24(n-3).                     (2)

Indeed every physical letter increases U by at most a factor2. Therefore,
for any prime p>24(n-3), a field solution decodes a literal word whose
ordinary count G agrees with n modulo p and has 0<G<p; it has G=n.
Conversely (1) supplies a word of length<=L for every root-three balanced
representation of n, and that integer word solves the system modulo every
such prime. No per-state binary range bits are required for this restricted
family. All raw states automatically belong to the same literal word.

Consequently the following is an explicit sufficient finishing lemma:

> **Logarithmic three-injection modular constructor (OPEN).** For every odd
> n>27,781, construct a word with at most L=floor(log2(n-3))+2 letters,
> at most three J letters, and G congruent to n modulo one prime
> p>24(n-3).

If established, (2) makes every constructed congruence exact, and the
previously verified finite repairs would yield the {7,21} coverage result.
The shared word, root budget, and prime scaling are part of the premise.

This is a sharpened constructor specification for the existing four-block
coverage problem, not an independent proof that its solutions exist. A
finite-field solver for individual n only decides the same bounded problem.
What remains useful is an all-n rule choosing its Boolean words.

## 3. A cheap hostile to an overly tight digit window

One tempting simplification is to use only the two adjacent near-maximal digit
budgets, k=floor(log2 n)-2 or floor(log2 n)-1. A complete target-capped
O/Z/J enumeration with r<=3 through50,001 refutes that restriction: it misses
1,498 odd targets in3<=n<=50,001 (including the separate singleton3),
with failures as large
as49,937. The unconstrained four-block class is already known to cover every
odd target above27,781 within this finite range. Thus one must retain shorter
retained-root words as well. For example,49,937 has the actual word
OOOZOZOJOOZOZJ, with profile(16593,6588,16672,6667),12 retained-root
letters and2 J letters. Its exact balanced count is16593+2*16672=49937.
No word with13 or14 retained-root letters and at most3 J realizes this
target. Thus a last-bit lift that blindly selects only those
two digit lengths cannot be the universal constructor.

`probe.py` enumerates triples (physical response, root count, retained-letter
count), not just scalar counts. It keeps all cost pairs because merging
only by response would lose the restricted-length predicate. U and G
strictly increase, so the target cap gives a complete finite universe.
`PROBE.json` records the641,253 state/cost triples and23,502 represented
odd values for the restricted two-length class. This census does not
establish infinitely many failures of that restriction.

## 4. Relation to the attached-paper lane

The dyadic saturation construction emits k about the modulus, whereas (1)
allows only logarithmically many retained-root emissions for an exact target
of comparable height. Each J is an additive injection, not a height reset.
Thus adding a few J letters after an exponentially tall residue witness
cannot shrink it. A successful constructor must choose the short word from
the outset, or replace the entire representation with a new short word.
That replacement must preserve the scalar congruence and all physical bits;
ordinary gcd division and independent CRT witnesses are known failures.
