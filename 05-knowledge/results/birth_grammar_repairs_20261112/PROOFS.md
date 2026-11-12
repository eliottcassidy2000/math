# Synchronized births and arithmetic limits of endpoint repair

**PROVED scoped mechanisms; FINITE-EXACT bank classifications. Chan–Pak OPEN.**
Local labels below are packet labels, not newly reserved global theorem IDs.
The inherited definitions and tail construction are in
[the previous packet](../joint_stability_forbidden_20261112/README.md) and
[the full run-tensor proof](../session288_twin_lengths_20261112/stable/PROOF.md).
Here H is the number of directed Hamilton paths, delta=H-W(O), and U is the
selected-start Hamilton/two-path-cover count. A unit seed has delta=1 and
produces H=1+2^t U for every integer t>=1. A stable-clone birth is strong,
has a stable unselected port, and has no stable unselected twin contraction.

## B1. A finite synchronized grammar giving infinitely many births

The 24 explicit records in `synchrony_kernels.json.gz` have two originally
unstable singleton ports a,b and an independent stable singleton port z.
For their positive transitive lengths the complete signed tensor is

    delta(a,b,z)=1+epsilon*(2^a-2^b), epsilon in {1,-1}.

This is an identity of full skeleton tensors, computed in two independent
ways, not an inference from the length-two and length-three discovery tests.
Hence its unit locus is exactly a=b, with z unrestricted. Twelve records
have each sign. Swapping the coupled ports reverses the sign.

Set a=b=m and z=1. Both growing blocks are unstable: increasing only one
of their lengths changes delta by plus or minus 2^m. Every old nonsingleton
unselected block has, when one member is duplicated, defect minus one equal
to either

    4+4*4^m,  or  12-12*2^m+8*4^m.

The first is positive. The second is 8x^2-12x+12 at x=2^m, and is positive
because its discriminant is negative and its leading coefficient positive.
Thus no nonsingleton block is stable. Positive inflation does not merge
distinct maximal marked transitive blocks; this is the inherited weighted
twin-quotient theorem. The free singleton z remains stable. Strongness is
preserved by substituting nonempty modules in a strong quotient. Therefore
every m>=1 produces a strong stable-clone birth of order 8+2m.

For ordinary marked isomorphism, compress the maximal unselected twin blocks
and retain their lengths and all selected singleton marks. Each of these
24 quotients is rigid. There are 12 pairs sharing a quotient. Equality of
total orders forces the same m in both charts; the exact coordinate equations
force m=1. Thus the 24 charts give 12 old order-ten births, and **24 distinct
births at every even order >=12**, with generating functions

    B(z)=12*z^10+24*z^12/(1-z^2),
    B_new(z)=24*z^12/(1-z^2).

The word grammar is finite: choose one of 24 starts at m=2 and repeat the
operation G that adds a twin simultaneously to each coupled block. A further
operation F grows the free port and preserves the seed; birth status above
is asserted only on the z=1 slice. This is a genuine infinite birth grammar,
with fixed quotients, not a complete grammar of all births or all counts.
Its order-twelve states are disjoint from the 834 births obtained by the
declared one-vertex extension grammar. All-length separation and literal
order10,12,14,16 controls are in `synchrony_grammar.json.gz`.

The compressed defect signature is the earlier two-port synchronization
gate multiplied by the free module's H. In arbitrary modules A,B,C it is

    H(C)*(H(A)H(B)+C2(A)H(B)-H(A)C2(B))

after choosing the positive orientation; C2 is the ordered two-path-cover
count. For transitive modules, H=1 and C2(TT_m)=2^m-2. Thus the bridge from
the earlier Farey-type gate to these new births preserves the signed defect,
but discards the frozen-block geometry and U. Those are indispensable sidecars
for birth status and actual Hamilton counts respectively.

One explicit record has rows
`[32,5,33,7,15,474,159,287,607,255]`, selected vertices {1,3},
coupled ports 6,9 and free port4. Its diagonal slope is

    U(m,m,z)=12*4^m+Q(m)*2^z,
    Q(m)=-50+100*2^m+100*4^m+100*6^m.

Here v2(Q)=1. For z>=2m+2, v2(U)=2m+2. Writing z=2m+1+k,
k>=1, its normalized odd coefficient is

    c_k(m)=A(m)*2^k+3,
    A(m)=-25+50*2^m+50*4^m+50*6^m.

Thus an extremely simple synchronized birth rule creates an **unbounded
lower-head invoice j=2m+2**. Finite growth rules do not imply a bounded number
of unpaid exponents. This is the first obstruction the repair language must
address. The displayed coefficient family is not claimed fully repaired.

## V1. Complete valuation classification of the 478 old even higher-mode axes

Write U(a)=sum e_b*b^a in its unique finite exponential form. In this bank
the only odd bases are 1,3,5. Let g=min v2(e_b) over nonzero odd-mode
coefficients, and f(a)=sum_(b odd) (e_b/2^g)*b^a. Even modes vanish modulo
any prescribed 2-power once a is sufficiently large.

If f(a) is odd, v2(U(a))=g eventually. This gives 352 families. Otherwise,
the exact bank coefficients have odd coefficient of3 and even coefficient
of5. At a=2t+p, p=0 or1, write

    f(2t+p)=C1+C3*3^p*9^t+C5*5^p*25^t.

For every r>=0 and integer t>=0, the difference at t+2^r and t has exact
valuation r+3: v2(9^(2^r)-1)=r+3, whereas the 25 term has at least r+4
because C5 is even. Thus a parity branch either has fixed valuation1 or2
modulo8, or has one uniquely lifting zero at every successive precision.
The latter yields arbitrarily high v2(U), since one can choose an arbitrarily
large representative of each residue class, making all even modes vanish
at the requested precision. This proves the classification:

| Eventual behavior | Families |
|---|---:|
| Constant valuation |352|
| Two bounded parity branches |60|
| Unbounded valuation on a parity branch |66|

Explicit bounded cutoffs are at most5. There are 52,736 literal bounded
controls and compatible lifting controls through precision g+32, in addition
to the symbolic argument. Of the 66 unbounded families, 65 use odd bases1,3;
one also uses5. At birth519, port4, the full coefficients are
`[26,-860,1798,1698,240]`; its normalized odd part is
`13+899*3^a+120*5^a`. The even branch has the unique lifting zero.

These are classifications of a named finite collection of infinite axes.
They neither classify all future seeds nor supply their missing heads.

## V2. Why a finite bank of affine transitive templates cannot close 126 families

An **affine TT chart** is a fixed finite tournament template with some vertices
replaced by transitive modules of lengths alpha_i*m+beta_i, alpha_i positive
integers, with all lengths positive for sufficiently large integer m. Fixed
modules may be arbitrary finite tournaments, absorbed into the fixed exterior.
Different charts may have different templates, offsets and slopes.

**Claim.** If the H-count of such a chart grows at most exponentially in m,
its unique exponential expansion has 2-integral coefficients at every odd
base. Here 2-integral means rational with odd denominator after reduction.

If the growing ports contain a directed triangle, take t vertices in each
of its three growing blocks, with t proportional to m. Cycling through the
three blocks allows independent permutations in each, producing at least
(t!)^3 Hamilton paths on that induced subtournament. Hamilton counts cannot
decrease when another vertex is added: insert that vertex in a Hamilton path
at its first permissible position, recovering the old path by deletion.
Thus a cyclic set of growing ports gives superexponential growth.

Otherwise the growing ports are acyclic. The inherited finite run expansion
has integer coefficients in products of F_r(alpha_i*m+beta_i), where

    F_r(l)=sum_(b=1..r) (-1)^(r-b)*binom(r,b)*b^l.

After substitution, each exponential base is a product of b_i^alpha_i and
each coefficient an integer times a product of b_i^beta_i. For an odd base
all b_i are odd, so even negative offsets yield only odd denominators.
Collecting terms preserves 2-integrality. Distinct positive exponential
bases are linearly independent, by Vandermonde.

Now fix an infinite arithmetic progression on which v2(U(a))=j. Its lowest
unpaid head is `1+U(a)/2^(j-1)`. If j>=g+2, some odd-mode coefficient has
negative 2-adic valuation. Restricting to a progression a=Lm+r changes the
bases to b^L and multiplies coefficients by odd b^r, so the defect persists.
Adding the constant1 cannot remove a negative valuation at base1 either.
Consequently no affine TT chart equals this head throughout that progression.

Nor can a **finite union** of charts, with an arbitrary choice of chart for
each parameter, cover the whole progression. Cyclic charts eventually exceed
the target. For an acyclic chart a nonzero difference of finite sums of
distinct positive real exponentials has finitely many real zeros: divide
by the smallest exponential and apply Rolle's theorem inductively. Unless
the profiles are identical, only finitely many integer matches are possible.
The identical case has just been excluded by 2-integrality.

Each of the 60 parity families has a branch j=g+2. Each of the 66 unbounded
families has a residue class modulo4 with j=g+3. The receipt supplies these
126 exact branches and 8,064 independent literal valuation controls. Thus
**126 of the 478 families cannot be fully repaired by any finite bank of
affine TT charts**, even allowing parameter-dependent chart selection.

This does not exclude finite recursive rules that change the core, nonlinear
clock lengths, nontransitive growing modules, or a richer observer. The
obstruction is to this precise representation class, not to Chan–Pak and
not to finite-rule grammars in general. The other 352 families are outside
this eventual obstruction; that is not a proof they admit such repairs.

## V3. A different obstruction to periodically repeating old seam words

The inherited 65-letter seam grammar acts on a three-coordinate state by
integer 3x3 matrices. The exact determinant list has absolute values1 or
at least5, and none equal6. Hence the multiplicative monoid generated by
those determinants does not contain6: two nonunit factors have absolute
product at least25, while a single one cannot be6.

Suppose a fixed word W is repeated m times between any fixed prefix and
suffix, and its scalar output is a sum with all three nonzero modes1,2,3.
Its minimal scalar recurrence has characteristic polynomial
`(x-1)(x-2)(x-3)`. Cayley–Hamilton bounds it by the characteristic polynomial
of the 3x3 word matrix; equality of degrees forces determinant6, a contradiction.

Therefore this particular periodic repair architecture cannot generate a
generic three-mode head. This is independent of V2: even wholly integral
three-mode profiles meet it. Arbitrary nonperiodic words and higher-dimensional
recursive constructions are not covered. `repair_obstructions.py` checks
the pinned inherited atlas, rather than assuming its matrices unimodular.
