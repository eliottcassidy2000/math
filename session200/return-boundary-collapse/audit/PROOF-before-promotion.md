# Physical boundary returns erase joint-fiber branching after two blocks

Session 200 continuation. **PROOF CANDIDATE.** A tempting amplification of
Carmichael-style same-boundary collisions fails for a structural reason:
two nonempty returns to the triangle ray force the auxiliary boundary onto
the same ray. All later returns preserve the ratio of the two readouts.
Consequently this entire repeated-return construction misses every prime
target. The two-block threshold is sharp. Universal odd-target coverage
remains OPEN.

## 1. Objects and attempted existence construction

Use the physical source grammar and pairing from the synchronized inverse:

    t=(1,0,2,1), b=(3,1,3,1), z=(1,0,2,0),
    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D-E+U,d,2U-d,d),
    J(D,E,U,d)=(D-E+U,D-E,D+U,D),
    K(v,w)=v^T Qw,
    Q=[[0,0,0,1],[0,1,0,-1],[0,0,1,-1],[1,-1,-1,1]].

Put G(v)=D_v+2U_v=K(v,b). A synchronized boundary is `(w,q)=(b,z).beta`.
It realizes the joint readouts `(m,h)=(K(v,w),K(v,q))` on one actual source.
The Mersenne target is `n=m+L_k h`, where `L_k=2(2^k-1)` and k>=1.

The attempted construction is to use different words returning b to the
same scalar multiple, thereby holding m fixed while varying h. The smallest
collision is OJOJ versus ZZOJ; both send b to9b, but send z to different
vectors. The question is whether concatenating these choices thickens each
joint fiber. The answer below is exact and negative, even for arbitrary
nonempty return words, not just these two examples.

A **return word** is a nonempty chronological O/Z/J word R whose matrix
satisfies `M_R b=lambda_R b`. Since its output is physical,
`lambda_R` is a positive odd integer; growth excludes1, so lambda_R>=3.
The word is part of the data. This definition does not claim a source for
an arbitrary prescribed lambda.

Let L be the full balanced linear plane

    L={(x,y,x,y): x,y in R}.

This includes formal vectors; the physical profiles in L form a proper
subset. No formal vector in L is assumed to be a source.

## 2. Every return has one balanced output and one rank-one restriction

**Lemma 1.** Every return matrix M_R has

    M_R(R^4) subset L,       M_R(L) subset span(b).             (1)

**Proof.** First show R is globally matrix-equivalent to a chronological
word `P O J`. This is a matrix statement, not only equality at b.
The physical positive coordinates are

    (a,f,e,p)=(D-E-d,U-D+E-d,E,d-E)>=0,

and balance D=U is exactly f=p=0. A final O cannot produce balance. The
coordinate f after O is strictly positive on any nontrivial physical source,
whereas Z preserves f. Thus if R has a nonempty terminal Z-run and some
letter before it, the last non-Z letter must be J. The global identity
`J Z=O J` in chronological notation rewrites `J Z^k` to `O^k J`, giving a
suffix OJ. If the final letter is already J, balance says its predecessor
has E=0. Starting from b, only O can make E zero: Z sets E=d>=1 and J sets
E=D-E>=d>=1. Therefore that word also ends OJ.

The remaining possibility is R=Z^k. But

    b.Z^k=(2^(k+1)+1,1,2^(k+1)+1,1),

which is a multiple of b only for k=0, contrary to the nonempty requirement.
This proves the asserted factorization in every case.

On every raw vector,

    (D,E,U,d).OJ=(D+2U,D,D+2U,D),

so M_R maps the full vector space into L. Since `b.P.OJ=lambda b`, the
physical predecessor `b.P` has D=lambda and U=lambda. Thus `b.P` is balanced.

In positive coordinates b=(1,0,1,0), with both its a and e coordinates
strictly positive. Every source generator, and hence P, has a nonnegative
matrix in these coordinates. The f and p outputs of P on b vanish. Their
a and e coefficients must therefore each vanish, so P maps the entire
linear plane f=p=0 into itself. This is exactly P(L) subset L.
Finally OJ sends `(x,y,x,y)` to `x b`, so M_R maps L into span(b). QED.

The nonnegative-matrix step is essential: balance at one input would not
imply preservation of a linear plane for an arbitrary signed matrix.

## 3. Two returns freeze the normalized sidecar

**Theorem 2.** Let R_1,...,R_r be return words, r>=2, with multipliers
lambda_1,...,lambda_r. There is a nonnegative integer mu, depending only
on R_1 and R_2, such that the synchronized result is

    w=(product_(i=1)^r lambda_i)b,
    q=mu (product_(i=3)^r lambda_i)b.                         (2)

Consequently, for every actual nontrivial source v,

    h/m = mu/(lambda_1 lambda_2),                            (3)

independent of all later choices and of v. The physical joint cone further
gives `1/2<mu/(lambda_1 lambda_2)<1`.

**Proof.** The first return sends z into L by Lemma1. The second sends it
into span(b), say mu b. Its coordinates are integers, and its E coordinate
is mu; the extended positive-source cone containing z proves mu>=0. Every
later return multiplies b by its lambda. The physical boundary has been
multiplied at every step, giving (2). Pairing with v and using
`K(v,b)=G(v)>0` gives (3). This is an actual synchronized source boundary,
so the already proved joint cone gives the final strict inequalities. QED.

In particular, if a fixed finite menu of a return words has the same
multiplier lambda, then for each r>=2 it produces at most a^2 distinct
normalized sidecars, not an exponentially expanding family. The first two
choices determine everything except the common scale. Repetition alone
cannot turn a collision into an interval of joint realizations.

**Corollary 3 (prime obstruction).** For r>=2 and k>=1, every resulting
Mersenne target is the integer product

    n=G(v) (product_(i=3)^r lambda_i)
              [lambda_1 lambda_2 + L_k mu].                 (4)

Here G(v)>=5 and the bracket is at least9. Thus no prime target is realized
by this repeated-return scheme. The balanced readout m is also composite.
This is a restriction of this construction, not a nonexistence assertion
for tournaments or the unrestricted synchronized inverse.

## 4. Exact collision family and sharpness at one return

Set A=OJOJ, B=ZZOJ. Direct calculation gives

    A(b)=B(b)=9b,
    A(z)=5b, B(z)=(23,7,23,7).

After two choices the auxiliary states are precisely

    45b (AA or AB), 69b (BA), 71b (BB).

After r>=2 choices these are multiplied by9^(r-2), and no later choice
changes the three normalized slopes

    h/m in {5/9,23/27,71/81}.

Thus this natural same-boundary collision yields exactly three slopes
at every depth r>=2, not2^r independently adjustable values.

The two-return threshold cannot be improved. The single return C=ZOJ has

    C(b)=5b, C(z)=(11,3,11,3).

With v=t it gives m=25,h=19. Taking k=2 gives the prime target
`n=25+6*19=139`. The synchronized compiler yields the literal source
`t.JZOOZZ`, followed by its closing J. The resulting nine-vertex graph has
139 Hamilton paths and one Hamilton cycle, checked directly below.

## 5. Meaning for all-target coverage and verification

The inherited inverse correctly preserves the joint sidecar. This result
shows a further constraint on trying to populate it: repeated scalar-ray
returns erase that sidecar after two blocks. A successful prime-target
construction must avoid this two-return restriction, for example by keeping
a varying non-ray physical boundary. This is a concrete stopping reason for
the collision-amplification lane, not a general impossibility result.

`verify_return_collapse.py` checks exact return matrices for all words
through length ten, proves the indicated linear-plane inclusions by their
columns, checks the displayed family for all words in A,B through depth ten,
and literally counts the139 graph. These finite checks validate the code
and examples; Lemma1 and Theorem2 supply the all-order argument. No
all-target existence conclusion follows from the finite search.
