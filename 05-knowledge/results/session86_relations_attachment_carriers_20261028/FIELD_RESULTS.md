# One marked line reconstructs the binary pairing in F64 and F256

2026-10-28. **PROVED classification; FINITE-EXACT marked-field tests.**
This extends source P3 and answers source programme O2 at m=3 and m=4.
It does not claim a new classical trace-form construction or an all-m
one-line theorem.

## Result first

Let E=F_(q^2), K=F_q, q=2^m. Keep the actual marked field and its absolute
Frobenius F:x->x^2. Among nondegenerate alternating binary forms, exactly
**one marked K-line other than K itself** suffices, in each tested field,
to recover

    B_1(x,y)=Tr_(E/F2)(x y^q).

Zero marked lines do not suffice. In fact **every** non-K line suffices in
these two fields. The K-line K itself supplies no restriction at all on
Frobenius-invariant alternating forms.

| Field | Polynomial for alpha | All F-invariant alternating forms | Nondegenerate ones | Minimum marked K-lines | Unique trace mask |
|---|---|---:|---:|---:|---:|
| F64 | alpha^6+alpha+1=0 | 8 | 2 | 1 | 27760 |
| F256 | alpha^8+alpha^4+alpha^3+alpha^2+1=0 | 16 | 8 | 1 | 124407440 |

Integers encode polynomial vectors least-significant coefficient first.
A form's mask lists its upper-triangle entries in lexicographic pair order
(0,1),(0,2),...,(n-2,n-1). These masks are marking-dependent, not invariants.
Both displayed polynomials are checked irreducible and primitive by exact
polynomial gcds and multiplicative-order prime-factor tests.

Explicit sufficient marked lines are alpha K, with nonzero vectors

* F64: `{2,28,30,44,46,48,50}`;
* F256: `{2,20,22,45,47,57,59,136,138,156,158,165,167,177,179}`.

A marked line is not a single bit: specifying total isotropy requires the
vanishing of all binom(m,2) basis-pair entries (three at m=3, six at m=4).
On the invariant form space their independent ranks are 2 for six F64
non-K lines, 1 for its two exceptional lines, and 3 for every F256 non-K
line; K has rank 0. Before the nondegeneracy restriction, these remove
exactly those many binary degrees of freedom. The proof of the minimum
is exact: enumerate the entire invariant linear space, discard singular
forms, and intersect its survivor set with each line's isotropy equations. All minimum singleton choices and all survivor
sets are retained in `results.json`.

## 1. Re-derivation of the all-size spread classification

Suppose C is alternating and vanishes on every K-line. From
C(a(x+y),x+y)=0 for a in K, the same-line terms vanish and give
C(ax,y)=C(x,ay). On a K-basis e,f write ell(t)=C(e,tf). Then

    C(ae+bf,ce+df)=ell(ad+bc).

Conversely any binary linear ell:K->F2 gives such an alternating form.
Consequently the space of all spread-compatible forms has dimension m and
exactly q elements, **including zero**. Put

    delta=e f^q+e^q f in K^*.

The forms

    B_c(x,y)=Tr_(E/F2)(c x y^q),  c in K,

have coordinate expression Tr_(K/F2)(c delta(ad+bc)), so run through the
whole space. Trace-pairing nondegeneracy makes every c!=0 nondegenerate.
Thus the number of nondegenerate forms is q-1, not the size of a linear
space. Frobenius transports B_c to B_(sqrt(c)); invariance forces c in F2.
After excluding zero, c=1 is unique. This independently reproduces P3.

## 2. All Frobenius-invariant alternating forms, without spread constraints

Every binary bilinear form can uniquely be written

    C(x,y)=Tr_(E/F2)(x L(y)),
    L(y)=sum_(i=0)^(n-1) a_i y^(2^i),  a_i in E, n=2m.

Indeed the trace pairing identifies linear maps with bilinear forms, and
linearized polynomials of degree below 2^n give all linear maps: injectivity
follows from the polynomial root bound and both spaces have binary
dimension n^2. Frobenius invariance is equivalent to L(y^2)=L(y)^2,
forcing every a_i in F2 by uniqueness. Symmetry forces a_i=a_(n-i).
In characteristic two, its diagonal is then a_0 Tr(x); the paired terms
cancel and the middle term Tr(x^(q+1)) vanishes. Alternation therefore
means a_0=0. The complete classification is

    L=sum_(i=1)^(m-1) a_i(F^i+F^(n-i)) + a_m F^m,
    a_1,...,a_m in F2.                                      (1)

This gives dimension m, independently of the finite enumerator.
The Frobenius operator has minimal polynomial t^n-1: no nonzero polynomial
in F of degree less than n can vanish, again by the linearized-polynomial
root bound. Thus L is invertible iff its coefficient polynomial is coprime
to t^n-1.

* For n=6, t^6-1=(t+1)^2(t^2+t+1)^2. Evaluating (1) at 1 and at a primitive
  cube root shows nondegeneracy iff a_3=1 and a_1=a_2. Exactly two survive.
* For n=8, t^8-1=(t+1)^8. Evaluating at 1 shows nondegeneracy iff a_4=1.
  Exactly eight survive.

These a_i are linearized-trace coefficients; they should not be confused
with the separately recorded circulant coefficients of a chosen normal
basis. The program independently checks the latter representation too.

There is also a structural reason K never helps. Since
K=ker(1+F^m)=im(1+F^m), Frobenius invariance implies

    C((1+F^m)x,(1+F^m)y)=0.

The first and last terms cancel, as do the two middle terms. Hence every
Frobenius-invariant bilinear form in this discussion vanishes on K.

## 3. Exact line constraints, including degenerate controls

For F64, the two nondegenerate masks are 656 and 27760. Form 656 is
isotropic only on K; form 27760 is isotropic on every spread line. Of the
eight non-K lines, six cut the entire invariant space to dimension 1.
The other two leave dimension 2, with both additional nonzero survivors
singular: masks 9568 and 18704 have ranks 4 and 2 respectively. These exceptional lines are

    {3,17,18,40,43,57,58},
    {5,8,13,51,54,59,62}.

They form a Frobenius orbit of length 2. It is important that the question
requires **nondegenerate** forms: these two lines do not isolate the trace
form from all nonzero alternating forms.

For F256, the eight nondegenerate masks, in ascending order, are

    29009424, 66200080, 86166160, 124407440,
    161571216, 199878032, 221842704, 259098896.

Every non-K line cuts the entire invariant space from dimension 4 to 1.
Thus each leaves only zero and the trace form even before nondegeneracy is
imposed. The two nonzero Frobenius orbits of spread lines both have length 8.

A useful algebraic control is

    C_alt(x,y)=Tr(xy)+Tr(x)Tr(y).                            (2)

For even n it is alternating, Frobenius-invariant and nondegenerate:
its trace-pairing linear map is y->y+Tr(y), which has trivial kernel since
Tr(1)=0. On a line tK, put s=t+t^q. The restriction becomes

    Tr_K(s^2 uv)+Tr_K(su)Tr_K(sv).

It vanishes for t in K. Otherwise s!=0 and this is equivalent, after
rescaling, to Tr_K(ab)+Tr_K(a)Tr_K(b). For m>=2 that form is nonzero:
its trace-pairing map b->b+Tr_K(b) has kernel of dimension at most 1,
strictly less than m. Therefore (2) vanishes on no other K-line.
At n=6 it is exactly the sole alternative mask 656; this proves the
one-non-K-line assertion algebraically, not just computationally.

## 4. Hostile marginals and an exact attachment condition

All q-1 nonzero B_c have **identical zero restrictions on each separate
K-line**. Even retaining the unlabelled bit histogram for every transverse
pair of lines does not distinguish them. If L and M are distinct K-lines,
then the restriction B_c:L x M->F2 is a nondegenerate binary pairing. For
uniform independent x in L,y in M, it has exactly

    q(q-1)/2 ones among q^2 marked pairs.

Thus these histograms agree, although the labelled joint tables differ.
This is a precise marginal-data model, not a claim that all possible
coordinate marginals are equal.

| Field | Hostile coefficient c | Hostile mask | B_1(1,4) | B_c(1,4) | Ones on every transverse L x M |
|---|---:|---:|---:|---:|---:|
| F64 | 14 | 6062 | 0 | 1 | 28 of 64 |
| F256 | 10 | 232908522 | 0 | 1 | 120 of 256 |

These hostile B_c are spread-compatible but **not Frobenius-invariant**.
That is exactly the constraint which removes their ambiguity. Across all
ordered pairs of E, the two forms differ at 2016 or 32640 pairs respectively,
because their difference B_(1+c) is nondegenerate. The program checks each
pair, not merely the formula.

For two transverse marked ports eK and fK, consider K-linear attachment
rescalings e->ae, f->bf with a,b in K^*. Their mixed pairing changes by the
coefficient ab. Consequently the exact compatibility condition is

    ab=1.                                                   (3)

The individual same-port restrictions remain zero regardless. More
generally a K-linear attachment A preserves B_1 iff det_K(A)=1, since
B_1(Ax,Ay)=B_(det A)(x,y). Independent per-port identifications cannot be
silently replaced by a compatible joint identification. The missing data
are the actual marked attachment maps, not another copy of the same
marginal histogram.

Even after **the same B_1 has been recovered on both sides**, local port
phase can change the joint value. In F64, F^2 preserves the exceptional
line 3K setwise and sends 3 to 17, yet

    B_1(14,3)=0,  B_1(14,17)=1.

The other exceptional line gives F^2(5)=13 and
B_1(14,5)=1, B_1(14,13)=0. Simultaneously rotating both arguments preserves
B_1; independently rotating only one port need not. Thus recovery of a
form is not recovery of arbitrary relative port alignment.

## 5. Exact stabilizers and what is not reconstructed

Here the group convention matters: we calculate **all binary linear
bijections commuting with the distinguished absolute Frobenius**, then
impose preservation of the spread as a set and/or preservation of B_1.
Commuting with F is stronger than merely normalizing its cyclic group.
Polynomial coordinates name the vectors for calculation; this group does
not fix all those labels pointwise. If alpha itself is required to be a
fixed named field element, the field-automorphism stabilizer is already
trivial, a different and stronger marking convention.

| Field | Invertible centralizer of F | Also preserve B_1 | Also preserve spread | Preserve both |
|---|---:|---:|---:|---:|
| F64 | 24 | 12 | 6 | 6 |
| F256 | 128 | 32 | 8 | 8 |

In both cases the last group is **exactly the Frobenius powers** C_6 or C_8.
This is an exhaustive finite group assertion, not an all-m classification.
The enumeration is small but complete: in a normal basis beta,F beta,...,
a commuting map is determined by u=g(beta), with g(F^i beta)=F^i u. Test
all 2^n values of u and retain those whose conjugates are independent.
Every resulting map and its spread-line permutation is in the receipt.

Adding a marked line means its **setwise** stabilizer. At F64 its order is
6 for K, 3 for either exceptional short-orbit line, and 1 for the other six
lines. At F256 its order is 8 for K and 1 for every other line. Thus even
when one marked short-orbit line recovers the form in F64, a threefold
attachment/transport ambiguity survives. For a primitive field generator
alpha, the formula behind these orbit sizes is: on E^*/K^* of order q+1,
Frobenius sends exponent j to 2j modulo q+1.

For comparison, if the full spread is **not** required to be preserved,
but B_1, F and the selected line are retained, the stabilizer orders are
12 on K and the same 1/3 values on the other F64 lines; at F256 they are
32 on K and 1 on every other line. Thus the stated non-K attachment
conclusions do not depend on covertly adding full-spread preservation.

The bare reconstructed form plus spread plus Frobenius does **not** select
an absolute phase of a normal basis: all Frobenius rotations remain. Nor
does it select one normal generator among different normal orbits. A
sufficient extra marked line can kill the field-automorphism stabilizer,
as the table shows, but that is additional marking; it should not be
misreported as phase recovered from the original three structures.
Prime labels, arithmetic sizes, and independent external port maps remain
outside the recovered form unless their actual maps are explicitly given.

## 6. Relation-to-object scope and reproducibility

The source relation construction sends the 15 nonzero elements of the
four-dimensional Hamming[7,4,3] kernel to the 15 points of F2^4 after a
basis choice. Iterating that same construction from Hamming[15,11,3]
would produce 2^11-1=2047 nonzero relation-objects, **not** 63 or 255.
The present fields test the pairing-reconstruction mechanism in binary
dimensions 6 and 8; they are not falsely identified with the immediate
next Hamming-kernel iteration. An arithmetic labelling of those carriers
would require a separately supplied realization and marking.

Source objects were read at pinned commit
`5ed1835235b36a57b8cd11eb705aa259b6bc0f08`:

* [Free-bits THEORY, P3](https://github.com/eliottcassidy2000/math/blob/5ed1835235b36a57b8cd11eb705aa259b6bc0f08/05-knowledge/results/free_bits_signed_incidence_20261028/THEORY.md), P3;
* [Free-bits PROGRAMME, O2](https://github.com/eliottcassidy2000/math/blob/5ed1835235b36a57b8cd11eb705aa259b6bc0f08/05-knowledge/results/free_bits_signed_incidence_20261028/PROGRAMME.md), O2;
* [Relation gluing THEOREMS, T3c](https://github.com/eliottcassidy2000/math/blob/5ed1835235b36a57b8cd11eb705aa259b6bc0f08/05-knowledge/results/relation_gluing_7_15_20261028/THEOREMS.md), T3c.

An initial read used predecessor pin
`2f97b2fee76e27a78d6979d44a0a2aac027a0314`; the final source pin above was
explicitly supplied by the coordinating session. The two THEORY/THEOREMS
source blobs are unchanged across those pins (respectively
`bb88c509bc19b03bbad6daff46786dba6caa326d` and
`242fe6e7f35daa18aaec8e92f3a8c51b99128bff`). No publisher files or git
state were changed by this subtask.

Run from the directory containing these source files:

    python check_field_forms.py > results.json
    python -O check_field_forms.py > results-optimized.json
    cmp results.json results-optimized.json

The stdlib-only checker uses explicit exceptions rather than removable
asserts. It exhausts all solutions of the exact Frobenius linear equations,
not all 2^28 ambient masks at n=8; nullspace completeness makes this an
exhaustive invariant-form enumeration. It additionally checks every spread
line, every minimal marked-line choice, the entire Frobenius centralizer,
the trace families, and all hostile ordered pairs. The independent audit
is supplied separately by the coordinating session.
