# Two local product shortcuts fail, but the global matching refinement survives

Status: exact counterexamples to two proposed *conditional* margin factors.
Neither example refutes the weighted stability conjecture, the accepted
Gram/Cauchy product bound, or an absolute weight-sensitive matching margin
bound. The rational checks are in `verify_local_factor_hostiles.py`.

Write f_+=u+v, where u=det(I+A), v=(-1)^n det A. Vary one skew entry
x=2p-1, put q=p(1-p), and let L be the endpoint interpolation of f_+.
The accepted local Gram factor is1/(1+4q). The following stronger factors
cannot be used without additional hypotheses or a compensating reserve.

## 1. The factor1-2q fails with other fractional cross edges, already at n5

Let the three-vertex binary skew core be

    C = [[ 0, 1, 1],
         [-1, 0, 1],
         [-1,-1, 0]],

and let the two cross rows be

    U = [[ 1/2, 1/2, 0],
         [-1/4, 1/4, 1/2]].

Complete the skew matrix as K=[[J_x,U],[-U^T,C]], with
J_x=[[0,x],[-x,0]], and set A=(J-I+K)/2. Every probability is in[0,1].
A direct Schur-complement calculation gives

    det(I+K)/det(I+C)=x^2-3x/8+81/64,
    f_+(x)=(1/4)(x^2-3x/8+81/64).

At x=3/4, p=7/8 and q=7/64. The exact values are

    u=165/256,          det A=33/128,
    f_+=u-|det A|=99/256,
    f_+(-1)=169/256,    f_+(1)=121/256,
    L=127/256,         1-2q=25/32,
    f_+-(1-2q)L=-7/8192<0.

At both endpoints the determinant sign is also positive, so the same
example refutes replacing f_+ by the actual minimum margin u-|det A|
in this conditional factor. The cross edges here are fractional; this
example does not address a matching-only uncertainty graph.

## 2. The factor1-q fails even for matching-only uncertainty, at n8

Take the ordinal sum of a weighted pair and two directed3-cycles. All
cross edges are binary. Only the pair is uncertain. Block triangularity
makes the determinant calculation immediate:

    u=4(1-q),       det A=-q,       f_+=u-|det A|=4-5q.

Both binary endpoint margins equal4, so L=4. Therefore

    f_+-(1-q)L=-q<0                 whenever q>0.

For a fair pair, q=1/4, the values are f_+=11/4 and (1-q)L=3.
This disproves a naive conditional1-q factor even with a binary remainder
and a matching uncertainty graph.

The absolute proposed matching floor1-q is NOT refuted: here4-5q>=1-q
for0<=q<=1/4. The nontransitive binary remainder carries substantial
excess margin. Thus an absolute matching-weight refinement could still
be true, but its proof must transfer that excess rather than multiply
a stronger conditional factor independently at every step.

For comparison, replacing both3-cycles by an ordinal binary transitive
remainder gives the sharp absolute margin1-q whenever there is a singleton.
The accepted support-only matching theorem and the full Gram/Cauchy proof
are unaffected by either hostile.

## Implication for the remaining finite range

The n>=526 dual certificate uses the accepted dimension-free gap bootstrap
only through c<1/4. Its local quadratic dual, after paying the deleted-edge
mass budget, gives h=d/16 with d<c. Thus this particular certificate has
an ideal ceiling below1/64, corresponding to dimension around512. Merely
optimizing its rational constants or geometric tail cannot settle small
orders10..525. A qualitative improvement must strengthen an input or use
more of the actual tournament structure. Weight-sensitive absolute matching
bounds and binary/Pfaffian defect compensation remain plausible routes,
not results proved by these examples.
