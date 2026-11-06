# A sparse-family bounded-source obstruction

Status: RESERVED / UNREVIEWED proof candidate. The twelve-template occupancy
lemma remains OPEN. This is not an occupancy proof or a counterexample to it.

## Statement

Use the twelve matrices of the reviewed Session 201 compact central arithmetic
packet, always on actual tau(u)=(g(u),u mod2,u,1), u>=2. If

    2^k-1 = Psi_W(u,v),    g(u),g(v)<=G,    G>=1,

then

    2^k-1 <= 18630 G^8.                                     (1)

Therefore every fixed bound on both g coordinates misses all sufficiently
large Mersenne targets. In particular any hypothetical uniform occupancy
construction must have max(g(u),g(v))>=(n/18630)^(1/8) along n=2^k-1.
This is stronger on this sparse target sequence than the inherited bounded-g
counting barrier, which alone does not exclude infinitely many Mersenne values.
There is no claim about bounding just one source coordinate, nor optimality of
the exponent or constant.

## Inheritance

The factorization candidate's elementary odd-core lemma, restated here:
write u=2^r a with a odd. Then g(u)=g(a), and 1<=g(a)<=a; for a>=3,
g(a)>=(a+3)/2. Hence g(u)<=G implies a<=2G (also for a=1).
Additionally g(2^j+1)=2^j+1 for j>=1, by the defining recurrence.
These facts have elementary direct proofs and require no occupancy premise.
The reviewed physical product bound is Psi_W(u,v)>=uv/2.

## Nonempty templates

Allow the transpose when exchanging the two source indices. The frozen matrix
list gives, writing A=g(u), B=g(v), p=u mod2,q=v mod2,

    n=lambda*A*B+c*A*v+d*u*B+b*u*v
      -e*A*q-f*p*B-h*p*v-i*u*q+j*p*q,

where lambda is positive odd, lambda<=5; 0<=c,d<=7; 0<=b<=8;
e,f,h,i>=0; and -1<=j<=1. The fixed-source g coefficient is strictly positive
at every actual source. This also follows directly from the matrix list.

Write u=2^r a,v=2^s b0, with odd a,b0<=2G, and exchange indices/matrix so r>=s.
If r=s=0 then u,v<=2G and n<=5G^2+28G^2+32G^2+1<=66G^2.

If r>=1 and s=0, then p=0, q=1 and

    n = A*(lambda*B+c*v-e) + u*(d*B+b*v-i) = C+uD.

Here C>0 by the fixed-source coefficient positivity, and C<=19G^2.
The product bound and v>=3 give n>=u, hence k>r. Reduction of n+1=2^k
modulo 2^r yields 2^r | C+1, whence 2^r<=20G^2. Thus u<=40G^3,v<=2G,
and n<=940G^4 using the displayed coefficient bounds.

If r>=s>=1, both parities vanish. The product bound gives n>=u and hence k>r.
Reduction modulo 2^s gives

    2^s | lambda*A*B+1,     2^s<=5G^2+1<=6G^2,
    v<=12G^3.

Now n=C+uD with C=A*(lambda*B+c*v)>0. A second reduction modulo 2^r gives

    2^r<=C+1<=5G^2+84G^4+1<=90G^4,
    u<=180G^5.

Consequently

    n<=5G^2+7G*(180G^5+12G^3)+8*(180G^5)*(12G^3)
      <=18629G^8.

The stated 18630G^8 bound safely covers all three cases.

## The empty template and its zero-constant boundary

For the empty form,

    n=A+B+uv-u-v+pq-p-q+1.

If both indices are odd they are at most 2G and n<=4G^2+2G+1.
Suppose r>=1; again exchange indices so r>=s. The product bound implies k>r
(the even/even case has v>=2 and the even/odd case v>=3).

If s>=1, reduction of n+1 modulo 2^s yields

    2^s | A+B+2,     2^s<=2G+2,     v<=8G^2.

Write n+1=(v-1)u+Z, where Z=A+B-v+2. If Z=0, then v=A+B+2>=4,
so the odd number v-1>=3 divides 2^k, impossible. Thus Z!=0 and

    2^r<=|Z|<=2G+v+2<=12G^2,
    u<=24G^3,
    n<=uv+A+B+1<=195G^5.

If s=0, then v is odd, v<=2G, and

    n+1=(v-1)u+Z,    Z=A+B-v+1.

If Z=0 then (v-1)u=2^k forces u and v-1 both to be powers of two.
Therefore A=1 and v=2^j+1 with j>=1, so B=v. But Z=A+B-v+1=2,
a contradiction. Hence Z!=0, so

    2^r<=|Z|<=2G+v+1<=5G,
    u<=10G^2,
    n<=uv+A+B+1<=23G^3.

These bounds are all less than 18630G^8. This handles the branch which a mere
sign argument would miss. The zero-constant case is ruled out arithmetically,
not dropped by a nonzero-divisor assumption.

## Limit

This argument identifies a necessary growing source coordinate on an infinite
specific target sequence. It does not prove that those targets are missing:
G is allowed to grow with the target, and the finite-menu occupancy problem
remains exactly where it was. No literature novelty claim is made.
