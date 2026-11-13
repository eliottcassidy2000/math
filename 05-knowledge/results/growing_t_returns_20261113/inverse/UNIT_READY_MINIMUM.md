# The exact minimum at every strong unit-ready order

**PROVED IN NOTE, with inherited OCF and the complete unique-start grammar.**
Let mu(M) be the smallest full Hamilton-path count of a **strong**
M-vertex tournament with a marked unit-ready root. Such sources exist
exactly for M>=4. Put r=M-4. Then

    mu(4)=5, mu(5)=9, mu(6)=17, mu(7)=33,

and for M>=8, with s=floor((r-1)/2),

    mu(M)=3*2^r+3*2^s+2^(r-s)-1.                       (1)

Every minimizer has a singleton core in its canonical funnel decomposition.
There is one rooted isomorphism type at each order except M7, where there
are two; those two also have different unrooted degree multisets.
This classifies an infinite extremal slice, not all tournaments or all
Hamilton counts. Chan--Pak and the target-preserving405 return remain OPEN.

## 1. Inheritance and the physical reduction

Unit-ready means F-H=1, using the ordered two-path-cover root convention
of [READY](../../two_entry_funnels_20261103/READY.md). That classification
proves its strong source is exactly a newborn y, beaten only by a root x
which has one Hamilton path starting at it in the strong deletion parent.
The [complete unique-start grammar F1--F2](../../polarized_path_calculus_20261103/FUNNELS.md)
then writes that parent uniquely as a funnel over a strong marked core C,
or a singleton, with a transitive binary tail and the root x.

Let c=|C|, let the tail length be L, and let its low-bit-first integer
address be n. The full source order is

    M=c+L+2.

If L>0, its final bit is1, so 2^(L-1)<=n<2^L. The singleton core requires
L>=1; for a nontrivial strong core c>=3, the empty tail n0 is also allowed.
The canonical grammar is a bijection of rooted isomorphism classes, not
merely a sufficient family.

Write h=H(C), e=E_C(1), p=P_C(2), f=E_C(2), and S=S_C(2). Here p counts
two ordered spanning paths with empty bins allowed; E and S require the
first bin to end or start at the marked core vertex. The inherited full
count is

    N(n)=(S+p)n+h g(n)+e(n mod2)+(S+p-e-f-1).            (2)

Put N0=S+p+h-e-f-1. Since g(n)>=1,

    N(n)=N0+n(S+p)+h(g(n)-1)+e(n mod2) >= N0+n p.       (3)

We use two elementary cover injections:

* S>=h: cut a Hamilton path immediately before the root, put its suffix
  in the distinguished first bin and its prefix in the second;
* f<=p/2: bin exchange pairs covers having the root in the first bin
  with those having it in the second; ending there is a further restriction.

Together with e<=h, these give

    N0>=p/2+h-1.                                      (4)

The closest prior inverse used the general strong-order bound. The new
mechanism is a stronger two-cover floor, which separates every nontrivial
core from the optimal singleton arrangements.

## 2. A strong c-vertex core has P(2)>=c*2^(c-1)

Let z(C) be the number of directed triangles. First,

    z(C)>=c-2.                                        (5)

For completeness, sort the outdegrees increasingly. Strongness makes
the sum of the smallest k outdegrees at least binom(k,2)+1 for every
1<=k<c: any proper nonempty vertex set has at least one outgoing arc.
For c>=4 the sorted vector

    (1,1,2,3,...,c-3,c-2,c-2)

has exactly these prefix sums and the same total. It therefore majorizes
the actual score vector. Convexity of binom(x,2) shows that the number
of transitive triangles, sum_v binom(outdegree(v),2), is at most
binom(c,3)-(c-2). This proves(5); c3 is the directed triangle itself.

The inherited **THM-002, Odd-Cycle Collection Formula**
([canonical file](../../../../01-canon/theorems/THM-002-ocf.md), with its
Grinberg--Stanley/Irving--Omar attribution) implies for every induced
subtournament U that H(U)>=1+2z(U). Summing over the first-bin vertex
set A gives

    p=sum_(A subset C) H(A)H(C-A)
      >=2^c+2 sum_A[z(A)+z(C-A)]
      =(z(C)+2)2^(c-1)
      >=c*2^(c-1).                                    (6)

Each triangle is wholly contained in either specified bin for exactly
2^(c-3) subsets, which supplies the middle identity. This uses only
positive OCF terms and keeps all path-cover counts ordinary and integral.
The score inequality and OCF are inherited mathematics; their use in
the extremal reduction is the contribution here.

## 3. Every non-singleton core loses to a singleton construction

Suppose c>=4. If L0, equations(4),(6) and h>=3 imply

    N>=c*2^(M-4)+2 >=4*2^(M-4)+2.

If L>=1, equation(3) and n>=2^(L-1) imply

    N>=c*2^(M-4)+c*2^(c-2)+2 >4*2^(M-4)+1.

At the same order M, a singleton core with address n=2^(M-4) has
full count4*2^(M-4)+1. Thus c>=4 cannot minimize.

The only strong core of order3 is the directed triangle. Its actual
fields are(h,e,S,p,f)=(3,1,3,12,3), so

    N(n)=15n+3g(n)+(n mod2)+10.

With L0 this is13 at M5, larger than the singleton count9. For L>=1,
it is at least15*2^(L-1)+13. At the same order M=L+5, the singleton
address2^(L+1)+1 gives14*2^(L-1)+5, which is strictly smaller.
Hence **every global minimizer uses the singleton core**.

That singleton has(h,e,S,p,f)=(1,1,1,2,1), and its complete count is

    N(n)=3n+g(n)+(n mod2),       bitlength(n)=M-3.       (7)

This also proves the useful preliminary bound

    H>=3*2^(M-4)+1,
    M<=4+floor(log2((H-1)/3))                          (8)

for every strong unit-ready source. The exact formula(1) is sharper.
There is no analogous bound on the total order of weak ready sources:
appending transitive sinks can preserve H indefinitely.

## 4. The scalar minimum and all equality addresses

Let r=M-4>=1 and write n=2^r+x, 0<=x<2^r. The inherited suffix law in
[central compiler C1](../../central_modular_chanpak_20261110/compiler/PROOF.md#c1-a-four-point-detector-for-the-deletion-core),
or induction in the r low bits of the recurrence, gives

    g(2^r+x)=2^r-rev_r(x)+g(x)-top_r(x).

Writing x_i for its r bits yields

    N(2^r+x)=4*2^r
      +sum_(i=0)^(r-1)(3*2^i-2^(r-i-1))x_i
      +g(x)+(x mod2)-top_r(x).                        (9)

For0<x<2^(r-1), the last line is at least2; if x is even it is at
least3. It equals2 exactly when x is a string of low ones. This uses
the recurrence fact that g(x)=1 exactly for x0 or x=2^a-1. For the
upper half x>=2^(r-1), the last line is at least1.

The linear coefficients in(9) are strictly increasing and none is0.
For r>=3 their negative entries are exactly the initial
s=floor((r-1)/2) entries. Thus the linear part is minimized by

    x*=2^s-1.

This lies in the lower half and attains residual2. Every upper-half x
has, in addition, the positive top coefficient3*2^(r-1)-1, so its
possible one-unit reduction of the residual cannot win. Every other
nonzero lower-half x has strictly larger linear part and residual at
least2. The separate x0 value is4*2^r+1.

At r3, the sum of negative coefficients is-1, so x0 and x1 tie.
At r>=4 it is at most-5, so x* is the unique minimizer. Substituting
x* into(9) proves(1). Directly checking the boundary r0,1,2 gives
the earlier values5,9,17 and the respective unique addresses1,2,4.

For r>=4 the actual extremal tail is, low bit first,

    1^s 0^(r-s) 1,       s=floor((r-1)/2).             (10)

Thus the optimal arrangement is an explicit infinite family with two
transitive runs and one closing bit. Its complete marked graph is
obtained by the inherited singleton funnel and newborn birth.
Canonical unique-start recognition and the unique newborn predecessor
make these equality addresses distinct rooted isomorphism types.
At M7 the two degree multisets are

    (1,1,2,3,4,5,5),       (1,1,2,4,4,4,5),

so even their unrooted tournaments are different.

## 5. Enumeration consequences and the unresolved k17 inverse

The exact minima begin, at M4,

    5,9,17,33,61,115,219,423,823,1615,3183,6303,12511,...

In particular mu(M)~3*2^(M-4), with an explicit error of order2^(M/2).
This is an all-order extremal theorem, not an extrapolation of that
finite prefix. It supplies a strictly smaller finite head for the
complete ready-value inverse. No new A000568 term is claimed.

For N_k=405*2^k+71, equation(8) gives source order at most k+11.
The targeted k17 inverse now excludes every canonical inner core
through order18 at every outer-word length. Equations(3),(4),(6) give
the additional necessary inequality

    N>=c*2^(c-2)(2n+1)+2,       c>=3.                 (11)

At N=53084231, only these twenty core-size/address choices remain:

|inner core size c|possible outer addresses n|possible tail lengths|
|--:|---|---|
|19|0 through10|0 through4|
|20|0 through4|0 through3|
|21|0 or1|0 or1|
|22|0|0|
|23|0|0|

Every c>=24 is excluded by(11). Counting M=c+L+2 now restricts any
undiscovered strong unit-ready witness to source orders21..25. The
table restricts constructor parameters; it does not enumerate all
large cores in those cells and does not assert that any cell is nonempty.
The earlier defect-three witness at the same full count remains valid
and is not relabelled as unit-ready.

## 6. Exact controls and reproduction

Run `python -B 05-knowledge/results/growing_t_returns_20261113/inverse/unit_ready_minimum.py`
and repeat with `-O -B`. The checker independently tests:

* every scalar address of bitlength at most19, including all ties;
* the score/triangle/P2 inequalities on all6008 strong order-eight
  representatives from the inherited complete6880-class census;
* every196417 unique-start rooted parent class through order15,
  comparing the actual minimum and number of equality classes;
* literal adjacency matrices for every extremizer of orders4..15,
  including full H, predecessor start count1, and actual readiness.

The census is inherited and pinned, not regenerated. The all-order
claim is the proof above. `joint_probe.py --order 18` separately records
the2,178,309-class next-core exclusion for the k17 target; both normal
and optimized receipts agree. These computations retain explicit
runtime checks rather than removable assertions.
