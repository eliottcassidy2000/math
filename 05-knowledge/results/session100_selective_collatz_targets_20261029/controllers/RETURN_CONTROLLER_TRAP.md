# A return-controller trap separates 3x+1 from 5x+1

**Status: PROVED** elementary classification of the stated auxiliary controllers;
**FINITE-EXACT** replay below. Not a convergence theorem or a novelty claim.

Inheritance: Session97 `controllers/INTEGER_TERMINATION.md` and
`SOURCE_PRESERVING_RETRIES.md`. The maintained source is Session99 verified-final,
HEAD `1e742862a5a442c87a6d7271cbc5e8305e2c4482`.
Concept board: native source / ternary queue / shifted potential / integer guard /
cycle obstruction / inverse cylinder. The crucial sidecar is the queue at a native
cycle entrance, not merely its contracting affine multiplier.

## 1. Exact controller and natural shift

For odd multiplier b>=3, let T_b(n)=(bn+1)/2^a, a=v_2(bn+1).
Fix r>=1, start Q=1, add 2 when r divides a, subtract 1 otherwise,
and stop at the first Q=0. This is always a ternary queue; r is the
valuation modulus, **not** its branching degree. For k branches a completed
word has J=3k+1 and D>=1+(r+2)k, so its affine slope is at most

    (b/2) [b^3/2^(r+2)]^k.

This bound describes completed words only. It cannot establish completion.

For r>=2, valuation 1 is a leaf, and the unique shift h making every
valuation-1 leaf preserve P=(n+h)(b/2)^Q is h=1/(b-2). For r=1
there are no leaves; the same shift still gives the formulas below, but
the leaf-based uniqueness characterization does not apply. Every actual
leaf weakly decreases P, since

    P'/P = (2/b)(bn+1+2^a h)/(2^a(n+h)) <= 1.

The inequality follows by subtracting the numerator from the denominator;
equality holds at a=1. At a branch,

    P'/P = (b^2/4)(bn+1+2^a h)/(2^a(n+h)).       (1)

For fixed a>=1 the last fraction is nonincreasing in n; for fixed n it is
strictly decreasing in a. These facts follow directly by differentiation,
or by cross multiplication. All potential calculations use exact rationals.

## 2. A finite integral-guard theorem

Suppose Delta=2^(r+2)-b^3>0 and define

    H=[b^2+(b+2)2^r]/Delta.

At a branch of exact valuation r, cross multiplication in (1) gives

    P'/P <1  iff n>H;  P'/P=1 iff n=H.

At larger branch valuations the ratio is smaller. This threshold is not an
unrelated estimate: the balanced block (r,1,1) has exactly

    C=b^2+(b+2)2^r, D=r+2, J=3,
    F(n)=H+[b^3/2^(r+2)](n-H).

Its queue increments are (+2,-1,-1). Thus a legal integral fixedpoint H of
this exact word is automatically an infinite active queue trap. A rational
fixedpoint or an integral point not realizing the exact valuations is not.

**Finite-guard theorem.** Put E={positive odd n<=H: r divides v_2(bn+1)}.
Every infinite active run visits E infinitely often. Consequently, if every
n in E has an independently certified terminating active run for **every**
positive initial queue, the controller terminates for all positive odd inputs.
In particular an empty E suffices.

Proof: let m be the least positive odd integer strictly above H. On all branch
inputs outside E, (1) is bounded by its value at (n,a)=(m,r), a fixed rho<1.
This bound is valid even if m does not itself realize valuation r. Once an
infinite run has made its last visit to E, its potential is bounded above by
a constant times rho to the number of subsequent branches. Active states have
P>=(1+h)b/2>0, while their queue identity forces infinitely many branches.
This is a contradiction. If every E-entry universally terminates, not even
one E-entry can belong to an infinite run. QED.

This is a genuinely finite arithmetic guard reduction for *auxiliary first
returns*, not for repeated descent. A finite E need not be easy to certify,
and its certification may fail: the5-case below fails precisely at83.

## 3. The (3,3) integral guard, and its failure at (5,5)

At (b,r)=(3,3), h=1 and the smallest branch input is n=13.
Formula (1) is at most27/28, sharply at n=13,a=3. Thus every branch
uniformly contracts, every leaf weakly contracts, and every active state has
P>=3. An infinite active run would require infinitely many branches because
Q=1+3k-J>=1. This proves first-return termination, as in Session97.

At (b,r)=(5,5), h=1/3, the exact branch classification is:

|branch input|valuation|P'/P|
|---|---|---|
|19|5|125/116 >1|
|83|5|1|
|any other input|multiple of5|<=875/884 <1|

Proof: a=5 iff n=19 mod64. After19 and83 the smallest input is147,
where (1) is875/884. For a>=10, n>=3 (n=1 has valuation1), and
monotonicity bounds (1) by its value at n=3,a=10, namely335/512,
which is less than875/884. The first two rows are direct evaluations.
With the noncanonical shift n+1, even83 is expanding, by25/24.

## 4. Exact classification of infinite active runs

**Theorem.** For the (5,5) controller, from any positive odd n and any
positive integer queue Q, the active run is infinite iff it visits83
before stopping.

**Proof.** At83 the native orbit is

    83 --a=5--> 13 --a=1--> 33 --a=1--> 83.

The queue runs Q,Q+2,Q+1,Q, so any positive Q survives forever.
If the run instead visits19 while active, its next states are3,1,3,1,...:
the valuations are5 followed by4,1,4,1,... . After that first branch every
step is a leaf, hence the queue reaches0 in finitely many steps.
An infinite run avoiding83 therefore also avoids19. Every branch contracts
P by at most rho=875/884 and every leaf weakly contracts. The active lower
bound P>=(1+1/3)(5/2)=10/3, together with Q=Q_initial+3k-J,
contradicts infinitely many branches. This proves both directions. QED.

This gives a terminating three-way diagnostic: stop at Q=0 (completed
return), at19 (guaranteed eventual completion), or at83 (infinite trap).
It is not merely a finite search for cycles. On trajectories avoiding both
special inputs, (10/3)<=P_initial*rho^k bounds the branch count explicitly;
the queue identity then bounds the accelerated clock.

## 5. Hostile source cylinders and the missed fixed-point coordinate

The block w=(5,1,1) has J=3,D=7,C=249 and native affine map

    F_w(n)=(125n+249)/128 = 83+(125/128)(n-83).

For every k>=1 the exact positive odd sources realizing w^k are

    n == 83 (mod 2^(7k+1)).                       (2)

Indeed composition gives F_w^k(n)=83+(125/128)^k(n-83).
The exact valuation cylinder has modulus2^(7k+1), with a unique odd residue;
83 realizes it, proving (2). Alternatively use the exact inverse guard below.
Every such prefix leaves the queue at its initial value and never lowers it
inside a block. Source83 realizes all prefixes, hence never completes.
Notice that the block's slope125/128<1 is perfectly compatible with a
positive integer cycle; the carry249 fixes its height83.

The analogous (3,3) block v=(3,1,1) has J=3,D=5,C=49:

    F_v(n)=(27n+49)/32 = 49/5+(27/32)(n-49/5).

Its repeated-word source cylinder is

    n == 49 * inverse(5) (mod 2^(5k+1)).          (3)

The common 2-adic point is49/5, not a positive integer. Thus this particular
balanced queue cycle is arithmetically excluded for3, not excluded by slope.
The integer termination theorem excludes *all* infinite active (3,3) words;
(3) illustrates its integral obstruction, not an alternative proof of it.

For either multiplier, retain exact native source guards and carry:

    C_next=b*C+2^D; D_next=D+a; J_next=J+1;
    2^D n_out=b^J n_source+C;
    n_source == (2^D-C)*inverse(b^J) (mod 2^(D+1)).

The inverse guard determines the exact word, not just total divisibility.
It is verified by chronological replay. The corresponding symbolic/native
map preserves source, chronological word, endpoint and clocks; throwing away
C destroys the cycle height and would erase the decisive distinction above.

## 6. Consequence and remaining debt

There **is** a useful bijective reframe: exact valuation word plus its inverse
source cylinder, with queue and carry retained. The natural shifted potential
makes the3-versus5 distinction an explicit integer-guard question. There is
**not** a universal generalized termination theorem based only on the two
completed-tree slope bounds3/2*(27/32)^k and5/2*(125/128)^k.

Even for3, restarting Q=1 at every completed return multiplies the potential
by3/2. First-return termination does not repay these resets or prove descent
below the original marked source. Any stronger rank must account for the
reset, affine carry and source. For5, no strictly decreasing well-founded rank
on the active83-cycle can exist, even with source-preserving finite sidecars
that return after each cycle. An unbounded sidecar cannot supply a well-founded
strict decrease along an infinite active trajectory either.

Reproduce: `python3 -B check_return_controller_trap.py`, and repeat with `-O`.
The checker explicitly states finite ranges and includes rational-hostile,
positive-integer trap, inverse-guard and first-return controls. Its computations
check the formulas; the infinite statements depend on the proofs above.
