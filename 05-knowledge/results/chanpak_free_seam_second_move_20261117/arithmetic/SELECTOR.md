# Complementary carries improve the free +2 seam to O(N^(1/3))

**PROVED from the inherited FREE_SEAM.md physical count law.** This is an
absolute, constructive, two-sided placement theorem, not exact universal
count coverage. Reflection still changes the selected count by EXACTLY2.

## 1. Retain the two free addresses, but freeze their sum's low bits

Use the physical seven-vertex middle core and reflected-flank construction
of `wheel_square_selector/FREE_SEAM.md`. For positive opposite-parity a,b,
put the ODD member first as u and the EVEN member second as w. The symmetric
lower count is exactly

    h(a,b)=256ab+142(a+b)+81+5(F(a)+F(b)).                 (1)

It belongs to an actual strong same-order pair(h,h+2). Interchanging which
parameter is called a does not change this lower read; the odd/even rule fixes
the actual physical orientation. We do not consume that switch to approximate N.

The original2/5 argument suppressed F by lengthening each address's suffix.
Here no outer suffix is required. Instead choose a+b=m with m≡−1 modulo
a power of two L. The low address words are complementary, so their large
carry contributions cancel. The remaining carry cost is m/L, not m. Since
m is odd, every pair along this line has the required opposite parities.

## 2. Complementary-carry bound (PROVED)

For x≥0 the recurrence is

    F(2x+1)=F(x),   F(2x)=F(x)+x+[x odd],

with F(0)=1, and 1≤F(x)≤x+1. The bound follows inductively; for odd x the
stronger F(x)≤(x+1)/2 when x≥1 supplies the even successor's parity case.
For L=2^t and0≤r<L, let rev(r) reverse its t bits. The suffix formula is

    F(Lx+r)=F(x)+(L−1−rev(r))x
              +[leading bit of r is0][x odd]+F(r)−1.       (2)

This includes x=0 and leading-zero suffixes. It follows by summing contributions
of appended zero bits, as in the inherited reflection proof.

**Lemma.** Let L≥2 be dyadic and m≡−1(mod L). For integers0≤a<a+1≤m,
put S_m(a)=F(a)+F(m−a). Then

    |S_m(a+1)−S_m(a)|
       ≤ |m−2a|+2L+m/L+5.                                (3)

**Proof.** Write a=Lx+r and m−a=Ly+(L−1−r), so x+y=(m−L+1)/L.
If r<L−1, the prefixes x,y do not change in the next step. In the sum of
(2), the variable reversal term is −rev(r)(x−y). Its change has magnitude
at most L|x−y|. The two suffix F terms together lie in[2,L+1], giving a
change≤L−1; the two parity-boundary terms change by at most2. Finally
L|x−y|≤|m−2a|+L−1. This proves the stated bound in this case.

If r=L−1, the next prefixes are x+1,y−1, with y≥1. Using the all-one and
all-zero suffixes directly, the change is

    F(x+1)−F(x)+F(y−1)−F(y)
      +(L−1)(x+1−y)+[x+1 odd]−[y odd].

The first two differences have total magnitude at most x+y+1≤m/L+1.
Since L(x+1−y)=2a−m+1, the linear term is at most|m−2a|+1.
The parity term is at most1. This gives (3), including this otherwise
potentially dangerous carry boundary. QED.

Equation(1) and(3) give the adjacent-count estimate

    |h(a+1,m−a−1)−h(a,m−a)|
      ≤261|m−2a|+10L+5m/L+281.                           (4)

No monotonicity of either F or h is assumed. Both endpoints are physically
available: if a changes parity, interchange its assignment to the two flank roles.

## 3. Universal constructive1/3 theorem, including an explicit threshold

**Theorem.** For every integer N≥2^18 there are two actual free-seam lower
counts h_low≤N≤h_high, each belonging to an exact+2 reflected pair, with

    max(N−h_low,h_high−N) ≤1000 N^(1/3).                  (5)

The addresses can be selected with O(log N) exact F evaluations, hence
O((log N)^2) elementary integer-arithmetic steps under the binary recurrence.
This operation count does not assume constant-time arbitrary-size integers.

**Construction.** Let M be the least integer with64 M²≥N. Thus M≥64.
Take L the largest power of two with L³≤M. Choose the least m≥M with
m≡−1(mod L), and set

    a0=(m−1)/2,  d=ceil(sqrt(2mL)),  a_left=a0−d.

All a,b used are positive: m≥64, L≤m^(1/3), and d≤(3/2)m^(2/3)+1≤3m/8+1.
Also m<M+L, L>M^(1/3)/2, and m≤2M.

At a0, the polynomial part of(1) is64m²+142m+17, already greater than N.
At a0−d its product part falls by256d²+256d. The F contribution is at most
5m+10. Since m−L≤M−1<sqrt(N/64),

    64m²−N <128mL.

Consequently the full value at a_left is below N: its excess is less than
128mL+147m+27−256d²−256d<0, because d²≥mL and L≥4.
This proves actual sign bracketing, including the full F term.

Use binary subdivision of the integer interval[a_left,a0], retaining one
endpoint below N and one endpoint at least N. This does NOT require h to be
monotone: at every midpoint keep the half whose endpoints bracket N. It ends
with adjacent addresses, or an exact hit. Throughout this interval
|m−2a|≤3m^(2/3)+3. Furthermore m/L≤3m^(2/3). Equation(4) is therefore
less than1000m^(2/3) for m≥64. Finally m≤sqrt(N)/2 for the stated threshold,
so this is certainly less than1000N^(1/3). The two adjacent endpoint values
prove(5). QED.

The constant is deliberately loose. The exponent improvement, not numerical
constant optimization, is the point. Initial cases below2^18 are not covered
by this asymptotic theorem; a finite search could handle a desired bounded
target list without being part of the all-order argument.

## 4. Adaptive carry sectors reach the quarter exponent

The proof can exploit how close the target's scale M lies to a dyadic sum.
For ANY dyadic L≥2, choose m≥M minimal with m≡−1(mod L) and put delta=m−M.
Replace d above by ceil(sqrt(2m(delta+1))). Whenever this leaves a_left≥1,
the same center and left-end calculations bracketN, and the final error is

    O(sqrt(m(delta+1))+L+m/L).                            (6)

The replacement uses64m²−N<128m(delta+1). All constants are absolute.
Testing O(log N) choices of L and retaining the best bracket is an exact
adaptive selector; it need not change the physical core or the+2 switch.

In particular, if delta is bounded and L is comparable to sqrt(M), (6)
is O(N^(1/4)). Such sectors occur infinitely often: take
M=2^(2t)−1, L=2^t and every N with64(M−1)²<N≤64M². Then delta=0.
This is a genuine stronger bound on explicitly specified target bands,
NOT a universal quarter-exponent theorem. For arbitrary M, balancing the
worst-case terms sqrt(mL) and m/L gives the universal choice L≈m^(1/3).

## 5. What remains

Complementary carries supply a new absolute approximation controller without
spending the signed reflection move. They do not correct the remaining residual
exactly. No density, primality, or Chan–Pak completeness statement follows.
The sharper next question is whether another physical move can correct the
O(N^(1/3)) residual while preserving this low-sum-bit cancellation and the
odd/even reflected sidecar.

`selector.py` uses only integer comparisons and integer square roots. `check.py`
checks the suffix law, every carry-boundary case in a finite box, exact physical
formula agreement using an independent direct prefix-sum read, bracketing across many
scales, and quarter-exponent target bands. These tests support the proof;
they are not its universal premise.
