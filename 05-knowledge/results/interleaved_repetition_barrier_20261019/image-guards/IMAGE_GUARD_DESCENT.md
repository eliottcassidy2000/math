# Two-sided guards and 3-adic image-type descent

19 October 2026. **PROVED, using the separately audited matched-block height
theorem and its cited logarithmic-form input.** This packet supplies a new
actual endpoint exclusion with unbounded episodes. Arbitrary bank entry,
Collatz convergence, and the complete realized S_d boundary remain open.

## 1. The bank-wide image-type theorem

For a positive integer z, write nu(z)=v_3(z-1), with nu(1)=infinity.
For every fully lawful bank transition the output types are exactly

| Letter | Odd-bit count h | v_3(b-M) | Output type |
| --- | ---: | ---: | --- |
| K | 19 | 3 | 3 |
| V | 77 | 4 | 4 |
| H | 4 | 1 | 1 |
| F | 1 | 1 | nu(input)+1 |

For K,V,H, the affine formula y=(3^h x+b)/M gives

    y-1=(3^h x+b-M)/M.

The denominator M is a power of two. In each of these three cases,
h>v_3(b-M), so the first numerator term cannot cancel the displayed lower
3-adic valuation of b-M. This proves the exact output type for every
integer input whose branch is lawful. The table is not inferred from
sampled values.

For F, h=v_3(b-M), so that shortcut would be invalid. Instead its exact
formula is

    F(x)-1=3(x-1)/4.

Thus F increments every finite type by one; if x=1, its output is again 1
and the type remains infinity. This exceptional fixed point is included
explicitly, not assigned a finite valuation.

Consequently, after the first K,V or H, a bank checkpoint's type is obtained
by resetting to 3,4 or 1 at those letters and incrementing at each F. The
height may be unbounded; this is a typed-counter rule, not a claim of a
finite-state quotient carrying all arithmetic.

## 2. General past-image and future-cylinder guards

Let A be any lawful past parity word ending at checkpoint y, with length n,
h odd bits and branch (3^h z+b)/2^n. Then

    y = b*2^(-n) mod 3^h.                        (1)

Indeed 2^n y=3^h z+b, and two is invertible modulo 3^h. This is an image
restriction, distinct from the familiar source restriction modulo 2^n.
A past suffix is enough: the integer state entering that suffix need not
be the original source of the whole word.

Suppose two checkpoints y_i,y_j support the same following r-bit parity
word. Its unique lawful source residue gives

    y_i-y_j=2^r k for some integer k.             (2)

If their past-image restrictions give unequal residues rho_i,rho_j modulo
3^s, then y_i!=y_j. Therefore, in an expansion-first first-crossing bank
word with total expansion-bit length L,

    100r>=L and rho_i!=rho_j mod 3^s              (3)

imply strict endpoint descent. This is the audited arbitrary-start
matched-block theorem with its exact-return alternative removed by (1).
The earlier context, all source guards, and the actual integer checkpoints
are retained. A symbolic factor comparison without these past guards is
insufficient.

In particular different finite types from Section 1 suffice. If the types
are t_i!=t_j, then v_3(y_i-y_j)=min(t_i,t_j), so their equality is impossible.
This applies even if their separator contracts and their observed future
continuations agree indefinitely.

## 3. Optional quantitative CRT refinement

Let a be the least residue of (rho_i-rho_j)*2^(-r) modulo 3^s. If a!=0,
(1)--(2) imply

    |y_i-y_j| >= 2^r * min(a,3^s-a).              (4)

This is the least absolute integer in the residue class of k. It requires
no equidistribution assumption. Put D_3=min(a,3^s-a). Combining (4) with the
first-crossing height bounds gives the necessary condition for non-descent

    r log2+log D_3
      < (L/30)log lambda+B+(B+1)log N+122 log2,   (5)

where lambda=3^19/2^30 and B=10^10. At macroblock starts the last term can
be omitted. One may also retain the actual prefix expansion length in
place of the full L. These are quantitative necessary conditions; (3)
already supplies a convenient all-length endpoint exclusion.

## 4. The V versus VF obstruction

Exact output residues modulo 243 are

    V output: 82 mod 243,
    VF output: 1 mod 243.                        (6)

For V, b_V=103 and M_V=22 modulo 243, and 22^(-1)=232, giving 82.
Appending F transforms that image residue to (3*82+1)/4=1 modulo 243.
Equivalently, the two output types are 4 and 5.

Their residue difference is 81 modulo 243. For every r, the multiplier
2^(-r) leaves it either 81 or 162 modulo 243, so (4) gives the stronger
separation |y_i-y_j|>=81*2^r. In particular they can never be an exact return.

## 5. An all-length greedy family excluded from S_d

Use the greedy first-crossing word W_m built from V and F in the companion
notes. Put theta=log(s_V)/log(4/3). Suppose integers 0<p<q, q>=2, m>2q obey

    gcd(p,q)=1, mp=1 mod q,
    1/(mq)<p/q-theta<1/((m-q)q).                 (7)

The residue proof in PERIODIC_TAIL_MISMATCH_HOSTILE.md shows

    W_m = V B^k R,

where B is a bank word of q iteration images, ends in VF, has length
P=122q+2p, and R is a prefix of B. In particular the continuation following
the first V is periodic as a formal string throughout the rest of W_m.

Set L=122m and r=ceil(L/100). There are equal r-bit factors starting just
after V and just after VB. They fit wholly in W_m: r<=2m, P<124q,
N>=122m, and m>2q with q>=2 give

    122+P+r < 122+124q+2m < 122m <= N.

The first checkpoint follows V and has type 4. The second follows the end
VF of B and has type 5. Section 2 therefore proves

    T^N(x)<x for every positive lawful source of W_m. (8)

No internal-cycle alternative remains. Thus no such word can witness
membership in S_d for any d>=0. This is an endpoint-realizability exclusion,
not merely the assertion that the formal word has no cycle.

## 6. The family has unbounded episodes without a partial-quotient hypothesis

Theta is irrational: a rational value p/q would imply
3^(77q+p)=2^(122q+2p), contradicting unique prime factorization.

Take any upper simple-continued-fraction convergent p/q of theta, and let
p'/q' be the preceding convergent. If t>1 is the next complete quotient,
the elementary convergent identities are

    p q'-p' q=1,
    theta=(p t+p')/(q t+q'),
    p/q-theta=1/[q(qt+q')].

They follow by multiplying the two-by-two matrices [[a,1],[1,0]] for the
partial quotients; the determinant alternates sign. Write a=floor(t)>=1
and choose

    m=(a+1)q+q'.

Then mp=1 modulo q, m>2q, and

    m-q=aq+q' < qt+q' < (a+1)q+q'=m,

which is exactly (7). Irrational continued fractions have infinitely many
upper convergents with denominators tending to infinity: the recurrence
q_next=a_next q+q_prev makes them strictly increase after the initial
terms, and a finite expansion would be rational. The two-by-two determinant
identity also makes each convergent reduced.

Hence (8) gives an unbounded sequence of fully lawful first-crossing word
families with unbounded explicit-reset counts floor(m theta)+1. For all
sufficiently large members, the greedy note proves that every consecutive
pure-expansion power is below the previous one-percent threshold, even
under K^4F<->V rewriting. The new coverage comes from the past 3-adic image
guard combined with the repeated future, not from increasing an episode cap
or reusing a one-episode repair theorem.

## 7. Exact controls and scope

The companion certificate checks two periodic-tail instances. The large
one has N=1833237826965784108483587883633758044204176, well beyond 2^52.
All verification uses rational logarithm enclosures, integer arithmetic,
short bank words and compressed periodic descriptions. Its giant word,
branch numerator and orbit are not materialized.

For that large instance every long equal factor has a contracting separator
and identical future continuation, as proved in the hostile note. The new
image guard still excludes an actual non-descending endpoint. This sharply
separates a failed purely formal shortcut from a successful arithmetic test.

The general S_d boundary may still contain other realized repetition-poor
sources. No theorem here supplies arbitrary entry into this bank language,
excludes every possible positive cycle, or proves the Collatz conjecture.

## 8. Precise next target: typed left-special factors in every greedy prefix

The proved unbounded subsequence does not settle all greedy W_m. A sharper
open target is:

> Is there an m_0 such that every m>=m_0 has two macroblock checkpoints
> carrying different incoming 3-adic types, with identical following r-bit
> parity factors wholly inside W_m, where r>=ceil(122m/100)?

A positive answer, together with the finite short-range closure, would give
strict endpoint descent for the entire greedy family. More generally one
can replace the different types by distinct certified past-image residues
modulo a power of three, and use the quantitative CRT budget if the match
is shorter than the displayed scale.

For this coding, the main checkpoint types are 4 after V and 5 after a
following explicit F (a final consecutive F can give type 6). Finding a
shared future with distinct such pasts is a typed left-special-factor
problem for the mechanical reset sequence and its variable-length literal
encoding, with the actual initial V and final forced F retained.

This is a finite-prefix recurrence question. The infinite irrational
rotation may supply left-special factors eventually without placing both
required occurrences, and enough of their futures, inside the present
word. Large partial quotients make that distinction important. Neither
irrationality, the global factor-complexity estimate, nor a formal infinite
extension proves the stated finite occurrence hypothesis. The target is
recorded as open, not used in any theorem above.
