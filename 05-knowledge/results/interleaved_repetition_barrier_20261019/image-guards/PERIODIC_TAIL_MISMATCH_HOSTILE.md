# A certified periodic-tail hostile to the mismatch shortcut

19 October 2026. **PROVED finite structural hostile and general conditional
construction.** The companion IMAGE_GUARD_DESCENT.md now excludes non-descending endpoints
for this whole conditional family, using an additional 3-adic image guard.
The smaller control is already excluded arithmetically by the primary
short-range theorem; Section 5 certifies a second instance far beyond that
range. Their formal long factors cannot remove the matched-block theorem's
exact-cycle alternative.

Use the greedy family W_m and notation v=s_V, mu=3/4,
theta=log(v)/log(4/3) from the separate extension notes.

## 1. Exact rational certificate

Set

    p=1003, q=9885, m=1116187,
    delta=p/q-theta.

Exact rational atanh-series enclosures for log2 and log3 certify

    gcd(p,q)=1, mp=1 mod q, m>2q,
    1/(mq)<delta<1/((m-q)q).                     (1)

These are certified inequalities about the specific logarithmic ratio;
no floating-point continued-fraction claim is used. The companion checker
uses 512 terms for each logarithm, sufficient for both controls, and verifies
the inequalities by integer cross multiplication. It also certifies floor(m theta)=113255.

The following proof applies to any 0<p<q with q>=2 satisfying (1).

## 2. The final forced reset restores a periodic tail

For n>=0 put

    H_n=floor(np/q)-1_{q divides n},
    g_n=H_n-H_(n-1) for n>=1.

Then H_0=-1, H_(n+q)=H_n+p, and g_n is a q-periodic binary sequence with
exactly p ones per period. It is the lower rational mechanical sequence
with the boundary pair 10 replaced by 01. In particular g_1=1.

For 1<=n<m, the true irrational floors obey

    floor(n theta)=H_n.                          (2)

At n=m there is precisely one extra loss:

    floor(m theta)=H_m-1.                        (3)

Here is the complete residue argument. Write np modulo q as k. Since
n delta<2/q, residues k>=2 never cross their integer boundary. At k=0,
the positive displacement n delta lowers the integer floor by exactly
one, as recorded in H_n. Residue k=1 occurs at n=m modulo q. Its last
occurrence before m is m-q, and (1) puts n delta<1/q there and at every
earlier occurrence. At m, (1) gives m delta>1/q. This proves (2)--(3).

Let f_n=floor(n theta)-floor((n-1)theta) be the greedy reset indicators.
The discrepancy H_0=-1 instead of floor(0 theta)=0 gives f_1=0 instead
of g_1=1. For 2<=n<m, f_n=g_n. Equations (2)--(3) give f_m=g_m-1;
since f_m is a nonnegative binary value and g_m is binary, g_m=1 and
f_m=0. The extra final F in W_m restores g_m=1.

Under the literal encoding 0->V, 1->VF, therefore,

    W_m is a prefix of A B B B ...,
    A=V,
    B=encode(g_2 g_3 ... g_q g_1).                (4)

This is a statement about the entire finite word, including the final F,
not an inference from the irrationality of an infinite rotation.

The preperiod has a=122 bits. A period has

    P=122q+2p=1207976 bits,
    h_B=77q+p=762148 odd bits.

Its affine multiplier is strictly contracting:

    rho=v^q mu^p,
    log rho=-q delta log(4/3)<0.                  (5)

For the concrete m,

    L=122m=136174814,
    N=L+2*(113255+1)=136401326,
    r=ceil(L/100)=1361749 > a+P=1208098.          (6)

## 3. No long matched factors have a divergent continuation

General elementary lemma: if a finite word is a prefix of A B-infinity,
with |A|=a and |B|=P, and r>=a+P, then any two equal length-r factors have
identical entire future continuations in that infinite extension.

Indeed, after at most a aligned positions, both continuations are in the
P-periodic tail. Their common length-r prefix contains at least P aligned
positions in that tail. Agreement on one full period implies agreement at
every later aligned position. All earlier positions were already included
in their matched prefixes.

Thus the matched factors in (6) cannot have a later witnessed parity
mismatch anywhere in W_m. The stronger noncontracting-separator test also
fails. If their starts are i<j, equality of their infinite continuations
makes j-i a period of the tail from i onward. That period has the same
odd-bit frequency h_B/P as B. The intervening multiplier consequently is

    exp(((j-i)/P) log rho)<1.                     (7)

This argument does not assume B is a primitive period. In particular,
all equal r-factor pairs have contracting separators, so neither a future
mismatch nor a separator of multiplier at least one removes the exact
return alternative.

There are plenty of matching pairs: (4) provides starts a and a+P with
identical r-bit factors wholly inside W_m, by (6). The periodic-tail
representation is not merely a claim that a factor-count test failed to
find a match.

## 4. What remains arithmetic, and what this example does not prove

Every positive integer lift of the unique lawful residue of W_m follows
the whole displayed first-crossing word. A nonnegative endpoint displacement
still requires, with the full branch (A_m x+B_m)/M_m,

    x=(B_m-M_m d)/(M_m-A_m) in the positive integers.

For this concrete word N<2^52, so the primary exact short-range theorem
already proves strict endpoint descent and rules out all such d>=0.
There is no actual cycle or counterexample hidden in this structural
hostile. It proves that a proposed *purely formal* mismatch argument cannot
cover every greedy word.

The conditional construction in Sections 1--3 would yield arbitrarily
large hostile words if one could prove infinitely many rational brackets
(1) with r>=122+122q+2p. That infinitude is not established here. One large
explicit word does not refute an eventual mismatch theorem beyond an
unspecified finite threshold. Any such theorem, however, needs a genuine
argument about this logarithmic rotation or about endpoint realizability;
irrationality alone cannot justify deleting the internal-cycle alternative.

## 5. A certified instance beyond the short-range theorem

The same exact certificate can be scaled far beyond the proved short range.
Use

    p=8116910408926339217391896275524023351,
    q=79995672446647253540458867489349070728,
    m=15001586057656515185428419394326237781535.

A 512-term rational logarithm enclosure verifies every inequality in (1),
without constructing the word, its affine numerator, or an orbit. It also
certifies the exact final floor

    floor(m theta)=1522163965844627930660358762978517428452.

The resulting dimensions are

    L=1830193499034094852622267166107801009347270,
    N=1833237826965784108483587883633758044204176,
    P=9775705859308817610370765626251634675518,
    r=18301934990340948526222671661078010093473.

In particular N>2^52 and r>122+P. The periodic tail contains 187 complete
periods, followed by a partial period. All conclusions of Sections 2--3
hold at these actual finite lengths: equal length-r factors have identical
future continuations and contracting separators. The earlier short-range
lemma cannot dispose of this instance merely from its length.

This is still a *formal* first-crossing word with fully realizable positive
lawful cylinders. The image-guard theorem now excludes a positive-displacement or zero-return
endpoint for this instance, independently of the short-range cutoff. The distinction is now genuinely beyond the existing
short-range exclusion, rather than a large-looking word that falls inside
that already-proved interval.

The rational candidate was found by interval-certified continued fractions,
but the certificate does not depend on its convergent status: its only
inputs are the displayed integers and the exact inequalities (1). This is
also why no unsupported asymptotic statement about the partial quotients
of theta enters the proof.

## 6. Exact reduction of endpoint realizability to a cycle equation

For either certified instance, write the compressed word as

    W_m=V B^k R,

where R is the appropriate literal prefix of B. The two starts immediately
after V and after VB support matching r-bit blocks wholly inside W_m.
Consequently the matched-block theorem implies that any non-descending
lawful source must satisfy

    y=T^122(x)=T^(122+P)(x).

Write the branch of B as (a_B y+b_B)/M_B. The exact necessary cycle equation
is therefore

    y=b_B/(M_B-a_B) in the positive integers.    (8)

In fact, for any fixed d>=0 the endpoint condition for this word is
*equivalent* to the following three tests:

    y=b_B/(M_B-a_B) is a positive integer;
    x=(m_V y-b_V)/a_V is a positive integer;
    T_R(y)=x+d.                                  (9)

Necessity is the matched-block theorem plus the actual V prefix and final
phase R. For sufficiency, (8) authenticates B's complete lawful cylinder
by endpoint integrality. The integer equation for V authenticates its
cylinder as well. Hence x follows V, then as many exact B cycles as needed,
then R, and the final equality in (9) is the required endpoint.

Thus the remaining barrier on this explicitly described, very long formal
family has a precise arithmetic form: a positive integral periodic point
for B, an integral V preimage, and the prescribed phase displacement.
These tests are in fact incompatible: the first checkpoint is a V output
(type 4), while returning around B would make it a VF output (type 5). The
companion IMAGE_GUARD_DESCENT.md proves this arithmetic exclusion. The purely
symbolic complexity and mismatch tests by themselves leave this cycle side
condition; the 3-adic past-image guard supplies the missing information.
