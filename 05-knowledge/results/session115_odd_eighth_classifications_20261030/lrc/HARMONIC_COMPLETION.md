# Harmonic completion: exact families, owner counts, and odd-eighth phases

Session 115, 30 October 2026. **PROVED scoped.** The full Lonely Runner Conjecture remains open. No literature-priority claim.

## Inherited mechanism and the new operation

[Session114's complete core theorem](../../session114_continuation_symmetry_20261030/lrc/GENERAL_CORE.md) proves the complete top-window description for

 C_q={1,...,q−1,q+1}, q>=3,

above radius1/(q+2). At each primitive a/q there are two directional owners. Their multiset contains each element of

 E_q={e:2<=e<=q−1, gcd(e,q)=1} union {q+1}

twice, once on each reflected side. Put s=s(q)=min E_q, the least prime not dividing q. At radiusr the side with owner e extends distance(1/q−r)/e from its center. The complete window statement, not only the maximizing grid, is inherited.

We now append a harmonic block of resonant tails:

 V(q,m,k)=C_q union {qm,2qm,...,kqm},  m,k>=1.

This preserves the named core and changes the continuation operation from one tail to a whole AP of tails. The least-owner side supplies the widest legal displacement. All tails vanish at every a/q, so their new observable depends only on that displacement.

## 1. Exact harmonic-block spectrum in the top chamber

**Theorem1.** If 1<=k<=q and 2m>s(q), then

 M(V(q,m,k))=min{m/(qm+s(q)), 1/(k+1)}.                 (1)

In particular, every block of at mostq−1 tails has the same exact maximum as its first tail. For the complete q-tail block the formula is m/(qm+s) when s/2<m<s, and1/(q+1) when m>=s.

**Proof.** Write t=a/q+x in a surviving core window, set y=q|x| and let e be its directional owner. Core safety at radiusr>1/(q+2) implies

 r <= (1−e y)/q <= (1−s y)/q.

The first tail gives r<=||m y||<=m y, so r<=m/(qm+s). The AP of k tails gives the other universal bound1/(k+1), by the usual k+1-point gap argument.

Let r be the right side of(1). The hypotheses make r>1/(q+2). Choose the side with owner s and displacement y=r/m. Since r<=m/(qm+s), this lies in its core window. Every tailj has phase distance||j r||>=r because r<=1/(k+1) and r<=jr<=kr<=1−r. Thus r is attained, proving both upper and lower bounds. This is a direct source-preserving construction; it uses no independent-phase model. QED.

When the first bound is strict, equality in the two inequalities forces the two reflected s-owner endpoints. Thus there are exactly two maximizers, with active speeds s andqm, and terminal mass slope

 B=2(1/s+1/(qm)).                                      (2)

At the collision between the two bounds, (2) is generally false: another tail becomes a faster binding owner. The next theorem retains this collision exactly.

## 2. Complete classification of harmonic-top preservation

Set k=q and n=q+1. The q-tail block alone has maximum1/n.

**Theorem2.** For every q>=3 and every m>=1,

 M(V(q,m,q))=1/(q+1)  if and only if  m>=s(q).          (3)

For the entire preserving family m>=s(q), all maximizers are precisely

 t_{e,j}=a_e/q − j/[q(q+1)m],  and 1−t_{e,j},           (4)

where e ranges over E_q, a_e is the inverse of e moduloq in{1,...,q−1}, and

 1<=j<=floor(m/e),   gcd(j,q+1)=1.

All listed points are distinct. Their number is exactly

 N(q,m)=2 sum_(e in E_q) sum_(d | q+1) mu(d) floor(m/(ed)).       (5)

**Proof.** The q-tail AP can attain1/n only when its base phase is a primitive multiple of1/n. At a core window around a/q this is exactly

 qm x = ±j/n,  gcd(j,n)=1.

Its nonzero displacement has magnitude at least1/(qnm). The largest core side at radius1/n has length1/(qns), so necessarily m>=s. Conversely the e=s,j=1 choice is safe if m>=s, proving(3). A side with owner e admits precisely those positive j with ej<=m, giving(4). The top windows are disjoint and no listed point is a center, so no duplicates occur. Every directional owner e occurs once as a left owner and once as its reflection. Counting and ordinary Möbius inversion give(5). QED.

This completely classifies preservation of the harmonic block's optimal value for all m, and completely classifies the maximizing times within the preserving family. For m<=s/2, the theorem says only that the maximum is below1/(q+1), and in fact at most1/(q+2). It does not supply a complete spectrum in that lower chamber.

## 3. Exact physical boundary weights

For a unit j modulo n=q+1 let ell_j be its inverse in{1,...,n−1}, and define

 w_n(j)=1/ell_j + 1/(n−ell_j).

**Theorem3.** Throughout the preserving family m>=s(q), the terminal mass slope is

 B(q,m)= [2/(qm)] sum_(e in E_q) sum_(1<=j<=floor(m/e), gcd(j,n)=1) w_n(j).       (6)

Here mu S(M−epsilon)=epsilon B for sufficiently small positive epsilon, with closed maximizing points retained.

**Proof.** At each listed point the two binding tails have indices ell_j and n−ell_j, up to reflection. Their speeds areqm ell_j andqm(n−ell_j). A core constraint can also bind when ej=m, but its speed e is smaller thanqm, hence smaller than either tail speed. Physical relaxation resolves a same-side tie by the maximum owner speed, so the core does not change this point's reciprocal-owner weight. Its contribution is w_n(j)/(qm); reflection doubles it. Summing proves(6). QED.

In particular, at m=s only the reflected e=s,j=1 pair survives, and

 B(q,s)=2(q+1)/(q^2 s).                                 (7)

This is the exact critical-owner correction to the invalid use of(2) at equality. For q=7,m=2, (2) would give8/7 whereas the correct answer is8/49: the newly active last tail is seven times faster than the old core owner.

The count is a quasipolynomial of degree1. More precisely, let

 L=(q+1) lcm(E_q), H_E=sum_(e in E_q)1/e,
 H_U=sum_(1<=a<=q, gcd(a,q+1)=1)1/a.

Then for every m>=s,

 N(q,m+L)−N(q,m)=2L phi(q+1) H_E/(q+1),                 (8)
 (m+L) B(q,m+L)−m B(q,m)=4L H_E H_U/[q(q+1)].           (9)

These identities follow by adding complete unit-residue periods in each inner sum. Consequently

 N(q,m)/m -> 2 phi(q+1) H_E/(q+1),
 B(q,m) -> 4 H_E H_U/[q(q+1)].                          (10)

The error in the first limit after subtracting its linear term is bounded; the second error is O_q(1/m). These are exact weighted incidence formulas, not assumptions of uniform sampling.

## 4. The odd-eighth family is a real phase carrier

Take q=7. Then n=8, s=2 and E_7={2,3,4,5,6,8}. For every m>=2, the complete seven-tail family

 {1,2,3,4,5,6,8} union {7m,14m,...,49m}

has maximum1/8. Its maximizing tail-base phases are drawn from

 1/8, 3/8, 5/8, 7/8,

because these are exactly the primitive8-grid phases of the AP7 tail block. All four occur once m>=6; for2<=m<=5 only the1/8 and7/8 classes occur. This threshold follows from the first allowed positive j=3 needing m>=3 min E_7=6.

The endpoint weights differ:

 phase1/8 or7/8: 8/(49m),
 phase3/8 or5/8: 8/(105m).                              (11)

Thus multiplication by a unit modulo8 acts transitively on the abstract phase set but does not preserve its physical relaxation weights. The owner inverse is the necessary sidecar. There is no map here identifying these numbers with a Collatz content rate merely because its coefficient is5/8.

The weight distortion has an exact quadratic-character decomposition. Define the real primitive character modulo8 by

 chi_2(j)=(-1)^((j^2−1)/8), for j odd.

It is +1 on residues1,7 and−1 on3,5 (the usual quadratic character of2). On the actual odd-eighth phase labels t=j/8, it also equals sqrt(2) cos(2 pi t). This identifies the correction with a real eighth-root character, not with an unrelated numerical occurrence of the fraction5/8. Directly from the four inverse residues,

 w_8(j)=88/105+(32/105) chi_2(j).                        (12)

Let A(m)=sum_(e in E_7) sum_(odd j<=floor(m/e)) chi_2(j). Then

 B(7,m)=88 N(7,m)/(735m)+64 A(m)/(735m).                 (13)

This separates the unweighted count from the missing signed-character tally. The character sum on each initial odd interval is0,+1 or−1, so |A(m)|<=6; A is periodic with period960. Consequently (13) is also an explicit O(1/m) correction to a count-only formula. Multiplication by3 swaps the two weight classes, while reflection preserves them. A single per-phase sign chi_2 is necessary and sufficient to recover its weight; it does not replace the complete displacement index needed to recover a physical point. This is a genuine map from the odd-eighth phase carrier to a quadratic-character observable, not an analogy based only on numerators.

The exact count specializes to

 N(7,m)=2 sum_(e in {2,3,4,5,6,8}) ceil(floor(m/e)/2),

and the asymptotic constants are

 N(7,m)/m -> 63/40,   B(7,m) -> 33/175.

This is a concrete transfer: the AP7 unit-phase carrier is pulled back through the widest complete core windows, with direction, displacement index, and reciprocal tail-owner weights retained. A bare four-state quotient loses both the multiplicities in(5) and the weights in(6).

## Inheritance, tests, and limits

The underlying AP maximum, Farey windows, max-owner derivative, and Goddyn–Wong tightness lineage are inherited from Session114 and its references. The new operation is harmonic completion; the new classification is its threshold, complete maximizing-time census, exact weighted slope and affine-periodic floor laws. No external priority claim is made.

The hostile cases are: m<s (no harmonic-top preservation), m=s (new fastest owner), a nonunit j (not an AP maximizer), and modular phase relabeling without derivative transport. Exact checks are in `check_harmonic.py` and its receipt. Arbitrary cores, arbitrary unrelated tail sets, and the full low-clearance chamber remain outside the theorem.
