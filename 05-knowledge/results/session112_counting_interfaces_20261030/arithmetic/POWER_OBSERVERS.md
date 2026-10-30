# Weighted aliquot observers: recover a power, then account for its roots

**Status: PROVED** elementary observer and conjugacy statements; the unit-root count uses the classical unit-group structures stated below. Direct checks are in `verify_power_observers.py`. These results extend `ALIQUOT_MODULAR.md`; no factoring speedup or literature-priority claim is made.

## 1. The weighted observer

For an integer weight `w≥1`, define

    sigma_w(n)=sum_(d|n) d^w,       s_w(n)=sigma_w(n)−n^w.

Fix a prime `p`, step `t≥1`, precision `k≥1`, and put `q=p^w` and `S_j(q)=1+q+...+q^(j−1)`. For `n=p^e u`, `p∤u`, let `a=s_w(n)` and `b=s_w(p^t n)`. Multiplicativity gives

    b−q^t a=S_t(q) sigma_w(u),
    b−a=q^e S_t(q)[u^w+q s_w(u)],
    v_p(b−a)=we.

The bracket and `S_t(q)` are p-adic units. Consequently

    n^w=S_(e+1)(q) S_t(q)^(-1)(b−q^t a)−a.          (1)

This reconstructs `n^w`, not necessarily `n`.

For residues `A,B mod p^k`, put

    h=floor(min(v_p(B−A),k−1)/w),

where a zero residue is assigned valuation at least `k`. Then

    F_w(A,B)=S_(h+1)(q) S_t(q)^(-1)(B−q^t A)−A mod p^k    (2)

satisfies `F_w(s_w(n),s_w(p^t n))=n^w mod p^k` for every positive integer `n`. The proof is the same stabilization argument as the ordinary observer: `q^(h+1)` is divisible by `p^k` at the top block, and the geometric sum has stabilized there.

## 2. All future observations retain the same power coordinate

Formula (2) also defines an operation on every formal pair of residues, including pairs that cannot arise from weighted divisor sums. On blocks `w h≤v_p(B−A)<w(h+1)`, it is

    F_w(A,B)=U_h(B−A)−q^(h+1)A,
    U_h=S_(h+1)(q)/S_t(q),

with unit slope and a translation of strictly higher valuation than the difference. The same shell argument as in `ALIQUOT_MODULAR.md` proves it is an isometry and permutation in either observation separately.

For the state `(x,A)=(n^w,s_w(n)) mod p^k`, set `h=floor(min(v_p(x),k−1)/w)` and

    C=(A+x)/(1−q^(h+1)).                            (3)

The denominator is a unit; (3) is a bijective coordinate change on the entire formal residue plane. Substitution in the inverse observer gives the exact conjugacy

    (x,C) ↦ (q^t x,C).

On genuine sources, `C=sigma_w(u)/(1−q)`. Thus every future same-prime weighted observation is determined by the current pair `(n^w,s_w(n)) mod p^k`. Additional later observations cannot resolve a root ambiguity when two genuine sources already agree in both these coordinates.

## 3. The root information that remains after decoding

Let `K(p,w,l)` denote the number of unit solutions of `z^w=1 mod p^l`. The standard unit-group structures give

    K(p,w,l)=gcd(w,(p−1)p^(l−1))                if p is odd;
    K(2,w,1)=1;
    K(2,w,2)=gcd(w,2);
    K(2,w,l)=gcd(w,2) gcd(w,2^(l−2))           if l≥3.

The odd-prime unit group is cyclic of order `(p−1)p^(l−1)`; the dyadic group for `l≥3` is `C2 × C_(2^(l−2))`. On a cyclic group of order `m`, the kernel of taking the wth power has `gcd(w,m)` elements. These classical structural inputs are explained in Keith Conrad's [Prime-power units and finite subgroups of GL_n(Q)](https://kconrad.math.uconn.edu/blurbs/gradnumthy/primepowerunitsandGLnQ.pdf). For the dyadic factorization one may also check directly that −1 and 5 generate the two factors: repeated squaring gives `v_2(5^(2^j)−1)=j+2`, and powers of 5 are precisely the units congruent to 1 modulo 4.

If a nonzero decoded power `y` lies in the image and `v_p(y)=we<k`, its fibre under `n↦n^w mod p^k` has exactly

    p^((w−1)e) K(p,w,k−we)                         (4)

source residues. To prove this, write `n=p^e u`. The unit equation for `u` is at precision `p^(k−we)` and has `K` roots; each lifts freely to `u mod p^(k−e)` in `p^((w−1)e)` ways. The zero-power fibre has exactly `p^(k−ceil(k/w))` elements, since its members are the multiples of `p^ceil(k/w)`.

**Scope:** these are fibres of the decoded power, not a claim that all roots yield the same full genuine aliquot-observation pair. The retained first observation may distinguish roots. The next examples show that sometimes the full ambiguity really survives.

## 4. Genuine two-root and four-root ambiguities

For a prime source `n≠p`, `s_w(n)=1`, and

    s_w(p^t n)=S_(t+1)(q)+S_t(q)n^w.

Hence prime sources with equal wth powers modulo `p^k` have identical paired observations, and indeed identical observations at every position on the same multiplication ray.

- **Two roots:** take `w=2,p=5,k=2`. The distinct prime sources `2` and `23` represent opposite residues modulo 25 and both square to 4. Their observations for `t=1` are `(1,5)` modulo 25. Every later ray observation also agrees.
- **Four roots:** take `w=2,p=2,k=5`. Prime sources `97,47,17,31` occupy the four distinct residues `1,15,17,31` modulo 32, all squaring to 1. Their observations for `t=1` are all `(1,6)` modulo 32. Every later ray observation also agrees.

For unit squares, the exact root kernel therefore has size two at odd primes and four at dyadic precision at least three. Selecting a root in the latter example needs two binary labels. This is an arithmetic root-selection ambiguity; it is **not** the same group action as Session 111's two lost orientation signs, and no isomorphism between those mechanisms is claimed.

The Gaussian square/doubling work likewise needs representation or root data beyond a norm. The precise common lesson is the need to retain the fibre coordinate of the actual power map, rather than identifying objects from equal scalar observations.

## 5. Two coprime weights resolve the unit root gauge

Suppose `p∤n`, and paired observers at two weights `w,z` recover `x=n^w` and `y=n^z` modulo `p^k`. If `gcd(w,z)=1`, choose integers `alpha,beta` with `alpha w+beta z=1`. Since `x,y` are units,

    n=x^alpha y^beta mod p^k.

Negative exponents mean modular inverses. More generally the joint unit-power map has kernel `K(p,gcd(w,z),k)`, because the intersection of the two power kernels equals the kernel of their gcd by Bézout's identity.

In particular weights 2 and 3 give

    n=(decoded cube)/(decoded square) mod p^k.

The unit guard is essential. For instance modulo 25, sources 5 and 10 both have decoded square and cube zero, while their source residues differ. This last example concerns the decoded powers alone; it does not assert equality of all four original weighted aliquot observations. No unguarded division by a nonunit is allowed.

## Verification and limits

The checker uses direct divisor enumeration for weighted sums, exhaustive root fibres, the explicit prime-source witnesses, and Bézout recovery on units. Universal claims have the proofs above. The residue operation does not make divisor sums inexpensive to evaluate, nor does an equality of powers identify their underlying Gaussian representations. The root counts describe precisely where power decoding stops and what an additional coprime-weight observation supplies.
