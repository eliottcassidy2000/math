# Carry-length fingerprints and an explicit negative-context family

**Status: proved, scoped to the inherited balanced two-apex chart.** The input is `apex_wedge_constructor_20261103/{THEORY,PROFILE,CARRIER}.md` on current main (commit `17b4acad2ee77e9bcec299f9fcba816a7b2a0cf7`); the earlier path-cover identity is used as a theorem, not reasserted here as new. No universal absolute-count realization is claimed.

## The carry-length theorem

For integer P≥1 let R=P+2 and suppose P,R both have k binary digits. Form the actual labelled strong tournament pair T−,T+ from the inherited balanced chart p01s,p10s. Let t be the number of consecutive 1 bits of P beginning with bit 1 (the second-least-significant bit). The leading bit cannot be consumed by the carry, so 0≤t≤k−3. Let Δ_r=F_r(T+)−F_r(T−) for ordered r-path covers, and Δ_0=0. Then

    Δ_r = 2 r^(2k−t−1) Q_t(r),
    Q_t(r) = r^(t+1) + (r²−1)((r−1)^t−(r+1)^t)/2.             (1)

The displayed half-expression is integral as a polynomial. For all such P the ordinary Hamiltonian jump is Δ_1=2, regardless of t. After the inherited strong triangle substitution G, the same labelled edge reversal has Hamiltonian jump

    Δ'_1 = Δ_1+Δ_2
         = 2+2^(2k−t) [2^(t+1)+(3−3^(t+1))/2].               (2)

For t=0,1 this is positive. For every t≥2 it is strictly negative. Every t≥0 actually occurs: take k≥t+3 and P=2^(k−1)+2^(t+1)−2. Thus arbitrarily long binary carries yield infinitely many literal, strong pairs with the *same +2 original jump*, whose jump under one fixed strong attachment changes sign. This sharpens the single sign-change example from the input and explains its cause. This does not assert the original sources have equal H.

**Proof.** The inherited PROFILE formula states Δ_r=2r(r^k C_r(R)−r^k C_r(P)), with C_r(N)=Σ n_i r^i L_(k−1−i), L_h=((r+1)^h+(r−1)^h)/2. Addition of 2 clears exactly the t ones in positions 1,...,t and sets the preceding zero at position t+1. Therefore

    C_r(R)−C_r(P)
      = r^(k−t−2) [L_(t+1)−Σ_(i=1)^t r^(t+1−i)L_i].

Summing the two geometric series using the roots r±1 gives the bracket Q_t, proving (1). At r=1, Q_t=1. At r=2, Q_t(2)=2^(t+1)+(3−3^(t+1))/2, giving (2) through Δ'_1=Δ_1+Δ_2. The quantity Q_2(2)=−4, and Q_(t+1)(2)−Q_t(2)=2^(t+1)−3^(t+1)<0 for t≥1, proving the sign. The construction has t carry ones, keeps the leading bit, and satisfies the inherited chart's strongness guard. □

## A polynomial fingerprint of carry length

For fixed k the entire path-cover response determines t exactly. Indeed ord_0 Δ=2k for t=0. If t=2s−1 or t=2s with s≥1, then ord_0 Δ=2k−2s in both cases, but its coefficient of that lowest power is respectively 2 or 4s. Thus t is recovered without knowing P beyond k. Conversely (1) shows that, for fixed k and t, every other prefix bit and the last bit are invisible to the entire edge-response polynomial. This is equality of differences only, not equality of source counts or of tournaments. The two suffix-probe language from CARRIER can recover the normalized word because it allows *further appended columns*; a single full-profile response at the current source has a strictly coarser fingerprint. For instance P=8 and P=9 have the same k,t, hence the same Δ_r, although their literal words differ; their lower/upper counts are 4939/4941 and 10345/10347, respectively.

**Proof.** For odd t the constant term of Q_t is 1. For positive even t its constant term is 0 and its linear term is t; for t=0, Q_0=r. Substituting these lowest terms into (1) gives the stated valuations and coefficients. □

## Exact starting counts for a new infinite family

For every j≥4 choose P=2^j+6, R=2^j+8. Both have k=j+1 bits and t=2. The source has 2j+6 vertices and

    H(T−) = (26665/64) 4^j + (31479/16) 2^j + 1721 − 191 [j=4],
    H(T+) = H(T−)+2,                                            (3)
    H(G(T+))−H(G(T−)) = 2−2^(2j+2).                           (4)

Thus (3) specifies the starting counts explicitly, rather than prescribing an unattained arbitrary H. The attachment retains the named old vertices and apex edge. At j=4 the literal words are 101100111100 and 101101011100, the counts are 139669 and 139671, and the attached counts are 11830709 and 11829687. They give 2 and −1022, respectively. At j=5 the lower starting count is 491319, and the post-attachment jump is −4094.

**Exact recurrence certificate for (3).** The incoming 27-state gap automaton counts H for each apex order. Encode each paired bit column (c,d) by c+2d. For j≥4 the exact literal word-pair template is

    3 0^(j−3) 330213330 3^(j−4) 0.                           (5)

Using the input automaton, with A and B the transition matrices for letters 0 and 3, respectively, the two apex-order counts have the forms α A^n M B^n β, n=j−4. Each left and right Krylov space has dimensions 2 and 3 for both orders. Their direct-sum product representation has dimension 12. The integer/rational code `verify.py` builds this representation, evaluates p(E) with p(x)=x(x−1)(x−2)(x−4), and checks the exact 12 residual moments α E^q p(E)β=0 for q=0,...,11. Cayley–Hamilton then proves the annihilator for the sum of the two counts for all n≥0. Their difference is the constant 2, which is also annihilated by p(E); hence each individual count obeys p(E). Both sides of (3) therefore obey p(E); equality at j=4,5,6,7 proves equality for all j≥4. The subtraction 191 is the necessary 0-root transient at n=0. Reproduction additionally checks (5), the actual gap counts through j=12, and literal graph subset DP at j=4, including both attached graphs. This is a finite exact linear certificate, not a fit inferred from samples.

The same script checks (1) symbolically for all 494 same-length balanced pairs with P<512. These finite computations test the symbolic construction but are not needed to make its quantifiers all-order. The script imports unmodified incoming `mechanisms.py`, `profile.py`, and `inherited_recurrences.py` and uses Python 3 and SymPy.
