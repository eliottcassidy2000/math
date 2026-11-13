# Random colors give predictable symmetry growth without altering counts

In the [binary wreath atlas](WREATH_ATLAS.md), color each leaf independently B or D with probability 1/2. Every outcome has the same full path-cover polynomial. If a height-j class X has group exponent e(X), its exact probability is

    p_j(X)=3^((3^j−1)/2−e(X))/2^(3^j),

by orbit–stabilizer for the recursively rotating tree. In particular the probability that three independent height-j child classes all agree is

    s_j=Σ_X p_j(X)^3
       = 3^(3(3^j−1)/2) · 2^(−3·3^j) · A_j(3^(−3)),

with the inventory polynomial A_j of the atlas. Thus if E_d is the random exponent of the automorphism order at height d, linearity gives

    E[E_d]/3^d = Σ_(j=0)^(d−1) s_j / 3^(j+1).

The first terms are s_0=1/4, s_1=7/64; higher terms follow from the exact inventory. The limit κ exists and is positive. In fact p_j(X)≤ 3^((3^j−1)/2)/2^(3^j), so s_j≤max_X p_j(X)^2; the tail in κ decays doubly exponentially in the depth j.

There is also concentration. Flipping one leaf changes the exponent by at most d, proved in the atlas. The bounded-differences inequality gives for every ε>0

    Pr(|E_d−E[E_d]|≥ε·3^d) ≤ 2 exp(−2 ε²·3^d/d²).

For a short derivation, reveal the 3^d independent leaf colors one by one and take conditional expectations of E_d. The two possible next conditional expectations differ by at most d, by coupling all unrevealed colors. Hence each centered martingale increment has conditional exponential moment at most exp(t²d²/8) (the elementary exponential-moment bound for an interval of length d, by convexity), and multiplying, then optimizing the resulting exponential Markov bound in t, gives the display. In particular E_d/3^d converges in probability to κ. Random independent colorings do not tend to rigidity: already at the bottom layer, there are independently 3^(d−1) monochromatic-triple tests, each succeeding with probability 1/4. On the other hand a neutral repainting can force either rigidity or the full group. This compares a typical construction with both extremes without conflating uniform sampling of labeled leaves and uniform sampling of isomorphism classes.
