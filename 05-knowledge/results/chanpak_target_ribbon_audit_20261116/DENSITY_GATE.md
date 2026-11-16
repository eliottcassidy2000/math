# What would certify a power saving for missing span-3 counts?

Let S be the spectrum of literal span-3 ribbons, and in a dyadic interval (X,2X] let A be the admissible odd targets and L=|A| (excluding the finitely many fixed obstructions if they occur there). Write r(n) for the number of *certified physical lifts* of n in any selected finite multiset of ribbon words, and R=Σ r(n).

The exact missing-value gate is:
```
E := Σ_{n∈A} (r(n)−R/L)^2 = Σ r(n)^2 − R^2/L.
#{n∈A:r(n)=0} ≤ E/(R/L)^2   (R>0).
```
Every missing target contributes the full squared mean to E, so the inequality is immediate. If exact absolute collisions of the selected words could be bounded by `Σ r(n)^2 ≤ (1+O(X^(−δ))) R^2/L`, the needed dyadic power saving would follow. Summing intervals would yield an O(X^(1−δ)) exceptional set for 0<δ<1. A constant-factor collision bound gives only a positive proportion of targets, not density one. Residue collisions alone are not the absolute collision estimate here unless the carry or interval has separately been controlled.

This clarifies the role of a rank-dependent lower-growth certificate. Its conclusion is **length control on physical words**; for example an H≥Cρ^n estimate bounds the lengths of small-count representatives where it applies. It does not imply that different words yield different counts or that the values are equidistributed. Even an exponential-growth lower bound for an entire selected language can coexist with a single geometric sequence of outputs. Nor does an upper growth bound H(wa)≤K H(w) remove this gap: it gives a multiplicative bracket, not arbitrary exact counts. The modular-to-absolute carry loss from the seam construction is one literal example of this logical difference.

The power saving requested is not established by the results in this packet. There is a verified lower-growth certificate in the separate language of thirteen floor syllables; passing to all literal span-3 words, proving exact collision dispersion, and selecting every individual target are three distinct tasks.

The [positive frontier envelope](POSITIVE_ENVELOPE.md) does give one rigorous corrected counting input: every word of length n satisfies H≤5(39/10)^n. There are 8^(n−2) distinct in-range span-3 words at order n≥3 (the first three positions have 1, 2, and 4 effective choices). Taking n=floor(log(X/5)/log(39/10)), for large X these alone provide Ω(X^(log 8/log(39/10))) distinct labelled natural-order ribbons, a word exponent about 1.528. This is not a bound on distinct counts, and does not amplify the restricted lower-growth certificate. [The source audit](DENSITY_AUDIT.md) explains why a finite sampled lower bound on the joint spectral radius could not be used instead of a certified upper one.
