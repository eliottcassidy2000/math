# Changing the primitive projector: three more shield interfaces

The [unit shield](SHIELD.md) has only two right-hand rays for its
fixed projector. The projector itself can change. Here are three
additional all-height physical interfaces, including two with
dyadic denominator 16384 rather than 128.

Choose the suffix y and the two left words from this table. Use the
same reverse bridge, and use the right word v=rc(x y), where x is
either 11 0^t or 101 0^t, t>=1.

| y | k | p | u0 | u1 | physical edits |
|--|--:|--:|--|--|--:|
| 1010010 | 7 | 91 | 1101001111100 | 1101100110010 | 5 |
| 11010010011100 | 14 | 12725 | 110101001001001000 | 111011010110110010 | 12 |
| 11101010011000 | 14 | 14761 | 101010011110000100 | 110001101011010110 | 10 |

For each y write a=A(y). Each suffix starts in 1, has
p=2^k−rev_k(a), and satisfies p(a+1)−2^k(h(y)−2)=1.
Thus for either prefix of defect 2, with A(x)=m2^t (m=3 or 5),

    h(v) = p m 2^t + h(y)−2,
    Q = 2^k m 2^t + a+1,
    −2^k h(v)+pQ = 1.

For each left pair the recurrence directly gives
(Δh,ΔP)=2(−2^k,p), so exact +2 follows for all t. The right
words are strong. Deleting both apices still leaves a transitive
spine, and only the chosen left apex contacts change. The two
baseline rays for each row, with D=2^t, are respectively:

| y | H0, m=3 | H0, m=5 |
|--|--|--|
| 1010010 | 4004667 D+865551 | 6674445 D+865551 |
| 11010010011100 | 18415108821 D+5046246339 | 30691848035 D+5046246339 |
| 11101010011000 | 15233834805 D+4649307361 | 25389724675 D+4649307361 |

These are constructive choices, not an exhaustive search beyond the
specified projectors. Together with the main pair, they produce eight
rays of exact adjacent odd counts. A finite union remains sparse;
getting dense absolute windows requires a growing bank of compatible
left profiles and right projectors, or a seam with more freedom.

The same structural test also provides a finite left-profile selector:
precompute h(A) using the simple one-apex recurrence, and ask whether
for the prescribed suffix there is an even normalized address A with
h(A+2p)−h(A)=−2^(k+1). It is a single exact profile-difference problem,
not selection by residue or by the unconstrained ratio Q/h.
