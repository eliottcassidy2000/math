# A sevenfold target directly from a norm reader

The [earlier norm forms](NORM_READERS.md) are not limited to transitive incoming blocks, and can produce a prime's sevenfold target without applying an external ×7 reader. Let u=B−h for the incoming block. The two indefinite forms Q=(2h+3u)^2−2u^2 and Q'=6h²+16hu+7u² in that note give:

| Incoming masks | (h,B) | Reader | Literal H |
|---|---|---|---:|
| (30,12,24,16,2) | (5,64) | Q | 28,007 = 7·4,001 |
| (126,92,24,112,96,70,4) | (29,640) | Q | 2,829,239 = 7·404,177 |
| (126,76,56,112,98,66,4) | (37,752) | Q' | 4,010,069 = 7·572,867 |

All three quotients are prime by trial division. For the first target the ordinary side has the separate order-10 literal witness (460,25,426,400,997,331,14,98,194,495). No minimality is claimed.

The construction has a reusable modular component. On inputs with 7∤h, divisibility by 7 is equivalent to u≡2h for Q, or u≡4h for Q'. Indeed Q≡h(4h+5u) and Q'≡2h(3h+u) modulo 7. The three transitive ratios u/h=2^k−1 run through 0,1,3 modulo 7 and miss both branches. This is one reason to change the source rather than merely change the width. Once a nontransitive source realizes either branch with valuation exactly one, append any multiple of 21 ordinal singletons: h is unchanged and B is multiplied by 2^21≡1 modulo 49 per step. Thus both congruence conditions survive at arbitrarily large actual outputs, each having 7-adic valuation exactly one and partial quotient fixed modulo 7. This is an infinite exact-value family, not a claim of infinitely many prime quotients and not an inverse for every 7p.

For a prescribed p one can instead solve Q=7p or Q'=7p in the physical cone and query the inverse profile (h,B). At p=4,001, x=187 and u=59 solve x²−2u²=28,007, requiring h=5,B=64, realized by the marked input above. The same norm target equation tells us which correlated profile to seek, without requiring p itself to have the same representation.

In general a two-copy, two-visit quotient reads x h² + y h v + z v², where v=B−2h is the ordered nonempty two-cover count, x counts neither block visited twice, y counts exactly one visited twice, and z counts both visited twice. The discriminant y²−4xz is entirely quotient data; the linear change from v to u changes no discriminant. Indefiniteness says the single-double-visit patterns dominate the geometric mean of the other two coefficients, while actual physical evaluations stay positive. This separates asking for a new form from asking for a new input point on it.

Replay: `python3 norm_seven.py`.
