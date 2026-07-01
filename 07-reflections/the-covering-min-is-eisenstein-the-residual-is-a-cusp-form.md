# The covering-min is Eisenstein; the residual is a cusp form

*klein-2026-06-30-S56. A reflection on HYP-3768 — bringing B2/E2/Dedekind to the hard residual, and finding the residual already had a name.*

For weeks the covering-min construction `{1,…,n-2, n(n-1)}` has felt like the *answer* to the wrong
question — the clean object that is nearly extremal but not quite (HYP-3764: it is the loose Farey rung).
This session, asked to think with Dedekind and E2, I computed its Dedekind sum, and it fell out in closed
form:
```
    s(n, Phi6(n)) = -(Phi6-1)/(12 Phi6)  ->  -1/12 = zeta(-1) = -B2/2.
```
The proof is one line of reciprocity, and it works for one reason: `Phi6 = n^2-n+1 ≡ 1 (mod n)`. That is
the same Farey-neighbor congruence that puts the covering-min on the Stern-Brocot ladder (HYP-3732).
And the limit, `-1/12`, is not a random constant — it is `zeta(-1)`, the constant term of the weight-2
Eisenstein series `E2`, the anomaly that makes `E2` fail to be modular. The covering-min construction
*is* the Eisenstein object. Its arithmetic is the E2 anomaly.

That reframes the whole covering-min story in the language the modular-forms side has been speaking all
along (HYP-3587, HYP-3041). On `X_0(2p)` the space of weight-2 forms splits into an **Eisenstein bulk**
(dimension `#cusps - 1 = 3`, constant across the family) and a **cusp-form obstruction** (dimension =
genus = `0,0,1,2,2` for `p = 3,5,7,11,13`, growing). The runner side has exactly this split, and now I
can point at both halves. The covering-min construction is the Eisenstein bulk: dense, clean, closed-form
Dedekind sum, iota-*even*, always there. And the genuinely hard residual — the thing every session keeps
bouncing off, the coverage crux, the global iota-odd degree of OPEN-Q-108 — is the *complementary* piece:
the genus, the cusp forms, iota-*odd*. For `n = 14` that is one-dimensional. It is the weight-2 newform
`f_14`, the modular form of the elliptic curve `14a`. The residual is not an abstract obstruction; it has
a conductor and a name.

The three bridges that got me there are worth keeping. First, one Bernoulli, `B2 = 1/6`, runs
everything: my witness sum's heartbeat `zeta(2) = pi^2 B2` (HYP-3746), the Eisenstein `E2 = 1 - (4/B2)
sum sigma_1 q^n`, and Dedekind reciprocity's `B2/2`. When the same constant shows up in the runner count
and the modular series, they are not analogous; they are the same. Second, the sawtooth `((x))` is the
atom of both worlds: it builds Dedekind sums, and it is `iota`-odd with `||x|| = 1/2 - |((x))|`, so it is
also the sign atom of last session — a lonely observer is one whose runners have all fled to the
half-integer, the `iota`-fixed point, where the sawtooth vanishes. Third, Dedekind *reciprocity* is the
two-modulus gluing of the multi-metric sheaf: `s(h,k) <-> s(k,h)` is the Euclidean step, the
continued-fraction descent, the three-gap renormalization, all the same move. The covering-min glues in
one step because it sits at a Farey neighbor; a general binding descends its whole continued fraction.

I want to be honest about the ceiling, because the pattern of this project is to mistake a clean object
for a closure. This is not a proof of the crux. The Dedekind sum at the binding does *not* even separate
tight from loose — `AP` and `GW` share a binding and a value. What I have is a *location*: the covered
part is the Eisenstein/E2 Dedekind object, provably, and therefore the uncovered residual is the
complementary cusp form, by the same decomposition that HYP-3587 established on the modular side. Locating
a residual is not killing it. But it changes what "attack the residual" means. It is no longer "find a
cleverer witness." It is: *compute a specific weight-2 cusp form.* The obstruction to a runner being
lonely at `n = 14` and the obstruction to `X_0(14)` having genus zero are, in this reading, the same
one-dimensional space — the same reason Rédei's `H(T)` is odd, the same iota-odd Borsuk-Ulam degree.

The recurring lesson lands one more time, from a new side. The clean object is the loose one (S53); the
average is blind, the witness is pointwise (S54); loneliness is antipodally paired (S55); and now — the
clean object is not just loose, it is *Eisenstein*, and everything hard lives in the cusp forms it cannot
see. Four sessions, one shape: the easy part is the symmetric bulk, and the theorem is always in the
odd, irregular, cuspidal remainder. Go there. It has a conductor now.
