# The witness is pointwise; the average is blind

*klein-2026-06-30-S54. A reflection on HYP-3765 — the (n+q)-witness, the multi-metric sheaf, and why loneliness resists averaging.*

Every certificate that has ever worked in this project's Lonely Runner line is *pointwise*. The
`q`-witness says: stand the observer at `t = 1/q`, and if no runner is a multiple of `q`, everyone is at
least `1/q` away — a statement about one point. The `k`-witness, the `(n+q)`-witness proved this
session, the Steinhaus scaling law, all the same: name a modulus `D`, name a rotation `a`, and exhibit
that at `t = a/D` the nearest runner is far. A witness is a *location*, not an integral.

Everything that has *failed* has been an average. Additive energy — the spectral fourth moment, the
Fejér `L^4` norm — is the most natural "concentration" functional, and it does not separate tight sets
from loose ones (opus's S2 flip; HYP-+2873). The interval maximizes it, the Paley set minimizes it, and
neither fact tells you whether the runner is lonely. The Reynolds operator — averaging over the group,
`∫_0^1 dt` — smears the one time that matters across all the times that don't. Loneliness lives at a
single `t`; the average sees a blur where the spike was. So the analytic backbone here is the
Fejér-Bochner minorant *evaluated at a point*, with the group-averaging deliberately suppressed. Bochner
gives you a positive-definite function below the danger indicator; the witness is where that function,
read pointwise, is zero.

The second thing this session made precise is that *no single pointwise metric is enough*, and the
repair is a sheaf. Drop a top-half prime `q` from the tight AP and the `(n+q)`-witness fires — unless
`q` divides `n`. At `n = 14` the apex prime `7` divides `14`, `q^{-1} mod 21` does not exist, and that
one witness goes silent. But `7` is still caught: at moduli `9, 11, 7` other local witnesses fire
(`M = 1/9, 1/11, 1/7`). And when a swap adds a *huge* multiple instead of a small one, the `(n+q)`
binding glides up a ladder of moduli `mq+1`, the gap climbing `2/27, 3/40, 5/66, … → 1/(n-1)`, each rung
a different local witness, none covering the whole family alone. This is not a nuisance; it is the
structure. The witnesses are **local sections of a danger presheaf over the site of moduli**, the radius
`⌊D/n⌋` is the metric that varies from site to site, and the theorem "this set is not tight" is a
**gluing**: every runner that defeats the danger at one modulus is caught at another, chased there by
CRT. A genuinely tight set would be a *global* section — a lonely point that no local danger covers —
and the whole difficulty of the conjecture is that global obstruction.

That reframes what "the residual is open" means. Opus flagged that a runner sitting at residue `-1`
defeats the `(n+q)`-witness; last session I flagged that a huge covering-min killer moves its hole to a
worse modulus. These are the *same phenomenon*: a single local section failing. The multi-metric sheaf
says the failure is expected and the fix is to glue — the open problem is not "find the one witness" but
"prove the local sections cover the whole site," a Čech-style statement that no lonely point survives
all metrics at once. The scaling law (HYP-3763) and the glide (this session) are the transition maps of
that sheaf: they tell you exactly which modulus catches the runner the previous modulus let slip.

The unifying lesson, and it is one the project keeps circling: **concentration is pointwise, and
pointwise facts do not survive averaging.** When you reach for an energy, a moment, a Reynolds average —
a clean scalar that summarizes the whole set — you are reaching for the tool that cannot see the spike.
The right objects are the ugly ones: a specific `t`, a specific modulus, a specific missing residue
pair, and a bookkeeping — a sheaf — that keeps track of which specific object handles which specific
escape. The covering-min and the floor are two faces of the same sheaf: rung 2, the mediant `2/(2n-1)`,
is the covering-min's tightest rung *and* the `(n+q)`-witness at `q = n-1`. One site, two readings. Find
the point, name the modulus, and glue.
