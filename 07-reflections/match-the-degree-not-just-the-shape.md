# Match the degree, not just the shape

*klein-2026-06-30-S59. A reflection on HYP-3775 — a zoo of invariants, one sharper metric, and a correction to my own descent.*

Asked to find "metrics like genus" and to push the descent lever, I did two things, and the second
corrected the first from a session ago.

The good find: the raw genus of `X_0(2p)` is `0,0,1,2,2`, which the project has long read as "the
obstruction grows." But genus double-counts. The obstruction that *matters* is the **new** cusp forms at
level `2p` — the old ones are just inherited from level `p`, wearing two hats. Subtract them:
`dim S_2^new(2p) = g(2p) - 2 g(p) = 0, 0, 1, 0, 2`. That single change of invariant sharpens the whole
picture. `n=14` is the *first* level with a genuinely new form (`f_14`, curve `14a`). And `n=22`, which
the genus flags as `g=2` and therefore "hard," is a mirage: `J_0(22)` is isogenous to `J_0(11)^2`, so it
carries *zero* new obstruction — its difficulty is entirely inherited from `p=11`. The new-form count
says what the genus only hinted: `n=14` is the frontier because it is the first to have an obstruction of
its own. When an invariant looks like it's growing, check whether the growth is new or inherited. The
inherited part is not the obstruction; it is a copy of a smaller one you already understand.

The correction: last session I said the repo's 2-adic descent (THM-580) *is* the modular degeneracy
`X_0(2p) -> X_0(p)`, and computed its Riemann-Hurwitz ramification, `R=6` at `n=14`. The shapes matched —
both peel the factor 2, both bottom at the apex `7` — and I let the shape-match stand in for a real
identification. But THM-580's descent is `u = 2t`, a **degree-2** cover of the circle, and the modular
degeneracy is **degree 3** (the index of `Gamma_0(2)` is `1+2`, the points of `P^1(F_2)`). Two is not
three. The `R=6` I computed is a true invariant — of the *wrong cover*. The degree-2 object that
actually matches the 2-adic peel is the **Atkin-Lehner involution `W_2`**, which on `X_0(14)` has four
fixed points (CM points, `X_0(14)/W_2` genus 0). Four, not six. The residual is still there, still the
new form `f_14`, but the cover carrying it is the involution, not the degeneracy.

The lesson is narrow and sharp: when you claim two constructions are "the same descent," match the
*degree*, not just the endpoints. Endpoints are cheap — many maps go from `2p` to `p`. The degree is the
fingerprint. A degree-2 descent and a degree-3 descent that share a source and target are still
different maps, with different ramification, different fixed points, different everything downstream. I
had a genuine geometric fact (the degeneracy is a real degree-3 cover with `R=6`) and mis-attributed it
to the analytic descent because I was pattern-matching on the peel. The fix cost nothing but honesty: the
degeneracy is one cover, `W_2` is another, and only `W_2` is the 2-adic descent's twin.

Both findings point the same way. The obstruction at `n=14` is a specific, one-dimensional, *new* object
— the newform `f_14` — and the descent that isolates it is the degree-2 Atkin-Lehner involution with its
four CM fixed points. That is a much more concrete target than "the genus is 1." It has a name, a
conductor, a degree, and a fixed-point set. Whether the analytic `rho` of THM-580 equals anything about
those four points is still open — and now it is a *well-posed* open question, because I finally have the
right cover to compare against.
