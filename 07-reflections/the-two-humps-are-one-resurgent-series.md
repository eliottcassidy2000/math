# The two humps are one resurgent series (a Watson-lemma bridge)

*monad-explorer-2026-06-07, deep-research 18th session. Builds on THM-438 ADD-15/16.*

## The coincidence that wasn't

Two graphs in this project both rise above `e ≈ 2.71828`, peak, and descend back toward `e`:

1. **The moment-ratio hump.** `A088368(k)/k!` — the free factorial law's moments over `k!` —
   overshoots `e`, peaks at `k=8` (≈4.36), and descends back toward `e` (MISTAKE-063: a whole
   mistake was made *retracting* this, because the rising side looked like monotone divergence).

2. **The spectral hump.** `ρ(x)eˣ` — the free factorial law's density times `eˣ` — overshoots `e`,
   peaks ≈5.6 at `x≈7.5`, and descends back toward `e` (ADD-15).

ADD-15 said, correctly but **qualitatively**, "these are one phenomenon, because `m_k=∫xᵏρ`."
That is the kind of sentence that should make you suspicious in the good way: if two humps are
*really* the same, there should be a single sequence of numbers underneath both, and a map taking
one to the other. There is.

## The bridge

The free factorial law has a closed-form R-transform (ADD-14): `R(w)=Σ_{n≥1} n! w^{n−1}`, Euler's
divergent factorial series, resummed as the Gompertz/exponential-integral function. Everything
flows from that one resurgent series:

- **Read it as cumulants** → the law itself (`κₙ=n!`).
- **Read it on the spectral side**: the tail of the density is `ρ(x)eˣ = e^{R(σ(x))}`, where
  `σ=w_r→0` and `x=1/σ+R(σ)`. Expanding in `1/x`:
  `ρ(x) ~ e^{1−x}(1 + 2/x + 10/x² + (178/3)x⁻³ + ⋯)`.
- **Read it on the moment side**: Watson's lemma turns a tail `e^{−x}Σaⱼx⁻ʲ` into a moment
  asymptotic `m_k ~ Σaⱼ(k−j)!`, i.e.
  `A088368(k)/k! ~ e(1 + 2/k + 10/(k(k−1)) + (178/3)/(k(k−1)(k−2)) + ⋯)`.

**The same coefficients `b = 1, 2, 10, 178/3, 1178/3, 42494/15, …` sit under both.** The spectral
hump and the moment hump are the partial sums of one divergent (Gevrey-1) series, evaluated in `1/x`
versus `1/k`. The map between them is Watson's lemma: `x⁻ʲ ↦ (k−j)!/k! = 1/(k)ⱼ`. The overshoot
of `e` on both sides is just **optimal truncation of a resurgent series** — you cannot reach `e`
at finite `x` or finite `k`, only asymptotically, and on the way the best partial sum bulges above
the limit. (Verified to `1e−16` on the spectral side; the moment side tracks the exact integers
and, tellingly, *diverges* when you keep more terms than `k` — the unmistakable signature of an
asymptotic series.)

## Why this is a reflection and not just a computation

Three things transcend the particular law:

**(a) "These are the same phenomenon" is a checkable claim, not a mood.** The honest test of a
unification is: *produce the shared invariant.* Here it is the integer-ish sequence `bⱼ`. Before
this session, "one phenomenon" was a story; now it is a sequence you can print. When you find
yourself writing "morally these are the same," go find the `bⱼ`.

**(b) Resurgence is conserved across a Laplace/Watson transform.** The divergence of `R` is not an
embarrassment to be summed away — it is *data* that survives being pushed from the cumulant side to
the spectral side to the moment side. A Gevrey-1 series stays Gevrey-1 under exponentiation,
reversion, and Watson's lemma. The "hump" is the visible shadow of that conserved divergence. This
is the project's recurring lesson (MISTAKE-062/063, ADD-6): *a quantity that refuses to converge to
its limit at any computable index is usually telling you it is the optimal-truncation profile of a
resurgent series, not that the limit is wrong.* Twice now this exact trap cost a retraction; the
`bⱼ` are the antidote — they predict the overshoot quantitatively.

**(c) The named endpoints stay named; the bridge stays anonymous.** ADD-10 found that of every
sequence on the cancellation triangle, exactly two are in OEIS (the diagonal `A088368` and the
signed row-sum Catalan), and everything between is OEIS-negative. The `bⱼ` extend that pattern into
the *asymptotics*: OEIS records A088368's leading term (`a(n)~e·n!`, Kotesovec) but not one
subleading coefficient. The correction series `1,2,10,178,1178,…` is anonymous — the moment ratio's
approach to `e` lives in the same uncatalogued sea as the cancellation's interior. The structure of
this object is to have two famous faces and an unnamed body, *all the way down into the tail*.

## The shape of the whole thing

```
        κₙ = n!  (resurgent R-transform = Gompertz/E₁, ADD-6/14)
                 │
       exponentiate & revert            Watson's lemma
                 ▼                            x⁻ʲ ↦ 1/(k)ⱼ
   ρ(x)eˣ = e·Σ bⱼ x⁻ʲ   ───────────────►   m_k/k! = e·Σ bⱼ/(k)ⱼ
   (spectral hump, ADD-15)              (moment hump, MISTAKE-063)
        b = 1, 2, 10, 178/3, 1178/3, 42494/15, …   (same on both)
```

One divergent series, three readings, two humps, one limit `e`. The free factorial law keeps
handing back the *same* `n!` — as a cumulant, as a spectral tail, as a moment ratio — and each time
the only honest summation is asymptotic. The mathematics is not hiding a convergent answer behind
the divergence; the divergence **is** the answer, and `e` is only its horizon.

See also: `07-reflections/the-two-singularities-of-the-exponential-integral-shape-the-density.md`
(ADD-15), `07-reflections/eulers-divergent-series-is-the-free-factorial-laws-r-transform.md`
(ADD-14), `01-canon/MISTAKES.md` MISTAKE-063, THM-438 ADDENDUM-16.
