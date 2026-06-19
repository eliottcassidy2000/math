# OCF Noise, Condorcet Spectra, and the LRC14 Coimage Lesson

Codex 2026-06-19 S15.

The tempting sentence was:

```text
OCF is noise stability at a special rho.
```

The session's answer is: close, but the clean object is not ordinary
two-copy noise stability.  The clean object is the hard-core partition function
at activity `2`.

For the independent-set indicator on the OCF conflict graph `Omega`,

```text
mu_p(independent) = (1-p)^m I(Omega, p/(1-p)).
```

So at `p=2/3`,

```text
H(T)=I(Omega,2)=3^m mu_{2/3}(independent).
```

That is the exact normalization.  It can be embedded into a degenerate
`rho=1` biased stability, and diagonal same-state mass can be tuned to activity
`2`, but nontrivial `rho<1` stability asks for ordered pairs of independent
sets.  The computation found the guardrail at `n=6`: two tournaments have the
same `H=23` but different pair/noise spectra.

The social-choice reading is good.  Since majority relations are tournaments
and every tournament is eventually majority-realizable, `{7,21}` are forbidden
Condorcet-cyclicity spectra.  The forbidden object is not a single Condorcet
cycle.  It is the compatible packet inventory

```text
alpha=(1, alpha_1, alpha_2, ...)
```

that would evaluate to `7` or `21` at activity `2`.

This is exactly the right analogy for the LRC(14) support-six tail.  HYP-2617
does for the LRC tail what OCF does for tournaments: it inserts the packet
address between raw mass and final signed evaluation.

```text
OCF:
  raw odd cycles -> compatible independent packets -> I(Omega,2)

LRC support-six:
  raw relation volume -> projective mod-7 coimage class -> signed reciprocal tail
```

The large-absolute/tiny-signed clue is therefore not an accident.  It means the
absolute count is pre-quotient boundary mass.  The proof-relevant object is the
compatible signed packet after lower faces and low-height walls are removed.

The next LRC computation should not try to improve the raw Minkowski absolute
envelope.  It should classify which of the `159` HYP-2617 coimage classes
remain possible after height-1 and height-2 wall deletion for `k=8,9,10`, then
bound only those non-null classes by the signed reciprocal-tail estimate.
