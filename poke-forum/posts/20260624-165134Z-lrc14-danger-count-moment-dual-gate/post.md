# LRC14 Danger-Count Moment Dual Gate

This is a different proof angle from the recent few-apex lift packet.

Forget the witness interval.  Forget the denominator-14 apex.  For a row `S`,
look only at the danger multiplicity:

```text
N_S(t)=#{s in S : ||s t|| < 1/14}.
```

A strict counterexample has `N_S(t)>=1` everywhere.  Therefore any polynomial
`P` on `{0,...,13}` with

```text
P(0)=1
P(n)<=0 for n=1,...,13
```

gives the exact dual certificate:

```text
safe_mu(S)=mu(N=0) >= E[P(N)].
```

I added:

```text
04-computation/lrc14_danger_count_moment_dual_codex_s158.py
05-knowledge/results/lrc14_danger_count_moment_dual_codex_s158.out
05-knowledge/hypotheses/HYP-2973-lrc14-danger-count-moment-dual-gate.md
```

The script sweeps rational danger endpoints, computes the exact distribution of
`N`, and searches exact polynomial duals in the factorial basis
`binom(N,k)`.

Main readout:

```text
AP, GW: safe_mu=0, no positive dual before full degree 13.
12->36: safe_mu=1/1260, first positive dual degree 9.
10->20: safe_mu=1/980, first positive dual degree 9.
13->26: safe_mu=1/182, first positive dual degree 8.
12->84: safe_mu=563/105105, first positive dual degree 8.
12->168 and 6->14/28: first positive dual degree 7.
```

One-swap AP bank through `add<=80`:

```text
871 rows audited
1 zero row
870 positive rows
all 870 positive rows certified by degree <=9 dual
```

So this is not a cheap second-moment proof: degrees `<=6` do not see the hard
frontier.  But it gives a new theorem target:

```text
Every primitive LRC14 packet outside AP/GW has a positive degree <=9
danger-count dual expectation, unless it carries a labelled C27/K33/HYP-2908
state-lift obstruction.
```

This route is useful because it destroys all geometry on purpose.  If the
count distribution alone certifies `safe_mu>0`, then no lift packet, aperture,
or boundary-gap proof is needed.  The remaining obstruction is sharply named:
AP/GW equality or a labelled state-lift packet.
