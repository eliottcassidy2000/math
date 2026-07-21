# HYP-8771 -- finite zero recurrence for primitive toral trinomials

**Status: RESERVED / OPEN.** For coprime positive `p0,q0`, put `r=p0+q0` and

```text
A_m(kappa)=sum_{0<=k<=m/r}
  m!/((q0*k)!(p0*k)!(m-r*k)!) * kappa^k.
```

Conjecture: for every fixed `kappa`, the set `{m>=1:A_m(kappa)=0}` is finite.
This is the exact toral factor governing THM-2018/THM-2021's proportional
central resonance. The symmetric `p0=q0=1` case is Legendre finite zero recurrence,
announced in 2026 by Mangoubi--Kadets--Weller Weiser; the higher-charge case is
open here. Full motivation, evidence, and Tournament Analysis will replace
this reservation stub in the next checkpoint.
