# Bounded ten-body parity-free clock-three completion

Status: **PROVED ELEMENTARY + FINITE-EXACT + INDEPENDENTLY AUDITED
(THM-4442); LRC(14) OPEN.**

## Result and mechanism

For every ten-subset \(C\subseteq[13]\) and every three distinct positive
ternary units \(T=\{a<b<c\}\),
\[
                         G_{3C\cup T}\ne\varnothing.
\]
No parity or tail-primitivity assumption is used.

If a closed component of \(G_C\) has length \(L\), three-sheet counting plus
connectedness proves
\[
                       7cL\ge3\Longrightarrow G_{3C\cup T}\ne\varnothing.
\]
For each body point \(y\), a ternary-unit tail can kill at most one of the
three lifts \((y+j)/3\). If all three tails killed all sheets, each would kill
exactly one. The largest tail's three relatively open sheet-label regions
would disjointly cover the connected component, so one fixed lift interval
would lie in one open \(c\)-danger tooth. Its closed length \(L/3\) must then
be strictly below the tooth length \(1/(7c)\).

## Exact finite residue

Two exact body constructions agree for all \(286\) bodies. The smallest
longest component is
\[
 {9\over1232}
\]
at \(C=(1,2,3,5,6,7,8,9,11,13)\). Thus the remaining largest tail is at most
58. Body-specific cutoffs leave exactly 174,045 rows.

```text
positive rows                         174,045
isolated-only rows                          0
direct positive components          5,171,992
quotient mask pieces                3,643,548
safe quotient pieces                3,559,274
parity masks 0..7                    all present
low circuits 111/112/122            7,092 / 14,249 / 11,765
```

Quotient lift-mask measure equals direct literal danger-union measure on
every row. The minimum completed measure is
\[
 {10517879\over643242600}
\]
at body \((1,2,3,5,6,7,8,9,11,13)\), tail \((8,34,50)\), with twenty
positive components.

## Reproduction

```powershell
python -B 04-computation/lrc14_bounded_ten_body_parity_free_thm4442.py
python -B -O 04-computation/lrc14_bounded_ten_body_parity_free_thm4442.py
python -B 04-computation/lrc14_bounded_ten_body_parity_free_thm4442_independent.py
python -B -O 04-computation/lrc14_bounded_ten_body_parity_free_thm4442_independent.py
```

The primary engine uses exact endpoint geometry and labelled quotient masks.
The clean-room engine uses an integer \(14\operatorname{lcm}\) grid and
literal physical rows. See
[THM-4442](../../01-canon/theorems/THM-4442-lrc14-bounded-ten-body-parity-free-scale-three-completion.md)
for the proof and scope.
