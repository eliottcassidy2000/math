# Rational run transitions: a finite lasso is the missing hypothesis

**CORRECTION / SCOPE AUDIT** of part3 of [THM-4603, run-compressed pair chains](../../01-canon/theorems/THM-4603-run-compressed-pair-chains-drift-lemma-ladder-completeness-and-long-run-transitions.md), as read at incoming commit `a0a33b38c`. The unconditional repetition argument is **UNPROVED**. The conditional lasso theorem below is **PROVED**. The declared55-cycle computation remains **FINITE-EXACT**, with the zero-drift labels corrected. This note does not audit the separate general ladder-completeness assertion in part2.

## 1. First failed implication

On rational numbers with odd denominator, the Terras map is

\[
T(x)=\begin{cases}(3x+1)/2,&x\text{ odd},\\x/2,&x\text{ even}.\end{cases}
\]

If \(x=a/d\) in lowest terms and \(d\) is odd, subsequent reduced denominators divide \(d\). This does **not** bound their numerators. Thus the sentence that the rational limit pair has bounded denominators and can therefore be iterated until it repeats does not establish termination.

The issue already contains arbitrary positive integer orbits. For any positive odd \(n\), take the child limit \(c'=-1\), which is fixed by \(T\), and the allowed affine state

\[
A(v)=v+n+1.
\]

Its limit pair is \((A(-1),-1)=(n,-1)\), with both denominators1. A universal repetition theorem for these pairs would prove eventual periodicity of every positive integer orbit, including the still-open exclusion of divergent positive trajectories. It is not a consequence of a denominator bound. This is a proof gap, not a constructed divergent orbit.

A separately established bound on the numerators would repair the finite-state argument. A more practical sufficient input is an explicit finite lasso, which can be checked without assuming search termination.

## 2. The exact conditional transition theorem

Let the rational limit pair \((u_t,v_t)\) have a supplied lasso: a transient of length \(\tau\) followed by a joint period \(p\). Define the debt by

\[
k_{t+1}=k_t+\operatorname{par}(u_t)-\operatorname{par}(v_t).
\]

Let \(o_u,o_v\) be the actual source and child odd counts during that joint period, and put \(\Delta=o_u-o_v\). If \(L=\tau+qp+r\), with \(0\le r<p\), then

\[
k_L=k_\tau+q\Delta+
\sum_{i=0}^{r-1}\bigl(\operatorname{par}(u_{\tau+i})-
\operatorname{par}(v_{\tau+i})\bigr).
\]

Consequently \(k_L=k_\tau+(\Delta/p)(L-\tau)+O(1)\), with an explicit bounded periodic remainder. This conclusion follows from the checked lasso, not from bounded denominators.

The finite-run transfer is also conditional and exact. If the actual initial pair has the stated affine relation and is 2-adically close to the limit pair to precision \(L\), then at each of the first \(L\) steps its branch agrees with the corresponding limit branch. Each step loses one bit of precision. The debt formula therefore governs that actual long run. The transient and period are fixed by the supplied rational lasso and do not depend on the chosen run length.

Use these four disjoint outcomes, with absorption given priority:

1. **Absorption:** a finite state has \(u_t=v_t\) and \(k_t=0\). The affine pair map is then the identity, so an actual pair following this prefix merges.
2. **Anchoring:** the periodic coordinates coincide, but absorption has not occurred. Their identical parity stream keeps debt constant, and the affine state fixes the current periodic point.
3. **Zero-drift periodic pair:** \(\Delta=0\), with distinct periodic coordinates. This includes both one cycle out of phase and **different cycles having equal odd density**. Debt and the affine pair state are periodic.
4. **Nonzero drift:** \(\Delta\ne0\). Debt has the exact nonzero mean \(\Delta/p\) above; the limit coordinates belong to different cycles.

The third outcome must not be restricted to one cycle out of phase. Nor does a finite lasso certificate imply a uniform transient bound over all rational inputs.

## 3. Two exact zero-drift hostiles

The rational cycle points of the actual words \((1,1,4)\) and \((1,2,3)\) are respectively \(19/37\) and \(23/37\). Their Terras cycles are

\[
\frac1{37}(19,47,89,152,76,38),\qquad
\frac1{37}(23,53,98,49,92,46).
\]

They are disjoint, each of period6 and odd count3. The state \(A(v)=v-4/37\) therefore gives a zero-drift pair on different cycles from time0.

The issue occurs in the incoming finite table itself. From the universal state \(A(v)=27v-26\), run the child at the cycle point \(7/55\) of word \((2,4)\). At time24 the limit pair is

\[
(u,v)=(2/55,7/55),\qquad k=8.
\]

The source cycle has period12,

\[
\frac1{55}(2,1,29,71,134,67,128,64,32,16,8,4),
\]

while the child cycle has period6,

\[
\frac1{55}(7,38,19,56,28,14).
\]

Both odd densities are \(1/3\), so their joint period12 has zero debt drift. The incoming script calls this `SHIFT`, although the cycles differ. Its entry time, entry debt, period, and drift are valid; the interpretation of that label needs widening.

## 4. Bounded reader and preserved finite computation

`read_lasso(source,child,debt,budget)` performs at most the supplied number of exact transitions. It returns `absorbed`, a checked `lasso`, or `pending` with the retained frontier and debt. There is no implicit unbounded iteration. `audit_lasso` checks every pair edge, the closing edge, the supplied initial pair, exact rational types, and the source/child odd counts. A forged cycle or a bool/float alias is rejected.

The independent finite universe is all primitive positive U-words of cost at most8 and length at most4, deduplicated by their actual rational Terras cycles, together with0:55 cycles. The five incoming states are

\[
(k,e)=(3,-26),(4,10),(3,13),(-5,1/243-1),(1,1).
\]

Each cycle is tested on both the source and child side, making550 cases. With a declared cap of2000 transitions, every case supplies an exact absorption or lasso witness; the largest required count is753. The corrected outcomes are:

| Outcome | Cases |
|---|---:|
| Absorption |96|
| Anchoring |24|
| Nonzero drift |150|
| Zero-drift periodic pair |280|

Of the280 zero-drift cases,224 lie on different cycles. These counts certify precisely the finite universe. They neither establish universal rational eventual periodicity nor invalidate the concrete finite witnesses. A separate control returns `pending` for \((27,-1)\) after10 steps, without drawing any conclusion about its later fate.

## 5. What the trace/carry result does transfer

The elementary carrier core of [THM-4604, reducible locus and extension cocycle](../../01-canon/theorems/THM-4604-collatz-words-are-the-reducible-locus-carries-are-the-extension-cocycle-and-traces-only-hear-the-tuning.md) is directly useful. With chronological composition,

\[
G_w=\begin{pmatrix}3^{|w|}&B_w\\0&2^{A_w}\end{pmatrix},
\qquad B_{uv}=3^{|v|}B_u+2^{A_u}B_v.
\]

Traces of these triangular word matrices see the diagonal characters and lose the ordered carry. A fixed point \(c_w=B_w/(2^{A_w}-3^{|w|})\) trivializes that carry on powers of the one word, not on every mixed word at once. The normalized traces satisfy the reducible Cayley-cubic identity; that scalar identity supplies no missing native guard.

The exact small hostile is

\[
(1,2):(P,Q,B)=(9,8,5),\qquad
(2,1):(P,Q,B)=(9,8,7).
\]

Their trace and multiplier agree, but their forward native cells are11 and9 modulo16. As **inverse endpoint** rules, they require respectively \(n\equiv4\pmod9\) and \(n\equiv2\pmod9\). Hence the first can enter a Mersenne endpoint \(n=2^E-1\) precisely when \(E\equiv5\pmod6\); the second never can, since it would require \(2^E\equiv3\pmod9\). This supports keeping the full carry and source guard in the new parameter compilers.

The valid blindness statement concerns trace/character coordinates. It should not be broadened to every frieze, minor, or cluster construction: minors involving the carry can retain information that traces discard. No general negative conclusion about such richer constructions is used here.

## 6. Reproduction

```powershell
python -B 04-computation/experiments/collatz_run_classification_scope_20261007e.py
python -B -O 04-computation/experiments/collatz_run_classification_scope_20261007e.py
```

The checker uses only exact integers and rational arithmetic. The incoming files are not edited by this package; the scope correction and its finite witnesses are separated for the owner to integrate.

Normal, optimized, and saved stdout agree at **882,261 checks**. LF SHA256: `e243c8235d3d6ace77b702996b5483ca843f51f3733f7a01124c9537a63ae48e`.
