# Chasing the predictions: the Pell supersymmetry is a real GMC(2) moment identity, and α(E_n)=2^{n−4}

**death-star-2026-07-21-S81** (HYP-8696). Owner: chase the S80 predictions; think the A+B+C−D−E−F=G /
A+B−C+D−E−F+G / A+B−C recursion modes as literally subtournaments — which iso classes come from smaller
subtournament classes. Two predictions **confirmed exactly**; the recursion modes resolve to the order-join.

## Prediction 1 — CONFIRMED: the Pell/Chebyshev identity is a GMC(2) moment identity
The a/b Pell identity $E_n^2-O_n^2=(x^2-1)^n$ (kps THM-1880) is the **polynomial lift of an exact GMC(2)
moment identity**. With $a\leftrightarrow Z$, $\bar a\leftrightarrow W$, $x^2-1=a\bar a\leftrightarrow ZW=s$
(radial), and $E$ the Gaussian functional ($E[Z^aW^b]=a!\delta_{ab}$), set the bosonic/fermionic pair
$\mathrm{sym}_n=(Z^n+W^n)/2$, $\mathrm{alt}_n=(Z^n-W^n)/2$. Then (exact Wick, verified $n=1..7$):
$$\boxed{\ E[\mathrm{sym}_n^2]-E[\mathrm{alt}_n^2]\;=\;E[(ZW)^n]\;=\;n!\ }$$
and generally, for **any** $P$ with charge-conjugate $\tilde P$ ($Z\!\leftrightarrow\!W$),
$$E[\mathrm{sym}(P)^2-\mathrm{alt}(P)^2]\;=\;E[P\cdot\tilde P]$$
(a difference of squares, $\mathrm{sym}^2-\mathrm{alt}^2=P\tilde P$; verified on several $P$). So the
"supersymmetry" prediction is real: **bosonic$^2-$fermionic$^2=$ the charge-conjugate norm**, and on the pure
mode $P=Z^n$ it localizes on the radial power $(ZW)^n$ with $E$-value $n!$. This is exactly the Pell relation
$E_n^2-O_n^2=(x^2-1)^n$ pushed through $E$ (which sends $s^n\mapsto n!$). The a/b Chebyshev pair *is* the
polynomial shadow of the GMC radial moment, made precise. `pell_supersymmetry_..._S81.py`.

**Why it matters.** $P\tilde P$ is manifestly charge-symmetric (charge-0-supported under $E$), so the identity
recasts the even/odd split as a *difference of squares whose value is a genuine (radial) moment* — a positivity-
flavored handle: $E[\mathrm{sym}^2]-E[\mathrm{alt}^2]$ is not a cancellation but the honest radial mass $n!>0$.
It sits on the toral/parity (DvdK-proved) side (S80), so it is a *provable* companion to GMC(2)'s open radial gap,
not the gap itself.

## Prediction 2 — CONFIRMED through n=7: α(E_n)=2^{n−4}
The even-graph metagraph independence number (S80: $\alpha(E_n)=1,1,2,4$ for $n=3..6$) continues:
**$\alpha(E_7)=8$** (built $E_7$, $V=54=$ A002854 exactly via a degree+spectrum+triangle hash-canonical; max
independent set by Bron–Kerbosch on the sparse complement). So
$$\alpha(E_n)=2^{\,n-4}\quad(n\ge4):\qquad 1,2,4,8\ \text{at}\ n=4,5,6,7,$$
a clean new law for the dual metagraph — while its $G_n$-dual $\alpha(G_n)=2,5,18$ has no such closed form (an
asymmetry between the tournament metagraph and its even-graph dual worth explaining). Prediction for $n=8$:
$\alpha(E_8)=16$.

## The recursion modes — resolved to the order-join
Read "A+B+C−D−E−F=G", "A+B−C+D−E−F+G", "A+B−C" as **subtournaments combined into a larger class**. Two readings:
- **Signed-additive deck** (my first interpretation: signed sums of the vertex-deleted subtournament invariants
  with the sign patterns ++−, +++−−−, ++−+−−+). Tested on $c_3$ and the signed Rédei $R$ at $n=3,6$: **no clean
  universal law** ($c_3$ 1280/32768, $R$ 128/32768 at $n=6$). The additive-signed reading is inconclusive.
- **Multiplicative order-join** (the correct "smaller→larger" mechanism, and where the fleet already is): every
  tournament decomposes uniquely into strongly-connected **atoms** joined in a linear order (the condensation),
  and the invariants **compose over the join** — $H$ and $|R|$ **multiply** (mac-mini THM-1936: $R$ join-
  multiplicative), $\mathrm{disc}$ **super-multiplies** (klein THM-1950's SCC law $\mathrm{disc}=\prod\mathrm{disc}(C_i)\cdot
  [\prod(1+s_i)+\prod(1-s_i)]/2^r$, with $s$ obeying the SL$_2$ velocity-addition $s(C_1{\Rightarrow}C_2)=(s_1{+}s_2)/(1{+}s_1s_2)$).
  So **the iso classes that "come from smaller subtournament classes" are exactly the reducible ones** (order-
  joins of strong atoms); the strong tournaments are the atoms that do *not*. The ± sign patterns are, in this
  light, the *permutation signs inside $R=\sum\mathrm{sgn}(\pi)$* — which is why $R$ (signed) still multiplies
  cleanly over the join: the cross-arc block contributes zero inversions (mac-mini's 2-line proof). The honest
  answer: the recursion is **multiplicative over the condensation, not additive-signed over the deck**.

## Fleet convergence (credited)
- **klein THM-1950** reduced my $H\ge\mathrm{disc}$ (HYP-8636) to the strong base via exactly this composition
  algebra — my S78 conjecture is now one strong-base inequality from proof.
- **mac-mini THM-1936** ($R$ join-multiplicative) is the signed "smaller→larger" law the recursion modes want.
- **kps THM-1885** (a/b $=$ BS(1,2)) + **boxeph THM-1926** (zeta) frame the a/b side the Pell identity lives on.

## Honest scope
Two exact confirmations (Pell moment identity, verified Wick $n\le7$ + general; $\alpha(E_7)=8$, exact) and a
resolution of the recursion modes to the (fleet-owned) order-join composition. The recursion modes' *additive-
signed* reading is a recorded dead end. No open problem closed; a confirmed prediction + a new metagraph law +
a clarified mechanism. Cross-links: S80 (the a/b shadow + $E_n$ sweep), kps THM-1880/1885, klein THM-1950,
mac-mini THM-1936, boxeph THM-1926, THM-1620 (Legendre/Hermite). Scripts
`pell_supersymmetry_and_deck_recursions_deathstar_S81.py` (+out), inline $E_7$ $\alpha$-check. HYP-8696.
