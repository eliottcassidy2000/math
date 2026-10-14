# Signed Collatz pairs: complete inverse fibres and authenticated terminal lifts

**Status: PROVED elementary identities and certificate transport;
FINITE-EXACT supplied terminal words and bounded controls.** These results
do not classify every integer cycle or prove that every positive source
reaches a known cycle. No priority claim is made for inverse-fibre identities.

The [checker](../../04-computation/experiments/signed_collatz_pair_lifts_20261014d.py)
and [exact output](signed_collatz_pair_lifts_20261014d.json) retain the source,
sign parameter, valuation, endpoint and terminal cycle separately.

## 1. Inheritance and the two clocks

The closest proved mechanism is the signed-parameter conjugacy and inverse
fibre in [the signed-cycle note](arithmetic_braids2_20260917_signed_cycles.md),
sections 1, 3 and 4. Its hostile is confusing parameter reflection with
time reversal. [THM-4472, four-vertex reading](../../01-canon/theorems/THM-4472-four-vertex-reading-of-3n-plus-minus-1.md)
keeps those operations distinct. The least-used coordinate here is the
**sheet of the continuation**, even after the integer endpoints coincide.
The live objects are same-source pairs, inverse fibres, affine carries,
signed terminal cycles and source-bound first-entry certificates.

For $b\in\{-1,1\}$ define the shortcut and fully odd maps by

$$
T_b(n)=\begin{cases}n/2&n\text{ even},\\(3n+b)/2&n\text{ odd},\end{cases}
\qquad U_b(n)=\frac{3n+b}{2^{a_b(n)}},\quad
a_b(n)=v_2(3n+b).
$$

The latter uses nonzero odd signed integers and the valuation of the
absolute numerator. It is always defined on this domain. One $U_b$ edge
contains exactly $a_b(n)$ shortcut steps. In particular, $1$ is fixed by
$U_+$ and $U_-$, but its shortcut cycles are respectively $(1,2)$ and $(1)$.
Zero is a separate shortcut fixed point, outside the fully odd domain.

For any odd signed $n$, the two shortcut odd outputs satisfy

$$
F_+(n)-F_-(n)=1,
\qquad F_+(n)\equiv2\pmod3,
\qquad F_-(n)\equiv1\pmod3.
$$

Exactly one of $a_+(n),a_-(n)$ is $1$, and the other is at least $2$.
This follows because $3n+1$ and $3n-1$ differ by $2$. Thus

$$
2^{a_+(n)}U_+(n)-2^{a_-(n)}U_-(n)=2,                 \tag{1}
$$

while the odd endpoints have residues
$U_+(n)\equiv(-1)^{a_+(n)}$ and
$U_-(n)\equiv-(-1)^{a_-(n)}\pmod3$.
The fixed residues $2/1$ belong to the shortcut outputs, not generally to
the fully odd endpoints.

Negation is an exact forward conjugacy:
$T_{-b}(-n)=-T_b(n)$ and $U_{-b}(-n)=-U_b(n)$, with the same valuations.
It transports a terminal to its negative on the other sheet, not to the
same positive terminal.

## 2. Complete inverse fibre over both sheets

**Theorem.** Fix a nonzero odd $q$ not divisible by $3$. For each $a\ge1$,
let $b_a\in\{-1,1\}$ be the representative of $2^a q\pmod3$. Then

$$
n_a=\frac{2^a q-b_a}{3},\qquad U_{b_a}(n_a)=q,
\qquad a_{b_a}(n_a)=a.                              \tag{2}
$$

These are exactly all integer source edges over both sheets with endpoint
$q$, with no repetitions of the marked pair $(b,n)$.
If $3\mid q$, there is no such edge. Every source in (2) has the sign of
$q$. In particular, the whole fibre is positive when $q>0$.

Indeed the numerator in (2) is an odd multiple of $3$; the equation
$3n_a+b_a=2^a q$ proves exact valuation. Conversely every edge has that
equation, which uniquely determines $b_a$ and $n_a$. Modulo $3$, its
nonzero right side is necessary. The valuation is itself unique, proving
the nonrepetition claim.

The enumeration has the exact recursion

$$
b_{a+1}=-b_a,\qquad n_{a+1}=2n_a+b_a.               \tag{3}
$$

Starting from any authenticated port $(b,n,a,q)$, $m\ge0$ lifts give

$$
L_{b,m}(n)=2^m n+b\frac{2^m-(-1)^m}{3},\quad
b'=(-1)^m b,\quad
3L_{b,m}(n)+b'=2^m(3n+b).                            \tag{4}
$$

Consequently the endpoint is unchanged and the valuation is $a+m$.
For $m=1$ this is $U_{-b}(2n+b)=U_b(n)$; for $m=2$ it is the inherited
same-sheet sibling identity $U_b(4n+b)=U_b(n)$.
No orbit conjugacy is asserted by (3): the common endpoint still has two
different possible next maps.

## 3. A valid seam and an infinite false-seam family

An exported terminal certificate is $(b,n,w,C)$: $w$ is the exact positive
valuation word from source $n$ to the first visit to a supplied simple
cycle $C$ under $U_b$. The cycle is replayed as well. Every preterminal
state must lie outside $C$; empty words are permitted exactly at $C$.

Given a port $(b,n,a,q)$ and a certificate $(b,q,w,C)$, prepend $a$ to $w$.
If $n\in C$, instead return its empty first-entry certificate. These are
the only two possibilities, since the first edge reaches $q$.
The API rejects both endpoint mismatches and sheet mismatches. All tuple
members and integers are reauthenticated before equality comparisons.

For example, $U_+(3)=U_-(7)=5$, but $U_+(5)=1$ and $U_-(5)=7$.
There is an infinite version: for every $j\ge0$,

$$
n_j=\frac{10\cdot4^j-1}{3},\qquad
m_j=2n_j+1=\frac{20\cdot4^j+1}{3}.                  \tag{5}
$$

The plus source $n_j$ has first-ROOT word $(2j+1,4)$.
The minus source $m_j$ enters the non-ROOT cycle $(5,7)$: its first-entry
word is empty for $j=0$, and $(2j+2)$ for $j\ge1$.
Both first edges have integer endpoint $5$. Thus even an exact common
integer endpoint does not authorize copying the other sheet's ROOT proof.

There is also a sharp word obstruction. A nonempty positive valuation word
$w=(a_1,\ldots,a_\ell)$ has

$$
F_{b,w}(n)=\frac{P n+bB}{Q},\quad
P=3^\ell,\ Q=2^{\sum a_i},\quad
B\equiv2^{a_1+\cdots+a_{\ell-1}}\not\equiv0\pmod3.
$$

If the same word on plus source $n$ and minus source $m$ had the same
integer endpoint, then $P(m-n)=2B$, impossible modulo $3$.
Changing sheets while retaining an integer common endpoint necessarily
changes the valuation word; (4) records precisely the one-edge cost change.
This is not an obstruction to negation, which negates the endpoint too.

## 4. Actual terminal records and a grounded infinite family

The four displayed words in the JSON were found using a declared cap of
10,000 odd steps for these four marked inputs only. Production functions
replay supplied words; they never perform orbit discovery.

| Source | Sheet | First-entry odd steps | Shortcut steps | Terminal |
|---:|---:|---:|---:|---|
| 1473 | $+1$ | 33 | 63 | $1$ |
| 1473 | $-1$ | 11 | 22 | enters at $61$ in the cycle through $17$ |
| 26861 | $+1$ | 48 | 91 | $1$ |
| 26861 | $-1$ | 45 | 86 | $1$ |

The numerical tokens originally have a separate prime-race role: 1473 is
a residue-class prime count at cutoff 26861, recovered in the parent
session. Treating either token as a Collatz source here is an explicitly
new finite control, not an identification of those two mathematical roles.

The two independently authenticated words at the unchanged endpoint
$q=26861$ now ground an **entire combined inverse fibre**. For every $a\ge1$,

$$
b_a=(-1)^{a+1},\qquad
n_a=\frac{26861\cdot2^a-(-1)^{a+1}}{3}              \tag{6}
$$

reaches $1$ under its single fixed map $U_{b_a}$, with no sign switching
along that orbit. The exact odd rank is $49$ for odd $a$ and $46$ for even
$a$; the shortcut rank is respectively $a+91$ and $a+86$.
The first two sources are $17907$ on the plus sheet and $35815$ on the
minus sheet. All these sources exceed $1$, so prepending the edge introduces
no ROOT padding. This is a constructive infinite grounded family, not a
claim that the inverse fibre contains arbitrary sources.

**All-exponent Mersenne boundary (PROVED).** This particular grounded
fibre contains **no** positive Mersenne source $M_E=2^E-1$, $E\ge1$.
More generally,

$$
U_+(M_E)=3\cdot2^{E-1}-1,\quad a_+(M_E)=1\quad(E\ge2),
$$
$$
U_-(M_E)=3\cdot2^{E-2}-1,\quad a_-(M_E)=2\quad(E\ge3). \tag{7}
$$

These follow by factoring $3\cdot2^E-2$ and $3\cdot2^E-4$ and checking
that the displayed remaining factor is odd. The small exceptions are
$U_+(M_1)=1$ with valuation $2$, and
$U_-(M_1)=U_-(M_2)=1$ with valuations $1,3$.
Thus for $q>1$ the complete fibre has a Mersenne source iff
$(q+1)/3=2^r$ for an integer $r\ge1$. In that case its only such ports
are $M_{r+1}$ on plus with valuation $1$, and $M_{r+2}$ on minus with
valuation $2$. Here $(26861+1)/3=8954=2\cdot4477$ is not a power of two,
and $q\ne1$, proving the claimed empty intersection for every $E$.
The `mersenne_ports(q)` API implements this exact classification.
Consequently (6) does not ground any of the outstanding positive Mersenne
sources; it is a different, explicitly proved infinite domain.

## 5. Known integer cycles and the denominator boundary

[THM-4641, Wilson dictionary, row 2b](../../01-canon/theorems/THM-4641-wilson-two-four-eight-fermat-multiplier-dictionary-structural-rows-and-coincidences.md)
records unit-denominator cycles and an explicitly bounded census. Our
controls verify the four known nonzero $U_+$ integer cycles:

$$
(1),\quad(-1),\quad(-5,-7),\quad
(-17,-25,-37,-55,-41,-61,-91).
$$

Their odd lengths/halving costs are $(1,2),(1,1),(2,3),(7,11)$.
The cycle equation is $n=B/(2^A-3^\ell)$, with positive ordered carry $B$.
The first three displayed denominators have absolute value $1$; the last
has denominator $-139$ and carry $2363$, giving $-17$.
The JSON also expands every cycle into its exact shortcut cycle.
Negation gives the corresponding cycles on the minus sheet.

This calculation does not prove that no other integer cycles exist.
Nor does a known-cycle list prove that every source enters a listed cycle.
The complete theorem in section 2 concerns a single specified inverse
fibre, not the union of every forward basin.

## 6. Reproduction and remaining obligation

Run both:

```text
python -B 04-computation/experiments/signed_collatz_pair_lifts_20261014d.py
python -B -O 04-computation/experiments/signed_collatz_pair_lifts_20261014d.py
```

Controls cover all odd signed sources from $-999$ to $999$, both sheets,
lift counts $0$ through $8$; odd endpoint units of magnitude at most $99$
and valuations $1$ through $16$; every word of lengths $1$ through $5$
over valuations $1,2,3,4$; four supplied terminal words via an independent
division-loop replay; the four known cycles; 32 instances of (5), 64 of
(6), the complete Mersenne classifier checked against all exponents
$1\le E\le128$ on both sheets and an independent scan for odd endpoints
$1\le q\le999$, and malformed-sheet/source/word/terminal controls.

Normal and optimized runs agree on **33,036 checks** and the saved LF JSON,
SHA256 `5c2e834e95a1a0f8aeff2764ee5189e1c77db8bee71d81f0c93dc1ddb7dfe047`.
Independent peer proof/API audits of both the core and the Mersenne
classifier delta pass; the independent optimized replay matches this
saved output. No correction remains.

The next source-specific obligation is explicit: a port endpoint needs an
authenticated terminal proof **on its own sheet**. Neither consecutive
shortcut outputs, reflection, a cycle denominator, nor a shared endpoint
discharges that obligation. Equation (6) shows how two supplied proofs
can discharge infinitely many ports while leaving general basin coverage
open.
