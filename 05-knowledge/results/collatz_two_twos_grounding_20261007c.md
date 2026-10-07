# Grounding emitted two-twos Mersenne children

**Status: FINITE-EXACT** for the retained first-hit ROOT certificates and their
explicit transports. **PROVED** for the fixed-head exponent phase, paid
common-future family, and completed seed-lift corollary. **CONDITIONAL** for
ROOT completion of an arbitrary member of the paid family: its final smaller
child remains an obligation.

The inherited entry mechanism is the
[general head-to-Mersenne compiler](collatz_general_head_phases_20261007b.md),
with the audited clearing states of
[THM-4601, two-anchor reduction](../../01-canon/theorems/THM-4601-two-anchor-reduction-of-the-residual-first-reset-two-branch.md).
The new objects are retained ROOT evidence for selected emitted children and
an exact further route from one such child. An exponent congruence, a paid
dependency, and a completed ROOT certificate are kept as different types.

## 1. The emitted obligation and its corrected short type

Write \(M_E=2^E-1\), and
\[
U(n)=\frac{3n+1}{2^{v_2(3n+1)}}.
\]
The inherited depth-four exits emit \(E=K-4\equiv13\pmod {16}\).
Consequently \(v_2(E-1)=2\). After the \(E-1\) initial valuation-one steps
there are exactly two valuation-two steps. The source checkpoint before
those twos is \(x=2\cdot3^{E-1}-1\).

The short clearing state differs from the long-run state:
\[
X=\frac{9x+7}{16},\qquad
Y=\frac{x-17}{144},\qquad X=81Y+10. \tag{1}
\]
The depth-three child clears with \((4)\), and the depth-four child clears
with \((2,2)\). Reusing the long-run relation \(X=27Y-26\) here would be
incorrect. The retained sidecars are the actual run length, affine anchor,
word costs, and native source residue.

The natural first attempts illustrate why matching only a short type fails.
The changed-anchor heads \((10)\) and \((1,12)\) produce respectively
\(E=1289\bmod4096\) and \(E=27861\bmod32768\), which are \(9\) and \(5\)
modulo \(16\). They miss every emitted child. The head \((1,1,14)\), with
partner \((3,1,1,1,4,3,1)\), does produce
\(E=177949\bmod262144\), but its shifted parent phase \(K=E+4\) is
disjoint from every one of the 225 fixed \(J=3\) entry phases. The script
checks all 225 dyadic compatibility tests. These are guard-intersection
obstructions, not failures of the underlying affine identities.

## 2. Four finite obligations are now discharged

The data file
[retained certificates](collatz_two_twos_grounding_20261007c.cert.json)
stores the following strict first-hit valuation words. Each word has its
initial \(E-1\) ones represented once and its remaining positive byte
letters compressed with zlib/base85. Metadata records rank, valuation cost,
peak odd-state bit length, and SHA256 of the fully decoded byte word.

| Child exponent \(E\) | Parent exponent \(K\) | Odd rank | Child valuation cost | Parent valuation cost | Child peak bits |
|---:|---:|---:|---:|---:|---:|
|6125|6129|28762|51712|51716|9709|
|18269|18273|90248|161309|161313|28956|
|46941|46945|226498|405932|405936|74402|
|98477|98481|475425|852008|852012|156086|

Discovery was explicitly bounded: start after the proved initial one-run
and permit at most 600000 further odd steps. All four reached ROOT within
that cap. The default production run never discovers an orbit. It decodes
the supplied data, checks its source-owned metadata and hash, then replays
every valuation with constant state memory and rejects any step after ROOT.
An independent reverse replay reconstructs the source from ROOT using
\((2^a y-1)/3\), checking integrality, oddness and absence of an earlier ROOT.

The four parent packets are, in row order,
\[
\begin{array}{c|c|c}
u&v&r\\\hline
(1,10)&(1,1,2,2,3)&1\\
(10)&(3,1,1,3)&1\\
(9,2)&(3,1,1,3,1)&1\\
(1,2,9)&(1,1,1,1,1,1)&3 .
\end{array}
\]
Each authenticated child prefix is replaced by the corresponding source
prefix at their exact common future. The resulting four parent certificates
are replayed independently. Their ranks agree with the child ranks and
their valuation costs are larger by four. Thus these eight specific
Mersenne sources are completed, not merely reduced to named obligations.

The public functions root_word, verify_word, parent_receipt, discharge and
reverse_discharge consume explicit sources and supplied proofs. There is
no ROOT oracle inside a production receipt. Large orbit-state lists are
avoided; the stored word itself remains explicit evidence.

## 3. A compact further exit on an actual emitted-child phase

A bounded comparison of deletion depths \(1\) through \(32\), and source
heads of length at most \(2048\) after the initial one-run, found the
following compact collision at \(E=18269\):
\[
\begin{split}
u={}&(2,2,1,1,3,1,1,3,3,1,3,1,2,1,2,1,4,7,1,1,2,1,2,4),\\
v={}&(5,1,4,3,2,3,1,1,1,2,1,1,1,1,2,1,4,1,2,2,2,1,2,1,3).
\end{split} \tag{2}
\]
This finite discovery is only provenance. The following algebra proves the
whole guarded family.

For \(F_w(z)=(3^{|w|}z+B_w)/2^{A_w}\), the exact carriers are
\[
\begin{array}{c|r|r|r}
 &3^{|w|}&2^{A_w}&B_w\\\hline
u&282429536481&1125899906842624&494713926107839\\
v&847288609443&281474976710656&776753761891457 .
\end{array}
\]
Thus \(|u|=24,\ |v|=25,\ A_u=50,\ A_v=48\), and
\[
F_v(-1)=4F_u(-1)+1. \tag{3}
\]
The source checkpoint is \(x=2\cdot3^{E-1}-1\); the child checkpoint is
\(y=2\cdot3^{E-2}-1\), so \(x=3y+2\). Equality of slopes together with
(3) gives the full identity
\[
F_v(y)=4F_u(x)+1. \tag{4}
\]

Let \(x_u=(2^{50}-B_u)3^{-24}\bmod2^{51}\) be the native odd-source
residue. Its exact pullback is
\[
3^{E-1}\equiv(x_u+1)/2\pmod {2^{50}},
\qquad
\boxed{E\equiv18269\pmod {2^{48}}}. \tag{5}
\]
The order of \(3\) modulo \(2^{50}\) is \(2^{48}\), so this is the
minimal exponent period and an iff for the native head, rather than a
sampled condition. Modular controls use parameters \(0,1,9,10^{30}\)
without materializing the enormous Mersenne integers.

Put \(Z=F_u(x)\) and \(c=v_2(3Z+1)\). The exact common-future words are
\[
\begin{array}{c|c}
M_E&1^{E-1}\,u\,c\\
M_{E-1}&1^{E-2}\,v\,(c+2).
\end{array} \tag{6}
\]
For every positive phase member, \(E\ge18269\).
During the source head its multiplicative coefficient relative to the
original source is at least
\((3/2)^{E-1}2^{-50}>1\); additive carries are positive.
Hence no source prefix through the head reaches ROOT, and \(Z>M_E\).
Equation (4) gives a positive odd child-head endpoint \(4Z+1>1\).
Reverse integrality and oddness therefore force every child valuation to
be exact; a prior ROOT would remain at ROOT under legal positive words and
cannot lead to \(4Z+1>1\). The terminal sibling identity then proves (6),
including first-hit safety.

The dependency is strictly paid because \(M_{E-1}<M_E\).
This does not assert ROOT for every phase member. At the retained point
\(E=18269\), reverse substitution of the already authenticated source
ROOT suffix instead proves ROOT for \(M_{18268}\), with rank \(90248\)
and valuation cost \(161308\). No further orbit search is used.

## 4. An infinite two-stage paid family and a completed subfamily

The phase in (5) intersects the old \((10)\), \(J=3\) parent packet:
\[
K=18273+2^{48}t,\quad
M_K\longrightarrow M_{K-4}\longrightarrow M_{K-5}. \tag{7}
\]
Each arrow means an authenticated smaller-child common-future receipt.
The narrower ray
\[
\boxed{K=18273+729\cdot2^{48}t,\qquad t\ge0} \tag{8}
\]
stays outside the declared \(G_1/G_5/G_{17}\) entry bank. Nonempty Mersenne
entry occurs when either \(K\equiv5\bmod6\) or
\(K\equiv733\bmod1458\); \(G_1\) never enters. On (8), the actual fixed
residues are \(K\equiv3\bmod6\) and \(K\equiv777\bmod1458\), since
\(1458\mid729\cdot2^{48}\). Both entry predicates are therefore false
for all parameters; 64 controls exercise the implementation. This is a
comparison with that entry bank, not with every prior adaptive controller.
ROOT for arbitrary members of (7) or (8) still requires the final child.

There is also an unconditional infinite completed family from each newly
authenticated finite seed, by the inherited
[completed terminal-lift mechanism](collatz_terminal_lifts_20261007.md).
For a seed \(a>1\) with retained ROOT word \(w\), let
\(L=|w|,\ P=3^L,\ Q=2^{A_w}\). For every integer \(s\ge0\), set
\[
N_s=a+\frac{4Q(4^{Ps}-1)}{3P},\qquad J=1+Ps. \tag{9}
\]
The quotient is integral since \(3P\mid4^{Ps}-1\).
The supplied word takes \(N_s\) to \((4^J-1)/3\). At \(s=0\) keep
the seed word alone; at \(s>0\) append the single valuation \(2J\).
Every proper prefix is positive and increases with the source parameter,
so the seed's strict first-hit property excludes ROOT padding. This gives
infinitely many completed positive sources from each retained seed, though
they are not asserted to be Mersenne numbers. Formula (9) is symbolic
proof storage; it does not justify expanding its enormous parameters.
The construction is credited to the earlier finite-seed and terminal-lift
grammar, not claimed as a new inverse-basin theorem.

## 5. Reproduction and boundaries

Run from the repository root:

    python -B 04-computation/experiments/collatz_two_twos_grounding_20261007c.py
    python -B -O 04-computation/experiments/collatz_two_twos_grounding_20261007c.py

The optional flag --rediscover explicitly reruns the four bounded orbit
discoveries and compares them with the frozen words. It is not used by
ROOT consumers. Default controls include independent reverse replay,
actual/affine receipt agreement, exponent-phase guards, all 225 target
compatibility tests, malformed exact types, a corrupted ROOT suffix,
wrong native phase and explicit materialization limits.

The finite completed points, all-height paid phase, and infinite constructed
completed seed-lifts have separate quantifiers. No finite certificate bank
or modular phase test is used to claim universal Collatz coverage.

Normal and optimized execution agree on **11,862,707 checks** and the
saved output. The LF-normalized stdout SHA256 is
`41fa847ba19dcf449eef219c9b4adc063a89dfb48acd245f7d97f2ed2973f084`.
The measured authentication times were approximately 61 and 63 seconds;
these timings are not included in the reproducible output.
