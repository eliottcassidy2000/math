# Parameter-complement guards and an infinite run insertion

**Status: PROVED** for the exact guard and common-future statements; **FINITE-EXACT** for the declared proposal search. ROOT completion of the resulting children is **OPEN** unless an independent child certificate is supplied.

The source chart is

\[
E(t)=924745897+2^{32}t,\qquad n(t)=2^{E(t)}-1,\qquad t\ge0.
\]

The inherited deep rule covered only \(t\equiv0\pmod{2^{47}}\). This package adds guards in its complement. The frozen finite bank has exact parameter density about \(0.19826862705329032\). A proved infinite run insertion adds exactly \(1/1024\), giving about \(0.19924518955329032\). These are densities of smaller-child receipts, not densities of completed ROOT certificates.

## 1. Inheritance and the retained coordinates

The closest mechanisms are the actual-word carrier and fixed-source compiler in [ordered context cancellation](collatz_context_cancellation_20261007d.md), the [eight-child phase](collatz_eight_child_routes_20261007d.md), and the endpoint coordinate in [marked periodicity](collatz_bott_marked_periodicity_20261007d.md). The old guard is retained in the new finite union. The cheapest remaining hostile is the entire cell \(t\equiv4\pmod{512}\), which neither the finite positive-gap bank nor its infinite run family covers.

For a valuation word \(w\), retain its ordered triple

\[
F_w(x)=\frac{P_wx+B_w}{Q_w},\qquad P_w=3^{|w|},\quad Q_w=2^{A_w}.
\]

The proposal reader also uses the exact integer

\[
I_w=3B_w-3P_w+Q_w,\qquad
I_{w a}=3I_w+2^{A_w+a}.
\]

Keeping the carry, actual source, valuation cost, and odd-clock difference is essential. A numerical collision at a proposed parameter is only a candidate until the full affine identity and native cylinder have been checked.

## 2. Every retained row proves an infinite guarded family

A row retains a source head \(u\), partner head \(v\), deletion \(D>0\), and sibling gap \(g>0\), satisfying

\[
|v|-|u|=D,\qquad A_u-A_v=2g,\qquad
F_v(-1)=4^gF_u(-1)+\frac{4^g-1}{3}.
\]

Let \(e=E-1\). After their respective initial valuation-one runs, the source and child have values

\[
x=2\cdot3^e-1,\qquad y=2\cdot3^{e-D}-1.
\]

The displayed identities imply \(F_v(y)=S^g(F_u(x))\), where \(S(z)=4z+1\). If the source head is native and its next valuation is \(c\), the partner's next valuation is exactly \(c+2g\). Thus the actual words

\[
1^e u c,\qquad 1^{e-D}v(c+2g)
\]

reach the same endpoint from \(n=2^E-1\) and the strictly smaller child \(2^{E-D}-1\). This is an actual common-future receipt. It does not assume the child has reached ROOT.

The source head's native cell is

\[
x\equiv (Q_u-B_u)P_u^{-1}\pmod{2Q_u}.
\]

The inherited exact power-of-three phase decoder transports this cell into a unique compatible dyadic cell in \(t\). Its generic sufficient growth cutoff is

\[
e\ge\max(2A_u,D+1).
\]

Indeed every source prefix has coefficient at least \((3/2)^e2^{-A_u}>1\), apart from the already growing initial run; its positive carry gives strict growth relative to the immutable source. The partner is positive and odd by the native endpoint identity and exact backward division. The first-hit receipt checks prevent a padded ROOT suffix. Every finite retained row has cutoff zero in the \(t\) coordinate.

`compile_rule` recomputes the full carrier identity, native phase, and cutoff. `receipt` consumes the supplied ordinary source; it never replaces it by a canonical sample. `symbolic_source` and `head_endpoint_residue` retain full denominator precision before modular division, so their arbitrary-modulus queries do not materialize \(2^{E(t)}\).

## 3. Declared finite discovery and exact coverage

The proposal universe is \(0\le t<1024\), \(1\le D\le32\), and source-head depth at most128, with1024 dyadic bits in the modular reader. At each still-uncovered sample, candidate rows are ordered by source cost, source depth, decreasing deletion, decreasing gap, then the literal words. Already covered samples are skipped. This is a deterministic bounded policy, not a complete classification of all possible source/child joins.

The saved JSON records all429 retained rows, the860 attempted parameters, and431 attempted failures. The resulting guard union covers593 of the1024 sample parameters. That sample fraction is not used as a density estimate.

The429 authenticated dyadic cells are pairwise disjoint. Their exact union mass is

\[
\beta=
\frac{87577585436686537861655933644679889301594718852088371532346146423787809}
{441711766194596082395824375185729628956870974218904739530401550323154944}.
\]

For scale, the11 cells of precision at most8 already have mass \(21/128\); the50 cells of precision at most12 have mass \(403/2048\). The largest retained precision is238. The cell \(t\equiv1\pmod{16}\) is a coarse new example, with head equal to the inherited15-letter partner followed by6 and deletion4. The old \(0\pmod{2^{47}}\) cell is included, so the finite increment is exactly \(\beta-2^{-47}\).

## 4. A reusable infinite operation

Write \(F_a(x)=(3x+1)/2^a\). The affine map

\[
S\circ F_2(x)=3x+2
\]

commutes with \(F_1\). Consequently, any valid gap-one pair with source head \((p,2)\) and partner \(v\) gives, for every \(k\ge0\), the pair

\[
(p,1^k,2),\qquad (v,1^k)
\]

with the same deletion and gap. This is a full ordered affine identity; neither a carry nor a source guard is discarded.

Apply it to the frozen row at \(t=29\). Put

\[
H=(2,2,4,1,3,3,1,3,3,1,1,1,2,2,3),
\]

\[
V=(2,2,1,2,2,1,1,2,1,1,2,1,2,2,2,3,2,4,2,1,1).
\]

Then

\[
u_k=(H,5,1^k,2),\qquad v_k=(V,1^k),\qquad D=4,\quad g=1.
\]

The costs are \(39+k\) and \(37+k\). The exact parameter cell has precision \(k+5\). For \(k=0,\ldots,8\), the residues are respectively

\[
29,5,117,21,213,853,1621,2133,5205.
\]

The cells are disjoint because the first valuation after the inserted ones is2 at different positions. Their total mass is \(\sum_{k\ge0}2^{-k-5}=1/16\). Every cell lies in the coarse parent \(t\equiv5\pmod8\).

For this family a uniform growth cutoff \(e\ge78\) suffices, independent of \(k\): the inserted ones only increase the relevant coefficient, and the fixed prefix plus final2 costs39. Hence the same lower bound \((3/2)^e2^{-39}>1\) applies to every source prefix. All chart exponents satisfy it. The general exported `run_rule(k)` deliberately retains its inherited conservative cutoff \(e\ge78+2k\); callers must honor that stated cutoff. The stronger uniform statement follows from the separate grammar proof, not from silently weakening the generic receipt API.

The finite bank meets this infinite family exactly in its six exits \(k=0,\ldots,5\). This is checked both by literal prefix classification and independently by dyadic-cell containment. Their mass is \(63/1024\), so the infinite extension adds exactly \(1/1024\).

## 5. Infinite-union density and finite overlap certificates

For every \(M\ge0\), all exits with \(k\ge M\) lie in the single native parent

\[
(H,5,1^M),
\]

whose parameter precision is \(M+3\). This shrinking parent bounds the entire tail, including any row-specific finite height cuts. Finite initial unions have their exact dyadic natural densities; the parent mass tends to zero. The sandwich proves existence of the full natural density and the stated geometric sum. It does not appeal to countable additivity of natural density.

For overlap with any finite dyadic bank of maximum precision \(B\), use \(M=\max(0,B-3)\). Test the individual exits \(k<M\). The remaining parent has precision at least \(B\), so it is either contained in one finite-bank cell or disjoint from the whole union. The remaining overlap is therefore either \(2^{-M-4}\) or zero. `run_overlap` implements this finite exact calculation; cell tuples are `(residue,bits)`.

Thus the finite bank plus the grammar has exact density \(\beta+1/1024\). The whole residual cell \(t\equiv4\pmod{512}\) remains outside: it misses the finite bank and also the grammar's \(5\pmod8\) parent. Negative-gap and ternary guards in concurrent packages are separate additions, not silently included here.

## 6. Reproduction and boundaries

The production default reads retained data, authenticates every row, independently checks modular actual prefixes and terminal reserves, replays429 modest ordinary-source receipts, and checks arbitrary-height modular endpoint identities. It also checks the run operation, tail overlap, residual cylinder, and malformed types/changed guards. It does no ROOT search.

From the repository root:

```powershell
python -B 04-computation/experiments/collatz_parameter_cover_20261007e.py
python -B -O 04-computation/experiments/collatz_parameter_cover_20261007e.py
python -B 04-computation/experiments/collatz_parameter_cover_20261007e.py --rediscover
```

`--rediscover` replays only the declared bounded proposal experiment and compares it with the frozen JSON. The explicit `--discover` option regenerates that data; normal production never writes it. No sampled frequency, successful proposal count, or paid child is promoted to universal or completed coverage. The next obligation is an additional guarded rule on the residual cells, followed by an independently authenticated ROOT proof at the selected smaller child.

Normal, optimized, and saved default output agree at **1,081,001 checks**; LF SHA256 `043c43dd4dfa97b6a915569515710d744fd9bc372ef379879ea040beadf4d18b`. Independent bounded rediscovery reproduces every retained row, attempted parameter, and miss, and passes **1,209,109 checks** including its discovery checks.
