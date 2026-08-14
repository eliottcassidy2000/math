# Disconnected-low affine tail and canonical reflected-residue branch closure

**Status.** PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Canonical theorem: [THM-3372](../01-canon/theorems/THM-3372-independent-superunit-affine-tail-and-reflected-residue-closure.md).  This is an independent alternative to THM-3355's weighted-horn proof; THM-3360 later sharpens its pair floor to `1/294`.  LRC(14) remains OPEN.  This note records the centered superunit repair, the exact finite dovetail, and the redundant canonical reflected-residue branch closure.

## 1. Scope and the missed chamber

On a nonzero affine resonance ray write

```text
D=d+a,                 1<=d<=8, 0<=a<=7d,
p=p0+dn, q=q0+Dn,      dq-Dp=c,
z=Lp-e, w=Lq-f,
A=dw-Dz=Lc+eD-df,
r=p mod d, N=(p-r)/d,
rho=|A|/z, T=N rho,    T_infinity=|A|/(Ld).
```

The refined turn-band theorem closes every finite `T>=1` point in the
analytic tail.  If `T<1`, it gives `1<=|c|<=13`; `c=0` is a separate fixed
ratio lane.  The earlier strict/equality split

```text
T_infinity<1, T_infinity=1
```

was not exhaustive.  For example

```text
d=2, a=0, c=2, p=781, q=782,
(L,j,e,f)=(784,420,2,1), A=1570, r=1
```

has

```text
T=390*1570/(784*781-2)=612300/612302<1,
T_infinity=1570/1568>1.
```

This is a real primitive physical point (`gcd(p,q)=1`), not a vacuous
parameter row.

## 2. Finite `T<1` forces only a slightly superunit limit

From `T<1`, with `A!=0`,

```text
|A|(p-r) < d(Lp-e).
```

Therefore

```text
T_infinity=|A|/(Ld)
 < (p-e/L)/(p-r)
 < p/(p-r).                                             (1)
```

For `p>=679`, `0<=r<=d-1`, and `d<=8`, the last expression is maximized
at `p=679,r=7`, so

```text
T_infinity < 679/672 = 97/96.                           (2)
```

The inequalities are strict even on the maximizing integer corner.  The
finite-exact census below harmlessly includes the closed cap
`T_infinity<=679/672`.

## 3. Expanded centered-continuum floor

For the limiting scaled tooth put

```text
Z0=L, W0=L(D/d),
F_(Z0,W0)(x)=(L/W0) int_0^(W0/(7Z0)) 1_C(x+u) du,
C=[-1/14,1/14] mod 1.
```

The convolution interval is centered.  In the integer compiler the phase is
`X-U/2`, not the left endpoint `X`.  For `A!=0`, integrate the exact
periodic piecewise-linear `d`-block path for length `T_infinity`, in the sign
of `A`; for `A=0`, evaluate the repeated zero-step block directly.

The exact universe is:

```text
all 2,530 small-ruler ordered contexts (L,j,e,f);
1<=d<=8, 0<=a<=7d;
nonzero |c|<=13 with gcd(a,d)|c;
the two endpoint-cone sign restrictions;
|A|/(Ld)<=679/672.
```

This deliberately overcounts rows that have no admissible raw residue or
`gcd(p,q)<=3` point.  Integer cross-products check `5,053,047`
triple-context rows, of which `362,926` have `T_infinity>1`.  The result is

```text
J0 >= 709/48048.                                        (3)
```

There are two equality rows; one witness is

```text
(d,a,c;L,j,e,f)=(3,8,-1;168,90,12,1), A=-39.
```

All `11,314` rows on `T_infinity=1` equal the mean `1/49`.  There are 29
zero-step triple-context rows.  Optimized and unoptimized C++ builds give
byte-identical output.

## 4. One uniform finite-to-continuum error

On the analytic superunit tail `p>=679`, put

```text
Z=L-e/p, W=L(D/d)+(Lc/d-f)/p,
zeta=3167/3168, omega=3155/3168,
k_f=6864/22085,
delta_d=d-1+1/12,
tau_d=1+delta_d/(264 zeta),
KP_d=d k_f+(97/96)[d(d-1)/(2 zeta)+d/zeta].                      (4)
```

The factor `97/96` on both drift terms is load-bearing: the exact
coefficients are `|A|/(dL zeta)=T_infinity/zeta`, once for the
within-block phase drift and once for the start-phase drift.  It cannot be
replaced by one on the superunit face.  The `d k_f` tooth-shape term is
unchanged.

The centered one-tooth comparison gives `||F-F0||_infinity<=k_f/p`;
phase drift across a `d`-block and start-phase drift give `KP_d/p`.  Also

```text
|T-T_infinity| <= T_infinity*delta_d/(p zeta),
max(T,T_infinity) <= tau_d T_infinity.                  (5)
```

For `p>=679`, the comparison path is shorter than two turns.  Indeed both
the cap in (1), refined by `d`, and `tau_d` increase with `d`, and at `d=8`

```text
(679/672) tau_8 = 26287/25336 < 2.                      (6)
```

Apply the signed composite-trapezoid identity to the finite block path and
compare it with the limiting path.  A safe five-term envelope is

```text
path-function perturbation        tau_d KP_d/d,
scaling/path-length displacement  2 tau_d delta_d/(7 zeta omega),
path endpoint displacement        tau_d d/(42 omega),
endpoint half-difference           d/(14 omega),
Peano slope-jump loss              d^2/(264 zeta omega). (7)
```

The second term deliberately retains an extra factor `tau_d`; cancellation
would permit its removal, but the inflated version avoids a face-dependent
derivation.  The Peano term does **not** require another `97/96` factor.  It
is controlled directly from the finite hypothesis `0<T<1`, rather than from
the length of the comparison path.  The physical signed Peano error is at most

```text
|A| L d /(2 Z W p^2).
```

Equivalently, after extracting the global factor `1/p`, its coefficient is
at most `|A|Ld/(2ZWp)`.  Since `T=(p-r)|A|/(dpZ)<1`, this coefficient is strictly less than
`d^2/[2 omega (p-r)]`.  For `p>=679` and `r<=7`, it is at most
`d^2/(1344 omega)`, hence certainly at most the displayed
`d^2/(264 zeta omega)`.  That term is therefore a deliberately inflated cap
independent of `T_infinity`.  Thus every nonzero-step finite `T<1` point in
this tail satisfies

```text
I(p) >= J0-K_d^*/p,                                    (8)

K_d^*=tau_d KP_d/d
      +2 tau_d delta_d/(7 zeta omega)
      +tau_d d/(42 omega)+d/(14 omega)
      +d^2/(264 zeta omega).
```

Exact arithmetic shows `K_d^*` increases for `1<=d<=8`, with

```text
K_8^*=1792138785426/221510098565,
K_8^*/(709/48048-Dmax/5)
 =144679594655915839062972000/
   207178827309738451742041
 in (698,699).                                          (9)
```

Consequently every nonzero-step `T<1` point with `p>=699` has
`I(p)>Dmax/5`.  The exact margin at `p=699` is

```text
46135211197112901571553/
4166631528336272997191190000 > 0.                       (10)
```

## 5. Zero step and zero resonance

If `A=0`, division by `|A|` is unnecessary.  Each complete `d`-block repeats
one value.  Comparing that value with its limit and discarding at most `d-1`
final teeth gives

```text
I(p)>=J0-K0_d/p,
KP0_d=d k_f+d(d-1)/(2 zeta)+d/zeta,
K0_d=KP0_d/d+(d-1)*2/(7 omega).                        (11)
```

Here `T_infinity=0`, so the superunit `97/96` factor in `KP_d` is neither
needed nor used; `KP0_d` is the zero-step specialization of the underlying
finite-to-limit comparison.

These constants also increase with `d`, and

```text
K0_8=477044832/69943195,
K0_8/(709/48048-Dmax/5) in (588,589).                   (12)
```

Hence the `A=0` face closes for `p>=589`.

If the Dirichlet resonance itself has `c=0`, reduction to coprime `(P,Q)`
gives `P|d`, hence `P<=8`.  In the remaining common-scale bank
`g=gcd(p,q)<=3`, this forces raw `p=gP<=24`; therefore no `c=0` point occurs
in the raw tail `p>=264`.  It belongs to the already certified finite head.

## 6. Dovetail

For the small-ruler, moderate-ratio, non-`3:5`, raw-gcd-`<=3` bank:

1. the finite exact scans cover every raw pair through `p=698`;
2. for `p>=699`, the refined turn bands close every `T>=1` Dirichlet
   witness;
3. (8)--(10) close every `A!=0,T<1` witness;
4. (11)--(12) close every `A=0,T<1` witness; and
5. `c=0` cannot occur there.

All three bridge segments have now received all-argmin reference audits and full-context hostile controls.  Therefore no affine-tail chamber remains.  This proves the strict high-pair floor used by THM-3372.

## 7. Reproduction

```bash
clang++ -std=c++17 -O3 -Wall -Wextra -Wconversion -Wshadow -Werror \
  04-computation/lrc14_affine_Tlt1_superunit_continuum_20260812.cpp \
  -o /tmp/lrc14-affine-expanded
git show HEAD:04-computation/lrc14_disconnected_head263_contexts_20260812.txt \
  | /tmp/lrc14-affine-expanded

python3 04-computation/lrc14_affine_Tlt1_superunit_unified_K_20260812.py
python3 -O 04-computation/lrc14_affine_Tlt1_superunit_unified_K_20260812.py
python3 04-computation/lrc14_affine_Tlt1_superunit_continuum_independent_audit_20260812.py
python3 -O 04-computation/lrc14_affine_Tlt1_superunit_continuum_independent_audit_20260812.py
```

An independent Python implementation constructs the periodized convolution
directly from its breakpoints (it does not call the C++ floor-area formulas),
then compares exact rational values on a deterministic reservoir sample of
`256` ordinary, `768` superunit, and all `29` zero-step rows.  All `1,053`
comparisons agree; ordinary and optimized outputs are byte-identical.

Frozen artifact hashes (LF-normalized bytes where relevant):

| artifact | SHA-256 |
|---|---|
| superunit continuum C++ | `2db0af3e0e2b678ee6b13d79dbaf4f6b8f91c8e8a83ee0a1d201de344abaf06a` |
| continuum output (O3/O0/UBSan identical) | `af1c6c944ea1e5997ed58809e4f8ddccde16326d4b90b39e3b777f4200498661` |
| independent PeriodicPL audit | `3d10014f093b2cbee89ba524fadd2e9a5b4b1dee92d28bc550da0869a18a2aa0` |
| independent output (normal/optimized identical) | `62a5eae88a76216fece5de8844f84a147b65392d8a77474e4cc88174096ab7c6` |
| unified exact-constant verifier | `db8575364318a66a8635ea9acebc80d1537520eaedfb880045a906017af576fb` |
| verifier output (normal/optimized identical) | `427cec0fe83780fb3115fa0c4ccfc8d21d82c6d645d692122bc9abc97375bda7` |

## 8. Inheritance, tree synthesis, and global consequence

The closest inherited objects are THM-2941's upper-median body bank and debt, THM-3350's all-channel tree transport, and THM-3352's exact arbitrary-ratio mass engine.  The disconnected low graph has a complete-multipartite high complement.  Choosing one vertex in each of two components gives an explicit five-edge cross-component spanning tree.  Since every high edge now has physical overlap strictly greater than `Dmax/5`, its credit strictly exceeds the whole singleton debt and Hunter closes the disconnected-low case.

Six distinct levels have connected or disconnected low graph.  THM-3352 closes the connected case and THM-3372 independently closes the disconnected case.  Repetitions are handled separately: THM-2941's same-level graph is complete on every one of the active 649 bodies.  Its only two chromatic exceptions have all fifteen robust edges and are disjoint from the robust-edge-`<=10` active bank.  Hence every positive canonical reflected-residue level assignment closes on those 649 bodies.  Together with the complementary bodies already closed by THM-2941, this exhausts all 3,003 bodies.  THM-3355 and THM-3360 give independent canonical proofs of the same branch conclusion.

## 9. Exact bridge provenance

The integer-only bridge covers `264<=p<=698`, `p<q<8p`, `gcd(p,q)<=3`, and all 2,530 small-ruler ordered contexts.  It evaluates `1,211,966` channels and `3,066,273,980` physical masses, with zero failures.  The weakest row is `(p,q;L,j,e,f)=(698,2559;168,90,12,1)`, with mass `20682154/1400220127>Dmax/5`.

The combined ledger SHA-256 is `932057abf10f674e4bb31f334c1ea94f39e4e17627c6a950bd5d727f8e595186`.  It is byte-identical to the concatenation of the separately generated `264..454`, `455..678`, and `679..698` ledgers.  Every one of the `1,211,966` reported argmins is replayed by the independently implemented slower THM-3352 reference engine; deterministic channels also receive full 2,530-context controls.

## 10. Corrected near-miss lineage

Three discarded compressions shaped the final theorem.  The first continuum compiler used a left endpoint instead of the centered convolution.  The second inferred `T_infinity<=1` from finite `T<1`; the explicit primitive witness in Section 1 refutes that implication.  The first superunit repair then put `97/96` on only one of the two phase-drift terms.  Both require it, moving the rigorous tail start from the provisional 695 to 699.  The exact bridge was extended through 698 before promotion.  See MISTAKE-378.

## 11. Scope and next concept board

This closes `z_e=q_eL-e` for positive levels, including the former 561-body residual.  It does not type arbitrary six-drift `k=1` packets into that residue form; other residues, projected `k=2,3`, the rung, physical entry, and LRC(14) remain open.

The next board has three live maps: (i) residue perturbation measured by symmetric-difference transport of a frozen Hunter tree; (ii) projected transfer, retaining denominator, unit, ray height, and local fibre coordinate; and (iii) the 14,168-ray primitive quotient as a finite-state object rather than a raw 22,890-row witness ledger.  Each map must name which physical coordinate it preserves and give a hostile zero-overlap probe before importing the pair floor.
