# Independent hostile audit of THM-2791

## Verdict

**PASS.  PROMOTION HAS LANDED; ADD THE SAME-ANCESTRY REFINEMENT.**

No mathematical or scope correction remains after the already-landed
Fourier-origin and raw-unit wording repairs.  The full address scan,
two-point physical chains, periods and conductor claims, THM-2788
lower-central word, quotient pushforward, normalized group-ring unit, and
transfer-versus-descent boundary all survive an independent reconstruction.

The correct status is:

```text
PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED
```

within the theorem's stated physical clock-one coefficient universe.  It is
not an LRC row exclusion or the missing endpoint-origin allocation.

After this audit, `origin/main` promoted THM-2791 and landed the
whole-cylinder endpoint check.  At fetched commit `8e57679b48`, the promoted
LF-normalized bytes are

```text
cebde0deb9dd1658303339f883b41aed1b578e31a32f66fb9b29f8a30fabb572
  01-canon/theorems/THM-2791-full-arm-orbit-transfer-and-lower-central-chord.md

2824d62c237fd9ac831d23236e6987ecabe96bebd68ba37a9abd0bb685ad0716
  04-computation/lrc14_full_arm_orbit_lower_central_chord_thm2791.py

9f2b8e69b9de430f201adb7758f98fff7bf505c5bd792b03f40b3ee7c9f46edd
  05-knowledge/results/lrc14_full_arm_orbit_lower_central_chord_thm2791.out
```

The changes are status/audit prose, a sharper open-boundary statement, and
an explicit whole-cylinder endpoint test; they do not alter the physical
construction audited below.

## Audited bytes and replay

LF-normalized SHA-256:

```text
7092c53aa1ac567dfe397f600c6cfd8cfe23395d5630da0dd2dbe6c7f096a4dd
  01-canon/theorems/THM-2791-full-arm-orbit-transfer-and-lower-central-chord.md

fc557d31a82ea52adf0abc7b26bfdacf1961facc99050061132d72dda231a9db
  04-computation/lrc14_full_arm_orbit_lower_central_chord_thm2791.py

9f2b8e69b9de430f201adb7758f98fff7bf505c5bd792b03f40b3ee7c9f46edd
  05-knowledge/results/lrc14_full_arm_orbit_lower_central_chord_thm2791.out
```

Fresh ordinary and `python3 -O` executions byte-match each other and the
stored output.  The source has zero AST `assert` nodes.

Independent audit artifacts:

```text
59d297e0fc79c046b281ce81d340e0476018bf1fd1f11c5283f9a0f5eea3794f
  .scratch/thm2791_hostile_audit/audit.py

5a21a1aaf614c473a9c7bdaf3d94d76bb312151c6710f1b2576ccf442bb33d6b
  .scratch/thm2791_hostile_audit/audit.out

5a21a1aaf614c473a9c7bdaf3d94d76bb312151c6710f1b2576ccf442bb33d6b
  .scratch/thm2791_hostile_audit/audit.opt.out

8bec0f61e17c0daa5e460a4b63328b22ba3c7b68f050cea84443af36020ee492
  .scratch/thm2791_hostile_audit/path_sheet_probe.py

2325aa126dc4a97af4ba4bde0348389bfdfe5720dec6f88792b9a06baa40afd3
  .scratch/thm2791_hostile_audit/path_sheet_probe.out

2325aa126dc4a97af4ba4bde0348389bfdfe5720dec6f88792b9a06baa40afd3
  .scratch/thm2791_hostile_audit/path_sheet_probe.opt.out

b2457337c00e92cda715a5a03daae789e0ea62060a3fcd8a8d3e2fe48fac11d4
  .scratch/thm2791_hostile_audit/PROPOSED-PATH-SHEET-ADDENDUM.md
```

The independent script does not import the THM-2791 companion.  It has zero
AST `assert` nodes, compiles, and its ordinary and optimized transcripts
byte-match.  The same statements hold for the separate path-sheet probe.

## Dependency-hash gate

All declared direct dependency/control hashes match current LF bytes:

```text
THM-2782 script
  7fbc6bb1ec303ded98eaad6e5d8205eb3d247258ada32b6f9904fc439ebb11fb
THM-2782 output
  13f570d63f212808171cecdb4d8f9aa41884fbdc7ed571dbfe27122b412fadc4

THM-2779 primary script/output
  4c6a58c80ddd4be0fd9bdd297b310df054bbc08996eb223d519d3cce6b8ed13a
  f7c96259777a3ab4a5e46cac8666181ae77a3be2e440cee8785997507706791a
THM-2779 secondary script/output
  004e06c617f9305e2f0bc30871926e3faa7843f47dcf63af1fd8a892e63101e4
  a89c00a3830ee9ff282cc5e4557d41293af5d6f0e7feabd5d3c7e7808591e754
THM-2779 tertiary script/output
  ef6e9f9bcb4f11152d291342a11ae215245d1d19b96c49940a01ba9ea850cbd9
  1feb463864015035ab8d7fcfcddf9cfe8b0ec0a3ed36481f2f66d6a9149182e6
THM-2779 independent script/output
  5019f87b24500a5a13825d3be01908ea983a08b360a384fd614107f476201f46
  1b9ad37b35e92a14dd90d0db8c1d0cf225761e2c37ca8e2fe2120bd0f64c47d4

THM-2788 script/output convention control
  d414bf2afb6aa3e40de9378ae20f03db1cb7bff75f59f13a60ac96e56cb95a89
  99ad33904617d45d76a285de5467b96408dc164839cb4168905c7fe678db8f66
```

The three direct source imports under the THM-2782 reconstruction were also
frozen and matched:

```text
d67c852c52f88feaadb2fcaa0a9a07a212f2e47018040b455855df886200595e
  lrc14_fully_marked_root_zero_target_profile_thm2749.py
25cbed38026d61891173c687006250a69fe38aea56d67439406bd8bb60fa2552
  lrc14_root_zero_clutch_mayer_vietoris_wing_shear_thm2751.py
f16754bd38ae0dfa0d7d91cc404b4447dbf359635101aa7b4223363f8064352f
  lrc14_relative_present_semantic_lift_probe_20260728.py
```

## 1. Complete physical universe and two-point chains

The independent audit rebuilds the physical row from THM-2782 operations
and uses a separate exact interval-overlap engine.  It checks all

```text
13^6 = 4,826,809
```

full address offsets for the common target profile.  The recovered support
is exactly

```text
k in {0,689364}
```

for each raw target `tau=3,...,11`; targets `0,1,2` are empty.  Equality of
the nine underlying weighted carriers makes the one full scan apply to all
nine columns, not only to `tau=3`.

At both support points:

```text
coefficient = c = 790161473087466480,
weighted mass = c/13^6
              = 60781651775958960/371293.
```

The absolute addresses and depth-two digits are

```text
3454614 -> (v,w)=(7,6),
4143978 -> (v,w)=(7,7).
```

Both whole open cylinders retain semantic record `(3,(1,2))`, predecessor
carry six, `sigma=0`, every target label `3,...,11`, and the sharp semantic
stability radius.  Their single weighted pieces are the entire selected
cylinders rather than truncated intersections.

The exceptional target `tau=12` was independently scanned through its full
`13^5=371293` central orbit.  It has exactly `121` positive cells, all with
the same coefficient and mass.

## 2. Periods, characters, and conductors

The translation stabilizers are trivial for:

```text
{0,689364} in C_(13^6),
{0,53028} in C_(13^5),
the exceptional 121-cell profile in C_(13^5).
```

Thus the claimed minimal periods `13^6` and `13^5` are exact.

For every full-address character,

```text
Fhat(chi)=c(1+chi(689364)^(-1)).
```

The value `chi(689364)` has order dividing the odd number `13^6`, so it
cannot equal `-1`.  The audit checked this exact order statement for all
`13^6` character indices, with no floating roots of unity.

Independent rational cyclotomic-block tests show survival at every primitive
central conductor

```text
13,13^2,13^3,13^4,13^5
```

for the main chord, the exceptional profile, and every one of the thirteen
decoded target profiles.

For `Z_r=O^(13^r)`, the five nonidentity finite differences have four signed
cells and normalized `l1=l2^2=4`; `Z_6` is the identity.  Restoring the
coefficient gives the theorem's `4c` and `4c^2` norms.

## 3. THM-2788 word and digit conventions

With `X:n->14n`, `O:n->n+1`, and
`[A,B]=ABA^(-1)B^(-1)`, direct affine calculation gives

```text
[X,O^(13^r)]=O^(13^(r+1)).
```

This matches THM-2788 exactly.

The full gap is

```text
689364=13*53028,
53028=(1,10,1,11,1)_13
     =1+10*13+13^2+11*13^3+13^4.
```

Since `Z_r=O^(13^r)`, the little-endian convention gives

```text
O^689364 = Z1 Z2^10 Z3 Z4^11 Z5.
```

Relative to the pure `Z1` step,

```text
53028-1=13*4079,       4079=10 mod13,
```

so the first discrepancy is exactly `10` in `Z2/Z3`.  It is divisible by
the first `u^13-1` factor and not by `u^169-1`.

## 4. Partial germ and quotient pushforward

The physical translation is exactly

```text
delta=7*689364/13^6=371196/371293.
```

It maps the entire first weighted cylinder to the second, with common weight
`27581135604`.  The exact integrality checks

```text
13^6 delta=4825548,
13^5 delta=371196
```

preserve the delayed phase, carry digit, and selected root-half data on this
restricted germ.  This is a genuine partial transfer with gain
`Z1^53028`; it is not global packet covariance.

The later literal-sheet challenge gives a strictly stronger local statement.
Rail eight is the route-two product

```text
U(x)V(x-1/13),
U=P_(13^5)(1_Q P_169 1_E),       V=P_(13^5)1_E.
```

Its pre-marginalization Boolean labels are `(a,b)` for `U` and `e'` for
`V`.  Independently enumerating their defining interval memberships at the
two translated cylinders gives

```text
#(a,b)=966606,       #e'=28534,
966606*28534=27581135604.
```

More strongly, the two ordered label sets are literally equal, not merely
equinumerous.  Both cylinders lie strictly inside the same raw contributor
chamber

```text
[140890500190440,144190879112280),
```

with no `Q`, `E`, or rotated-`E` path wall between them.  Thus translation
acts by the identity on every one of the `27,581,135,604` product labels.
The fixed collision labels are `(u,v,s,t)=(5,6,1,12)`, and the supplied
positive copy

```text
(a,b,e')=(59162,26,56658)
```

is active on both cylinders.  Its common ordered label-set digest is

```text
15c804c7cea9f61feab3b641eccdc035d937142b446d1cc14e059210eb1534fd.
```

The separately found same-atom interval
`[142004992589460,142005019034340)` lies in this same chamber as well.
Consequently THM-2791's partial germ may be sharpened to a
**same-ancestry partial germ on the THM-2471 rail sheet**.  No label is lost
through `(u,v,a,b,e')`.  This still does not construct the later
THM-2779/THM-2625 endpoint-origin or endpoint-atom labels.

Reduction modulo `169` sends the two absolute addresses to

```text
85=7+13*6,       98=7+13*7,
```

and fibre summation gives

```text
q_!F=c z^6(1+z).
```

The group-ring calculation is correct:

```text
(1+z)(1-z+z^2-...+z^12)=2
```

in every `K[C13]` with `char(K)!=2`.  Hence the scalar-normalized chain is a
unit.  The inherited coefficient Bockstein is also valid:

```text
v_13(c)=1,       (c/13) mod13=2,
```

so the characteristic-thirteen normalized chain is a unit as well.

The raw integral chain is **not** a unit in `Z[C13]`; THM-2791 now says this
explicitly.  The earlier repaired wording is therefore sharp.

The exceptional pushforward is exactly

```text
9N13+(1+z+z^2+z^3).
```

The integral `K_beta` convention agrees with THM-2779:

```text
(A_q,B_q)=
((58,5),(64,5),(59,5),(55,7),(48,12),(51,2),(48,8),
 (53,2),(56,9),(50,8),(47,11),(49,0),(55,3)).
```

Every decoded quotient row is nonconstant.  Over the rational prime-cycle
group ring, nonconstancy is equivalent to nonvanishing at every nontrivial
character, because a single primitive zero would force a multiple of
`Phi_13=1+z+...+z^12`.

## 5. Transfer is not descent or endpoint allocation

All three boundaries are stated and computationally witnessed correctly:

1. **Not descent.**  In offset coordinates,
   `F_tau(0)=c` but `F_tau(169)=0`, although the addresses agree modulo
   `169`.  The quotient object is a fibre **sum**, not a fibrewise value.
2. **Not pure `Z1` covariance.**  The remote physical lift is
   `Z1 Z2^10 Z3 Z4^11 Z5`, while the physical adjacent lift `j=1` is empty.
3. **Not endpoint-origin allocation.**  No map identifies these physical
   rail-sheet cylinders with THM-2779/THM-2625 endpoint origins or endpoint
   atoms.  The new identity on `(u,v,a,b,e')` removes the earlier
   rail-ancestry uncertainty, but it does not supply that endpoint sidecar.

The theorem also correctly withholds global packet covariance, positivity
of the signed inverse, current equivariance, THM-2542 root transition, row
exclusion, and an LRC(14) conclusion.

## Promotion recommendation

The promotion and independent-audit metadata have landed on `origin/main`.
Retain status

```text
PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
```

The theorem may additionally be strengthened, without changing its
dependencies, to call the graded-chord germ same-ancestry on the inherited
THM-2471 rail sheet.  Retain the endpoint boundary precisely: the missing
prize is a natural map from that now-fixed rail sheet to the endpoint-origin
central edge, not recovery of the rail ancestry itself.
