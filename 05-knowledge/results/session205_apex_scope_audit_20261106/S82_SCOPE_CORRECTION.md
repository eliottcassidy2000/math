# S82 one-apex join lead: scope audit and repaired quantitative obligation

Status: REVIEWED SCOPE CORRECTION (independent review by audit_threshold_leads, 2026-11-06). The fixed-parameter first-moment mechanism is checked below; no all-order coverage theorem is claimed. The reported finite coverage was inspected, not independently rerun here.

Source inspected on 2026-11-06: main/05-knowledge/results/one_apex_join_representation_mass_pi2_over_12_20261106.md, blob ce92ff847f40759f079a5cc01a9da455952c203d. Original checkpoint 5bde18fade89dc103d809c1e9df8462a2c0c5d50. Its audit and 2^35 coverage sections still contain placeholders. Do not report the 2^35 run as complete.

## 1. What survives

Write g(0)=1, g(2v+1)=g(v), g(2v)=g(v)+v+(v mod 2), and e(v)=[v=0 or v even]. For fixed a,c>=1 and b>=0, the count of v>=0 with a g(v)-b e(v)+c(1+v)<=X is

    (X/a) log(1+a/c) + O_{a,b,c}(sqrt X).

The bit-reversal grid proof is valid in mechanism: within a dyadic block the high half of the bits controls position and the reversed low half controls rho, with one sample in each product cell. A monotone boundary intersects O(sqrt(block size)) cells. The pointwise error g(v)-v rho(v)=O(log(v+2)) is harmless after summing blocks and sandwiching.

For a fixed valid rooted bank Y, summing this law gives the stated mean multiplicity, not a pointwise lower bound. The F1 leading coefficient pi^2/24 follows by summing log(1+t)/t over the bit-reversal grid. To justify the asserted linear constant for all real X, supplement the proof's complete dyadic-block estimate with its identical clipped-block estimate; summable complete-block errors alone would not justify arbitrary cutoffs. The triangular subtraction also needs a joint-grid limiting argument. This note does not substitute a complete independent proof of the second-order asymptotic for the pending S82 audit.

## 2. Uniformity correction

The displayed error O(sqrt(X/c)+(a/c)log X+1) is not uniform for arbitrary b>=0. With a=c=1, X=1 and b tending to infinity, every positive even v<=b/2-1 satisfies g(v)+v+1-b<=1 because g(v)<=v+1. Thus the count is unbounded while the displayed bound is constant.

A safe uniform replacement, for X>=1, is

    O(sqrt((X+b)/c) + (a/c) log(2+(X+b)/c) + b/c + 1).

Here is a uniform justification. Let L=(X+b)/c and D=a(2 log_2(2+L)+1)+b+c. Every counted v is at most L. For 1<=v<=L, the approximation error between Phi(v) and v(c+a rho(v)) is at most D. Let C0(Y) count positive v with v(c+a rho(v))<=Y, setting C0(Y)=0 for Y<=0. The preceding grid argument gives, with an absolute implied constant,

    C0(Y)=J(Y)+O(sqrt(max(Y,0)/c)+1),
    J(Y)=integral_0^1 (Y/(c+a s)-1)_+ ds.

The lower sandwich C0(X-D)<=C(X) is valid because X-D<=X+b, so its candidate v lie in the horizon L. The upper sandwich C(X)<=C0(X+D)+1 allows the extra v=0. The function J is 1/c-Lipschitz on the real line, and its difference from (X/a)log(1+a/c) is at most 1 for X>=0. Thus the total error is O(sqrt((X+D)/c)+D/c+1), which is absorbed by the displayed replacement because sqrt(D/c)<=D/c+1. All implied constants are absolute, and the logarithm base affects only a constant factor. This correction does not invalidate the F1 application, where b=e(u)<=1 and a/c<=1.

## 3. The actual covering obligation

Let Phi_y(v)=y_h g(v)-y_e e(v)+y_R(1+v). For a finite bank Y_N, truncate v to 0<=v<=floor((N+y_e)/y_R), which includes every representation of N. Put

    F_N(theta)=sum_{y in Y_N} sum_v exp(2 pi i theta Phi_y(v)).

Fourier inversion gives the exact identity

    r_{Y_N}(N)=integral_0^1 F_N(theta) exp(-2 pi i N theta) dtheta.

For a chosen measurable major-arc set A_N, an adequate all-order theorem would prove that the REAL PART of its integral is at least L_N>0, and establish the bound

    integral over [0,1)\A_N of |F_N(theta)| dtheta < L_N,

with major-arc errors already included in L_N. This suffices directly for r_{Y_N}(N)>0.

The S82 individual bound |sum_v exp(2 pi i theta(a g(v)+cv))| << V^(1-eta), even if available uniformly over all required coefficients, plus M_{Y_N}->infinity does not imply this inequality. A triangle bound sums errors proportional to lengths raised to 1-eta, whereas each bank member contributes only a constant-order mean representation count; increasing bank size alone does not repair that scale mismatch. One needs cancellation across bank members, stronger integrated estimates, or a different exact covering mechanism. The parity term and all local congruence factors must also be retained. Finite tests modulo 3,9,27,81 are not a uniform lower bound for a singular series.

## 4. Scope boundaries

- The source explicitly reports F2 coverage through 2^30, excluding 7 and 21, but its 2^35 result is pending.
- Neither the first moment nor a decreasing observed miss fraction proves eventual coverage.
- The Zaremba comparison is a heuristic analogy. No reduction establishing comparative hardness is supplied.
- Failure of bounded-defect or append-only compilers does not exclude a new all-order constructive family.

No change to the eventual odd-spectrum threshold is established by this audit.

## 5. Review and finite controls

Independent reviewer audit_threshold_leads checked the counterexample and repaired uniform-bound algebra, conditional on the stated absolute grid estimate, and the Fourier scope diagnosis. The review requested explicit real-part language for a nonsymmetric major-arc set; this is incorporated above. No independent rerun of the historical 2^30 coverage is claimed.

Reproduce the small controls with `python3 check_scope.py` from this directory. The script checks the rho approximation for 1<=v<4096, five exact b-growth controls, and 100 exact no-alias Fourier coefficient controls. These controls test the formulas and counterexample only; they are not evidence for an asymptotic coverage claim.
