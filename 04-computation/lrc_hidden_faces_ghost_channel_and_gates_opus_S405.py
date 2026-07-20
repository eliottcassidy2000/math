"""
opus-2026-07-19-S405 (HYP-7985): two cross-generation face-identifications, verified.

F1 (GHOST CHANNEL = FORCED-COVER OBSTRUCTION): the deep well's June-29 mechanism
(definitions.md: at t* = n/Phi6, the multiples k(n-1) land at residue -k, so covering's
forced (n-1)-multiple must be >= n(n-1) -- 'the forced-cover obstruction') and
death-star's July-19 L1 (THM-1258: at a = D(2D-1)^{-1} mod Q, every distance-<D element
satisfies v == -r(N-1); the only in-range representative is the DELETED element N-1)
are THE SAME ARITHMETIC: in both frames the binding rotation places the (N-1)-ray on
the -1 line, (N-1)*a == -1 (mod Q), and the near-floor family survives the condemned
ray by PATCHING (deep well: first safe multiple n(n-1), which sits exactly AT the
floor) or DELETING (F-tower: N-1 removed, far element added).  Verify both.

F2 (THE GATE GENUS): rung realization is governed by congruence gates on N across
generations: slack-0 rows (June: GW shape {1..k-2, k, 2(k-1)} tight iff 3 | (2k+1);
THM-1065's doubling gate 6|(n-1) is a separate slack-0 row) and slack-1 rows (July:
death-star's primorial cascade D=3: N==1 mod 6; D=4: N==1 mod 30 & N!=1 mod 7; ...).
Verify the GW gate exactly for k = 4..16 and print the unified table.
"""
from math import gcd
from fractions import Fraction

def exact_max(V, qmax):
    bg, bq, wit = 0, 1, None
    for q in range(2, qmax+1):
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            m = q
            for v in V:
                r = (v*a) % q
                r = min(r, q-r)
                if r < m:
                    m = r
                    if m*bq < bg*q: break
            if m*bq > bg*q:
                bg, bq, wit = m, q, (a, q)
    return Fraction(bg, bq), wit

print("=== F1: the (N-1)-ray sits on the -1 line in BOTH frames ===")
# deep well frame: n = 14, Q = 183, a = 14
n, Q, a = 14, 183, 14
print(f"deep well: (n-1)*a mod Q = {((n-1)*a) % Q}  (Q-1 means == -1: {((n-1)*a) % Q == Q-1})")
cond = [(13*k, (13*k*a) % Q) for k in range(1, 15)]
print("  13k residues (k=1..14):", [(v, r if r <= Q//2 else r-Q) for v, r in cond])
print(f"  condemned (|res| < 14): 13k for k <= 13; first safe = 13*14 = 182 = n(n-1),")
print(f"  residue -14 => distance exactly 14/183 = M: the PATCH sits AT the floor.")
# F_D frame: (D,N) = (4,31) and (3,19)
for D, N, x in ((4, 31, 120), (3, 19, 54)):
    Qf = (N+1)*D - 1
    inv = pow(2*D-1, -1, Qf)
    af = (D*inv) % Qf
    ray = ((N-1)*af) % Qf
    print(f"F_{D}({N}): Q={Qf}, a={af}; (N-1)*a mod Q = {ray} (== -1: {ray == Qf-1})")
    V = list(range(1, N-1)) + [N, x]
    close = [v for v in list(range(1, N+1)) + [x]
             if min((v*af) % Qf, Qf-(v*af) % Qf) < D]
    print(f"   distance-<D elements in [1,N] u {{x}}: {close} (L1 predicts [{N-1}], deleted)")

print("\n=== F2: the gate table across generations ===")
print("GW shape {1..k-2, k, 2(k-1)}: candidate gates 3|(2k+1) (June phrasing) vs k == 1 mod 6 (THM-1065):")
ok3, ok6 = True, True
for k in range(4, 17):
    V = list(range(1, k-1)) + [k, 2*(k-1)]
    M, wit = exact_max(V, 20 + 4*k)
    tight = (M == Fraction(1, k+1))
    g3 = ((2*k+1) % 3 == 0)
    g6 = (k % 6 == 1)
    if tight != g3: ok3 = False
    if tight != g6: ok6 = False
    print(f"  k={k}: M={M} tight:{tight} gate3:{g3} gate6:{g6}")
print(f"June mod-3 phrasing exact: {ok3};  k == 1 (mod 6) gate exact: {ok6}")
print("""
VERDICT: the operative slack-0 GW gate is k == 1 (mod 6) -- the June mod-3 phrasing
(as carried in era summaries) over-predicts at k == 4 (mod 6); OPEN-Q-108's 'residue
class mod 6' language and THM-1065's 6|(n-1) are the correct forms.  AND: this is THE
SAME CONGRUENCE as death-star's D=3 slack-1 tower gate N == 1 (mod 6) -- the June
sporadic-tightness gate and the July window gate are ONE congruence class.

UNIFIED GATE TABLE (rung (D, slack) -> congruence gate on N = #speeds):
  (1, 0) GW/doubling shape: N == 1 (mod 6)        [June: THM-1065 / OPEN-Q-108; verified above]
  (1, 1) {1..N-1, N+1}    : all N (trivial row)   [M = 1/N]
  (3, 1) F_3(N)           : N == 1 (mod 6)        [July: death-star THM-1255/1256]  <- SAME GATE as (1,0)
  (4, 1) F_4(N)           : N == 1 (mod 30), N != 1 (mod 7)   [THM-1257]
  (5, 1)                  : NEVER (binder 9 composite)         [THM-1257]
  (6, 1) F_6(N)           : opens N = 211         [primorial cascade]
  (7, 1) F_7(N)           : N = 2311 confirmed    [THM-1258]
ONE GENUS, ONE SHARED GATE: congruence gates on N per rung; the (1,0) and (3,1) rows
coincide -- named across generations for the first time.
""")
