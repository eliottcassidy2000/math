"""aud_density.py -- audit S11: covering density D_k(t) of HYP-9471 from the audit's own matrix files.

D_k(t) = sum over the DISTINCT nonsingular 2x2 restrictions [[a, b], [c, d]] = [[M00, M02], [M20, M22]] of the
S-row-zero matrices of order <= k of (a + b t)^(-2).
Usage: python3 aud_density.py DIR KMAX     (reads DIR/aud_mats1.txt ... aud_matsKMAX.txt)
Also prints: the same sum taken over distinct 3x3 matrices instead of distinct restrictions (sensitivity check),
the uniform slope average of D_k over [2, T] for a few T, and the growth ratios D_k(2) / D_(k-1)(2).
"""
import sys, math
D = sys.argv[1]; K = int(sys.argv[2])
first = {}
for k in range(1, K + 1):
    for line in open(f'{D}/aud_mats{k}.txt'):
        m = tuple(map(int, line.split()))
        first.setdefault(m, k)
ts = (2.0, 2.5, 3.0, 4.0)
prev = None
for k in range(1, K + 1):
    s0 = [m for m, o in first.items() if o <= k and m[3] == 0 and m[4] == 0 and m[5] == 0]
    restr = {}
    for m in s0:
        A = (m[0], m[2], m[6], m[8])
        if A[0] * A[3] - A[1] * A[2] != 0:
            restr.setdefault(A, k)
    oddd = [A for A in restr if (A[0] * A[3] - A[1] * A[2]) % 2]
    Dk = {t: math.fsum((A[0] + A[1] * t) ** -2 for A in restr) for t in ts}
    D3 = math.fsum((m[0] + m[2] * 2.0) ** -2 for m in s0 if m[0] * m[8] - m[2] * m[6] != 0)
    ratio = f"{Dk[2.0] / prev:.2f}" if prev else "-"
    prev = Dk[2.0] if Dk[2.0] > 0 else None
    print(f"k<={k}: S-row-zero distinct {len(s0)}, distinct nonsingular 2x2 restrictions {len(restr)}, odd det {len(oddd)}; "
          + ", ".join(f"D({t})={v:.4f}" for t, v in Dk.items())
          + f"; growth D(2) {ratio}; [over distinct 3x3 instead: D(2)={D3:.4f}]")
    if k == K:
        for T in (2.5, 3.0, 4.0, 5.0):
            N = 2000
            avg = math.fsum(math.fsum((A[0] + A[1] * (2 + (T - 2) * (i + 0.5) / N)) ** -2 for A in restr) for i in range(N)) / N
            print(f"   uniform slope average of D_{k} over [2, {T}]: {avg:.4f}")
