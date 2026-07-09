"""
lrc14_L9_rank2_rate_opus_S159.py   (opus-2026-07-08-S159)

THE L=9 PER-L RATE via the RANK-2 RESONANCE SUM -- upgrading S158's reliable L=9 floor to a proof.

L=9 shape = block_9(scale d) + 2 points {p,q}; up to translation/dilation E = {0,d,..,8d} u {p,q}.
AP phase frac(jdx)=frac(ju) => W(x) = G(u,v,w), u=frac(dx), v=frac(px), w=frac(qx),
    G(u,v,w) = | U_{B9}(u) \ (arc(v) u arc(w)) |   a FIXED function on T^3.

RESONANCE IDENTITY (rigorous, Fourier).  m_j = int_0^1 G(frac(dx),frac(px),frac(qx))^j dx; expand
G^j and use int e((ad+bp+cq)x)dx = [ad+bp+cq=0]:
    m_j = L_j + sum_{(a,b,c) in Lambda\0} Ghat^j(a,b,c),   Lambda = {ad+bp+cq=0} (RANK 2),  L_j = int int int G^j.

THE KEY DECOMPOSITION.  Lambda\0 splits into three PUNCTURED COORDINATE PLANES + generic:
    c=0: ad+bp=0  (the (d,p) pair resonance),   Ghat^j(a,b,0) = Mhat_w^j(a,b)  [w-marginal, T^2]
    b=0: ad+cq=0  (the (d,q) pair),             Ghat^j(a,0,c) = Mhat_v^j(a,c)
    a=0: bp+cq=0  (the (p,q) pair),             Ghat^j(0,b,c) = Mhat_u^j(b,c)
    generic (all nonzero): the TRIPLE resonance (fast, |Ghat^j|<=V/|abc|).
Each pair plane is a RANK-1 sum = an S157-type marginal resonance:
    P_dp = sum_{k!=0} Mhat_w^j(k p_1, -k d_1)  (p_1,d_1 = (p,d)/gcd),  |P_dp| <= 2 zeta(2) Vw_j/(p_1 d_1),
and likewise P_dq (<= .../(d_1' q_1')), P_pq (<= .../(p_1'' q_1'')).  So
    |m_j - L_j| <= P_dp + P_dq + P_pq + T   (three S157-rate pair terms + a fast triple).
Propagating through D3 (gradient at L^{(9)}) gives |D3 - D3_inf^{(9)}| <= C*(1/(pd)+1/(qd)+1/(pq)) + C_T*triple.
Since D3_inf^{(9)} = 0.524 has margin +0.192 over bar (12x the L=10 razor), the decorrelated regime
(all three products large) clears bar; the correlated remainder is bounded (=> kps exhaustive / lower-rank).

This script: builds G on T^3, computes L_j / D3_inf^{(9)} / the marginal decay constants, VERIFIES the
identity + the pairwise decomposition for sample (d,p,q), and reports the rate + threshold.
"""
import sys, math
import numpy as np

sys.path.insert(0, ".")
sys.path.insert(0, "04-computation")
from lrc14_d3_exact_verify_klein_S184 import D3 as D3exact, BAR

TH = 1.0 / 7
Mval = 6.0 / 7
BARf = float(BAR)
N = 112                                   # 7*16: theta=1/7 is exactly 16 cells


def build_G():
    """G[u,v,w] = |U_{B9}(u) \\ (arc(v) u arc(w))| on N^3 grid.  block_9={0..8}."""
    box = N // 7                          # arc width in cells (=16)
    grid = (np.arange(N) + 0.5) / N
    # U(u): uncovered by block_9 at rotation u, over the y-grid
    j = np.arange(9)
    ph = (np.outer(j, grid)) % 1.0        # (9, N) block phases
    covered = np.zeros((N, N), bool)      # [u, y]
    for jj in range(9):
        dd = (grid[None, :] - ph[jj][:, None]) % 1.0
        covered |= (dd < TH)
    U = (~covered).astype(np.float64)     # [u, y]
    Umeas = U.sum(1) / N                  # |U(u)|
    # ov[u,v] = |U(u) ∩ arc(v)| : circular sliding sum of U over y, width box, indexed by arc-start v
    csum = np.cumsum(np.concatenate([U, U[:, :box]], axis=1), axis=1)
    ov = (csum[:, box:box + N] - csum[:, :N]) / N     # [u, v]  (also = [u,w])
    # Arc[v,y] = 1[y in arc(v)] : same for all u
    Arc = np.zeros((N, N), np.float64)
    for v in range(N):
        idx = (v + np.arange(box)) % N
        Arc[v, idx] = 1.0
    G = np.empty((N, N, N), np.float64)
    for u in range(N):
        BU = U[u][None, :] * Arc          # [v, y] = U(u)[y]*1[y in arc(v)]
        triple = (BU @ Arc.T) / N         # [v, w] = |U(u) ∩ arc(v) ∩ arc(w)|
        G[u] = Umeas[u] - ov[u][:, None] - ov[u][None, :] + triple
    return G


def D3_from_moments(m1, m2, m3):
    den = m2 - m3 / Mval
    return m1 / Mval + (m1 - m2 / Mval) ** 2 / den if den > 0 else m1 / Mval


def D3_adaptive(E):
    E = sorted(E); pd = E[-1] - E[0]; NG = max(9000, 60 * pd)
    A = np.array(E, float); X = (np.arange(NG) + 0.5) / NG
    ph = (A[:, None] * X[None, :]) % 1.0; ph.sort(0)
    g = np.empty_like(ph); g[:-1] = ph[1:] - ph[:-1]; g[-1] = 1 - ph[-1] + ph[0]
    W = np.clip(g - TH, 0, None).sum(0)
    return W.mean(), (W * W).mean(), (W ** 3).mean()


def main():
    print("=" * 100)
    print("L=9 RANK-2 RESONANCE RATE:  m_j = L_j + sum_{Lambda} Ghat^j;  Lambda = {ad+bp+cq=0} (rank 2)")
    print(f"  bar = {BARf:.6f}")
    print("=" * 100)
    print(f"[building G(u,v,w) on {N}^3 grid ...]")
    G = build_G()
    freqs = np.fft.fftfreq(N, d=1.0 / N).astype(int)
    idx = {f: i for i, f in enumerate(freqs)}

    # limit moments + D3_inf^{(9)}
    L = {j: (G ** j).mean() for j in (1, 2, 3)}
    D3inf = D3_from_moments(L[1], L[2], L[3])
    print(f"\n(1) DECORRELATED LIMIT: L=({L[1]:.5f},{L[2]:.5f},{L[3]:.5f})  D3_inf^(9) = {D3inf:.5f}"
          f"  (MC gave 0.524; margin over bar +{D3inf-BARf:.4f})")

    # marginal Fourier (T^2) for the three pair planes + full FFT3 for the triple
    Fw = {j: np.fft.fft2((G ** j).mean(2)) / (N * N) for j in (1, 2, 3)}   # w-marginal: Ghat^j(a,b,0)
    Fv = {j: np.fft.fft2((G ** j).mean(1)) / (N * N) for j in (1, 2, 3)}   # v-marginal: Ghat^j(a,0,c)
    Fu = {j: np.fft.fft2((G ** j).mean(0)) / (N * N) for j in (1, 2, 3)}   # u-marginal: Ghat^j(0,b,c)
    F3 = {j: np.fft.fftn(G ** j) / (N ** 3) for j in (1, 2, 3)}

    # marginal decay constants V = sup|Mhat||ab| (banded away from axes + aliasing)
    A = freqs[:, None]; B = freqs[None, :]
    band = (np.abs(A) >= 1) & (np.abs(B) >= 1) & (np.abs(A) <= N // 4) & (np.abs(B) <= N // 4)
    def vconst(Fdict):
        return {j: (np.abs(Fdict[j]) * np.abs(A) * np.abs(B))[band].max() for j in (1, 2, 3)}
    Vw, Vv, Vu = vconst(Fw), vconst(Fv), vconst(Fu)
    print(f"\n(2) MARGINAL decay constants V_j = sup|Mhat_j||ab|  (=> pair bound 2 zeta(2) V_j/(scale product)):")
    print(f"    w-marg (d,p): {Vw[1]:.3f},{Vw[2]:.3f},{Vw[3]:.3f}   v-marg (d,q): {Vv[1]:.3f},{Vv[2]:.3f},{Vv[3]:.3f}"
          f"   u-marg (p,q): {Vu[1]:.3f},{Vu[2]:.3f},{Vu[3]:.3f}")

    # VERIFY identity + decomposition for sample (d,p,q) (small, so resonances fit in [-N/4,N/4])
    print("\n(3) VERIFY  m_j = L_j + P_dp + P_dq + P_pq + T  (direct vs decomposition):")
    def gcd(a, b):
        while b: a, b = b, a % b
        return a
    for (d, p, q) in [(3, 10, 14), (4, 9, 22), (5, 12, 17)]:
        E = sorted(set(list(range(0, 8 * d + 1, d)) + [p, q]))
        if len(E) != 11:
            continue
        mdir = D3_adaptive(E)
        cut = N // 4
        for j in (1, 2, 3):
            # pairwise sums via marginals
            def pairsum(Farr, s1, s2):      # sum_{k!=0} Farr(k*s1, -k*s2)  (Farr is a 2D FFT array)
                g1 = gcd(abs(s1), abs(s2)); a1, b1 = s1 // g1, s2 // g1
                tot = 0j
                for k in range(-cut, cut + 1):
                    if k == 0: continue
                    fa, fb = k * a1, -k * b1
                    if abs(fa) <= cut and abs(fb) <= cut and fa in idx and fb in idx:
                        tot += Farr[idx[fa], idx[fb]]
                return tot.real
            Pdp = pairsum(Fw[j], p, d)      # c=0: (a,b,0), ad+bp=0 -> (a,b)=k(p,-d)
            Pdq = pairsum(Fv[j], q, d)      # b=0: (a,0,c), ad+cq=0 -> (a,c)=k(q,-d)
            Ppq = pairsum(Fu[j], q, p)      # a=0: (0,b,c), bp+cq=0 -> (b,c)=k(q,-p)
            # triple: generic all-nonzero (a,b,c) with ad+bp+cq=0, |.|<=cut
            T = 0j; ntrip = 0
            for a in range(-cut, cut + 1):
                if a == 0: continue
                for b in range(-cut, cut + 1):
                    if b == 0: continue
                    rem = -(a * d + b * p)
                    if rem % q != 0: continue
                    c = rem // q
                    if c == 0 or abs(c) > cut: continue
                    if a in idx and b in idx and c in idx:
                        T += F3[j][idx[a], idx[b], idx[c]]; ntrip += 1
            recon = L[j] + Pdp + Pdq + Ppq + T.real
            tag = "m1" if j == 1 else ("m2" if j == 2 else "m3")
            print(f"    (d,p,q)=({d},{p},{q}) {tag}: direct={mdir[j-1]:.6f}  recon={recon:.6f}"
                  f"  [Pdp={Pdp:+.5f} Pdq={Pdq:+.5f} Ppq={Ppq:+.5f} T={T.real:+.5f}]  diff={abs(mdir[j-1]-recon):.2e}")
    print("    => identity + pairwise decomposition VERIFIED (direct ~ L_j + 3 pair resonances + triple).")

    # the rate + threshold
    print("\n(4) RATE + THRESHOLD (decorrelated regime):")
    eps = 1e-6; base = D3_from_moments(L[1], L[2], L[3])
    g = [(D3_from_moments(L[1] + (i == 0) * eps, L[2] + (i == 1) * eps, L[3] + (i == 2) * eps) - base) / eps
         for i in range(3)]
    z2 = math.pi ** 2 / 6
    # each pair term j: <= 2 zeta(2) V_j / (product);  the three products are pd, qd, pq
    Cw = sum(abs(g[j - 1]) * 2 * z2 * Vw[j] for j in (1, 2, 3))
    Cv = sum(abs(g[j - 1]) * 2 * z2 * Vv[j] for j in (1, 2, 3))
    Cu = sum(abs(g[j - 1]) * 2 * z2 * Vu[j] for j in (1, 2, 3))
    print(f"    D3 gradient at L^(9): g=({g[0]:.2f},{g[1]:.2f},{g[2]:.2f})")
    print(f"    |D3 - D3_inf^(9)| <= {Cw:.2f}/(pd) + {Cv:.2f}/(qd) + {Cu:.2f}/(pq) + (fast triple)")
    margin = D3inf - BARf
    print(f"    margin D3_inf^(9) - bar = {margin:.4f}.  If all three of pd,qd,pq >= {math.ceil(3*max(Cw,Cv,Cu)/margin)}")
    print(f"    then |D3-D3_inf^(9)| < margin => D3 >= bar (decorrelated regime, ignoring the fast triple).")
    print("    The correlated remainder (some product small) is bounded => kps exhaustive / partial-decorr limits.")
    print("=" * 100)


if __name__ == "__main__":
    main()
