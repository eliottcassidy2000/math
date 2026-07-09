"""
lrc14_hembed_ruler_embedding_opus_S176.py   (opus-2026-07-09-S176)

Ground THM-527 Part A / klein-S203's hembed (the SHARED blocker of the good-period AND density routes):
  a good period j (cluster teeth {frac(e_i j/Vmax)} have circular maxgap > 1/7, phase phi clears them:
  nearInt(phi - tooth_i) > 1/14) should be realizable at a REAL time tau with min_i nearInt(v_i tau) >= 1/14,
  v_i = Vmax - e_i the actual speeds.

TWO-SCALE structure: slow x = tau, fast phi = frac(Vmax tau).  Naive embed tau = (j + phi)/Vmax has a
COUPLING error: frac(e_i tau) = tooth_i + e_i*phi/Vmax (mod 1), shift ~ spread/Vmax.  This script tests:
  (1) the naive embed tau=(j+phi)/Vmax: min_i nearInt(v_i tau) vs 1/14 -- does it clear, and margin?
  (2) the coupling error size (spread/Vmax) vs the gap slack (maxgap - 1/7);
  (3) the TRUE M(S)=max_tau min_i nearInt(v_i tau) >= 1/14 (LRC holds for the set) + the witness tau,
      and whether a BETTER tau (search near the ruler period) realizes the good period cleanly.
"""
import numpy as np

INV14 = 1.0 / 14


def nearInt(x):
    f = x - np.floor(x)
    return np.minimum(f, 1 - f)


def circ_maxgap(pts):
    p = np.sort(pts % 1.0)
    g = np.diff(np.concatenate([p, [p[0] + 1.0]]))
    return g.max()


def good_periods(e, Vmax):
    """j in [1,Vmax) with maxgap{frac(e_i j/Vmax)} > 1/7 (cluster good period). Returns list of (j, gap, phi*)."""
    out = []
    for j in range(1, Vmax):
        teeth = (np.array(e) * j % Vmax) / Vmax
        p = np.sort(teeth)
        gaps = np.diff(np.concatenate([p, [p[0] + 1.0]]))
        k = gaps.argmax()
        if gaps[k] > 1.0 / 7 + 1e-12:
            phistar = (p[k] + gaps[k] / 2) % 1.0   # gap center = best clearing phase
            out.append((j, gaps[k], phistar))
    return out


def M_of_set(v, N=200003):
    """true M(S) = max_tau min_i nearInt(v_i tau), tau in [0,1)."""
    tau = (np.arange(N) + 0.5) / N
    V = np.array(v)
    reach = np.min(nearInt(np.outer(tau, V)), axis=1)   # N x 13 -> N
    k = reach.argmax()
    return reach[k], tau[k]


print("=" * 100)
print("hembed / THM-527 Part A: does a good period embed into a lonely real time?  1/14 = 0.07143")
print("  setup: e = co-offset cluster (incl 0), Vmax >= max e, speeds v_i = Vmax - e_i")
print("=" * 100)
# e = dissociated co-offset cluster (has good periods); vary Vmax from >>spread to ~spread
E_CLUSTER = [0, 1, 3, 7, 8, 12, 15, 17, 20, 22, 25, 27, 30]   # 13 dissociated co-offsets, spread 30
tests = {
    "Vmax=1000 >> spread 30 (clean)": (E_CLUSTER, 1000),
    "Vmax=200  >  spread 30        ": (E_CLUSTER, 200),
    "Vmax=60   ~2x spread 30       ": (E_CLUSTER, 60),
    "Vmax=31   ~ spread 30 (hard)  ": (E_CLUSTER, 31),
}
for name, (elist, Vmax) in tests.items():
    e = sorted(set(elist))
    if len(e) != 13:
        print(f"  {name}: |e|={len(e)} skip"); continue
    v = [Vmax - ei for ei in e]                 # speeds (Vmax-runner has e=0)
    gps = good_periods(e, Vmax)
    Mtrue, tauM = M_of_set(v)
    print(f"\n  {name}:  Vmax={Vmax}, spread(e)={max(e)-min(e)}, #good-periods={len(gps)}")
    print(f"    TRUE M(S) = max_tau min nearInt(v_i tau) = {Mtrue:.5f}  ({'>=1/14 LRC-OK' if Mtrue>=INV14-1e-6 else '<1/14!'}) at tau={tauM:.5f}")
    if gps:
        # test naive embed on the widest good period
        j, gap, phi = max(gps, key=lambda t: t[1])
        tau_embed = (j + phi) / Vmax
        reach_embed = np.min(nearInt(np.array(v) * tau_embed))
        coupling = (max(e)) * phi / Vmax
        slack = gap - 1.0 / 7
        print(f"    widest good period j={j}, gap={gap:.4f} (slack over 1/7 = {slack:.4f}), phi*={phi:.4f}")
        print(f"    NAIVE embed tau=(j+phi)/Vmax={tau_embed:.5f}: min nearInt(v_i tau) = {reach_embed:.5f} "
              f"({'clears 1/14' if reach_embed>=INV14-1e-9 else 'FAILS'})")
        print(f"    coupling error e_max*phi/Vmax = {coupling:.4f}  vs gap slack/2 = {slack/2:.4f}  "
              f"({'slack absorbs error' if slack/2 > coupling else 'ERROR EXCEEDS SLACK'})")
        # search tau near j/Vmax (fine) for the best realization
        taus = j / Vmax + (np.arange(20000) / 20000) / Vmax   # one ruler cell
        reach = np.min(nearInt(np.outer(taus, np.array(v))), axis=1)
        kk = reach.argmax()
        print(f"    BEST tau in ruler cell [j/Vmax, (j+1)/Vmax): min nearInt = {reach[kk]:.5f} at tau={taus[kk]:.6f} "
              f"({'>=1/14' if reach[kk]>=INV14-1e-9 else '<1/14'})")
print()
print("  READING: if the naive embed clears 1/14 when slack/2 > coupling error (Vmax >> spread), the")
print("  embedding is clean there; if Vmax~spread the coupling error may exceed slack (the hard case).")
print("  BEST-tau-in-cell tells whether SOME tau in the ruler cell realizes loneliness (existence of embed).")
print("=" * 100)
