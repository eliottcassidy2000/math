#!/usr/bin/env python3
"""audit_mcraw.py -- independent audit of section 7 and of the Outcomes table.

(1) Recompute the summary numbers of Run A from the author's raw files ../out/mc11_*.out (and of the order-10 control from
    ../out/mc10_*.out): N, D, standard error, z, per-run values, rho, E|b|, rank fractions; apply the registered rules R1-R4.
(2) Recompute the Run C verdicts (R7, R8, control) from ../out/c35_*.out.
(3) My own estimator (audit_cube.c, mode mc; new seeds): out/mc_audit_11_*.out, out/mc_audit_10_*.out; consistency with
    the author's estimate at order 12 and with the exact D_11.
(4) The reference rank law quoted in R2 (random alternating 11 x 11 matrix over F_2 with an odd number of ones above the
    diagonal): exact, from the eigenvalues of the association scheme of alternating forms (Delsarte-Goethals), cross-checked
    by brute force for 4 x 4 and 5 x 5 matrices.
"""
import glob, os, re, math, itertools
from fractions import Fraction as F

here = os.path.dirname(os.path.abspath(__file__))
pk = os.path.dirname(here)

def parse_author(fn):
    t = open(fn).read().split()
    d = {"m": int(t[1]), "N": int(t[3])}
    i = t.index("odd,");
    j = t.index("edges)"); d["fail"] = [int(x) for x in t[j + 1:j + 5]]
    for key in ("S1", "S2", "SA", "SU"):
        d[key] = float(t[t.index(key) + 1])
    d["npos"] = int(t[t.index("npos") + 1])
    d["rank"] = {}
    for r in range(0, d["m"] + 1):
        k = t.index(f"r{r}"); d["rank"][r] = (int(t[k + 1]), int(t[k + 2]))
    return d

def summarize(files, label, exact=None):
    runs = [parse_author(f) for f in files]
    m = runs[0]["m"]; N = sum(r["N"] for r in runs)
    S1 = sum(r["S1"] for r in runs); S2 = sum(r["S2"] for r in runs); SA = sum(r["SA"] for r in runs); SU = sum(r["SU"] for r in runs)
    mean = S1 / N; var = S2 / N - mean ** 2; se = math.sqrt(var / N); unit = 4 ** (m + 1)
    print(f"{label}: m = {m}, runs {len(runs)}, N = {N}, control failures {[sum(r['fail'][i] for r in runs) for i in range(4)]}")
    print(f"  D_{m+1} = {mean:.6e} +- {se:.2e} ; times 4^{m+1}: {mean*unit:.1f} +- {se*unit:.1f} ; z = {mean/se:.2f}")
    per = [r["S1"] / r["N"] * unit for r in runs]
    print("  per-run values (times 4^n):", ", ".join(f"{x:.0f}" for x in per))
    se_run = math.sqrt(var / runs[0]["N"]) * unit
    chi2 = sum(((x - mean * unit) / se_run) ** 2 for x in per)
    print(f"  per-run standard error {se_run:.0f}; chi^2 of the runs about the pooled mean = {chi2:.2f} on {len(runs)-1} d.o.f.; "
          f"empirical s.e. of the mean from the run scatter = {math.sqrt(sum((x - mean*unit)**2 for x in per)/(len(per)-1)/len(per)):.1f}")
    rho = 2 ** m * S2 / N
    rhos = [2 ** m * r["S2"] / r["N"] for r in runs]
    se_rho = math.sqrt(sum((x - rho) ** 2 for x in rhos) / (len(rhos) - 1) / len(rhos))
    print(f"  rho_{m} = 2^m E[b^2] = {rho:.5f} (per run: {', '.join(f'{x:.5f}' for x in rhos)}; s.e. from the run scatter {se_rho:.5f})")
    print(f"  E|b_R| = {SA/N:.6f} ; E[b_R] (unsigned vertex-deletion correlation) = {SU/N:.6e} ; forced part 2^(1-m) = {2.0**(1-m):.6e}")
    tot = {r: [sum(run["rank"][r][0] for run in runs), sum(run["rank"][r][1] for run in runs)] for r in range(m + 1)}
    for r in range(m + 1):
        c = tot[r][0] + tot[r][1]
        if c:
            fr = c / N
            nz = tot[r][1] / c
            print(f"   rank {r:2d}: fraction {fr:.6f} +- {math.sqrt(fr*(1-fr)/N):.6f} ; non-zero-bias frequency {nz:.5f} +- {math.sqrt(nz*(1-nz)/c):.5f}")
    if exact is not None:
        print(f"  exact D_{m+1} 4^{m+1} = {float(exact*unit):.1f}; deviation of the estimate {((mean - float(exact))/se):+.2f} s.e.")
    return dict(mean=mean, se=se, unit=unit, rho=rho, se_rho=se_rho, tot=tot, N=N, SA=SA)

print("(1) Run A from the author's raw files")
A = summarize(sorted(glob.glob(os.path.join(pk, "out", "mc11_1*.out"))), "Run A")
est = A["mean"] * A["unit"]; se = A["se"] * A["unit"]
if 1500 <= est <= 2500 and A["mean"] / A["se"] > 3: v1 = "CONFIRMED"
elif (est < 1500 and (1500 - est) > 2 * se) or (est > 2500 and (est - 2500) > 2 * se): v1 = "FAILED"
else: v1 = "UNDECIDED"
print(f"  R1 rule: point estimate {est:.1f}, interval [1500, 2500], distance to the interval {(1500-est)/se:.2f} s.e. -> {v1}")
N = A["N"]; f10 = sum(A["tot"][10]) / N; f6 = sum(sum(A["tot"][r]) for r in range(0, 7)) / N
print(f"  R2 rule: rank-10 fraction {f10:.6f}, |diff from 0.83898| = {abs(f10-0.83898):.6f} (<= 0.006: {abs(f10-0.83898) <= 0.006}); rank <= 6 fraction {f6:.6f} in [0.0012, 0.0030]: {0.0012 <= f6 <= 0.0030}")
print(f"  R3 rule: rho_11 = {A['rho']:.5f} in [1.015, 1.045]: {1.015 <= A['rho'] <= 1.045} (distance to the lower end {(A['rho']-1.015)/A['se_rho']:.1f} run-scatter s.e.)")
nz10 = A["tot"][10][1] / sum(A["tot"][10])
print(f"  R4 rule: corank-1 non-zero-bias frequency {nz10:.5f} in [0.500, 0.512]: {0.500 <= nz10 <= 0.512}")
print("\n    order-10 control from the author's raw files")
summarize(sorted(glob.glob(os.path.join(pk, "out", "mc10_2*.out"))), "control", exact=F(260546929, 2 ** 39))

print("\n(2) Run C from the author's raw files")
rows = {}
for fn in sorted(glob.glob(os.path.join(pk, "out", "c35_1*.out"))):
    t = open(fn).read().split()
    if len(t) == 5: rows[int(t[0])] = (int(t[1]), int(t[2]), int(t[3]), int(t[4]))
def kron2(n): return 1 if n % 8 in (1, 7) else -1
def closed(n):
    if n % 2 == 1:
        e3 = F(1, 2 ** ((n - 1) // 2))
        return e3, e3 * (2 + kron2(n) * F(1, 2 ** ((n - 1) // 2)) - F(1, 2 ** (n - 2))), F(1, 2 ** (n - 1)) * (3 - kron2(n) * F(1, 2 ** ((n - 3) // 2)))
    if n % 4 == 2: return F(1, 2 ** ((n - 2) // 2)), F(1, 2 ** ((n - 2) // 2)), F(1, 2 ** (n - 2))
    return (F(1, 2 ** ((n - 2) // 2)),) + closed(n - 1)[1:]
zmax = 0; zmax3 = 0
for n, (Nn, s3, s5, s35) in sorted(rows.items()):
    se_ = 1 / math.sqrt(Nn); e = closed(n)
    z3 = (s3 / Nn - float(e[0])) / se_; z5 = (s5 / Nn - float(e[2])) / se_; z35 = (s35 / Nn - float(e[1])) / se_
    zmax = max(zmax, abs(z5), abs(z35)); zmax3 = max(zmax3, abs(z3))
    print(f"  n={n}: N={Nn} E[{{3}}]={s3/Nn:.6f} (z={z3:+.2f}) E[{{5}}]={s5/Nn:.7f} (z={z5:+.2f}) E[{{3,5}}]={s35/Nn:.7f} (z={z35:+.2f})")
se_ = 1 / math.sqrt(rows[14][0])
print(f"  R7: n=14: ({rows[14][3]/rows[14][0]:.6f} - 2^-6)/s.e. = {(rows[14][3]/rows[14][0] - 2**-6)/se_:+.2f} ; n=18: {(rows[18][3]/rows[18][0] - 2**-8)/se_:+.2f}  (rule: within 3)")
print(f"  R8: (E_12 - E_11)/s.e.(diff) = {(rows[12][3]-rows[11][3])/rows[11][0]/(se_*math.sqrt(2)):+.2f} ; (E_16 - E_15) = {(rows[16][3]-rows[15][3])/rows[15][0]/(se_*math.sqrt(2)):+.2f}  (rule: within 3)")
print(f"  control: largest |z| of the triangle column {zmax3:.2f} (rule: within 3); largest |z| of the 16 pentagon comparisons {zmax:.2f}")

print("\n(3) My own estimator (audit_cube.c mc, direct evaluation of H mod 4 on all out-sets; new seeds)")
def mine(pattern, exact=None, compare=None):
    fl = sorted(glob.glob(os.path.join(here, "out", pattern)))
    if not fl:
        print("  (no output yet for", pattern, ")"); return
    S1 = S2 = SA = 0.0; N = 0; per = []; m = None; h8 = [0, 0, 0, 0]
    for fn in fl:
        t = open(fn).read().split()
        if len(t) < 10: continue
        m = int(t[1]); n_ = int(t[3]); N += n_
        s1 = float(t[t.index("S1") + 1]); S1 += s1; S2 += float(t[t.index("S2") + 1]); SA += float(t[t.index("SA") + 1])
        per.append(s1 / n_ * 4 ** (m + 1))
        k = t.index("h8"); h8 = [a + int(b) for a, b in zip(h8, t[k + 1:k + 5])]
    if not N: return
    unit = 4 ** (m + 1); mean = S1 / N; se = math.sqrt((S2 / N - mean ** 2) / N)
    print(f"  m = {m}: N = {N} in {len(per)} runs; D_{m+1} = {mean:.6e} +- {se:.2e}; times 4^{m+1}: {mean*unit:.1f} +- {se*unit:.1f}; per run: {', '.join(f'{x:.0f}' for x in per)}")
    print(f"     rho_{m} = {2**m*S2/N:.5f}; E|b_R| = {SA/N:.6f}; relative masses of H mod 8 (1,3,5,7): " + ", ".join(f"{4*x/sum(h8):.5f}" for x in h8))
    if exact is not None:
        print(f"     exact value times 4^{m+1} = {float(exact*unit):.1f}; deviation {(mean-float(exact))/se:+.2f} s.e.")
    if compare is not None:
        a, sa = compare
        z = (mean * unit - a) / math.sqrt((se * unit) ** 2 + sa ** 2)
        w1 = 1 / (se * unit) ** 2; w2 = 1 / sa ** 2
        print(f"     author's estimate {a:.1f} +- {sa:.1f}: difference {mean*unit - a:+.1f} = {z:+.2f} combined s.e.; pooled estimate {(mean*unit*w1 + a*w2)/(w1+w2):.1f} +- {1/math.sqrt(w1+w2):.1f}")
mine("mc_audit_11_*.out", compare=(est, se))
mine("mc_audit_10_*.out", exact=F(260546929, 2 ** 39))

print("\n(4) Reference rank law of R2: random alternating n x n matrix over F_2 with an odd number of ones above the diagonal")
def gauss(a, b, q):
    if b < 0 or b > a: return 0
    num = 1; den = 1
    for i in range(b):
        num *= q ** (a - i) - 1; den *= q ** (i + 1) - 1
    return num // den
def counts(n):
    """N_r = number of alternating matrices of rank 2r; S_r = sum of (-1)^(weight) over them (eigenvalue at the dual class of J+I)."""
    m = n // 2; b = 4; c = 2 ** (2 * m + 1) if n % 2 else 2 ** (2 * m - 1)
    # rank of J + I over F_2: n-1 (n odd) or n (n even) -> index i = m in both cases
    def P(k, i):
        return sum((-1) ** (k - j) * b ** ((k - j) * (k - j - 1) // 2) * gauss(m - j, m - k, b) * gauss(m - i, j, b) * c ** j for j in range(k + 1))
    return [P(k, 0) for k in range(m + 1)], [P(k, m) for k in range(m + 1)]
def brute(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    Nr = {}; Sr = {}
    for t in range(1 << len(pairs)):
        rows = [0] * n
        for bit, (i, j) in enumerate(pairs):
            if (t >> bit) & 1: rows[i] |= 1 << j; rows[j] |= 1 << i
        rk = 0; rr = rows[:]
        for col in range(n):
            p = next((i for i in range(rk, n) if (rr[i] >> col) & 1), None)
            if p is None: continue
            rr[rk], rr[p] = rr[p], rr[rk]
            for i in range(n):
                if i != rk and (rr[i] >> col) & 1: rr[i] ^= rr[rk]
            rk += 1
        Nr[rk // 2] = Nr.get(rk // 2, 0) + 1; Sr[rk // 2] = Sr.get(rk // 2, 0) + (-1) ** bin(t).count("1")
    return [Nr.get(k, 0) for k in range(n // 2 + 1)], [Sr.get(k, 0) for k in range(n // 2 + 1)]
for n in (4, 5, 6):
    print(f"  n={n}: formula {counts(n)} ; brute force {brute(n)} : {'OK' if counts(n) == brute(n) else 'MISMATCH'}")
Nr, Sr = counts(11)
tot_odd = (sum(Nr) - sum(Sr)) // 2
print("  n=11, odd weight: P(rank 10), P(rank 8), P(rank <= 6) =", ", ".join(f"{x:.5f}" for x in
      (float(F((Nr[5] - Sr[5]) // 2, tot_odd)), float(F((Nr[4] - Sr[4]) // 2, tot_odd)), float(F(sum((Nr[k] - Sr[k]) // 2 for k in range(4)), tot_odd)))),
      " (quoted in R2: 0.83898 / 0.15965 / 0.00137)")
print("  n=11, unconditional:", ", ".join(f"{float(F(Nr[k], sum(Nr))):.5f}" for k in (5, 4)), "(rank 10, rank 8)")
p10 = float(F((Nr[5] - Sr[5]) // 2, tot_odd)); p8 = float(F((Nr[4] - Sr[4]) // 2, tot_odd))
print(f"  Run A against the exact reference: rank 10: {f10:.6f} - {p10:.6f} = {f10-p10:+.6f} = {(f10-p10)/math.sqrt(p10*(1-p10)/N):+.2f} s.e.; "
      f"rank 8: {sum(A['tot'][8])/N:.6f} - {p8:.6f} = {sum(A['tot'][8])/N-p8:+.6f} = {(sum(A['tot'][8])/N-p8)/math.sqrt(p8*(1-p8)/N):+.2f} s.e.")
