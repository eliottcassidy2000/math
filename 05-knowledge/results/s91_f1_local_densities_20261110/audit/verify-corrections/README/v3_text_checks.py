#!/usr/bin/env python3
"""Verifier of the corrections to the S91 README: mechanical checks of the text.

(1) phrases excluded by the session author's decisions D7, D9 and by the findings do not occur;
(2) every receipt path written in backticks exists (audit/... next to the README; bridges/, rooted_spectrum/,
    rates_exact/, blind_proof/ in the scripts directory; 04-computation/... from the repository root);
(3) every THM-nnnn / HYP-nnnn cited exists in canon (MISTAKE numbers are listed, not required: the ledger entry of
    the audit is written at close-out);
(4) numbers of the corrected passages are found in the receipts named for them;
(5) line width (tables, code lines and the title excepted) and ASCII.
"""
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RDIR = os.path.normpath(os.path.join(HERE, "..", "..", ".."))                 # the README's directory
ROOT = os.path.normpath(os.path.join(RDIR, "..", "..", ".."))                  # repository root
SDIR = os.path.join(ROOT, "04-computation", "experiments", "s91_f1_local_densities_20261110")
README = os.path.join(RDIR, "README.md")
text = open(README, encoding="utf-8").read()
lines = text.split("\n")
RES = []


def check(name, cond, detail=""):
    RES.append(bool(cond))
    print(("PASS  " if cond else "FAIL  ") + name + (("   [" + str(detail) + "]") if detail != "" else ""), flush=True)


# (1) excluded phrases
flat = re.sub(r"\s+", " ", text)
bad = ["level of distribution", "major-arc", "major arc", "local side is closed", "cannot work",
       "the reason single-base", "relative variance about", "exactly the dichotomy", "feed the first-order",
       "three of the four inputs", "PRE-AUDIT", "Not independently audited", "Windows fail",
       "Local structure belongs", "boosted class", "a third of the balanced", "rank has period 6 in",
       "partial attribution", "followed by sparse ones (`13`"]
hits = [b for b in bad if b in flat and not (b == "rank has period 6 in" and "corank has period 6 in" in flat
                                              and flat.count("rank has period 6 in") == flat.count("corank has period 6 in"))]
check("1  excluded phrases absent", not hits, hits)
ss = [m.start() for m in re.finditer("singular series", flat)]
ctx = [flat[max(0, i - 60): i + 40] for i in ss]
ok_ss = all(("not the singular series" in c) or ("singular series proper" in c) or ("the singular series of a sheet" in c)
            for c in ctx)
check("1b 'singular series' occurs only as 'not the singular series', 'the singular series proper would ...' or "
      "'would make ... the singular series of a sheet' (%d occurrences)" % len(ss), ok_ss, ctx if not ok_ss else "")

# (2) receipt paths
paths = set(re.findall(r"`([A-Za-z0-9_./<>*-]+\.(?:out|txt|py|md|npz|json))`", text))
missing = []
for p in sorted(paths):
    if "<" in p or "*" in p or p.startswith("OUT"):
        continue
    if p.startswith("audit/"):
        cand = [os.path.join(RDIR, p)]
    elif p.startswith("04-computation/") or p.startswith("05-knowledge/") or p.startswith("01-canon/"):
        cand = [os.path.join(ROOT, p)]
    elif p.split("/")[0] in ("bridges", "rooted_spectrum", "rates_exact", "blind_proof"):
        cand = [os.path.join(SDIR, p)]
    else:
        # bare file names: scripts directory, its sub-directories, or the audit directory quoted just before
        cand = [os.path.join(SDIR, p)] + glob.glob(os.path.join(SDIR, "*", p)) + glob.glob(os.path.join(RDIR, "audit", "*", p))
    if not any(os.path.exists(c) for c in cand):
        missing.append(p)
check("2  the %d receipt paths in backticks exist" % len(paths), not missing, missing)
bare = ["twisted_spectrum.py", "twisted_spectrum_2adic.py", "mod4_symbol.py", "box_identity.py", "divisible_states.py",
        "constants.py", "first_order_exact.py", "ap_constants.py", "regions_check.py", "window_split.py",
        "first_order_exact_28.out", "ap_constants.out", "regions_check.out", "prereg_window_27_24.json",
        "rooted_spectrum/census.py", "rooted_spectrum/source_supply.py", "rooted_spectrum/source_supply_levels78.py",
        "rooted_spectrum/line_density_classes.py", "rooted_spectrum/head_parity_check.py",
        "rooted_spectrum/head_values_check.py", "rooted_spectrum/source_witnesses.py",
        "rooted_spectrum/line_density_order10.py", "rooted_spectrum/line_density_order9.out",
        "rooted_spectrum/line_density_order10.out", "rooted_spectrum/head_values_check.out",
        "rates_exact/scan_003_101.npz", "rates_exact/scan_103_521.npz", "rates_exact/scan_103_521.out",
        "rates_exact/scan_1023.npz", "rates_exact/scan_1143.npz", "PREREGISTRATION.md"]
mb = [b for b in bare if not os.path.exists(os.path.join(SDIR, b))]
check("2b the scripts and outputs of section 7 exist", not mb, mb)
notes = ["05-knowledge/results/odd_c_seed_symmetry_20261110/EXTENSIONS.md",
         "05-knowledge/results/odd_c_seed_symmetry_20261110/PRIMES_AND_PAPERS.md",
         "05-knowledge/results/kps_s22_sign_flip_towers_20261110/README.md",
         "05-knowledge/results/kps_s22_sign_flip_towers_20261110/SOURCES_AND_DOMINATION.md",
         "05-knowledge/results/session257_permutohedron_squares_20261111/history/sources/INCOMING_EXTENSION_CONTENT.md",
         "05-knowledge/results/spectrum_analogies_nearjoin_20261102/README.md",
         "05-knowledge/results/seam_gadget_thickness_spectrum_20261104.md",
         "05-knowledge/results/seven_twentyone_mersenne_openai_math_20261006.md",
         "05-knowledge/results/one_apex_join_representation_mass_pi2_over_12_20261106.md",
         "05-knowledge/results/old_mark_bound_lorentzian_gadgets_one_apex_join_coverage_20261104.md",
         "05-knowledge/results/session193_triangulated_proof_20261105/modular_sync/PROOF.md",
         "05-knowledge/results/session192_surprising_approaches_20261105/dyadic/PROOF.md",
         "05-knowledge/results/session25_spectrum_wiring_20261021/SYNTHESIS.md",
         "05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md",
         "05-knowledge/results/session223_root_changes_and_local_factors_20261107",
         "05-knowledge/results/session261_creative_chanpak_20261111",
         "05-knowledge/results/session262_interface_supply_20261111",
         "session206/exact-modular-image/DYADIC_SURJECTIVITY.md",
         "session206/weighted-parameter-capacity/WEIGHTED_CAPACITY.md",
         "04-computation/experiments/chanpak_mass_exploratory_20261106/local_factor_mod3_prediction.py",
         "04-computation/experiments/one_apex_join_mass_20261106_audit_S82/aud82_localfactor_28.out"]
mn = [n for n in notes if not os.path.exists(os.path.join(ROOT, n))]
check("2c the notes of other sessions that the README credits exist (%d files)" % len(notes), not mn, mn)

# (3) canon ids
thm = sorted(set(re.findall(r"THM-(\d+)", text)), key=int)
hyp = sorted(set(re.findall(r"HYP-(\d+)", text)), key=int)
mt = [t for t in thm if not glob.glob(os.path.join(ROOT, "01-canon", "theorems", "THM-%s-*.md" % t))
      and not glob.glob(os.path.join(ROOT, "01-canon", "theorems", "THM-%s-*.md" % t.zfill(3)))]
mh = [h for h in hyp if not glob.glob(os.path.join(ROOT, "05-knowledge", "hypotheses", "HYP-%s-*.md" % h))]
check("3  the %d theorems and %d hypotheses cited exist" % (len(thm), len(hyp)), not mt and not mh, (mt, mh))
led = open(os.path.join(ROOT, "01-canon", "MISTAKES.md"), encoding="utf-8").read()
mk = sorted(set(re.findall(r"MISTAKE-(\d+)", text)), key=int)
print("      MISTAKE numbers cited: %s; present in 01-canon/MISTAKES.md now: %s"
      % (mk, [m for m in mk if ("## MISTAKE-%s " % m) in led]))


# (4) numbers against receipts
def rd(*parts):
    return open(os.path.join(*parts), encoding="utf-8", errors="replace").read()


a4 = rd(RDIR, "audit", "proof-5128-CF", "a4_moduli_26.out")
blk = a4[a4.index("m = 9 = 3^2"): a4.index("m = 25 = 5^2")]
s1 = blk[blk.index("sheet 1"):]
row26 = [ln for ln in s1.split("\n") if ln.strip().startswith("26 |")][0]
nums = re.findall(r"\d\.\d{4}", row26)
check("4a modulo 9, sheet 1, X = 2^26: measured 0.848, 1.229, 1.456; first-order 0.852, 1.216, 1.457; limits "
      "0.889, 1.185, 1.296", ["%.3f" % float(v) for v in nums[:6]] == ["0.848", "0.852", "1.229", "1.216", "1.456", "1.457"]
      and "v=0: 0.8889  v=1: 1.1852  v>=2: 1.2963" in s1, nums[:6])
ap = rd(RDIR, "audit", "recompute-5128", "a5_Ap.out")
check("4b A_13 = -0.0038, -0.0049 and A_19 = -0.00014, -0.00032 (means of the last six partial sums)",
      "mean of last 6: -0.00377, -0.00488" in ap and "mean of last 6: -0.00014, -0.00032" in ap)
dv = rd(RDIR, "audit", "recompute-5128", "a4_divisible.out")
check("4c divisible states: 13: 32, 721, 877; 19: 218, 256, 294, 313, 351, 959",
      "13  [32, 721, 877, 1111, 1462, 1618]" in dv and "19  [218, 256, 294, 313, 351, 959]" in dv)
ms = rd(RDIR, "audit", "recompute-5128", "a7_misses.out")
check("4d dispersion and misses: Var/mean^2 0.5639 / 0.5066 and 0.4459 / 0.3900; Var/mean 2.76-2.79 and 3.06-3.11; "
      "miss rates 1.381, 0.823, 0.615, 0.347 %; ratio 2.52 to 3.69; difference +0.67",
      all(t in ms for t in ("sheet 0 0.5639, sheet 1 0.5066", "sheet 0 0.4459, sheet 1 0.3900", "= 1.381 %", "= 0.823 %",
                            "= 0.615 %", "= 0.347 %", "= 2.52", "= 3.69", "difference +0.6655", "sheet 0 2.760, sheet 1 3.108",
                            "sheet 0 2.782, sheet 1 3.082", "F1 misses: 65")))
pf = rd(RDIR, "audit", "prior-art-wording", "a1_prime_power_factor.out")
check("4e p = 3, sheet 1, [2^22, 2^24): 1.2228 for 3 || N, 1.4497 for 9 | N", "p | . exactly 1.2228; p^2 | . 1.4497" in pf)
b3 = rd(SDIR, "bridges", "b3_f1_data.out")
check("4f 1039: rate 0.763636; S_1 = 0.81135, S_0 = 1.05523, S = 0.91667; E = -0.00121, +0.00189, predicted -0.00025, "
      "-0.00031, se 0.00364; 27: +0.49121, +0.48898; 58 divisible states, 62.2 expected",
      all(t in b3 for t in ("max_k rho(T_k)/2 = 0.763636", "N = 1039: S_1 = 0.81135  S_0 = 1.05523  S = 0.91667",
                            "E measured = -0.00121 (s = 0), +0.00189 (s >= 1)   predicted -0.00025, -0.00031", "se = 0.00364",
                            "p =   27: E measured = +0.49121 (s = 0), +0.48898", "count 58, expected 62.2")))
rt = {}
for ln in rd(SDIR, "rates_exact", "rates_by_modulus.txt").split("\n"):
    if ln and not ln.startswith("#"):
        f = ln.split()
        rt[int(f[0])] = float(f[1])
want = {3: "0.8660", 5: "0.7418", 7: "0.9320", 11: "0.8587", 13: "0.8272", 17: "0.8289", 19: "0.8053", 29: "0.7894",
        31: "0.9002", 43: "0.8342", 257: "0.8551"}
check("4g rates of the section 1 table; 262 moduli = every odd m <= 521 and 1023, 1143; maximum 0.932039 at the "
      "multiples of 7 only; the list 93, 31, 381, 15, 219, 171, 127, 51, 3, 73; 5 the fastest",
      all("%.4f" % rt[m] == s for m, s in want.items()) and len(rt) == 262
      and sorted(rt) == list(range(3, 522, 2)) + [1023, 1143] and max(rt.values()) == 0.932039
      and all((v == 0.932039) == (m % 7 == 0) for m, v in rt.items())
      and ["%.3f" % rt[m] for m in (93, 31, 381, 15, 219, 171, 127, 51, 3, 73)]
      == ["0.914", "0.900", "0.900", "0.892", "0.879", "0.878", "0.875", "0.869", "0.866", "0.852"]
      and min(rt, key=rt.get) == 5 and "%.3f" % rt[5] == "0.742")
top = sorted({v for m, v in rt.items() if m % 7}, reverse=True)[:9]
check("4h the nine largest rates at moduli prime to 7 are those of 93, 31, 381, 15, 219, 171, 127, 51, 3",
      top == [rt[m] for m in (93, 31, 381, 15, 219, 171, 127, 51, 3)], top)
ss_ = rd(SDIR, "rooted_spectrum", "source_supply.out")
check("4i source_supply.out (6): correlations 0.68 (26 offsets), 0.66, 0.41 on 3 <= h <= 73 (34 offsets); poorest "
      "63, 47, 67 (six classes), 35, 39",
      all(t in ss_ for t in ("23 <= h <= 73 (26 attained offsets): 0.68", "15 <= h <= 73 (29 attained offsets): 0.66",
                             "3 <= h <= 73 (34 attained offsets): 0.41",
                             "(0.156, 63, 2), (0.208, 47, 2), (0.265, 67, 6), (0.278, 35, 2), (0.350, 39, 2)")))
lc = rd(SDIR, "rooted_spectrum", "line_density_classes.out")
check("4j line_density_classes.out: 0.5784 against 0.5574 at 3 (ratio 1.038); 5, 7, 11, 13 within 0.004 of 1",
      " 3 | 0.5784 | 0.5574" in lc
      and all(abs(float(re.search(r"^ ?%d \| [^|]*\| [^|]*\| (\d\.\d{3}) \|" % q, lc, re.M).group(1)) - 1) <= 0.0041
              for q in (5, 7, 11, 13))
      and re.search(r"^ ?3 \| [^|]*\| [^|]*\| (\d\.\d{3}) \|", lc, re.M).group(1) == "1.038")
pr = rd(RDIR, "audit", "recompute-5128", "a10_prereg_full.out")
check("4k all 279 items inside at 4, 3, 2 measured standard errors; 43691 and 241 intervals contain 0",
      pr.count("total items: 279") == 4 and "outside -> REFUTED" not in pr and pr.count("(contains 0)") == 4)
lt = rd(RDIR, "audit", "recompute-5128", "a13_lightest.out")
check("4l lightest integers: share 0.0873; 0.7273, 0.7244, 0.7319 of the mean; 22.2, 25.2, 28.0 % of the misses; "
      "word pairs 0.8016, 0.8097, 0.8087", all(t in lt for t in ("0.0873", "= 0.7273", "= 0.7244", "= 0.7319", "(22.2 %",
                                                              "(25.2 %", "(28.0 %", "= 0.8016", "= 0.8097", "= 0.8087")))

# (5) form
wide = [i + 1 for i, ln in enumerate(lines) if len(ln) > 118 and not ln.startswith("|") and not ln.startswith("    ")
        and not ln.startswith("# ")]
check("5  no text line wider than 118 characters; ASCII only", not wide and all(ord(ch) < 128 for ch in text), wide)

print("ALL CHECKS PASSED (%d checks)" % len(RES) if all(RES) else "SOME CHECKS FAILED")
sys.exit(0 if all(RES) else 1)
