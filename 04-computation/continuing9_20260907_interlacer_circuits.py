"""Exact C-only/full-C,D circuit sidecar at S=13,Q=59.

No inference of an excluded full-model circuit word is made. All27 C-only
words are explicit rational witnesses; the neighborhood theorem is continuity
at a strictly positive Hankel center. No producer module is imported.
"""
from pathlib import Path
from hashlib import sha256
from itertools import product
from math import comb
import json
import sys
import sympy as sp

sys.stdout.reconfigure(newline="\n")
HERE = Path(__file__).resolve().parent
DEST = HERE.parent/"05-knowledge/results" if HERE.name == "04-computation" else HERE
STEM = Path(__file__).stem
t = sp.Symbol("t")
gates = 0


def need(test, message):
    global gates
    gates += 1
    if not bool(test):
        raise RuntimeError(message)


def moments(x, y, z, which):
    den = [1, -13, 55, -x, y, -z]
    num = [1, -12, 45, -2*x/3, 3*y/7] if which == "C" else [1, -11, 36, -5*x/12, y/7]
    out = []
    for k in range(9):
        out.append((num[k] if k < len(num) else 0) -
                   sum(den[j]*out[k-j] for j in range(1, min(k, 5)+1)))
    return out


def data(x, y, z):
    x, y, z = map(sp.Rational, (x, y, z))
    need(min(x, y, z) > 0, "positive coefficient coordinates")
    B = sp.Poly(t**5-13*t**4+55*t**3-x*t**2+y*t-z, t)
    need(B.count_roots(0, sp.oo) == 5 and sp.degree(sp.gcd(B, B.diff())) == 0,
         "literal positive simple beta roots")
    ratios = [sp.Rational(831875, 8788)/x, 13*x**3/(166375*y), 44*y**3/(x**3*z)]
    row = {"xyz": [x, y, z], "C": ratios, "word": [int(sp.sign(c-1)) for c in ratios]}
    for A in ("C", "D"):
        mm = moments(x, y, z, A)
        minors = [sp.Matrix([[mm[i+j] for j in range(k)] for i in range(k)]).det()
                  for k in range(1, 6)]
        row[A+"_moments"] = mm
        row[A+"_leading_minors"] = minors
    return row


q = sp.Rational(275, 338)
mu = sp.Rational(13, 5)


def chart(cc):
    out = []
    for k in range(6):
        value = sp.binomial(5, k)*mu**k*q**comb(k, 2)
        for j in range(2, k):
            value /= cc[j-2]**comb(k-j+1, 2)
        out.append(sp.cancel(value))
    need(out[:3] == [1, 13, 55], "exact fixed S13 Q59")
    return out[3:]


center = data(*chart([1, 1, 1]))
need(center["word"] == [0, 0, 0], "all-one circuit center")
need(all(x > 0 for x in center["C_leading_minors"]), "strict C-only center")
need(all(x > 0 for x in center["D_leading_minors"][:3]), "D ordinary degree4 packet survives")
need(center["D_leading_minors"][3] < 0, "D ordinary degree6 packet rejects center")
need(center["D_leading_minors"][4] < 0, "full D packet rejects center")

main_cert = HERE/"continuing9_20260907_fixed_moment_circuits_certificate.json"
if not main_cert.exists():
    main_cert = DEST/main_cert.name
raw = main_cert.read_bytes()
need(sha256(raw).hexdigest() == "291b62a5638f8a057ccf514851cb0b5d100a2a61a3da8eff9151b49040be11b9",
     "frozen fullmoment coefficient certificate")
main = json.loads(raw)
reference = {tuple(rec["word"]): rec for rec in main["cases"][2]["targets"]}
bank = []
for word in product((-1, 0, 1), repeat=3):
    cc = [1+sp.Rational(w, 2048) for w in word]
    xyz = chart(cc)
    inherited = [sp.Rational(ee)*mu**k for k, ee in enumerate(reference[word]["coefficients"])]
    need(xyz == inherited[3:], "same literal main target")
    record = data(*xyz)
    need(record["word"] == list(word) and record["C"] == cc, "exact requested word")
    need(all(x > 0 for x in record["C_leading_minors"]), "explicit strict C-only witness")
    need(record["D_leading_minors"][3] < 0, "each displayed C-only witness fails D")
    bank.append(record)

controls = [("+++", ("77.454", "8.902", ".02558")),
            ("++-", ("77.3613", "8.6001", ".17694")),
            ("+-+", ("86.2333", "51.3919", "8.6469")),
            ("-++", ("97.7028", "70.6021", "14.5020"))]
full = []
for label, xyz in controls:
    record = data(*map(sp.Rational, xyz))
    expected = [1 if c == "+" else -1 for c in label]
    need(record["word"] == expected, "full-model exact circuit word")
    need(all(x > 0 for x in record["C_leading_minors"]+record["D_leading_minors"]),
         "strict simultaneous full-model point")
    record["label"] = label
    full.append(record)


def encode(value):
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    if isinstance(value, list):
        return [encode(v) for v in value]
    if isinstance(value, sp.Basic):
        return str(value)
    return value


certificate = {"status": "Exact positive C-only and full-model examples; omitted full words OPEN",
               "center": center, "C_only_27": bank, "full_model_four": full,
               "main_certificate_sha256": sha256(raw).hexdigest(), "gates": gates}
payload = (json.dumps(encode(certificate), sort_keys=True, indent=2)+"\n").encode()
(DEST/(STEM+"_certificate.json")).write_bytes(payload)
print("CENTER", center["xyz"])
print("CENTER C leading minors:+++++; D leading minors:+++--; circuits000")
print("C_ONLY: all27 exact radius1/2048 target words have positive simple beta roots and strict C interlacing")
print("D_BOUNDARY: all27 displayed C-only targets fail the ordinary degree6 D Hankel")
print("FULL_MODEL: four exact strict words +++,++-,+-+,-++ with both Hankels positive definite")
print("OPEN: no excluded full-model word follows from this finite positive bank")
print("CERTIFICATE_SHA256", sha256(payload).hexdigest())
print("PASS", gates, "always-active exact gates; raw LF")
