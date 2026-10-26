# Audit of HYP-9391 (C): the O(q^2) closed formula for TT7(P_q) against the auditor's own TT7 counts.
# Inputs (all produced by the auditor's code):
#   T5 counts from tt7_count.c (TT7 = q(q-1)/2 * T5(D)),  Y sums from ysums.c,  modular coefficients from forms_dump.gp.
# Checks:
#   1. sympy: the HYP's printed formula = A7/(q(q-1)) + B7 S_W4 + C7 S_(K6-P6) + g(q) with the audit's reduced form
#      (audit_tt7_explicit.out), S_W4 = -4 - a_q, S_(K6-P6) = |b|^2 - q^2 - 7, and the lane's g_formula.json.
#   2. 2^21 TT7(P_q) == q(q-1) [formula] exactly at every prime q = 3 mod 4 counted.
#   3. the auditor's TT7 counts against the lane's 71 counted values (paley_tt7_data*.out).
# Usage: tt7_formula_check.py T5file ysumsfile formsfile ; writes tt7_formula_check.out
import sys, re, json, os
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
T5f, Yf, Ff = sys.argv[1], sys.argv[2], sys.argv[3]
out = []
q, a, c, d, b2, Y02, Y05, Y09, Y18 = sp.symbols('q a c d b2 Y02 Y05 Y09 Y18')
hyp = (q**5 - 90*q**4 + 3116*q**3 - 49252*q**2 + 322374*q - 673743 - 2*(14*q**2 - 701*q + 7254)*a
       + 46*c - 83*d - 64*(q - 31)*b2 - 64*Y02 - 48*Y05 - 64*Y09 + 96*Y18)
A7 = q*(q-1)*(q-3)*(q**4 - 87*q**3 + 2785*q**2 - 38827*q + 198842)
B7 = 28*q**2 - 1288*q + 13238
C7 = 1912 - 64*q
gj = json.load(open(os.path.join(HERE, '..', 'g_formula.json')))
names = {'1': 1, 'q': q, 'q^2': q**2, 'q^3': q**3, 'a84': a, 'q a84': q*a, 'a86': c, 'a46': d, '|b|^2': b2,
         'Y02': Y02, 'Y05': Y05, 'Y09': Y09, 'Y18': Y18}
g = sum(sp.Rational(v) * names[k] for k, v in gj.items())
recon = sp.expand(sp.cancel(A7 / (q*(q-1))) + B7*(-4 - a) + C7*(b2 - q**2 - 7) + g)
out.append(f"1. printed formula == A7/(q(q-1)) + B7 S_W4 + C7 S_K6P6 + g(q) (g_formula.json): {sp.expand(hyp - recon) == 0}")
# also check the audit's A7 against audit_tt7_explicit.out's monomial list
A7_list = 596526*q - 911849*q**2 + 362505*q**3 - 50228*q**4 + 3136*q**5 - 91*q**6 + q**7
out.append(f"   A7 = q(q-1)(q-3)(q^4-87q^3+2785q^2-38827q+198842) equals audit_tt7_explicit.out's monomials: {sp.expand(A7 - A7_list) == 0}")

forms = {}
for line in open(Ff):
    p, x84, x86, x46, nb = map(int, line.split()); forms[p] = (x84, x86, x46, nb)
T5 = {}
for line in open(T5f):
    p, t = map(int, line.split()); T5[p] = t
Ys = {}
for line in open(Yf):
    r = list(map(int, line.split())); Ys[r[0]] = r[1:]
TT7 = {p: p * (p - 1) // 2 * t for p, t in T5.items()}
fh = sp.lambdify((q, a, c, d, b2, Y02, Y05, Y09, Y18), hyp, 'math')
good, bad, badlist = 0, 0, []
for p in sorted(TT7):
    if p not in Ys or p not in forms or p < 7: continue
    y = Ys[p]
    x84, x86, x46, nb = forms[p]
    val = int(sp.Integer(p**5 - 90*p**4 + 3116*p**3 - 49252*p**2 + 322374*p - 673743) - 2*(14*p*p - 701*p + 7254)*x84
              + 46*x86 - 83*x46 - 64*(p - 31)*nb - 64*y[0] - 48*y[2] - 64*y[4] + 96*y[5])
    if (1 << 21) * TT7[p] == p * (p - 1) * val: good += 1
    else: bad += 1; badlist.append(p)
out.append(f"2. 2^21 TT7(P_q) == q(q-1)[formula] with the auditor's TT7 counts and Y sums: {good} primes agree, {bad} disagree "
           f"(q = 3 mod 4, {min(TT7)}..{max(TT7)}){'; failures: ' + str(badlist[:20]) if bad else ''}")
lane = {}
for fn in ('paley_tt7_data.out', 'paley_tt7_data2.out'):
    for line in open(os.path.join(HERE, '..', '..', fn)):
        m = re.match(r'p=(\d+) TT7=(\d+)', line)
        if m: lane[int(m.group(1))] = int(m.group(2))
common = sorted(set(lane) & set(TT7))
mism = [p for p in common if lane[p] != TT7[p]]
out.append(f"3. lane's counted TT7 values: {len(lane)} primes {min(lane)}..{max(lane)}; auditor's counts agree at {len(common) - len(mism)}/{len(common)}"
           f"{'; mismatches ' + str(mism) if mism else ''}")
out.append(f"   zeros of TT7 among counted primes: {[p for p in sorted(TT7) if TT7[p] == 0]}")
out.append(f"   sample: TT7(P_31) = {TT7.get(31)}, TT7(P_43) = {TT7.get(43)}, TT7(P_811) = {TT7.get(811)}, TT7(P_{max(TT7)}) = {TT7[max(TT7)]}")
print('\n'.join(out))
open(os.path.join(HERE, 'tt7_formula_check.out'), 'w').write('\n'.join(out) + '\n')
