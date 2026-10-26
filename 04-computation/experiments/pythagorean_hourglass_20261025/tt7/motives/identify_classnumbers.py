# Do class numbers fill the new directions?  Columns q^j h(D) for D in {-q, -4q, -8q, -3q, -12q, -24q, -2q?} (as
# discriminants; h = qfbclassno, the form class number) and Hurwitz H(4q), H(q) (via gp), j <= 2.
import re, os, subprocess, tempfile
exec(open('identify_new.py').read().split('results = []')[0].replace('print(f"{len(qs)} primes', 'pass  # print(f"{len(qs)} primes'))
gpf = os.path.join(tempfile.gettempdir(), 'cls.gp')
Ds = {'h(-q)': '-q', 'h(-4q)': '-4*q', 'h(-8q)': '-8*q', 'h(-3q)': '-3*q', 'h(-12q)': '-12*q', 'h(-24q)': '-24*q', 'h(-16q)': '-16*q'}
prog = 'default(parisize,100000000);\n' + ''.join(
    f'print("{k} ", vector({len(qs)}, i, my(q = {qs}[i]); qfbclassno({v})));\n' for k, v in Ds.items()
) + f'print("H(4q) ", vector({len(qs)}, i, my(q = {qs}[i]); qfbhclassno(4*q)));\nprint("H(q) ", vector({len(qs)}, i, my(q = {qs}[i]); qfbhclassno(q)));\nquit;\n'
open(gpf, 'w').write(prog)
out = subprocess.run(['gp', '-q', gpf], capture_output=True, text=True).stdout.strip().splitlines()
cls = {}
for l in out:
    k, v = l.split(' ', 1); cls[k] = [eval(x) for x in v.strip('[]').split(',')]
print({k: v[:6] for k, v in cls.items()})
from fractions import Fraction
for k, v in cls.items():
    cols = [[Fraction(x) * q ** j for x, q in zip(v, qs)] for j in range(3)]
    rBC = rank(base + cols); d = rank(base + cols + newcols) - rBC
    print(f"{k}: new dimension {rBN - rB} -> {d}")
allc = []
for k, v in cls.items(): allc += [[Fraction(x) * q ** j for x, q in zip(v, qs)] for j in range(3)]
rBC = rank(base + allc); print("all class-number columns:", rBN - rB, "->", rank(base + allc + newcols) - rBC)
