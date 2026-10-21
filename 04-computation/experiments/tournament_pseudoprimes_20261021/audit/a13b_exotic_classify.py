# AUDIT a13b: classify the 60 exotic squarefree Camion pseudoprimes <= 20000 (from a13_list20000.out) by residues used.
from sympy import factorint
import re, ast
txt = open("a13_list20000.out").read()
odd = ast.literal_eval(re.search(r"exotic odd: (\[.*\])", txt).group(1))
even = ast.literal_eval(re.search(r"exotic even: (\[.*\])", txt).group(1))
for n in odd + even:
    res = {p: (n - 1) % (p - 1) + 1 for p in factorint(n) if p > 2}
    spor = {p: r for p, r in res.items() if r >= 4}
    print(n, factorint(n), res, "SPORADIC " + str(spor) if spor else "")
