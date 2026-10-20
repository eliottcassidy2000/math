#!/usr/bin/env python3
"""Exact finite base cases and monotone lower bound for N <= K_n^2."""
from fractions import Fraction
from math import factorial
import json
from pathlib import Path
p=Path(__file__).resolve().parent
small=[(3,2),(4,3),(5,10),(6,34)]
rows=[]
for n,K in small:
    N=1<<((n-1)*(n-2)//2)
    assert N<=K*K
    rows.append({'n':n,'N':N,'K':K,'K_squared':K*K,'pass':True})
n=7
R=Fraction(1<<((n-1)*(n+2)//2),4*factorial(n)**2)
assert R>=1
# For all n >= 7, 2^(n+1)/(n+1)^2 >= 4. Its successive ratio is
# 2*((n+1)/(n+2))^2 > 1, since 2*(n+1)^2 > (n+2)^2 for n >= 2.
assert Fraction(1<<8,8**2)==4
report={'status':'PASS','finite_base_cases':rows,'n7_lower_bound_ratio':str(R),'lower_bound':'K_n >= 2^(n(n-1)/2)/(2 n!)','ratio':'K_n^2/N >= R_n = 2^((n-1)(n+2)/2)/(4(n!)^2)','growth':'R_(n+1)/R_n = 2^(n+1)/(n+1)^2 >= 4 for n >= 7','all_order_conclusion':'N <= K_n^2 for every n >= 3','certificate_consequence':'O(N n^2) = O(K_n^2 log K_n) bits, polynomial in explicit matrix/support size'}
(p/'certificate_size_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
