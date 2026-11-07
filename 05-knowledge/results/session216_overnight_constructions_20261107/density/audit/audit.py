from fractions import Fraction as F
from itertools import groupby
import json
from pathlib import Path

def log_bounds(x, terms=30):
    z=(x-1)/(x+1)
    s=2*sum((z**(2*j+1)/F(2*j+1) for j in range(terms)),F(0))
    return s,s+2*z**(2*terms+1)/(F(2*terms+1)*(1-z*z))
a,b=log_bounds(F(2)); c,d=log_bounds(F(3))
assert F(6931,10000)<a<b<F(6932,10000)
assert d+2*b<F(2485,1000)
assert F(524,100)**3<144
first=2+2*F(6932,10000)
second=(F(3374,100)-12*F(6931,10000)+8*F(2485,1000))/12
third=(2+9*F(6932,10000))/(16*F(524,100))
assert first+second+third<F(73,10)
assert F(292,10)**3*4 < 60**3*F(6931,10000)**2
# Independently read each word from its high bit and append bits, without a g-table.
def direct(v):
    ans,prefix=1,0
    for ch in bin(v)[2:]:
        if ch=='0': ans+=prefix+(prefix%2)
        prefix=2*prefix+int(ch)
    return ans
original={r['n']:r for r in json.loads((Path(__file__).resolve().parents[1] / 'check.json').read_text())}
rows=[]
for n in range(12,19):
    m=max(j for j in range(1,n) if n*2**(3*(j-1))<=2**n)
    C=2**m; P=2**(n-m-1)
    assert P>=8*C
    values=[direct(v) for v in range(2**(n-1),2**n)]
    sizes=[len(list(group)) for _,group in groupby(sorted(values))]
    energy=sum(s*s for s in sizes)
    assert energy==original[n]['energy'] and len(sizes)==original[n]['distinct']
    parity_mutant_failures=0
    for v,g in zip(range(2**(n-1),2**n),values):
        h,l=divmod(v,C)
        rev=int(format(l,f'0{m}b')[::-1],2)
        c=C-1-rev
        uncorrected=direct(h)+direct(c)+h*c-1
        assert g==uncorrected+(h%2)*(c%2)
        parity_mutant_failures+=(g!=uncorrected)
        assert 1<=g<=v+1
    assert parity_mutant_failures>0
    # Literal ordered collision enumeration at boundary n=12, no histogram identity.
    if n==12: assert sum(x==y for x in values for y in values)==energy
    rows.append(dict(n=n,m=m,energy=energy,support=len(sizes),parity_mutant_failures=parity_mutant_failures))
print(json.dumps(dict(rational_upper_bound=str(first+second+third),rational_upper_bound_decimal=float(first+second+third),log_certificates='PASS',spectrum_constant_certificate='PASS',independent_rows=rows),indent=2))
