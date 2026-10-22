"""Independent modular enumeration and exact proof-side checks; no producer imports."""
from fractions import Fraction as Q
from math import comb
from itertools import product
import hashlib,json,pathlib

def check(condition,label):
    if not condition: raise RuntimeError(label)
def val2(n):
    a=0
    while n%2==0: a+=1;n//=2
    return a,n

def bad_residue(n,h,s):
    # Once total valuation reaches s, the bad event is impossible.
    total=0;first=None
    for _ in range(h):
        a,n=val2(3*n+1); total+=a
        if first is None:first=a
        if total>=s:return False,first
    return True,first

rows=[]; modular_checks=0
register_expected=[(1,1),(7,7),(1,26),(7,22),(1,21),(7,12)]
for h in range(1,9):
    s=(3**h).bit_length();M=2**s
    # k has period lcm(12,M)=3M; enumerate complete odd-exponent periods.
    p=3;count=[0]*6;bad=[0]*6;images=set()
    for k in range(1,3*M,2):
        n=(p+1)//4
        images.add(n)
        flag,first=bad_residue(n,h,s)
        index=(k%12)//2
        count[index]+=1;bad[index]+=int(flag)
        check((first==1)==(k%4==3),'first valuation guard')
        if h<=5:
            u=(3**k+1)//4
            check((u%15,u%35)==register_expected[index],'literal source register')
        p=p*9%(4*M);modular_checks+=1
    check(images==set(range(1,M,2)),'exponent/source residue bijection')
    expected=sum((Q(comb(j-1,h-1),2**j) for j in range(h,s)),Q())
    check(Q(sum(bad),sum(count))==expected,'full modular frequency vs composition law')
    conditionals=[Q(a,b) for a,b in zip(bad,count)]
    check(all(q<=2*expected for q in conditionals),'all six guarded bounds')
    check(conditionals[0]==conditionals[2]==conditionals[4],'k=1mod4 guards')
    check(conditionals[1]==conditionals[3]==conditionals[5],'k=3mod4 guards')
    rows.append({'h':h,'valuation_threshold':s,'period':3*M,'bad_density':str(expected),'six_guarded_bad_densities':list(map(str,conditionals))})

# Independent finite-binomial identity and Chernoff fifth-power certificate.
for h in range(1,301):
    s=(3**h).bit_length()
    q=sum((Q(comb(j-1,h-1),2**j) for j in range(h,s)),Q())
    tail=Q(sum(comb(s-1,j) for j in range(h,s)),2**(s-1))
    check(q==tail,'negative-binomial identity')
    check(q**5*84375**h<=65536**h,'rational exponential bound')

# Actual trajectories: use the affine carry recurrence in rational arithmetic.
trajectory_checks=0
for k in range(1,602,2):
    n=(3**k+1)//4;x=n;carry=Q();A=0
    for h in range(1,41):
        a,x=val2(3*x+1); A+=a;carry=(3*carry+1)/2**a
        check(x==Q(3**h*n,2**A)+carry,'affine endpoint identity')
        B=Q(3,2)**h-1;s=(3**h).bit_length();T=B/(1-Q(3**h,2**s))
        check(0<carry<=B,'uniform carry bound')
        if n>T:check((x<n)==(A>=s),'finite-source-exception threshold')
        trajectory_checks+=1
# Worst-case uniform bound reached by all valuation letters one.
for h in range(1,301):
    check(sum((Q(3,2)**j/2 for j in range(h)),Q())==Q(3,2)**h-1,'sharp carry bound')

src=pathlib.Path(__file__).resolve().parent/'frozen_inputs'
inputs={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(src.glob('*')) if p.is_file()}
print(json.dumps({'status':'ACCEPTED','input_sha256':inputs,'modular_source_cases':modular_checks,'literal_endpoint_cases':trajectory_checks,'exact_tail_depths':300,'modular_frequency_rows':rows,'scope':'Natural odd-exponent density of some finite descent, including each source register; no pointwise convergence or original-source credit.'},indent=2,sort_keys=True))
