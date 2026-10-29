"""Valuation-adjusted, source-marked reverse decoder for growing returns."""
from fractions import Fraction as F
from pathlib import Path
import json

# In m=(n+1)/2 coordinates: 3^e (a*m+b)/2^d, e+v3(a)=J.
B={1:(4,1,1,6,4,63,64),2:(6,3,11,11,7,679,2048),
   3:(5,9,25,11,7,1135,2048),4:(4,27,59,11,7,1439,2048)}
SIGNATURES={61:1,13:2,70:3,34:4}
def valuation(n,p):
    if n<=0:raise ValueError('positive integer required')
    v=0
    while n%p==0:n//=p;v+=1
    return v
def signature(m):
    v=valuation(m,3)
    return (m//3**v)*pow(2,v,81)%81
def forward(m,i):
    if i==0:
        if m%2:raise ValueError('zero-letter guard')
        return 3*m//2
    e,a,b,d,J,r,mod=B[i]
    if m%mod!=r:raise ValueError('nonzero-letter guard')
    y,rem=divmod(3**e*(a*m+b),2**d)
    if rem:raise ArithmeticError('inexact forward return')
    return y
def decode(m):
    """Unique candidate B_i followed by z zero letters from positive 81Z.

    None means there is no such suffix. This does not choose a stopping
    source: a marked source can occur inside the decoded zero run.
    """
    i=SIGNATURES.get(signature(m))
    if i is None:return None
    e,a,b,d,J,r,mod=B[i];z=valuation(m,3)-e
    if z<0:return None
    y=(2**z)*(m//3**z)
    numerator=(2**d)*y-3**e*b
    x,rem=divmod(numerator,3**e*a)
    if rem or x<=0 or x%81 or x%mod!=r:return None
    if forward(x,i)!=y:return None
    return dict(source_m=x,letter=i,zeros=z,odd_clock=J+z,division_clock=d+z)

def native(n):
    q=1;word=[]
    while q:
        n=3*n+1;a=0
        while n%2==0:n//=2;a+=1
        word.append(a);q+=-1+3*(a%3==0)
    return n,word

def verify():
    checks=0;examples=0;depths=0
    def need(x,msg):
        nonlocal checks
        checks+=1
        if not x:raise ValueError(msg)
    words={1:[3,1,1,1],2:[3,3,1,1,1,1,1],3:[3,1,3,1,1,1,1],4:[3,1,1,3,1,1,1]}
    residues=[]
    for t in range(1,5):
        mod=3**t
        residues.append([pow(2,e,mod)*b*pow(2**d,-1,mod)%mod for e,a,b,d,J,r,q in B.values()])
    need(residues==[[1,1,1,1],[7,4,7,7],[7,13,16,7],[61,13,70,34]],'least tested signature precision')
    for i,(e,a,b,d,J,r,mod) in B.items():
        # All stable inputs up to a million; every legal zero suffix of each.
        first=81*(r*pow(81,-1,mod)%mod)
        for x in range(first,1000001,81*mod):
            y=forward(x,i);h,w=native(2*x-1)
            need(h==2*y-1 and w==words[i],'native word/endpoint');examples+=1
            need(valuation(y,3)==e and signature(y) in SIGNATURES,'valuation reset')
            for z in range(valuation(y,2)+1):
                endpoint=3**z*y//2**z
                rec=decode(endpoint)
                need(rec==dict(source_m=x,letter=i,zeros=z,odd_clock=J+z,division_clock=d+z),'unique suffix and both clocks')
                need(signature(endpoint)==signature(y),'zero-run invariant');depths+=1
                # Exact marked-source debt before and after undoing this suffix.
                n0=2*x-1;n1=2*endpoint-1
                need(n1>n0,'reverse descent is forward growth, not repayment')
    # Decode for an independent finite endpoint universe and replay all accepted results.
    accepted=0
    for endpoint in range(1,100001):
        rec=decode(endpoint)
        if rec is None:continue
        y=forward(rec['source_m'],rec['letter'])
        for _ in range(rec['zeros']):y=forward(y,0)
        need(y==endpoint,'endpoint enumeration replay');accepted+=1
    need(forward(63,1)==81 and decode(81) is None,'stable-domain hypothesis is necessary')
    need(forward(4671,1)==forward(5535,4)==5913,'27-divisible sources have a genuine last-letter collision')
    need(valuation(4671,3)==valuation(5535,3)==3,'81 layer is sharp among divisibility layers')
    need(forward(6399,1)==8100 and forward(8100,0)==12150,'source stopping hostile')
    need(decode(12150)['source_m']==6399,'decoder would jump past marked source 8100')
    N=[1,2,3,4,6,9,13,18,24,31,40]
    for k in range(11,75):N.append(N[k-1]+N[k-6]+3*N[k-11])
    r=json.loads((Path(__file__).parent/'bound74.json').read_text())
    lower=r['exclusive_bound'];ell=((lower+1)//2).bit_length()-1
    need(r['k']==74 and not r['counterexample'] and lower==49*(1090945607561904127+1)-1,'two-stage consequence from exact s60')
    bound=F(N[74],4**74)+sum((F(N[j+1],2*4**j) for j in range(ell,74)),F())
    need(bound==F(641060333941629,2**148),'all-original-source mass bound')
    return dict(status='PASS',explicit_checks=checks,stable_sources=examples,legal_suffixes=depths,
                endpoint_universe=100000,accepted_endpoints=accepted,signature_residues=residues,
                bound74=dict(least_source_lower_bound=lower,ell_lower_bound=ell,mass_upper_bound=str(bound)),
                source_hostile=dict(endpoint_m=12150,macro_source_m=6399,marked_source_inside_suffix_m=8100),
                weaker_layer_hostile=dict(source_m=[4671,5535],letters=[1,4],endpoint_m=5913))

if __name__=='__main__':print(json.dumps(verify(),indent=2,sort_keys=True))
