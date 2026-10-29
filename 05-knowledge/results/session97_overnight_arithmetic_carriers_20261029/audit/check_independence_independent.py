"""Independent square-set and dense F2 matrix verifier, no producer imports."""
import json, pathlib
from fractions import Fraction as F
from math import isqrt

def require(x,msg):
    if not x:raise ArithmeticError(msg)

def prime(p):return p>2 and all(p%d for d in range(2,isqrt(p)+1))
def dense_rank(rows,n):
    a=[r[:] for r in rows];k=0
    for col in range(n):
        idx=next((i for i in range(k,len(a)) if a[i][col]),None)
        if idx is None:continue
        a[k],a[idx]=a[idx],a[k]
        for i in range(len(a)):
            if i!=k and a[i][col]:a[i]=[x^y for x,y in zip(a[i],a[k])]
        k+=1
    return k

def check(c,z):
    a1,a2,a3,a4,a6=map(int,c['ainvs'])
    b2=a1*a1+4*a2;b4=a1*a3+2*a4;b6=a3*a3+4*a6
    b8=a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4
    disc=-b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6
    require(disc==int(c['discriminant']),'recomputed discriminant')
    pts=[tuple(map(F,p)) for p in c['points']]
    for x,y in pts:require(y*y+a1*x*y+a3*y==x**3+a2*x*x+a4*x+a6,'membership')
    rows=[]
    for row in z['independent_rows']:
        p,r=row['p'],row['root']
        require(prime(p) and disc%p!=0,'odd good prime')
        require((4*r**3+b2*r*r+2*b4*r+b6)%p==0,'root')
        squares={(i*i)%p for i in range(1,p)}
        vector=[]
        for x,y in pts:
            require(x.denominator%p!=0,'integral reduction')
            inv=next(i for i in range(1,p) if x.denominator*i%p==1)
            v=(x.numerator*inv-r)%p
            require(v!=0,'ordinary character chart')
            vector.append(int(v not in squares))
        require(sum(b*2**i for i,b in enumerate(vector))==row['bits'],'named row bits')
        rows.append(vector)
    rk=dense_rank(rows,len(pts))
    require(rk==len(pts)==z['character_rank'],'full rank')
    p=z['two_torsion_exclusion_prime']
    require(prime(p) and disc%p!=0,'torsion exclusion good prime')
    require(not any((4*x**3+b2*x*x+2*b4*x+b6)%p==0 for x in range(p)), 'no rational 2 torsion')
    return {'curve':c['id'],'certified_rank_lower_bound':rk,'two_torsion_exclusion_prime':p}

def main():
    root=pathlib.Path(__file__).resolve().parents[1]/'ranks'
    certs=json.loads((root/'character_certificates.json').read_text())
    print(json.dumps([check(json.loads((root/f"curve{z['curve_id']}.json").read_text()),z) for z in certs],sort_keys=True))
if __name__=='__main__':main()
