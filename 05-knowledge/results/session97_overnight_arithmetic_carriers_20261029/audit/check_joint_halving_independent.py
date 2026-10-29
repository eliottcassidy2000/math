"""Independent local characters, full binary kernels and actual joint halves."""
import json
from pathlib import Path
from fractions import Fraction as F
from collections import Counter

def need(b,m):
    if not b:raise ArithmeticError(m)
def rank(rows):
    basis={}
    for x in rows:
        while x:
            i=x.bit_length()-1
            if i in basis:x^=basis[i]
            else:basis[i]=x;break
    return len(basis)
def plus(P,Q,A,B):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    if x==u and y+v==0:return None
    slope=(3*x*x+2*A*x+B)/(2*y) if P==Q else (v-y)/(u-x)
    z=slope*slope-A-x-u
    return z,slope*(x-z)-y

def main():
    root=Path(__file__).resolve().parents[1]/'halving'
    packet=json.loads((root/'halving_certificate.json').read_text())
    kernels={r['curve_id']:r for r in json.loads((root/'binary_kernel_certificate.json').read_text())}
    out=[]
    for entry in packet['curves']:
        curve=entry['curve'];ident=curve['id'];a1,a2,a3,a4,a6=map(F,curve['ainvs'])
        A=a2+a1*a1/4;B=a4+a1*a3/2;C=a6+a3*a3/4
        pts=[]
        for x,y in [tuple(map(F,P)) for P in curve['points']]:
            Y=y+(a1*x+a3)/2
            need(Y*Y==x**3+A*x*x+B*x+C,'point membership');pts.append((x,Y))
        disc=16*(A*A*B*B-4*B**3-4*A**3*C-27*C*C+18*A*B*C)
        bits=[]
        for row in entry['character_rows']:
            p,r=row['p'],row['root']
            need(p>2 and all(p%d for d in range(2,int(p**0.5)+1)),'odd prime')
            red=lambda q:q.numerator*pow(q.denominator,-1,p)%p
            need(red(disc)!=0,'good reduction');need((r**3+red(A)*r*r+red(B)*r+red(C))%p==0,'2torsion root')
            squares={i*i%p for i in range(1,p)};mask=0
            for j,(x,Y) in enumerate(pts):
                val=(red(x)-r)%p;need(val!=0,'ordinary character chart')
                if val not in squares:mask|=1<<j
            need(mask==row['bits'],'recomputed row');bits.append(mask)
        need(rank(bits)==17,'rank17 bank');need(all(any(b>>j&1 for b in bits) for j in range(30)),'no individual doubles')
        kc=kernels[ident];basis=kc['kernel_basis_masks']
        need(rank(basis)==13 and all((b&r).bit_count()%2==0 for b in basis for r in bits),'full kernel basis')
        span={0}
        for b in basis:span|={x^b for x in span}
        need(len(span)==8192,'all kernel words')
        enum=Counter(v.bit_count() for v in span)
        need(dict(enum)=={int(k):v for k,v in kc['weight_enumerator'].items()},'complete weight enumeration')
        minimum=min(v.bit_count() for v in span if v)
        support={tuple(i for i in range(30) if v>>i&1) for v in span if v.bit_count()==minimum}
        need(support==set(map(tuple,kc['minimum_supports_zero_based'])),'all minimum supports')
        joint=json.loads((root/f'joint_candidate_{ident}.json').read_text())
        terms=[pts[i] if s==1 else (pts[i][0],-pts[i][1]) for i,s in zip(joint['indices_zero_based'],joint['signs'])]
        # Balanced summation, separately from producer's chronological addition.
        while len(terms)>1:terms=[plus(terms[k],terms[k+1],A,B) if k+1<len(terms) else terms[k] for k in range(0,len(terms),2)]
        target=terms[0];need(target==tuple(map(F,joint['completed_square_sum'])),'joint sum')
        x,y=map(F,joint['native_half']);Y=y+(a1*x+a3)/2
        need(Y*Y==x**3+A*x*x+B*x+C,'native half membership')
        need(plus((x,Y),(x,Y),A,B)==target,'actual rational halving identity')
        need(plus((x,-Y),(x,-Y),A,B)!=target,'opposite sign hostile')
        out.append({'curve':ident,'local_rank':17,'kernel_dimension':13,'minimum_weight':minimum,'minimum_support_count':len(support),'joint_half_verified':True})
    print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
