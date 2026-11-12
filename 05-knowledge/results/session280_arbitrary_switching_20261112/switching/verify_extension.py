"""Exact check of the151-letter supercritical-candidate switching metric."""
from verify import *
def main():
    raw=(ROOT/'metric151.json').read_bytes();data=json.loads(raw);words=data['words'];tau=F(*data['tau']);D=data['Q_den']
    check(len(words)==151 and len(set(words))==151,'alphabet151 unique')
    check(all(2<=len(w)<=9 and w[0]=='1' and w[-1]=='0' and sum(a!=b for a,b in zip(w,w[1:]))<=3 for w in words),'active lowoscillation')
    excluded={'10000010','11100000','11110000','11111000','11111110','100000010','111000000','111100000','111110000','111111000','111111110'}
    check(words==[w for w in active(9) if sum(a!=b for a,b in zip(w,w[1:]))<=3 and w not in excluded],'exact151 language')
    Q=[((F(a,D),F(b,D)),(F(b,D),F(c,D))) for a,b,c in data['Q_num']];M=list(map(block,words));ell=F(1,4)
    for q in Q:check(psd(sub(q,((ell,0),(0,ell))),True),'uniform norm lower bound')
    minimum=None;tests=0
    for u,U in enumerate(M):
        for v,V in enumerate(M):
            q=V[0][1];t=V[1][1]
            for r in R:
                a,w,d=mv(V,mv(U,r));rr=F(w,a);s=F(d,a);K=((t-q*rr,rr),(-q*s,s))
                loss=sub(scale(tau,Q[u]),mm(tuple(zip(*K)),mm(Q[v],K)))
                check(psd(loss,True),('metric inequality',words[u],words[v],r))
                det=loss[0][0]*loss[1][1]-loss[0][1]**2;minimum=det if minimum is None else min(minimum,det);tests+=1
    initial=F(0)
    for u,U in enumerate(M):
        rays=[mv(U,r) for r in R];positive=[(a,F(w+d,2),F(w-d,2)) for a,w,d in rays]
        for j in range(3):check(max(r[j] for r in positive)<=2*min(r[j] for r in positive),'ratio2')
        for i,j in product(range(3),repeat=2):
            if i==j:continue
            x,y=rays[i],rays[j];z=tuple(F(y[t])-F(y[0],x[0])*x[t] for t in (1,2));initial=max(initial,quad(Q[u],z))
    B=96;k=1;bound=4*B*B*initial/ell
    while bound>=1:k+=1;bound*=tau
    out={'status':'PROVED rational switched metric; not scalar coverage','alphabet_size':151,'metric_PSD_tests':tests,'tau':str(tau),'ell':str(ell),'initial_energy':str(initial),'B':B,'uniform_threshold':k,'squared_width_bound_decimal_display':float(bound),'minimum_PSD_determinant':str(minimum),'certificate_sha256':hashlib.sha256(raw).hexdigest()}
    (ROOT/'extension_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
