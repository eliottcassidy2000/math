"""Structural clone and marked-halving controls for the actual 6325 return."""
from return6325 import HERE,need,junction,phases,values
import json
def clone(T,q,m):
    ad=[(i,j) for i in range(len(T)) for j in range(m if i==q else 1)]
    rows=tuple(sum(1<<b for b,(v,k) in enumerate(ad) if (i==v==q and j<k) or (i!=v and T[i]>>v&1)) for i,j in ad)
    return rows,ad
def main():
    d=json.loads((HERE/'return6325.json').read_text());d2=json.loads((HERE/'return6325_plus2.json').read_text());q=d['root'];out=[]
    for side,T in (('source',tuple(d['source'])),('target_plus4',tuple(d['target'])),('target_plus2',tuple(d2['target']))):
        sample=[];keep=[i for i in range(len(T)) if i!=q];A=junction.sub(T,tuple(keep));p=next(i for i in range(len(T)) if T[i]>>q&1);p=keep.index(p)
        a,s,dd,F=junction.profile(A,p)
        for m in range(1,5):
            U,ad=clone(T,q,m);r=ad.index((q,0));p=junction.profile(U,r)
            need(values(U,r)[1]==p,'independent marked clone profile')
            need(p[0]==a+F-s+(2**m-2)*(F-dd),'all-size visit law at literal control')
            need(p[3]-p[0]==(s if m==1 else dd),'first clone changes readiness to parent deletion count')
            sample.append(dict(m=m,profile=p));print(side,m,p,flush=True)
        H=[z['profile'][0] for z in sample];A=H[1]-H[0];B=2*H[0]-H[1]
        need(all(h==A*2**i+B for i,h in enumerate(H)),'capacity-one scalar recurrence')
        out.append(dict(side=side,A=A,B=B,parent_profile=[a,s,dd,F],samples=sample))
    (HERE/'clone_probe.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
