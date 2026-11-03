"""Exact all-exponent intersections of positive shifted dyadic rays."""
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
def need(ok,message):
    if not ok:raise ValueError(message)
def v2(n):
    n=abs(n);need(n>0,'nonzero valuation');return (n&-n).bit_length()-1
def pow2(n):return n>0 and n&(n-1)==0
def intersections(U,V,A,B):
    need(U>0 and A>0,'positive slopes')
    u,v=v2(U),v2(A);a,b=U>>u,A>>v;D=B-V;out=set()
    if not D:
        return dict(infinite=(a==b),shift=u-v if a==b else None,points=[])
    def add(M,N):
        if M>=u and N>=v:
            m,n=M-u,N-v;need(U*2**m+V==A*2**n+B,'exact solution');out.add((m,n))
    t=v2(D)
    if a!=b and D%(a-b)==0 and pow2(D//(a-b)):
        e=(D//(a-b)).bit_length()-1;add(e,e)
    z=a*2**t-D
    if z>0 and z%b==0 and pow2(z//b):
        N=(z//b).bit_length()-1
        if N>t:add(t,N)
    z=D+b*2**t
    if z>0 and z%a==0 and pow2(z//a):
        M=(z//a).bit_length()-1
        if M>t:add(M,t)
    need(len(out)<=2,'at most two isolated intersections')
    return dict(infinite=False,shift=None,points=sorted(out))

def run():
    checks=0
    for U in range(1,13):
      for A in range(1,13):
       for V in range(-4,5):
        for B in range(-4,5):
            r=intersections(U,V,A,B);actual={(m,n) for m in range(9) for n in range(9) if U*2**m+V==A*2**n+B}
            predicted={(m,n) for m in range(9) for n in range(9) if n==m+r['shift']} if r['infinite'] else {tuple(p) for p in r['points'] if max(p)<9}
            need(actual==predicted,'independent bounded grid');checks+=1
    need(intersections(3,0,1,2)['points']==[(0,0),(1,2)],'sharp two-intersection control')
    inputs={1457:(1188,269),1459:(1846,-387),1463:(1240,223),1465:(1076,389)}
    outputs=[(h,c,3*inputs[h][0],3*inputs[h][1]+c) for h,c in ((1457,2),(1457,4),(1459,2),(1459,4),(1463,4),(1465,4))]
    data=[]
    for h,c,U,V in outputs:
      for g,(A,B) in inputs.items():
        data.append(dict(source=h,offset=c,target_input=g,result=intersections(U,V,A,B)))
    out=dict(checks=checks,inputs=inputs,outputs=outputs,intersections=data)
    (HERE/'ray_gate.json').write_text(json.dumps(out,indent=2)+'\n')
    print('exact dyadic gate',checks,'nonempty',[(d['source'],d['offset'],d['target_input'],d['result']) for d in data if d['result']['points'] or d['result']['infinite']],flush=True)
    return out
if __name__=='__main__':run()
