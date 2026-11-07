"""Two exact one-arc escapes from the Session213 comb; standard library only."""
import json
from counting import adj,count,dp

def need(test,message):
    if not test:raise RuntimeError(message)

def inherited_comb_formula(m,q,tau,reversed_apex):
    # Literal Session213 cycle-to-path/formula.py coefficients, separately pinned.
    v=int(reversed_apex);p=1<<m
    base=p+1+(1-v)*(4*p-2)+tau*(2*p+2)-tau*(1-v)*p
    linear=[]
    for i in range(m):
        u=(1<<(m-1-i))+(1<<i);first=int(i==0);last=int(i==m-1)
        d=(1<<(m+1-i))*3**i+(1<<(i+2))-2*last
        linear.append(d-(1-v)*(u-first+last)-tau*(u+first-last))
    pair={d:(1<<(m-d))*(3**(d-1)+1)+(1<<(d-1))+int(d==1) for d in range(1,m)}
    selected=[i for i in range(m) if q>>i&1]
    return base+sum(linear[i] for i in selected)-sum(pair[j-i] for i in selected for j in selected if i<j)

def run():
    records=[]
    for m,q,tau,ab,r,oldh,newh,expected in [
      (8,14,0,0,0,3187,7477,[1022,1020,2040,4080,4064,4032,1920,1792,1536,1024,3,1991]),
      (9,179,1,1,5,24693,27371,[4094,6140,8184,8176,4064,1984,8064,7936,3584,7168,2048,4130,1329])]:
        L=m+2;a=[int(i==1) for i in range(L)];b=[1,1-tau]+[1-(q>>i&1) for i in range(m)]
        parent=adj(a,b,ab);child_a=a.copy();child_a[r]=1;child=adj(child_a,b,ab)
        changed=[(i,parent[i]^child[i]) for i in range(L+2) if parent[i]!=child[i]]
        need(changed==[(r,1<<L),(L,1<<r)],('single flip',changed))
        need(parent[r]>>L&1 and child[L]>>r&1,'flip direction')
        need(child==expected,'witness rows')
        need(count(a,b,ab)==dp(parent)==inherited_comb_formula(m,q,tau,not ab)==oldh,'parent count')
        need(count(child_a,b,ab)==dp(child)==newh,'child count')
        records.append({'m':m,'q':q,'tau':tau,'apex_orientation':'+' if ab else '-',
          'flip_chain_index':r,'parent_H':oldh,'child_H':newh,'parent_rows':parent,'child_rows':child,
          'methods':['subset dynamic programming','three-run transfer'],'parent_extra_control':'inherited closed comb formula'})
    return {'status':'PASS','scope':'Two finite physical one-arc constructions; no all-target recurrence','cases':records}
if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
