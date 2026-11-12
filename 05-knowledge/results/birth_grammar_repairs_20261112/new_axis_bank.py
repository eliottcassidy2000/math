"""Exact selected-start/two-cover run DP for the new strong birth layers."""
from common import *
from hamilton_clocks import trim
from numba import njit
from fractions import Fraction
from collections import Counter
import time

@njit(cache=True)
def slope_skeleton(rows,O,p):
    n=len(rows);verts=np.array([v for v in range(n) if v!=p],dtype=np.int64);m=n-1
    degree=0
    for v in range(n):degree+=(rows[p]>>v)&1
    R=min(n-degree,degree+2);N=1<<m
    # channels: ordinary uncut, selected-start uncut, one selected-start cut.
    dp=np.zeros((N,m+1,R+1,3),dtype=np.int64)
    dp[0,m,1,0]=1
    for j in range(m):
        dp[1<<j,j,0,0]=1
        if O>>verts[j]&1:dp[1<<j,j,0,1]=1
    for mask in range(N):
        for r in range(R+1):
            for j in range(m+1):
                v=p if j==m else verts[j]
                if j!=m and r<R and rows[v]>>p&1:
                    for q in range(3):dp[mask,m,r+1,q]+=dp[mask,j,r,q]
                for k in range(m):
                    if mask>>k&1:continue
                    mm=mask|(1<<k)
                    if rows[v]>>verts[k]&1:
                        for q in range(3):dp[mm,k,r,q]+=dp[mask,j,r,q]
                    if O>>verts[k]&1:dp[mm,k,r,2]+=dp[mask,j,r,0]
    out=np.zeros(R,dtype=np.int64)
    for r in range(1,R+1):
        for j in range(m+1):out[r-1]+=dp[N-1,j,r,1]+dp[N-1,j,r,2]
    return out

def run(top=12):
    start=time.time();out=[];controls=0;hist=Counter()
    for n in range(11,top+1):
        for idx,r in enumerate(read(f'birth_extensions_{n}.json.gz')['records']):
            if r['U']%2==0:continue
            hist[n]+=1
            for p in r['stable_ports']:
                u=[int(x) for x in slope_skeleton(np.array(r['rows'],dtype=np.int64),r['selected'],p)]
                need(u[0]==r['U'] and all(c>=0 for c in u),'positive U DP and base control')
                if idx%17==0:
                    for a in [2,3]:
                        labs,rr=ib.inflate(r['rows'],p,a);h,uu,w,_=tables(np.array(rr,dtype=np.int64));O=[v for v,q in enumerate(labs) if r['selected']>>q&1]
                        need(sum(int(uu[v]) for v in O)==value(expcoeff(u),a),'independent literal slope');controls+=1
                out.append({'birth_id':r['birth_id'],'rows':r['rows'],'selected':r['selected'],'port':p,'Stirling_coefficients':u,'exponential_coefficients':list(trim(expcoeff(u)))})
        print(json.dumps({'order':n,'axes':len(out),'seconds':round(time.time()-start,2)}),flush=True)
    lookup={tuple(r['exponential_coefficients']):i for i,r in enumerate(out)};matches=[]
    for r in read('even_valuations.json.gz')['records']:
        if r['type']!='bounded_constant':continue
        j=r['eventual_valuations_by_parity'][0]
        for shift in range(-8,9):
            e=[Fraction(c,2**j)*Fraction(b+1)**shift for b,c in enumerate(r['exponential_coefficients'])]
            if any(c.denominator!=1 for c in e):continue
            ix=lookup.get(trim([int(c) for c in e]))
            if ix is not None:matches.append({'birth_id':r['birth_id'],'port':r['port'],'valuation':j,'shift':shift,'new_axis_index':ix})
    result={'scope':f'Odd U axes of the retained strong birth extension layers11..{top}; exact positive repeated-port DP, not a complete tournament-size census.','odd_sources_by_order':dict(hist),'axes':len(out),'distinct_profiles':len(lookup),'literal_controls':controls,'matches':matches,'records':out}
    save(f'new_axis_bank_{top}.json.gz',result);print(json.dumps({k:v for k,v in result.items() if k!='records'}),flush=True)
if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 12)
