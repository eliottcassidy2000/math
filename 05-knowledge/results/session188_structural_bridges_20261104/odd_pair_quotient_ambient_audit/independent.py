"""Independent physical replay of every retained ambient pair-cube value."""
from pathlib import Path
from functools import lru_cache
import hashlib,json
ROOT=Path(__file__).resolve().parent.parent
S=ROOT/'odd_pair_quotient'
def need(ok,why):
    if not ok: raise RuntimeError(why)
def valid(T):
    n=len(T)
    need(all(0<=r<(1<<n) and not r>>i&1 for i,r in enumerate(T)), 'invalid row')
    need(all(((T[i]>>j)&1)+((T[j]>>i)&1)==1 for i in range(n) for j in range(i)), 'invalid pair')
@lru_cache(None)
def H(T):
    n=len(T);pred=[sum(1<<u for u in range(n) if T[u]>>v&1) for v in range(n)]
    @lru_cache(None)
    def end(mask,v):
        rest=mask^(1<<v)
        if not rest:return 1
        eligible=rest&pred[v];z=0
        while eligible:
            bit=eligible&-eligible;eligible-=bit;z+=end(rest,bit.bit_length()-1)
        return z
    return sum(end((1<<n)-1,v) for v in range(n))
def F(T,a,c):
    n=len(T);R=list(T)+[1<<(n+1),0]
    need(0<=a<(1<<n) and 0<=c<(1<<n),'port mask')
    for v in range(n):
        for port,mask in ((n,a),(n+1,c)):
            u,w=(v,port) if mask>>v&1 else (port,v)
            R[u]|=1<<w
    valid(R);return H(tuple(R))
def substitute(K,r,A):
    outside=[i for i in range(len(K)) if i!=r];n=len(A);R=list(A)+[0]*len(outside)
    for i,u in enumerate(outside):
        for j,v in enumerate(outside):
            if K[u]>>v&1:R[n+i]|=1<<(n+j)
        for v in range(n):
            x,y=(n+i,v) if K[u]>>r&1 else (v,n+i)
            R[x]|=1<<y
    return tuple(R)
def main():
    raw=(S/'ambient.normal.json').read_bytes();need(raw==(S/'ambient.optimized.json').read_bytes(),'producer modes')
    data=json.loads(raw);values=0;traces=[]
    A=(14,12,56,48,35,3);B=(50,48,11,3,44,12)
    # These physical TT2 inflations use 0->1,2->3,4->5 internally.
    valid(A);valid(B)
    for rec in data['contexts']:
        K=tuple(rec['outer_quotient']);r=rec['root'];valid(K)
        T,U=tuple(rec['source']),tuple(rec['target']);valid(T);valid(U)
        need(T==substitute(K,r,A) and U==substitute(K,r,B),'physical root substitution')
        incoming=sum(1<<v for v in range(6,len(T)) if T[v]&1)
        d={};unit=H(substitute(K,r,(2,0)));need(unit==rec['unit'] and unit%2,'odd multiplier')
        for coarse,query,expected in zip(rec['coarse_queries'],rec['physical_queries'],rec['source_target_values']):
            a,c=coarse;want=[incoming|sum(3<<(2*i) for i in range(3) if m>>i&1) for m in (a,c)]
            need(query==want,'physical exterior incidence lift')
            actual=[F(T,*query),F(U,*query)];need(actual==expected,'literal context count');values+=2
            local=[sum(3<<(2*i) for i in range(3) if m>>i&1) for m in (a,c)]
            gap=actual[0]-actual[1];need((gap-unit*(F(A,*local)-F(B,*local)))%8==0,'fixed multiplier per query');d[tuple(coarse)]=gap
        L=[];Q=(2,4,1)
        for j in range(3):
            av=1<<j;z=0
            for i in range(3):
                if i==j:continue
                au=1<<i;z+=d[au,av]-d[au,0]-d[0,av]+d[0,0]
                if Q[j]>>i&1:z+=d[au|av,av]-d[au,av]-d[av,av]+d[0,av]-d[au|av,0]+d[au,0]+d[av,0]-d[0,0]
            L.append(z)
        need(L==rec['L'] and sum(L)%8==4,'complete finite-difference trace');traces.append(sum(L))
    cube=data['ambient_cube'];seen=set()
    for state in cube['states']:
        theta,bits=state['theta'],state['bits'];base=data['contexts'][3]['target' if theta else 'source'];T=list(base)
        for i in range(3):
            if bits>>i&1:T[2*i]^=1<<(2*i+1);T[2*i+1]^=1<<(2*i)
        need(T==state['physical_rows'],'cube physical XOR target')
        actual=[F(tuple(T),*q) for q in cube['queries']];need(actual==state['query_values'],'all cube query values');values+=len(actual)
        need(tuple(actual) not in seen,'distinct complete retained signatures');seen.add(tuple(actual))
        for dec in cube['decoders']:
            if dec['theta']==theta:need(actual[dec['query_index']]==dec['values'][bits>>dec['child']&1],'individual child decoder')
    need(len(seen)==16 and len(cube['queries'])==28 and cube['adaptive_bound']==25,'state and query census')
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [S/'AMBIENT.md',S/'verify_ambient.py',S/'ambient.normal.json',S/'ambient.optimized.json',Path(__file__)]}
    return dict(status='PASS',method='No producer imports. Backward-ending Hamilton recursion, physical reconstruction, literal full queries.',contexts=len(traces),traces=traces,cube_states=len(seen),literal_retained_values=values,source_sha256=pins)
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
