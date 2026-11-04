"""Exact named-port controls; Python standard library only."""
from pathlib import Path
from itertools import combinations, permutations
from functools import lru_cache
import json, hashlib
HERE=Path(__file__).resolve().parent

def need(x,msg):
    if not x: raise RuntimeError(msg)
def valid(T):
    n=len(T)
    need(all(not(T[i]>>i&1) for i in range(n)), 'loop')
    need(all(((T[i]>>j&1)+(T[j]>>i&1))==1 for i in range(n) for j in range(i)), 'tournament')
def rev(T):
    m=(1<<len(T))-1
    return tuple(m^(1<<i)^r for i,r in enumerate(T))
def birth(T,cut):
    n=len(T);return tuple(r|(((cut>>i)&1)<<n) for i,r in enumerate(T))+(((1<<n)-1)^cut,)
def counter(T):
    n=len(T);R=rev(T)
    @lru_cache(None)
    def end(S,v):
        rem=S^(1<<v)
        if not rem:return 1
        cand=rem&R[v];ans=0
        while cand:
            bit=cand&-cand;cand-=bit
            ans+=end(rem,bit.bit_length()-1)
        return ans
    @lru_cache(None)
    def h(S):
        if not S:return 1
        ans=0;u=S
        while u:
            bit=u&-u;u-=bit;ans+=end(S,bit.bit_length()-1)
        return ans
    return end,h

def fields(T,w,other=None):
    end,h=counter(T);start,_=counter(rev(T));full=(1<<len(T))-1
    F=sum(end(S,w)*h(full^S) for S in range(1,full+1) if S>>w&1)
    ans={'H':h(full),'F':F}
    if other is not None:
        # Two nonempty ordered slots: first ending w, second starting other.
        J=sum(end(S,w)*start(full^S,other) for S in range(1,full) if S>>w&1 and not(S>>other&1))
        ans['J_end_start']=J
    return ans

def polyadd(P,Q,scale=1):
    R=P.copy()
    for m,c in Q.items():
        R[m]=R.get(m,0)+scale*c
        if not R[m]:del R[m]
    return R

def mul(P,Q):
    R={}
    for a,c in P.items():
        for b,d in Q.items():
            if not a&b:R[a|b]=R.get(a|b,0)+c*d
    return {m:c for m,c in R.items() if c}

def matrix(T):
    n=len(T);G=[[{} for _ in T] for _ in T]
    def walk(s,v,mask):
        G[s][v][mask]=G[s][v].get(mask,0)+1
        u=T[v]&~mask
        while u:
            bit=u&-u;u-=bit;walk(s,bit.bit_length()-1,mask|bit)
    for s in range(n):walk(s,s,1<<s)
    return G

def update_birth(G,cut):
    n=len(G);left=[{} for _ in G];right=[{} for _ in G]
    for i in range(n):
        for j in range(n):
            if cut>>j&1:left[i]=polyadd(left[i],G[i][j])
            else:right[i]=polyadd(right[i],G[j][i])
    z={1<<n:1};R=[[{} for _ in range(n+1)] for _ in range(n+1)]
    for i in range(n):
        for j in range(n):R[i][j]=polyadd(G[i][j],mul(z,mul(left[i],right[j])))
        R[i][n]=mul(z,left[i]);R[n][i]=mul(z,right[i])
    R[n][n]=z
    return R

def flip(T,i,j):
    R=list(T);R[i]^=1<<j;R[j]^=1<<i;return tuple(R)
def update_flip(G,i,j):
    # Existing i->j becomes j->i.
    return [[polyadd(polyadd(G[a][b],mul(G[a][j],G[i][b])),mul(G[a][i],G[j][b]),-1) for b in range(len(G))] for a in range(len(G))]

def main():
    audit=json.loads((HERE/'input_audit.json').read_text());zero=json.loads((HERE/'input_zero_ready.json').read_text())
    report={'source_ref':'4ad839b20dd280cdcfbadf86cb62b2be3c9615ad','literal':[], 'zero':[], 'hostile':[]}
    for z in audit['states']:
        T=tuple(z['rows']);valid(T);_,h=counter(T);H=h((1<<len(T))-1)
        need(H==z['H'],'literal count')
        report['literal'].append({'stage':z['stage'],'H':H})
    for z in zero:
        T=tuple(z['suffix_rows']);w=z['suffix_names'].index(z['inport']);f=fields(T,w,z['suffix_names'].index(6))
        need((f['H'],f['F'])==(z['h_B'],z['F_w']),'suffix evaluator')
        need(f['H']+(1<<len(z['prefix_names']))*f['F']==z['H'],'zero clock')
        report['zero'].append(f)
    for T,expected in zip(audit['hostile']['rows'],audit['hostile']['responses']):
        _,h=counter(tuple(T));_,h2=counter(birth(tuple(T),256))
        base=h((1<<len(T))-1);out=h2((1<<(len(T)+1))-1)
        need(base==462915 and out==expected,'81/27 distinguished cut')
        report['hostile'].append({'H':base,'cut256_H':out})
    births=flips=0
    for n in range(1,5):
        edges=list(combinations(range(n),2))
        for code in range(1<<len(edges)):
            T=[0]*n
            for k,(i,j) in enumerate(edges):
                if code>>k&1:T[i]|=1<<j
                else:T[j]|=1<<i
            T=tuple(T);G=matrix(T)
            for cut in range(1<<n):
                need(update_birth(G,cut)==matrix(birth(T,cut)),'whole marked birth matrix');births+=1
            for i,j in edges:
                if not(T[i]>>j&1):i,j=j,i
                need(update_flip(G,i,j)==matrix(flip(T,i,j)),'whole marked flip matrix');flips+=1
    report['exhaustive']={'base_orders':[1,2,3,4],'birth_full_matrix_checks':births,'flip_full_matrix_checks':flips}
    # Distinct branch-word counts and rigorous finite-bank depth obstruction.
    seed=56935;k=2
    for depth in range(1,30):
        upper=3**depth*(seed+2)-2
        cap=k*(upper.bit_length()-1) # k floor(log2 upper)
        if 2**depth>cap:break
    report['bank_bound']={'seed':seed,'bank_size':k,'depth':depth,'distinct_ternary_values':2**depth,'upper_count':upper,'bank_capacity_upper':cap}
    vals={seed}
    for d in range(1,depth+1):
        vals={3*x+e for x in vals for e in (2,4)}
        need(len(vals)==2**d,'branch distinctness')
    report['status']='PASS exact arithmetic and full squarefree-matrix identities'
    (HERE/'receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
