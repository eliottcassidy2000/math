"""Actual named returns; old deletions are explicit and both live ports survive."""
from pathlib import Path
from functools import lru_cache
import sys,json,random
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'return6325_completion_20261103'))
from return6325 import junction,need,flip
from audit import hp,profile
def relabel(T,p):return tuple(sum(1<<j for j in range(len(p)) if T[p[i]]>>p[j]&1) for i in range(len(p)))
def guard(T,a,b,q):
    n=len(T);full=(1<<n)-1
    @lru_cache(None)
    def dfs(mask,last):
        if mask==full:return () if last not in (b,q) else None
        avail=T[last]&~mask if last>=0 else full
        if last==a:avail&=1<<b
        else:avail&=~(1<<b)
        while avail:
            bit=avail&-avail;avail-=bit;v=bit.bit_length()-1;tail=dfs(mask|bit,v)
            if tail is not None:return (v,)+tail
        return None
    return dfs(0,-1)
def main():
    targets=json.loads((HERE/'targets.json').read_text());out=[]
    for sourcefile in ('return6325_plus2.json','return6325.json'):
        old=json.loads((HERE.parent/'return6325_completion_20261103'/sourcefile).read_text())
        U=tuple(old['next_rows']);names=old['next_names'];r=old['next_root'];b=names.index(6);n=len(U);H=old['next_H'];p=next(i for i in range(n) if U[i]>>r&1)
        for budget in (1,2):
            candidates=[t for t in targets if t['H']==H and t['budget']==budget and len(t['rows'])<=n]
            # Prefer targets retaining the greatest number of old names.
            size=max(map(lambda t:len(t['rows']),candidates));candidates=[t for t in candidates if len(t['rows'])==size]
            rng=random.Random(H*10+budget);best=None;bestscore=10**9
            for tid,z in enumerate(candidates):
                Z=tuple(z['rows']);free=[i for i in range(n) if i not in (r,p)];values=[i for i in range(size) if i not in (z['root'],z['inport'])]+[-i-1 for i in range(n-size)]
                def score(m):
                    if m[b]<0:return 100000
                    return sum(((U[i]>>j)^(Z[m[i]]>>m[j]))&1 for i in range(n) if m[i]>=0 for j in range(i+1,n) if m[j]>=0)
                for start in range(100):
                    rng.shuffle(values);m=[None]*n;m[r]=z['root'];m[p]=z['inport']
                    for i,j in zip(free,values):m[i]=j
                    d=score(m)
                    while True:
                        low=d;moves=[]
                        for ix,i in enumerate(free):
                            for j in free[ix+1:]:
                                m[i],m[j]=m[j],m[i];v=score(m);m[i],m[j]=m[j],m[i]
                                if v<low:low=v;moves=[(i,j)]
                                elif v==low and v<d:moves.append((i,j))
                        if not moves:break
                        i,j=rng.choice(moves);m[i],m[j]=m[j],m[i];d=low
                    keep=[i for i in range(n) if m[i]>=0];S=relabel(U,keep)
                    if d<bestscore and m[b]>=0 and junction.strong(S):bestscore=d;best=(z,m[:],S,keep);print('alignment',H,budget,'arcs',d,'old deletions',n-size,flush=True)
            need(best is not None,'target alignment with strong deletion endpoint');z,m,T,keep=best;local=[names[i] for i in keep];rr=local.index(names[r]);sp=local.index(6);W=relabel(tuple(z['rows']),[m[i] for i in keep]);N=len(T)
            edits=[(i,j) for i in range(N) for j in range(i+1,N) if ((T[i]>>j)^(W[i]>>j))&1];startT=T
            @lru_cache(None)
            def schedule(mask):
                S=startT
                for j,(a,c) in enumerate(edits):
                    if mask>>j&1:S=flip(S,a,c)
                if not junction.strong(S):return None
                if mask==(1<<len(edits))-1:return ()
                for j,(a,c) in enumerate(edits):
                    if mask>>j&1 or not junction.strong(flip(S,a,c)):continue
                    if not S[a]>>c&1:a,c=c,a
                    P=guard(S,a,c,rr)
                    if P is None:continue
                    rest=schedule(mask|1<<j)
                    if rest is not None:return ((j,P),)+rest
                return None
            sched=schedule(0);need(sched is not None,'strong guarded repair ordering');trace=[];fresh=max(names)+1
            for j,P in sched:
                a,c=edits[j]
                if not T[a]>>c&1:a,c=c,a
                i=P.index(a);cut=((1<<N)-1)^(1<<a);X=junction.ear(T,cut);P1=(*P,N);P2=(*P[:i+1],N,*P[i+1:]);Y=junction.native(X,P1);V=junction.native(Y,P2);T2=tuple(x&((1<<N)-1) for x in V[:N])
                need(T2==flip(T,a,c),'literal isolated edge repair');ns=local+[fresh]
                trace.append(dict(arc=[local[a],local[c]],helper=fresh,cut=cut,paths=[[ns[v] for v in P1],[ns[v] for v in P2]],rows=T2));fresh+=1;T=T2
            pr,g=profile(T,rr);need(T==W and pr[0]==H and pr[3]-H==budget,'independent returned profile')
            V=junction.ear(junction.ear(T,1<<rr),1<<N);Hnext=hp(V);need(Hnext==3*H+2*budget and junction.strong(V),'literal next ternary branch')
            rec=dict(source_H=H,budget=budget,source=U,source_names=names,source_root=r,source_spare=b,deleted_old_names=[names[i] for i in range(n) if i not in keep],after_deletion=startT,
                names=local,root=rr,spare=sp,target=T,profile=pr,golden=g,target_grammar=z['grammar'],trace=trace,native_moves=2*len(trace),temporary_helpers=len(trace),peak_repair_order=N+1,
                next_rows=V,next_names=local+[fresh,fresh+1],next_root=N,next_spare=sp,next_H=Hnext)
            out.append(rec);(HERE/'returns.json').write_text(json.dumps(out,indent=2)+'\n');print('RETURN',H,budget,pr,'next',Hnext,'moves',2*len(trace),flush=True)
if __name__=='__main__':main()
