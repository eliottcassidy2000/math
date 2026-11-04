"""Finite source-specific returns with an explicit move of the live mark.

Alignment is only a search heuristic. Every realized edit is compiled into
literal Hamilton-guarded complements and independently replayed afterward.
"""
from pathlib import Path
from functools import lru_cache
import sys,json,random,subprocess,importlib.util
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('prior_named_returns',HERE.parent/'return18979_palindrome_20261103/returns.py')
prior=importlib.util.module_from_spec(spec);spec.loader.exec_module(prior)
guard,relabel,junction,need,flip=prior.guard,prior.relabel,prior.junction,prior.need,prior.flip

def observe(T,r,s,exe):
    v=list(map(int,subprocess.run([exe],input=' '.join(map(str,[len(T),r,s,*T]))+'\n',text=True,capture_output=True,check=True).stdout.split()))
    return [v[0],v[3],v[4],v[5]],[v[0],(v[6]-v[2])//2,(v[6]+v[2])//2]

def main():
    exe=sys.argv[1];targets=json.loads((HERE/'targets.json').read_text());out=[]
    sources=json.loads((HERE.parent/'return18979_palindrome_20261103/returns.json').read_text())
    for old in sources:
      U=tuple(old['next_rows']);oldnames=old['next_names'];oldroot=old['next_root'];H=old['next_H']
      for budget in (1,2):
        candidates=[t for t in targets if t['H']==H and t['budget']==budget]
        if not candidates:continue
        size=max(len(t['rows']) for t in candidates if len(t['rows'])<=len(U)+1);candidates=[t for t in candidates if len(t['rows'])==size]
        names=oldnames[:];S=U;permanent=None
        if size>len(U):
            S=junction.ear(U,1<<oldroot);r=len(U);names.append(max(names)+1)
            permanent=dict(name=names[-1],cut=1<<oldroot,rows=S)
        else:r=len(U)-1
        p=oldroot;n=len(S);spare=names.index(6)
        need([i for i in range(n) if S[i]>>r&1]==[p],'new live mark has old root as unique inport')
        rng=random.Random(H*10+budget);best=None;bestscore=10**9
        for z in candidates:
            Z=tuple(z['rows']);free=[i for i in range(n) if i not in (r,p)];values=[i for i in range(size) if i not in (z['root'],z['inport'])]+[-i-1 for i in range(n-size)]
            def score(m):
                if m[spare]<0:return 100000
                return sum(((S[i]>>j)^(Z[m[i]]>>m[j]))&1 for i in range(n) if m[i]>=0 for j in range(i+1,n) if m[j]>=0)
            for start in range(70):
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
                keep=[i for i in range(n) if m[i]>=0];T=relabel(S,keep)
                if d<bestscore and m[spare]>=0 and junction.strong(T):bestscore=d;best=(z,m[:],T,keep)
        need(best is not None,'strong named alignment');z,m,T,keep=best
        local=[names[i] for i in keep];rr=local.index(names[r]);sp=local.index(6);W=relabel(tuple(z['rows']),[m[i] for i in keep]);N=len(T)
        edits=[(i,j) for i in range(N) for j in range(i+1,N) if ((T[i]>>j)^(W[i]>>j))&1];startT=T
        print('align',H,budget,'arcs',len(edits),'deleted',n-size,'new mark',local[rr],flush=True)
        @lru_cache(None)
        def schedule(mask):
            V=startT
            for j,(a,c) in enumerate(edits):
                if mask>>j&1:V=flip(V,a,c)
            if not junction.strong(V):return None
            if mask==(1<<len(edits))-1:return ()
            for j,(a,c) in enumerate(edits):
                if mask>>j&1 or not junction.strong(flip(V,a,c)):continue
                if not V[a]>>c&1:a,c=c,a
                P=guard(V,a,c,rr)
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
        pr,g=observe(T,rr,sp,exe);need(T==W and pr[0]==H and pr[3]-H==budget,'independent returned profile')
        V=junction.ear(junction.ear(T,1<<rr),1<<N);Hnext=observe(V,N,sp,exe)[0][0];need(Hnext==3*H+2*budget and junction.strong(V),'literal next ternary branch')
        rec=dict(source_H=H,budget=budget,source=U,source_names=oldnames,source_root=oldroot,source_spare=old['next_spare'],permanent_birth=permanent,
            moved_mark_name=names[r],inport_name=names[p],deleted_old_names=[a for a in oldnames if a not in local],after_deletion=startT,
            names=local,root=rr,spare=sp,target=T,profile=pr,golden=g,target_grammar=z['grammar'],trace=trace,native_moves=2*len(trace),temporary_helpers=len(trace),peak_repair_order=N+1,
            next_rows=V,next_names=local+[fresh,fresh+1],next_root=N,next_spare=sp,next_H=Hnext)
        out.append(rec);(HERE/'moving_returns.json').write_text(json.dumps(out,indent=2)+'\n');print('RETURN',H,budget,pr,'next',Hnext,'moves',2*len(trace),flush=True)
if __name__=='__main__':main()
