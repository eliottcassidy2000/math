"""An actual same-root return at 6325, with all old names retained."""
from pathlib import Path
import sys,json,random
from functools import lru_cache
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'digit_cycle_ternary_20261103'))
from support import need,junction,phases,values,observer,paths,flip
from observer import records

def insertion(T,vertices):
    p=[]
    for v in vertices:
        i=0
        while i<len(p) and not (T[v]>>p[i]&1):i+=1
        p.insert(i,v)
    need(all(T[a]>>b&1 for a,b in zip(p,p[1:])),'insertion Hamilton path')
    return p
def relabel(T,p):return tuple(sum(1<<j for j in range(len(T)) if T[p[i]]>>p[j]&1) for i in range(len(T)))
def main():
    inherited=json.loads((HERE.parent/'digit_cycle_ternary_20261103/return2107.json').read_text())
    U=tuple(inherited['next_rows']);names=inherited['next_names'];r=inherited['next_root'];spare=names.index(6);n=len(U)
    a=next(i for i in range(n) if U[i]>>r&1)
    rec=json.loads((HERE/'ready_search.json').read_text())['hits'][0]
    need(rec['parent_order']==n-1 and rec['budget']==2,'same-order target')
    Z=junction.ear(tuple(rec['rows']),1);rr=n-1;aa=0
    slots=[i for i in range(n) if i not in (r,a)];other=[i for i in range(n) if i not in (rr,aa)]
    rng=random.Random(6325);best=None;distance=10**9
    def score(p):return sum(((U[i]>>j)^(Z[p[i]]>>p[j]))&1 for i in range(n) for j in range(i+1,n))
    for restart in range(180):
        rng.shuffle(other);p=[None]*n;p[r]=rr;p[a]=aa
        for i,j in zip(slots,other):p[i]=j
        d=score(p)
        while True:
            moves=[];low=d
            for ix,i in enumerate(slots):
              for j in slots[ix+1:]:
                p[i],p[j]=p[j],p[i];z=score(p);p[i],p[j]=p[j],p[i]
                if z<low:low=z;moves=[(i,j)]
                elif z==low and z<d:moves.append((i,j))
            if not moves:break
            i,j=rng.choice(moves);p[i],p[j]=p[j],p[i];d=low
        if d<distance:distance=d;best=p[:];print('candidate arc budget',distance,flush=True)
    W=relabel(Z,best);edits=[(i,j) for i in range(n) for j in range(i+1,n) if ((U[i]>>j)^(W[i]>>j))&1]
    need(all(r not in e for e in edits),'entire old root star retained')
    @lru_cache(None)
    def seek(mask):
        T=U
        for j,(a,b) in enumerate(edits):
            if mask>>j&1:T=flip(T,a,b)
        if not junction.strong(T):return None
        if mask==(1<<len(edits))-1:return ()
        for j,(a,b) in enumerate(edits):
            if mask>>j&1 or not junction.strong(flip(T,a,b)):continue
            if not T[a]>>b&1:a,b=b,a
            P=next((P for P in paths(T) if P[-1]!=r and any(P[i:i+2]==(a,b) for i in range(n-2))),None)
            if P is None:continue
            tail=seek(mask|1<<j)
            if tail is not None:return ((j,P),)+tail
        return None
    schedule=seek(0);need(schedule is not None,'strong one-helper schedule exists')
    T=U;trace=[];fresh=max(names)+1
    for edit_id,P0 in schedule:
        i,j=edits[edit_id]
        a,b=(i,j) if T[i]>>j&1 else (j,i)
        position=P0.index(a);cut=((1<<n)-1)^(1<<a)
        X=junction.ear(T,cut);P=(*P0,n);Q=(*P0[:position+1],n,*P0[position+1:])
        need(all(X[v]>>w&1 for v,w in zip(P,P[1:])),'one-helper first legal guard')
        Y=junction.native(X,P);need(all(Y[v]>>w&1 for v,w in zip(Q,Q[1:])),'one-helper second legal guard')
        V=junction.native(Y,Q);V=tuple(x&((1<<n)-1) for x in V[:n])
        expect=list(T);expect[i]^=1<<j;expect[j]^=1<<i
        need(V==tuple(expect),'exact single old-edge repair')
        def frame(path):
            j=path.index(r);return [path[j-1] if j else None,path[j+1] if j+1<len(path) else None]
        need(frame(P)==frame(Q),'common retained clone frame')
        localnames=names+[fresh]
        trace.append(dict(arc=[names[a],names[b]],helpers=[fresh],cuts=[cut],
            paths=[[localnames[v] for v in p] for p in (P,Q)],rows=V,root_frame=[localnames[v] if v is not None else None for v in frame(P)],
            H_after=junction.H(V)))
        T=V;fresh+=1
    need(T==W and junction.strong(T),'actual strong target')
    profile=junction.profile(T,r);need(profile==(6325,891,893,6327),'literal ready profile')
    need(values(T,r)[1]==profile,'independent endpoint root counts')
    V=junction.ear(junction.ear(T,1<<r),1<<n);need(junction.H(V)==18979,'actual next offset-four child')
    out=dict(status='FINITE_EXACT_SUPPLIED_SOURCE_RETURN_NOT_UNIFORM_CLOSURE',source=U,names=names,root=r,spare=spare,
        target=T,target_identification=best,arc_repairs=distance,search_scope='180 deterministic steepest-swap starts, no minimum claim',
        trace=trace,native_moves=2*len(trace),temporary_births=len(trace),temporary_deletions=len(trace),old_deletions=0,peak_order=n+1,
        profile=profile,source_phases=phases(U,r),target_phases=phases(T,r),spare_phases=phases(T,spare),two_port_observer=records(observer(T,(r,spare))),
        next_rows=V,next_names=names+[fresh,fresh+1],next_root=n,next_H=18979,next_phases=phases(V,n))
    (HERE/'return6325.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PASS6325',profile,'native',out['native_moves'],'peak',out['peak_order'],'next18979',flush=True)
if __name__=='__main__':main()
