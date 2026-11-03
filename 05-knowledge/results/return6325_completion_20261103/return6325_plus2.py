"""Retain the same 6325 source, add one permanent vertex, return to 1-ready."""
from return6325 import *
def main():
    inherited=json.loads((HERE.parent/'digit_cycle_ternary_20261103/return2107.json').read_text())
    U=tuple(inherited['next_rows']);oldnames=inherited['next_names'];r=inherited['next_root'];spare=oldnames.index(6);oldn=len(U);n=oldn+1
    a=next(i for i in range(oldn) if U[i]>>r&1);names=oldnames+[max(oldnames)+1]
    rng=random.Random(63252);best=None;distance=10**9;best_id=None
    candidates=json.loads((HERE/'one_ready13.json').read_text())['hits']
    slots=[i for i in range(n) if i not in (r,a)];other=[i for i in range(n) if i not in (n-1,0)]
    for target_id,rec in enumerate(candidates):
        Z=junction.ear(tuple(rec['rows']),1)
        def score(p):return sum(((U[i]>>j)^(Z[p[i]]>>p[j]))&1 for i in range(oldn) for j in range(i+1,oldn))
        for restart in range(120):
            rng.shuffle(other);p=[None]*n;p[r]=n-1;p[a]=0
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
            if d<distance:
                distance=d;best=p[:];best_id=target_id;W=relabel(Z,best);print('plus2 arc budget',distance,flush=True)
    cut=sum(1<<i for i in range(oldn) if W[i]>>oldn&1);S=junction.ear(U,cut)
    edits=[(i,j) for i in range(oldn) for j in range(i+1,oldn) if ((U[i]>>j)^(W[i]>>j))&1]
    need(all(r not in e for e in edits),'retained old root star')
    @lru_cache(None)
    def seek(mask):
        T=S
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
    schedule=seek(0);need(schedule is not None,'strong one-helper schedule for plus2')
    T=S;trace=[];fresh=max(names)+1
    for j,P0 in schedule:
        a,b=edits[j]
        if not T[a]>>b&1:a,b=b,a
        i=P0.index(a);cc=((1<<n)-1)^(1<<a);X=junction.ear(T,cc);P=(*P0,n);Q=(*P0[:i+1],n,*P0[i+1:])
        need(all(X[v]>>w&1 for v,w in zip(P,P[1:])),'first guard')
        Y=junction.native(X,P);need(all(Y[v]>>w&1 for v,w in zip(Q,Q[1:])),'second guard')
        Z=junction.native(Y,Q);V=tuple(z&((1<<n)-1) for z in Z[:n]);need(V==flip(T,a,b) and junction.strong(V),'single strong edge repair')
        def frame(p):
            i=p.index(r);return [p[i-1] if i else None,p[i+1] if i+1<len(p) else None]
        need(frame(P)==frame(Q),'same actual root clone frame');local=names+[fresh]
        trace.append(dict(arc=[names[a],names[b]],helpers=[fresh],cuts=[cc],paths=[[local[v] for v in p] for p in (P,Q)],rows=V,H_after=junction.H(V)))
        T=V;fresh+=1
    need(T==W,'actual target');profile=junction.profile(T,r)
    need(profile[0]==6325 and profile[3]==6326 and values(T,r)[1]==profile,'independent ready1 profile')
    V=junction.ear(junction.ear(T,1<<r),1<<n);need(junction.H(V)==18977,'next offset-two child')
    out=dict(status='FINITE_EXACT_SUPPLIED_SOURCE_RETURN',source=U,old_names=oldnames,names=names,root=r,spare=spare,permanent_birth=dict(name=names[-1],cut=cut,rows=S,H=junction.H(S)),
        target=T,target_id=best_id,target_identification=best,arc_repairs=distance,search_scope='120 steepest-swap starts on each of two ready targets; no arc minimum claim',
        trace=trace,native_moves=2*len(trace),temporary_births=len(trace),temporary_deletions=len(trace),old_deletions=0,peak_order=n+1,profile=profile,
        source_phases=phases(U,r),target_phases=phases(T,r),spare_phases=phases(T,spare),two_port_observer=records(observer(T,(r,spare))),
        next_rows=V,next_names=names+[fresh,fresh+1],next_root=n,next_H=18977,next_phases=phases(V,n))
    (HERE/'return6325_plus2.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PASS6325 plus2',profile,'native',out['native_moves'],'peak',out['peak_order'],'next18977',flush=True)
if __name__=='__main__':main()
