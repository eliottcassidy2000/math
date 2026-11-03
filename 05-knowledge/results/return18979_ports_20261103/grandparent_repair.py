"""Exact grandparent-defect obstruction and a uniform two-arc escape."""
from itertools import combinations,permutations,product
from support import HERE,PREV,need,save,barrier_source,junction,flip,phases,observer
from matching_carry import clone,matching
from observer import records

def chain(T):
    p=tuple(sorted(range(len(T)),key=lambda i:-T[i].bit_count()))
    return p if all(T[a]>>b&1 for i,a in enumerate(p) for b in p[i+1:]) else None

def path(T):
    p=[]
    for v in range(len(T)):
        i=next((i for i,w in enumerate(p) if T[v]>>w&1),len(p));p.insert(i,v)
    need(all(T[a]>>b&1 for a,b in zip(p,p[1:])),'insertion path');return tuple(p)

def minimum_feedback(T):
    size=1<<len(T);dp=[0]*size;last=[-1]*size
    for S in range(1,size):
        opts=[(dp[S^(1<<v)]+(T[v]&(S^(1<<v))).bit_count(),v) for v in range(len(T)) if S>>v&1]
        dp[S],last[S]=min(opts)
    S=size-1;p=[]
    while S:v=last[S];p.append(v);S^=1<<v
    return dp[-1],tuple(reversed(p))

def framed_flip(T,a,b,clone_ids,names,fresh):
    if not T[a]>>b&1:a,b=b,a
    # All members of the clone block must be consecutive in this fixed frame.
    keep=[i for i in range(len(T)) if i not in (a,b) and i not in clone_ids[1:]]
    small=junction.sub(T,tuple(keep));P0=[keep[i] for i in path(small)]
    j=P0.index(clone_ids[0]);P0[j:j+1]=clone_ids
    n=len(T);z,s,w=n,n+1,n+2
    X=junction.ear(T,1<<b);X=junction.ear(X,(1<<(n+1))-1);X=junction.ear(X,((1<<(n+2))-1)^(1<<a))
    P=(a,b,z)+tuple(P0)+(s,w);Q=(a,w,b,z)+tuple(P0)+(s,)
    need(sorted(P)==list(range(n+3)) and sorted(Q)==list(range(n+3)),'literal full guards')
    need(all(X[x]>>y&1 for x,y in zip(P,P[1:])),'first native guard')
    Y=junction.native(X,P);need(all(Y[x]>>y&1 for x,y in zip(Q,Q[1:])),'second native guard')
    Z=junction.native(Y,Q);U=junction.sub(Z,tuple(range(n)))
    need(U==flip(T,a,b),'isolated named old arc')
    labels=list(names)+[f'helper{fresh+j}' for j in range(3)]
    return U,{'arc':[names[a],names[b]],'guards':[[labels[i] for i in R] for R in (P,Q)],'rows':U,
              'birth_incoming_masks':[1<<b,(1<<(n+1))-1,((1<<(n+2))-1)^(1<<a)],'deleted':labels[-3:]}

def main():
    T,r,sp,names=barrier_source();n=len(T);incoming=[i for i in range(n) if T[i]>>r&1]
    need(len(incoming)==1,'one in-port');w=incoming[0];ids=[i for i in range(n) if i not in (r,w)]
    D=junction.sub(T,tuple(ids));need(junction.H(D)==27,'defect is the grandparent H')
    C=junction.sub(T,tuple(i for i in range(n) if i!=r));cw=w-(w>r)
    h,s,q,F=junction.profile(T,r);ch,cs,cq,cF=junction.profile(C,cw)
    need((s,h,F)==(ch-cs,ch-cs+cF,ch+cF),'three deletion identities')
    cost,order=minimum_feedback(D);need(cost==2,'minimum grandparent arc repair')
    flips=[(ids[a],ids[b]) for i,a in enumerate(order) for b in order[i+1:] if D[b]>>a&1]
    need({frozenset((names[a],names[b])) for a,b in flips}=={frozenset((2,32)),frozenset((2,4))},'literal two-arc witness')
    untouched=[i for i in ids if i!=sp];U0=junction.sub(T,tuple(untouched));need(junction.H(U0)==15,'root/spare-invisible cyclic remainder')
    extra=[v for v in ids if v!=sp and chain(junction.sub(T,tuple(i for i in ids if i not in (sp,v)))) is not None]
    need([names[v] for v in extra]==[2],'unique one-extra-old-port repair budget')
    # Now allow the entire root star to change, not just a fixed in-port.
    def possible_inports(extra_ports):
        return [v for v in range(n) if v!=r and chain(junction.sub(T,tuple(i for i in range(n) if i not in {r,sp,v,*extra_ports}))) is not None]
    need(not possible_inports(()),'all mixed root/spare corrections retain an uneditable cycle')
    completion={names[p]:[names[v] for v in possible_inports((p,))] for p in range(n) if p not in (r,sp) and possible_inports((p,))}
    need(completion=={2:[21],21:[2]},'complete one-extra-port menu with moving in-port')
    # Positive alternative: activate old21 and move the unique in-port to old2.
    neww=names.index(2);active=(sp,w);fixed=[i for i in range(n) if i not in (r,neww,*active)]
    fixed=[fixed[i] for i in chain(junction.sub(T,tuple(fixed)))];candidates=[]
    for positions in combinations(range(n-2),2):
        for placed in permutations(active):
            it=iter(fixed);order2=[placed[positions.index(i)] if i in positions else next(it) for i in range(n-2)]
            for bits in product((0,1),repeat=2):
                V=T
                for i in range(n):
                    if i!=r and bool(V[i]>>r&1)!=(i==neww):V=flip(V,i,r)
                for j,a in enumerate(order2):
                    for c in order2[j+1:]:
                        if not V[a]>>c&1:V=flip(V,a,c)
                for a,bit in zip(active,bits):
                    if bool(V[neww]>>a&1)!=bool(bit):V=flip(V,neww,a)
                edits=[(a,c) for a,c in combinations(range(n),2) if (T[a]^V[a])>>c&1]
                need(all(a in (r,*active) or c in (r,*active) for a,c in edits),'alternate permitted support')
                if junction.strong(V):candidates.append((len(edits),V,edits,order2))
    altcost,Alt,altedits,altorder=min(candidates)
    Alt2,altids=clone(Alt,r,2);altprofile=junction.profile(Alt2,altids.index((r,0)))
    need(altprofile[3]-altprofile[0]==1,'positive alternate support completion')
    alternate={'extra_port':21,'new_inport':2,'minimum_arcs_in_declared_support':altcost,'strong_candidates':len(candidates),
      'rows':Alt,'flips':[[names[a],names[c]] for a,c in altedits],'grandparent_order':[names[i] for i in altorder],
      'base_profile':junction.profile(Alt,r),'two_clone_profile':altprofile,'root_phases':phases(Alt,r),'spare_phases':phases(Alt,sp)}
    mixed_witnesses=[]
    for v in range(n):
        if v==r:continue
        keep=[i for i in range(n) if i not in (r,sp,v)]
        tri=next(c for c in combinations(keep,3) if sum(bool(T[a]>>b&1) for a,b in ((c[0],c[1]),(c[1],c[2]),(c[2],c[0]))) in (0,3))
        mixed_witnesses.append({'prospective_inport':names[v],'untouched_directed_triangle':[names[i] for i in tri]})
    U=T
    for a,b in flips:U=flip(U,a,b)
    need(junction.strong(U) and junction.profile(U,r)==(1425,181,181,1425),'actual strong repaired source')
    need(chain(junction.sub(U,tuple(ids))) is not None,'grandparent made transitive')
    samples=[];programs=[]
    for m in range(1,5):
        A,addresses=clone(T,r,m);B,_=clone(U,r,m);rr=addresses.index((r,0))
        z=junction.profile(B,rr);za=junction.profile(A,rr)
        need(z[0]==1243*2**m-1061 and za[0]==5407*2**m-4489,'two matching-one count formulas')
        need(len(matching(A,rr)[0])==1 and len(matching(B,rr)[0])==1,'matching-one upper certificate')
        if m>=2:
            need(za[3]-za[0]==27 and z==(1243*2**m-1061,1243*2**(m-1)-1062,1243*2**(m-1)-1061,1243*2**m-1060),'uniform grandparent response')
        if m in (1,2,4):
            labels=[f'{names[v]}:{j}' if v==r else str(names[v]) for v,j in addresses];cur=A;trace=[]
            for k,(a,b) in enumerate(flips):
                cur,tr=framed_flip(cur,addresses.index((a,0)),addresses.index((b,0)),[addresses.index((r,j)) for j in range(m)],labels,100+3*k);trace.append(tr)
            need(cur==B,'uniform framed compiler endpoint')
            programs.append({'m':m,'trace':trace,'native_moves':4,'births':6,'deletions':6,'old_deletions':0,'peak_order':len(A)+3,'expanded_arc_reversals':4*((len(A)+2)*(len(A)+1)//2)})
        samples.append({'m':m,'original_profile':za,'repaired_profile':z})
        junction.endpoints.cache_clear()
    # Independent full small-tournament check of the algebraic deletion identity.
    checked=0
    for nn in range(2,5):
        pairs=list(combinations(range(nn),2))
        for bits in range(1<<len(pairs)):
            A=[0]*nn
            for k,(a,b) in enumerate(pairs):
                if bits>>k&1:A[a]|=1<<b
                else:A[b]|=1<<a
            A=tuple(A)
            for v in range(nn):
                incoming=[i for i in range(nn) if A[i]>>v&1]
                if len(incoming)!=1:continue
                ii=incoming[0];B,adr=clone(A,v,2);z=junction.profile(B,adr.index((v,0)))
                R=junction.sub(A,tuple(i for i in range(nn) if i not in (v,ii)))
                need(z[3]-z[0]==junction.H(R),'all small unique-in-port sources');checked+=1
                junction.endpoints.cache_clear()
    save('grandparent_repair.json',{'status':'EXACT_ALL_SIZE_MECHANISM_WITH_FINITE_CERTIFICATES','source_rows':T,'source_names':names,'root':r,'in_port':w,'spare':sp,
      'grandparent_rows':D,'grandparent_H':27,'root_spare_untouched_H':15,'unique_extra_port':2,'feedback_minimum':cost,
      'arbitrary_mixed_root_spare_obstruction':mixed_witnesses,'one_extra_port_completion_moving_inport':completion,'alternate_completion':alternate,
      'transitive_order':[names[ids[i]] for i in order],'flips':[[names[a],names[b]] for a,b in flips],'repaired_rows':U,'samples':samples,'programs':programs,
      'original_formula':'5407*2^m-4489','repaired_formula':'1243*2^m-1061','count_debt':'-4164*2^m+3428','uniform_ready_scope':'m>=2, first clone, defect1; next reset gives 3H+2',
      'phases':{str(names[v]):{'before':phases(T,v),'after':phases(U,v)} for v in (r,sp,names.index(2))},
      'two_port_before':records(observer(T,(r,sp))),'two_port_after':records(observer(U,(r,sp))),'extra_port_after':records(observer(U,(r,names.index(2)))),
      'small_identity_checks':checked,'scope':'root star and its unique in-port fixed; permitted corrections touch declared old ports; no neutral count return or complete ternary closure claimed'})
    print('PASS grandparent repair: defect27 -> 1; unique extra port2; two arcs; scalar debt explicit',flush=True)
if __name__=='__main__':main()
