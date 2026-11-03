"""Induced regular-five obstruction and a source-changing apex readiness compiler."""
from itertools import combinations,permutations,product
from collections import Counter
import json
from effects import HERE,old,observer,need,junction
from observer import records

def cycles(T):
    return sum(all(T[a]>>b&1 for a,b in zip(p,p[1:]+p[:1])) for tail in permutations(range(1,len(T))) for p in [(0,)+tail])
def regular5(T):
    return [S for S in combinations(range(len(T)),5) if all(sum(T[a]>>b&1 for b in S if b!=a)==2 for a in S)]

def exceptional(tail):
    word=[(1,1)]+list(tail);n=len(word)+3;rows=[0]*n
    def arc(a,b):rows[a]|=1<<b
    arc(1,2);arc(1,0);arc(2,0)
    for i,(A,B) in enumerate(word,3):
        arc(0,i)
        for port,bit in ((1,A),(2,B)):arc(i,port) if bit else arc(port,i)
        for j in range(3,i):arc(j,i)
    return tuple(rows)
def normalize(w,u,forward=True):
    T=old.source(w,u,forward);m=len(w);a,b=(m,m+1) if forward else (m+1,m)
    selected=w if forward else u;trace=[]
    need(w[0]==u[0]==1 and w[-1]==u[-1]==0,'outer apex guards')
    for i in range(m-1):
        if selected[i]:continue
        P=(i,a,b)+tuple(j for j in range(m) if j!=i)
        need(all(T[x]>>y&1 for x,y in zip(P,P[1:])),'uniform traced Hamilton path')
        n=len(T);incoming=((1<<n)-1)^(1<<i);X=junction.ear(T,incoming)
        p1=P+(n,);p2=P[:1]+(n,)+P[1:]
        V=junction.native(junction.native(X,p1),p2)
        V=tuple(row&((1<<n)-1) for row in V[:n]);expected=list(T);expected[i]^=1<<a;expected[a]^=1<<i
        need(V==tuple(expected) and junction.strong(V),'literal star flip and strong endpoint')
        trace.append({'arc':[i,a],'ports':[a,b],'temporary_incoming':incoming,'paths':[p1,p2]});T=V
    profile=junction.profile(T,a)
    need(profile[3]-profile[0]==1,'actual one-ready root')
    return T,a,trace,profile

def main():
    R5=tuple(sum(1<<((i+d)%5) for d in (1,2)) for i in range(5))
    need(cycles(R5)==2 and all(junction.profile(R5,r)[1]==3 for r in range(5)),'two cycles but three starts: the stronger rooted obstruction')
    regular_count=0
    for T in junction.all_tournaments(5):
        if all(row.bit_count()==2 for row in T):
            need(all(junction.profile(T,r)[1]==3 for r in range(5)),'every labelled regular-five tournament has three starts per vertex')
            regular_count+=1
    need(regular_count==24,'complete regular-five universe')
    need(all(junction.H(junction.sub(R5,tuple(v for v in range(5) if v!=r)))>1 for r in range(5)),'every deletion nontransitive')
    exceptional_controls=[]
    for tail in ((),((1,0),),((0,1),),((1,0),(0,1)),((0,1),(1,0))):
        T=exceptional(tail);need(not regular5(T),'five collapsed exceptional ready cores are R5-free')
        exceptional_controls.append({'tail':tail,'rows':T,'profile':junction.profile(T,0)})
    controls=[];w,u=old.chart(3,4)
    for forward in (True,False):
        T=old.source(w,u,forward);S,r,trace,z=normalize(w,u,forward)
        V=observer(S,(7,8));slot=(7,8).index(r)
        e=sum(c for (sk,degree),c in V.items() if len(sk)+degree==1 and any(b==slot for a,b in sk))
        u2=sum(c for (sk,degree),c in V.items() if len(sk)+degree in (1,2) and any(a==slot for a,b in sk))
        controls.append({'forward':forward,'source_rows':T,'regular5_copies':regular5(T),'target_rows':S,'root':r,'native_trace':trace,
                         'profile':z,'target_regular5_copies':regular5(S),'ternary_child_H':z[0]+2*z[3],
                         'retained_ports':[7,8],'golden':[z[0],(u2-e)//2,(u2+e)//2],'full_observer':records(V)})
    tables=[]
    for m in range(2,10):
        vals=[]
        for mid in product((0,1),repeat=m-2):
            u=(1,)+mid+(0,);w=(1,)*(m-1)+(0,);T=old.source(w,u)
            h=old.gap_h(w,u);need(h==old.hp(T) if m<=5 else True,'independent small ready counts')
            vals.append(h)
        tables.append({'m':m,'H':vals,'distinct':len(set(vals)),'minimum':min(vals),'maximum':max(vals)})
        print('m',m,'counts',len(vals),'distinct',len(set(vals)),'range',min(vals),max(vals),'head',vals[:12],flush=True)
    normalization_checks=0
    for m in range(2,6):
        words=[(1,)+mid+(0,) for mid in product((0,1),repeat=m-2)]
        for w,u in product(words,repeat=2):
            for forward in (True,False):
                T,r,trace,z=normalize(w,u,forward)
                need(not regular5(T),'star normalization removes every regular-five witness')
                normalization_checks+=1
    (HERE/'ready.json').write_text(json.dumps({'status':'FINITE_EXACT','regular5':R5,'R5_Hamilton_cycles':2,'R5_start_counts':[3]*5,
        'collapsed_exceptional_ready_controls':exceptional_controls,'all_small_star_repairs':normalization_checks,'H497_controls':controls,'ready_tables':tables},indent=2)+'\n',encoding='utf8',newline='\n')
    print('H497 repairs',[(r['profile'],r['ternary_child_H']) for r in controls])
if __name__=='__main__':main()
