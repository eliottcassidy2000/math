#!/usr/bin/env python3
"""Independently reconstruct both families and audit their literal named cones."""
import json,hashlib,pathlib
from audit_checks import *
H=pathlib.Path(__file__).resolve().parent

def build(a,odd):
    m=a+4;r=m-1;n=2*m+int(odd)
    def u(i,j):
        if i==j:return False
        if i<a and j<a:return i<j
        if a<=i<a+3 and a<=j<a+3:return (j-i)%3==1
        return (i<a and a<=j<a+3) or (a<=i<a+3 and j==r) or (i==r and j<a)
    block=tuple(tuple(u(i,j) for j in range(m)) for i in range(m))
    def edge(i,j):
        if i==j:return False
        if i<m and j<m:return u(j,i)
        if m<=i<2*m and m<=j<2*m:return u(i-m,j-m)
        if i<2*m and j<2*m:
            if not odd and {i,j}=={r,m+r}:return i==m+r
            return i<m
        if i==2*m:return (j>=m and j!=m+r) or j==r
        return not edge(j,i)
    return tuple(tuple(edge(i,j) for j in range(n)) for i in range(n)),block

def check(a,odd):
    t,u=build(a,odd);m=a+4;n=len(t);r=m-1;w=2*m;shore=sum(1<<v for v in range(m,2*m));expected=[shore]+([shore+(1<<w)] if odd else [])
    need(strong(t) and strong(u),'initial strongness');lam,delta,shores,vals=cut_certificate(t);need((lam,delta,shores)==(1,2,expected),'minimum cut family')
    E=[tuple(range(a))+tuple(a+(i+k)%3 for k in range(3))+(r,) for i in range(3)]
    for e in E:legal(u,e)
    if a<=8:need(set(paths(u,end=r))==set(E),'rooted Hamilton orbit incomplete')
    for k in range(3):
        phi=[v if not a<=v<a+3 else a+(v-a+k)%3 for v in range(m)]
        need(phi[r]==r and all(u[i][j]==u[phi[i]][phi[j]] for i in range(m) for j in range(m)),'rooted rotation not automorphism')
    dec=[tuple(v+m for v in x)+((w,) if odd else ())+tuple(reversed(y)) for x in E for y in E]
    for p in dec:legal(t,p);need(anti_failure(t,p) is None,'decreasing path not full anti-automorphism')
    p=tuple(range(a+1,-1,-1))+tuple(range(m,2*m))+((w,) if odd else ())+(r,a+2)
    q=tuple(range(2*m-1,m-1,-1))+tuple(range(a))+(a,a+2,r,a+1)+((w,) if odd else ())
    terminal,moves=replay(t,[p,q],detail=True)
    need(cut(terminal,set(range(m,2*m)))==0,'terminal zero cut missing')
    mid=reverse(t,p);U2=tuple(tuple(mid[i][j] for j in range(m,2*m)) for i in range(m,2*m));V2=tuple(tuple(mid[i][j] for j in range(m)) for i in range(m))
    sq=[sum(d*d for d in degrees(U2)),sum(d*d for d in degrees(opposite(V2)))]
    need(abs(sq[0]-sq[1])==2,'degree-square anti-isomorphism obstruction')
    if odd:need(degrees(mid)[q[0]]+degrees(mid)[q[-1]]==n-4,'odd endpoint obstruction')
    z={'a':a,'n':n,'odd':odd,'rows':rows(t),'minimum_shores':shores,'lambda':lam,'delta':delta,'decreasing_words':dec,'moves':moves,'block_square_pair':sq,'terminal_core_ports':core(terminal),'exact_black_nonstrong_distance':2}
    if a<=2:
        pp=list(paths(t));black=[p for p in pp if anti_failure(t,p) is not None]
        for s in shores:
            cone=[p for p in pp if s>>p[0]&1 and not s>>p[-1]&1];need(set(cone)==set(dec),'complete cone mismatch')
        need(black and all(strong(reverse(t,p)) for p in black),'black one-step exit exists')
        direct=[(s,cut(t,{v for v in range(n) if s>>v&1})) for s in range(1,(1<<n)-1)]
        need(min(c for s,c in direct)==1 and [s for s,c in direct if c==1]==shores,'direct cuts mismatch')
        z['exhaustive']={'all_hamilton_words':len(pp),'black_hamilton_words':len(black),'all_black_neighbors_strong':True,'cuts':len(direct)}
    return z

def main():
    cases=[check(a,odd) for a in list(range(1,26))+[100] for odd in (False,True)]
    # Preserve exact producer first-member named starts, no canon relabelling.
    for odd,name in [(False,'normal.json'),(True,'odd-normal.json')]:
        src=json.loads((H/'producer'/name).read_text())
        for c in src['cases']:
            own=next(x for x in cases if x['a']==c['a'] and x['odd']==odd)
            need(own['rows']==c['initial_rows'],'producer source labels differ')
            for mine,theirs in zip(own['moves'],c['route']):
                need(list(mine['R_path'])==theirs['R_word'] and list(mine['K_path'])==theirs['K_word'],'producer word labels differ')
                need(mine['K_source_rows']==theirs['K_source'] and mine['K_target_rows']==theirs['K_target'],'producer native endpoint differs')
    # Fixed-minimum-shore failure is not global distance failure.
    t=from_rows((8,1,3,6));sh=1<<1;pp=list(paths(t));cone=[p for p in pp if sh>>p[0]&1 and not sh>>p[-1]&1]
    need(cone==[(1,0,3,2)] and anti_failure(t,cone[0]) is None,'fixed-shore hostile')
    exits=[p for p in pp if anti_failure(t,p) is not None and not strong(reverse(t,p))];need(exits,'fixed-shore confused with global obstruction')
    print(json.dumps({'status':'PASS','representation':'independent adjacency matrices; source built by pair predicate','scope':'a=1..25 and a=100; complete Hamilton/cut universes n10,11,12,13; structural all-order proof audited separately','cases':cases,'fixed_shore_n4':{'shore':sh,'blue_cone':cone,'black_other_shore_exit':exits[0]}},sort_keys=True,indent=2))
if __name__=='__main__':main()
