"""Independent graph checks of the new constructors and source-preserving moves."""
from itertools import product,permutations
from pathlib import Path
import json
from mechanisms import *
HERE=Path(__file__).resolve().parent

def literal(rows):
    n=len(rows)
    return sum(all(rows[a]>>b&1 for a,b in zip(p,p[1:])) for p in permutations(range(n)))

def strong(rows):
    n=len(rows)
    for v in range(n):
        seen=1<<v
        while True:
            old=seen
            for i in range(n):
                if seen>>i&1:seen|=rows[i]
            if seen==old:break
        if seen!=(1<<n)-1:return False
    return True

def replay(record):
    rows=record['source'];n=len(rows);ports=record['ports']
    M=[[bool(rows[i]>>j&1) for j in range(n)] for i in range(n)]
    for i in range(n):M[i].append(not bool(record['birth_out']>>i&1))
    M.append([bool(record['birth_out']>>i&1) for i in range(n)]+[False])
    for path in record['paths']:
        need(sorted(path)==list(range(n+1)),'independent path permutation')
        need(all(M[a][b] for a,b in zip(path,path[1:])),'independent legal native move')
        fixed={frozenset((a,b)) for a,b in zip(path,path[1:])}
        for i in range(n+1):
            for j in range(i+1,n+1):
                if frozenset((i,j)) not in fixed:M[i][j],M[j][i]=M[j][i],M[i][j]
    result=tuple(sum(1<<j for j in range(n) if M[i][j]) for i in range(n))
    need(list(result)==record['target'],'independent exact endpoint')
    diff=[(i,j) for i in range(n) for j in range(i+1,n) if (rows[i]>>j&1)!=(result[i]>>j&1)]
    need(diff==[tuple(ports)],'exactly the named port edge changed')
    need(record['vertex_clock']==[n,n+1,n+1,n+1,n],'explicit birth/move/delete clock')

def main():
    checks=0;permchecks=0;resetchecks=0;strongcharts=0
    for m in range(6):
        for w in product((0,1),repeat=m):
            for u in product((0,1),repeat=m):
                counts=[]
                for forward in (True,False):
                    rows=source(w,u,forward);h=hp(rows)
                    need(h==gap_h(w,u,forward),'subset DP versus gap partitions')
                    if m<=3:
                        need(h==literal(rows),'independent permutation count');permchecks+=1
                    counts.append(h);checks+=1
                need(counts[0]-counts[1]==jump(w,u)==2*pair(wedge(m),(1,)+w,(1,)+u),'edge contraction and wedge')
                if m>=2 and w[0]==u[0]==1 and w[-1]==u[-1]==0:
                    need(strong(source(w,u)) and strong(source(w,u,False)),'strong endpoint guard')
                    replay(reset(w,u));resetchecks+=1;strongcharts+=1
    jumps=[]
    for d in list(range(1,257))+[1023,1024,65535,65536,10**20+1]:
        w,u=build_jump(d);n=len(w)+2
        need(n==d.bit_length()+6 and jump(w,u)==2*d,'logarithmic every-even-jump constructor')
        hp0=gap_h(w,u);hm=gap_h(w,u,False)
        need(hp0-hm==2*d and strong(source(w,u)) and strong(source(w,u,False)),'literal graph endpoints')
        record=reset(w,u);replay(record);resetchecks+=1
        if d<=32:
            need(hp(source(w,u))==hp0 and hp(source(w,u,False))==hm,'jump graph DP')
        if d<=16 or d>256:jumps.append({'d':d,'order':n,'plus':hp0,'minus':hm,'words':[w,u]})
    for P in range(1,33):
        for R in range(1,33):
            w,u=chart(P,R)
            need(jump(w,u)==2*((R+1)//2-(P+1)//2),'balanced chart formula')
    families=[]
    for off,ranges in ((-2,range(2,5)),(0,range(2,5)),(3,range(3,6))):
        for r in ranges:
            P=(1<<r)+off;w,u=chart(P,P+2)
            a,b=hp(source(w,u)),hp(source(w,u,False))
            need(a==gap_h(w,u) and b==gap_h(w,u,False) and a-b==2,'family source DP')
            rec=reset(w,u);replay(rec);resetchecks+=1
            families.append({'offset':off,'r':r,'order':len(w)+2,'upper':a,'lower':b,'program':rec})
    w,u=chart(1,3);v,z=chart(3,5)
    need(jump(w,u)==jump(v,z)==2 and gap_h(w,u)!=gap_h(v,z),'source sidecar is necessary')
    controls={'same_jump_different_counts':[[gap_h(w,u,False),gap_h(w,u)],[gap_h(v,z,False),gap_h(v,z)]],
              'zero_jump':jump(w,w),'reversed_sign':jump(u,w)}
    need(controls['zero_jump']==0 and controls['reversed_sign']==-2,'zero/sign controls')
    out={'status':'FINITE-EXACT','universe':'all ordered Boolean word pairs of equal length 0 through 5; both apex-edge orientations',
         'subset_DP_graphs':checks,'independent_permutation_graphs':permchecks,'strong_small_charts':strongcharts,
         'independent_native_replays':resetchecks,'jump_examples':jumps,'families':families,'controls':controls}
    (HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
    print('PASS',checks,'graph DP checks;',resetchecks,'independent native replays')
if __name__=='__main__':main()
