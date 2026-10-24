#!/usr/bin/env python3
"""Same-agent independent audits: BFS, literal owner decoding, and subset counts.
This is not represented as review by a separate agent.
"""
from itertools import permutations,combinations,product
from collections import deque
import json
import one_bit

CHECKS=0
def need(x,m):
    global CHECKS
    CHECKS+=1
    if not x:raise RuntimeError(m)
def bfs_audit():
    faces=[]
    for pins in product((-1,0,1),repeat=2):
        faces.append({x for x in range(4) if all(v==-1 or (x>>i)&1==v for i,v in enumerate(pins))})
    menus=0;states=0;boundaries=0
    for s in range(1,4):
        for E in product(faces,repeat=s):
            ends=[p for p in permutations(range(4),s) if all(x in e for x,e in zip(p,E))]
            if not ends:continue
            start=ends[0];seen={start};todo=deque([start])
            while todo:
                x=todo.popleft()
                for i in range(s):
                    for bit in range(3):
                        v=x[i]^(1<<bit)
                        if v in x or v%4 not in E[i]:continue
                        y=x[:i]+(v,)+x[i+1:]
                        if y not in seen:seen.add(y);todo.append(y)
            need(all(p in seen for p in ends),'independent full transition BFS covers all boundary assignments')
            menus+=1;states+=len(seen);boundaries+=len(ends)
    return {'menus':menus,'reachable_states':states,'boundary_states':boundaries}
def native_audit():
    code=29558785;B=[[0]*8 for _ in range(8)]
    for i,(u,v) in enumerate(combinations(range(8),2)):B[u][v]=(code>>(27-i))&1;B[v][u]=1-B[u][v]
    paths=((3,0,6,7,1,5,4,2),(3,1,6,7,0,5,4,2));pin=(0,2,3,5,6,7);g=(4,0,2,3,1,5,6,7)
    def partner(A,p):
        out=[row[:] for row in A];consecutive={(p[i],p[i+1]) for i in range(len(p)-1)}
        for u,v in combinations(range(len(A)),2):
            if (u,v) not in consecutive and (v,u) not in consecutive:out[u][v],out[v][u]=out[v][u],out[u][v]
        return out
    reports=[]
    for r in range(1,6):
        d=r+1;n=8*d;N=1<<r
        A=[[0 if u==v else B[u%8][v%8] if u//8==v//8 else int(u//8<v//8) for v in range(n)] for u in range(n)]
        Ps=[tuple(8*i+v for i in range(d) for v in paths[(x>>i)&1]) for x in range(1<<d)]
        Ks=[];resources=[];closed=set();chords=[(i,j) for i,j in combinations(range(n),2) if j>i+1];idx={z:i for i,z in enumerate(chords)};full=(1<<len(chords))-1
        for p in Ps:
            need(len(set(p))==n and all(A[u][v] for u,v in zip(p,p[1:])),'literal full Hamiltonian path')
            K=partner(A,p);Ks.append(K)
            need(partner(K,p)==A,'inverse path square')
            C=list(map(list,zip(*A)))
            need(partner(C,p[::-1])==list(map(list,zip(*K))),'converse with exact reversed integer clocks')
            w=sum(A[p[i]][p[j]]<<k for k,(i,j) in enumerate(chords))
            rev=sum(((w>>k)&1)<<idx[n-1-j,n-1-i] for k,(i,j) in enumerate(chords))
            orb={w,w^full,rev,rev^full}
            need(len(orb)==4 and not closed&orb,'normalized full-orbit resources are actual distinct capacity')
            closed.update(orb);resources.append(min(orb))
        cube={x:{x^(1<<j) for j in range(r)} for x in range(N)}
        h,e,L=one_bit.complete([cube]*N,tuple(range(N)),tuple(x^(N-1) for x in range(N)))
        for old,new in zip(h,h[1:]):
            changed=[i for i,(x,y) in enumerate(zip(old,new)) if x!=y]
            need(len(changed)==1,'one literal slot per clock')
            words=[x+N*l for x,l in new]
            need(len({resources[w] for w in words})==N,'no reused normalized port')
            i=changed[0];x=old[i][0]+N*old[i][1];y=new[i][0]+N*new[i][1]
            diff=x^y;need(diff.bit_count()==1,'one actual block bit')
            block=diff.bit_length()-1
            need([j for j in range(n) if Ps[x][j]!=Ps[y][j]]==[8*block+1,8*block+4],'two exact owner clocks')
            need(all(Ps[y][8*t+j]==8*t+paths[0][j] for t in range(d) for j in pin),'every original pinned owner retained')
            # g carries the zero-bit partner to the one-bit partner; reverse events need g^-1.
            transport=g if not (x>>block)&1 else tuple(g.index(j) for j in range(8))
            t=tuple(8*block+transport[u%8] if u//8==block else u for u in range(n))
            need(all(Ks[x][u][v]==Ks[y][t[u]][t[v]] for u,v in combinations(range(n),2)),'literal target transporter for each primitive')
        reports.append({'base_dimension':r,'enlarged_owners':n,'slots':N,'integer_events':len(e),'closed_reservoir_ports':len(closed)})
    return reports
def main():
    bfs=bfs_audit();native=native_audit()
    return {'status':'PASS; same-agent independent computational paths','independent_BFS':bfs,'literal_native_routes':native,'checks':CHECKS}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
