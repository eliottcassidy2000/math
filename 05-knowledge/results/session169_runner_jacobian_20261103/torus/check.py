#!/usr/bin/env python3
"""Exact controls for flat signed-unit covers; standard library only."""
from collections import Counter, defaultdict, deque
from fractions import Fraction
from functools import lru_cache
from itertools import permutations, product
from math import gcd, isqrt, lcm
from pathlib import Path
import json


def need(ok, label):
    if not ok:
        raise RuntimeError(label)


@lru_cache(None)
def sign(n, a, b):
    if not n:
        return 0
    h = 1 << (n-1)
    u, r, v, t = a % h, a // h, b % h, b // h
    if not r and not t:
        return sign(n-1, u, v)
    if not r and t:
        return sign(n-1, v, u)
    if r and not t:
        return sign(n-1, u, v) ^ bool(v)
    return sign(n-1, v, u) ^ bool(v) ^ 1


def compose(a, b):
    return tuple(a[b[x]] for x in range(len(a)))


def power(p, k):
    r = tuple(range(len(p)))
    for _ in range(k):
        r = compose(p, r)
    return r


def full_group(n):
    N=1 << n
    gens=[tuple((x^a)+N*(s^sign(n,a,x)) for s in range(2) for x in range(N)) for a in range(1,N)]
    identity=tuple(range(2*N)); found={identity}; queue=deque([identity])
    while queue:
        p=queue.popleft()
        for g in gens:
            q=compose(g,p)
            if q not in found:
                found.add(q);queue.append(q)
    return found


def algebra_controls():
    full3=full_group(3); need(len(full3)==128,'full octonion group')
    reports=[]
    for n in range(3,8):
        h=1<<(n-1); N=2*h; points=tuple(range(2*N)); identity=points
        def decode(v):
            s,x=divmod(v,N); t,y=divmod(x,h)
            return s,y,t
        def enc(s,y,t):return s*N+y+h*t
        Q=tuple(enc(s^(t and bool(y)),y,t) for s,y,t in map(decode,points))
        Z=tuple(v^N for v in points)
        H=tuple(enc(s^bool(y)^t,y,t^1) for s,y,t in map(decode,points))
        J=tuple(enc(s^bool(y),y,t^1) for s,y,t in map(decode,points))
        Bs=[]
        for u in range(h):
            B=tuple(enc(s^(t and (bool(y)^bool(y^u))),y^u,t) for s,y,t in map(decode,points))
            Bs.append(B)
            need(sum(sign(n-1,u,y) for y in range(h))%2==0,'even lower sign row')
            if u:
                D=tuple(enc(s^sign(n-1,u,y)^t,y,t) for s,y,t in map(decode,points))
                L=tuple(enc(s^sign(n,u,y+h*t),y^u,t) for s,y,t in map(decode,points))
                need(compose(L,D)==B,'B belongs to kernel-left coset')
            translated=tuple(enc(s,y^u,t) for s,y,t in map(decode,points))
            need(compose(Q,compose(B,Q))==translated,'B coordinate conjugacy')
            need(power(B,2)==identity,'B involution')
            need(compose(B,H)==compose(H,B),'B H commute')
            need(compose(B,J)==compose(J,B),'B J commute')
        need(H==tuple((x^h)+N*(s^sign(n,h,x)) for s in range(2) for x in range(N)),'H actual left multiplication')
        need(power(H,2)==Z,'H square sign')
        need(power(J,2)==identity,'J involution')
        need(compose(Q,Q)==identity,'Q involution')
        for u in range(h):
            for v in range(h):
                need(compose(Bs[u],Bs[v])==Bs[u^v],'all B composition')
        AA={compose(B,power(H,k)) for B in Bs for k in range(4)}
        EE={compose(power(Z,z),compose(B,power(J,t))) for B in Bs for z,t in product(range(2),repeat=2)}
        for name,GG in [('A',AA),('E',EE)]:
            need(len(GG)==2*N,'group size '+name)
            need({p[0] for p in GG}==set(points),'regular orbit '+name)
            need(all(p==identity or all(p[x]!=x for x in points) for p in GG),'semiregular '+name)
        if n==3:
            need(AA<=full3 and EE<=full3,'enumerated full-group membership')
        ordersA=Counter(next(k for k in (1,2,4) if power(p,k)==identity) for p in AA)
        ordersE=Counter(next(k for k in (1,2) if power(p,k)==identity) for p in EE)
        # Characters of C2^(n-1) x C4 factor through k modulo two.
        a_char_values_at_Z=[(((mask>>(n-1))&1)*(2%2))%2 for mask in range(1<<n)]
        e_char_values_at_Z=[mask&1 for mask in range(1<<(n+1))]
        need(set(a_char_values_at_Z)=={0},'A binary quotients kill Z')
        need(sum(e_char_values_at_Z)==1<<n,'E binary quotients detect Z')
        reports.append({'level':n,'states':2*N,'A_order_counts':dict(sorted(ordersA.items())),'E_order_counts':dict(sorted(ordersE.items()))})
    GG=full_group(2); need(len(GG)==8,'quaternion order')
    need(any(compose(p,q)!=compose(q,p) for p in GG for q in GG),'quaternion nonabelian')
    need(len({p[0] for p in GG})==8,'quaternion regular')
    return {'levels':reports,'full_octonion_group':len(full3),'quaternion_group_order':len(GG)}


def honey_step(vertex,c):
    u,v,e=vertex
    return u+((1 if e else -1) if c=='b' else 0),v+((1 if e else -1) if c=='c' else 0),e^1


def reg_step(z,c,e):
    s,y,t=z
    if c=='b':return s^t,y^1,t
    if c=='c':return s^y^t^(0 if e else 1),y,t^1
    return z


def register_label(z):
    s,y,t=z
    return ('-' if s else '+')+'e_'+str(y+4*t)


def position_register(u,v):
    y=u%2;k=v%4;t=k%2;q=k//2
    return q^(t*y),y,t


def cycle_family(m):
    a,b,c=(2*m+1)**2,(2*m+3)**2,(2*m+5)**2
    p,q=b-a,c-a;g=gcd(p,q);r,s=q//g,p//g;D=a+b-c+g;L=D//2
    need((2*L-a)%g!=0,'disjoint arithmetic layer residues')
    xs=[L+g*((j*s)%r) for j in range(r)]
    ys=[a-x if x<L+q-p else c-x for x in xs]
    vertices=[z for pair in zip(xs,ys) for z in pair]
    colors=[]
    for i,x in enumerate(vertices):
        target=x+vertices[(i+1)%len(vertices)]
        need(target in (a,b,c),'correct target')
        need(isqrt(target)**2==target,'square sum')
        colors.append('abc'[(a,b,c).index(target)])
    need(len(set(vertices))==len(vertices),'simple cycle')
    need(min(vertices)>0,'positive cycle')
    need(max(vertices)==2*m*m+14*m+13,'cycle height')
    cur=(0,0,0);z=(0,0,0)
    for color in colors:
        z=reg_step(z,color,cur[2]);cur=honey_step(cur,color)
        need(z==position_register(cur[0],cur[1]),'position unit transport')
    need(cur==(r,-s,0),'primitive winding')
    ell=4 if s%2 else 2
    lapstates=[(0,0,0)];z=(0,0,0)
    for lap in range(1,ell+1):
        for i,color in enumerate(colors):z=reg_step(z,color,i%2)
        lapstates.append(z)
        need((z==(0,0,0))==(lap==ell),'exact lap return')
    need(ell==(2 if m%2 else 4),'parity family')
    return {'m':m,'length':len(vertices),'height':max(vertices),'winding':[r,-s],'return_multiplier':ell,'period':ell*len(vertices),'lap_states':[register_label(z) for z in lapstates]},vertices,colors


def honey_controls():
    vertices=list(product(range(2),range(4),range(2)))
    graph={x:[] for x in vertices};edges=[]
    for u,v,e in vertices:
        for color in 'abc':
            dest=honey_step((u,v,e),color);dest=(dest[0]%2,dest[1]%4,dest[2])
            graph[u,v,e].append(dest)
            z=reg_step(position_register(u,v),color,e)
            need(z==position_register(dest[0],dest[1]),'finite honey transport')
            if not e:edges.append(((u,v,e),dest,color))
    seen={vertices[0]};queue=deque(seen)
    while queue:
        u=queue.popleft()
        for v in graph[u]:
            if v not in seen:seen.add(v);queue.append(v)
    need(len(seen)==16 and len(edges)==24,'honey size connected')
    for start in vertices:
        cur=start
        for c in 'abcabc':
            nxt=honey_step(cur,c);cur=(nxt[0]%2,nxt[1]%4,nxt[2])
        need(cur==start,'hexagonal cells fill')
    # Hostile unmarked example: a triangle changes the retained bipartite layer.
    cur=(0,0,0);labels=[6];layers=[0]
    for color in 'acbacb':
        cur=honey_step(cur,color)
        u,v,e=cur;value=6+11*u+24*v
        labels.append(25-value if e else value);layers.append(e)
    need(labels==[6,19,30,6,19,30,6],'unmarked arithmetic triangle')
    need(layers[3]==1 and cur==(0,0,0),'triangle returns only after two marked laps')
    fam=[cycle_family(m)[0] for m in range(3,101)]
    return {'honeycomb_vertices':16,'honeycomb_edges':24,'family_members':len(fam),'unmarked_triangle_control':{'targets':[25,36,49],'labels':labels,'layers':layers},'examples':[f for f in fam if f['m'] in (3,4,5,6,100)]}


def physical_q1():
    # Direct independent order-six tournament enumeration. Labels are retained
    # only as canonical mask names, never read from an external result file.
    n=6;pairs=list((i,j) for i in range(n) for j in range(i+1,n));idx={e:i for i,e in enumerate(pairs)}
    chords=[(i,j) for i in range(n) for j in range(i+2,n)];ci={e:i for i,e in enumerate(chords)}
    maps=[]
    for p in permutations(range(n)):
        source=[idx[tuple(sorted((p[i],p[j])))] for i,j in pairs]
        rev=sum((p[i]>p[j])<<k for k,(i,j) in enumerate(pairs));maps.append((source,rev))
    cache={}
    def canon(mask):
        if mask not in cache:
            orbit=set()
            for source,rev in maps:
                value=rev
                for k,j in enumerate(source):value^=((mask>>j)&1)<<k
                orbit.add(value)
            c=min(orbit)
            cache.update((v,c) for v in orbit)
        return cache[mask]
    C=[]
    for t in range(1<<len(chords)):
        full=sum(((t>>k)&1)<<idx[e] for k,e in enumerate(chords))
        C.append(canon(full))
    classes=sorted(set(C));op={c:canon(c^((1<<len(pairs))-1)) for c in classes}
    free=sorted({min(c,op[c]) for c in classes if c!=op[c]})
    sidx=[ci[n-1-j,n-1-i] for i,j in chords]
    S=[sum(((t>>k)&1)<<j for k,j in enumerate(sidx)) for t in range(1<<len(chords))]
    allmask=(1<<len(chords))-1;seen=set();edges=[]
    for t in range(1<<(len(chords)-1)):
        tc=t^allmask
        if S[t]==t:continue
        key=tuple(sorted((t,tc,S[t],S[tc])))
        if key in seen:continue
        seen.add(key)
        a,b=C[t],C[tc];u,v=min(a,op[a]),min(b,op[b])
        if u not in free or v not in free:continue
        sig=int((a==u)!=(b==v));edges.append((u,v,sig))
    need(len(free)==22 and len(edges)==157,'physical Q1 sizes')
    adj=defaultdict(list)
    for i,(u,v,sig) in enumerate(edges):adj[u].append((v,i,sig));adj[v].append((u,i,sig))
    gauge={free[0]:0};tree=set();queue=deque(gauge)
    while queue:
        u=queue.popleft()
        for v,e,sig in adj[u]:
            if v not in gauge:gauge[v]=gauge[u]^sig;tree.add(e);queue.append(v)
    need(len(gauge)==22,'physical Q1 connected')
    cotree=[i for i in range(len(edges)) if i not in tree]
    norm=[s^gauge[u]^gauge[v] for u,v,s in edges]
    need(len(cotree)==136,'cycle rank')
    need(all(norm[i]==0 for i in tree),'tree voltage zero')
    negative_loops=sum(u==v and s for u,v,s in edges);need(negative_loops==9,'negative physical loops')
    pivot=next(i for i in cotree if norm[i]);others=[i for i in cotree if i!=pivot][:3]
    voltage=[0]*len(edges)
    for i in cotree:voltage[i]=norm[i]
    for j,i in enumerate(others):voltage[i]^=1<<(j+1)
    liftadj=defaultdict(list);projection_edge_counts=Counter()
    for i,(u,v,sig) in enumerate(edges):
        for x in range(16):
            y=x^voltage[i]
            liftadj[u,x].append((v,y));liftadj[v,y].append((u,x))
            # bit zero is corrected sign q. Undo vertex switching to physical sheet.
            pu=(x&1)^gauge[u];pv=(y&1)^gauge[v]
            need(pu^pv==sig,'physical edge projection')
            projection_edge_counts[i,pu]+=1
    seen={(free[0],0)};queue=deque(seen)
    while queue:
        u=queue.popleft()
        for v in liftadj[u]:
            if v not in seen:seen.add(v);queue.append(v)
    need(len(seen)==352,'extended cover connected')
    need(len(projection_edge_counts)==314 and set(projection_edge_counts.values())=={8},'eight lifts per physical edge')
    need(all(len(liftadj[v,x])==len(adj[v]) for v in free for x in range(16)),'local cover degrees with loop flags')
    return {'Q1_vertices':22,'Q1_edges':157,'cycle_rank':136,'negative_loops':9,'regular_extension_vertices':352,'regular_extension_edges':2512,'physical_lift_vertices':44,'physical_lift_edges':314,'extension_degree_over_physical':8,'voltage_pivot_edge':pivot,'extra_basis_edges':others}


def seed_profiles():
    n=6;pairs=[(i,j) for i in range(n) for j in range(i+1,n)];full=(1<<n)-1
    reports=[];masks=(153,185)
    def edge(mask,i,j):
        a,b=sorted((i,j));k=pairs.index((a,b));return bool((mask>>k)&1)==(i<j)
    for mask in masks:
        paths=[0]*(1<<n)
        for S in range(1,1<<n):
            @lru_cache(None)
            def ending(T,j):
                if T==(1<<j):return 1
                return sum(ending(T^(1<<j),i) for i in range(n) if i!=j and T>>i&1 and edge(mask,i,j))
            paths[S]=sum(ending(S,j) for j in range(n) if S>>j&1)
        @lru_cache(None)
        def covers(S):
            if not S:return (1,)+(0,)*n
            out=[0]*(n+1);bit=S&-S;T=S
            while T:
                if T&bit:
                    prev=covers(S^T)
                    for k in range(1,n+1):out[k]+=paths[T]*prev[k-1]
                T=(T-1)&S
            return tuple(out)
        profile=covers(full)
        need(profile==(0,29,131,162,77,15,1),'exact seed profile')
        for r in range(7):
            falling=1;value=0
            for j in range(1,7):
                falling*=r-j+1;value+=profile[j]*falling
            need(value==r**6+12*r**4+16*r**2,'seed slot polynomial')
        reports.append({'mask':mask,'profile':profile})
    need((153^185).bit_count()==1,'literal single arc')
    iso=False
    for p in permutations(range(n)):
        if all(edge(153,i,j)==edge(185,p[i],p[j]) for i,j in pairs):iso=True;break
    need(not iso,'nonisomorphic seeds')
    for mask in masks:
        for src in range(n):
            seen={src}
            while True:
                new=seen|{j for i in seen for j in range(n) if i!=j and edge(mask,i,j)}
                if new==seen:break
                seen=new
            need(len(seen)==n,'strong seed')
    return reports


def circuit_control():
    data,verts,colors=cycle_family(4);N=len(verts);bits=(N-1).bit_length()+3;size=1<<bits
    def unpack(z):return (z>>2)&1,(z>>1)&1,z&1
    def pack(z):s,y,t=z;return (s<<2)|(y<<1)|t
    F=[]
    for state in range(size):
        i,z=divmod(state,8)
        F.append(8*((i+1)%N)+pack(reg_step(unpack(z),colors[i],i%2)) if i<N else state)
    need(len(set(F))==size,'clock map permutation')
    cycles=[];seen=set()
    for start in range(size):
        if start in seen:continue
        cyc=[];x=start
        while x not in seen:seen.add(x);cyc.append(x);x=F[x]
        need(x==start,'cycle closes');cycles.append(cyc)
    counts=Counter(map(len,cycles));need(counts=={1:80,88:2},'complete clock cycle type')
    gates=[]
    def endpoint_swap(a,b):
        path=[a];cur=a
        for bit in range(bits):
            if (a^b)>>bit&1:cur^=1<<bit;path.append(cur)
        path_edges=list(zip(path,path[1:]));gates.extend(path_edges+path_edges[-2::-1])
    for cyc in cycles:
        # With rightmost-first functional notation, chronological endpoint swaps
        # (c0,c1),(c0,c2),... implement c0 -> c1 -> ... -> c0.
        for x in cyc[1:]:endpoint_swap(cyc[0],x)
    result=list(range(size))
    for a,b in gates:
        need((a^b).bit_count()==1,'single-bit controlled gate')
        result=[b if x==a else a if x==b else x for x in result]
    need(result==F,'full reversible schedule equals clock permutation')
    artifact={'bits':bits,'macro_permutation':F,'adjacent_transpositions':gates,'arithmetic_vertices':verts,'colors':''.join(colors)}
    Path(__file__).with_name('circuit_m4.json').write_text(json.dumps(artifact,indent=2)+'\n')
    return {'m':4,'bits':bits,'tournament_vertices':6*bits,'active_tournaments':N*8,'cycle_counts':dict(sorted(counts.items())),'scheduled_controlled_arc_gates':len(gates),'Hamilton_count':29**bits,'arithmetic_vertices':verts,'color_word':''.join(colors)}


def runner_controls():
    out=[]
    for n in range(3,13):
        v=list(range(1,n+1));ell=lcm(*[2//gcd(2,x) for x in v[:-1]],4//gcd(4,v[-1]))
        need(ell==(4 if n%2 else 2),'runner primitive period')
        t=Fraction(1,n+1);dist=[]
        for x in v:
            q=(x*t)%1;dist.append(min(q,1-q))
        need(min(dist)==t,'consecutive speed witness')
        out.append({'n':n,'attained_loneliness':str(t),'lift_return_time':ell})
    laps=[position_register(i,3*i) for i in range(5)]
    need([register_label(z) for z in laps]==['+e_0','+e_5','-e_0','-e_5','+e_0'],'runner sheet sequence')
    return out


def main():
    report={'algebra':algebra_controls(),'honeycomb':honey_controls(),'converse_extension':physical_q1(),'seed_profiles':seed_profiles(),'clock_circuit':circuit_control(),'runner_examples':runner_controls()}
    data=json.dumps(report,indent=2,sort_keys=True)+'\n'
    Path(__file__).with_name('results.json').write_text(data)
    print(data,end='')


if __name__=='__main__':main()
