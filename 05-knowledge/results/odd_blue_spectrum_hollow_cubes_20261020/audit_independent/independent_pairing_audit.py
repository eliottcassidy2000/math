#!/usr/bin/env python3
"""Independent marked-skeleton audit; imports no producer/census modules.

Universe: polarized anti-involutions on 2m+1 vertices. X={0,...,m-1},
Y={m,...,2m-1}, root=2m, r(i)=i+m mod 2m. Each skeleton has a
m-vertex tournament A and a symmetric off-diagonal 0/1 matrix B.
Matched-arc bits b complete it. Quotient only by pair-index permutations.
Full vertex-permutation HP/automorphism checks are independent at n<=7.
"""
import itertools as it
import collections as co
import math,json,time
from pathlib import Path
OUT=Path(__file__).resolve().parent


def make_e(m,e):
    A=[[0]*m for _ in range(m)];B=[[0]*m for _ in range(m)]
    pairs=list(it.combinations(range(m),2));k=len(pairs)
    for t,(i,j) in enumerate(pairs):
        A[i][j]=(e>>t)&1;A[j][i]=1-A[i][j]
        B[i][j]=B[j][i]=(e>>(t+k))&1
    return A,B


def ecode(A,B,p):
    m=len(A);pairs=list(it.combinations(range(m),2));k=len(pairs)
    return sum(A[p[i]][p[j]]<<t for t,(i,j) in enumerate(pairs))+sum(B[p[i]][p[j]]<<(t+k) for t,(i,j) in enumerate(pairs))


def perm_b(b,p):return sum(((b>>v)&1)<<i for i,v in enumerate(p))


def graph(A,B,b):
    m=len(A);n=2*m+1;f=2*m;G=[[0]*n for _ in range(n)]
    for i in range(m):
        G[i][f]=G[f][m+i]=1
        G[i][m+i]=(b>>i)&1;G[m+i][i]=1-G[i][m+i]
        for j in range(m):
            if i!=j:
                G[i][j]=A[i][j];G[m+i][m+j]=A[j][i]
                G[i][m+j]=B[i][j];G[m+j][i]=1-B[i][j]
    return G


def c3(G):
    n=len(G);return math.comb(n,3)-sum(math.comb(sum(r),2) for r in G)


def sym_paths(G):
    m=(len(G)-1)//2;f=2*m;r=lambda v:v+m if v<m else v-m
    answer=[]
    for p in it.permutations(range(m)):
        for s in range(1<<m):
            q=tuple(v+(((s>>v)&1)*m) for v in p)
            P=q+(f,)+tuple(r(v) for v in q[::-1])
            if all(G[u][v] for u,v in zip(P,P[1:])):answer.append(P)
    return set(answer)


def direct(G):
    """Full S_n path-orbit computation, not polarized quotient counting."""
    n=len(G);pairs=list(it.combinations(range(n),2));aut=[];anti=[];fibre=set();hp=0
    for p in it.permutations(range(n)):
        if all(G[p[i]][p[j]]==G[i][j] for i,j in pairs):aut.append(p)
        if all(G[p[i]][p[j]]!=G[i][j] for i,j in pairs):anti.append(p)
        if all(G[p[i]][p[i+1]] for i in range(n-1)):
            hp+=1;fibre.add(tuple(G[p[i]][p[j]] for i,j in pairs))
    blue=sum(all(w[t]==w[pairs.index((n-1-j,n-1-i))] for t,(i,j) in enumerate(pairs)) for w in fibre)
    invol=[p for p in anti if all(p[p[i]]==i for i in range(n))]
    r=tuple(list(range(n//2,n-1))+list(range(n//2))+[n-1])
    central=[p for p in aut if all(p[r[i]]==r[p[i]] for i in range(n))]
    # Every anti-involution conjugate to r by an actual automorphism.
    conj=set()
    for a in aut:
        ai=[a.index(i) for i in range(n)]
        conj.add(tuple(a[r[ai[i]]] for i in range(n)))
    assert set(invol)==conj
    assert len(fibre)*len(aut)==hp
    return {'H':hp,'aut':len(aut),'fibre':len(fibre),'blue':blue,'anti_involutions':len(invol),'centralizer':len(central)}


def run(m):
    start=time.time();n=2*m+1;perms=list(it.permutations(range(m)));reps={};pairing={}
    skeleton_groups=co.Counter();orbit_sizes=co.Counter();raw_checks=0
    for e in range(1<<(m*(m-1))):
        A,B=make_e(m,e);codes=[ecode(A,B,p) for p in perms]
        H=[p for p,c in zip(perms,codes) if c==e]
        skeleton_groups[len(H)]+=1;assert len(H)%2
        # Canonical rooted isomorphism code, root forced into slot 0.
        roots=[min(c for p,c in zip(perms,codes) if p[0]==i) for i in range(m)]
        O=[i for i in range(m) if roots[i]==min(roots)]
        assert set(O)=={h[O[0]] for h in H};assert len(O)%2
        flip=sum(1<<i for i in O);orbit_sizes[len(O)]+=1
        P=sym_paths(graph(A,B,0))
        for b in range(1<<m):
            d=b^flip;stabilizer=[p for p in H if perm_b(b,p)==b]
            assert stabilizer==[p for p in H if perm_b(d,p)==d]
            G=graph(A,B,b);Gd=graph(A,B,d)
            assert (c3(G)-c3(Gd))%2
            assert sym_paths(G)==P==sym_paths(Gd)
            # Equivariance including all markings and every pair relabeling.
            key=min((c,perm_b(b,p)) for p,c in zip(perms,codes))
            target=min((c,perm_b(d,p)) for p,c in zip(perms,codes))
            assert len(P)%len(stabilizer)==0
            data={'e':e,'b':b,'c3':c3(G),'blue':len(P)//len(stabilizer),'centralizer':len(stabilizer),'P_r':len(P),'orbit':O,'target':target}
            if key in reps:
                assert reps[key]['blue']==data['blue'] and reps[key]['c3']==data['c3']
                assert pairing[key]==target
            else:reps[key]=data;pairing[key]=target
            raw_checks+=1
    for key,target in pairing.items():
        assert pairing[target]==key and key!=target
        assert reps[key]['blue']==reps[target]['blue']
        assert (reps[key]['c3']-reps[target]['c3'])%2
    spectra=[dict(sorted(co.Counter(r['blue'] for r in reps.values() if r['c3']%2==parity).items())) for parity in range(2)]
    assert spectra[0]==spectra[1]
    direct_rows=[]
    if n<=7:
        for key,row in reps.items():
            A,B=make_e(m,row['e']);G=graph(A,B,row['b']);exact=direct(G)
            assert exact['blue']==row['blue']
            assert exact['centralizer']==row['centralizer']
            assert exact['blue']*exact['centralizer']==row['P_r']
            row.update(exact)
        for key,row in reps.items():
            other=reps[pairing[key]]
            if key<pairing[key]:direct_rows.append({'source':list(key),'target':list(pairing[key]),'blue':row['blue'],'source_c3':row['c3'],'target_c3':other['c3'],'source_H':row['H'],'target_H':other['H'],'source_aut':row['aut'],'target_aut':other['aut'],'source_fibre':row['fibre'],'target_fibre':other['fibre'],'centralizer':row['centralizer'],'P_r':row['P_r']})
    report={'n':n,'raw_marked_completions':raw_checks,'SC_classes':len(reps),'pairs':len(reps)//2,'skeleton_group_orders':dict(skeleton_groups),'selected_orbit_sizes':dict(orbit_sizes),'each_parity_spectrum':spectra[0],'full_Sn_checks':len(reps) if n<=7 else 0,'pair_examples':direct_rows,'seconds':round(time.time()-start,3)}
    print(json.dumps(report),flush=True)
    return report

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--max-m',type=int,default=4);a=p.parse_args()
    results=[run(m) for m in range(1,a.max_m+1)]
    (OUT/'independent_pairing_audit.json').write_text(json.dumps(results,indent=2)+'\n')
