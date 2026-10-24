#!/usr/bin/env python3
"""Independent domain, full-group, and hostile checks of the frozen decoder."""
import argparse, importlib.util, itertools as it, json, random, sys
sys.dont_write_bytecode = True
from pathlib import Path
from collections import Counter
from reference import need,cyc,relabel,stable,strong,autos,error,graph
checks=0
def require(x,m):
    global checks
    checks+=1;need(x,m)
def load(path):
    spec=importlib.util.spec_from_file_location('candidate',path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def construct(sizes,crossings):
    cells=[];start=0
    for p in sizes:cells.append(tuple(range(start,start+p)));start+=p
    rows=[0]*start
    for p,C in zip(sizes,cells):
        for i,r in enumerate(cyc(p)):rows[C[i]]=sum(1<<C[j] for j in range(p) if r>>j&1)
    for i,j in it.combinations(range(len(cells)),2):
        B=crossings[i,j]
        for x,a in enumerate(cells[i]):
            for y,b in enumerate(cells[j]):
                u,v=(a,b) if B[x][y] else (b,a);rows[u]|=1<<v
    return tuple(rows)
def mat(p,s=1,c=0):return tuple(tuple(int(y==(s*x+c)%p) for y in range(p)) for x in range(p))
def permat(p,perm):return tuple(tuple(int(y==perm[x]) for y in range(p)) for x in range(p))
def uni(a,b,v):return tuple((v,)*b for _ in range(a))
def closure(n,gens):
    identity=tuple(range(n));seen={identity};todo=[identity]
    while todo:
        g=todo.pop()
        for h in gens:
            v=tuple(h[g[i]] for i in range(n))
            if v not in seen:seen.add(v);todo.append(v)
    return seen
def check(D,rows):
    out=D.decode(rows);group=autos(rows);e,orbits=error(rows,group)
    require({frozenset(c) for c in out['cells']}=={frozenset(c) for c in stable(rows)},'actual stable partition')
    require(out['group_order']==len(group),'full automorphism order');require(closure(len(rows),out['generators'])==set(group),'actual generated permutation group')
    size={v:len(o) for o in orbits for v in o};require(out['orbit_sizes']==[size[v] for v in range(len(rows))],'actual vertex orbits');require(out['error']==e,'literal refinement error');return out
def regulars(n):
    target=(n-1)//2;rows=[0]*n;deg=[0]*n
    def rec(i):
        if i==n:
            if all(d==target for d in deg):yield tuple(rows)
            return
        k=target-deg[i];tail=tuple(range(i+1,n))
        if not 0<=k<=len(tail):return
        for out in it.combinations(tail,k):
            outgoing=set(out);old=list(deg);ri=rows[i]
            for j in tail:
                if j in outgoing:rows[i]|=1<<j;deg[i]+=1
                else:rows[j]|=1<<i;deg[j]+=1
            if all(deg[j]<=target and deg[j]+n-i-2>=target for j in range(i+1,n)):yield from rec(i+1)
            rows[i]=ri
            for j in tail:
                if j not in outgoing:rows[j]^=1<<i
            deg[:]=old
    yield from rec(0)
def hostile_and_multicell(D):
    result={};rng=random.Random(202651)
    for p in (3,5,7,11,13):
        result[f'internal_{p}']=check(D,cyc(p))['group_order'];pair=construct([p,p],{(0,1):mat(p,2,1)})
        require(check(D,pair)['group_order']==p,'nonunit gain pair')
        good={(0,1):mat(p,2,1),(0,2):mat(p,4%p,2),(1,2):mat(p,2,2)}
        balanced=construct([p]*3,good);require(check(D,balanced)['group_order']==p,'balanced gain cycle')
        bad=dict(good);bad[0,2]=mat(p,1 if 4%p!=1 else 2,2);frustrated=construct([p]*3,bad);require(check(D,frustrated)['group_order']==1,'unbalanced gain cycle')
        for rows in (pair,balanced,frustrated):
            pmap=list(range(len(rows)));rng.shuffle(pmap);out=check(D,relabel(rows,pmap));require(out['error']==D.decode(rows)['error'],'label equivariance')
            mask=(1<<len(rows))-1;op=tuple(mask^row^(1<<i) for i,row in enumerate(rows));require(check(D,op)['error']==out['error'],'converse compatibility')
        if p>=5:
            perm=list(range(p));perm[-1],perm[-2]=perm[-2],perm[-1];pin={(0,1):permat(p,perm),(0,2):mat(p,2),(1,2):mat(p,3)}
            rows=construct([p]*3,pin);require(check(D,rows)['group_order']==1,'rigid pin propagation')
            two=construct([p,p],{(0,1):permat(p,perm)});require(check(D,two)['error']==2*p*(p-1),'both rigid endpoints pinned')
    mixed=construct([3,5,7],{(0,1):uni(3,5,1),(0,2):uni(3,7,1),(1,2):uni(5,7,1)});require(check(D,mixed)['group_order']==105,'mixed prime product')
    B=mat(3);comp=tuple(tuple(1-v for v in row) for row in B);rows=construct([3,3,3],{(0,1):B,(0,2):comp,(1,2):B});require(list(map(len,stable(rows)))==[9],'hostile equitable blocks collapse')
    try:D.decode(rows)
    except ValueError:require(True,'coarsened composite cell rejected')
    else:raise RuntimeError('accepted arbitrary equitable blocks instead of stable cells')
    paley=tuple(sum(1<<((i+d)%7) for d in (1,2,4)) for i in range(7))
    for name,rows in [('paley7',paley),('bad7',(74,84,97,22,37,11,56))]:
        try:D.decode(rows)
        except ValueError:require(True,'non-Rp prime cell rejected');result[name]='rejected'
        else:raise RuntimeError('accepted non-Rp prime regular cell '+name)
    for name,rows in [('loop',(1,)),('missing_arc',(0,0)),('double_arc',(2,1)),('out_of_range',(2,)),('negative',(-1,))]:
        try:D.decode(rows)
        except (ValueError,TypeError):require(True,'malformed tournament rejected')
        else:raise RuntimeError('accepted malformed tournament '+name)
    return result
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--decoder',default='source/decoder.py');ap.add_argument('--quick',action='store_true');ap.add_argument('--output',default='receipt.json');a=ap.parse_args();D=load(a.decoder);result={'hostile_and_multicell':hostile_and_multicell(D)}
    if not a.quick:
        universe=Counter()
        for n in range(7):
            for bits in range(1<<(n*(n-1)//2)):check(D,graph(n,bits));universe[n]+=1
        result['all_labelled_through_six']={str(k):v for k,v in sorted(universe.items())};reg=Counter()
        for rows in regulars(7):
            group=autos(rows)
            try:out=D.decode(rows)
            except ValueError:require(len(group) in (3,21),'regular-seven rejection class');reg['rejected_'+str(len(group))]+=1
            else:require(len(group)==7,'accepted non-cyclic regular-seven class');check(D,rows);reg['accepted_7']+=1
        require(dict(reg)=={'accepted_7':720,'rejected_3':1680,'rejected_21':240},'regular-seven full recognition census');result['all_regular_seven']=dict(reg)
    result['checks']=checks;result['status']='PASS';Path(a.output).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result,sort_keys=True))
if __name__=='__main__':main()
