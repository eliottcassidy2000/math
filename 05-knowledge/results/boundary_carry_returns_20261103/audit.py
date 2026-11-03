#!/usr/bin/env python3
"""Independent skeleton/word audit. No imports from the producer modules."""
from itertools import permutations,combinations
from math import comb
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent


def need(ok,message):
    if not ok:raise RuntimeError(message)


def matrix(rows):
    n=len(rows);t=[[bool(rows[i]>>j&1)for j in range(n)]for i in range(n)]
    need(all(not t[i][i]for i in range(n)),'no loops')
    need(all(t[i][j]!=t[j][i]for i,j in combinations(range(n),2)),'tournament pairs')
    return t


def skeleton(rows,v):
    t=matrix(rows);n=len(t);c=[0]*(n+1)
    for p in permutations(i for i in range(n)if i!=v):
        required=0;optional=int(t[v][p[0]])+int(t[p[-1]][v]);ok=True
        for a,b in zip(p,p[1:]):
            no=t[a][b];yes=t[a][v]and t[v][b]
            if not(no or yes):ok=False;break
            if not no:required+=1
            elif yes:optional+=1
        if ok:
            for j in range(optional+1):c[required+j]+=comb(optional,j)
    return c


def matching(rows,v):
    t=matrix(rows);out=[i for i in range(len(t))if t[v][i]]
    inside=[i for i in range(len(t))if t[i][v]]
    def go(k,used):
        if k==len(out):return 0
        return max([go(k+1,used)]+[1+go(k+1,used|{i})for i in inside if i not in used and t[out[k]][i]])
    return go(0,set())


def polynomial(c,m):
    return sum(c[r]*sum((-1)**(r-j)*comb(r,j)*j**m for j in range(r+1))for r in range(1,len(c)))


def perform(rows,path):
    t=matrix(rows);n=len(t)
    need(sorted(path)==list(range(n)),'path uses exact source names')
    need(all(t[a][b]for a,b in zip(path,path[1:])),'directed path guard')
    spine={frozenset((a,b))for a,b in zip(path,path[1:])}
    return [sum(1<<j for j in range(n)if i!=j and (t[i][j]if frozenset((i,j))in spine else not t[i][j]))for i in range(n)]


def all_prefix_defects(rows):
    n=len(rows);full=(1<<n)-1;ends=[[0]*n for _ in range(full+1)]
    def walk(last,mask):
        ends[mask][last]+=1
        for j in range(n):
            if not(mask>>j&1)and rows[last]>>j&1:walk(j,mask|(1<<j))
    for i in range(n):walk(i,1<<i)
    h=[sum(e)for e in ends];h[0]=1
    return [sum(ends[s][v]*h[full^s]for s in range(1,full+1))-h[full]for v in range(n)]


def ordered_covers_literal(rows,k):
    n=len(rows);answer=0
    for p in permutations(range(n)):
        bad=sum(not(rows[a]>>b&1)for a,b in zip(p,p[1:]))
        if bad<=k-1<=n-1:answer+=comb(n-1-bad,k-1-bad)
    return answer


def add(rows,mask):
    n=len(rows);return [r+((1<<n)if mask>>i&1 else 0)for i,r in enumerate(rows)]+[((1<<n)-1)^mask]


def expand(rows,v,m):
    addresses=[(i,j)for i in range(len(rows))for j in range(m if i==v else 1)]
    rows2=[sum(1<<j for j,(b,y)in enumerate(addresses)if (a==b and x<y)or(a!=b and rows[a]>>b&1))for a,x in addresses]
    def path(p):return [addresses.index((i,j))for i in p for j in range(m if i==v else 1)]
    return rows2,path


def uniform_program(r):
    src=r['source'];v=r['clone'];target=r['target'];c=r['offset']
    cs=skeleton(src,v);ct=skeleton(target,v)
    for m in range(1,max(len(cs),len(ct))+2):need(polynomial(ct,m)==3*polynomial(cs,m)+c,'coefficient identity')
    for m in (1,2,4,7):
        t,p=expand(src,v,m)
        got=perform(perform(t,p(r['first_path'])),p(r['second_path']))
        mask=sum(1<<j for j,(i,k)in enumerate([(i,k)for i in range(len(src))for k in range(m if i==v else 1)])if r['mask']>>i&1)
        got=add(got,mask)
        need(got==expand(target,v,m)[0],'literal common-frame lift')


if __name__=='__main__':
    banks=json.loads((HERE/'small_banks.json').read_text());rooted=0
    for key,bank in banks.items():
        n=int(key);seen=set();pairs=list(combinations(range(n),2))
        for t in bank:
            orbit={sum(1<<k for k,(i,j)in enumerate(pairs)if t[p[i]]>>p[j]&1)for p in permutations(range(n))}
            need(not(seen&orbit),'small-bank orbits are disjoint');seen.update(orbit)
        need(len(seen)==1<<len(pairs),'independent small-bank orbit coverage')
        if int(key)<2:continue
        for t in bank:
            for v in range(len(t)):
                c=skeleton(t,v);r=max(i for i,x in enumerate(c)if i and x)
                need(r==matching(t,v)+1,'matching equals highest run count');rooted+=1
                need(not(c[1]in(287,289,293,295)and c[2:4]==[270,72]and not any(c[4:])),
                     'no grandchild family has base order at most six')
    base=[10,36,1,22,7,29];mate=[26,36,1,22,6,29]
    need(skeleton(base,3)[1:4]==[31,30,8],'31 source coefficients')
    need(skeleton(mate,3)[1:4]==[33,26,6],'33 mate coefficients')
    atlas_checked=0
    for file in sorted(HERE.glob('catalog*.jsonl')):
        for line in file.read_text().splitlines():
            row=json.loads(line);t=row['rows'];c=skeleton(t,0);scale=3**row['level']
            need(c[1:4]==[row['H'],30*scale,8*scale]and not any(c[4:]),'independent atlas polynomial')
            for start in range(len(t)):
                seen={start}
                while True:
                    nxt=seen|{j for i in seen for j in range(len(t))if t[i]>>j&1}
                    if nxt==seen:break
                    seen=nxt
                need(len(seen)==len(t),'independent atlas strongness')
            atlas_checked+=1
    progs=json.loads((HERE/'scout31.json').read_text())['programmes']
    for r in progs:uniform_program(r)
    # Audit every saved next-level positive, not only the chosen 97->293 route.
    next_count=0
    for path in sorted(HERE.glob('next_*.json')):
        data=json.loads(path.read_text())
        if 'programmes'not in data:continue
        for r in data['programmes']:uniform_program(r);next_count+=1
    results=json.loads((HERE/'returns.json').read_text())
    hypothesis_controls=[progs[0],progs[1],json.loads((HERE/'next_0.json').read_text())['programmes'][1],*results]
    for r in hypothesis_controls:
        old=len(r['source']);v=r['clone'];source=r['source'];target=r['target']
        need(len(target)-old<=3,'hypothesis permanent-birth control')
        need(sum(bool((source[a]^target[a])>>b&1)for a,b in combinations(range(old),2))<=4,
             'hypothesis exterior-edit control')
        need(all(bool(source[v]>>j&1)==bool(target[v]>>j&1)for j in range(old)),'hypothesis root-star control')
    for r in results:
        t=r['source'];v=r['clone'];cs=skeleton(t,v)
        for mask in r['permanent_birth_masks']:t=add(t,mask)
        for recipe in r['exterior_flips']:
            ext=t
            for mask in recipe['masks']:ext=add(ext,mask)
            flipped=list(t);a,b=recipe['a'],recipe['b'];flipped[a]^=1<<b;flipped[b]^=1<<a
            for m in (1,2,4,7):
                expanded,p=expand(ext,v,m)
                got=perform(perform(expanded,p(recipe['first_path'])),p(recipe['second_path']))
                old=len(t)+m-1;mask=(1<<old)-1
                need([row&mask for row in got[:old]]==expand(flipped,v,m)[0],'independent exterior compiler')
            t=flipped
        need(t==r['target'],'all physical labels reach target')
        ct=skeleton(t,v)
        for m in range(1,12):need(polynomial(ct,m)==3*polynomial(cs,m)+r['offset'],'repaired uniform return coefficients')
        need(ct[1:4]==[r['values'][0][2],270,72]and not any(ct[4:]),'grandchild coefficients')
    # Exact marked isomorphisms justify only census deduplication, never a free move.
    isos=[]
    for left,right in ((7,0),(11,10)):
        a,b=progs[left]['target'],progs[right]['target'];fixed={1:1,5:5,3:3};other=[0,2,4,6]
        found=None
        for p in permutations(other):
            f=dict(fixed);f.update(zip(other,p))
            if all(bool(a[i]>>j&1)==bool(b[f[i]]>>f[j]&1)for i in range(7)for j in range(7)):
                found=[f[i]for i in range(7)];break
        need(found is not None,'marked census identification');isos.append({'from':left,'to':right,'map':found})
    for r in json.loads((HERE/'boundary.json').read_text()):
        for m in range(1,5):
            d=all_prefix_defects(expand(r['source'],r['clone'],m)[0])
            for q,value in r['values'][m-1].items():
                q=int(q);index=q+(m-1 if q>r['clone']else 0)
                need(d[index]==value,'independent all-prefix boundary count')
    moments=json.loads((HERE/'moment_budget.json').read_text())
    for bank in moments:
        for r in bank['witnesses']:
            need(matching(r['rows'],0)==r['matching'],'actual matching witness')
            for m in range(1,4):
                need(ordered_covers_literal(expand(r['rows'],0,m)[0],r['components'])==r['values'][m-1],
                     'independent permutation/cut moment control')
    report={'rooted_small_bank_checks':rooted,'first_programmes':len(progs),'second_programmes':next_count,
            'repaired_returns':len(results),'marked_census_identifications':isos,
            'boundary_samples_independently_checked':16,
            'atlas_positive_presentations_checked':atlas_checked,
            'moment_budget_dimensions':[r['dimension']for r in moments],
            'method':'same agent, independently implemented skeleton and boolean-word algorithms; not an independent agent audit'}
    (HERE/'audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
