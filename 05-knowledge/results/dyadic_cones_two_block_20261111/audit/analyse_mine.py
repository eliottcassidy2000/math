# Auditor: union of own census (orders 2..10 roots) and comparison with the packet's claims and files.
import sys
def load(fn):
    d={}; W={}
    for l in open(fn):
        if l.startswith('#'): continue
        if l.startswith('W'):
            p=l.split(); c,r,j=int(p[1]),int(p[2]),int(p[3]); W[(c,r)]=(j,list(map(int,p[4:]))); continue
        c,r,j=map(int,l.split()); d[(c,r)]=j
    return d,W
U={}; per={}; Wall={}
for o in range(2,11):
    d,W=load(f'ord{o}.txt'); per[o]=d
    for k,j in d.items():
        if k not in U or j<U[k]: U[k]=j
    for k,w in W.items():
        Wall.setdefault(k,{})[o]=w
win={k:j for k,j in U.items() if abs(k[1])<=2045}
print('types c<=255 |r|<=2045 (orders<=10):',len(win))
# packet's files: union of order9 and order10 blocks
th={}
for fn in ['../out/blocks_roots_order9.txt','../out/blocks_roots_order10.txt']:
    for l in open(fn):
        if l.startswith('#'): continue
        c,r,j=map(int,l.split()); k=(c,r)
        if k not in th or j<th[k]: th[k]=j
print('packet union order9+10:',len(th))
mine910={}
for o in (9,10):
    for k,j in per[o].items():
        if abs(k[1])<=2045 and (k not in mine910 or j<mine910[k]): mine910[k]=j
print('mine order 9+10 == packet order 9+10:', mine910==th)
d10={k:j for k,j in per[10].items() if abs(k[1])<=2045}
th10={}
for l in open('../out/blocks_roots_order10.txt'):
    if l.startswith('#'): continue
    c,r,j=map(int,l.split()); th10[(c,r)]=j
print('order-10 table identical:', d10==th10, len(d10), len(th10))
# 1) cones
cones={r:j for (c,r),j in U.items() if c==1}
print('cone gaps |r|<1000:',[r for r in range(-999,1000,2) if r not in cones])
print('max j cones in [-791,647]:',max(cones[r] for r in range(-791,648,2)))
# boxes
def box(cmax,rmax):
    miss=[(c,r) for c in range(1,cmax+1,2) for r in range(-rmax,rmax+1,2) if (c,r) not in U]
    mj=max([U[(c,r)] for c in range(1,cmax+1,2) for r in range(-rmax,rmax+1,2) if (c,r) in U])
    return miss[:10],len(miss),mj
print('box c<=31 |r|<=275:',box(31,275))
print('box c<=31 |r|<=277:',box(31,277))
print('box c<=15 |r|<=295:',box(15,295))
print('box c<=15 |r|<=297:',box(15,297))
print('box c<=63 |r|<=9:',box(63,9), 'box c<=63 |r|<=11', box(63,11))
print('(245,1) present:',(245,1) in U)
# orders needed: types in box c<=31 |r|<=275 whose ONLY witnesses (any j) are order 10, and with c>7
only10=[(c,r) for c in range(1,32,2) for r in range(-275,276,2) if (c,r) in per[10] and all((c,r) not in per[o] for o in range(2,10))]
print('box types present only at order 10:',len(only10),' with c>7:',len([k for k in only10 if k[0]>7]))
# types whose LEAST exponent is attained only at order 10 (min j over order<=9 larger)
# witness-file coverage: c=3..31 order 9, c=3,5,7 order 10, c=1 cones
import re
cov=set()
for fn in ['../out/types_c3_to_c31_roots_order9.txt','../out/types_c3_c5_c7_roots_order10.txt']:
    for l in open(fn):
        m=re.match(r'r=(-?\d+) j=(\d+) m=(\d+) O=(\d+) c=(\d+)',l)
        if m: cov.add((int(m.group(5)),int(m.group(1))))
for o in (8,9,10):
    for l in open(f'../out/cones_roots_order{o}.txt'):
        m=re.match(r'r=(-?\d+) j=',l)
        if m: cov.add((1,int(m.group(1))))
nowit=[(c,r) for c in range(1,32,2) for r in range(-295 if c<=15 else -275,(295 if c<=15 else 275)+1,2) if (c,r) not in cov]
print('types in the two boxes without a witness line on file:',len(nowit), nowit[:20])
import json; json.dump(nowit,open('nowitness_types.json','w'))
print('--- breakdown of types without witness line')
at9=[k for k in nowit if any(k in per[o] for o in range(2,10))]
print('of these present at order<=9 in own census:',len(at9),at9[:10])
for k in at9[:5]:
    print(k,{o:per[o].get(k) for o in range(2,11) if k in per[o]})
from collections import Counter
print('by c:',sorted(Counter(c for c,r in nowit).items()))
print('--- Theorem 3 offsets from own census (least j over orders <= 10)')
print('1+2^a:',[(a,U.get((1,1+2**a))) for a in range(1,12)])
print('2^a-1:',[(a,U.get((1,2**a-1))) for a in range(1,12)])
