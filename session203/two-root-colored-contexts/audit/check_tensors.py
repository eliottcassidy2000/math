#!/usr/bin/env python3
"""Independent ordered-gap choice enumeration, then commutative aggregation."""
from itertools import permutations, combinations
from collections import Counter
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
ZERO=(0,0,0,0)
CASES=[
 ('A',[2,24,3,37,13,23],[15,63],[2,8,67,101,15,87,19],[79,127],(1,4)),
 ('B',[32,5,1,39,47,70,31],[0,5],[0,5,129,39,47,199,159,27],[0,5],(0,5)),
]

def need(b,s):
    if not b:raise RuntimeError(s)

def valid(a):
    n=len(a)
    for i,row in enumerate(a):
        need(0<=row<1<<n and not row&(1<<i),'row')
        for j in range(i):need(((row>>j)&1)+((a[j]>>i)&1)==1,'pair')

def plus(*polys):
    ans=Counter()
    for poly in polys:
        for key,value in poly.items():ans[key]+=value
    return {key:value for key,value in ans.items() if value}

def scale(poly,c):return {key:c*value for key,value in poly.items() if c*value}

def times(a,b):
    result=Counter()
    for x,c in a.items():
        for y,d in b.items():result[tuple(x[i]+y[i] for i in range(4))]+=c*d
    return {key:value for key,value in result.items() if value}

ONE={ZERO:1}
VARIABLES=[{tuple(int(i==j) for i in range(4)):1} for j in range(4)]
x,y,z,w=VARIABLES
H=plus(ONE,x,y,z,w)
W=plus(y,scale(z,-1),times(x,y),times(x,w),times(y,y),times(y,w))
F=plus(y,times(x,y),times(y,y))
J=plus(y,times(y,y),times(y,w))
FACTORS={
 'A':plus(times(plus(scale(ONE,4),x),W),times(plus(x,z),F),times(y,plus(scale(ONE,4),scale(x,3),scale(y,4),scale(z,3),scale(w,4)))),
 'B':plus(times(plus(scale(ONE,21),scale(y,42),scale(w,20)),W),scale(times(plus(z,w),J),20),scale(times(y,H),22)),
}

def tensor(core,masks):
    answer=Counter();order_count=0;legal_orders=0
    n=len(core)
    for order in permutations(range(n)):
        order_count+=1;choices=[]
        for gap in range(n+1):
            left=order[gap-1] if gap else None
            right=order[gap] if gap<n else None
            options=[]
            if left is None or right is None or core[left]&(1<<right):options.append(None)
            for kind in range(4):
                first,last=divmod(kind,2)
                if (left is None or not masks[first]&(1<<left)) and (right is None or masks[last]&(1<<right)):
                    options.append(kind)
            if not options:break
            choices.append(options)
        if len(choices)!=n+1:continue
        legal_orders+=1
        def choose(g,word):
            if g==len(choices):answer[word]+=1;return
            for item in choices[g]:choose(g+1,word if item is None else word+(item,))
        choose(0,())
    aggregated=Counter()
    for word,value in answer.items():aggregated[tuple(word.count(i) for i in range(4))]+=value
    return dict(answer),dict(aggregated),order_count,legal_orders

def encode(poly):return [{'powers':list(key),'coefficient':value} for key,value in sorted(poly.items())]
def decode(rows):return {tuple(row['powers']):row['coefficient'] for row in rows}

def main():
    published={item['family']:item for item in json.loads((HERE/'PRODUCER_TENSORS.json').read_text())['tensors']}
    cpp={};current=None
    for line in (HERE/'PRODUCER_ORDERED_TENSORS.txt').read_text().splitlines():
        parts=line.split()
        if not parts:continue
        if parts[0]=='FAMILY':current=parts[1][0];cpp[current]={}
        elif parts[0]=='SIGNS':continue
        else:cpp[current][tuple(int(t,2) for t in parts[1].split(','))]=int(parts[2])
    out=[];all_orders=0
    for name,p,pm,c,cm,flip in CASES:
        valid(p);valid(c);n=len(p)
        need(len(c)==n+1,'one birth')
        need([m&((1<<n)-1) for m in cm]==pm,'retained body incidences')
        changed=[(i,j) for i,j in combinations(range(n),2) if bool(p[i]&(1<<j))!=bool(c[i]&(1<<j))]
        need(changed==[flip] and p[flip[0]]&(1<<flip[1]),'old arc reversal')
        po,pa,pcount,plegal=tensor(p,pm);co,ca,ccount,clegal=tensor(c,cm);all_orders+=pcount+ccount
        excess=plus(ca,scale(pa,-2),scale(H,-1))
        need(excess==FACTORS[name],'exact factorization')
        need(pa==decode(published[name]['parent']) and ca==decode(published[name]['child']) and excess==decode(published[name]['excess']),'producer full tensor disagreement')
        ordered={}
        for key in set(po)|set(co)|{(),(0,),(1,),(2,),(3,)}:
            value=co.get(key,0)-2*po.get(key,0)-int(len(key)<=1)
            if value:ordered[key]=value
        need(ordered==cpp[name],'ordered C++ coefficient disagreement')
        out.append({'family':name,'parent':encode(pa),'child':encode(ca),'excess':encode(excess),
                    'ordered_excess':[{'types':list(key),'coefficient':value} for key,value in sorted(ordered.items())],
                    'parent_core_orders':pcount,'child_core_orders':ccount,
                    'parent_legal_orders':plegal,'child_legal_orders':clegal,
                    'parent_ordered_terms':len(po),'child_ordered_terms':len(co)})
    need(all_orders==51120,'core permutation universe')
    print(json.dumps({'status':'PASS','core_orders_examined':all_orders,'families':out,
                      'method':'Enumerate every permitted ordered choice of occupied core gaps before aggregating; no body-size interpolation.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
