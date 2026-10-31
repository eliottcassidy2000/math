#!/usr/bin/env python3
"""Inspect inert integer arrays only; never compile, import, or execute external C."""
import hashlib,json,re,sys
from pathlib import Path
from math import isqrt
checks=0

def check(ok,msg):
    global checks
    checks+=1
    if not ok:raise RuntimeError(msg)

def arrays(path):
    raw=Path(path).read_bytes()
    check(hashlib.sha256(raw).hexdigest()=='602c1919ee51101519a2f19f2c1817407006824f6f09c85a3b772a818a0b736e','pinned source digest')
    text=raw.decode('ascii')
    out={}
    for name in ('smallseq','data','S'):
        m=re.search(r'int\s+'+name+r'\[[^;=]+?=\s*(\{.*?\});',text,re.S)
        check(m is not None,'missing array '+name)
        body=m.group(1)
        check(re.fullmatch(r'[\s\d{},\-]+',body) is not None,'non-integer data syntax')
        out[name]=json.loads(body.replace('{','[').replace('}',']'))
    return out,hashlib.sha256(raw).hexdigest()

def pathcheck(a,n,cycle):
    check(len(a)==n and sorted(a)==list(range(1,n+1)),'coverage')
    for u,v in zip(a,a[1:]+(a[:1] if cycle else [])):
        check(isqrt(u+v)**2==u+v,'square edge')

def main():
    data,sha=arrays(sys.argv[1]);small=data['smallseq'];S=data['S'];words=data['data']
    check(len(small)==41,'small table length')
    small_ns={15,16,17,23}|set(range(25,41))
    for n,a in enumerate(small):
        if n in small_ns:pathcheck(a,n,n>=32)
        else:check(not a,'unexpected small entry')
    check(len(S)==3984,'base count')
    for n in range(41,2033):
        a,b=S[2*(n-41):2*(n-41)+2]
        for j,w in enumerate((a,b)):
            pathcheck(w,n+j,True)
            check(w[0]==1 and w[-1]==(3 if (n+j)%2 else 8),'base endpoint')
        pa={v:i%2 for i,v in enumerate(a,1)}
        pb={v:i%2 for i,v in enumerate(b,1)}
        check(all(pa[v]==pb[v] for v in a),'base pair parity')
    check(len(words)==196,'word count')
    profiles={};seen=set()
    for row_index,row in enumerate(words):
        q,d,p,ext=row[0];state=(d,p,ext)
        check(row_index==4*(d-24)+2*p+ext,'physical state addressing')
        check(q==49 and 24<=d<=72 and p in (0,1) and ext in (0,1),'state')
        check(state not in seen,'duplicate state');seen.add(state)
        check(len(row)==74,'atom count')
        offsets=set();singletons={};shifts={};end=None;first=None
        # Keep accumulated position as alpha*n+beta, evaluated mod2 only.
        posparity=1;blocks=0;extra=0
        for atom in row[1:]:
            check(len(atom)==3,'atom schema');c,orientation,e=atom
            if orientation==0:
                check(e==0 and 1<=c<=24 and c not in singletons,'singleton')
                singletons[c]=posparity;a=b=c;lengthparity=1
            else:
                check(orientation in (-1,1) and e in (0,1),'block schema')
                check(-24<=c<=24 and c not in offsets,'offset coverage')
                offsets.add(c);blocks+=1;extra+=e
                check(e==int(49+(-1)**p*c<=d+ext),'forced block choice')
                lengthparity=(p+e)%2
                a=49+c;b=49*(3 if lengthparity else 8)+(-1)**(lengthparity+1)*c
                shifts[c]=(posparity-1)%2 if orientation==1 else (posparity+lengthparity)%2
                if orientation==-1:a,b=b,a
            if first is None:first=a
            if end is not None:check(isqrt(end+a)**2==end+a,'connector')
            end=b;posparity=(posparity+lengthparity)%2
        check(offsets==set(range(-24,25)) and set(singletons)==set(range(1,25)),'atom coverage')
        check(blocks==49 and extra+24==d+ext,'symbolic exact length')
        check(first==1 and end==(3 if (p+d+ext)%2 else 8),'output endpoints')
        check(isqrt(first+end)**2==first+end,'closing edge')
        profiles[state]=(singletons,shifts)
    for d in range(24,73):
        for p in range(2):check(profiles[d,p,0]==profiles[d,p,1],'universal output pair parity')
    print(json.dumps({'status':'PASS','source_sha256':sha,'base_pairs':1992,'base_range':[41,2032],'small_cycle_range':[32,40],'transition_words':196,'paired_transitions':98,'induction_start':2033,'checks':checks,'external_code_executed':False,'scope':'Inert integer table verification plus stated affine induction proves cycles for every n>=32; small paths15,16,17,23,25..31'},sort_keys=True,indent=2))
if __name__=='__main__':main()
