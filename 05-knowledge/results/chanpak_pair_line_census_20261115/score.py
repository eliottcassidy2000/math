"""Compute finite line hits; no claimed tournament nonexistence outside the census."""
import argparse, hashlib, json
from collections import defaultdict

def mark_w(mask):
    """TT host: first-path end and second-path start both marked."""
    bits=[i for i in range(mask.bit_length()) if mask>>i&1]
    return sum((1 << (j-i-1))+(j-i==1) for j in bits for i in bits if i<j)

def band(A,lo,hi):
    a=sorted(A.intersection(range(lo,hi+1,2)))
    gaps=[i for i in range(lo,hi+1,2) if i not in A][:10]
    longest=(0,None,None)
    if a:
        current=a[0]
        for i,b in enumerate(a):
            if i==0 or b==a[i-1]+2:
                longest=max(longest,((b-current)//2+1,current,b))
            else:
                current=b
                longest=max(longest,(1,b,b))
    return {"hit":len(a),"possible":(hi-lo)//2+1,"first_misses":gaps,"longest_run":longest}

def read(path):
    inventories=defaultdict(list)
    for line in open(path):
        n,h,p,x,z=line.split()
        inventories[int(n)].append((int(h),int(p)+2*int(h),int(x),z))
    return inventories

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('profiles')
    ap.add_argument('--hi',type=int,default=59999)
    ap.add_argument('--marks',type=int,default=6)
    ap.add_argument('--lo',type=int,default=69)
    ap.add_argument('--orders',nargs='+',type=int,default=[8,9,10])
    args=ap.parse_args()
    result={}
    inventories=read(args.profiles)
    readers={'H':(1,0),'Q+H':(1,1),'2Q+H':(1,2),'4Q+H':(1,4),'8Q+H':(1,8),'3Q-H':(-1,3),'6Q-H':(-1,6),'1024Q-H':(-1,1024)}
    marks={m: (1-mark_w(m),m) for m in range(1,1<<args.marks)}
    for n in args.orders:
        if n not in inventories:
            continue
        profiles=inventories[n]
        feasible=[p for p in profiles if p[1]<=args.hi]
        def image(r,m):
            return {t for h,q,_,_ in (profiles if m==0 else feasible) if args.lo <= (t:=r*h+m*q) <= args.hi}
        scored={name:image(r,m) for name,(r,m) in readers.items()}
        scored['TT reader bank']=set.union(*(image(r,m) for r,m in marks.values()))
        scored['all lines'] = set.union(*(A for name,A in scored.items() if name!='H'))
        scored['all lines + direct H'] = scored['all lines']|scored['H']
        result[n]={'profiles':len(profiles),'primitive':sum(__import__('math').gcd(h,q)==1 for h,q,_,_ in profiles),'strong_profiles':sum(z!='-' for _,_,_,z in profiles),'scores':{name:band(A,args.lo,args.hi) for name,A in scored.items()}}
    print(json.dumps(result,indent=2,sort_keys=True))
