"""Audit profile counts and gaps in the cap sweep."""
from score import read, mark_w, band
import json, pathlib, sys
path=pathlib.Path(sys.argv[1])
lo=int(sys.argv[2]) if len(sys.argv)>2 else 69
hi=59999
reader_bank={m:(1-mark_w(m),m) for m in range(1,64)}
families = {"Q+H":(1,1),"2Q+H":(1,2),"4Q+H":(1,4),"8Q+H":(1,8),"3Q-H":(-1,3),"6Q-H":(-1,6),"1024Q-H":(-1,1024)}
fixed_scc={}
ext_scc={}
for n in range(1,16):
    f=path/f'pairline_1_{n}.tsv'
    if f.is_file():
        x=read(f)[n]
        print(n,len(x),*[sum(lo<=a*h+m*q<=hi for h,q,_,_ in x) for a,m in families.values()])
        for a,m in [*families.values(),*reader_bank.values()]:
            akey=(a,m)
            ext_scc.setdefault(akey,set()).update(a*h+m*q for h,q,_,_ in x if lo<=a*h+m*q<=hi)
alllow=read(path/'pairline10.tsv')
for x in alllow.values():
    x=[(h,q) for h,q,_,_ in x if q < 59999]
    for a,m in [*families.values(),*reader_bank.values()]:
        fixed_scc.setdefault((a,m),set()).update(a*h+m*q for h,q in x if lo<=a*h+m*q<=hi)
u=set.union(*fixed_scc.values(),*ext_scc.values())
old=set.union(*fixed_scc.values())
print(json.dumps({'fixed':{k:band(fixed_scc[v],lo,hi) for k,v in families.items()},
    'all':{k:band(fixed_scc[v]|ext_scc.get(v,set()),lo,hi) for k,v in families.items()},
    'expanded_union':band(u,lo,hi),'orders_up_to_10_union':band(old,lo,hi),'new_hits':sorted(u-old)},indent=2))
