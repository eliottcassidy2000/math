import sys,json,random
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'session275_exact_selection_20261111/switches'))
from carrier import *
D={}
for T in tournaments(6):
 p=profile(T); hp.cache_clear();key=tuple(p[k] for k in ('h','Q','V','e'))
 D.setdefault(key,{})[tuple(p[k] for k in ('H','U','d'))]=T
rank=sorted(D.items(),key=lambda kv:len(kv[1]),reverse=True)
print(json.dumps([dict(key=k,values=[dict(p=p,T=t) for p,t in v.items()]) for k,v in rank[:40]],indent=2))
