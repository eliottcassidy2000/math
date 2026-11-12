"""Query exact all-order scalar fibers for certified core universes."""
import argparse,json
from receiver_fiber import TRI,S4,A5,E6,BC,C_SC_T,fiber
p=argparse.ArgumentParser();p.add_argument('receiver',type=int);p.add_argument('--universe',choices=['h3','h5','nine-bank-cores'],default='h5');args=p.parse_args()
cores={'h3':[('triangle',TRI)],'h5':[('strong4',S4)],'nine-bank-cores':[('A',A5),('E',E6),('BC',BC),('C-singleton-B',C_SC_T)]}[args.universe]
rows=[]
for name,D in cores:
 for w in fiber(D,args.receiver):rows.append(dict(core=name,**w))
support=sorted({w['profile']['H'] for w in rows});runs=[]
for H in support:
 if H-2 in support:continue
 end=H
 while end+2 in support:end+=2
 runs.append(dict(first=H,last=end,length=(end-H)//2+1))
print(json.dumps(dict(V=args.receiver,universe=args.universe,scope='GLOBAL at h3/h5 by certified core classification; otherwise exactly the named core/padding universe',support=support,runs=runs,witnesses=rows),indent=2))
