from module_tools import *
import json
seeds=json.loads(Path(__file__).with_name('even_births.json').read_text())['seeds'];D={}
for T,O,q in seeds:D.setdefault(q,[]).append((tuple(T),O))
out=[]
for q,ls in D.items():
 if len(ls)<2:continue
 views=[]
 for T,O in ls:
  menu=[]
  for v in range(6):
   for mark in (['first','last'] if O>>v&1 else ['none']):
    Z,P=m.inflate(T,O,v,2,mark);h,w,u=m.data(Z,P)
    if h-w==1:
     Z,P=m.inflate(T,O,v,3,mark);h3,w3,u3=m.data(Z,P)
     menu.append((v,mark,u,h3-w3,u3))
  views.append({'rows':T,'O':O,'H':sum(m.endpoints(T)[-1]),'menu':menu})
 if len({str(sorted((x[2],x[3],x[4]) for x in r['menu'])) for r in views})>1:
  z={'q':q,'views':views};out.append(z);print(z,flush=True)
Path(__file__).with_name('seed_returns.json').write_text(json.dumps(out,indent=2))
