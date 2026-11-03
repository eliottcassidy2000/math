import itertools,json,re,pathlib
checks=0
def require(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
def partition(labels):
 return tuple(tuple(i for i,x in enumerate(labels) if x==v) for v in dict.fromkeys(labels))
def refine(labels,maps):
 sig=[(labels[i],)+tuple(labels[f[i]] for f in maps) for i in range(len(labels))]
 names={v:j for j,v in enumerate(dict.fromkeys(sig))}
 return tuple(names[v] for v in sig)
cases=0
for n in range(1,4):
 maps=list(itertools.product(range(n),repeat=n))
 for fs in itertools.product(maps,repeat=2):
  for obs in itertools.product(range(2),repeat=n):
   current=obs;hist=[]
   for k in range(n+1):
    nxt=refine(current,fs);q=len(set(current));r=len(set(nxt));hist.append(q)
    require(r>=q,'monotonicity')
    if r==q:
     require(partition(current)==partition(nxt),'equal counts equal relation')
     for i,j in itertools.combinations(range(n),2):
      if current[i]==current[j]:
       require(all(current[f[i]]==current[f[j]] for f in fs),'congruence')
     require(partition(nxt)==partition(refine(nxt,fs)),'permanent stability')
    current=nxt
   require(len(set(hist[-2:]))==1,'finite stabilization')
   cases+=1
text=pathlib.Path(__file__).with_name('NOTE.md').read_text()
links=re.findall(r'\]\((\.\./[^)]+)\)',text)
for link in links:
 part=link.split('/')[1]
 require((pathlib.Path(__file__).parent / link).is_file(),'declared sibling target')
print(json.dumps({'status':'pass','finite_systems':cases,'checks':checks,'relative_link_occurrences':len(links),'scope':'all 1-3-state two-map systems with binary observations; links verified against intended sibling packaging'},sort_keys=True,indent=2))
