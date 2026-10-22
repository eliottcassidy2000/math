import collections,json,math

def machine(m,v=(1,2),invertible=True,observe=True,run_cap=1):
 start=(0,0,0); sink=None; states=[start]; ids={start:0}; trans=[]; outputs=[]
 for s in states:
  if s is None: ns=[None,None];o='invalid'
  else:
   x,y,last=s;o=x if observe else 0
   ns=[]
   for d in [0,1]:
    if d and last==run_cap: ns.append(None)
    else:
     a,b=(y,x+y) if invertible else (0,0)
     ns.append(((a+d*v[0])%m,(b+d*v[1])%m,(last+1 if d else 0)))
  row=[]
  for t in ns:
   if t not in ids:ids[t]=len(states);states.append(t)
   row.append(ids[t])
  trans.append(row); outputs.append(o)
 # independent partition refinement by behavioral signatures
 groups={};part=[]
 for o in outputs:
  if o not in groups:groups[o]=len(groups)
  part.append(groups[o])
 while True:
  groups={};new=[]
  for i,o in enumerate(outputs):
   key=(o,part[trans[i][0]],part[trans[i][1]])
   if key not in groups:groups[key]=len(groups)
   new.append(groups[key])
  if new==part:break
  part=new
 return {'modulus':m,'reachable':len(states),'minimal':len(set(part))}
res={'fibonacci':[],'lucas':[],'hostiles':{}}
for m in range(1,41):
 for name,v,formula in [('fibonacci',(1,2),2*m*m+1),('lucas',(2,1),2*m*m//math.gcd(m,5)+1)]:
  r=machine(m,v);assert r['minimal']==r['reachable']==formula,(name,r,formula);res[name].append(r)
res['hostiles']['noninvertible_M_zero']=machine(5,(1,0),False)
res['hostiles']['unobservable_output']=machine(5,(1,2),True,False)
# first-round scalar compression hostile from inherited 2 vs 8
res['hostiles']['scalar_state_2_vs_8']={'modulus':3,'prefixes':['10','10000'],'value_residues':[2,2],'shift_residues':[0,1],'after_zero_values':[3,13]}
res['bounded_runs']=[]
for cap in range(2,5):
 for m in range(1,21):
  for name,v,index in [('fibonacci',(1,2),1),('lucas',(2,1),math.gcd(m,5))]:
   r=machine(m,v,run_cap=cap);assert r['minimal']==r['reachable']==(cap+1)*m*m//index+1
   res['bounded_runs'].append(dict(r,run_cap=cap,family=name))
print(json.dumps(res,indent=2))
