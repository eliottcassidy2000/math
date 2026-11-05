from pathlib import Path
import importlib.util,json
P=Path(__file__).resolve().parent; S=P/'source/session197-lowcycle';D=P/'source/session197/lowcycle-two-port'
def mod(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def need(b,m):
 if not b:raise RuntimeError(m)
c=mod('compiler',D/'compiler.py');t=mod('core',D/'twoport_core.py')
def constants(L):
 if L==1:return 1,1,1
 bits=[1-((L-1)>>i&1)for i in range((L-1).bit_length())];g=1+bits[0];w=0
 for i,b in enumerate(bits):
  if not b:g+=w
  if i and not bits[i-1] and b:g+=1
  w=2*w+b
 return g,1+sum(2**(len(bits)-1-i)for i,b in enumerate(bits)if b),1<<len(bits)
results=[];arcs=0
for x in json.loads((S/'arithmetic-hard-witnesses.json').read_text()):
 rows=[0];root=0;h,e,R=1,1,1;trace=[[h,e,R]]
 for L in x['recipe']:
  old=root;root=len(rows);bits=[1-((L-1)>>i&1)for i in range((L-1).bit_length())]
  n=len(rows);rows.extend([0]*(len(bits)+1))
  for a in range(n):
   for b in range(n,len(rows)):
    if b==root and a==old:rows[b]|=1<<a
    else:rows[a]|=1<<b
  for a in range(n,len(rows)):
   for b in range(a+1,len(rows)):
    if a==root and bits[b-root-1]==0:rows[b]|=1<<a
    else:rows[a]|=1<<b
  g,Q,B=constants(L);h,e,R=g*h-(L%2)*e+L*R,h-(L==1)*e,Q*h+B*R;trace.append([h,e,R])
 need(rows==x['rows']and root==x['root'],'independent physical adjacency')
 need(trace==x['profile_trace']and trace[-2]==x['prefix_state'],'exact trace and prefix')
 need(R>=h and h==x['target']==x['Hamilton_readout'] and len(rows)==x['order'],'target metadata')
 comp=c.compile_one(rows,{root});word=[[s['mask'],s['drop']]for s in comp['steps'][2:]]
 need(word==x['compiled_word'],'stored word differs from physical compiler')
 # Replay stored word into its literal birth-labelled graph, not just a scalar.
 out=[2,0];ports=[0,1];state=t.brute(out,ports);cycles=0
 for mask,drop in word:
  v=len(out);out.append(0)
  for u in range(v):
   if u in ports and mask>>(ports.index(u))&1:out[v]|=1<<u
   else:out[u]|=1<<v
  ports.append(v);del ports[drop]
  cycles=(0,state[1]+state[2],state[3]+state[5],state[2]+state[5])[mask]
  state=t.step(state,2,mask,drop)
 order=comp['birth_order']
 for i in range(len(rows)):
  for j in range(i+1,len(rows)):
   need((out[i]>>j&1)==(rows[order[i]]>>order[j]&1),'stored word physical arc');arcs+=1
 need(t.H(state,2)==h and cycles==1==x['cycle_readout'],'uncapped audited interface readouts')
 results.append(dict(target=h,order=len(rows),root=root,profile=trace[-1],word_length=len(word),cycles=cycles))
print(json.dumps(dict(status='PASS',hard_witnesses=results,independent_arc_checks=arcs),indent=2,sort_keys=True))
