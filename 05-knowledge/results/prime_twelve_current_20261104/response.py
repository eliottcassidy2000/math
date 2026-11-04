"""Retain complete named polynomials for four prime twelve-core controls."""
import argparse,json,subprocess,sys,tempfile
from pathlib import Path
from hashlib import sha256
from collections import Counter
HERE=Path(__file__).resolve().parent
OLDER=HERE.parent/'homogeneous_pair_response_20261103'
sys.path.insert(0,str(OLDER))
import check
import response as engine
need=check.need

def delete(rows,v):
 names=[u for u in range(len(rows)) if u!=v]
 return [sum(1<<j for j,w in enumerate(names) if rows[u]>>w&1) for u in names]

def boundary(rows):
 h,s,e,c=check.sidecar(rows)
 return [2*(e[u]-s[u])+sum(c[u])-sum(row[u] for row in c) for u in range(len(rows))]

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=HERE/'response.json');parser.add_argument('--tables',type=Path,default=HERE/'response_tables.json');args=parser.parse_args()
 paths={'response.py':HERE/'response.py','../homogeneous_pair_response_20261103/response.cpp':OLDER/'response.cpp',
        '../homogeneous_pair_response_20261103/response.py':OLDER/'response.py','../homogeneous_pair_response_20261103/check.py':OLDER/'check.py'}
 pins={k:sha256(p.read_bytes()).hexdigest() for k,p in paths.items()}
 orbit=json.loads((HERE/'switch.json').read_text());base=orbit['base_rows']
 cores=[]
 for flip in [1,3,11]:
  rows=[r^((flip^(4095 if flip>>u&1 else 0))&~(1<<u)) for u,r in enumerate(base)]
  need(orbit['records'][flip]['prime'],'prime orbit control')
  cores.append(('conference-switch-'+str(flip),rows))
 cores.append(('inherited-endpoint-balanced-prime',[3706,680,3603,1508,3914,1300,294,1397,519,1256,2370,1002]))
 ms=engine.masks(12);pos={m:i for i,m in enumerate(ms)};size=len(ms)**2
 high=[(4095,4095),(4095,0),(0,4095),(4094,4093),(1365,2730),(4063,4031)]
 qs=[(a,c) for a in ms for c in ms]+high
 inp=[str(2*len(cores))]
 for name,rows in cores:
  for a in [rows,check.opposite(rows)]:
   inp += [f'12 {len(qs)}',' '.join(map(str,a))]+[f'{x} {y}' for x,y in qs]
 with tempfile.TemporaryDirectory(prefix='prime-full-response-') as folder:
  exe=Path(folder)/'counter.exe'
  subprocess.run(['g++','-std=c++17','-O3',str(OLDER/'response.cpp'),'-o',str(exe)],check=True)
  out=subprocess.run([str(exe)],input='\n'.join(inp)+'\n',capture_output=True,text=True,check=True)
 vals=list(map(int,out.stdout.split()));need(len(vals)==len(cores)*2*len(qs),'query stream')
 tables=[];records=[]
 for j,(name,rows) in enumerate(cores):
  f=vals[2*j*len(qs):(2*j+1)*len(qs)];g=vals[(2*j+1)*len(qs):(2*j+2)*len(qs)]
  fc,gc=engine.coefficients(f[:size],ms),engine.coefficients(g[:size],ms)
  delta=[x-y for x,y in zip(fc,gc)];degrees=Counter()
  for k,x in enumerate(delta):
   a,c=ms[k//len(ms)],ms[k%len(ms)];d=a.bit_count()+c.bit_count()
   need(not x or d<=3,'whole-converse quartic')
   if d==2 and (not a or not c):need(x==0,'pure quadratic')
   if x:degrees[d]+=1
  cards=[];shared=[]
  for v in range(12):
   names=[u for u in range(12) if u!=v];card=delete(rows,v);b=boundary(card)
   hc,sc,ec,cc=check.sidecar(card)
   cards.append(dict(deleted=v,names=names,one_port_current=b,endpoint_boundary_parity=[(x-y)%2 for x,y in zip(ec,sc)]))
   for i,u in enumerate(names):
    actual=delta[pos[(1<<u)|(1<<v)]*len(ms)+pos[1<<v]]
    dual=delta[pos[1<<v]*len(ms)+pos[(1<<u)|(1<<v)]]
    need(actual==dual==-b[i],'exact deletion-current cubic')
    shared.append(dict(u=u,v=v,coefficient=actual))
  ht,st,et,ct=check.sidecar(rows)
  for v,card in enumerate(cards):
   rhs=1
   for i,u in enumerate(card['names']):
    if rows[v]>>u&1:rhs^=card['endpoint_boundary_parity'][i]
   need(rhs==(st[v]+et[v])%2,'apex boundary identity')
  for k,(a,c) in enumerate(high):
   need(engine.evaluate(fc,ms,a,c)==f[size+k] and engine.evaluate(gc,ms,a,c)==g[size+k],'full high-weight reconstruction')
  for a,c in [(1,0),(0,1),(3,2),(4094,4093)]:
   need(check.h(check.extend(rows,a,c))==engine.evaluate(fc,ms,a,c),'independent forward query')
  witness=next((dict(a=ms[k//len(ms)],c=ms[k%len(ms)],F=f[k],F_op=g[k]) for k in range(size) if (f[k]-g[k])%4),None)
  need(witness is not None,'actual mod4 separator')
  # Every proper subset tested: no isomorphism or score substitute for primality.
  modules=[]
  for m in range(1,4095):
   if m.bit_count()<2:continue
   if all(not(rows[u]&m) or (rows[u]&m)==m for u in range(12) if not m>>u&1):modules.append(m)
  need(not modules,'prime core')
  tables.append(dict(name=name,rows=rows,masks=ms,order='row-major (a,c)',values=f[:size],converse_values=g[:size],
                     coefficients=fc,converse_coefficients=gc,contrast_coefficients=delta,high_queries=[
                       dict(a=a,c=c,F=f[size+k],F_op=g[size+k]) for k,(a,c) in enumerate(high)]))
  records.append(dict(name=name,rows=rows,H=ht,nonzero_contrast_by_degree=dict(degrees),current=boundary(rows),
                      endpoint_boundary_parity=[(x-y)%2 for x,y in zip(et,st)],cards=cards,shared_cubics=shared,
                      mod4_query=witness,proper_subsets_checked=4082,full_battery=6241,cubic_battery=1753))
 tabletext=json.dumps(dict(status='FINITE-EXACT full named responses of four specified prime12 cores',cases=tables),sort_keys=True,separators=(',',':'))+'\n'
 args.tables.write_text(tabletext,encoding='utf8',newline='\n')
 need(pins=={k:sha256(p.read_bytes()).hexdigest() for k,p in paths.items()},'source changed')
 receipt=dict(status='PASS: finite prime-core controls, not universal exclusion',source_sha256=pins,cases=records,
              actual_compiled_queries=len(vals),independent_forward_queries=16,exact_shared_cubic_checks=1056,
              tables_sha256=sha256(tabletext.encode()).hexdigest())
 args.output.write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n',encoding='utf8',newline='\n')
 print(json.dumps(dict(actual_queries=len(vals),cores=len(cores),shared_checks=1056),indent=2))
if __name__=='__main__':main()

