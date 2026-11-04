"""Integer covers, exact card currents, literal queries, and full named cubics."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'homogeneous_pair_response_20261103'
sys.path.insert(0,str(OLD))
import check
need=check.need


def delete(a,v):
 names=[u for u in range(len(a)) if u!=v]
 return names,[sum(1<<j for j,w in enumerate(names) if a[u]>>w&1) for u in names]


def current(a):
 _,s,e,c=check.sidecar(a)
 return [2*(e[u]-s[u])+sum(c[u])-sum(row[u] for row in c) for u in range(len(a))]


def star(a,v):
 d=a[v].bit_count()
 return [((1<<u)|((1<<v) if a[v]>>u&1 else 0),1<<v)
         for u in range(len(a)) if u!=v]+[((1<<v) if d%2 else 0,1<<v)]


def exact(a):
 n=len(a);h,s,e,c=check.sidecar(a)
 need(all((c[u][v]-c[v][u])%2==0 for u in range(n) for v in range(n)), 'half-cover integrality')
 d=[[(c[u][v]-c[v][u])//2%2 for v in range(n)] for u in range(n)]
 boundary=[(x+y)%2 for x,y in zip(s,e)]
 cards=[]
 if n>=3:
  for v in range(n):
   names,q=delete(a,v);g=current(q)
   need(all(x%2==0 for x in g),'half-current integrality')
   parity=(sum(d[u][v] for u in range(n))+sum(g[i]//2 for i,u in enumerate(names) if a[v]>>u&1))%2
   need(parity==1,'integer cover/current obstruction')
   cards.append(dict(deleted=v,names=names,current=g,parity=parity))
 if n%2 and n>=3:
  need(sum(d[u][v] for u in range(n) for v in range(u+1,n))%2==1,'integer odd cover obstruction')
 return dict(rows=a,H=h,S=s,E=e,C=c,D=d,endpoint_boundary=boundary,cards=cards)


def evaluate(coeff,a,c):
 return sum(v for x,y,v in coeff if x&a==x and y&c==y)


def cyclic_quotient(a):
 n=len(a);out=[r|(1<<n) for r in a]+[1<<(n+1),(1<<n)-1]
 return out


def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=HERE/'certificate.json')
 p.add_argument('--tables',type=Path,default=HERE/'full_named_cubic.json');args=p.parse_args()
 sources=['certificate.py','../homogeneous_pair_response_20261103/check.py',
          '../homogeneous_pair_response_20261103/response.cpp',
          '../prime_twelve_current_20261104/response_tables.json',
          '../session187_named_response_20261104/prime12/cubic.json']
 pins={f:sha256((HERE/f).read_bytes()).hexdigest() for f in sources}
 # Independent small integer and literal-word controls, with no class-bank input.
 small=0;literal=0
 for n in range(1,6):
  for code in range(1<<(n*(n-1)//2)):
   a=check.tournament(n,code);rec=exact(a);small+=1
   if n>4:continue
   need(check.literal_sidecar(a)==(rec['H'],rec['S'],rec['E'],rec['C']), 'literal cover counts')
   if n==1:
    need(a==check.opposite(a),'singleton positive equality control');continue
   if n==2:
    queries=[(1,0)]
   elif n%2==0:
    queries=star(a,0)+[(1,0)]
   else:
    u,v=next((u,v) for u in range(n) for v in range(u+1,n) if rec['D'][u][v])
    queries=[(0,0),(1<<u,0),(0,1<<v),(1<<u,1<<v)]
   deltas=[]
   for x,y in queries:
    f=check.literal_h(check.extend(a,x,y));g=check.literal_h(check.extend(check.opposite(a),x,y))
    need(f==check.h(check.extend(a,x,y)) and g==check.h(check.extend(check.opposite(a),x,y)), 'literal vs forward DP')
    literal+=2;deltas.append(f-g)
   if n>=4:
    need(sum(deltas[:-1])%4==2*(1+rec['endpoint_boundary'][0])%4,'literal star identity')
    need(deltas[-1]%4==2*rec['endpoint_boundary'][0],'literal singleton identity')
   need(any(x%4==2 for x in deltas),'literal actual separator')
 # Preserve every allowed integer coefficient, not merely the new certificate.
 inherited=json.loads((HERE/sources[3]).read_text())['cases']
 tb=next(t for t in inherited if t['name']=='inherited-endpoint-balanced-prime')
 ms=tb['masks']
 bc=[[x,y,tb['contrast_coefficients'][i*len(ms)+j]] for i,x in enumerate(ms) for j,y in enumerate(ms)
     if x.bit_count()+y.bit_count()<=3 and not ((x==0 and y.bit_count()==2) or (y==0 and x.bit_count()==2))]
 tp=json.loads((HERE/sources[4]).read_text())
 pc=[[x,y,v] for x,y,v in tp['coefficients'] if not ((x==0 and y.bit_count()==2) or (y==0 and x.bit_count()==2))]
 need(len(bc)==len(pc)==1753,'complete allowed named cubic')
 for x,y,v in tp['coefficients']:
  if (x==0 and y.bit_count()==2) or (y==0 and x.bit_count()==2):need(v==0,'inherited pure quadratic zero')
 polynomials=[dict(name='prime-B',rows=tb['rows'],coefficients=bc),dict(name='prime-P',rows=tp['rows'],coefficients=pc)]
 base=[sum(1<<v for v in range(11) if (v-u)%11 in [1,3,4,5,9]) for u in range(11)]+[2047]
 c3pairs=[sum(1<<v for v in range(6) if v//2==(u//2+1)%3 or (v//2==u//2 and u<v)) for u in range(6)]
 cores=[('TT2',[2,0]),('C3',[2,4,1]),('TT4',[14,12,8,0]),('C3-pairs',c3pairs),
        ('balanced-eight',[4,17,42,131,45,75,31,119]),('Paley11-source',base),
        ('prime-B',tb['rows']),('prime-P',tp['rows'])]
 cases=[];jobs=[];query_sets=[]
 for name,a in cores:
  rec=exact(a);n=len(a);rec['name']=name
  if n%2==0 and n>=4:
   queries=sorted(set(q for v in range(n) for q in star(a,v)+[(1<<v,0)]))
  elif n==2:queries=[(1,0)]
  else:queries=[(0,0),(1,0),(0,2),(1,2)]
  query_sets.append(queries);cases.append(rec)
  jobs += [(a,queries),(check.opposite(a),queries)]
 # Whole core versus internal module converse in the actual strong C3 quotient.
 ambient_specs=[]
 for name,a in cores[-2:]:
  n=len(a);queries=star(a,0);lift=[(x|(1<<(n+1)),y|(1<<(n+1))) for x,y in queries]
  t,u=cyclic_quotient(a),cyclic_quotient(check.opposite(a))
  ambient_specs.append(dict(name=name,rows=t,reversed_module_rows=u,module=list(range(n)),queries=lift,core_queries=queries))
  jobs += [(t,lift),(u,lift)]
 # A positive same-module control is actually counted twice in the same context.
 positive=cyclic_quotient(c3pairs);positive_queries=[(0,0),(7,5),(255,255)]
 jobs += [(positive,positive_queries),(positive,positive_queries)]
 inp=[str(len(jobs))]
 for a,qs in jobs:inp += [f'{len(a)} {len(qs)}',' '.join(map(str,a))]+[f'{x} {y}' for x,y in qs]
 with tempfile.TemporaryDirectory(prefix='parity-certificates-') as folder:
  exe=Path(folder)/'counter.exe'
  subprocess.run(['g++','-O3','-std=c++17',str(OLD/'response.cpp'),'-o',str(exe)],check=True)
  out=subprocess.run([str(exe)],input='\n'.join(inp)+'\n',text=True,capture_output=True,check=True)
 values=list(map(int,out.stdout.split()));need(len(values)==sum(len(q) for _,q in jobs),'actual query stream')
 responses=[];at=0
 for _,qs in jobs:responses.append(values[at:at+len(qs)]);at+=len(qs)
 full_polynomial_queries=0;forward_queries=0
 for i,(name,a) in enumerate(cores):
  n=len(a);rec=cases[i];qs=query_sets[i];f,g=responses[2*i:2*i+2]
  vals={q:(x,y) for q,x,y in zip(qs,f,g)}
  rec['queries']=[dict(a=q[0],c=q[1],F=x,F_op=y) for q,x,y in zip(qs,f,g)]
  rec['stars']=[]
  for v in range(n) if n%2==0 and n>=4 else []:
   total=sum(vals[q][0]-vals[q][1] for q in star(a,v))
   need(total%4==2*(1+rec['endpoint_boundary'][v])%4,'compiled actual star identity')
   single=vals[(1<<v,0)];need((single[0]-single[1])%4==2*rec['endpoint_boundary'][v],'compiled singleton identity')
   rec['stars'].append(dict(anchor=v,queries=[list(q) for q in star(a,v)],exact_contrast_sum=total,residue=total%4))
  witness=next((dict(a=q[0],c=q[1],F=x,F_op=y) for q,x,y in zip(qs,f,g) if (x-y)%4==2),None)
  need(witness is not None,'core actual separator');rec['witness']=witness
  for q in qs[:2]:
   need(check.h(check.extend(a,*q))==vals[q][0] and check.h(check.extend(check.opposite(a),*q))==vals[q][1], 'forward DP core query')
   forward_queries+=2
  if name in ('prime-B','prime-P'):
   poly=next(t for t in polynomials if t['name']==name);coeff=poly['coefficients'];cd={(x,y):z for x,y,z in coeff}
   need(not any(rec['endpoint_boundary']),'endpoint balanced prime control')
   for q,(x,y) in vals.items():need(evaluate(coeff,*q)==x-y,'complete named cubic evaluated on actual query');full_polynomial_queries+=1
   for card in rec['cards']:
    v=card['deleted']
    for u,gvalue in zip(card['names'],card['current']):
     need(cd[(1<<u)|(1<<v),1<<v]==cd[1<<v,(1<<u)|(1<<v)]==-gvalue,'exact saved deletion cubic')
   need(all(cd[1<<u,1<<v]//2%2==rec['D'][u][v] for u in range(n) for v in range(n)), 'saved mixed quadratic')
 # The two hostile mechanisms fail separately, exactly as the proof requires.
 cp=next(r for r in cases if r['name']=='C3-pairs')
 need(not any(x for row in cp['D'] for x in row) and any(x%4 for c in cp['cards'] for x in c['current']), 'cover-only hostile')
 ps=next(r for r in cases if r['name']=='Paley11-source')
 need(not any(x%4 for c in ps['cards'] for x in c['current']) and any(x for row in ps['D'] for x in row), 'shared-cubic-only hostile')
 for i,spec in enumerate(ambient_specs):
  f,g=responses[2*len(cores)+2*i:2*len(cores)+2*i+2]
  core=next(r for r in cases if r['name']==spec['name']);lookup={(q['a'],q['c']):q['F']-q['F_op'] for q in core['queries']}
  spec['values']=[dict(a=q[0],c=q[1],F=x,F_reversed=y) for q,x,y in zip(spec['queries'],f,g)]
  need(all((x-y-3*lookup[q])%4==0 for q,x,y in zip(spec['core_queries'],f,g)),'actual strong quotient residue transfer')
  need(sum(x-y for x,y in zip(f,g))%4==2,'ambient star certificate')
 need(responses[-1]==responses[-2],'same-module positive equality control')
 tabletext=json.dumps(dict(status='FINITE-EXACT full named integer converse cubics; all allowed zeros retained',cases=polynomials),sort_keys=True,separators=(',',':'))+'\n'
 args.tables.write_text(tabletext,encoding='utf8',newline='\n')
 need(pins=={f:sha256((HERE/f).read_bytes()).hexdigest() for f in pins},'sources changed')
 result=dict(status='PASS integer and actual-query controls for E1-E4',source_sha256=pins,
             small_labelled_integer_controls=small,literal_hamilton_counts=literal,
             actual_compiled_counts=len(values),independent_forward_counts=forward_queries,
             full_cubic_actual_queries=full_polynomial_queries,cases=cases,ambient=ambient_specs,
             positive_control_counts=responses[-1],tables_sha256=sha256(tabletext.encode()).hexdigest())
 args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf8',newline='\n')
 print(json.dumps({k:v for k,v in result.items() if k not in ('cases','ambient','source_sha256')},indent=2))


if __name__=='__main__':main()
