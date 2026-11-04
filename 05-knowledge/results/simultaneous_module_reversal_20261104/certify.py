"""Complete named bidegree-(2,2) reversal cubes and exact isolation queries."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'homogeneous_pair_response_20261103'
sys.path.insert(0,str(OLD))
import check
import response as poly
need=check.need


def blow(q,cores):
 names=[];at=0
 for c in cores:names.append(list(range(at,at+len(c))));at+=len(c)
 a=[0]*at
 for i,c in enumerate(cores):
  for u,row in enumerate(c):a[names[i][u]]=sum(1<<names[i][v] for v in range(len(c)) if row>>v&1)
 for i in range(len(q)):
  for j in range(len(q)):
   if q[i]>>j&1:
    for u in names[i]:a[u]|=sum(1<<v for v in names[j])
 return a,names


def qprofile(a):
 n=len(a);full=(1<<n)-1
 h=[1]+[check.h([sum(1<<j for j,w in enumerate([i for i in range(n) if m>>i&1]) if a[u]>>w&1)
                  for u in range(n) if m>>u&1]) for m in range(1,full+1)]
 d=[[0]*(full+1) for _ in range(n+1)];d[0][0]=1
 for k in range(1,n+1):
  for m in range(1,full+1):
   b=m
   while b:d[k][m]+=h[b]*d[k-1][m^b];b=(b-1)&m
 return [d[k][full] for k in range(1,n+1)]


def separator(a):
 for x in poly.masks(len(a)):
  for y in poly.masks(len(a)):
   f=check.h(check.extend(a,x,y));g=check.h(check.extend(check.opposite(a),x,y))
   if (f-g)%4==2:return x,y,f,g
 raise RuntimeError('core separator')


def lift(q,names,i,x,y):
 a=sum(1<<names[i][u] for u in range(len(names[i])) if x>>u&1)
 c=sum(1<<names[i][u] for u in range(len(names[i])) if y>>u&1)
 for j,block in enumerate(names):
  if q[j]>>i&1:a|=sum(1<<v for v in block);c|=sum(1<<v for v in block)
 return a,c


def uniform(names,x,y):
 return (sum(1<<v for i,b in enumerate(names) if x>>i&1 for v in b),
         sum(1<<v for i,b in enumerate(names) if y>>i&1 for v in b))


def word_polynomial(a):
 n=len(a);out=Counter()
 for word in permutations(range(n+2)):
  terms={(0,0):1};valid=True
  for u,v in zip(word,word[1:]):
   if u<n and v<n:
    if not a[u]>>v&1:valid=False;break
   elif u>=n and v>=n:
    if (u,v)!=(n,n+1):valid=False;break
   else:
    old=u if u<n else v;port=v if u<n else u
    key=(1<<old,0) if port==n else (0,1<<old)
    if u<n:terms={(x|key[0],y|key[1]):z for (x,y),z in terms.items()}
    else:
     new=Counter(terms)
     for (x,y),z in terms.items():new[x|key[0],y|key[1]]-=z
     terms=dict(new)
  if valid:
   for key,value in terms.items():out[key]+=value
 return out


def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=HERE/'certificate.json')
 p.add_argument('--tables',type=Path,default=HERE/'full_response.json');args=p.parse_args()
 pins={f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['certify.py','../homogeneous_pair_response_20261103/check.py',
       '../homogeneous_pair_response_20261103/response.py','../homogeneous_pair_response_20261103/response.cpp']}
 families=[('two-pairs',[2,0],[[2,0],[2,0]]),
           ('cyclic-three-pairs',[2,4,1],[[2,0],[2,0],[2,0]]),
           ('cyclic-mixed',[2,4,1],[[2,0],[2,4,1],[6,4,0]])]
 jobs=[];states=[];fam_records=[]
 for name,q,cores in families:
  a,names=blow(q,cores);n=len(a);r=len(q);ms=poly.masks(n)
  battery=[(x,y) for x in ms for y in ms]
  unifs=[uniform(names,x,y) for x in range(1<<r) for y in range(1<<r)]
  seps=[separator(c) for c in cores]
  profiles=[qprofile(c) for c in cores]
  need(all(qprofile(check.opposite(c))==v for c,v in zip(cores,profiles)),'complete profile converse invariance')
  fr=dict(name=name,quotient=q,cores=cores,blocks=names,core_separators=seps,profiles=profiles,states=[])
  for top in range(2):
   qq=q if not top else check.opposite(q)
   probes=[lift(qq,names,i,*s[:2]) for i,s in enumerate(seps)]
   high=sorted(set(unifs+probes+[((1<<n)-1,(1<<n)-1),((1<<n)-1,0),(0,(1<<n)-1)]))
   qs=battery+high
   for eps in range(1<<r):
    cc=[check.opposite(c) if eps>>i&1 else c for i,c in enumerate(cores)]
    aa,_=blow(qq,cc)
    for block in names:
     mask=sum(1<<v for v in block)
     need(all((aa[u]&mask) in (0,mask) for u in range(n) if u not in block),'actual disjoint modules')
    index=len(states);fr['states'].append(index)
    states.append(dict(family=name,top=top,epsilon=eps,rows=aa,masks=ms,queries=qs,
                       battery=len(battery),probes=probes,uniform_queries=unifs))
    jobs.append((aa,qs))
  fam_records.append(fr)
 # Complete labelled quotient universe through four; every vertex becomes TT2.
 small=[]
 for r in range(1,5):
  for code in range(1<<(r*(r-1)//2)):
   q=check.tournament(r,code);base,names=blow(q,[[2,0]]*r)
   sep=separator([2,0]);probes=[lift(q,names,i,*sep[:2]) for i in range(r)]
   for eps in range(1<<r):
    a,_=blow(q,[check.opposite([2,0]) if eps>>i&1 else [2,0] for i in range(r)])
    small.append(dict(r=r,code=code,epsilon=eps,job=len(jobs),probes=probes));jobs.append((a,probes))
 inp=[str(len(jobs))]
 for a,qs in jobs:inp += [f'{len(a)} {len(qs)}',' '.join(map(str,a))]+[f'{x} {y}' for x,y in qs]
 with tempfile.TemporaryDirectory(prefix='simultaneous-response-') as folder:
  exe=Path(folder)/'counter.exe'
  subprocess.run(['g++','-O3','-std=c++17',str(OLD/'response.cpp'),'-o',str(exe)],check=True)
  out=subprocess.run([str(exe)],input='\n'.join(inp)+'\n',text=True,capture_output=True,check=True)
 stream=list(map(int,out.stdout.split()));need(len(stream)==sum(len(qs) for a,qs in jobs),'query stream')
 values=[];at=0
 for a,qs in jobs:values.append(stream[at:at+len(qs)]);at+=len(qs)
 tables=[];lookups=[];literal_poly=0;forward_counts=0
 for i,s in enumerate(states):
  ms=s['masks'];vals=values[i];co=poly.coefficients(vals[:s['battery']],ms)
  lookup=dict(zip(s['queries'],vals));lookups.append(lookup)
  for query,value in zip(s['queries'][s['battery']:],vals[s['battery']:]):need(poly.evaluate(co,ms,*query)==value,'full quartic reconstruction')
  if len(s['rows'])==4:
   direct=word_polynomial(s['rows'])
   need(all(co[j*len(ms)+k]==direct[x,y] for j,x in enumerate(ms) for k,y in enumerate(ms)),'literal full polynomial')
   literal_poly+=1
  for query in s['probes'][:1]:need(check.h(check.extend(s['rows'],*query))==lookup[query],'independent forward count');forward_counts+=1
  tables.append(dict(family=s['family'],top=s['top'],epsilon=s['epsilon'],rows=s['rows'],masks=ms,
                     order='row-major (a,c); every bidegree-(2,2) coefficient retained',
                     values=vals[:s['battery']],coefficients=co,
                     supplementary=[dict(a=q[0],c=q[1],F=v) for q,v in zip(s['queries'][s['battery']:],vals[s['battery']:])]))
 for fr in fam_records:
  r=len(fr['quotient']);ids=fr['states'];fr['decoders']=[]
  for top in range(2):
   base=next(i for i in ids if states[i]['top']==top and not states[i]['epsilon'])
   probes=states[base]['probes'];dec=[]
   for j,probe in enumerate(probes):
    unit=next(i for i in ids if states[i]['top']==top and states[i]['epsilon']==1<<j)
    first,second=lookups[base][probe],lookups[unit][probe]
    need((second-first)%4==2,'actual isolated nonzero residue')
    for i in ids:
     if states[i]['top']==top:need(lookups[i][probe]==(second if states[i]['epsilon']>>j&1 else first),'exact isolation')
    dec.append(dict(block=j,a=probe[0],c=probe[1],F_zero=first,F_one=second,difference=second-first))
   fr['decoders'].append(dict(top=top,coordinates=dec))
  need(len({tuple(tables[i]['coefficients']) for i in ids})==len(ids),'finite two-level full response injectivity')
  base=ids[0];fr['contrast_degrees']=[]
  for i in ids:
   contrast=[x-y for x,y in zip(tables[i]['coefficients'],tables[base]['coefficients'])]
   tables[i]['contrast_from_top0_epsilon0']=contrast
   counts=Counter();modcounts=Counter();ms=states[i]['masks']
   for k,z in enumerate(contrast):
    degree=ms[k//len(ms)].bit_count()+ms[k%len(ms)].bit_count()
    if z:counts[degree]+=1
    if z%4:modcounts[degree]+=1
   fr['contrast_degrees'].append(dict(top=states[i]['top'],epsilon=states[i]['epsilon'],nonzero_by_degree=dict(counts),mod4_nonzero_by_degree=dict(modcounts)))
  top1=next(i for i in ids if states[i]['top']==1 and states[i]['epsilon']==0)
  uq=states[base]['uniform_queries'];diff=[lookups[top1][q]-lookups[base][q] for q in uq]
  fr['quotient_only_uniform']=dict(assignments=len(uq),nonzero_exact=sum(x!=0 for x in diff),nonzero_mod4=sum(x%4!=0 for x in diff),differences=diff)
  # A complete mixed flip difference, retaining all higher-degree interaction.
  ids0={states[i]['epsilon']:i for i in ids if not states[i]['top']}
  interaction=[a-b-c+d for a,b,c,d in zip(tables[ids0[0]]['coefficients'],tables[ids0[1]]['coefficients'],tables[ids0[2]]['coefficients'],tables[ids0[3]]['coefficients'])]
  fr['two_flip_interaction_coefficients']=interaction
  fr['interaction_nonzero']=sum(x!=0 for x in interaction)
  for probe in states[base]['probes']:need(poly.evaluate(interaction,states[base]['masks'],*probe)==0,'interaction vanishes on decoder')
 grouped={}
 for rec in small:grouped.setdefault((rec['r'],rec['code']),[]).append(rec)
 for group in grouped.values():
  r=group[0]['r'];by={x['epsilon']:values[x['job']] for x in group}
  for j in range(r):
   zero,one=by[0][j],by[1<<j][j];need((one-zero)%4==2,'complete quotient decoder residue')
   need(all(vals[j]==(one if eps>>j&1 else zero) for eps,vals in by.items()),'complete quotient isolation')
 for fr in fam_records:
  if fr['name']=='cyclic-three-pairs':need(fr['quotient_only_uniform']['nonzero_mod4']==0 and fr['quotient_only_uniform']['nonzero_exact']>0,'uniform mod4 hostile')
 text=json.dumps(dict(status='FINITE-EXACT complete named responses; no cubic or quartic coordinates omitted',states=tables),sort_keys=True,separators=(',',':'))+'\n'
 args.tables.write_text(text,encoding='utf8',newline='\n')
 need(pins=={f:sha256((HERE/f).read_bytes()).hexdigest() for f in pins},'source changed')
 result=dict(status='PASS disjoint cube isolation; nested examples finite only',source_sha256=pins,
             actual_compiled_queries=len(stream),full_response_states=len(states),literal_complete_polynomials=literal_poly,
             independent_forward_counts=forward_counts,complete_labelled_quotients=len(grouped),small_reversal_states=len(small),
             small_actual_queries=sum(rec['r'] for rec in small),families=fam_records,tables_sha256=sha256(text.encode()).hexdigest())
 args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf8',newline='\n')
 print(json.dumps({k:v for k,v in result.items() if k not in ('source_sha256','families')},indent=2))


if __name__=='__main__':main()
