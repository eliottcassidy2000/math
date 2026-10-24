#!/usr/bin/env python3
import importlib.util,json,random,sys
from pathlib import Path
from reference import need,stable,autos,error,relabel,strong
sys.dont_write_bytecode=True
spec=importlib.util.spec_from_file_location('decoder','source/decoder.py');D=importlib.util.module_from_spec(spec);spec.loader.exec_module(D)
ROWS=((102,108,408,657,771,220,412,787,555,103),(102,108,664,401,771,220,412,787,551,107))

def subset_paths(rows):
 n=len(rows);full=(1<<n)-1;endcounts=[0]*n
 # A scalar DP for each distinguished start, hence retaining both endpoints.
 matrix=[]
 for start in range(n):
  levels={(1<<start,start):1}
  for size in range(1,n):
   nxt={}
   for (mask,end),count in levels.items():
    choices=rows[end]&~mask
    while choices:
     bit=choices&-choices;choices^=bit;key=(mask|bit,bit.bit_length()-1);nxt[key]=nxt.get(key,0)+count
   levels=nxt
  values=[levels.get((full,v),0) for v in range(n)];matrix.append(values)
  endcounts=[a+b for a,b in zip(endcounts,values)]
 return sum(endcounts),matrix,endcounts

def literal_paths(rows):
 n=len(rows);full=(1<<n)-1;ends=[0]*n;nodes=0
 def visit(last,used):
  nonlocal nodes
  nodes+=1
  if used==full:ends[last]+=1;return
  rest=rows[last]&~used
  while rest:
   bit=rest&-rest;rest^=bit;visit(bit.bit_length()-1,used|bit)
 for start in range(n):visit(start,1<<start)
 return sum(ends),ends,nodes

def main():
 records=[];decoded=[];rng=random.Random(51)
 for rows in ROWS:
  need(all(not(r>>i&1) and r<1<<10 for i,r in enumerate(rows)),'tournament bounds')
  need(all((rows[i]>>j&1)+(rows[j]>>i&1)==1 for i in range(10) for j in range(i)),'tournament arcs')
  dec=D.decode(rows);group=autos(rows);e,orbits=error(rows,group);need(group==(tuple(range(10)),),'full ordinary rigidity');need(e==40==dec['error'],'refinement error');need(strong(rows),'strong witnesses')
  need(stable(rows)==(tuple(range(5)),tuple(range(5,10))),'actual stable cells');need(all(sum(rows[x]>>y&1 for y in range(5,10))==2 for x in range(5)),'degree two crossing')
  H,pairs,ends=subset_paths(rows);literal,literal_ends,nodes=literal_paths(rows);need(H==literal and ends==literal_ends,'DP and literal path enumeration')
  converse=tuple(1023^r^(1<<i) for i,r in enumerate(rows));need(subset_paths(converse)[0]==H,'converse invariance')
  perm=list(range(10));rng.shuffle(perm);need(subset_paths(relabel(rows,perm))[0]==H,'relabel invariance')
  decoded.append(dec);records.append({'rows':rows,'hamiltonian_paths':H,'endpoint_matrix':pairs,'terminal_counts':ends,'literal_prefixes_enumerated':nodes,'automorphism_order':len(group),'vertex_orbits':orbits,'error':e})
 need(decoded[0]==decoded[1],'entire decoder outputs equal');need([r['hamiltonian_paths'] for r in records]==[15465,15549],'literal witness totals')
 print(json.dumps({'status':'PASS','whole_decoded_outputs_equal':True,'common_decoder_output':decoded[0],'path_mass_difference':records[1]['hamiltonian_paths']-records[0]['hamiltonian_paths'],'witnesses':records},indent=2,sort_keys=True))
if __name__=='__main__':main()
