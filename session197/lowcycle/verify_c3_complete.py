#!/usr/bin/env python3
"""Exact low-cycle census wrapper. Requires a standard C++17 compiler (g++).
All graph enumeration/counting code is in the adjacent auditable C++ source.
No package, network, repository checkout or external graph bank is used.
"""
from pathlib import Path
import hashlib,json,subprocess,tempfile,collections,itertools
P=Path(__file__).resolve().parent

def need(ok,why):
 if not ok:raise RuntimeError(why)

def main():
 expected1=[1,3,8,21,55,144,377,987,2584]
 expected2=[0,2,5,18,60,191,589,1775,5257]
 expected3=[0,1,11,38,124,390,1193,3574,10537]
 with tempfile.TemporaryDirectory(prefix='c3-lowcycle-') as d:
  d=Path(d);exe=d/'verify';records=d/'records.csv'
  subprocess.run(['g++','-O2','-std=c++17',str(P/'verify_c3_complete.cpp'),'-o',str(exe)],check=True)
  run=subprocess.run([str(exe),str(records)],check=True,text=True,capture_output=True)
  need(not run.stderr,'unexpected census stderr')
  inventory=[json.loads(line)for line in run.stdout.splitlines()];need(len(inventory)==9,'missing census order')
  byorder=collections.defaultdict(list);labelled=collections.Counter()
  for line in records.read_bytes().splitlines(keepends=True):
   z=list(map(int,line.split(b',')));n,c=z[:2];r=z[2:]
   need(len(r)==n and 1<=c<=3,'record shape')
   need(all(not(r[i]>>i&1)and 0<=r[i]<(1<<n)for i in range(n)),'diagonal/range')
   need(all(((r[i]>>j&1)+(r[j]>>i&1))==1 for i in range(n)for j in range(i+1,n)),'tournament record')
   if c==3:need(any(x.bit_count() in (1,n-2)for x in r),'C3 degree-one gate')
   byorder[n].append(line);labelled[n,c]+=1
  for i,entry in enumerate(inventory):
   n=i+4;need(entry['order']==n,'order sequence')
   for c,expected in ((1,expected1),(2,expected2),(3,expected3)):
    need(entry['C'+str(c)]==expected[i]==labelled[n,c],('class inventory',n,c))
   need(entry['no_degree_one_C3']==0,'C3 atom')
   entry['canonical_graph_sha256']=hashlib.sha256(b''.join(byorder[n])).hexdigest()
  need(sum(len(x)for x in byorder.values())==sum(sum(e['C'+str(c)]for c in (1,2,3))for e in inventory),'total records')
  # Exact degree-cover arithmetic, including the empty low-degree set.
  bounds=[l*(7-l)//2 for l in range(6)]
  need(bounds==[0,3,5,6,6,5] and max(bounds)*2==12,'two-sided atom bound')
  hostile=[14,52,56,50,33,1]
  hc=sum(all(hostile[u]>>v&1 for u,v in zip((0,)+w,w+(0,))) for w in itertools.permutations(range(1,6)))
  degrees=[x.bit_count() for x in hostile]
  need(hc==3 and degrees.count(1)==1 and 4 not in degrees,'one-sided degree hostile')
  payload={'status':'PASS','one_sided_hostile':{'rows':hostile,'C':hc,'outdegrees':degrees},'maximum_order':12,'inventory':inventory,'total_canonical_graphs':sum(len(x)for x in byorder.values()),'all_records_sha256':hashlib.sha256(records.read_bytes()).hexdigest(),'two_sided_cover_bounds':bounds,'minimal_bridge_free_C3_order_bound':12,'all_order_theorem':'Every strong tournament with exactly three Hamiltonian cycles has a strong bridge, and hence an indegree-one or outdegree-one vertex.'}
 data=json.dumps(payload,sort_keys=True,indent=2)+'\n';(P/'c3_complete_certificate.json').write_text(data);print(data,end='')
if __name__=='__main__':main()
