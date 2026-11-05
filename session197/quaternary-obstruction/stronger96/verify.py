from pathlib import Path
import hashlib,json,subprocess,tempfile
P=Path(__file__).resolve().parent

def need(v,s):
 if not v:raise RuntimeError(s)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def run():
 expected='fd502c3c82507fe72991f2db9bcd56f3c399fc4d39b005df557cb627562a7a25';need(sha(P/'strong-bank.jsonl')==expected,'frozen target bank')
 with tempfile.TemporaryDirectory(prefix='quaternary96-')as d:
  d=Path(d)
  for name in ['generate','clone_counts','gap_counts']:subprocess.run(['g++','-O2','-std=c++17',str(P/(name+'.cpp')),'-o',str(d/name)],check=True)
  bank=d/'bank.jsonl';generation=subprocess.check_output([str(d/'generate'),str(bank)],text=True);need(sha(bank)==expected,'entire regenerated complete bank')
  records=[json.loads(s)for s in bank.read_text().splitlines()];hist={}
  for x in records:hist[x['n']]=hist.get(x['n'],0)+1
  need([hist[n]for n in range(3,12)]==[1,1,6,35,353,3343,3153,1517,331],'complete order inventory')
  csv=d/'input.csv';csv.write_text(''.join(','.join(map(str,[x['n'],x['H']]+x['rows']))+'\n'for x in records))
  literal=subprocess.check_output([str(d/'clone_counts'),str(csv)],text=True).strip();gaps=subprocess.check_output([str(d/'gap_counts'),str(csv)],text=True).strip()
 need(literal=='PASS graphs 8740 eligible_graph_signs 1398 roots 12263 m2_survivors 6 m3_survivors 0 m4_survivors 0','literal root inventory')
 need(gaps=='DONE rank3_roots 13456 hits 0','gap root inventory')
 need(361**3==47045881 and 5**11==48828125 and 361**3<5**11,'strict order11 cutoff')
 need(max(A for A in range(1,96)if A%6==0)==90,'factorial amplitude bound')
 print(json.dumps({'status':'PASS','bank_sha256':expected,'bank_classes':len(records),'order_counts':hist,'generation':generation.splitlines(),'literal_clone_counter':literal,'exterior_gap_counter':gaps,'cutoff':{'H_max':361,'maximum_strong_order':11,'H_cubed':361**3,'five_to_eleven':5**11},'excluded_amplitudes':list(range(6,91,6)),'offsets':[-1,1],'necessary_amplitude':96,'source_hashes':{name:sha(P/name)for name in ['generate.cpp','clone_counts.cpp','gap_counts.cpp']}},sort_keys=True,indent=2))
if __name__=='__main__':run()
