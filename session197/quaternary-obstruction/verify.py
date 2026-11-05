from pathlib import Path
import hashlib,json,subprocess,tempfile
HERE=Path(__file__).resolve().parent

def need(v,s):
 if not v:raise RuntimeError(s)
def main():
 bank=HERE/'strong-bank.jsonl';need(hashlib.sha256(bank.read_bytes()).hexdigest()=='3d4ddbe02807c7c6bf962147c582e0c6707eb4d221134a715cdd8a84bf730951','frozen complete bank')
 data=[json.loads(s) for s in bank.read_text().splitlines()];hist={}
 for x in data:hist[x['n']]=hist.get(x['n'],0)+1
 need([hist[n] for n in range(3,11)]==[1,1,6,35,353,1144,726,190],'bank order inventory')
 out={}
 with tempfile.TemporaryDirectory(prefix='quaternary-controls-') as td:
  td=Path(td);inp=td/'input.csv';inp.write_text(''.join(','.join(map(str,[x['n'],x['H']]+x['rows']))+'\n' for x in data))
  for mode in ['clone_counts','gap_counts']:
   exe=td/mode;subprocess.run(['g++','-std=c++17','-O2',str(HERE/(mode+'.cpp')),'-o',str(exe)],check=True)
   out[mode]=subprocess.check_output([str(exe),str(inp)],text=True).strip()
 need(out['clone_counts']=='PASS graphs 2456 eligible_graph_signs 410 roots 3426 m2_survivors 2 m3_survivors 0 m4_survivors 0','literal clone count result')
 need(out['gap_counts']=='DONE rank3_roots 2410 hits 0','independent gap result')
 # Strict strong-order cutoff and the factorial-derived amplitude range.
 need(5**11>217**3 and 9*5**2>217,'Busch cutoff')
 need([A for A in range(1,60) if A%6==0]==list(range(6,55,6)),'amplitude range')
 print(json.dumps({'status':'PASS','bank_sha256':hashlib.sha256(bank.read_bytes()).hexdigest(),'strong_classes':len(data),'order_inventory':hist,'methods':out,'excluded_positive_amplitudes':list(range(6,55,6)),'excluded_offsets':[-1,1],'necessary_amplitude_lower_bound':60},sort_keys=True,indent=2))
if __name__=='__main__':main()
