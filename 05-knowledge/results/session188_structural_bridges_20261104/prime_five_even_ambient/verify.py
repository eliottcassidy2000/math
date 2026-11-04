from pathlib import Path
import hashlib,json,subprocess,tempfile,sys
ROOT=Path(__file__).resolve().parent
def need(c,m):
 if not c:raise RuntimeError(m)
def main():
 inp='1\n24 17 3 6 12\n2\n0 1\n2 18\n'
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'covers';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'covers.cpp'),'-o',str(exe)],check=True);ls=[list(map(int,s.split())) for s in subprocess.check_output([str(exe)],input=inp.encode()).decode().splitlines()]
 need(len(ls)==4 and all(len(v)==48 for v in ls),'all first-three profile channels');records=[]
 for j,query in enumerate([[0,1],[2,18]]):
  T,U=ls[2*j:2*j+2];base=[(T[k]-U[k])%16 for k in range(3)];derivs=[[(T[3+3*i+k]-U[3+3*i+k])%2 for k in range(3)] for i in range(15)];need(base==[8*j,0,0],'two baseline full cover gaps');expect=[[int(i in[6,9]),0,0] for i in range(15)];need(derivs==expect,'all fifteen perturbation differences, including q3/q5');records.append({'query':query,'source_channels':T,'target_channels':U,'base_contrast':base,'perturbation_contrasts_mod2':derivs})
 histogram={};allstates=[]
 for bits in range(1<<15):
  vv=[[ (r['base_contrast'][k]+8*sum(((bits>>i)&1)*r['perturbation_contrasts_mod2'][i][k] for i in range(15)))%16 for k in range(3)] for r in records];need(sorted(vv)==[[0,0,0],[8,0,0]],'every formal state has full first-three separator');key=tuple(tuple(v) for v in vv);histogram[key]=histogram.get(key,0)+1;allstates.extend(v for row in vv for v in row)
 need(set(histogram.values())=={16384},'complete balanced state count')
 physical=json.loads(subprocess.check_output([sys.executable,str(ROOT/'physical_controls.py')]))
 need(physical['contrast_mod16']==8 and physical['physical_order']==24,'actual strong-ambient separator')
 # For exterior C3 with two singleton children, exact module slots have b1=3,b2=1.
 expected=[(3*ls[2+side][0]+ls[2+side][1])%16 for side in [0,1]];need([r['value_mod16'] for r in physical['records']]==expected,'literal lift matches independently generated core profiles')
 pins={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','covers.cpp','physical.cpp','physical_controls.py']};deps={f:hashlib.sha256((ROOT.parent/f).read_bytes()).hexdigest() for f in ['prime_five_even_profiles/PROOF.md','even_cyclic_mod16/PROOF.md','mod8_transfer/PROOF.md','module_decomposition/PROOF.md','monotone_dual_queries/PROOF.md']}
 print(json.dumps({'status':'PASS','scope':'Prime-five arbitrary-even>=4 child quotient converse in every fixed ambient; special full-cover mod16 certificate, not generic mod16 transfer.','channel_order':['base_q1','base_q2','base_q3']+[f'child{i}_run{k}_q{h}' for i in range(5) for k in [2,3,5] for h in [1,2,3]],'ordered_cover_records':records,'formal_states':32768,'state_histogram':[{'profile_pair':k,'count':v} for k,v in histogram.items()],'complete_state_vector_sha256':hashlib.sha256(bytes(allstates)).hexdigest(),'physical_controls':physical,'source_pins':pins,'dependency_pins':deps},sort_keys=True,indent=2))
if __name__=='__main__':main()
