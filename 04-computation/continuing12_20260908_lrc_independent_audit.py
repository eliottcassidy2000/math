"""Independent unpruned full-profile consumer and native depth-wedge audit.

No producer import or execution. Both full native banks are scanned separately
without incumbent shortcuts; the word engine visits every seven-multiset.
"""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
from math import gcd
import argparse,json,os,shutil,subprocess,sys,tempfile
sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve().parent
STEM=Path(__file__).stem
GATES=0
def need(ok,why):
 global GATES
 GATES+=1
 if not ok:raise ArithmeticError(why)
def canonical(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def main():
 ap=argparse.ArgumentParser();filed=HERE.name=='04-computation'
 ap.add_argument('--root',type=Path,default=HERE.parent if filed else Path('C:/w/s0905'))
 ap.add_argument('--producer',type=Path,default=HERE.parent/'05-knowledge/results' if filed else Path('C:/w/continuing12_20260908_lrc'))
 ap.add_argument('--work-dir',type=Path,default=Path(tempfile.gettempdir())/STEM)
 args=ap.parse_args()
 raw=(args.root/'04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json').read_bytes()
 need(sha256(raw).hexdigest()=='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f','full inherited profile pin')
 P=json.loads(raw)['levels'];states=sorted({row[0] for row in P['6']['profiles']})
 need(states==P['6']['gcds'] and len(states)==42,'whole inherited alphabet')
 oldraw=(args.root/'05-knowledge/results/continuing11_20260908_lrc_composite_clocks_certificate.json').read_bytes()
 need(sha256(oldraw).hexdigest()=='2c9b4fa1a1103d3309094d06a422ec2f5e56e9b804961037cab9d5a05a352b73','actual predecessor certificate')
 old=json.loads(oldraw)['new_scales']
 need(len(old)==7622 and max(old)==11935 and sha256(canonical(old)).hexdigest()=='832714c179ab4f425a76172c41d7d8eecf05dc15adfdd3c5e08bc7eb5a13e48b','whole actual predecessor array')
 pr=(args.producer/'continuing12_20260908_lrc_zero_clock_probe_certificate.json').read_bytes()
 wr=(args.producer/'continuing12_20260908_lrc_depth_wedge_certificate.json').read_bytes()
 need(sha256(pr).hexdigest()=='e8caa081f3d9e8bffb2c7e4a003256ae4b39e1116799fd97294ceac014febf2b','frozen complete probe certificate')
 need(sha256(wr).hexdigest()=='c56ea56a87c98cae1a16efe186696d4424942a9325d4c4189c5ce2bedbe3abea','frozen whole depth-wedge certificate')
 J=json.loads(pr);Z=json.loads(wr);records=J['clocks']
 need(all(t in old for t in [6930,7560,9240,11088]),'entire ordered declared universe in predecessor')
 need([R['t'] for R in records]==[6930,7560] and J['removed_scales']==[6930],'exact stop after first complete residual')
 need(not records[0]['survivors'] and len(records[1]['survivors'])==6,'first residual is complete and follows one closure')
 need(J['profile_sha256']==sha256(raw).hexdigest() and J['baseline_semantic_sha256']==sha256(canonical(old)).hexdigest(),'supplier identities retained')
 middle=[t for t in old if t!=6930]
 need(J['new_scales']==middle and len(middle)==7621 and sha256(canonical(middle)).hexdigest()==J['new_scales_sha256']=='567d76d5b776e805a0606298399e7ecedc14657c6910f03321ca7e64e6c4c092','minimum-tree intermediate array')
 final=[t for t in middle if t!=7560]
 need(Z['new_scales']==final and Z['removed_scales']==[7560] and len(final)==7620 and max(final)==11935,'depth theorem removes only 7560')
 need(sha256(canonical(final)).hexdigest()==Z['new_scales_sha256']=='66f903430fd804074e6c8adb330aca4f6bb6b8f97d955e014a73309c39d17a8d','entire final semantic array')
 need(Z['probe_sha256']==sha256(pr).hexdigest() and Z['profile_word_count']==37209,'actual full residual supplier')
 need(all(gcd(d,7)==1 for d in states),'all inherited sheet states coprime to seven')
 domains=sorted({tuple(d for d in states if R['t']%d==0) for R in records})
 bankpins={tuple(D['domain']):D['words_sha256'] for D in J['domains']}
 need(set(domains)==set(bankpins) and sorted(map(len,domains))==[16,23],'complete exact divisor domains')
 profiles=[(int(k),c,w) for k,row in P.items() for c,w in row['profiles']]
 lines=[str(len(profiles))]+[' '.join(map(str,[k,c]+w)) for k,c,w in profiles]
 lines+=[str(len(domains))]+[' '.join(map(str,[len(D)]+list(D))) for D in domains]+[str(len(records))]
 for R in records:
  D=tuple(d for d in states if R['t']%d==0)
  need(list(D)==R['domain'],'clock divisibility determines every domain')
  need(all(R['t']%(7*d)==0 for d in D),'zero budget on the whole native domain')
  need(sha256(canonical(R['words'])).hexdigest()==bankpins[D] and len(R['words'])==R['word_count'],'recorded whole word array hash')
  table={tuple(ab):v[0] for ab,v in R['weights']}
  need(set(table)==set(combinations(D,2))|{(d,d) for d in D},'complete pair table including repeated sheets')
  parent=list(range(7))
  def find(i):
   while parent[i]!=i:i=parent[i]
   return i
  total=0;need(len(R['owner_edges'])==6,'complete owner edge list')
  for i,j,cost in R['owner_edges']:
   need(0<=i<j<7 and find(i)!=find(j),'owner positional tree')
   parent[find(i)]=find(j)
   need(cost==table[tuple(sorted((R['owner'][i],R['owner'][j])))],'owner uses actual pair costs')
   total+=cost
  need(total==R['owner_tree'],'full owner cost')
  lines+=[' '.join(map(str,[R['t'],domains.index(D),R['word_count'],R['minimum_margin'],R['owner_E'],R['owner_tree'],R['event_evaluations'],R['compatible_atlas_edges']]+R['owner']))]
  lines+=[str(len(R['survivors']))]+[' '.join(map(str,w+[e,m])) for w,e,m in R['survivors']]
  lines+=[str(len(R['weights']))]+[' '.join(map(str,ab+[len(v)]+v)) for ab,v in R['weights']]
 # Derive the entire positional zero graph; enumerate every connected subgraph.
 R=records[1];table={tuple(ab):v[0] for ab,v in R['weights']}
 need([T['word'] for T in Z['topology']]==[w for w,e,m in R['survivors']],'complete topology corresponds to all six residuals')
 connected_controls=0
 for T,(word,e,m) in zip(Z['topology'],R['survivors']):
  need(e==m==0 and word.count(24)==word.count(18)==1 and 4 in word,'entire residual type and zero cost')
  edges=[(i,j) for i,j in combinations(range(7),2) if table[tuple(sorted((word[i],word[j]))) ]==0]
  need([list(e) for e in edges]==T['zero_edges'],'every positional zero edge, not a selected tree')
  mid=word.index(24);other=word.index(18);four=[i for i,d in enumerate(word) if d==4]
  need(all({j if i==v else i for i,j in edges if v in (i,j)}=={mid} for v in four),'each actual four endpoint forces unique middle')
  part={i for i,d in enumerate(word) if d in (9,18,27)}
  cross=[e for e in edges if (e[0] in part)!=(e[1] in part)]
  need(cross==[tuple(sorted((mid,other)))],'entire partition forces actual 18-24 bridge')
  for mask in range(1<<len(edges)):
   picked=[edge for j,edge in enumerate(edges) if mask>>j&1]
   seen={0}
   while True:
    nxt=seen|{j for i,j in picked if i in seen}|{i for i,j in picked if j in seen}
    if nxt==seen:break
    seen=nxt
   if len(seen)!=7:continue
   connected_controls+=1
   need(tuple(sorted((mid,other))) in picked and all(tuple(sorted((v,mid))) in picked for v in four),'every connected actual subgraph contains the depth wedge')
 # Supplied separate realizations and phases are checked on literal finite grids.
 depthsets=[]
 def vp(n):
  v=0
  while n%2==0:v+=1;n//=2
  return v
 for bank in Z['zero_banks']:
  depths=set()
  need(bank['e']==gcd(bank['end'],24) and bank['n']==7560//bank['e'],'actual changed quotient clock')
  for row in bank['zero']:
   p,q=row['primitive_pair'];num,den=row['ratio'];alpha,beta=row['zero_phase'];n=bank['n']
   need(gcd(p,q)==1 and sorted([p,q])==sorted([num,den]),'literal primitive directed ratio')
   need(row['separate_realization']==[bank['e']*num,bank['e']*den],'separate actual arm realization')
   need([gcd(7560,s) for s in row['separate_realization']]==[bank['end'],24],'separate arm has requested margins')
   count=0
   for j in range(n):
    xx=j*beta+alpha;D=n*beta
    hit=lambda v:14*min((v*xx)%D,D-(v*xx)%D)<D
    count+=hit(p) and hit(q)
   need(count==0,'literal complete quotient grid at the recorded separate zero phase')
   need(vp(bank['end'])<vp(7560),'unclipped endpoint depth')
   depth=vp(bank['end'])-vp(num)+vp(den);depths.add(depth)
   need(depth==row['middle_depth'],'same physical middle depth determined by native ratio')
  depthsets.append(sorted(depths))
 need(depthsets==Z['middle_depth_sets']==[[6],[4,5,8]] and not set(depthsets[0])&set(depthsets[1]),'entire two depth sets are incompatible')
 # Saturation hostile: local margin compatibility is insufficient without endpoint unsaturation.
 need([gcd(2,s) for s in [1,4]]==[1,2] and [gcd(2,s) for s in [2,3]]==[2,1],'both hostile arms separately satisfy their proposed local margins')
 need(1*4==4*1 and 4*3==6*2 and gcd(gcd(1,4),6)==1,'hostile joint primitive realization retains both native ratios')
 need([gcd(2,s) for s in [1,4,6]]==[1,2,2],'saturation hostile retains its actual changed third margin')
 work=args.work_dir/('optimized' if sys.flags.optimize else 'normal');work.mkdir(parents=True,exist_ok=True)
 inp=work/'input.txt';inp.write_bytes(('\n'.join(lines)+'\n').encode())
 compiler=shutil.which('g++') or 'C:/Users/Eliott/scoop/apps/gcc/current/bin/g++.exe'
 exe=work/('referee.exe' if os.name=='nt' else 'referee')
 options=['-std=c++17','-O3','-DNDEBUG'] if sys.flags.optimize else ['-std=c++17','-O2']
 subprocess.run([compiler,*options,str(HERE/(STEM+'.cpp')),'-o',str(exe)],capture_output=True,check=True)
 env=os.environ.copy();env['PATH']=str(Path(compiler).resolve().parent)+os.pathsep+env.get('PATH','')
 with (work/'progress.txt').open('wb') as err:proc=subprocess.run([str(exe),str(inp),str(work)],stdout=subprocess.PIPE,stderr=err,env=env)
 if proc.returncode:raise ArithmeticError((work/'progress.txt').read_text()[-4000:])
 for j,D in enumerate(domains):need(sha256((work/('words_'+str(j)+'.json')).read_bytes()).hexdigest()==bankpins[D],'complete independently unpruned word-list digest')
 print('INDEPENDENT: full unpruned126-profile domains; raw spatial cells; signed all-phase sweep; Kruskal; complete native banks.')
 print(proc.stdout.replace(b'\r\n',b'\n').decode().strip())
 print('CONNECTED_ZERO_SUBGRAPHS',connected_controls,'ALL_FORCE_THE_SAME_MIDDLE')
 print('FINAL_NECESSARY_ARRAY',len(final),'MAXIMUM',max(final),'SEMANTIC_SHA256',sha256(canonical(final)).hexdigest())
 print('PASS',GATES,'always-active Python completeness/topology/native-depth gates; output LF')
if __name__=='__main__':main()
