"""Independent unpruned full-profile consumer and complete native zero-cut audit.

No producer import or execution. All native ratio banks are scanned separately
without incumbent shortcuts; the word engine visits every seven-multiset.
"""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
from collections import Counter
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
 ap.add_argument('--producer',type=Path,default=HERE.parent/'05-knowledge/results' if filed else Path('C:/w/continuing13_20260908_lrc'))
 ap.add_argument('--work-dir',type=Path,default=Path(tempfile.gettempdir())/STEM)
 args=ap.parse_args()
 raw=(args.root/'04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json').read_bytes()
 need(sha256(raw).hexdigest()=='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f','full inherited profile pin')
 P=json.loads(raw)['levels'];states=sorted({row[0] for row in P['6']['profiles']})
 need(states==P['6']['gcds'] and len(states)==42,'whole inherited alphabet')
 oldraw=(args.root/'05-knowledge/results/continuing12_20260908_lrc_depth_wedge_certificate.json').read_bytes()
 need(sha256(oldraw).hexdigest()=='c56ea56a87c98cae1a16efe186696d4424942a9325d4c4189c5ce2bedbe3abea','actual predecessor certificate')
 old=json.loads(oldraw)['new_scales']
 need(len(old)==7620 and max(old)==11935 and sha256(canonical(old)).hexdigest()=='66f903430fd804074e6c8adb330aca4f6bb6b8f97d955e014a73309c39d17a8d','whole predecessor array')
 pr=(args.producer/'continuing13_20260908_lrc_zero_clock_probe_certificate.json').read_bytes()
 wr=(args.producer/'continuing13_20260908_lrc_zero_cuts_certificate.json').read_bytes()
 need(sha256(pr).hexdigest()=='a56f8d8bdee84f642d7f837a6c9e6398008eadbc0318b0778d6d0e435c0b0201','frozen complete probe certificate')
 need(sha256(wr).hexdigest()=='23ae81ba6cb39bc14aecf2d93fc3fcf6d4e3c587b5240db35ea1d8f3992551f1','frozen zero-cut certificate')
 J=json.loads(pr);Z=json.loads(wr);records=J['clocks']
 need(J['declared_clocks']==Z['declared_clocks']==[9240,11088],'exact complete declared universe')
 need([R['t'] for R in records]==[9240,11088] and all(not R['survivors'] for R in records),'both complete clocks have no survivors')
 need(J['profile_sha256']==sha256(raw).hexdigest() and J['baseline_semantic_sha256']==sha256(canonical(old)).hexdigest(),'actual inherited profiles and array')
 final=[t for t in old if t not in [9240,11088]]
 need(J['new_scales']==Z['new_scales']==final and len(final)==7618 and max(final)==11935,'entire actual final array')
 need(J['removed_scales']==Z['removed_scales']==[9240,11088],'remove only the two proved clocks')
 need(sha256(canonical(final)).hexdigest()==J['new_scales_sha256']==Z['new_scales_sha256']=='89a3e7545cc77467c5f85fbe0bdaab1f71071cff65184bbc1212eb86c9f07177','final semantic array pin')
 need(Z['probe_sha256']==sha256(pr).hexdigest(),'sidecar uses actual complete probe')
 need(all(gcd(d,7)==1 for d in states),'all inherited states coprime to seven')
 domains=sorted({tuple(d for d in states if R['t']%d==0) for R in records})
 bankpins={tuple(D['domain']):D['words_sha256'] for D in J['domains']}
 need(set(domains)==set(bankpins) and sorted(map(len,domains))==[20,21],'complete exact divisor domains')
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
 # Compare the complete independent native zero sets, then every positional cut.
 exceptional=0;words_checked=0
 for R,T in zip(records,Z['records']):
  clock=R['t'];need(T['clock']==clock and T['domain']==R['domain'],'exact sidecar clock and domain')
  actual={tuple(map(int,line.split())) for line in (work/('zeros_'+str(clock)+'.txt')).read_text().splitlines()}
  expected=set()
  for pair,rows in T['zero_banks']:
   for p,q,dp,dq,num,den in rows:
    e=gcd(dp,dq)
    need(sorted([dp,dq])==pair,'declared zero ratio retains actual margin roles')
    expected.add((e,p,q,dp,dq))
  need(actual==expected and len(actual)==T['zero_ratio_count'],'entire independently unpruned native zero bank')
  zero_pairs={tuple(sorted((row[3],row[4]))) for row in actual}
  table={tuple(ab):v[0] for ab,v in R['weights']}
  need(zero_pairs=={ab for ab,v in table.items() if v==0},'entire zero relation agrees with pair minima')
  def components(word):
   remaining=set(range(len(word)));parts=[]
   while remaining:
    comp={min(remaining)}
    while True:
     nxt=comp|{j for j in remaining if any(tuple(sorted((word[i],word[j]))) in zero_pairs for i in comp if i!=j)}
     if nxt==comp:break
     comp=nxt
    parts.append(sorted(comp));remaining-=comp
   return parts
  shapes=Counter();isolated=0;outliers=[]
  for word in R['words']:
   parts=components(word);shape=tuple(sorted(map(len,parts)));shapes[shape]+=1;words_checked+=1
   need(len(parts)>1,'every full-profile positional zero graph disconnects')
   if 1 in shape:isolated+=1
   else:
    exceptional+=1;outliers.append(word)
    report=next(B for B in T['no_isolation'] if B['word']==word)
    need(sorted(parts)==sorted(report['components']),'every non-isolation component retained')
    need(sorted(sorted(word[i] for i in part) for part in parts)==sorted(report['component_margins']),'component margins retain repeated positions')
    cut=set(report['cut']);need(sorted(cut) in parts,'declared cut is a full actual zero component')
    crossing=[table[tuple(sorted((word[i],word[j])))] for i in cut for j in range(7) if j not in cut]
    need(min(crossing)==report['minimum_crossing_credit']>0,'complete positive crossing minimum on each exceptional cut')
  need(isolated==T['isolated_word_count'] and outliers==[B['word'] for B in T['no_isolation']],'entire isolated versus exceptional partition')
  need([[list(k),v] for k,v in sorted(shapes.items())]==T['shapes'],'complete component-size histogram')
  alphabet=components(R['domain'])
  need(sorted(sorted(R['domain'][i] for i in p) for p in alphabet)==sorted(T['alphabet_components']),'full alphabet sidecar retains nontrivial connected components')
 need(exceptional==4 and words_checked==41249,'entire exact positional cut universe')
 print('INDEPENDENT: all126 profiles on every unpruned multiset; raw spatial geometry; signed phase sweep; Kruskal; every native ratio and positional cut.')
 print(proc.stdout.replace(b'\r\n',b'\n').decode().strip())
 print('POSITIONAL_CUTS',words_checked,'DISCONNECTED; NON_ISOLATION_EXCEPTIONS',exceptional)
 print('FINAL_NECESSARY_ARRAY',len(final),'MAXIMUM',max(final),'SEMANTIC_SHA256',sha256(canonical(final)).hexdigest())
 print('PASS',GATES,'always-active Python completeness/native-bank/whole-cut gates; output LF')
if __name__=='__main__':main()
