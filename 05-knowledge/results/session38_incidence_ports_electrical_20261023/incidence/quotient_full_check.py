"""Reproduce the complete downward E33 quotient-forest search, plus Python checks.
Compiles the bundled standard C++17 helper locally. No external packages or solver.
"""
from pathlib import Path
from collections import Counter
import argparse,csv,hashlib,io,json,subprocess,tempfile
import converse_check as C
HERE=Path(__file__).resolve().parent

def canon(part):
 d={};return ''.join(str(d.setdefault(x,len(d))) for x in part)
def children(part):
 for a in sorted(set(part)):
  for b in sorted(set(part)):
   if a<b:yield canon(part.replace(b,a))
def recolor(rows,part):return b''.join(sorted(bytes([len(row)]+sorted(int(part[x]) for x in row)) for row in rows))
def need(ok,msg):
 if not ok:raise AssertionError(msg)
def main():
 p=argparse.ArgumentParser();p.add_argument('--classes',type=Path,default=HERE.parent/'wiring/classes_input_n8.tsv');p.add_argument('--read-existing',action='store_true',help='validate existing search file without recompiling or rerunning enumeration');p.add_argument('--write',action='store_true');args=p.parse_args()
 if not args.classes.exists() and args.classes==HERE.parent/'wiring/classes_input_n8.tsv':args.classes=HERE.parent/'session38-wiring/classes_input_n8.tsv'
 R=sorted(csv.DictReader(args.classes.open(),delimiter='\t'),key=lambda r:int(r['canon']));need(len(R)==6880,'class count');look={int(r['class']):i for i,r in enumerate(R)};conv=[look[int(r['sigmaClass'])] for r in R];cc=[int(r['canon']) for r in R];base=[];lines=['6880']
 for i,r in enumerate(R):
  A,c,e=C.incidence(cc[i]);rows=[[c[v]%3*3+e[u][v]%3 for v in range(8) if A[u][v]] for u in range(8)];base.append(rows);s=C.sig(A,c,e,3,3);lines.append(str(conv[i])+' '+''.join(str(x) for x in s))
 data=('\n'.join(lines)+'\n').encode();target=HERE/'quotient_search.tsv'
 if not args.read_existing:
  with tempfile.TemporaryDirectory(prefix='session38-quotient-') as tmp:
   tmp=Path(tmp);(tmp/'input.txt').write_bytes(data)
   subprocess.run(['g++','-O2','-std=c++17',str(HERE/'quotient_forest.cpp'),'-o',str(tmp/'search')],check=True)
   r=subprocess.run([str(tmp/'search'),str(tmp/'input.txt')],check=True,capture_output=True)
   if args.write:target.write_bytes(r.stdout)
   else:need(target.read_bytes()==r.stdout,'search transcript differs')
 records=list(csv.DictReader(io.StringIO(target.read_text()),delimiter='\t'));by={}
 for row in records:
  k=int(row['level']);by.setdefault(k,[]).append(row)
 need(sorted(by,reverse=True)==list(range(8,min(by)-1,-1)),'missing level');parents={'012345678'};summary=[];last_success=None;S=[[0]*10 for _ in range(10)];S[0][0]=1
 for n in range(1,10):
  for k in range(1,n+1):S[n][k]=S[n-1][k-1]+k*S[n-1][k]
 for k in sorted(by,reverse=True):
  expected={q for p in parents for q in children(p)};actual={r['partition'] for r in by[k]};need(len(actual)==len(by[k]),'duplicate partition');need(actual==expected,'incomplete child enumeration');good={r['partition'] for r in by[k] if r['forest']=='1'}
  need(all(len(set(r['partition']))==k for r in by[k]),'wrong block count');summary.append({'symbols':k,'total_partitions':S[9][k],'tested':len(actual),'pruned_by_failed_ancestor':S[9][k]-len(actual),'forest_successes':len(good)})
  if good:last_success=k
  else:need(k==min(by),'continued after zero-success layer')
  parents=good
 need(not parents,'search stopped before first empty success layer');need(last_success==min(by)+1,'minimum boundary')
 successes=[r for r in by[last_success] if r['forest']=='1'];selected='012120201';need(any(r['partition']==selected for r in successes),'ternary sum rule is not a minimum-level success');colors=[recolor(rows,selected) for rows in base];info,G,pc=C.forest_info(colors,conv,cc);need(info['forest'],'selected C++ success fails Python forest check');need(info['paired_signatures']==6880,'forest must imply distinct edge pairs')
 # Independently replay all 36 immediate quotients in Python; this cross-checks the compiled carrier.
 py8=[]
 for row in by[8]:
  part=row['partition'];inf,_,_=C.forest_info([recolor(rows,part) for rows in base],conv,cc);need(inf['forest']==(row['forest']=='1'),'Python/C++ eight-symbol mismatch');py8.append(inf['forest'])
 O={'scope':'FINITE-EXACT minimum forest-certified alphabet within all set partitions of the fixed nine E33 symbols; not absolute wiring impossibility','classes_sha256':hashlib.sha256(args.classes.read_bytes()).hexdigest(),'prepared_input_sha256':hashlib.sha256(data).hexdigest(),'helper_sha256':hashlib.sha256((HERE/'quotient_forest.cpp').read_bytes()).hexdigest(),'search_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'universe_partitions':sum(S[9]),'levels':summary,'minimum_forest_certified_symbols':last_success,'selected_partition':selected,'selected_symbol_blocks':[[[x//3,x%3] for x in range(9) if selected[x]==b] for b in sorted(set(selected))],'selected_python_verification':info,'independent_eight_symbol_replays':len(py8),'all_lower_levels_excluded_by_cycle_monotonicity':True}
 O=json.loads(json.dumps(O));out=HERE/'quotient_full_receipt.json'
 if args.write:out.write_text(json.dumps(O,indent=2)+'\n')
 else:need(json.loads(out.read_text())==O,'frozen full quotient receipt differs')
 print(json.dumps(O,indent=2))
if __name__=='__main__':main()
