"""Single fixed-query original-family controls, with full physical counts."""
from pathlib import Path
import hashlib,importlib.util,itertools,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('v',ROOT.parent/'covered_response'/'verify.py');v=importlib.util.module_from_spec(s);s.loader.exec_module(v)
def modules(rows):
 n=len(rows);out=[]
 for m in range(1,(1<<n)-1):
  if m.bit_count()<2:continue
  if all(not(m>>u&1) and (rows[u]&m) in (0,m) or m>>u&1 for u in range(n)):out.append(m)
 return out
def flip(rows,family):
 a=rows[:]
 for m in family:
  for u in range(len(a)):
   if m>>u&1:a[u]^=m^(1<<u)
 return a
def ambient(core,cyclic):
 if not cyclic:return core[:],7,56,1
 n=len(core);rows=[r|(1<<n) for r in core]+[1<<(n+1),(1<<n)-1];return rows,7|(1<<(n+1)),56|(1<<(n+1)),3
def main():
 tasks=[];records=[];sources=[]
 for label,child in [('TT3',v.tt(3)),('C3',[2,4,1])]:
  T=v.blowup(child);mods=modules(T);v.need(len(mods)==(9 if label=='TT3' else 3),'complete proper module list');q=v.profile(child);P=-v.gap(q);H=v.hp(T);profile_cache={}
  for cyclic in [False,True]:
   old,a,c,outerH=ambient(T,cyclic);source_index=len(tasks);tasks.append((old,a,c));source={'child':label,'cyclic_exterior':cyclic,'old_rows':old,'query':[a,c],'count_index':source_index,'P':P,'outerH':outerH};sources.append(source)
   for bits in range(1<<len(mods)):
    family=[511]+[m for i,m in enumerate(mods) if bits>>i&1];U=flip(T,family);target_profiles=[]
    for block in range(3):
     local=tuple((U[3*block+i]>>(3*block))&7 for i in range(3))
     if local not in profile_cache:profile_cache[local]=v.profile(list(local))
     target_profiles.append(profile_cache[local]);v.need(all(x>=y for x,y in zip(profile_cache[local],q)),'hereditary child domination')
    old,a,c,outerH=ambient(U,cyclic);index=len(tasks);tasks.append((old,a,c));records.append({'child':label,'cyclic_exterior':cyclic,'family':family,'source_index':source_index,'target_index':index,'target_old_rows':old,'query':[a,c],'target_child_profiles':target_profiles,'scalar_strict':v.hp(U)>H,'minimum_gap':outerH*P})
 payload=str(len(tasks))+'\n'+''.join(f'{len(rows)} {a} {c} '+' '.join(map(str,rows))+'\n' for rows,a,c in tasks)
 with tempfile.TemporaryDirectory() as td:
  binary=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(binary)],check=True);values=list(map(int,subprocess.check_output([str(binary)],input=payload.encode()).split()))
 v.need(len(values)==1044 and len(records)==1040,'full physical workload')
 for r in records:
  r['source_value']=values[r.pop('source_index')];r['target_value']=values[r.pop('target_index')];r['gap']=r['target_value']-r['source_value'];v.need(r['gap']>=r['minimum_gap']>0,'one fixed query always separates')
 for r in sources:r['value']=values[r.pop('count_index')]
 v.need(sum(r['scalar_strict'] for r in records if r['child']=='TT3' and not r['cyclic_exterior'])==296,'crossing scalar-strict controls')
 deps=[ROOT.parent/'module_decomposition'/'PROOF.md',ROOT.parent/'covered_response'/'PROOF.md',ROOT.parent/'covered_response'/'verify.py',ROOT.parent/'equal_small_ambient'/'PROOF.md',ROOT.parent/'equal_nine_ambient'/'PROOF.md',ROOT.parent/'equal_ten_ambient'/'PROOF.md'];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
 return {'status':'PASS','scope':'One fixed query for all root-containing original-module families at equal-profile cyclic source roots; arbitrary ambient only within separately certified baseline domains.','family_context_records':len(records),'physical_counts':len(values),'sources':sources,'records':records,'source_pins':{p.name:sha(p) for p in [ROOT/'PROOF.md',ROOT/'verify.py',ROOT/'count.cpp']},'dependency_pins':{str(p.relative_to(ROOT.parent)):sha(p) for p in deps}}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
