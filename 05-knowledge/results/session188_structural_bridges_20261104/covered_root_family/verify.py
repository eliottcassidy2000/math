"""Exact root-containing family controls; independent literal predecessor query counts."""
from pathlib import Path
import importlib.util,itertools,json,hashlib
ROOT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('covered',ROOT.parent/'covered_response'/'verify.py');v=importlib.util.module_from_spec(s);s.loader.exec_module(v)
def modules(a):
 n=len(a);return [m for m in range(1,1<<n) if m.bit_count()>=2 and all(a[u]&m in (0,m) for u in range(n) if not m>>u&1)]
def xor(a,fs):
 r=list(a)
 for m in fs:
  for u,w in itertools.combinations(range(len(a)),2):
   if m>>u&1 and m>>w&1:r[u]^=1<<w;r[w]^=1<<u
 return r
T=v.blowup(v.tt(3));M=modules(T);root=(1<<9)-1;proper=[m for m in M if m!=root];v.need(len(proper)==9,'all nine proper interval modules');H=v.hp(T);baseline=v.hp(v.extend(T,7,56));records=[];strict=equal=0
for z in range(1<<len(proper)):
 fs=[root]+[m for i,m in enumerate(proper) if z>>i&1];U=xor(T,fs);h=v.hp(U);v.need(h>=H,'decomposition scalar monotonicity control')
 if h>H:query=[0,0];values=[H,h];strict+=1
 else:
  f=v.hp(v.extend(U,7,56));v.need(f-baseline==-v.gap(v.profile(v.tt(3))),'equal-profile root physical separator');query=[7,56];values=[baseline,f];equal+=1
 v.need(values[0]!=values[1],'actual full response differs');records.append({'family':fs,'query':query,'values':values,'target_rows':U})
deep=[];T=v.blowup(v.tt(4));root=(1<<12)-1;child=[15<<(4*i) for i in range(3)];pairs=[3<<(2*i) for i in range(6)];baseline=v.hp(v.extend(T,15,240));H=v.hp(T)
for fs in ([root]+child+pairs,[root]+pairs,[root]+child+[7,7<<4,7<<8]+[3,3<<4,3<<8]):
 v.need(all(m in modules(T) for m in fs),'deep original modules');U=xor(T,fs);f=v.hp(v.extend(U,15,240));v.need(v.hp(U)==H and f-baseline==-v.gap(v.profile(v.tt(4))),'deep laminar literal separator');deep.append({'family':fs,'query':[15,240],'values':[baseline,f],'target_rows':U})
# Unequal-size corollary: all original-module families at TT2,TT3,TT3 root.
sizes=[2,3,3];off=[0,2,5];T=[0]*8
for i,m in enumerate(sizes):
 for u,row in enumerate(v.tt(m)):
  T[off[i]+u]=(row<<off[i])|(((1<<sizes[(i+1)%3])-1)<<off[(i+1)%3])
root=255;proper=[m for m in modules(T) if m!=root];v.need(len(proper)==7,'unequal proper module count');H=v.hp(T);battery=[(1<<u,0) for u in range(2,5)]+[(0,1<<u) for u in range(2,5)];baseline={q:v.hp(v.extend(T,*q)) for q in battery};mixed=[];mixed_strict=mixed_equal=0
for z in range(1<<len(proper)):
 fs=[root]+[m for i,m in enumerate(proper) if z>>i&1];U=xor(T,fs);h=v.hp(U);v.need(h>=H,'unequal scalar monotonicity control')
 if h>H:mixed_strict+=1;mixed.append({'family':fs,'query':[0,0],'values':[H,h],'target_rows':U});continue
 mixed_equal+=1;values=[v.hp(v.extend(U,*q)) for q in battery];contrast=sum(baseline[battery[j]]-values[j] for j in range(3))-2*sum(baseline[battery[j]]-values[j] for j in range(3,6));v.need(contrast==192,'unequal profile-only start battery');j=next(j for j,q in enumerate(battery) if baseline[q]!=values[j]);mixed.append({'family':fs,'query':battery[j],'values':[baseline[battery[j]],values[j]],'target_rows':U,'endpoint_battery':[[*q,baseline[q],values[j]] for j,q in enumerate(battery)],'battery_contrast':contrast})

print(json.dumps({'status':'PASS','scope':'All root-containing original-module XOR families at an equal-profile C3 whole root, and at a cyclic root with transitive children of sizes at least two; no arbitrary strong ambient lift.','nine_name_source':v.blowup(v.tt(3)),'root_containing_families':len(records),'strict_scalar_families':strict,'scalar_equal_families':equal,'family_records':records,'deep_twelve_name_controls':deep,'unequal_transitive_source':T,'unequal_root_families':len(mixed),'unequal_scalar_strict':mixed_strict,'unequal_scalar_equal':mixed_equal,'unequal_records':mixed,'source_sha256':{str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'PROOF.md',ROOT/'verify.py',ROOT.parent/'covered_response'/'PROOF.md',ROOT.parent/'module_decomposition'/'PROOF.md',ROOT.parent/'covered_quotient'/'PROOF.md']}},sort_keys=True,indent=2))
