"""Exact physical metadata, affine battery and known-obstruction transports."""
from physical_audit import D,P,model,op,need
import hashlib,json

def main():
 R=json.loads((P/'physical_records.json').read_text());A=json.loads((P/'anchor_exploration.json').read_text());checks=0
 for kind in ['anchor_probes','ambient_probes']:
  for p in R[kind]:
   for side in range(2):
    Q=R['records'][p['ri']]['rows'];Q=Q if side==0 else op(Q);T,a,c,B=model(Q,p['V'],*p['query'],kind=='ambient_probes')
    need(T==p['augmented_rows'][side],'retained complete augmented physical rows');need(B==p['child_sets'],'retained actual child names');need(p['port_names']==[len(T)-2,len(T)-1],'retained port names');need(p['physical_queries'][side]==[a,c],'actual full masks')
    if kind=='ambient_probes':need(p['exterior_names']==[len(T)-4,len(T)-3],'retained exterior names')
    checks+=1
 for ri,entry in enumerate(A):
  for row in entry['patterns']:
   V=row['V'];v=lambda j:(V>>j)&1
   if ri==0:
    actual=[row['first'][5]//2,row['first'][7]//2,row['second'][0]//4];expected=[v(2)^v(4),1^v(3)^v(5)^v(6),v(2)^v(3)^v(4)^v(5)^v(6)]
   else:
    actual=[row['first'][0]//2,row['first'][2]//2,row['second'][7]//4];expected=[1^v(1)^v(2)^v(4),v(3)^v(5),v(1)^v(2)^v(3)^v(4)^v(5)]
   need(actual==expected and (actual[0]^actual[1]^actual[2])==1,'three affine anchor bits sum to one')
 known=[238,4,240,118,195,146,34,74];transports=[(known,A[1]['rows'],[1,6,5,7,2,4,3,0]),(op(known),A[0]['rows'],[0,3,4,2,7,5,6,1]),(op(A[0]['rows']),A[1]['rows'],list(reversed(range(8))))]
 for T,U,p in transports:need(sorted(p)==list(range(8)) and all(((T[p[i]]>>p[j])&1)==((U[i]>>j)&1) for i in range(8) for j in range(8)),'known obstruction common physical relabeling')
 return {'status':'PASS','physical_side_metadata_checks':checks,'three_affine_battery_patterns':512,'explicit_shared_transports':3,'proof_sha256':hashlib.sha256((P/'PROOF.md').read_bytes()).hexdigest(),'source_sha256':hashlib.sha256((D/'label_audit.py').read_bytes()).hexdigest(),'physical_input_sha256':hashlib.sha256((P/'physical_records.json').read_bytes()).hexdigest(),'anchor_input_sha256':hashlib.sha256((P/'anchor_exploration.json').read_bytes()).hexdigest()}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
