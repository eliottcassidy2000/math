"""Independent literal-path consumer and source/target capacity replay.
Default compares frozen receipt. --write intentionally regenerates it.
"""
from pathlib import Path
from collections import Counter,defaultdict
from itertools import combinations
import argparse,csv,gzip,hashlib,heapq,json,importlib.util
import numpy as np
P=Path(__file__).resolve().parent;ROOT=P.parents[2];C=6880;N=1<<21;ALL=N-1
checks=0
def need(test,why):
 global checks
 checks+=1
 if not test:raise AssertionError(why)
def data(name):return gzip.decompress((P/(name+'.gz')).read_bytes())
def matrix(code):
 A=[[0]*8 for _ in range(8)]
 for k,(u,v) in enumerate(combinations(range(8),2)):A[u][v]=(code>>(27-k))&1;A[v][u]=1-A[u][v]
 return A
def sigs(A):
 v4=[0]*8
 for B in combinations(range(8),4):
  if sorted(sum(A[u][v] for v in B) for u in B)==[1,1,2,2]:
   for v in B:v4[v]+=1
 records=[(x,) for x in v4];out=[tuple(sorted(records))]
 for r in (1,2):
  records=[(records[u],tuple(sorted(records[v] for v in range(8) if A[u][v]))) for u in range(8)];out.append(tuple(sorted(records)))
 return out

def peel(rowcap,colcap,edges):
 r=dict(rowcap);c=dict(colcap);rn=defaultdict(set);cn=defaultdict(set);q=[];fixed={}
 for a,b in edges:
  if r[a] and c[b]:rn[a].add(b);cn[b].add(a)
 for a in r:
  if len(rn[a])==1:heapq.heappush(q,(0,a))
 for b in c:
  if len(cn[b])==1:heapq.heappush(q,(1,b))
 def remove(a,b):
  rn[a].discard(b);cn[b].discard(a)
  if r[a] and len(rn[a])==1:heapq.heappush(q,(0,a))
  if c[b] and len(cn[b])==1:heapq.heappush(q,(1,b))
 while q:
  side,v=heapq.heappop(q)
  if side==0:
   if not r[v] or len(rn[v])!=1:continue
   a=v;b=next(iter(rn[v]));k=r[a]
  else:
   if not c[v] or len(cn[v])!=1:continue
   b=v;a=next(iter(cn[v]));k=c[b]
  need(k<=r[a] and k<=c[b],'peeling infeasible');r[a]-=k;c[b]-=k;fixed[a,b]=fixed.get((a,b),0)+k
  if not r[a]:
   for bb in list(rn[a]):remove(a,bb)
  if not c[b]:
   for aa in list(cn[b]):remove(aa,b)
 need(all(not k or rn[a] for a,k in r.items()),'source residual unsupported');need(all(not k or cn[b] for b,k in c.items()),'target residual unsupported')
 return fixed,{a:k for a,k in r.items() if k},{b:k for b,k in c.items() if k}

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
 source=ROOT/'04-computation/experiments/tiling_line_metagraph_census_20261020_classes_n8.tsv';rows=list(csv.DictReader(source.open(),delimiter='\t'));need(len(rows)==C,'classes');mats=[matrix(int(r['canon'])) for r in rows];auts=np.array([int(r['Aut']) for r in rows]);Hs=np.array([int(r['H']) for r in rows]);hs=np.array([int(r['tc']) for r in rows]);bs=np.array([int(r['bm']) for r in rows]);conv=np.array([int(r['sigmaClass']) for r in rows]);owner=np.frombuffer(data('owner_map.u16'),dtype='<u2');need(len(owner)==N and int(owner.max())<C,'owner carrier');seen=bytearray(N);ts=list((i,j) for i in range(8) for j in range(i+2,8));ix={e:k for k,e in enumerate(ts)}
 # Independent path construction: append each new vertex and form only its new chord bits.
 total_paths=0
 for ci,A in enumerate(mats):
  out=[sum(b<<v for v,b in enumerate(r)) for r in A];word=[];count=0
  def dfs(used,mask):
   nonlocal count
   k=len(word)
   if k==8:
    need(int(owner[mask])==ci,'literal source certificate owner mismatch');seen[mask]+=1;count+=1;return
   avail=(out[word[-1]] if word else 255)&~used
   while avail:
    bit=avail&-avail;avail-=bit;v=bit.bit_length()-1;new=mask
    for i,u in enumerate(word[:-1]):
     if out[u]>>v&1:new|=1<<ix[i,k]
    word.append(v);dfs(used|bit,new);word.pop()
  dfs(0,0);need(count==int(Hs[ci]),'independent Hamilton count');total_paths+=count
 need(np.array_equal(np.frombuffer(seen,dtype=np.uint8),auts[owner]),'Aut orbit multiplicity');need(np.array_equal(np.bincount(owner,minlength=C),hs),'source fibre counts');ham_checksum=int(Hs[owner].sum());need(ham_checksum==int((Hs*Hs//auts).sum())==815136768,'Session37 Hamilton polynomial checksum');print('independent path coverage',total_paths,flush=True)
 masks=np.arange(N,dtype=np.uint32);sigma=np.zeros(N,dtype=np.uint32)
 for k,(u,v) in enumerate(ts):sigma|=((masks>>k)&1)<<ix[7-v,7-u]
 need(np.array_equal(owner[sigma],conv[owner]),'all literal converse targets');blue=masks==sigma;need(np.array_equal(np.bincount(owner[blue],minlength=C),bs),'blue counts')
 closed=masks[(masks&(1<<ix[0,7]))==0];src=owner[closed].astype(np.uint64);dst=owner[closed^ALL].astype(np.uint64);sa=1+np.bitwise_count(closed&sum(1<<ix[0,j] for j in range(2,8)));sb=6-np.bitwise_count(closed&sum(1<<ix[i,7] for i in range(6)));col=blue[closed].astype(np.uint64)
 keys=(((src*C+dst)*8+sa)*8+sb)*2+col;uk,ct=np.unique(keys,return_counts=True);stored=np.frombuffer(data('ports.bin'),dtype=[('key','<u8'),('count','<u4')]);need(np.array_equal(uk,stored['key']) and np.array_equal(ct,stored['count']),'all source/target ports')
 # Independent target endpoint capacities from literal OPEN presentations.
 opened=closed^ALL;oa=1+np.bitwise_count(opened&sum(1<<ix[0,j] for j in range(2,8)));ob=6-np.bitwise_count(opened&sum(1<<ix[i,7] for i in range(6)));need(np.array_equal(oa,8-sa) and np.array_equal(ob,6-sb),'literal endpoint transport')
 node_rows=list(csv.DictReader((P/'merged_nodes.tsv').open(),delimiter='\t'));need(len(node_rows)==3528,'merged nodes');mid=np.zeros(C,dtype=np.uint16)
 for r in node_rows:
  j,a,b=map(int,(r['merged_id'],r['class_a'],r['class_b']));need(b==int(conv[a]),'node converse pair');mid[a]=mid[b]=j;need(int(r['fibre'])==(int(hs[a]) if a==b else int(hs[a]+hs[b])) and int(r['blue'])==int(bs[a]),'node capacities')
 line=masks[:N//2];aa=owner[line];bb=owner[line^ALL];edge_summaries={}
 for name,uu,vv in [('ordinary_edges.tsv',aa,bb),('merged_edges.tsv',mid[aa],mid[bb])]:
  a=np.minimum(uu,vv).astype(np.uint64);b=np.maximum(uu,vv).astype(np.uint64);ek=((a<<16)|b)*2+blue[line];ekey,ecount=np.unique(ek,return_counts=True);pair=ekey>>1;uniq,starts=np.unique(pair,return_index=True);bc=np.zeros(len(uniq),dtype=np.uint32);kc=bc.copy();which=np.searchsorted(uniq,pair);bc[which[(ekey&1)==1]]=ecount[(ekey&1)==1];kc[which[(ekey&1)==0]]=ecount[(ekey&1)==0]
  text=data(name).decode().splitlines();need(text[0]=='source\ttarget\tblue\tblack','edge header');need(len(text)-1==len(uniq),'edge support size')
  for k,line_text in enumerate(text[1:]):
   a0,b0,x,y=map(int,line_text.split('\t'));need((a0,b0,x,y)==(int(uniq[k]>>16),int(uniq[k]&65535),int(bc[k]),int(kc[k])),'edge mismatch')
  if name.startswith('merged'):need(np.all(kc%2==0),'black pair multiplicities')
  edge_summaries[name]={'support_pairs':len(uniq),'blue_lines':int(bc.sum()),'black_lines':int(kc.sum()),'blue_self_lines':int(bc[(uniq>>16)==(uniq&65535)].sum()),'black_self_lines':int(kc[(uniq>>16)==(uniq&65535)].sum()),'blue_multiplicity_histogram':dict(sorted(Counter(map(int,bc[bc>0])).items())),'black_multiplicity_histogram':dict(sorted(Counter(map(int,kc[kc>0])).items()))}
 print('complete edge tables matched',edge_summaries['merged_edges.tsv']['support_pairs'],flush=True)
 # Decode native port groups once. Every target signature is computed on its certified actual tournament.
 z=uk.copy();colors=(z%2).astype(np.uint8);z//=2;bvals=(z%8).astype(np.uint8);z//=8;avals=(z%8).astype(np.uint8);z//=8;targets=(z%C).astype(np.uint16);sources=(z//C).astype(np.uint16);all_sigs=[sigs(A) for A in mats];levels={}
 api_path=P.parent/'metagraph_minimal_signatures_20261023/check.py';spec=importlib.util.spec_from_file_location('minimal_signatures',api_path);api=importlib.util.module_from_spec(spec);spec.loader.exec_module(api)
 for i,A in enumerate(mats):all_sigs[i].append(api.adaptive_wedge_signature(A))
 need(sum(v[3][1]!=0 for v in all_sigs)==2,'adaptive wedge bit support')
 for level in range(4):
  bank={};class_sig=[];members=defaultdict(list)
  for ci,v in enumerate(all_sigs):
   signature=v[level]
   if signature not in bank:bank[signature]=len(bank)
   sid=bank[signature];class_sig.append(sid);members[sid].append(ci)
  by_port=defaultdict(list)
  for i,cj in enumerate(targets):by_port[class_sig[int(cj)]].append(i)
  stats=Counter();residuals=[];row_total=0
  for sid,indices in by_port.items():
   if len(members[sid])==1:
    stats['rows']+=len(indices);stats['allowed_entries']+=len(indices);stats['forced_entries']+=len(indices);stats['forced_mass']+=int(ct[indices].sum());continue
   rowids={};colids={};rk=[];ck=[];rowcap=Counter();colcap=Counter();actual=Counter()
   for i in indices:
    row=(int(sources[i]),int(avals[i]),int(bvals[i]),int(colors[i]));column=(int(targets[i]),8-int(avals[i]),6-int(bvals[i]),int(colors[i]));mass=int(ct[i])
    if row not in rowids:rowids[row]=len(rk);rk.append(row)
    if column not in colids:colids[column]=len(ck);ck.append(column)
    r=rowids[row];c=colids[column];rowcap[r]+=mass;colcap[c]+=mass;actual[r,c]+=mass
   allowed=[]
   for r,(ci,a,b,color) in enumerate(rk):
    for cj in members[sid]:
     column=(cj,8-a,6-b,color)
     if column in colids:allowed.append((r,colids[column]))
   fixed,rr,cc=peel(rowcap,colcap,allowed);need(all(actual[e]==v for e,v in fixed.items()),'peeling differs from native actual')
   stats.update({'rows':len(rk),'allowed_entries':len(allowed),'forced_entries':len(fixed),'forced_mass':sum(fixed.values()),'residual_mass':sum(rr.values()),'residual_rows':len(rr),'residual_columns':len(cc)})
   if rr:residuals.append({'target_class_ids':members[sid],'target_canonical_words':[int(rows[c]['canon']) for c in members[sid]],'mass':sum(rr.values()),'rows':len(rr),'columns':len(cc)})
  need(stats['forced_mass']+stats['residual_mass']==1<<20,'capacity conservation');levels['adaptive_wedge' if level==3 else str(level)]={'class_signatures':len(bank),**dict(stats),'residual_groups':residuals};print('v4 refinement',level,len(bank),dict(stats),flush=True)
 # A reflected2x2 counterfeit proves the one-round residual is genuinely nonunique.
 ambiguous={3490,3676,4239,5173};native_amb={};examples={}
 for i,cj in enumerate(targets):
  if int(cj) in ambiguous:native_amb[int(sources[i]),int(avals[i]),int(bvals[i]),int(colors[i]),int(cj)]=int(ct[i])
 rows_by_type=defaultdict(list)
 for (ci,a,b,color,cj),mass in native_amb.items():
  if cj in (3490,3676) and color==0:rows_by_type[a,b,color].append((ci,cj))
 counterfeit=None
 for typ,values in sorted(rows_by_type.items()):
  a,b,color=typ
  left=sorted({ci for ci,cj in values if cj==3490});right=sorted({ci for ci,cj in values if cj==3676})
  for ci in left:
   if counterfeit:break
   for di in right:
    if ci==di or int(mid[ci])==int(mid[di]):continue
    delta=Counter({(ci,a,b,color,3490):-1,(ci,a,b,color,3676):1,(di,a,b,color,3676):-1,(di,a,b,color,3490):1})
    for (ei,ea,eb,ec,ej),value in list(delta.items()):delta[int(conv[ei]),7-eb,7-ea,ec,int(conv[ej])]+=value
    if any(native_amb.get(e,0)+v<0 for e,v in delta.items()):continue
    need(all(all_sigs[e[4]][1]==all_sigs[3490 if e[4] in (3490,3676) else 4239][1] for e in delta),'fake target signature unsupported')
    row_delta=Counter();column_delta=Counter();edge_delta=Counter()
    for (ei,ea,eb,ec,ej),value in delta.items():
     row_delta[ei,ea,eb,ec,all_sigs[ej][1]]+=value;column_delta[ej,8-ea,6-eb,ec]+=value;x,y=sorted((int(mid[ei]),int(mid[ej])));edge_delta[x,y]+=value
    need(not any(row_delta.values()) and not any(column_delta.values()),'counterfeit margin failure');need(any(edge_delta.values()) and all(v%2==0 for v in edge_delta.values()),'counterfeit fails merged twin changes')
    witness=[]
    for e,value in sorted(delta.items()):
     if not value:continue
     ei,ea,eb,ec,ej=e
     selected=closed[(src==ei)&(dst==ej)&(sa==ea)&(sb==eb)&(col==ec)]
     witness.append({'source_class':ei,'source_canon':int(rows[ei]['canon']),'source_endpoint_scores':[ea,eb],'target_class':ej,'target_canon':int(rows[ej]['canon']),'delta':value,'native_count':native_amb.get(e,0),'native_example_mask':int(selected[0]) if len(selected) else None})
    counterfeit={'edits':witness,'merged_black_delta':[[list(e),v] for e,v in sorted(edge_delta.items()) if v]};break
  if counterfeit:break
 need(counterfeit is not None,'one-round hostile switch not found')
 artifact_hashes=json.loads((P/'artifact_hashes.json').read_text())
 for name,h in artifact_hashes.items():need(hashlib.sha256(data(name)).hexdigest()==h['sha256'],'artifact content hash')
 out={'status':'FINITE-EXACT complete n8 native ownership and wiring','class_source_lf_sha256':hashlib.sha256(source.read_text().encode()).hexdigest(),'signature_api_lf_sha256':hashlib.sha256(api_path.read_text(encoding='utf-8-sig').encode()).hexdigest(),'artifact_hashes':artifact_hashes,'producer_summary':json.loads((P/'producer_summary.json').read_text()),'independent_paths':total_paths,'hamilton_tiling_checksum':ham_checksum,'tilings':N,'closed_ports':1<<20,'port_groups':len(uk),'edges':edge_summaries,'v4_refinement':levels,'one_round_counterfeit':counterfeit,'checks':checks}
 if args.write:(P/'receipt.json').write_text(json.dumps(out,indent=2)+'\n')
 else:need(json.loads((P/'receipt.json').read_text())==json.loads(json.dumps(out)),'frozen receipt mismatch')
 print('checks',checks,'written' if args.write else 'matched',flush=True)
if __name__=='__main__':main()
