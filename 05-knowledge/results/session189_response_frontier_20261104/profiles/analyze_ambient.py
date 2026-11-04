from analyze_models import *
M=readmodels(json.load(open(D/'cycle_models.json')),(D/'cycle_values.txt').read_text());json.dump(M,open(D/'cycle_banks.json','w'),separators=(',',':'))
S=json.load(open(D/'affine_banks.json'));out=[];bs=M[0]['block_masks'];pos={a:i for i,a in enumerate(bs)};witness=[]
def tovalues(coeff):return [sum(coeff[pos[s]*29+pos[t]] for s in subs(a) for t in subs(c)) for a in bs for c in bs]
for k,case in enumerate(S):
 A,B=case['full_run_coefficients'];K={key:[[a-b for a,b in zip(x,y)] for x,y in zip(A[key],B[key])] if key!='K0' else [a-b for a,b in zip(A[key],B[key])] for key in A}
 base=tovalues(case['all_pairs_contrast']);cols=[tovalues([4*x for x in v]) for v in K['K2']]+[tovalues([2*x for x in v]) for v in K['K3']]
 z=[(a-b)%2 for a,b in zip(flatten(M[2*k]['named_values']),flatten(M[2*k+1]['named_values']))];safe=[i for i,v in enumerate(z) if v==0]
 rank,sol=solve(bitvec([base[i] for i in safe]),[bitvec([c[i] for i in safe]) for c in cols]);vectors=[bitvec(c) for c in cols];bb=bitvec(base);sf=sum(1<<i for i in safe);blind=[];safe_blind=[];records=[]
 for mask in range(1<<14):
  bits=bb
  for i,col in enumerate(vectors):
   if mask>>i&1:bits^=col
  if not bits:blind.append(mask)
  if not bits&sf:safe_blind.append(mask)
  odd=(mask>>7).bit_count()%2
  usable=bits&sf if odd else bits
  ix=(usable&-usable).bit_length()-1 if usable else -1
  if ix<0 and not(k==1 and mask==127):raise RuntimeError(('unexpected uncovered type',k,mask))
  records.append([mask,odd,ix,None if ix<0 else bs[ix//29],None if ix<0 else bs[ix%29],0 if ix<0 else 4])
 out.append(dict(Q=case['Q'],cycle_parity_contrast=z,safe_indices=safe,safe_rank=rank,safe_blind_solution=sol,safe_blind_count=len(safe_blind),safe_blind_masks=safe_blind,local_blind_masks=blind))
 witness.append(dict(Q=case['Q'],columns=['type_code_b_then_c','total_order_odd','basis_index','a','c','normalized_gap_mod8'],records=records))
 print(case['Q'],'safe',len(safe),'rank',rank,'safe_blind',safe_blind,'local_blind',blind,flush=True)
json.dump(out,open(D/'ambient_analysis.json','w'),separators=(',',':'));json.dump(witness,open(D/'type_witnesses.json','w'),separators=(',',':'))
