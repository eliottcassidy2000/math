#!/usr/bin/env python3
from verify import *

def brute(r):
 return sum(all(r[a]>>b&1 for a,b in zip(p,p[1:]))for p in permutations(range(len(r))))
def main():
 r=[24,121,235,48,4,145,185,27];H,C,h=boundary(r)
 need(strong(r) and H==brute(r)==241,'parent')
 need(literal_boundary(r)==(H,C),'literal pair capacities')
 states=[]
 for S in range(256):
  child=extend(r,S);d=hp(child)-H;t=sum((r[i]&~S).bit_count()for i in bits(S))
  need(d==cut(r,C,S),'every child independent count')
  need(strong(child)==(0<S<255),'all strong cuts')
  if 0<S<255:states.append((t,d,S))
 need(min(t for t,d,S in states)==1,'one triangle exists')
 need(min(d for t,d,S in states if t==1)==224,'one triangle floor')
 need([(t,d,S)for t,d,S in states if d==min(x[1]for x in states)]==[(2,220,8),(2,220,253)],'all minimizers')
 for S,value in [(8,461),(253,461),(16,465)]:need(brute(extend(r,S))==value,'literal child')
 best,records=minstrong(r,C);need(best[0]==220,'all pair minimum')
 A={5,6,8,9,10,12,15,16};D={0};holes={}
 for degree in range(1,5):
  D={x+y for x in D for y in A};holes[degree]=sorted(set(range(5*degree,16*degree+1))-D)
 need(holes=={1:[7,11,13,14],2:[29],3:[],4:[]},'toric holes')
 out={'status':'PASS','rows':r,'H':H,'capacities':C,'all_256_children_checked':True,'all_56_terminal_pairs_certified':len(records)==56,'flow_records':records,'one_triangle_states':[list(x)for x in states if x[0]==1],'minimum_states':[list(x)for x in states if x[1]==220],'literal_child_counts':{'8':461,'253':461,'16':465},'toric_holes':holes}
 raw=json.dumps(out,indent=2,sort_keys=True)+'\n';Path(__file__).with_name('hostile_receipt.json').write_text(raw)
 print(json.dumps({'status':'PASS','minimum_increment':220,'minimum_one_triangle_increment':224,'terminal_pairs':len(records),'receipt_sha256':hashlib.sha256(raw.encode()).hexdigest()},sort_keys=True))
if __name__=='__main__':main()
