from pathlib import Path
import json,sys,hashlib
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P/'source/session197-lowcycle'))
import verify_arithmetic_inverse as S

def need(b,m):
 if not b:raise RuntimeError(m)
def direct(L):
 if L==1:return 1,1,1
 bits=[1-((L-1)>>i&1)for i in range((L-1).bit_length())];g=1+bits[0];W=0
 for i,w in enumerate(bits):
  if not w:g+=W
  if i and bits[i-1]==0 and w:g+=1
  W=2*W+w
 Q=1+sum(2**(len(bits)-1-i)for i,w in enumerate(bits)if w)
 return g,Q,2**len(bits)
# Enumerate the actual left/right base assignment of every apex Hamilton word.
assignments=0
for L in range(1,2049):
 m=(L-1).bit_length();zeros=L-1;full=(1<<m)-1;ones=full^zeros;g=Q=Lv=0
 for left in range(1<<m):
  right=full^left;endok=not left or bool(zeros>>(left.bit_length()-1)&1);startok=not right or bool(ones&(right&-right))
  g+=endok and startok;Q+=startok;Lv+=endok;assignments+=1
 need((g,Q,1<<m)==direct(L)==S.gadget(L),'literal compiled gadget constants')
 need(Lv==L,'literal ending two-cover code')
# Direct finite-L evaluation versus the source modular inverse, including
# even virtual h and even targets, every e boundary, and modulus one.
controls=0;cone=0;profiles=0
for h in range(1,16):
 for e in range(h+1):
  for R in range(h,2*h+1):
   profiles+=1;expected={}
   for L in range(1,401//R+1):
    g,Q,B=direct(L);N=g*h-(L&1)*e+L*R
    if N<=401:expected.setdefault(N,[]).append(L)
    hp=N;ep=h-(L==1)*e;rp=Q*h+B*R
    need(hp<=rp<=2*(hp-ep),'formal response cone');cone+=1
   for N in range(1,402):
    need(S.inverse(N,(h,e,R))==expected.get(N,[]),('exact inverse',N,h,e,R));controls+=1
# Compare all sorted full triples, not just the source's prefix counts.
prefixes=[]
for B in (500,1000):
 ours={tuple(map(int,l.split(',')))for l in (P/f'image-normal-{B}-prefixes.csv').read_text().splitlines()}
 expected=set(S.prefixes(B));need(ours==expected,'complete prefix triples')
 prefixes.append({'B':B,'states':len(ours),'csv_sha256':hashlib.sha256((P/f'image-normal-{B}-prefixes.csv').read_bytes()).hexdigest()})
print(json.dumps({'status':'PASS','literal_gadgets':2048,'literal_base_assignments':assignments,'virtual_profiles':profiles,'arithmetic_inverse_controls':controls,'cone_controls':cone,'prefix_banks':prefixes},sort_keys=True,indent=2))
