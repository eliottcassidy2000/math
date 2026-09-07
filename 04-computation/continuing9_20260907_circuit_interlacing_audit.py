"""Independent fixed-moment C-only/full-C-D referee; no producer import."""
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
from math import comb
from pathlib import Path
import json
import sys
sys.stdout.reconfigure(newline="\n")
HERE=Path(__file__).resolve().parent
G=0
def gate(ok,label):
    global G
    G+=1
    if not ok: raise RuntimeError(label)
def trim(f):
    f=list(map(Q,f))
    while len(f)>1 and not f[-1]: f.pop()
    return f
def ev(f,x):
    out=Q(0)
    for a in reversed(f): out=out*x+a
    return out
def deriv(f): return trim([i*f[i] for i in range(1,len(f))] or [0])
def rem(f,g):
    f=trim(f); g=trim(g)
    while f!=[0] and len(f)>=len(g):
        n=len(f)-len(g); factor=f[-1]/g[-1]
        for j,a in enumerate(g): f[n+j]-=factor*a
        f=trim(f)
    return f
def sturm(f):
    a=trim(f); b=deriv(a); out=[a,b]
    while b!=[0]:
        c=[-q for q in rem(a,b)]
        if c==[0]: break
        c=[q/abs(c[-1]) for q in c]
        out.append(c); a,b=b,c
    return out
def variation(values):
    signs=[1 if v>0 else -1 for v in values if v]
    return sum(a!=b for a,b in zip(signs,signs[1:]))
def var_at(chain,x):
    return variation([ev(f,x) for f in chain])
def count(chain,a,b):
    va=var_at(chain,a)
    vb=variation([f[-1] for f in chain]) if b is None else var_at(chain,b)
    return va-vb
def reflected(f): return [(-1)**k*a for k,a in enumerate(f)]
def center(d,q): return [Q(comb(d,k))*q**(k*(k-1)//2) for k in range(d+1)]
def ratios(e):
    d=len(e)-1; h=[e[k]/comb(d,k) for k in range(d+1)]
    R=[h[k]*h[k]/(h[k-1]*h[k+1]) for k in range(1,d)]
    C=[R[k]/R[k-1] for k in range(1,len(R))]
    return R,C
def reconstruct(d,q,cs):
    h=[Q(1),Q(1)]; rk=1/q
    for k in range(1,d):
        if k>=2: rk*=cs[k-2]
        h.append(h[-1]*h[-1]/(rk*h[-2]))
    return [comb(d,k)*h[k] for k in range(d+1)]
def refine(f,chain,interval):
    lo,hi=interval
    if ev(f,hi)==0: return hi,hi
    mid=(lo+hi)/2
    return (lo,mid) if count(chain,lo,mid) else (mid,hi)
def isolate_positive(f):
    chain=sturm(f); degree=len(f)-1
    gate(count(chain,Q(0),None)==degree,'fresh exact positive root count')
    upper=Q(1)
    while count(chain,Q(0),upper)<degree: upper*=2
    pending=[(Q(0),upper)]; intervals=[]
    while pending:
        lo,hi=pending.pop(); n=count(chain,lo,hi)
        if n==0: continue
        if n==1:
            for _ in range(34): lo,hi=refine(f,chain,(lo,hi))
            intervals.append((lo,hi)); continue
        mid=(lo+hi)/2; pending.extend(((mid,hi),(lo,mid)))
    return sorted(intervals)

def model(x,y,z):
    return ([-z,y,-x,Q(55),Q(-13),Q(1)],
            [Q(3,7)*y,-Q(2,3)*x,Q(45),Q(-12),Q(1)],
            [Q(1,7)*y,-Q(5,12)*x,Q(36),Q(-11),Q(1)])

def imul(a,b):
    values=[x*y for x in a for y in b]
    return min(values),max(values)

def iev(f,interval):
    result=(Q(0),Q(0))
    for c in reversed(f):
        lo,hi=imul(result,interval)
        result=(lo+c,hi+c)
    return result

def sign_interval(f,interval):
    lo,hi=iev(f,interval)
    return 1 if lo>0 else -1 if hi<0 else 0

def verify_C_geometry(x,y,z):
    B,C,D=model(x,y,z)
    br=isolate_positive(B); cr=isolate_positive(C)
    gate(len(br)==5 and len(cr)==4,'five B/four C positive roots')
    gate(all(br[i][1]<cr[i][0] and cr[i][1]<br[i+1][0] for i in range(4)),
         'strict C interlacing by independently isolated polynomial roots')
    gate(all(sign_interval(C,br[i])==(-1)**i for i in range(5)),
         'positive C/B residues at every native root')
    return B,C,D,br,cr

q=Q(275,338); mu=Q(13,5)
anchor=[Q(comb(5,k))*mu**k*q**(k*(k-1)//2) for k in range(6)]
gate(anchor[:3]==[1,13,55],'actual fixed first two moments')
gate(ratios(anchor)==([Q(338,275)]*4,[Q(1)]*3),'actual all-one circuit center')
B,C,D,anchor_roots,_=verify_C_geometry(*anchor[3:])
gate(sign_interval(D,anchor_roots[0])==-1,'all-one center strictly fails D/B at smallest B root')

def locate(name):
    for p in (HERE/name,HERE.parent/'05-knowledge/results'/name,
              Path('C:/w/continuing9_20260907_moments')/name):
        if p.is_file(): return p
    raise FileNotFoundError(name)

STEM='continuing9_20260907_interlacer_circuits'
for suffix,pin in (('.py','a2d27b5925c3992d23089405bf41911936eac97c31b346a80664129b8fd3c369'),
                   ('.out','0b3cc97abe03805dd757b21b2e88b364ae83c820bcd16e51a065ecffcd0b3841')):
    raw=locate(STEM+suffix).read_bytes()
    gate(sha256(raw).hexdigest()==pin,'frozen sidecar '+suffix)
    if suffix=='.out': gate(b'\r' not in raw,'actual producer raw LF')
data_raw=locate(STEM+'_certificate.json').read_bytes()
gate(sha256(data_raw).hexdigest()=='dfbae81c3078b21042cd4a33cc0717860caec3d9b858fdca461b136225a84657',
     'frozen sidecar certificate')
data=json.loads(data_raw)
main_raw=locate('continuing9_20260907_fixed_moment_circuits_certificate.json').read_bytes()
gate(sha256(main_raw).hexdigest()==data['main_certificate_sha256']=='291b62a5638f8a057ccf514851cb0b5d100a2a61a3da8eff9151b49040be11b9',
     'retained fixed-moment carrier certificate')

def quotient_moments(B,A):
    den=list(reversed(B)); num=list(reversed(A))
    moments=[]
    for k in range(9):
        value=num[k] if k<len(num) else Q(0)
        value-=sum(den[j]*moments[k-j] for j in range(1,min(k,len(den)-1)+1))
        moments.append(value)
    # Literal coefficient multiplication, including all five recurrence terms.
    for k in range(9):
        gate(sum(den[j]*moments[k-j] for j in range(min(k,len(den)-1)+1))==
             (num[k] if k<len(num) else 0),'literal rational-function normalization')
    return moments

def leading_determinants(moments):
    # Symmetric rational LDL elimination, no symbolic determinant routines.
    A=[[moments[i+j] for j in range(5)] for i in range(5)]
    answer=[]; determinant=Q(1)
    for k in range(5):
        pivot=A[k][k]
        gate(pivot!=0,'nonzero exact symmetric elimination pivot')
        determinant*=pivot; answer.append(determinant)
        for i in range(k+1,5):
            for j in range(i,5):
                A[i][j]-=A[i][k]*A[k][j]/pivot
                A[j][i]=A[i][j]
    return answer

def certificate_row(row,kind):
    xyz=list(map(Q,row['xyz'])); es=[Q(1),Q(13),Q(55)]+xyz
    gate(all(v>0 for v in xyz),'strictly positive native coefficient packet')
    rr,cs=ratios(es)
    gate(cs==list(map(Q,row['C'])),'literal actual circuit ratios')
    gate([1 if c>1 else -1 if c<1 else 0 for c in cs]==row['word'],
         'literal circuit signs including exact ties')
    gate(rr[0]==Q(338,275),'fixed first-ratio carrier')
    B,C,D,br,cr=verify_C_geometry(*xyz)
    for label,A in (('C',C),('D',D)):
        moments=quotient_moments(B,A)
        dets=leading_determinants(moments)
        gate(moments==list(map(Q,row[label+'_moments'])),'full nine '+label+' moments')
        gate(dets==list(map(Q,row[label+'_leading_minors'])),'all five '+label+' leading minors')
        if label=='C' or kind=='full':
            gate(all(d>0 for d in dets),'strict '+label+' Gram positivity')
        else:
            gate(dets[3]<0,'explicit negative D four-by-four principal minor')
    if kind=='full':
        dr=isolate_positive(D)
        gate(all(br[i][1]<dr[i][0] and dr[i][1]<br[i+1][0] for i in range(4)),
             'full D interlacing by independently isolated polynomial roots')
        gate(all(sign_interval(D,br[i])==(-1)**i for i in range(5)),
             'strict positive D/B residues at all five native roots')
    else:
        gate(sign_interval(D,br[0])==-1,'strict D failure at the first native B root')
    return xyz

gate(certificate_row(data['center'],'C-only')==anchor[3:],
     'sidecar center equals independently reconstructed center')
seen=set()
for row in data['C_only_27']:
    word=tuple(row['word']); seen.add(word)
    actual=certificate_row(row,'C-only')
    coefficients=reconstruct(5,q,[1+Q(w,2048) for w in word])
    expected=[coefficients[k]*mu**k for k in range(3,6)]
    gate(actual==expected,'exact earlier ternary radius targets, with actual scaling')
gate(seen==set(product((-1,0,1),repeat=3)) and len(data['C_only_27'])==27,
     'complete ternary word universe exactly once')
full_words=set()
for row in data['full_model_four']:
    certificate_row(row,'full'); full_words.add(tuple(row['word']))
gate(full_words=={(1,1,1),(-1,1,1),(1,-1,1),(1,1,-1)} and len(data['full_model_four'])==4,
     'exact stated four full-model words only')
print('PASS independent C-only/full-interlacer circuit audit')
print('Native rational Sturm root ordering: center, all 27 ternary targets, and four full-model points.')
print('All 576 moments and 320 leading determinants match independent convolution/LDL.')
print('Every C-only control has an explicit negative D/B residue at its smallest B root.')
print('A C-only neighborhood exists; no entire radius cube or omitted full-word exclusion is inferred.')
print('Always-active exact gates:',G)
