"""Exact central-word arithmetic. No third-party modules; explicit checks.

python verify.py --cap 5001
python verify.py --header central_forms.h
g++ -O3 -std=c++17 inventory.cpp -o /tmp/central-inventory
/tmp/central-inventory 10000001
"""
import argparse
import hashlib
import heapq
import json
import math
from functools import lru_cache
from itertools import product
from pathlib import Path

T=(1,0,2,1)
Q=((0,0,0,1),(0,1,0,-1),(0,0,1,-1),(1,-1,-1,1))
P=((1,0,2,1),(0,0,1,0),(1,1,2,2),(0,0,1,1))
U=(1,1,2,2)
I=tuple(tuple(int(i==j) for j in range(4)) for i in range(4))
GAPS=[7,13,21,41,43,55,67,77,83,163,173,231,233,247,251,277,
      299,421,449,503,533,569,583,637,821,1067,1339,1861,2647]

def need(ok,message):
    if not ok: raise RuntimeError(message)

def transpose(A): return tuple(zip(*A))
def mm(A,B):
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(4))
                       for j in range(4)) for i in range(4))
def mv(A,v): return tuple(sum(x*y for x,y in zip(row,v)) for row in A)
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def step(s,c):
    D,E,V,d=s
    if c=='O': return D,0,2*V,d
    if c=='Z': return D-E+V,d,2*V-d,d
    if c=='J': return D-E+V,D-E,D+V,D
    raise ValueError(c)

GENERATORS={c:transpose(tuple(step(v,c) for v in I)) for c in 'OZJ'}
def matrix(w):
    M=I
    for c in w: M=mm(GENERATORS[c],M)
    return M
def state(w,s=T):
    for c in w: s=step(s,c)
    return s
def dagger(w): return ''.join({'O':'Z','Z':'O','J':'J'}[c] for c in w[::-1])
def K(v,w): return dot(v,mv(Q,w))
def F(s): return s[0]-s[1]+s[2]
def positive(s):
    D,E,V,d=s
    return D-E-d,V-D+E-d,E,d-E

@lru_cache(None)
def g(k):
    if k==1: return 1
    if k%2==0: return g(k//2)
    h=(k+1)//2
    return g(h)+h-h%2
def tau(k): return g(k),k%2,k,1
def unit_word(k):
    need(k>=2,'nontrivial unit index')
    return ''.join('Z' if c=='0' else 'O' for c in bin(k-1)[3:])
def central_words(limit=8):
    yield ''
    if limit>=1: yield 'J'
    for length in range(2,limit+1):
        for middle in product('OZJ',repeat=length-2):
            yield 'J'+''.join(middle)+'J'

def form_data(limit=8):
    """Matrix deduplication retains a shortest representative, by construction."""
    answer=[];seen=set()
    for w in central_words(limit):
        A=mm(Q,matrix(w))
        if A in seen: continue
        seen.add(A)
        C=mm(mm(transpose(P),A),P)
        lower4=min(4*C[i][j]//(U[i]*U[j]) for i in range(4) for j in range(4))
        # All denominators divide four, so this is exactly four times lambda.
        need(all((4*C[i][j])%(U[i]*U[j])==0 for i in range(4) for j in range(4)),
             'exact rational coefficient scale')
        need(lower4>=2,'universal half-product bound')
        answer.append((w,A,lower4))
    return answer

def row_for(A,u):
    return tuple(sum(tau(u)[i]*A[i][j] for i in range(4)) for j in range(4))
def value(row,v): return row[0]*g(v)+row[1]*(v%2)+row[2]*v+row[3]

def inverse(n,limit=8,max_small_index=None):
    """Complete progression inverse in the symmetric whole central-word family."""
    found=[]
    for w,A,b in form_data(limit):
        top=math.isqrt(4*n//b)
        if max_small_index is not None: top=min(top,max_small_index)
        for u in range(2,top+1):
            alpha,gamma,beta,delta=row_for(A,u)
            need(alpha>0 and beta>0,'positive arithmetic inverse coefficients')
            bound=4*n//(b*u)
            for parity in (0,1):
                lo=(u-parity+1)//2;hi=(bound-parity)//2
                h=math.gcd(2*beta,alpha)
                rhs=n-delta-(beta+gamma)*parity
                if rhs%h: continue
                modulus=alpha//h
                residue=0 if modulus==1 else (rhs//h)*pow(2*beta//h,-1,modulus)%modulus
                first=lo+(residue-lo)%modulus
                for t in range(first,hi+1,modulus):
                    v=2*t+parity
                    if value((alpha,gamma,beta,delta),v)==n:
                        found.append((w,u,v))
    return sorted(found,key=lambda z:(len(z[0]),z[0],z[1],z[2]))

def inventory(cap,limit=8):
    best={};witness={};tests=0
    for w,A,b in form_data(limit):
        for u in range(2,math.isqrt(4*cap//b)+1):
            row=row_for(A,u)
            for v in range(u,4*cap//(b*u)+1):
                n=value(row,v);tests+=1
                need(n>=3 and n%2==1,'physical odd output')
                if n<=cap and len(w)<best.get(n,100):
                    best[n]=len(w);witness[n]=(w,u,v)
    return best,witness,tests

def source_bank(cap):
    out={T};heap=[T]
    while heap:
        s=heapq.heappop(heap)
        for c in 'OZJ':
            q=step(s,c)
            if q[2]<=cap and q not in out:
                out.add(q);heapq.heappush(heap,q)
    return out

def graph(w):
    rows=[2,0];root=0
    for c in w:
        n=len(rows);rows.append(0)
        for i in range(n):
            if c!='O' and i==root: rows[n]|=1<<i
            else: rows[i]|=1<<n
        if c=='J':root=n
    return rows

def graph_counts(rows):
    n=len(rows);full=(1<<n)-1
    need(all(not (rows[i]>>i&1) for i in range(n)), 'no graph loops')
    need(all(((rows[i]>>j&1)+(rows[j]>>i&1))==1 for i in range(n) for j in range(i)),
         'literal tournament')
    def count(starts,closed):
        dp=[{} for _ in range(1<<n)]
        for i in starts:dp[1<<i][i]=1
        for mask in range(1,full+1):
            for last,multiplicity in dp[mask].items():
                todo=rows[last]&~mask&full
                while todo:
                    bit=todo&-todo;todo-=bit;j=bit.bit_length()-1
                    out=dp[mask|bit];out[j]=out.get(j,0)+multiplicity
        return sum(c for j,c in dp[full].items() if not closed or rows[j]&1)
    return count(range(n),False),count([0],True)

def graph_word(w,u,v): return unit_word(v)+w+dagger(unit_word(u))+'J'

def local_component(w):
    seen={w};todo=[w]
    for x in todo:
        for i in range(len(x)-1):
            if x[i:i+2] in ('JZ','OJ'):
                y=x[:i]+('OJ' if x[i:i+2]=='JZ' else 'JZ')+x[i+2:]
                if y not in seen:seen.add(y);todo.append(y)
    return seen

def normal_menu():
    kept=[];seen=set();certificates=[]
    for w in central_words():
        if w in seen:continue
        component=local_component(w);seen.update(component)
        need(all(matrix(x)==matrix(w) for x in component),'exact rewrite equality')
        loose=[x for x in component if x and (x[0]!='J' or x[-1]!='J')]
        if loose:
            x=min(loose);core=x.strip('OZ')
            need(len(core)<len(w),'strict boundary absorption')
            certificates.append((w,x,core))
        else:kept.append(min(component))
    return kept,certificates

def checks(cap):
    need(cap>=5001,'finite cap includes all displayed gaps')
    for k in range(2,4097):
        need(state(unit_word(k))==tau(k),'unit physical compiler')
        need(len(unit_word(k))==(k-1).bit_length()-1,'unit word length')
    for c in 'OZJ':
        need(mm(transpose(GENERATORS[c]),Q)==mm(Q,GENERATORS[dagger(c)]),'adjoint')
    forms=form_data();need(len(forms)==592,'distinct finite matrix count')
    all_matrices={mm(Q,matrix(w)) for w in central_words()}
    need(all(transpose(A) in all_matrices for A in all_matrices),'transpose-closed central menu')
    minima=[]
    for length in range(9):
        minima.append(min(b for w,A,b in forms if len(w)==length))
    need(minima==[2,4,4,8,12,20,36,60,100],'finite exact lower coefficient minima')
    for w,A,b in forms:
        C=mm(transpose(P),A)
        for column in (0,2):
            need(all(C[i][column]>=0 for i in range(4)) and C[2][column]>0 and C[3][column]>0,
                 'all-index positivity of alpha and beta')
        for u,v in ((2,2),(3,5),(11,23),(41,86)):
            n=value(row_for(A,u),v)
            need(n==F(state(graph_word(w,u,v)[:-1])),'chronological physical compiler')
            need(len(graph_word(w,u,v))+2 <= n.bit_length()-1+7,'uniform logarithmic graph order')
    kept,certificates=normal_menu()
    need(len(kept)==235 and len(certificates)==376,'normal-menu finite counts')
    best,witness,tests=inventory(cap)
    gaps=[n for n in range(3,cap+1,2) if n not in best]
    need(gaps==GAPS,'finite arithmetic coverage gap list')
    bank=source_bank(cap-1)
    direct={F(s) for s in bank if F(s)<=cap}
    need(set(best)==direct,'finite equality with complete forward physical source inventory')
    critical={}
    for n in (1877,17383):
        found=inverse(n)
        need(found and min(len(w) for w,u,v in found)==8,'complete length-seven obstruction')
        w,u,v=found[0];word=graph_word(w,u,v);rows=graph(word)
        need(graph_counts(rows)==(n,1),'literal critical C1 graph')
        critical[n]={'word':w,'left_index':u,'right_index':v,'graph_word':word,
                     'vertices':len(rows),'representations':found}
    # An independent arithmetic path for every target in a modest interval.
    inverse_tests=0
    for n in list(range(3,301,2))+GAPS+[1877,17383]:
        hit=inverse(n)
        need(bool(hit)==(n not in GAPS),'progression inverse membership control')
        for w,u,v in hit:need(value(row_for(mm(Q,matrix(w)),u),v)==n,'inverse exact replay')
        inverse_tests+=1
    repairs=json.loads((Path(__file__).parent/'dependencies/REPAIR_GRAPHS.json').read_text())
    need(sorted(r['H'] for r in repairs)==[n for n in GAPS if n not in (7,21)],'all finite repair targets')
    for r in repairs:need(graph_counts(r['rows'])==(r['H'],r['C']),'literal finite repair graph')
    return {'status':'PASS','cap':cap,'raw_words':sum(1 for _ in central_words()),
            'distinct_matrices':len(forms),'normal_menu_words':len(kept),
            'normal_menu_dominated_components':len(certificates),
            'lower4_minima_by_length':minima,'arithmetic_pair_tests':tests,
            'physical_source_count':len(bank),'finite_gaps':gaps,
            'progression_inverse_controls':inverse_tests,'critical_witnesses':critical,
            'repair_graphs':len(repairs),'universal_occupancy':'OPEN'}

def header(path):
    forms=form_data()
    s='// Generated exactly by verify.py; order is (D,E,U,d).\n'
    s+='const int FORMCOUNT='+str(len(forms))+';\n'
    s+='const int lengths[]={'+','.join(str(len(w)) for w,A,b in forms)+'};\n'
    s+='const int lower4[]={'+','.join(str(b) for w,A,b in forms)+'};\n'
    s+='const long long forms[][4][4]='+json.dumps([A for w,A,b in forms]).replace('[','{').replace(']','}')+';\n'
    Path(path).write_text(s)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--cap',type=int,default=5001)
    parser.add_argument('--header');args=parser.parse_args()
    if args.header:header(args.header)
    else:print(json.dumps(checks(args.cap),indent=2,sort_keys=True))
