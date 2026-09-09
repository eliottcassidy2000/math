"""Independent exact referee of the real trace-square-zero cubic theorem.

No primary source imports or execution. General signed real matrices are
enumerated directly, with separate weighted-triple and trace contractions.
"""
from pathlib import Path
from collections import Counter
from itertools import product, combinations
from fractions import Fraction as Q
from hashlib import sha256
import json
import sys

sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve().parent; STEM=Path(__file__).stem
PRIMARY='continuing15_20260909_universal_cubic_bound'
GATES=Counter()
PINS={'.py':'27af22e5785f3a82f71930381d5dc3b85a6d963c6e7744e1ad95e9fc0d8af6ad',
      '.out':'e29c6fb53c742fe08c8062848d72d0b7279b20173b3d2256f8e50480474b6121',
      '_certificate.json':'3dc2f6540cc76290256321b558fdc73076a2577df832a77eaa719d31998e1263'}

def need(ok,label):
    GATES[label]+=1
    if not ok:raise ArithmeticError('always-active independent gate: '+label)

def transpose(A):return [list(row) for row in zip(*A)]
def multiply(A,B):return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def trace(A):return sum(A[i][i] for i in range(len(A)))
def energy(A):return sum(x*x for row in A for x in row)
def eye(n):return [[int(i==j) for j in range(n)] for i in range(n)]
def scale(A,c):return [[c*x for x in row] for row in A]
def add(A,B,c=1):return [[a+c*b for a,b in zip(row,other)] for row,other in zip(A,B)]
def rank(A):
    B=[[Q(x) for x in row] for row in A];n=len(A);r=0
    for j in range(n):
        pivot=next((i for i in range(r,n) if B[i][j]),None)
        if pivot is None:continue
        B[r],B[pivot]=B[pivot],B[r];q=B[r][j];B[r]=[x/q for x in B[r]]
        for i in range(n):
            if i!=r and B[i][j]:
                q=B[i][j];B[i]=[a-q*b for a,b in zip(B[i],B[r])]
        r+=1
    return r

EQUALITY=Counter()
def inspect(A,label,require_zero_square=True):
    square=multiply(A,A);moment2=trace(square);cube=multiply(square,A);moment3=trace(cube);E=energy(A)
    if require_zero_square:need(moment2==0,label+' exact trace-square hypothesis')
    else:
        if moment2:return None
    need(3*moment3*moment3<=E**3,label+' absolute cubic inequality without radicals')
    if not E:
        need(all(x==0 for row in A for x in row),'zero-energy equality has exactly the zero matrix')
        EQUALITY['zero']+=1
    elif 3*moment3*moment3==E**3:
        B=A if moment3>0 else scale(A,-1);r=Q(abs(moment3),E)
        need(multiply(A,transpose(A))==multiply(transpose(A),A),'every enumerated equality is real normal')
        need(rank(A)==3 and trace(A)==0,'every nonzero enumerated equality has rank three and zero first trace')
        need(multiply(multiply(multiply(B,B),B),B)==scale(B,r**3),'exact equality spectral polynomial with positive scale')
        EQUALITY['positive' if moment3>0 else 'negative']+=1
    return E,moment3

def triple_check(A):
    n=len(A);C=T=0
    for i,j,k in combinations(range(n),3):
        vertices=(i,j,k);edges=[];out={v:0 for v in vertices}
        for a,b in combinations(vertices,2):
            if A[a][b]:edges.append(A[a][b]);out[a]+=1
            elif A[b][a]:edges.append(A[b][a]);out[b]+=1
            else:break
        if len(edges)!=3:continue
        w=edges[0]*edges[1]*edges[2]
        if all(out[v]==1 for v in vertices):C+=w
        else:T+=w
    square=multiply(A,A);moment3=trace(multiply(square,A))
    plus=add(A,transpose(A));minus=add(A,transpose(A),-1)
    F=Q(trace(multiply(plus,multiply(minus,minus))),2)
    need(moment3==3*C,'direct trace versus signed cyclic triple products')
    need(F==3*C-T,'direct symmetric-skew contraction versus cyclic-transitive production')
    if all(x>=0 for row in A for x in row):
        E=energy(A)
        need(T>=0 and F<=3*C,'nonnegative orientation pays the production comparison')
        need(27*C*C<=E**3,'every tested nonnegative oriented support satisfies universal cyclic bound')
        if F>0:need(3*F*F<=E**3,'positive production bound squared only after its sign is paid')
    return C,T,F

def encode(A):return [[str(x) for x in row] for row in A]

def main():
    filed=HERE.name=='04-computation'
    producer=HERE.parent/'05-knowledge/results' if filed else Path('C:/w/continuing15_20260909_cubic')
    for suffix,pin in PINS.items():
        path=(HERE if filed and suffix=='.py' else producer)/(PRIMARY+suffix)
        need(sha256(path.read_bytes()).hexdigest()==pin,'frozen primary '+suffix)
    full_counts=[];eligible_counts=[]
    for n,alphabet in ((1,range(-2,3)),(2,range(-2,3)),(3,range(-1,2))):
        total=eligible=0
        for flat in product(alphabet,repeat=n*n):
            total+=1;A=[list(flat[i*n:(i+1)*n]) for i in range(n)]
            # Independent direct trace-square filter includes diagonal terms
            # and cancellations between opposite signed entries.
            moment2=sum(A[i][i]**2 for i in range(n))+2*sum(A[i][j]*A[j][i] for i,j in combinations(range(n),2))
            if moment2:continue
            eligible+=1;inspect(A,'complete arbitrary real integer matrix')
        full_counts.append(total);eligible_counts.append(eligible)
    need(full_counts==[5,625,19683],'complete unrestricted integer matrix cubes')
    oriented_counts=[]
    for n in range(1,5):
        edges=list(combinations(range(n),2));count=0
        for choices in product(range(5),repeat=len(edges)):
            A=[[0]*n for _ in range(n)]
            for (i,j),choice in zip(edges,choices):
                if choice in (1,2):A[i][j]=1 if choice==1 else -1
                elif choice in (3,4):A[j][i]=1 if choice==3 else -1
            inspect(A,'complete signed oriented matrix');triple_check(A);count+=1
        oriented_counts.append(count)
    need(oriented_counts==[1,5,125,15625],'complete distinct signed-oriented matrices through order four')
    five_counts=0
    for weights in product((1,2),repeat=10):
        A=[[0]*5 for _ in range(5)]
        for (i,j),w in zip(combinations(range(5),2),weights):
            if (j-i)%5 in (1,2):A[i][j]=w
            else:A[j][i]=w
        inspect(A,'all weights on prime cyclic five-support');triple_check(A);five_counts+=1
    need(five_counts==1024,'complete ten-edge two-amplitude prime cyclic support')
    cyclic_support=[]
    for vertices in combinations(range(5),3):
        if all(sum(bool(A[v][u]) for u in vertices if u!=v)==1 for v in vertices):
            cyclic_support.append({(i,j) for i,j in product(vertices,repeat=2) if A[i][j]})
    need(len(cyclic_support)==5 and not set.intersection(*cyclic_support),'five-support cyclic triangles have no common edge')
    for size in range(2,5):
        for chosen in combinations(range(5),size):
            is_module=all(len({bool(A[v][u]) for u in chosen})==1 for v in range(5) if v not in chosen)
            need(not is_module,'five-support has no nontrivial homogeneous substitution block')
    P=[[0,1,0],[0,0,1],[1,0,0]]
    controls=[]
    def control(A,name):
        E,t3=inspect(A,name)
        controls.append({'name':name,'matrix':encode(A),'energy':str(E),'trace_cube':str(t3),
                         'normal':multiply(A,transpose(A))==multiply(transpose(A),A),
                         'equality':3*t3*t3==E**3})
        return E,t3
    control(P,'positive equal cycle');control(scale(P,-1),'negative equal cycle')
    old=[[0]*6 for _ in range(6)];old[0][1]=2
    for i in range(2,6):old[1][i]=old[i][0]=1
    E,t3=control(old,'incoming six-vertex multi-return equality')
    need(E==12 and t3==24 and rank(old)==3,'old non-triangle equality satisfies the complete spectral criterion')
    need(triple_check(old)==(8,0,24),'incoming multi-return production equality has exact zero transitive cost')
    control([[0,7,1],[0,0,3],[0,0,0]],'nonzero nilpotent zero-spectrum strictness')
    for numerator in range(-12,13):
        t=Q(numerator,4);a=1-t*t/2;b=1+t*t/2
        A=[[2*t,0,0],[0,-a,-b],[0,b,-a]]
        E,t3=control(A,'rational real-complex spectrum '+str(t))
    for n in (3,4,6):
        base=[[0]*n for _ in range(n)]
        for i in range(3):base[i][(i+1)%3]=1
        for j in range(1,7):
            v=[Q(((j+1)*(i+2))%7-3) for i in range(n)];norm=sum(x*x for x in v)
            if not norm:continue
            H=[[Q(int(i==l))-2*v[i]*v[l]/norm for l in range(n)] for i in range(n)]
            need(multiply(H,transpose(H))==eye(n),'exact rational Householder gauge is orthogonal')
            A=multiply(multiply(H,base),transpose(H));E,t3=control(A,'rational orthogonal image '+str((n,j)))
            need(E==3 and t3==3,'orthogonal images preserve equality beyond the oriented coordinate cone')
        for j in range(1,7):
            t=Q(j,3);S=eye(n);Sinv=eye(n);S[0][1]=t;Sinv[0][1]=-t
            need(multiply(S,Sinv)==eye(n),'exact unipotent similarity inverse')
            A=multiply(multiply(S,base),Sinv);E,t3=control(A,'nonnormal unipotent similarity '+str((n,j)))
            need(t3==3 and E>3 and 3*t3*t3<E**3,'similarity preserves cubic moment but its Frobenius defect is strictly paid')
    two_cycles=[[0]*6 for _ in range(6)]
    for start in (0,3):
        for i in range(3):two_cycles[start+i][start+(i+1)%3]=1
    E,t3=control(two_cycles,'two positive normal cyclic blocks')
    need(E==6 and t3==6 and 3*t3*t3<E**3,'normality alone and multiple active spectral triples are insufficient for equality')
    signed=[[0,1,-1],[0,0,1],[0,0,0]]
    control(signed,'signed transitive orientation outside nonnegative production scope')
    C,T,F=triple_check(signed)
    need(C==0 and T==-1 and F==1,'outside-cone transitive sign refutes F less than or equal to trace cube')
    need(3*trace([[1]])**2>energy([[1]])**3,'omitting trace-square-zero hypothesis fails already on a one by one real matrix')
    certificate={'status':'INDEPENDENTLY ACCEPTED real trace-square-zero cubic inequality and oriented consequence',
        'reviewed_primary_report_sha256':'2ecfe0ee8fd7fe2aa65113a7e6c29df0f2214cefa45e739444cfcb92fb2c1bd6',
        'primary_executable_pins':PINS,'full_real_matrix_universe':full_counts,'trace_square_zero_counts':eligible_counts,
        'complete_signed_oriented_counts':oriented_counts,'prime_cyclic_five_amplitude_bank':five_counts,
        'named_exact_controls':controls,'all_tested_equality_counts':dict(EQUALITY),
        'gate_counts':dict(sorted(GATES.items())),'total_gates':sum(GATES.values()),
        'scope':'Every finite-dimensional real matrix analytically; exact declared finite controls. Oriented production needs nonnegative amplitudes and no opposing entries. No PDE trajectory or global priority claim.'}
    target=HERE.parent/'05-knowledge/results' if filed else HERE
    (target/(STEM+'_certificate.json')).write_bytes((json.dumps(certificate,sort_keys=True,separators=(',',':'))+'\n').encode())
    print('Independent universal real cubic and oriented-support referee: PASS')
    print('COMPLETE_REAL_MATRIX_UNIVERSES',full_counts,'TRACE_SQUARE_ZERO',eligible_counts)
    print('COMPLETE_SIGNED_ORIENTED',oriented_counts,'PRIME_FIVE_SUPPORT_AMPLITUDES',five_counts)
    print('NAMED_RATIONAL_CONTROLS',len(controls),'EQUALITY_COUNTS',dict(EQUALITY))
    print('TOTAL_ALWAYS_ACTIVE_GATES',sum(GATES.values()))

if __name__=='__main__':main()
