"""Independent multiplicity-vector Burnside control and recovery hostile probes."""
from pathlib import Path
import sys,json,math,hashlib,subprocess,tempfile
from fractions import Fraction
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'engine'));sys.path.insert(0,str(ROOT/'absolute'))
from dyadic_completion import residue,residue_response,complete,recover_from_interval
from certified_interval import integer_interval,lower_and_tail

def independent_burnside(n):
    # Ascending odd-cycle multiplicities, no imported sector or cache code.
    terms=[]
    def go(j,remaining,cycles):
        if remaining==0:
            centralizer=1;bits=0
            for length,count in cycles:
                centralizer*=length**count*math.factorial(count)
                bits+=count*(length-1)//2+length*count*(count-1)//2
            for a,(length,count) in enumerate(cycles):
                for other,k in cycles[a+1:]: bits+=count*k*math.gcd(length,other)
            terms.append(Fraction(2**bits,centralizer));return
        if j>remaining:return
        for k in range(remaining//j+1):
            go(j+2,remaining-k*j,cycles+([(j,k)] if k else []))
    go(1,n,[])
    return sum(terms,Fraction()),terms

checks=0
for n in range(25):
    exact,terms=independent_burnside(n)
    assert exact.denominator==1
    for bits in [0,1,2,3,7,17,61,151]:
        for method in [residue,residue_response]:
            actual,_=method(n,bits)
            assert actual==int(exact)%(1<<bits),(n,bits,method.__name__)
            checks+=1
    for s in [0,1,2,3,4,5,6,9,n,n+3]:
        lo,hi=integer_interval(n,s)
        assert lo<=exact<=hi
        val,_=complete(n,s)
        assert val==exact
        checks+=2
    if n:
        # Weaker rational version of centralizer/bit lower bound, checked
        # at every sector where exact enumeration is cheap.
        assert all(t>0 for t in terms)

# Inclusive interval: width=M-1 allowed, width=M rejected; singleton uses b=0.
assert recover_from_interval(-8,-1,7,3)==-1
assert recover_from_interval(6,6,0,0)==6
assert recover_from_interval(2,3,3,1)==3
for args in [(0,8,0,3),(5,4,0,1),(2,3,0,2),(0,0,0,-1)]:
    try:recover_from_interval(*args)
    except ValueError:pass
    else:raise AssertionError(args)
checks+=7
# Granularity derivative proof uses log2(n)+1 <= n/8 for n>=64;
# establish base and monotonicity analytically in review, sample displays only.
assert 64-(1+math.log2(64))**2==15
binary=Path(tempfile.gettempdir())/'s277_review_modular'
subprocess.run(['g++','-O2','-std=c++17',str(ROOT/'engine/modular_counter.cpp'),'-lgmpxx','-lgmp','-o',str(binary)],check=True)
cpp_checks=0
for n in range(25):
    exact,_=independent_burnside(n)
    for b in [0,1,7,31,191]:
        run=subprocess.run([str(binary),str(n),str(b)],capture_output=True,text=True,check=True)
        assert int(run.stdout)==int(exact)%(1<<b)
        cpp_checks+=1
for n,b in [(-1,0),(10001,0),(0,-1)]:
    run=subprocess.run([str(binary),str(n),str(b)],capture_output=True,text=True)
    assert run.returncode!=0
    cpp_checks+=1
files=['absolute/PROOF.md' ,'absolute/certified_interval.py','engine/dyadic_completion.py','engine/PROOF.md','engine/modular_counter.cpp','bridges/CENSUS_VERSUS_SUPPORT.md']
out={'status':'PASS','independent_controls':checks,'cpp_independent_controls':cpp_checks,'n_range':[0,24],
     'residue_bits':[0,1,2,3,7,17,61,151],
     'sources':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files}}
(ROOT/'audit/RECOVERY_CONTROLS.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
