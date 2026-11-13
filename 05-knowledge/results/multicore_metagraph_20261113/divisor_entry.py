"""Target-driven entry into the empty-prefix, one-clock funnel sublanguage."""
from pathlib import Path
import sys,json,gzip,time
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'berggren_target_return_20261111'/'reset'))
import funnel_six_state as f
from rooted_growth import save

def run(N=15):
    targets=set()
    for a,b in [(462,-5),(455,59),(45,7),(33,7),(99,23)]:
        for depth in range(9):targets.add((a,b));a,b=3*a,3*b+2
    hits=[];tested=0;t0=time.time()
    for n in [1]+list(range(3,N+1)):
        for rank in range(1 if n==1 else f.a.class_counts(n)[n]):
            state=f.SINGLE if n==1 else f.compile_rank(n,rank);h,e,s,p,ee,t=state;G=s+p
            for a,b in targets:
                if a%G:continue
                alpha=a//G
                if alpha<1:continue
                right=tuple((alpha-1)>>i&1 for i in range((alpha-1).bit_length()))
                beta=f.a.tail_coefficients((1,)+right,1)[2]
                tested+=1
                if h*beta-ee-1==b:
                    hits.append(dict(A=a,B=b,order=n,rank=rank,state=state,left=(),middle=1,right=right,fixed_order=n+len(right)+2))
        print(n,'hits',len(hits),'divisor candidates',tested,'seconds',round(time.time()-t0,2),flush=True)
    save('divisor_entry',dict(scope=f'All unique-start strong cores through{N}, plus singleton; empty fixed prefix; middle1; every possible right word determined by target A; displayed target set',max_core_order=N,targets=sorted(targets),tested=tested,hits=hits))
    print('new targets',sorted({(r['A'],r['B']) for r in hits}))
if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 15)
