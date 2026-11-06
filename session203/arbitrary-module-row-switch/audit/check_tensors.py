"""Independent core partitions into actual paths, without auxiliary directed cores."""
from itertools import permutations
from math import gcd,factorial
from pathlib import Path
import json

CASES=[('A',[102,0,10,83,135,30,182,47],5),('C',[142,0,34,22,71,155,47,94],3)]
KEYS=('h','s','e','L','R','d','beta','B','V')

def need(b,s):
    if not b:raise RuntimeError(s)

def poly_add(*polys):
    ans=[0]*12
    for factor,p in polys:
        for j,v in enumerate(p):ans[j]+=factor*v
    return ans

def multiply(p,q):
    ans=[0]*12
    for i,a in enumerate(p):
        for j,b in enumerate(q):
            if a and b:need(i+j<12,'degree');ans[i+j]+=a*b
    return ans

def gap_path(a,replaced,word,root_first=False,root_last=False):
    # A path with no fixed core is empty or consists of one module path.
    if not word:return [1,1]+[0]*10
    if root_first and word[0]!=0:return [0]*12
    if root_last and word[-1]!=0:return [0]*12
    ans=[1]+[0]*11
    for k in range(len(word)+1):
        l=word[k-1] if k else None;r=word[k] if k<len(word) else None
        empty=int(l is None or r is None or a[l]>>r&1)
        filled=int((l is None or a[l]>>replaced&1) and (r is None or a[replaced]>>r&1))
        if root_first and k==0:filled=0
        if root_last and k==len(word):filled=0
        ans=multiply(ans,[empty,filled]+[0]*10)
    return ans

def profile(a,v):
    for i in range(len(a)):
        need(not(a[i]>>i&1),'loop')
        for j in range(i):need((a[i]>>j&1)+(a[j]>>i&1)==1,'physical tournament')
    keep=[i for i in range(len(a)) if i!=v];out={k:[0]*12 for k in KEYS};orders=0;cuts=0
    for p in permutations(keep):
        orders+=1
        for kind,rf,rl in [('h',False,False),('s',True,False),('e',False,True)]:
            out[kind]=poly_add((1,out[kind]),(1,gap_path(a,v,p,rf,rl)))
        root_index=p.index(0)
        # R=P followed by Q, with Q starting at the retained root.
        out['R']=poly_add((1,out['R']),(1,multiply(gap_path(a,v,p[:root_index]),gap_path(a,v,p[root_index:],True))))
        # L=P followed by Q, with P ending at the retained root.
        out['L']=poly_add((1,out['L']),(1,multiply(gap_path(a,v,p[:root_index+1],False,True),gap_path(a,v,p[root_index+1:]))))
        cuts+=2
    for p in permutations([i for i in keep if i!=0]):
        orders+=1;out['d']=poly_add((1,out['d']),(1,gap_path(a,v,p)))
        for cut in range(len(p)+1):
            cuts+=1;out['beta']=poly_add((1,out['beta']),(1,multiply(gap_path(a,v,p[:cut]),gap_path(a,v,p[cut:]))))
    out['B']=poly_add((1,out['h']),(1,out['L']),(-1,out['s']));out['V']=poly_add((1,out['R']),(1,out['beta']))
    return out,orders,cuts

def determinant(a):
    if len(a)==1:return a[0][0]
    return sum((-1)**j*a[0][j]*determinant([r[:j]+r[j+1:] for r in a[1:]]) for j in range(len(a)))

def main():
    producer=json.loads((Path(__file__).with_name('PRODUCER_TENSORS.json')).read_text());ans={};orders=cuts=0
    for name,a,v in CASES:
        out,o,c=profile(a,v);orders+=o;cuts+=c;ans[name]=out
        archived=producer['A' if name=='A' else 'B']
        for key in KEYS:need(out[key]==(archived[key]+[0]*12)[:12],('full coefficient array',name,key))
    delta={key:poly_add((1,ans['C'][key]),(-1,ans['A'][key])) for key in KEYS}
    expected=[0,2,2]+[0]*9
    for key in KEYS:need(delta[key]==(expected if key in ('h','s') else [0]*12),('exact difference including empty boundary',key))
    mat=[[ans['A'][key][j] for j in range(1,5)] for key in ('B','V','d')]
    kernel=[(-1)**j*determinant([row[:j]+row[j+1:] for row in mat]) for j in range(4)]
    g=0
    for x in kernel:g=gcd(g,abs(x))
    need(g>0,'rank three');kernel=[x//g for x in kernel]
    if kernel[0]>0:kernel=[-x for x in kernel]
    need(kernel==[-1728,3420,-5427,7768],'primitive integer kernel')
    endpoint=[sum(ans[name]['s'][j+1]*kernel[j] for j in range(4)) for name in ('A','C')]
    need(endpoint==[1944,5328],'endpoint changes')
    allowed=[k for k in range(1,25) if all(k*kernel[j]%factorial(j+1)==0 for j in range(4))]
    need(allowed==[6,12,18,24],'factorial lattice restriction')
    unordered=[6*kernel[j]//factorial(j+1) for j in range(4)]
    need(unordered==[-10368,10260,-5427,1942],'unordered primitive kernel')
    primitive=0
    for value in unordered:primitive=gcd(primitive,abs(value))
    need(primitive==1,'primitive unordered direction')
    permutation=[0,1,2,5,4,3,7,6]
    A=CASES[0][1];C=CASES[1][1]
    relabeled=[sum(1<<j for j,v in enumerate(permutation) if C[u]>>v&1) for u in permutation]
    changed=[(i,j) for i in range(8) for j in range(i) if bool(A[i]>>j&1)!=bool(relabeled[i]>>j&1)]
    need(changed==[(5,3)] and A[5]>>3&1,'single module-bundle flip after relabeling')
    print(json.dumps({'status':'PASS','method':'Partition actual core vertices between labeled paths; no auxiliary directed cores.','core_orders':orders,'root_or_separator_cuts':cuts,'tensors':ans,'delta':delta,'primitive_moment_kernel':kernel,'endpoint_changes':endpoint,'factorial_lattice_multiplier':6,'unordered_moment_kernel':unordered,'C_relabeling':permutation,'relabeled_C_rows':relabeled,'single_reversed_arc_before_inflation':[5,3]},indent=2,sort_keys=True))

if __name__=='__main__':main()
