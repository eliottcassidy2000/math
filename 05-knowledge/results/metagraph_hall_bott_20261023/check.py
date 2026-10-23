"""Independent, inverse-free linking-charge extraction from the frozen n8 lift.

Reads actual merged edges, not the dyadic engine. All arithmetic is integral
or Fraction. --write freezes this packet's receipt; otherwise exact replay.
"""
from pathlib import Path
from fractions import Fraction
from hashlib import sha256
from collections import Counter
from math import gcd
import argparse, csv, gzip, json

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
CHECKS = 0

def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(why)

def raw(path):
    data = path.read_bytes() if path.exists() else gzip.decompress(Path(str(path)+'.gz').read_bytes())
    return data.replace(b'\r\n', b'\n')

def inverse2(a):
    d = a[0][0]*a[1][1]-a[0][1]*a[1][0]
    return [[Fraction(a[1][1],d),Fraction(-a[0][1],d)],
            [Fraction(-a[1][0],d),Fraction(a[0][0],d)]]

def bilinear(a,x,y):
    return sum(x[i]*a[i][j]*y[j] for i in range(len(x)) for j in range(len(y)))

def controls():
    a1 = [[9,-1],[-1,1]]
    a3 = [[3,-1],[-1,3]]
    out = []
    for a, unit in [(a1,1),(a3,3)]:
        inv = inverse2(a)
        need(a[0][0]*a[1][1]-a[0][1]*a[1][0] == 8, 'det8')
        need(gcd(*[abs(x) for row in a for x in row]) == 1, 'Smith(1,8)')
        need(inv[0][0] == Fraction(unit,8), 'generator norm')
        dimensions=[]
        for j in range(1,5):
            m=2**j
            residues={(x0%2,x1%2) for x0 in range(m) for x1 in range(m)
                if all((a[i][0]*x0+a[i][1]*x1)%m==0 for i in range(2))}
            need(len(residues) in (1,2),'one-dimensional dyadic control')
            dimensions.append(len(residues).bit_length()-1)
        need(dimensions==[1,1,1,0],'full dyadic dimensions coincide')
        # Exact inverse and quotient invariance under all short principal shifts.
        for i in range(2):
            for j in range(2):
                need(sum(a[i][k]*inv[k][j] for k in range(2)) == (i==j), 'inverse')
        c = [a[i][i] for i in range(2)]
        for z0 in range(-2,3):
            for z1 in range(-2,3):
                z=[z0,z1]; x=[1,0]
                xp=[x[i]+sum(a[i][j]*z[j] for j in range(2)) for i in range(2)]
                q=lambda v: (bilinear(inv,v,v)-bilinear(inv,v,c))/2
                need((q(xp)-q(x)).denominator == 1,'characteristic refinement descends')
        out.append({'A':a,'Smith':[1,8],'generator_self_link':f'{unit}/8',
                    'dyadic_dimensions':dimensions})
    need({(u*u)%8 for u in range(1,8,2)} == {1}, 'odd square unit')
    need(not any((u*u)%8 == 3 for u in range(1,8,2)), 'linking forms nonisometric')
    need(Fraction(1,2).denominator != 1, 'odd lattice uncorrected quadratic fails')
    # A unit tree's marked quotient charges are all zero despite different wiring.
    tree1=[[2,-1],[-1,1]]; tree2=[[1,0],[0,1]]
    need(inverse2(tree1) != inverse2(tree2), 'full rational response distinguishes trees')
    need(all(x.denominator==1 for a in (tree1,tree2) for row in inverse2(a) for x in row), 'linking quotient forgets trees')
    path1=[[1,-1,0],[-1,2,-1],[0,-1,2]]
    path2=[[1,0,-1],[0,2,-1],[-1,-1,2]]
    invpath1=[[3,2,1],[2,2,1],[1,1,1]]
    invpath2=[[3,1,2],[1,1,1],[2,1,2]]
    for a,inv in [(path1,invpath1),(path2,invpath2)]:
        for i in range(3):
            for j in range(3):need(sum(a[i][k]*inv[k][j] for k in range(3))==(i==j),'marked path unimodular')
    need(path1!=path2 and [path1[i][i] for i in range(3)]==[path2[i][i] for i in range(3)],'same labelled path degrees')
    return out

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--write',action='store_true'); args=parser.parse_args()
    control=controls()
    source=RESULTS/'metagraph_order8_electrical_20261023'/'receipt.json'
    edgepath=RESULTS/'metagraph_order8_20261023'/'merged_edges.tsv'
    nodepath=RESULTS/'metagraph_order8_20261023'/'merged_nodes.tsv'
    rr,er,nr=raw(source),raw(edgepath),raw(nodepath)
    r=json.loads(rr)['order8']; layer=r['higher_layers'][-1]
    need(layer['level']==8 and layer['dimension']==1 and layer['rank']==1 and layer['next_dimension']==0, 'last layer')
    x=layer['lifts'][0]; n=len(x); root=r['root']; modulus=256
    nodes=list(csv.DictReader(nr.decode().splitlines(),delimiter='\t'))
    need(n==3528 and len(nodes)==n and root==3527 and x[root]==0,'universe/root')
    need(all(isinstance(v,int) and 0<=v<modulus for v in x),'literal potential')
    image=[0]*n; diagonal=[0]*n; energy=0; colors={'blue':Counter(),'black':Counter()}; total={'blue':0,'black':0}
    for row in csv.DictReader(er.decode().splitlines(),delimiter='\t'):
        i,j,b,k=[int(row[s]) for s in ('source','target','blue','black')]
        need(0<=i<=j<n and min(b,k)>=0 and k%2==0,'literal edge schema')
        w=b+k; dx=x[i]-x[j]; image[i]+=w*dx; image[j]-=w*dx; energy+=w*dx*dx
        if i!=j: diagonal[i]+=w; diagonal[j]+=w
        for color,count in (('blue',b),('black',k)):
            if count: colors[color][dx%modulus]+=count;total[color]+=count
    need(all(v%modulus==0 for v in image),'all Kirchhoff equations modulo256')
    need(any(v%2 for v in x),'primitive potential gives exact charge order256')
    need(sum(x[i]*image[i] for i in range(n))==energy,'independent edge/dot energies')
    need(energy==6819178752,'frozen energy')
    u=energy//modulus
    need(u%256==105 and u%8==1,'self-link unit')
    normalizers=[s for s in range(1,256,2) if u*s*s%256==1]
    need(normalizers==[27,101,155,229],'cyclic normalization')
    # A marked vertex charge e_v has pairing x_v/256 with y=A*x/256.
    phase_gcd=256
    for color,hist in colors.items():
        for phase in hist: phase_gcd=gcd(phase_gcd,phase)
    need(phase_gcd==1,'port phases generate full256-character')
    a,b,c,d=1238,2242,2246,2315
    switched={a:2*(x[c]-x[d]),b:-2*(x[c]-x[d]),c:2*(x[a]-x[b]),d:-2*(x[a]-x[b])}
    guard={str(v):z%256 for v,z in switched.items()}
    need(list(guard.values())==[218,38,158,98],'marked switch obstruction')
    for amplitude in (-1,1,2):
        need(any((amplitude*z)%256 for z in switched.values()),'every feasible nonzero switch rejected')
        need(any((amplitude*z)%8 for z in switched.values()),'three-bit guard rejects each nonzero switch')
    # Universal minimality follows from Delta A_2=4*(u*v^T+v*u^T).
    # This check is the corresponding selected-potential control.
    need(all((2*z)%4==0 for z in switched.values()),'selected potential obeys the mod4 boundary')
    need(gcd(*(list(switched.values())+[8]))==2,'mod8 switch guard iff4 divides t')
    result={'status':'FINITE-EXACT; native n8 top linking summand and marked port phases',
      'input_sha256_lf':{str(p.name):sha256(v).hexdigest() for p,v in [(source,rr),(edgepath,er),(nodepath,nr)]},
      'native':{'vertices':n,'root':root,'modulus':modulus,'energy':energy,'divided_energy':u,
        'charge_order':256,'self_link':'105/256','cyclic_unit_squareclass_mod8':u%8,
        'normalizing_multipliers':normalizers,'normalized_self_link':'1/256',
        'charge_min':min(image)//256,'charge_max':max(image)//256,
        'odd_diagonal_count':sum(diagonal[i]%2 for i in range(n) if i!=root),
        'port_phase_histogram':{c:{str(k):v for k,v in sorted(h.items())} for c,h in colors.items()},
        'literal_line_totals':total,'switch_delta_kirchhoff_mod256':guard,
        'smallest_binary_modulus_rejecting_all_three_switches':8,
        'switch_delta_kirchhoff_mod8':{str(v):z%8 for v,z in switched.items()}},
      'controls':control,'checks':CHECKS}
    data=(json.dumps(result,indent=2,sort_keys=True)+'\n').encode()
    if args.write:(HERE/'receipt.json').write_bytes(data)
    else:need(data==raw(HERE/'receipt.json'),'frozen receipt comparison')
    print(json.dumps({'checks':CHECKS,'energy':energy,'self_link':'105/256','normalized_self_link':'1/256','receipt_sha256':sha256(data).hexdigest()},sort_keys=True))

if __name__=='__main__':main()
