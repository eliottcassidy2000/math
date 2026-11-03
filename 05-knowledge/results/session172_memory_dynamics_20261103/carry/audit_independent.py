import itertools, json
checks=0
for p in (2,3,5,7):
    for a,d in itertools.product(range(p),repeat=2):
        for entries in itertools.product(range(p), repeat=4):
            p11,p12,p21,p22=entries
            outputs=[]
            for b in range(p):
                # Direct matrix diagonal of P*G, not author mul routine.
                left=p11*a
                right=p21*b+p22*d
                outputs.append((left+right)%p)
                if entries == (0,p-1,1,0) and outputs[-1]!=b:
                    raise RuntimeError('J read failure')
            if (len(set(outputs))==p) != (p21!=0):
                raise RuntimeError('Probe injectivity failure')
            checks+=1
# Independently solve exact two-letter guards in each residue mod32.
guards={}
for e,f in ((1,3),(3,1)):
    valid=[]
    for r in range(1,32,2):
        first=3*r+1
        if first%(2**e) or first%(2**(e+1))==0: continue
        nxt=first//2**e
        second=3*nxt+1
        if second%(2**f) or second%(2**(f+1))==0: continue
        valid.append(r)
    guards[str((e,f))]=valid
if guards != {'(1, 3)':[19], '(3, 1)':[29]}:
    raise RuntimeError('Guards mismatch')
print(json.dumps({'status':'PASS','finite_field_probe_checks':checks,'guards':guards},indent=2,sort_keys=True))
