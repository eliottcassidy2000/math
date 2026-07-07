import numpy as np
GRID=300000; xs=(np.arange(GRID)+0.5)/GRID
def mu_17(E):
    ph=np.mod(np.outer(xs,np.array(E,float)),1.0); ph.sort(axis=1)
    allg=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1]).reshape(-1,1)],axis=1)
    return float((allg.max(axis=1)>1/7).mean())

ap13 = mu_17(list(range(13)))
print(f"=== k=13 STRUCTURED adversary test vs AP-minimality (opus-S131) ===")
print(f"mu_17(AP {{0..12}}) = {ap13:.5f} = 477/1078 = {477/1078:.5f}\n")
tests = {
  "AP {0..12}": list(range(13)),
  "12-AP + outlier@20": list(range(12))+[20],
  "12-AP + outlier@100": list(range(12))+[100],
  "12-AP + outlier@1000": list(range(12))+[1000],
  "11-AP + 2 out @{50,51}": list(range(11))+[50,51],
  "two 6-blocks {0..5}+{50..56}"[:40]: list(range(6))+list(range(50,57)),
  "two 6-blocks gap100": list(range(7))+list(range(100,106)),
  "AP step2 {0,2,..,24}": list(range(0,25,2)),
  "AP with 1 internal gap {0..5,7..13}": list(range(6))+list(range(7,14)),
  "AP with hole {0..12}\{6}+{13}": [i for i in range(13) if i!=6]+[13],
  "near-AP jitter": [0,1,2,3,4,5,6,7,8,9,10,11,13],
  "near-AP jitter2": [0,1,2,3,4,5,7,8,9,10,11,12,14],
  "3 blocks of 4-ish": list(range(4))+list(range(20,24))+list(range(40,45)),
  "geometric": [1,2,3,4,6,9,13,19,28,42,63,94,141],
}
below=[]
for name,E in tests.items():
    if len(set(E))!=13: 
        print(f"  [skip {name}: not 13 distinct]"); continue
    v=mu_17(E); flag = "  <-- BELOW AP!" if v<ap13-0.0005 else ("  = AP" if abs(v-ap13)<0.0005 else "")
    if v<ap13-0.0005: below.append((name,v))
    print(f"  {name:<34} mu_17={v:.5f}{flag}")
print(f"\nVERDICT: {'AP-minimality HOLDS (no structured adversary below AP).' if not below else 'BROKEN: '+str(below)}")
