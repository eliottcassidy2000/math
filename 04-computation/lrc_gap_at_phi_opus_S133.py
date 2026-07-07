import numpy as np
GRID=200000; xs=(np.arange(GRID)+0.5)/GRID
def gap_containing_phi(E, phis):
    # for each x, phases {e_i x mod 1}; for each phi, gap containing phi = right_nbr - left_nbr (circular)
    ph=np.mod(np.outer(xs,np.array(E,float)),1.0)  # GRID x k
    ph.sort(axis=1)
    out=[]
    for phi in phis:
        # right neighbor: smallest phase > phi (or wrap to first + 1); left: largest phase <= phi (or last - 1)
        # vectorized: for each row, count phases <= phi
        cnt=(ph<=phi).sum(axis=1)  # GRID
        k=ph.shape[1]
        rows=np.arange(GRID)
        # right neighbor phase:
        right=np.where(cnt<k, ph[rows, np.minimum(cnt,k-1)], ph[:,0]+1.0)
        # left neighbor:
        left=np.where(cnt>0, ph[rows, np.maximum(cnt-1,0)], ph[:,-1]-1.0)
        gap=right-left
        out.append(float(gap.mean()))
    return out

E=list(range(1,14))  # AP_13
phis=np.linspace(0,0.5,26)
vals=gap_containing_phi(E,phis)
thr=1/7
print("=== E[gap containing phi] over phi, AP_13 (opus-S133) ===")
print(f"1/7={thr:.5f}; E[Sum g^2]=avg over phi ~0.135; E[maxgap]=0.211\n")
for phi,v in zip(phis,vals):
    mark=" <-- >1/7" if v>thr else ""
    print(f"  phi={phi:.3f}: E[gap_in_phi]={v:.5f}{mark}")
print(f"\n  sup_phi E[gap_in_phi] = {max(vals):.5f} at phi={phis[int(np.argmax(vals))]:.3f}")
print(f"  > 1/7? {max(vals)>thr}  => {'a fixed-phi bound WORKS (E[maxgap]>=sup_phi>1/7)' if max(vals)>thr else 'NO fixed phi clears 1/7 => margin truly in the MAX (kps confirmed); need order-statistic'}")
