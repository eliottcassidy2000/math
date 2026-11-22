"""Exercise the split inverse and drive the remaining one-flank faces."""
import argparse
from random import Random
import subprocess
from fermat_receiver import A, D, EDGE, TARGET, flipped, join, p, receiver, tagged_data


def invert(program, slope, h, k0, k1, target, bits, shard=None):
    args = [program, slope, h, k0, k1, target, bits]
    if shard is not None:
        args.extend(shard)
    output = subprocess.check_output(list(map(str, args)), text=True)
    hits = [int(line.split()[1]) for line in output.splitlines() if line.startswith("address ")]
    assert all(slope*y+h*p(y)+(k1 if y%2 else k0) == target for y in hits)
    return sorted(hits), output


def controls(program):
    rng = Random(61005)
    for s, h, k0, k1 in ((11, 3, 5, 7), (21, 9, 9, 3), (31097, 2465, 235, 789)):
        for i in range(16):
            y = rng.randrange(20_000)
            target = s*y + h*p(y) + (k1 if y%2 else k0) + i%2
            wanted = [z for z in range(target//s+1)
                      if s*z+h*p(z)+(k1 if z%2 else k0) == target]
            result, _ = invert(program, s, h, k0, k1, target, 5)
            assert result == wanted, (s, h, target, result, wanted)
    print("All planted and adjacent full-width controls passed.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("program")
    ap.add_argument("--spine", type=int)
    ap.add_argument("--bits", type=int, default=14)
    ap.add_argument("--side", choices=("left", "right", "both"), default="both")
    ap.add_argument("--fixed", type=int, default=0)
    ap.add_argument("--through-fixed", type=int, help="cover fixed values from --fixed up to this bound")
    ap.add_argument("--shard", nargs=2, type=int)
    ar = ap.parse_args()
    if ar.spine is None:
        controls(ar.program)
    else:
        h,q,c,t,e,r,k = join(join(tagged_data(A,12),tagged_data(flipped(D,*EDGE),35)),receiver(ar.spine))
        for side in ("left", "right"):
            if ar.side not in (side, "both"):
                continue
            slope,odd,base_slope,base_odd = (t,r,c,e) if side=="left" else (c,e,t,r)
            stop = ar.fixed if ar.through_fixed is None else ar.through_fixed
            if not 0 <= ar.fixed <= stop:
                ap.error("invalid fixed-address interval")
            totals = [0,0,0]
            for fixed in range(ar.fixed,stop+1):
                base = k+base_slope*fixed+base_odd*(fixed%2)+h*p(fixed)
                hits, report = invert(ar.program,slope+q*fixed,h,base,base+odd,TARGET,ar.bits,ar.shard)
                for line in report.splitlines():
                    parts=line.split()
                    if parts[:1] == ["parity"]:
                        totals[0]+=int(parts[parts.index("blocks")+1])
                        totals[1]+=int(parts[parts.index("candidates")+1])
                totals[2]+=len(hits)
                if hits or stop==ar.fixed:
                    print("spine",ar.spine,"varying",side,"fixed",fixed,flush=True)
                    print(report,end="",flush=True)
            print("spine",ar.spine,"varying",side,"fixed range",ar.fixed,stop,
                  "blocks/candidates/hits",*totals,flush=True)
