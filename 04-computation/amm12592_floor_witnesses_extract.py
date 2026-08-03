"""Extract and exactly verify the THM-3029 gamma* floor-profile witnesses.

Produces 04-computation/amm12592_floor_witnesses_R8_R16_R32.json:
  R = 8, 16 : direct beam solutions of THM-3002 (*) at the gamma* floor
              profile d_i = floor(gamma*(R+i)), D0 = 0;
  R = 32    : the (gamma = 1/2, D0 = 3) beam solution LIFTED onto the gamma*
              floor profile via THM-3026 (L) (convolution with the constant-1
              block [binom(d'-d,k)]_k) -- THM-3029 (M).

Every witness is verified EXACTLY (integer arithmetic only) before saving:
  (a) each block admissible at its profile degree
      (|delta_k| <= binom(d_i,k), delta_k == binom(d_i,k) mod 2), and
  (b) the epoch identity of THM-3002 (*):
      sum_{i<R} x^i sum_k delta_{i,k} x^{d_i-k}(1-x)^k == (1-x)^{R-1} in Z[x].

gamma* is entered as the 16-digit truncation GS = 5979874356654402/10^16
(gamma* = log(phi)/log(sqrt 5) = 0.59798743566544014975...); at R <= 64 the
floor profile of GS coincides with that of gamma* itself (THM-3029 R2).
"""
import io, json, os, sys
from contextlib import redirect_stdout

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import amm12592_gamma35_beam_deathstar as beam
with redirect_stdout(io.StringIO()):           # suppress liftrate's import-time demo
    from liftrate import prof, lift_block, admissible, epoch_identity, eff

GS = (5979874356654402, 10**16)

def check(R, blocks, profile):
    """Exact verification: admissibility of every block + epoch identity."""
    assert len(blocks) == R == len(profile)
    adm = all(admissible(blocks[i], profile[i]) for i in range(R))
    idt = epoch_identity(R, blocks, profile)
    return adm and idt

witnesses = []

# --- R = 8, 16: direct solves at the gamma* floor profile, D0 = 0 ----------
for R in (8, 16):
    p = prof(R, *GS, 0)
    sol, msg = beam.solve(R, g1=GS[0], g2=GS[1], D0=0, beam=800, ctrl=2, span=2)
    assert sol is not None, f"R={R}: {msg}"
    assert beam.verify(R, sol, g1=GS[0], g2=GS[1], D0=0)
    ok = check(R, sol, p)
    assert ok, f"R={R}: exact re-verification failed"
    witnesses.append({"R": R, "profile": p, "blocks": sol, "verified": True})
    print(f"R={R:2d} direct: {msg}; profile={p}")
    print(f"      blocks admissible + epoch identity exact: {ok}; eff rate = {eff(R,p)} = {float(eff(R,p)):.6f}")

# --- R = 32: lift the (gamma=1/2, D0=3) solution onto the floor profile ----
R = 32
src = prof(R, 1, 2, 3)
sol, msg = beam.solve(R, g1=1, g2=2, D0=3, beam=250, ctrl=2, span=2)
assert sol is not None, f"R=32 source: {msg}"
assert beam.verify(R, sol, g1=1, g2=2, D0=3)
tgt = prof(R, *GS, 0)
assert all(tgt[i] >= src[i] for i in range(R)), "target not pointwise >= source"
lifted = [lift_block(sol[i], src[i], tgt[i]) for i in range(R)]
ok = check(R, lifted, tgt)
assert ok, "R=32: exact re-verification of lifted witness failed"
witnesses.append({"R": R, "profile": tgt, "blocks": lifted, "verified": True})
print(f"R=32 lifted from (gamma=1/2, D0=3) source ({msg}); profile={tgt}")
print(f"      blocks admissible + epoch identity exact: {ok}; eff rate = {eff(R,tgt)} = {float(eff(R,tgt)):.6f}")

out = os.path.join(HERE, "amm12592_floor_witnesses_R8_R16_R32.json")
with open(out, "w") as f:
    json.dump(witnesses, f, indent=1)
print(f"\nwrote {out}")
print("ALL WITNESSES VERIFIED EXACTLY (admissibility + THM-3002 epoch identity).")
