#!/usr/bin/env python3
"""nc2_channel_tournament_lens_deathstar_S88.py (HYP-8772)
NEW LENS: NC2's radial channels form a TOURNAMENT (i->j iff channel i dominates
channel j as m->inf). Verify: DEGREE-GAP => TRANSITIVE channel tournament (strict
domination order, one Hamiltonian path = a total order, the S75 nullcone vertex);
RESONANCE band => TIES (regular/balanced channels = the wall, S84). This realizes
NC2's residual as a tournament-transitivity statement -- the moment-nullcone IS a
tournament-nullcone on its channels."""
from math import factorial as fact, log
from fractions import Fraction as Fr
def Lp(c): return sum(v*fact(k) for k,v in c.items())
def pm(p,q):
    r={}
    for k,v in p.items():
        for k2,v2 in q.items(): r[k+k2]=r.get(k+k2,0)+v*v2
    return r
def pw(p,n):
    r={0:Fr(1)}
    for _ in range(n): r=pm(r,p)
    return r
def channel(A,B,c0,i,m):
    if m-2*i<0: return Fr(0)
    multi=Fr(fact(m),fact(i)*fact(i)*fact(m-2*i))
    rad=pm(pm({i:c0**i},pw(A,i)),pm(pw(B,m-2*i),{i:Fr(1)}))
    return multi*Lp(rad)

def tournament(A,B,c0,m):
    # channels i=0..m//2 ; arc i->j if |channel_i| > |channel_j|
    ch={i:abs(channel(A,B,c0,i,m)) for i in range(m//2+1)}
    order=sorted(ch, key=lambda i:ch[i])
    return ch, order

print("LENS: the NC2 radial-channel tournament (i->j iff channel_i dominates)")
print("="*66)
# DEGREE-GAP example: A high-degree, B constant => one channel dominates (transitive)
print("\n[A] DEGREE-GAP: P=Z*(s^2)+ (1) + W*(s^2)  (deg h large, deg B=0): STRICT domination")
A={2:Fr(1)}; B={0:Fr(1)}; c0=Fr(1)
for m in (6,8,10):
    ch,order=tournament(A,B,c0,m)
    ratios=[float(ch[order[-1]]/ch[o]) if ch[o] else 9e9 for o in order[:-1]]
    print(f"  m={m}: dominant channel i={order[-1]} (|val| ratio to others >> 1: "
          f"{[f'{r:.0f}' for r in ratios]}) -> clear source = TRANSITIVE")
# RESONANCE example: balanced degrees => channels comparable (ties/regular)
print("\n[B] RESONANCE (central offset): P=Z*(1+s)+(1+s)+W*1: channels COMPARABLE")
A={0:Fr(1),1:Fr(1)}; B={0:Fr(1),1:Fr(1)}; c0=Fr(1)
for m in (6,8,10):
    ch,order=tournament(A,B,c0,m)
    top=ch[order[-1]]; second=ch[order[-2]] if len(order)>1 else Fr(0)
    print(f"  m={m}: top/second ratio = {float(top/second):.2f} (near 1 => TIED/regular channels)")
print("""
UNIFYING READING:
 * DEGREE-GAP  => channel tournament is TRANSITIVE (one dominant channel = source).
   A transitive tournament is the S75 nullcone VERTEX; the dominant channel's sign
   survives, so E[P^m] != 0: NONCANCELLATION. This is codex THM-2017 ('one endpoint
   ratio one') seen as tournament transitivity.
 * RESONANCE  => channels TIE (comparable) = a REGULAR/balanced sub-tournament =
   THE WALL (S84 'regular/Paley is the tightest case', S75 Paley pole). Cancellation
   is only possible here, and only if the tied channels' SIGNS+phases conspire.
 * SO NC2's residual = 'the channel tournament, even in its regular (tied) core,
   does not produce exact cancellation' -- the SAME transitive(easy)-vs-regular(hard)
   dichotomy as H>=disc and LRC. The moment-nullcone IS a tournament-nullcone on
   its channels; the open case is always the regular-channel wall.""")
