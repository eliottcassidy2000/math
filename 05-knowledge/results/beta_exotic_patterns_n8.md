# Exotic Betti profiles at n=8

Discovered by kind-pasteur-S45 (2026-03-09).

## Profiles observed (from ~2000 random n=8 tournaments)

| Profile | chi | Freq | Description |
|---------|-----|------|-------------|
| [1,0,0,0,0,0,0,0] | 1 | ~83% | Contractible |
| [1,0,0,1,0,0,0,0] | 0 | ~14% | H_3 = Z |
| [1,1,0,0,0,0,0,0] | 0 | ~1% | H_1 = Z |
| [1,0,0,0,1,0,0,0] | 2 | ~0.8% | H_4 = Z (EVEN!) |
| [1,0,0,1,1,0,0,0] | 1 | ~0.15% | H_3 + H_4 (coexist!) |
| [1,1,0,0,0,1,0,0] | -1 | ~0.05% | H_1 + H_5 (coexist!) |
| [1,0,0,0,5,0,0,0] | 6 | rare | H_4 = Z^5 |

## Key observations

1. **beta_2 = 0 ALWAYS** (0 violations in 2000+ samples)
2. **beta_4 is the first even Betti to appear** (onset at n=8)
3. **beta_3 and beta_4 CAN coexist** (chi=1 when they do)
4. **beta_1 and beta_5 CAN coexist** (chi=-1, extremely rare)
5. **beta_1 and beta_3 NEVER coexist** (THM-095 seesaw, 0 violations in 2000+)
6. All exotic cases are strongly connected
7. c3 range for beta_4>0: 13-20 (not restricted to near-regular)

## Coexistence rules (empirical)

| Pair | Coexist? | Note |
|------|----------|------|
| beta_1 + beta_3 | NEVER | THM-095 (seesaw via beta_2=0) |
| beta_3 + beta_4 | YES | chi = 0 - 1 + 1 = 0 or 1 |
| beta_4 + beta_5 | Not observed | May occur at higher n |
| beta_1 + beta_5 | YES | chi = -1 |
| beta_1 + beta_4 | Not observed in 2000 | May be impossible? |

## Score sequences for exotic cases

- beta_1+beta_5: near-regular (3,3,3,3,4,4,4,4), c3=20
- beta_4 only: varied scores, c3=13-20
- beta_3+beta_4: c3=14-17, not near-regular
