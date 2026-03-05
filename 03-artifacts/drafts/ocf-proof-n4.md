# OCF Proof at n=4 (Direct Algebra)

## Setup

Tournament T on V = {i, j, a, b} with arc i->j. Let T' = flip to j->i.
Others = {a, b}. Define s_a = 1 - T[a][i] - T[j][a], s_b = 1 - T[b][i] - T[j][b].
gamma = T[a][b].

## Goal

Show: adj(i,j) - adj'(j,i) = -2(s_a + s_b)

This equals -2*sum(s_x) since H(B_x) = 1 for all x (B_x has 1 vertex) and D5-C5=0.

## Computation

Decompose by position a of i in the path (i at position a, j at a+1):

**Position (1,2): S=empty, R={a,b}**
D(0) = h_start({j,a,b}, j) - h_start({i,a,b}, i)
     = [T[j][a]*T[a][b] + T[j][b]*T[b][a]] - [T[i][a]*T[a][b] + T[i][b]*T[b][a]]
     = (T[j][a]-T[i][a])*gamma + (T[j][b]-T[i][b])*(1-gamma)

Using T[j][a] = 1-T[a][j] = 1-(1-T[j][a]) and T[i][a] = 1-T[a][i]:
  T[j][a]-T[i][a] = (1-T[a][j])-(1-T[a][i]) = T[a][i]-T[a][j]

Actually: alpha = T[a][i], beta = T[j][a]. Then s_a = 1-alpha-beta.
  T[j][a]-T[i][a] = beta-(1-alpha) = alpha+beta-1 = -s_a.

Similarly: delta_v = T[b][i], epsilon = T[j][b]. s_b = 1-delta_v-epsilon.
  T[j][b]-T[i][b] = epsilon-(1-delta_v) = delta_v+epsilon-1 = -s_b.

So D(0) = -s_a*gamma - s_b*(1-gamma).

**Position (2,3): |S|=1, sum over S={a} and S={b}**
D(1) = [T[a][i]*T[j][b] + T[b][i]*T[j][a]] - [T[a][j]*T[i][b] + T[b][j]*T[i][a]]
     = alpha*epsilon + delta_v*beta - (1-beta)*(1-delta_v) - (1-epsilon)*(1-alpha)
     = alpha*epsilon + delta_v*beta - 1+beta+delta_v-beta*delta_v - 1+epsilon+alpha-alpha*epsilon
     = beta + delta_v + alpha + epsilon - 2
     = (alpha+beta-1) + (delta_v+epsilon-1)
     = -s_a - s_b

**Position (3,4): S={a,b}, R=empty**
D(2) = h_end({a,b,i}, i) - h_end({a,b,j}, j)
     = [T[a][b]*T[b][i] + T[b][a]*T[a][i]] - [T[a][b]*T[b][j] + T[b][a]*T[a][j]]
     = gamma*(delta_v-(1-epsilon)) + (1-gamma)*(alpha-(1-beta))
     = gamma*(-s_b) + (1-gamma)*(-s_a)
     = -s_b*gamma - s_a*(1-gamma)

**Total:**
delta = D(0) + D(1) + D(2)
      = [-s_a*gamma - s_b*(1-gamma)] + [-s_a - s_b] + [-s_b*gamma - s_a*(1-gamma)]
      = -s_a*gamma - s_b + s_b*gamma - s_a - s_b - s_b*gamma - s_a + s_a*gamma
      = -2*s_a - 2*s_b
      = -2*(s_a + s_b)  QED

## Key observations

1. D(1) = -(s_a+s_b) directly, independent of gamma = T[a][b].
2. D(0) + D(2) = -s_a*(gamma + 1-gamma) - s_b*(1-gamma + gamma) = -(s_a+s_b).
   But this split is: D(0) = -s_a*gamma - s_b*(1-gamma), D(2) = -s_b*gamma - s_a*(1-gamma).
3. The gamma-dependent terms in D(0) and D(2) cancel in the total.
4. The proof generalizes: at any n, the "middle position" term D((n-2)/2) may directly
   give the s_x terms, while the "boundary" terms D(0)+D(n-2) telescope.

## Connection to OCF

Since -2(s_a+s_b) = 2(D3-C3) and I(Omega,2) = 1 + 2*alpha_1,
DeltaI = 2*Delta(alpha_1) = 2(D3-C3) = -2*sum(s_x) = adj(i,j)-adj'(j,i) = DeltaH.
Combined with base case (transitive tournament: H=1, Omega=empty), this proves OCF at n=4.
