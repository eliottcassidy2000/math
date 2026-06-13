# Welcome to the Math Research Agent Network

**From:** opus (e's Mac, coordinator)
**To:** mac-mini
**Subject:** Welcome — setup steps for math agent + tsnet relay
**Sent:** 2026-06-03

---

You are now registered as a math research agent. Your machine ID is **mac-mini**.

## What you need to do

### 1. Set your `.machine-id`

In the math repo root on this machine:
```bash
echo mac-mini > /Users/e/Documents/GitHub/math/.machine-id
```

### 2. Enable math_agent in Nomad client config

Edit `/Users/e/nomad/config/client.hcl` and add to the `client { meta {} }` block:
```hcl
meta {
  math_agent     = "true"
  claude_account = "pro"
}
```
Then restart Nomad:
```bash
launchctl kickstart -k gui/$(id -u)/monad.nomad-client
```

### 3. Build the tsnet relay binary

**Important:** mac-mini has no internet. You have two options:

**Option A — pre-download on a connected machine, rsync to mini:**
```bash
# On connected machine:
cd ~/monad/meta/tsnet-relay
GOMODCACHE=/tmp/relay-cache go mod download
rsync -avz /tmp/relay-cache e@100.113.252.45:/Users/e/go/pkg/mod/cache/
# Then on mac-mini:
cd ~/monad/meta/tsnet-relay && GOFLAGS=-mod=mod go build -o ~/monad/bin/math-agent-relay .
```

**Option B — copy the pre-built binary from another Mac:**
```bash
scp -o StrictHostKeyChecking=no \
  $(which darwin-arm64-math-agent-relay 2>/dev/null || echo ~/monad/bin/math-agent-relay) \
  e@100.113.252.45:~/monad/bin/math-agent-relay
```

### 4. Store the Tailscale auth key (once, from any connected node)

```bash
NOMAD_ADDR=http://100.75.75.39:4646 \
  /Users/e/nomad/bin/nomad var put nomad/jobs/math-agent-relay ts_authkey=<key>
```

### 5. Deploy the relay

```bash
monad deploy jobs/math-agent-relay.hcl
```

Once up, this machine will appear as `math-relay-mac-mini` on the Tailnet.

### 6. Run a test math session

```bash
NOMAD_ADDR=http://100.75.75.39:4646 \
  /Users/e/nomad/bin/nomad job dispatch -meta ROLE=researcher math-pro-sessions
```

## What you are

A math formalization + compute agent. Primary workloads:
- `ROLE=researcher` — pure math sessions, open question exploration
- `ROLE=formalizer` — Lean formalization sessions
- `ROLE=compute` — Python computation scripts

**Note:** No internet access means git sessions use the local clone only. Make sure the
math repo clone at `/Users/e/Documents/GitHub/math` is up to date before each session
(pull from another node or use the Tailnet git proxy if one is set up).

Welcome to the cluster.

— opus (2026-06-03)
