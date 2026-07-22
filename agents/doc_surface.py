"""Single source of truth for bounded agent-facing documentation and routing."""

FRONTIER_EPOCH = "2026-07-22"

ALWAYS_READ_ROUTES = (
    "00-navigation/START-HERE.md",
    "00-navigation/CURRENT-FRONTIER.md",
    "01-canon/ACTIVE-GUARDRAILS.md",
)

WORKFLOW_ROUTES = (
    "00-navigation/RESEARCH-PROTOCOL.md",
    "00-navigation/META-PATTERNS.md",
)

ON_DEMAND_ROUTES = (
    "05-knowledge/reference/CORE-PAPERS.md",
    "00-navigation/PROBLEM-LEDGER.md",
)

CANONICAL_ROUTES = ALWAYS_READ_ROUTES + WORKFLOW_ROUTES + ON_DEMAND_ROUTES

STARTUP_DOCS = (
    "AGENTS.md",
    "CLAUDE.md",
    "README.md",
    *ALWAYS_READ_ROUTES,
    *WORKFLOW_ROUTES,
    "05-knowledge/reference/CORE-PAPERS.md",
    "00-navigation/CONCURRENT-SESSIONS.md",
    "01-canon/README.md",
    "05-knowledge/README.md",
    "07-reflections/README.md",
    "agents/README.md",
)

TOPIC_SEARCH_ROUTES = (
    "00-navigation/CURRENT-FRONTIER.md",
    "01-canon/ACTIVE-GUARDRAILS.md",
    "05-knowledge/reference/CORE-PAPERS.md",
    "00-navigation/PROBLEM-LEDGER.md",
    "00-navigation/LRC14-PROOF-MAP.md",
    "00-navigation/META-PATTERNS.md",
)

# Only the maintained prefix of these historical ledgers is eligible for a
# "current route" match. The suffix remains searchable as explicitly labelled
# provenance through ordinary rg or the historical file group.
SEARCHABLE_PREFIX_END = {
    "00-navigation/PROBLEM-LEDGER.md": "## Legacy portfolio snapshot",
    "00-navigation/LRC14-PROOF-MAP.md": "> ## 2026-07-19",
}

TOPIC_SEARCH_TREES = (
    "01-canon/theorems",
    "05-knowledge/hypotheses",
    "07-reflections",
)

TOPIC_ARTIFACT_TREES = (
    "04-computation",
    "05-knowledge/results",
    "00-navigation/TANGENTS.md",
)

LINE_BUDGETS = {
    "AGENTS.md": 140,
    "CLAUDE.md": 70,
    "README.md": 120,
    "00-navigation/START-HERE.md": 155,
    "00-navigation/CURRENT-FRONTIER.md": 450,
    "01-canon/ACTIVE-GUARDRAILS.md": 190,
    "00-navigation/RESEARCH-PROTOCOL.md": 300,
    "00-navigation/META-PATTERNS.md": 400,
    "00-navigation/SESSION-LOG.md": 120,
    "05-knowledge/hypotheses/INDEX.md": 120,
    "05-knowledge/reference/CORE-PAPERS.md": 625,
    "00-navigation/CONCURRENT-SESSIONS.md": 160,
    "01-canon/README.md": 140,
    "05-knowledge/README.md": 120,
    "07-reflections/README.md": 100,
    "agents/README.md": 140,
}

PREFIX_LINE_BUDGETS = {
    "00-navigation/PROBLEM-LEDGER.md": 40,
    "00-navigation/LRC14-PROOF-MAP.md": 150,
}

# A runaway generated ledger must never silently become part of the bounded
# startup packet, even when its line count happens to remain small.
STARTUP_BYTE_BUDGET = 180_000
MAX_STARTUP_LINE_BYTES = 500
MAX_EMITTED_MATCH_BYTES = 360
MAX_STARTUP_PACKET_BYTES = 30_000
