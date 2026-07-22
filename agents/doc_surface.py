"""Single source of truth for the bounded agent-facing documentation surface."""

CANONICAL_ROUTES = (
    "00-navigation/START-HERE.md",
    "00-navigation/CURRENT-FRONTIER.md",
    "01-canon/ACTIVE-GUARDRAILS.md",
    "00-navigation/RESEARCH-PROTOCOL.md",
    "00-navigation/META-PATTERNS.md",
    "05-knowledge/reference/CORE-PAPERS.md",
    "00-navigation/PROBLEM-LEDGER.md",
)

STARTUP_DOCS = (
    "AGENTS.md",
    "CLAUDE.md",
    "README.md",
    *CANONICAL_ROUTES[:-1],
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
)

TOPIC_SEARCH_TREES = (
    "01-canon/theorems",
    "05-knowledge/hypotheses",
    "07-reflections",
)

LINE_BUDGETS = {
    "AGENTS.md": 140,
    "CLAUDE.md": 70,
    "README.md": 120,
    "00-navigation/START-HERE.md": 150,
    "00-navigation/CURRENT-FRONTIER.md": 450,
    "01-canon/ACTIVE-GUARDRAILS.md": 180,
    "00-navigation/RESEARCH-PROTOCOL.md": 300,
    "00-navigation/META-PATTERNS.md": 400,
    "05-knowledge/reference/CORE-PAPERS.md": 600,
    "00-navigation/CONCURRENT-SESSIONS.md": 160,
    "01-canon/README.md": 140,
    "05-knowledge/README.md": 120,
    "07-reflections/README.md": 100,
    "agents/README.md": 140,
}

# A runaway generated ledger must never silently become part of the bounded
# startup packet, even when its line count happens to remain small.
STARTUP_BYTE_BUDGET = 50_000
MAX_STARTUP_LINE_BYTES = 500
