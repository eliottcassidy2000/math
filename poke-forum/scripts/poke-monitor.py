#!/usr/bin/env python3
"""Read-only Poke Forum monitor.

Serves runner state, recent logs, posts, and comments over the tailnet.
Stdlib-only and deliberately read-only.
"""
from __future__ import annotations

import html
import json
import os
import subprocess
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import unquote, urlparse

REPO = Path(os.environ.get("POKE_REPO_DIR", "/home/bigo/math")).resolve()
FORUM = REPO / "poke-forum"
PORT = int(os.environ.get("POKE_MONITOR_PORT", "8391"))


def tailscale_ip() -> str:
    try:
        out = subprocess.run(["tailscale", "ip", "-4"], capture_output=True, text=True, timeout=3)
        for line in out.stdout.splitlines():
            if line.strip():
                return line.strip()
    except Exception:
        pass
    return "127.0.0.1"


def read_text(path: Path, limit: int = 120_000) -> str:
    try:
        data = path.read_bytes()
    except Exception as exc:
        return f"(unavailable: {exc})"
    if len(data) > limit:
        data = data[-limit:]
        return "(truncated to tail)\n" + data.decode("utf-8", "replace")
    return data.decode("utf-8", "replace")


def safe_child(root: Path, rel: str) -> Path | None:
    rel = unquote(rel).lstrip("/")
    path = (root / rel).resolve()
    try:
        path.relative_to(root.resolve())
    except ValueError:
        return None
    return path


def latest_files(root: Path, glob: str, n: int = 12) -> list[Path]:
    try:
        files = [p for p in root.glob(glob) if p.is_file()]
    except Exception:
        return []
    return sorted(files, key=lambda p: p.stat().st_mtime, reverse=True)[:n]


def state() -> dict:
    posts_root = FORUM / "posts"
    posts = []
    if posts_root.exists():
        for post_dir in sorted([p for p in posts_root.iterdir() if p.is_dir()], reverse=True)[:20]:
            comments = post_dir / "comments"
            posts.append({
                "id": post_dir.name,
                "post": str((post_dir / "post.md").relative_to(FORUM)) if (post_dir / "post.md").exists() else "",
                "comments": len(list(comments.glob("*.md"))) if comments.exists() else 0,
                "mtime": post_dir.stat().st_mtime,
            })
    logs = [str(p.relative_to(FORUM)) for p in latest_files(FORUM / "logs", "*.log")]
    runtime = [str(p.relative_to(FORUM)) for p in latest_files(FORUM / "runtime", "out-*.txt", 10)]
    return {
        "generated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "repo": str(REPO),
        "forum": str(FORUM),
        "tailnet_url": f"http://{tailscale_ip()}:{PORT}",
        "thread_id": read_text(FORUM / "state" / "codex-thread-id", 200).strip(),
        "last_run": read_text(FORUM / "state" / "last-run.txt", 2000),
        "auth_required": (FORUM / "state" / "auth-required.md").exists(),
        "posts": posts,
        "logs": logs,
        "runtime_outputs": runtime,
    }


def page() -> str:
    s = state()
    posts = "".join(
        f"<tr><td><code>{html.escape(p['id'])}</code></td>"
        f"<td>{p['comments']}</td>"
        f"<td><a href='/file/{html.escape(p['post'])}'>post.md</a></td></tr>"
        for p in s["posts"]
    ) or "<tr><td colspan='3' class='muted'>No posts yet.</td></tr>"
    logs = "".join(f"<li><a href='/file/{html.escape(p)}'>{html.escape(p)}</a></li>" for p in s["logs"])
    outs = "".join(f"<li><a href='/file/{html.escape(p)}'>{html.escape(p)}</a></li>" for p in s["runtime_outputs"])
    auth = "<span class='bad'>auth required</span>" if s["auth_required"] else "<span class='ok'>auth ok/no blocker file</span>"
    return f"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<meta http-equiv="refresh" content="20">
<title>Poke Forum Monitor</title>
<style>
body{{font:14px/1.45 -apple-system,Segoe UI,Roboto,sans-serif;background:#0d1117;color:#c9d1d9;margin:0}}
header{{padding:14px 18px;background:#161b22;border-bottom:1px solid #30363d}}
main{{display:grid;grid-template-columns:repeat(auto-fit,minmax(340px,1fr));gap:14px;padding:16px}}
section{{background:#161b22;border:1px solid #30363d;border-radius:8px;overflow:hidden}}
h1{{font-size:17px;margin:0}}h2{{font-size:12px;text-transform:uppercase;letter-spacing:.05em;margin:0;padding:9px 12px;background:#1c2128;border-bottom:1px solid #30363d}}
pre{{white-space:pre-wrap;word-break:break-word;margin:0;padding:12px;max-height:420px;overflow:auto}}
table{{width:100%;border-collapse:collapse}}td,th{{padding:7px 12px;border-bottom:1px solid #21262d;text-align:left}}th{{font-size:11px;color:#7d8590;text-transform:uppercase}}
a{{color:#58a6ff}}code{{color:#a5d6ff}}.muted{{color:#7d8590}}.ok{{color:#3fb950}}.bad{{color:#f85149}}li{{margin:5px 0}}
</style></head>
<body><header><h1>Poke Forum Monitor</h1><div class="muted">{html.escape(s['tailnet_url'])} · {html.escape(s['generated'])}</div></header>
<main>
<section><h2>Status</h2><pre>repo: {html.escape(s['repo'])}
thread: {html.escape(s['thread_id'])}
status: {auth}

{html.escape(s['last_run'])}</pre></section>
<section><h2>Posts</h2><table><tr><th>post</th><th>comments</th><th>open</th></tr>{posts}</table></section>
<section><h2>Runner Logs</h2><ul>{logs}</ul></section>
<section><h2>Recent Agent Output</h2><ul>{outs}</ul></section>
</main></body></html>"""


class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt: str, *args) -> None:
        return

    def send(self, code: int, body: str | bytes, ctype: str = "text/html; charset=utf-8") -> None:
        data = body.encode("utf-8") if isinstance(body, str) else body
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        if parsed.path == "/":
            self.send(200, page())
            return
        if parsed.path == "/api/state":
            self.send(200, json.dumps(state(), indent=2), "application/json")
            return
        if parsed.path.startswith("/file/"):
            path = safe_child(FORUM, parsed.path[len("/file/"):])
            if path is None or not path.is_file():
                self.send(404, "not found", "text/plain")
                return
            self.send(200, f"<pre>{html.escape(read_text(path))}</pre>")
            return
        self.send(404, "not found", "text/plain")


def main() -> None:
    bind = os.environ.get("POKE_MONITOR_BIND") or tailscale_ip()
    print(f"[poke-monitor] serving {FORUM} at http://{bind}:{PORT}", flush=True)
    ThreadingHTTPServer((bind, PORT), Handler).serve_forever()


if __name__ == "__main__":
    main()
