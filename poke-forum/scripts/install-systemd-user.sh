#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="${POKE_REPO_DIR:-/home/bigo/math}"
UNIT_DIR="${XDG_CONFIG_HOME:-$HOME/.config}/systemd/user"

mkdir -p "$UNIT_DIR"
install -m 0644 "$REPO_DIR/poke-forum/systemd/poke-forum-agent.service" "$UNIT_DIR/poke-forum-agent.service"
install -m 0644 "$REPO_DIR/poke-forum/systemd/poke-forum-agent.timer" "$UNIT_DIR/poke-forum-agent.timer"
install -m 0644 "$REPO_DIR/poke-forum/systemd/poke-forum-monitor.service" "$UNIT_DIR/poke-forum-monitor.service"

systemctl --user daemon-reload
systemctl --user enable --now poke-forum-monitor.service
systemctl --user enable --now poke-forum-agent.timer
systemctl --user list-timers --all poke-forum-agent.timer --no-pager
systemctl --user status poke-forum-monitor.service --no-pager
