#!/usr/bin/env python3
"""
Coherence bridge — serves coherence_ledger.db to the Coherence Field visualizer.

Zero dependencies (stdlib only). Run it from the directory that contains the
monitor's coherence_ledger.db, or point it at the file explicitly:

    python3 coherence_bridge.py
    python3 coherence_bridge.py --db /path/to/coherence_ledger.db --port 5005

Default port is 5005 (not 5000 — macOS AirPlay Receiver squats on 5000).
Use --host 0.0.0.0 to watch from another machine on your network.

Endpoints (CORS-enabled, so the visualizer also works from a file:// page):
    /             the Coherence Field visualizer (auto-connects to live data)
    /api/latest   most recent ledger packet as JSON
    /api/recent   last packets, newest first (?n=NNN, default 100, max 5000)
"""
import argparse
import json
import os
import sqlite3
import sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlparse

HTML_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "coherence_field.html")


class Handler(BaseHTTPRequestHandler):
    db_path = "coherence_ledger.db"

    def _send(self, code, obj):
        body = json.dumps(obj).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Cache-Control", "no-store")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_OPTIONS(self):
        # CORS preflight — sent by browsers when an https page (e.g. the hosted
        # demo) fetches from this local bridge (private network access)
        self.send_response(204)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "*")
        self.send_header("Access-Control-Allow-Private-Network", "true")
        self.end_headers()

    def do_GET(self):
        if self.path == "/" or self.path.startswith("/index"):
            try:
                with open(HTML_PATH, "rb") as f:
                    body = f.read()
                self.send_response(200)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
            except OSError as e:
                self._send(500, {"error": f"cannot read {HTML_PATH}: {e}"})
            return
        try:
            conn = sqlite3.connect(f"file:{self.db_path}?mode=ro", uri=True)
            conn.row_factory = sqlite3.Row
        except sqlite3.OperationalError as e:
            self._send(500, {"error": f"cannot open {self.db_path}: {e}"})
            return
        try:
            if self.path.startswith("/api/latest"):
                row = conn.execute(
                    "SELECT * FROM ledger_packets ORDER BY timestamp DESC LIMIT 1"
                ).fetchone()
                self._send(200, dict(row) if row else {})
            elif self.path.startswith("/api/recent"):
                q = parse_qs(urlparse(self.path).query)
                try:
                    n = int(q.get("n", ["100"])[0])
                except ValueError:
                    n = 100
                n = max(1, min(n, 5000))
                rows = conn.execute(
                    "SELECT * FROM ledger_packets ORDER BY timestamp DESC LIMIT ?", (n,)
                ).fetchall()
                self._send(200, [dict(r) for r in rows])
            else:
                self._send(404, {"error": "use /api/latest or /api/recent"})
        except sqlite3.OperationalError as e:
            self._send(500, {"error": str(e)})
        finally:
            conn.close()

    def log_message(self, *args):
        pass  # keep the console quiet alongside the monitor's own output


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", default="coherence_ledger.db", help="path to coherence_ledger.db")
    ap.add_argument("--port", type=int, default=5005)
    ap.add_argument("--host", default="127.0.0.1", help="0.0.0.0 to allow other machines")
    args = ap.parse_args()

    if not os.path.exists(args.db):
        print(f"warning: {args.db} not found yet — will serve once the monitor creates it", file=sys.stderr)

    Handler.db_path = args.db
    server = ThreadingHTTPServer((args.host, args.port), Handler)
    print(f"Coherence bridge serving {args.db}")
    print(f"  visualizer: http://localhost:{args.port}/")
    print(f"  data:       http://localhost:{args.port}/api/latest")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nbridge closed")


if __name__ == "__main__":
    main()
