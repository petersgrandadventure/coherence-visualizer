# Coherence Field

A particle-field visualizer for [coherence_monitor](https://codeberg.org/TaoishTechy/coherence_monitor) —
the state of the field is readable at a glance, the way the GCP Dot compresses
network variance into a single color.

**▶ [Live demo](https://petersgrandadventure.github.io/coherence-visualizer/)** —
runs in simulation mode; try the **Coherence event** and **Pazuzu surge** buttons.
(For live monitor data, run the bridge below and open `http://localhost:5005`
instead — or paste your bridge's address into the demo's endpoint field.)

Deep links jump straight to a state:
[`?demo=coherence`](https://petersgrandadventure.github.io/coherence-visualizer/coherence_field.html?demo=coherence) ·
[`?demo=pazuzu`](https://petersgrandadventure.github.io/coherence-visualizer/coherence_field.html?demo=pazuzu)

![A high-coherence event: Health 0.96, the field unified in gold, particles
migrating from the continuum cloud to the boundary ring](docs/coherence-event.png)

## Quick start

Open `coherence_field.html` in any browser. It starts in **simulation mode** with a
realistic packet generator; use the **Coherence event** / **Pazuzu surge** buttons to
tour the states, and **How to read** for the full legend.

## Live data

On the machine running the monitor:

```bash
python3 coherence_bridge.py            # from the folder containing coherence_ledger.db
```

Then open <http://localhost:5005> — the page loads already connected, polling
`/api/latest` every 2 s (the monitor logs a packet every 4.2 s). Alternatively open
the HTML file directly and press **Connect live**.

- Default port is **5005**, not 5000 — macOS AirPlay Receiver occupies 5000.
- To watch from another machine: `python3 coherence_bridge.py --host 0.0.0.0`,
  then point the endpoint field at `http://<monitor-ip>:5005/api/latest`.
- `--db /path/to/coherence_ledger.db` if the bridge doesn't run beside the database.

The bridge is stdlib-only (no Flask required) and opens the database read-only, so
it is safe to run alongside the monitor.

## History playback

The **History** button replays logged data through the field. With a bridge
connected it pulls the last 1200 packets (~84 minutes) from `/api/recent?n=1200`;
without one it replays whatever this page has seen during the current session.
The timeline strip charts Health across the loaded span — gold dots are
high-coherence events (Health > 0.9), red ticks are Pazuzu flags — so you can
scrub directly to the moments that matter. Playback runs at 1× (real time,
one packet per 4.2 s), 8×, or 30×.

## Visual grammar

| Signal | Encoding |
|---|---|
| Health / CI_C | Hue coherence: scattered cool blues → one unified gold band as Health passes 0.9; motion aligns from Brownian wander into shared orbital flow |
| σ (noise) | Per-particle jitter amplitude |
| ρ (spectral radius) | Rotation speed of the coherent flow |
| CI_B spikes | Particles migrate from the central continuum cloud to the outer boundary ring (Axiom H₁₃ conservation, made visible) |
| λmax/λMP > 1.15 | Particles lock onto a seven-pointed star — the heptagonal residue (Gate 7 / Marchenko–Pastur spike) |
| λmax/λMP > 1.30 | Crimson turbulence — the Pazuzu abort state |
| PSI < 0.3 | The structure disintegrates outward |
| Gates 1–7 | Seven particle families; each family's brightness follows its gate's normalized ledger index (G1–G7 tiles in the HUD) |

HUD: threshold ticks on each metric bar match the monitor's safety table
(σ 0.053, ρ 0.95, r/dₛ 0.93, λ 1.15/1.30), plus a 3-minute Health sparkline and the
latest payload hash. Press **H** to hide the panel for a pure field view.

## License

MIT — see [LICENSE](LICENSE).
