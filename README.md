# Coherence Field

An at-a-glance particle-field display for coherence in randomness — the way the
GCP Dot compresses network variance into a single color, spatialized.

The project has two generations. It began as a visualizer for
[coherence_monitor](https://codeberg.org/TaoishTechy/coherence_monitor); a
statistical review showed that monitor's Health metric is structurally constant
(a step function of a matrix rank), so the project grew its own measurement
layer: **a GCP-style random-event-generator network built from the Mac's
physical noise sources, with calibrated statistics and a live control channel.**
That instrument, in `reg/`, is the current heart of the project.

## The REG instrument

Three physical *eggs* — microphone ADC noise, the built-in accelerometer read
at its native ~800 Hz through IOKit, and the camera sensor's LSB plane behind a
translucent diffuser — each whitened to match its measured physics and formed
into GCP-convention 200-bit trials every second, beside a `/dev/random` CSPRNG
control run through byte-identical code. Every egg is scored against its own
empirically calibrated mean and variance (estimated only from history older
than the live window), the network is summarized by Stouffer Z and network
variance, and every claim of "unusual" is made against a bootstrap null from
the control egg. Full design, methodology, and honest-expectations notes:
[reg/README.md](reg/README.md).

**Chance, honestly displayed** — the live instrument on real data. Each ring is
an egg breathing with its current z; the grey outer ring is the CSPRNG control
ghost; the outer annulus is the GCP cumulative-deviation plot bent into a
circle with a bootstrap 95 % simultaneous band. Here the status chip is
quantifying a mild ten-minute fluctuation (Stouffer Z +2.03, 95.7th percentile
of the control's own history) — the kind chance produces routinely:

![The live REG field: four labelled rings of particles, cumulative traces
inside the translucent band, and an ATTENTION chip quantifying a mild
fluctuation against the control](docs/field-live.png)

**What a real signal would look like** — simulation mode with the "correlated
eggs" injection (a shared component across the physical eggs, the signature
network variance is built to detect). The families converge toward gold, the
cumulative trace erupts past the dashed 99 % ring, and the control ghost stays
home; the chip states the numbers that triggered EXCURSION:

![Simulated excursion: the gold cumulative trace far outside the band while
the control trace stays inside; netvar z +8.77](docs/field-excursion.png)

Deep links tour the states:
[`?sim=corr`](https://petersgrandadventure.github.io/coherence-visualizer/reg/field.html?sim=corr) ·
[`?sim=mean`](https://petersgrandadventure.github.io/coherence-visualizer/reg/field.html?sim=mean) ·
[`?sim=var`](https://petersgrandadventure.github.io/coherence-visualizer/reg/field.html?sim=var) ·
[`?sim=fault`](https://petersgrandadventure.github.io/coherence-visualizer/reg/field.html?sim=fault)
(the hosted page runs simulation mode; live mode needs the local bridge).

It also carries the experiment: **pre-registered windows** (declare start,
duration, statistic and hypothesised direction before the window opens; the
bridge evaluates it automatically afterwards and the ledger has no file
drawer) and **the walk** — a PEAR-style intention ball whose full lap is the
session's 1 % line — plus a **zen view** (press Z) that leaves nothing but the
rings. Details in [reg/README.md](reg/README.md).

### Running it

From Terminal.app on a Mac (camera access is granted per launching app):

```bash
cd reg && ./run_calibration.sh 0        # acquisition daemon, runs until stopped
../monitor/venv/bin/python reg_bridge.py &   # statistics bridge + field at http://localhost:5006
```

`./report.sh` generates the full calibration report (constants, whiteness,
stationarity, cross-egg structure, covariate coupling, Monte-Carlo bands) at
any time. Requirements: numpy, opencv-python, pyaudio in a venv (see
`monitor/`'s setup, or any venv — the daemon takes `--no-camera` etc. to run
with whatever hardware is available).

## v1 — visualizer for coherence_monitor

The original single-file visualizer for the upstream monitor's ledger packets
is still here and still fun:
**▶ [Live demo](https://petersgrandadventure.github.io/coherence-visualizer/)**
(simulation mode; try the **Coherence event** and **Pazuzu surge** buttons, or
jump straight in:
[`?demo=coherence`](https://petersgrandadventure.github.io/coherence-visualizer/coherence_field.html?demo=coherence) ·
[`?demo=pazuzu`](https://petersgrandadventure.github.io/coherence-visualizer/coherence_field.html?demo=pazuzu)).

![v1: a high-coherence event — the field unified in gold, particles migrating
from the continuum cloud to the boundary ring](docs/coherence-event.png)

Its grammar maps the monitor's schema directly: Health drives hue convergence
(scattered cool blues → unified gold), σ drives jitter, ρ drives rotation,
CI_B spikes migrate particles to the boundary ring, and a Marchenko–Pastur
eigenvalue spike locks the field into a seven-pointed star. To run it against
a real coherence_monitor ledger: `python3 coherence_bridge.py` beside
`coherence_ledger.db`, then open <http://localhost:5005> (5005 because macOS
AirPlay squats on 5000). The **History** button replays logged packets with a
scrubbable timeline. Note the review in the project history before treating
the upstream metrics as physically meaningful.

## License

MIT — see [LICENSE](LICENSE).
