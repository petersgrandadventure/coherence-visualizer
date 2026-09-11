# I built a random-event instrument out of my laptop. It found nothing, correctly.

This started as a toy. An upstream project called coherence_monitor emits a stream of numbers — Health, sigma, rho, something called CI_B — and I wanted to *see* them rather than read them. So, working with Claude (Anthropic's Claude Code, my collaborator across this whole build), I made a particle field: Health pulled the hues from scattered cool blues toward gold, sigma jittered the particles, rho spun them, CI_B spikes pushed them out to a boundary ring, and a Marchenko–Pastur eigenvalue spike locked the whole thing into a seven-pointed star. Press the "Coherence event" button and the screen goes gold. It was lovely.

![The v1 visualizer in simulation mode: particles unified in gold with a seven-pointed star. Pretty — but it depicts a metric the review found not physically meaningful.](https://raw.githubusercontent.com/petersgrandadventure/coherence-visualizer/main/docs/coherence-event.png)

Then I put the monitor itself under a statistical review, and the review was uncomfortable. Its Health metric is structurally constant: a step function of a matrix rank, so it sits at one value until the rank jumps. It doesn't measure the world. I had built a beautiful display for a number that meant nothing.

I could have stopped. Instead I asked the question underneath the toy. The Global Consciousness Project, which the visualizer was quietly gesturing at, runs a network of hardware random-event generators — "eggs" — and looks for departures from chance. Could a laptop be an egg? Could it be three?

## Three eggs and a ghost

Yes, with caveats that took three weeks to learn.

The microphone was the obvious first egg: ADC thermal noise plus room acoustics, captured at the device's native 44.1 kHz. Bit planes 0–3 (measured fair to bit 6), each XOR-folded 8:1, give about 22 kbit/s.

The surprise was the accelerometer. I hadn't known my MacBook Pro had one reachable from userland, but IOKit exposes it at about 795 Hz in 2⁻¹⁶ g steps, with roughly 50 LSB rms of noise riding on the signal. That noise is the product: LSB parity per axis, XOR-folded 4:1, gives about 600 bit/s. Tiny beside the mic, but bit rate is irrelevant beyond a few hundred bits per second — the GCP convention is one 200-bit trial per second per device. What matters is serial correlation below about 10⁻³, stationarity over hours, and independence between eggs. (CPU timing jitter auditioned as a fourth egg and was cut: it is a deterministic function of the 24 MHz clock phase and machine load.)

Beside the physical eggs sits a ghost: `/dev/random`, the kernel's Fortuna CSPRNG, run through byte-identical code. Zero physical content by construction. Every claim of "unusual" the instrument makes is made against a bootstrap null drawn from that control.

## The tape saga

The camera egg taught me the most, mostly by refusing to work.

The idea is sound: sensor noise, read off the least-significant bit of every pixel. But the raw LSB plane is a mess after the image processor has had its way — lag-1 temporal correlation of 0.65, and 49% of pixels showing a fixed pattern where the ideal is 5.8%. XOR-folding each 16×16 block into one bit fixes it (bias 0.5013, temporal r 0.005, fixed-pattern fraction 5.0%); 8×8 was marginal. And the whole sensor counts as one egg, because quadrants share the processor's gain and exposure and correlate up to −0.28.

Then I covered the lens, so the sensor would see a stable flat field rather than my room (no frames are ever stored — only trial sums and statistics). The first tape leaked light, and the sensor quietly tracked the daylight through it. So I switched to properly opaque tape, which was worse: below about 8 gray the processor black-clamps the sensor and the LSB plane simply freezes. The code now labels that regime, a little wearily, *invalid (too dim: sensor near black clamp — add light or a lighter diffuser)*.

What works is a translucent diffuser. A covered camera doesn't read black at all: auto-exposure maxes the gain and the sensor sees a flat mid-gray, and that flat-field signature is how the instrument recognizes the covered regime and calibrates it separately from uncovered. The dim-light lesson came from the calibration data as the room dimmed in the evenings: between 8 and 30 gray the camera isn't clamped, but it isn't clean either — its ten-minute Stouffer Z had a standard deviation of 1.1 to 3.6 instead of 1.0, from minute-scale drift. Measured over seven days, it is trustworthy only from about 30 gray up, so the instrument never scores the camera below a gray mean of 30, and never falls back to a pooled constant.

## Calibration, and a week of nothing

Each egg is scored against its own empirical mean and variance, estimated only from history older than the live window. Theory says 100 and √50 = 7.071; the calibrated constants land within a hair — mic 99.9985 / 7.0701 on 64.7 million trials, accelerometer 99.9978 / 7.0722 on 1.69 million, control 99.9988 / 7.0712 on 18.7 million. The daemon also logs, per second and per egg, bias, lag-1..10 serial correlation, a runs test and mask-phase bias on the pre-mask bits, each egg's health (mic level, camera gray, accelerometer rms, device identity), and host covariates (load, clock drift, loop lateness) — so any excursion can be explained or stand on its own.

Then I let it run. Over 20.6 wall-clock days it recorded 7.33 days (the daemon runs under `caffeinate`, which prevents idle sleep but not a closed lid). The network — summarized by the Stouffer Z across the physical eggs and by the network-variance statistic — was indistinguishable from the ghost. Rolling ten-minute Stouffer |Z| exceeded the control's own 99th percentile 0.57% of the time (100 episodes, 14.5 a day) against the control's 0.34% (94 episodes, 13.6 a day); on the one-hour variance statistic the order flips (27 episodes against the control's 38); and once the camera's dim-light seconds are excluded, the ten-minute rate lands at 0.35% against 0.34%. Largest single-second |Z| all week: 4.98, where the expected maximum for that many seconds is about 5.2; the control's was 4.95. Cross-egg per-second correlations all sit within |r| ≤ 0.0022.

The biggest ten-minute variance excursion of the week peaked at netvar z +4.54, one Monday night, lasting nine minutes. Another coincided with the mic at −12.88 dBFS and the accelerometer at 1,369 LSB rms — the room was loud and the machine was being jostled. That is what the covariate log is for.

![The live REG field on real data: four labelled rings, cumulative traces inside the bootstrap band, and an ATTENTION chip quantifying a mild ten-minute fluctuation (Stouffer Z +2.03, 95.7th percentile of the control's own history) — the kind chance produces routinely.](https://raw.githubusercontent.com/petersgrandadventure/coherence-visualizer/main/docs/field-live.png)

## What you see on screen

Four labelled rings of particles. Three breathe with their egg's current z — mic, accelerometer, camera, each with its own hue. The fourth is a grey ghost, the control, and it never changes character. As the ten-minute network statistic climbs, the physical families' drift aligns and their hues converge toward gold.

Around the outside runs the GCP cumulative-deviation plot bent into a circle, with the bootstrap 95% band as a translucent annulus, the 99% limits dashed, and the control's own cumulative trace ghosted alongside. The classic pointwise ±1.96√n parabola is exited about half the time under the null — a display that draws it cries wolf every other window — so the band here is *simultaneous*: chosen so that 95% of whole random paths stay inside it at every step. The thing to watch for is the real trace leaving the band while the ghost stays home. A status chip reads CHANCE-LIKE, ATTENTION or EXCURSION with the numbers that triggered it, and crimson INSTRUMENT if an egg goes stale or a constant is uncalibrated — never for the physics.

Press Z and everything but the rings, the ball and the outer trace disappears.

## The walks

A null instrument is only interesting if you can ask it a question, so the field carries a pre-registration ledger. You declare start, duration, one statistic and a hypothesized direction at least 60 seconds before the window opens; the record is hashed, appended, never edited; the bridge evaluates it automatically after the window closes, and every registration gets a result row, including "no data." No file drawer.

The friendliest way in is the walk, descended from PEAR's random-walk intention displays: a small ball on a thin ring. Each second's network Stouffer Z is one step; the ball's angle is the cumulative sum since the session opened. HI means clockwise; LO mirrors the ring so clockwise is still toward your intention; BL is PEAR's baseline — intend nothing, and the ball should stay near the top. Small ticks mark the 5% line for your session length, and one full lap is the 1% line: precisely "this session ended as a p < 0.01 result." Gold shows only while the ball is currently past the 5% line, never latched, and the only position that counts is where it sits when the window closes.

A session goes like this. You choose HI, LO or BL and a duration, click Start walk, and click again to confirm — the registration is permanent and cannot be cancelled, because a window you could withdraw after seeing the data wouldn't be a pre-registration. A sixty-second settle follows; the window opens; a gold arc appears on the outer ring with a countdown. You sit. The ball wanders. Thirty to ninety seconds after the window closes, the ring label switches to the ledger's numbers — Z, sidedness, p versus the control, n — and the result row behind it records every egg's value, the control over the same window, the band check and the covariates.

My ledger holds four windows so far: two smoke tests, one ten-minute HI session that Claude ran (Z −0.18, p 0.59 against the control), and my own five-minute HI walk — Stouffer Z +0.24, p 0.42; mic +0.54, accelerometer +0.42, camera −0.55; the control egg read −1.78 in the same window. The formal series combines to Z +0.105 (two-sided p 0.92), and the control's matched windows combine to −2.13. Four windows is far too few for either number to mean anything — that a null channel with zero physical content sits at two sigma after four windows is exactly what a handful of windows looks like, and exactly why the ghost is there.

## Honest expectations

The GCP's published effect is about 0.3σ per pre-registered event, pooled over roughly 60 devices. A three-egg home instrument cannot detect an effect that size, and this one has not detected anything. It does not measure consciousness, intention or psi. What it can do is show, with a defined null and logged covariates, whether *your* field ever departs from chance — and what the room was doing when it did. A striking walk is not evidence; the evaluation is. And an instrument that showed something in its first week would more likely be reporting a bug than a discovery.

![For contrast, what a real signal would look like — simulation mode with the correlated-eggs injection: the gold cumulative trace far outside the band while the control ghost stays inside, netvar z +8.77. Nothing in this image is data.](https://raw.githubusercontent.com/petersgrandadventure/coherence-visualizer/main/docs/field-excursion.png)

## Try it

No Mac handy? The field runs in simulation mode with nothing installed — [`field.html?sim=corr`](https://petersgrandadventure.github.io/coherence-visualizer/reg/field.html?sim=corr) (also `?sim=mean`, `?sim=var`, `?sim=fault`) — and the original v1 visualizer is at [petersgrandadventure.github.io/coherence-visualizer](https://petersgrandadventure.github.io/coherence-visualizer/). Both are simulations, not data.

Live mode needs a Mac (tested on an Apple Silicon MacBook Pro; the accelerometer comes through IOKit, and the daemon takes `--no-camera` and similar flags to run with whatever hardware you have), a Python venv at `monitor/venv` with numpy, opencv-python and pyaudio (the scripts call `../monitor/venv/bin/python`), something translucent over the camera, and ports 5005 and 5006 free. Then, from **Terminal.app** — macOS grants camera access per launching application:

1. `cd reg && ./run_calibration.sh 0` — the acquisition daemon, until `./run_calibration.sh stop`. It runs under `caffeinate -i`, which prevents idle sleep but not a closed lid, so leave the lid open.
2. `../monitor/venv/bin/python reg_bridge.py &` — the read-only statistics bridge, serving the field at `http://localhost:5006`.
3. Let it calibrate. An egg needs 50,000 trials before empirical constants replace the theoretical ones; the extra trials each second get you there in hours. `./report.sh` writes the calibration report at any time.
4. Open the field, choose an intention and a duration, click Start walk, confirm, and sit for a while. The bridge evaluates the window after it closes (or run `./reg_register.py evaluate`); `./reg_register.py list` shows the ledger.
5. `./reg_register.py summary` prints your formal series.

## Compare ledgers with me

My formal series is four windows long. Yours will start empty. Run a few sessions and some baselines, then post the output of `summary` in an issue on the repo, with the control's matched result beside it. I am curious whether anyone's home field ever leaves the band while its ghost stays inside — and at least as curious about the ledgers that never do. Both answers are the point.

Code, MIT: [github.com/petersgrandadventure/coherence-visualizer](https://github.com/petersgrandadventure/coherence-visualizer).

---

## Social versions

**Short (one post):**
I turned my MacBook into a GCP-style random-event instrument — mic noise, accelerometer, camera sensor — with a CSPRNG control and a pre-registration ledger. A week of data: indistinguishable from chance, as it should be. Story + code: https://github.com/petersgrandadventure/coherence-visualizer

**Thread:**
1. Started as a particle visualizer for a "coherence monitor." Looked great. Then a statistical review found its Health metric is a step function of a matrix rank — structurally constant. A beautiful display for a number that meant nothing. So: pivot. 1/7
2. The pivot: build a GCP-style random-event-generator network from the Mac's own physics. Mic ADC noise. Camera sensor noise. And, a surprise to me, the built-in accelerometer, readable via IOKit at ~795 Hz — its ~50 LSB rms of noise is the product. 2/7
3. Beside the three physical eggs: /dev/random, the kernel CSPRNG, run through byte-identical code. Zero physical content by construction. Every claim of "unusual" is made against a bootstrap null drawn from that control. 3/7
4. The tape saga: the first tape leaked light and the sensor tracked the daylight; opaque tape drove it into black clamp and froze the LSB plane. What works: a translucent diffuser. And below 30 gray the camera drifts, so it is never scored dim. 4/7
5. 7.33 days recorded. Rolling 10-min |Z| above the control's own 99th percentile: network 0.57% of the time, control 0.34%; other cuts flip the order. Max single-second |Z| 4.98 vs control 4.95. Cross-egg |r| ≤ 0.0022. Indistinguishable from chance — what a null instrument should show. 5/7
6. Then the walks: PEAR-style intention sessions that are pre-registered windows — hashed, never edited, evaluated automatically; one full lap of the ball is the session's 1% line. My 5-min HI walk: Z +0.24, p 0.42. Four windows combine to Z +0.105, p 0.92. 6/7
7. Honest expectations: the GCP's published effect is ~0.3σ per event pooled over ~60 devices. Three eggs at home can't see that, and this one hasn't seen anything. Built with Claude Code as collaborator. MIT, come compare ledgers: https://github.com/petersgrandadventure/coherence-visualizer 7/7

**LinkedIn:**
I turned my MacBook into a GCP/PEAR-style random-event instrument and spent three weeks learning what it takes for "nothing" to mean something.

Three physical noise sources — the microphone's ADC noise, the built-in accelerometer (readable via IOKit at ~795 Hz, which surprised me), and the camera sensor behind a translucent diffuser — each whitened to its measured physics and formed into GCP-convention 200-bit trials every second, beside a /dev/random control run through byte-identical code. Every egg is calibrated against its own history; every claim of "unusual" is judged against a bootstrap null from the control; every experiment is a pre-registered window in a ledger with no file drawer.

A week of recorded data: indistinguishable from the control in the frequency, duration and magnitude of excursions. Four pre-registered windows so far combine to Z +0.105 (two-sided p 0.92). The GCP's published effect is ~0.3σ per event over ~60 devices — far below what three eggs at home can detect — and I claim no detection of anything. What the instrument can do is show whether your field ever leaves chance, and what the room was doing when it did.

Built in collaboration with Claude Code. MIT-licensed; I'd love to compare ledgers: https://github.com/petersgrandadventure/coherence-visualizer
