# Part A output provenance — refreshed with Mac results

The current smoke and RDF data were generated on the user's M1 Mac, macOS
27.0.1, on 2026-10-08 and uploaded as A3_PartA 2.zip. The simulation C source
matches the earlier implementation. No user simulation data were regenerated
by the assistant during this update. Notebook plots now use the Mac data.

Mac-generated files:
- a_smoke_thermo.csv: 1000-step DPD run, 3000 monomers, box 10^3, all a=25.
- a_smoke_rdf.csv / a_smoke_rdf_blocks.csv: 30 frames from that run.
- a_smoke.pdb: five frames, A labelled C and B labelled O for display.
- a_smoke.restart: final native-format M1 Mac restart.
- a_smoke_run.log / a_smoke_stderr.log: run output and empty stderr. Shell timing
  was shown in Terminal, not saved to that stderr file.
- a_uniform_rdf.csv / a_uniform_rdf_blocks.csv: 100 independent uniform
  configurations, 500 beads, A:B=150:350, box 10^3. Normalisation check, not B2 MD.
- a_smallbox_rdf.csv / a_smallbox_rdf_blocks.csv: one-frame box-6 fallback check.

Transcripts made by the assistant from the user's messages:
- a_mac_checks_transcript.log: pasted output of sh tests/run_checks.sh; all passed.
- a_mac_run_transcript.txt: exact commands, 3.841-second timing, final thermo rows,
  and limited OVITO loading/playback observations.

Historical Linux checks (not Mac-generated logs):
- a_checks.log / a_checks_stderr.log: original focused checks and timing.
- a_sanitizer.log / a_sanitizer_stderr.log: original address/undefined-behaviour
  checks; leak detection disabled due to execution-environment restrictions.

Other:
- a_package_verification.txt: refreshed local packaging/notebook verification;
  this is not the official course checker.

No B-D production runs are included. See the notebook Build & run log.
