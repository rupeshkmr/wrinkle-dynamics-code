#!/usr/bin/env bash
# Full pipeline: simulate a config with checkpoints, then upsample the result.
#
#   scripts/run_pipeline.sh configs/<id>.json [upsample.py options, e.g. --every 20]
#
# Needs a build with both executables, e.g.
#   cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release -DARGUS_BUILD_UPSAMPLER=ON
#   cmake --build build/release -j
# Set BUILD_DIR to pick a build directory (default: build/, else the newest build/*/ with both binaries).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ $# -lt 1 ]]; then
    sed -n '2,9p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi
CONFIG="$1"; shift

if [[ -z "${BUILD_DIR:-}" ]]; then
    for d in "$ROOT/build" $(ls -td "$ROOT"/build/*/ 2>/dev/null); do
        if [[ -x "$d/bin/wrinkle_dynamics" && -x "$d/bin/upsample" ]]; then BUILD_DIR="$d"; break; fi
    done
fi
if [[ -z "${BUILD_DIR:-}" || ! -x "$BUILD_DIR/bin/wrinkle_dynamics" || ! -x "$BUILD_DIR/bin/upsample" ]]; then
    echo "No build with bin/wrinkle_dynamics and bin/upsample found; build with -DARGUS_BUILD_UPSAMPLER=ON or set BUILD_DIR." >&2
    exit 1
fi

# The simulator resolves the config path relative to the repository root.
CONFIG_ABS="$(cd "$(dirname "$CONFIG")" && pwd)/$(basename "$CONFIG")"
CONFIG_REL="${CONFIG_ABS#"$ROOT"/}"
if [[ "$CONFIG_REL" == "$CONFIG_ABS" ]]; then
    echo "Config must be inside the repository ($ROOT)." >&2
    exit 1
fi

mkdir -p "$ROOT/output"
LOG="$ROOT/output/simulate_$(basename "$CONFIG" .json)_$(date +%Y%m%d_%H%M%S).log"
echo "== simulate $CONFIG_REL (log: $LOG)"
"$BUILD_DIR/bin/wrinkle_dynamics" "$CONFIG_REL" --nogui --checkpoint 2>&1 | tee "$LOG"

CKPT="$(sed -n 's/^Initiated checkpointing to //p' "$LOG" | tail -1)"
if grep -q "no checkpoints are requested" "$LOG"; then
    echo "This build of wrinkle_dynamics cannot save checkpoints; reconfigure with -DARGUS_CHECKPOINT=ON (implied by -DARGUS_BUILD_UPSAMPLER=ON) and rebuild." >&2
    exit 1
fi
if [[ -z "$CKPT" || ! -d "$CKPT" ]]; then
    echo "Could not find the checkpoint directory in the simulation log." >&2
    exit 1
fi
CKPT="$(cd "$CKPT" && pwd)"

echo "== upsample $CKPT"
python3 "$ROOT/scripts/upsample.py" "$CKPT" --bin "$BUILD_DIR/bin/upsample" "$@"
echo "== done: checkpoint $CKPT (upsampled meshes are in the render/ folder of the run directory listed above)"
