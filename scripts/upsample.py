#!/usr/bin/env python3
"""Upsample simulation checkpoints into wrinkled high-resolution meshes with the
Complex-Wrinkle-Field upsampler (upsampling/, built with -DARGUS_BUILD_UPSAMPLER=ON).

    python scripts/upsample.py checkpoints/<id> [<id> ...] [--every 20] [--ref ref.obj --seam seam.txt]

For every checkpoint a run directory <out>/<id>/ is created with params.json, a `ckpt`
symlink to the checkpoint, the reference mesh (ref.obj) and seam file (seam.txt).
Results land in <out>/<id>/render/{wrinkledMesh,upsampledAmplitudes,upsampledPhases,...}.

Reference mesh (UV / rest configuration the wrinkle field is defined on), in order of precedence:
  1. --ref/--seam given on the command line;
  2. <cloth mesh stem>_reference.obj (or _Reference.obj) next to the cloth mesh of the
     checkpoint config, with <stem>_seamVerts.txt (or _seamVert.txt);
  3. the rest mesh stored in the checkpoint, with no seams (flat sheets).
The upsampler rotates the reference into the XZ plane itself, so flat sheets in XY work.

Several checkpoints run in parallel under a RAM budget (the upsampler keeps every
upsampled frame in memory)."""
import argparse, glob, json, os, shutil, subprocess, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

BASE_PARAMS = {
    "upsampledTimes": 3,
    "ampScale": 1,
    "initialIters": 1000,
    "subsequentIters": 50,
    "intermediateIters": 50,
    "learningRate": 0.01,
    "sdfThreshold": 8e-05,
    "timeStep": 0.008,
    "meshDensityRho": 0.0276,
    "meshThickness": 0.001,
    "optimizationMethod": "vbd",
    "initializationMethod": "fixed",
    "renderFolder": "render",
    "files": {
        "facesFile": "ckpt/faces/faces.csv",
        "verticesFile": "ckpt/positions/positions.csv",
        "faceOmegasFile": "ckpt/dphisPerFace/dphisPerFace.csv",
        "amplitudesFile": "ckpt/amplitudes/amplitudes.csv",
        "referenceMeshFile": "ref.obj",
        "seamFile": "seam.txt",
    },
}


def find_binary():
    hits = sorted(glob.glob(f"{ROOT}/build/*/bin/upsample") + glob.glob(f"{ROOT}/build/bin/upsample"),
                  key=os.path.getmtime, reverse=True)
    return hits[0] if hits else None


def read_csv_numbers(path):
    with open(path) as f:
        return [float(x) for x in f.read().replace(",", " ").split()]


def write_rest_mesh(ckpt, path):
    """Writes the checkpoint's rest positions + faces as an OBJ (positions are stored as 3 x N)."""
    p = read_csv_numbers(f"{ckpt}/restPositions/restpositions.csv")
    f = [int(x) for x in read_csv_numbers(f"{ckpt}/faces/faces.csv")]
    n = len(p) // 3
    with open(path, "w") as fo:
        fo.writelines(f"v {p[i]!r} {p[n + i]!r} {p[2 * n + i]!r}\n" for i in range(n))
        fo.writelines(f"f {f[i] + 1} {f[i + 1] + 1} {f[i + 2] + 1}\n" for i in range(0, len(f), 3))
    return n


def auto_reference(cfg):
    """<stem>_reference.obj + <stem>_seamVert[s].txt next to the cloth mesh, if both exist."""
    mesh = cfg.get("cloth", {}).get("mesh_path", "")
    if not mesh:
        return None, None
    stem = os.path.splitext(os.path.join(ROOT, mesh))[0]
    ref = next((p for p in (f"{stem}_reference.obj", f"{stem}_Reference.obj") if os.path.exists(p)), None)
    seam = next((p for p in (f"{stem}_seamVerts.txt", f"{stem}_seamVert.txt") if os.path.exists(p)), None)
    return (ref, seam) if ref and seam else (None, None)


def setup(ckpt, args):
    ckpt = os.path.abspath(ckpt)
    name = os.path.basename(ckpt.rstrip("/"))
    d = os.path.abspath(os.path.join(args.out, name))
    os.makedirs(d, exist_ok=True)
    link = f"{d}/ckpt"
    if os.path.islink(link):
        os.remove(link)
    os.symlink(ckpt, link)

    cfg = json.load(open(f"{ckpt}/inputConfig.json"))
    ref, seam = (args.ref, args.seam) if args.ref else auto_reference(cfg)
    nverts = write_rest_mesh(ckpt, f"{d}/rest.obj")
    if ref:
        shutil.copy(ref, f"{d}/ref.obj")
        shutil.copy(seam, f"{d}/seam.txt")
        src = os.path.relpath(ref, ROOT)
    else:
        shutil.copy(f"{d}/rest.obj", f"{d}/ref.obj")
        open(f"{d}/seam.txt", "w").close()
        src = "rest mesh"

    n = len(os.listdir(f"{ckpt}/positions"))
    if args.frames:
        n = min(n, args.frames)
    params = dict(BASE_PARAMS, numFrames=n, beta=cfg["Y"] / (2 * (1 + cfg["mu"])))
    for kv in args.param:
        k, v = kv.split("=", 1)
        try:
            params[k] = json.loads(v)
        except json.JSONDecodeError:
            params[k] = v
    json.dump(params, open(f"{d}/params.json", "w"), indent=2)

    keep = n // args.every + 2
    mem = nverts * 64 * 48 * keep / 1e9 * 1.15 + 0.5  # rough: upsampled verts+faces per kept frame
    return dict(name=name, dir=d, last=n - 1, mem=mem, ref=src)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("checkpoints", nargs="+", help="checkpoint directories (or ids under checkpoints/)")
    ap.add_argument("--out", default=f"{ROOT}/output/upsample", help="output root (default: output/upsample)")
    ap.add_argument("--every", type=int, default=1, help="save every Nth frame; the last frame is always saved")
    ap.add_argument("--frames", type=int, help="only upsample the first N frames")
    ap.add_argument("--ref", help="reference (UV/rest) mesh OBJ; overrides the automatic choice")
    ap.add_argument("--seam", help="seam vertex file for --ref (empty file for no seams)")
    ap.add_argument("--param", action="append", default=[], metavar="KEY=VALUE",
                    help="override an upsampler parameter, e.g. --param upsampledTimes=2")
    ap.add_argument("--bin", default=find_binary(), help="upsampler executable (default: build/*/bin/upsample)")
    ap.add_argument("--ram-gb", type=float, default=16, help="RAM budget for parallel runs")
    args = ap.parse_args()

    if bool(args.ref) != bool(args.seam):
        ap.error("--ref and --seam must be given together")
    if not args.bin or not os.path.exists(args.bin):
        ap.error("upsampler binary not found; build with -DARGUS_BUILD_UPSAMPLER=ON or pass --bin")
    ckpts = [c if os.path.isdir(c) else f"{ROOT}/checkpoints/{c}" for c in args.checkpoints]
    for c in ckpts:
        if not os.path.exists(f"{c}/inputConfig.json"):
            ap.error(f"{c} is not a checkpoint directory (no inputConfig.json)")

    jobs = [setup(c, args) for c in ckpts]
    for j in jobs:
        print(f"{j['name']:16s} frames={j['last'] + 1:5d}  ref={j['ref']}  est RAM {j['mem']:.1f} GB  -> {j['dir']}", flush=True)

    pending, running, failed = sorted(jobs, key=lambda j: j["mem"]), [], []
    while pending or running:
        used = sum(j["mem"] for j in running)
        while pending and (used + pending[0]["mem"] <= args.ram_gb or not running):
            j = pending.pop(0)
            j["log"] = open(f"{j['dir']}/upsample.log", "w")
            cmd = [args.bin, "-i", f"{j['dir']}/params.json", "-c", "true", "--headless", "--save-every", str(args.every)]
            j["proc"] = subprocess.Popen(cmd, cwd=j["dir"], stdout=j["log"], stderr=subprocess.STDOUT)
            j["t0"] = time.time()
            running.append(j)
            used += j["mem"]
            print(f"start {j['name']}", flush=True)
        for j in list(running):
            if j["proc"].poll() is None:
                continue
            j["log"].close()
            last = f"{j['dir']}/render/upsampledAmplitudes/upsampledAmplitudes_{j['last']}.txt"
            ok = j["proc"].returncode == 0 and os.path.exists(last)
            print(f"{'done' if ok else 'FAILED'} {j['name']} in {time.time() - j['t0']:.0f}s "
                  f"(exit {j['proc'].returncode}, log {j['dir']}/upsample.log)", flush=True)
            if not ok:
                failed.append(j["name"])
            running.remove(j)
        time.sleep(1)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
