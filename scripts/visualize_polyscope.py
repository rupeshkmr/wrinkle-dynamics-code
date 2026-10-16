#!/usr/bin/env python3
"""
Polyscope visual comparison: reference vs our re-run, for a given config.

Shows:
  - "reference" and "ours" meshes (amplitude, wrinkle frequency dphi, grad-phi vectors)
  - "diff" mesh colored by:
        * position diff  : per-vertex  |reference_pos - ours_pos|
        * amplitude diff : per-vertex  |reference_amp - ours_amp|
        * frequency diff : per-face    |reference_dphi - ours_dphi|
  - obstacles (spheres)

Usage:
    conda activate crom
    python visualize_polyscope.py <config_name>
"""
import os
import re
import sys
import json
import numpy as np
import polyscope as ps
import polyscope.imgui as psim
from glob import glob

ROOT = os.path.dirname(os.path.abspath(__file__))
config_name = sys.argv[1] if len(sys.argv) > 1 else "1773844743"

# ---- resolve our run dir from RUN_MAP.txt ----
mapping = {}
with open(os.path.join(ROOT, "checkpoints", "RUN_MAP.txt")) as f:
    for line in f:
        m = re.search(r"^(\d+)\s+exit=0\s+checkpoint=checkpoints/(\d+)", line)
        if m:
            mapping[m.group(1)] = m.group(2)

if config_name not in mapping:
    sys.exit(f"Config '{config_name}' not found in RUN_MAP.txt. Available: {sorted(mapping)}")

REF = os.path.join(ROOT, "checkpoints_transfer/checkpoints", config_name)
MINE = os.path.join(ROOT, "checkpoints", mapping[config_name])
print(f"config={config_name}  reference={REF}  ours={MINE}")

def _num(p):
    return int(re.search(r"(\d+)$", p).group(1))

def sorted_files(base, what):
    return sorted(glob(os.path.join(base, what, "*.csv.*")), key=_num)

def load_frames(base, what, transpose=False, ndmin=2):
    out = []
    for f in sorted_files(base, what):
        d = np.loadtxt(f, delimiter=",", ndmin=ndmin)
        out.append(d.T if transpose else d)
    return out

# ---- load data ----
faces = np.loadtxt(os.path.join(REF, "faces", "faces.csv"), delimiter=",").astype(np.int64)
pos_ref = load_frames(REF, "positions", transpose=True)
pos_mine = load_frames(MINE, "positions", transpose=True)
amp_ref = load_frames(REF, "amplitudes", ndmin=1)
amp_mine = load_frames(MINE, "amplitudes", ndmin=1)
dphi_ref = load_frames(REF, "dphisPerFace", transpose=True)   # list of (F,2)
dphi_mine = load_frames(MINE, "dphisPerFace", transpose=True)
N = min(len(pos_ref), len(pos_mine))

# shared amplitude / wrinkle-frequency color ranges
amp_all = np.concatenate([np.concatenate(amp_ref), np.concatenate(amp_mine)])
amp_lo, amp_hi = float(amp_all.min()), float(amp_all.max())
if dphi_ref and dphi_mine:
    wf = np.concatenate([np.linalg.norm(np.concatenate(dphi_ref), axis=1),
                         np.linalg.norm(np.concatenate(dphi_mine), axis=1)])
    wf_lo, wf_hi = float(wf.min()), float(wf.max())
else:
    wf_lo, wf_hi = 0.0, 1.0

# grad-phi (3D) per face, batched (same math as viz_helpers.computeFaceGradPhiFromDphiPerFace_batched)
def compute_grad_phi(dphi, faces, pos):
    vi = pos[faces[:, 0]]; vj = pos[faces[:, 1]]; vk = pos[faces[:, 2]]
    BX = vj - vi; BY = vk - vi
    a = np.einsum('fi,fi->f', BX, BX)
    b = np.einsum('fi,fi->f', BX, BY)
    c = np.einsum('fi,fi->f', BY, BY)
    det = a * c - b * b
    d0 = dphi[:, 0]; d1 = dphi[:, 1]
    g0 = (c * d0 - b * d1) / det
    g1 = (-b * d0 + a * d1) / det
    return BX * g0[:, None] + BY * g1[:, None]

# ---- obstacles (SPHERE types only) ----
def sphere_mesh(center, radius, res=40):
    th = np.linspace(0, np.pi, res); ph = np.linspace(0, 2*np.pi, res)
    th, ph = np.meshgrid(th, ph)
    x = center[0] + radius*np.sin(th)*np.cos(ph)
    y = center[1] + radius*np.sin(th)*np.sin(ph)
    z = center[2] + radius*np.cos(th)
    V = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)
    F = []
    for i in range(res-1):
        for j in range(res-1):
            p1 = i*res+j; p2 = p1+1; p3 = p1+res; p4 = p3+1
            F += [[p1, p2, p3], [p2, p4, p3]]
    return V, np.array(F, dtype=np.int64)

obstacles = []
cfg = json.load(open(os.path.join(ROOT, "configs", config_name + ".json")))
for o in cfg.get("obstacles", []):
    t = o.get("type", "")
    if t in ("SPHERE", "SPHERE_"):
        obstacles.append(sphere_mesh(o["center"], o["radius"]))
    else:
        print(f"Note: obstacle type {t!r} not drawn (only spheres are rendered).")

# ---- polyscope ----
ps.init()
ps.set_up_dir("y_up")
ps.set_ground_plane_mode("none")

mesh_ref = ps.register_surface_mesh("reference", pos_ref[0], faces)
mesh_mine = ps.register_surface_mesh("ours", pos_mine[0], faces)
mesh_diff = ps.register_surface_mesh("diff", pos_mine[0], faces)
for i, (sphV, sphF) in enumerate(obstacles):
    ps.register_surface_mesh(f"obstacle_{i}", sphV, sphF, transparency=0.6)

def add_quantities(mesh, pos, amp, dphi):
    mesh.add_scalar_quantity("amplitude", amp, defined_on="vertices", vminmax=(amp_lo, amp_hi))
    if dphi is not None:
        mesh.add_scalar_quantity("dphi_x", dphi[:, 0], defined_on="faces")
        mesh.add_scalar_quantity("dphi_y", dphi[:, 1], defined_on="faces")
        mesh.add_scalar_quantity("wrinkle_freq |dphi|", np.linalg.norm(dphi, axis=1),
                                 defined_on="faces", vminmax=(wf_lo, wf_hi))
        mesh.add_vector_quantity("grad_phi", compute_grad_phi(dphi, faces, pos), defined_on="faces")

def add_diff_quantities(mesh_diff, kr, km):
    # position diff (per vertex)
    mesh_diff.add_scalar_quantity("pos diff |Δx|",
        np.linalg.norm(pos_ref[kr] - pos_mine[km], axis=1), defined_on="vertices")
    # amplitude diff (per vertex)
    mesh_diff.add_scalar_quantity("amplitude diff |Δa|",
        np.abs(amp_ref[kr] - amp_mine[km]), defined_on="vertices")
    # frequency diff (per face)
    if dphi_ref and dphi_mine:
        mesh_diff.add_scalar_quantity("frequency diff |Δdphi|",
            np.linalg.norm(dphi_ref[kr] - dphi_mine[km], axis=1), defined_on="faces")

add_quantities(mesh_ref, pos_ref[0], amp_ref[0], dphi_ref[0] if dphi_ref else None)
add_quantities(mesh_mine, pos_mine[0], amp_mine[0], dphi_mine[0] if dphi_mine else None)
add_diff_quantities(mesh_diff, 0, 0)

state = {"frame": 0, "play": False, "align": True}

def callback():
    if psim.Button("Play/Pause"):
        state["play"] = not state["play"]
    changed, f = psim.SliderInt("Frame", state["frame"], v_min=0, v_max=N-1)
    if changed:
        state["frame"] = f
    k = state["frame"]
    kr = k
    km = max(0, k - 1) if state["align"] else k   # ours is 1 step ahead

    mesh_ref.update_vertex_positions(pos_ref[kr])
    mesh_mine.update_vertex_positions(pos_mine[km])
    mesh_diff.update_vertex_positions(pos_mine[km])

    add_quantities(mesh_ref, pos_ref[kr], amp_ref[kr], dphi_ref[kr] if dphi_ref else None)
    add_quantities(mesh_mine, pos_mine[km], amp_mine[km], dphi_mine[km] if dphi_mine else None)
    add_diff_quantities(mesh_diff, kr, km)

    if state["play"]:
        state["frame"] = (state["frame"] + 1) % N

ps.set_user_callback(callback)
print(f"Loaded {N} frames. ours[k] aligns with reference[k+1] (set state['align']=False to disable).")
ps.show()
