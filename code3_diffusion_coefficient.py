import os, traceback
import numpy as np
from scipy.stats import linregress
import MDAnalysis as mda
import warnings
warnings.filterwarnings("ignore")

AQP_FOLDERS = [

    "/media/maajid/dsk_4/mj/Paper_review_simulations/rep_2_tryp_aqp_2/charmm-gui-1057003267/gromacs",

]

GRO_FILE  = "water_oh.gro"
XTC_FILE  = "water_oh_last_500.xtc"
XVG_FILES = {
    "A": "npa_last_500_chain_A.xvg",
    "B": "npa_last_500_chain_B.xvg",
    "C": "npa_last_500_chain_C.xvg",
    "D": "npa_last_500_chain_D.xvg",
}

WATER_SEL = "resname TIP3 and name OH2"

Z_BELOW_NM = 0.5
Z_ABOVE_NM = 1.5

CYLINDER_RADIUS_NM = 0.5

USE_MONOMER_CENTER_Z = True

TIME_MIN_PS = 500000.0
TIME_MAX_PS = 1000000.0

MAX_LAG_PS   = 500.0
FIT_START_PS =  50.0
FIT_END_PS   = 400.0

try:
    _script_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    _script_dir = os.getcwd()
GLOBAL_SUMMARY = os.path.join(_script_dir, "all_aqp_dz_summary.csv")

def load_xvg(path):
    data = np.loadtxt(path, comments=("#", "@"))
    if data.ndim == 1:
        data = data[None, :]
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]

def min_image_1d(delta, boxlen):
    return delta - np.round(delta / boxlen) * boxlen

def compute_msd_1d_vectorised(resids_list, zvals_list, max_lag_frames, Lz_nm):
    n_frames = len(zvals_list)
    max_lag  = min(max_lag_frames, n_frames - 1)

    msd_acc    = np.zeros(max_lag + 1, dtype=np.float64)
    counts_acc = np.zeros(max_lag + 1, dtype=np.int64)

    tracks = {}
    for f in range(n_frames):
        rr = resids_list[f]
        zz = zvals_list[f]
        for r, z in zip(rr, zz):
            if r not in tracks:
                tracks[r] = ([], [])
            tracks[r][0].append(f)
            tracks[r][1].append(z)

    for r, (flist, zlist) in tracks.items():
        frames = np.asarray(flist, dtype=np.int32)
        zvals  = np.asarray(zlist, dtype=np.float64)
        n_obs  = len(frames)
        if n_obs < 2:
            continue

        fmap = {fr: idx for idx, fr in enumerate(frames)}

        for i in range(n_obs):
            f0  = frames[i]
            z0  = zvals[i]

            max_lag_here = min(max_lag, n_frames - 1 - int(f0))
            if max_lag_here < 1:
                continue

            target_frames = f0 + np.arange(1, max_lag_here + 1, dtype=np.int32)

            j_indices = np.array([fmap.get(tf, -1) for tf in target_frames])
            valid     = j_indices >= 0
            if not np.any(valid):
                continue
            lags_valid = np.where(valid)[0] + 1
            j_valid    = j_indices[valid]
            dz         = zvals[j_valid] - z0

            dz         = min_image_1d(dz, Lz_nm)
            np.add.at(msd_acc,    lags_valid, dz * dz)
            np.add.at(counts_acc, lags_valid, 1)

    valid_mask = counts_acc > 0
    msd_out = np.zeros_like(msd_acc)
    msd_out[valid_mask] = msd_acc[valid_mask] / counts_acc[valid_mask]

    return np.arange(max_lag + 1, dtype=np.float64), msd_out, counts_acc

def safe_linfit(time_ps, msd, counts, t0, t1):
    mask = (time_ps >= t0) & (time_ps <= t1) & (counts > 0) & np.isfinite(msd)
    if np.sum(mask) < 3:
        return None
    slope, intercept, r, _, stderr = linregress(time_ps[mask], msd[mask])
    return {
        "slope":     slope,
        "intercept": intercept,
        "r2":        r * r,
        "stderr":    stderr,
        "npts":      int(np.sum(mask)),
    }

def run_folder(folder, label, global_rows):
    print(f"\n{'='*60}")
    print(f"[AQP] {label}")
    print(f"      {folder}")
    print(f"{'='*60}")

    if not os.path.isdir(folder):
        print("  [SKIP] folder not found.")
        return

    topo = os.path.join(folder, GRO_FILE)
    traj = os.path.join(folder, XTC_FILE)
    missing = []
    if not os.path.isfile(topo): missing.append(GRO_FILE)
    if not os.path.isfile(traj): missing.append(XTC_FILE)
    for ch, xvg in XVG_FILES.items():
        if not os.path.isfile(os.path.join(folder, xvg)):
            missing.append(xvg)
    if missing:
        print(f"  [SKIP] missing: {', '.join(missing)}")
        return

    outdir = os.path.join(folder, "dz_outputs_latest")
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(os.path.join(outdir, "msd"), exist_ok=True)

    com_data = {}
    for ch, xvg in XVG_FILES.items():
        t_ps, x, y, z = load_xvg(os.path.join(folder, xvg))
        com_data[ch]  = {"t": t_ps, "x": x, "y": y, "z": z}
        print(f"  xvg [{ch}] {xvg}: {len(t_ps)} frames")

    try:
        U     = mda.Universe(topo, traj)
        water = U.select_atoms(WATER_SEL)
    except Exception as e:
        print(f"  [ERROR] MDAnalysis: {e}")
        return

    if water.n_atoms == 0:
        print(f"  [SKIP] WATER_SEL matched 0 atoms.")
        return

    resids_all = water.resids

    mon_resids = {ch: [] for ch in XVG_FILES}
    mon_zvals  = {ch: [] for ch in XVG_FILES}

    Lz_nm = None

    frame_idx = 0
    n_frames  = 0
    for ts in U.trajectory:
        t_ps = float(ts.time)
        if TIME_MIN_PS is not None and t_ps < TIME_MIN_PS:
            frame_idx += 1
            continue
        if TIME_MAX_PS is not None and t_ps > TIME_MAX_PS:
            frame_idx += 1
            continue

        if Lz_nm is None:
            Lz_nm = float(ts.dimensions[2]) / 10.0

        Lx_nm = float(ts.dimensions[0]) / 10.0
        Ly_nm = float(ts.dimensions[1]) / 10.0

        pos_nm = water.positions.astype(np.float64) / 10.0
        x_w    = pos_nm[:, 0]
        y_w    = pos_nm[:, 1]
        z_w    = pos_nm[:, 2]

        for ch in XVG_FILES:
            cd  = com_data[ch]

            fi  = frame_idx
            if fi >= len(cd["t"]):
                fi = len(cd["t"]) - 1
            cx, cy, cz = cd["x"][fi], cd["y"][fi], cd["z"][fi]

            dx = min_image_1d(x_w - cx, Lx_nm)
            dy = min_image_1d(y_w - cy, Ly_nm)

            z_lo = cz - Z_BELOW_NM
            z_hi = cz + Z_ABOVE_NM

            r2   = dx*dx + dy*dy
            mask = (r2 <= CYLINDER_RADIUS_NM**2) & (z_w >= z_lo) & (z_w <= z_hi)

            resids_in = resids_all[mask]
            if USE_MONOMER_CENTER_Z:
                z_store = z_w[mask] - cz
            else:
                z_store = z_w[mask]

            mon_resids[ch].append(resids_in)
            mon_zvals[ch].append(z_store)

        frame_idx += 1
        n_frames  += 1

    print(f"  frames processed: {n_frames}")
    if n_frames < 10:
        print("  [WARN] very few frames.")
        return

    if Lz_nm is None:
        print("  [SKIP] could not read box dimensions.")
        return

    first_ch = next(iter(XVG_FILES))
    t_arr    = com_data[first_ch]["t"]
    dt_ps    = float(np.median(np.diff(t_arr)))
    max_lag_frames = int(round(MAX_LAG_PS / dt_ps))
    print(f"  dt={dt_ps:.2f} ps  max_lag={max_lag_frames} frames  Lz={Lz_nm:.3f} nm")

    folder_rows = []
    for ch in XVG_FILES:
        print(f"  [MSD] monomer {ch} ...")

        lag_idx, msd, counts = compute_msd_1d_vectorised(
            mon_resids[ch], mon_zvals[ch], max_lag_frames, Lz_nm
        )
        time_lags_ps = lag_idx * dt_ps

        msd_csv = os.path.join(outdir, "msd", f"{ch}_msd.csv")
        with open(msd_csv, "w") as f:
            f.write("lag_ps,msd_nm2,count\n")
            for t, m, c in zip(time_lags_ps, msd, counts):
                f.write(f"{t:.4f},{m:.10e},{int(c)}\n")

        fit = safe_linfit(time_lags_ps, msd, counts, FIT_START_PS, FIT_END_PS)
        if fit is None:
            print(f"    [{ch}] not enough valid points in fit window.")
            Dz_cm2s = Dz_err = r2 = float("nan")
        else:
            Dz_cm2s  = fit["slope"] / 2.0 * 1e-2
            Dz_err   = fit["stderr"] / 2.0 * 1e-2
            r2       = fit["r2"]
            print(f"    [{ch}] Dz={Dz_cm2s:.3e} ?? {Dz_err:.3e} cm??/s  R??={r2:.3f}  n={fit['npts']}")

        row = [label, ch, n_frames, f"{dt_ps:.2f}", f"{MAX_LAG_PS:.0f}",
               f"{FIT_START_PS:.0f}", f"{FIT_END_PS:.0f}",
               f"{Dz_cm2s:.6e}", f"{Dz_err:.6e}", f"{r2:.4f}",
               "rel" if USE_MONOMER_CENTER_Z else "abs"]
        folder_rows.append(row)
        global_rows.append(row)

    hdr = ("monomer,frames,dt_ps,max_lag_ps,fit_start_ps,fit_end_ps,"
           "Dz_cm2_s,Dz_err_cm2_s,R2,z_mode")
    np.savetxt(
        os.path.join(outdir, "dz_summary.csv"),
        np.array([r[1:] for r in folder_rows], dtype=object),
        fmt="%s", delimiter=",", header=hdr, comments=""
    )

    summary_txt = os.path.join(outdir, "dz_summary.txt")
    dz_vals = [float(r[7]) for r in folder_rows if r[7] != "nan"]
    with open(summary_txt, "w") as f:
        f.write(f"Longitudinal Dz ??? {label}\n")
        f.write("="*50 + "\n")
        f.write(f"frames: {n_frames}  dt: {dt_ps:.2f} ps\n")
        f.write(f"Z-window: COM_z ?? [{Z_BELOW_NM}, {Z_ABOVE_NM}] nm  "
                f"r_cyl: {CYLINDER_RADIUS_NM} nm\n")
        f.write(f"MSD max lag: {MAX_LAG_PS} ps  "
                f"fit: [{FIT_START_PS}, {FIT_END_PS}] ps\n")
        f.write(f"Z mode: {'monomer-centred' if USE_MONOMER_CENTER_Z else 'absolute'}\n\n")
        for r in folder_rows:
            f.write(f"  {r[1]}: Dz={r[7]} ?? {r[8]} cm??/s  R??={r[9]}\n")
        if dz_vals:
            f.write(f"\n  Mean Dz = {np.mean(dz_vals):.3e} ?? {np.std(dz_vals):.3e} cm??/s\n")

    print(f"  saved ??? {outdir}/dz_summary.csv")

def main():
    global_rows = []

    for folder in AQP_FOLDERS:
        label = os.path.basename(folder.rstrip("/"))
        try:
            run_folder(folder, label, global_rows)
        except Exception:
            print(f"\n  [ERROR] unhandled exception in {folder}:")
            traceback.print_exc()
            print("  Continuing with next folder...\n")

    if not global_rows:
        print("\n[DONE] No results ??? check folder paths and filenames.")
        return

    hdr = ("aqp,monomer,frames,dt_ps,max_lag_ps,fit_start_ps,fit_end_ps,"
           "Dz_cm2_s,Dz_err_cm2_s,R2,z_mode")
    np.savetxt(
        GLOBAL_SUMMARY,
        np.array(global_rows, dtype=object),
        fmt="%s", delimiter=",", header=hdr, comments=""
    )
    print(f"\n{'='*60}")
    print(f"[DONE] Global summary => {GLOBAL_SUMMARY}")
    print(f"       {len(global_rows)} monomer rows written.")

if __name__ == "__main__":
    main()
