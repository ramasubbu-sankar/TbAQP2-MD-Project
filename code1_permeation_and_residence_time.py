import os, traceback
import numpy as np
import pandas as pd
import MDAnalysis as mda

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

WATER_SEL  = "resname TIP3 and name OH2"

Z_BELOW_NM = 0.5
Z_ABOVE_NM = 1.3

RADIUS_NM  = 0.5

MIN_INSIDE_FRAC         = 0.80
MIN_MIDDLE_FRAMES       = 3
MAX_CONSEC_OUTSIDE_GAP  = 15
REQUIRE_MONOTONIC_TREND = True

USE_MIN_IMAGE_XY = False

try:
    _script_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    _script_dir = os.getcwd()
GLOBAL_SUMMARY = os.path.join(_script_dir, "all_aqp_permeation_summary.csv")

ANG_TO_NM = 0.1

def parse_xvg(fn):
    data = []
    with open(fn) as fh:
        for ln in fh:
            if not ln or ln[0] in ("@", "#"):
                continue
            sp = ln.split()
            if len(sp) < 4:
                continue
            data.append((float(sp[0]), float(sp[1]),
                         float(sp[2]), float(sp[3])))
    arr = np.asarray(data, dtype=np.float64)
    if arr.size == 0:
        raise ValueError(f"No numeric data in {fn}")
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]

def min_image_delta(d, L):
    return d - L * np.round(d / L)

def zones_from_z(z_nm, z_lo, z_hi):
    zc = np.full(z_nm.shape, -1, dtype=np.int8)
    zc[z_nm < z_lo] = 0
    zc[(z_nm >= z_lo) & (z_nm <= z_hi)] = 1
    zc[z_nm > z_hi] = 2
    return zc

def run_monomer(monomer_id, xvg_path, u, outdir):
    t_com, x_com, y_com, z_com = parse_xvg(xvg_path)

    water  = u.select_atoms(WATER_SEL)
    nW     = water.n_atoms
    resids = water.resids.astype(np.int64)

    origin_zone   = np.full(nW, -1,       dtype=np.int8)
    origin_time   = np.full(nW, np.nan,   dtype=np.float64)
    origin_z      = np.full(nW, np.nan,   dtype=np.float64)
    armed         = np.zeros(nW,          dtype=np.int8)
    t_start       = np.full(nW, np.nan,   dtype=np.float64)
    z_start       = np.full(nW, np.nan,   dtype=np.float64)
    z_arm0        = np.full(nW, np.nan,   dtype=np.float64)
    inside_frames = np.zeros(nW,          dtype=np.int32)
    total_frames  = np.zeros(nW,          dtype=np.int32)
    middle_frames = np.zeros(nW,          dtype=np.int32)
    outside_run   = np.zeros(nW,          dtype=np.int32)
    last_zone     = np.full(nW, -1,       dtype=np.int8)
    last_time     = np.full(nW, np.nan,   dtype=np.float64)
    last_z_cont   = np.full(nW, np.nan,   dtype=np.float64)
    last_z_raw    = np.full(nW, np.nan,   dtype=np.float64)
    z_offset      = np.zeros(nW,          dtype=np.float64)

    events = []

    for frm_idx, ts in enumerate(u.trajectory):
        t_ps  = float(ts.time)
        Lx    = float(ts.dimensions[0]) * ANG_TO_NM
        Ly    = float(ts.dimensions[1]) * ANG_TO_NM
        Lz    = float(ts.dimensions[2]) * ANG_TO_NM

        pos   = water.positions.astype(np.float64) * ANG_TO_NM
        x_w   = pos[:, 0]
        y_w   = pos[:, 1]
        z_raw = pos[:, 2]

        cx = float(np.interp(t_ps, t_com, x_com))
        cy = float(np.interp(t_ps, t_com, y_com))
        cz = float(np.interp(t_ps, t_com, z_com))

        if USE_MIN_IMAGE_XY:
            dx = min_image_delta(x_w - cx, Lx)
            dy = min_image_delta(y_w - cy, Ly)
        else:
            dx = x_w - cx
            dy = y_w - cy
        inside = (dx*dx + dy*dy) <= (RADIUS_NM * RADIUS_NM)

        if frm_idx == 0:
            last_z_raw[:] = z_raw
        else:
            dz_raw = z_raw - last_z_raw
            z_offset[dz_raw >  (Lz / 2)] -= Lz
            z_offset[dz_raw < -(Lz / 2)] += Lz
            last_z_raw[:] = z_raw
        z_cont = z_raw + z_offset

        z_lo  = cz - Z_BELOW_NM
        z_hi  = cz + Z_ABOVE_NM
        zcode = zones_from_z(z_raw, z_lo, z_hi)

        efb = (last_zone == 0) & (zcode == 1)
        efa = (last_zone == 2) & (zcode == 1)
        if np.any(efb):
            idx = np.where(efb)[0]
            origin_zone[idx] = 0
            origin_time[idx] = last_time[idx]
            origin_z[idx]    = last_z_cont[idx]
        if np.any(efa):
            idx = np.where(efa)[0]
            origin_zone[idx] = 2
            origin_time[idx] = last_time[idx]
            origin_z[idx]    = last_z_cont[idx]

        can_arm = (armed == 0) & (zcode == 1) & inside & (origin_zone != -1)
        if np.any(can_arm):
            up = can_arm & (origin_zone == 0)
            dn = can_arm & (origin_zone == 2)
            armed[up] = +1;  armed[dn] = -1
            t_start[can_arm]       = origin_time[can_arm]
            z_start[can_arm]       = origin_z[can_arm]
            z_arm0[can_arm]        = origin_z[can_arm]
            inside_frames[can_arm] = 0
            total_frames[can_arm]  = 0
            middle_frames[can_arm] = 0
            outside_run[can_arm]   = 0

        changed = zcode != last_zone
        if np.any(changed):
            last_zone[changed]   = zcode[changed]
            last_time[changed]   = t_ps
            last_z_cont[changed] = z_cont[changed]

        left_mid = (zcode != 1) & (origin_zone != -1) & (armed == 0)
        origin_zone[left_mid] = -1

        active = armed != 0
        if np.any(active):
            total_frames[active]               += 1
            inside_frames[active & inside]     += 1
            middle_frames[active & (zcode==1)] += 1
            outside_run[active & (~inside)]    += 1
            outside_run[active & inside]        = 0
            drop = active & (outside_run > MAX_CONSEC_OUTSIDE_GAP)
            armed[drop] = 0;  origin_zone[drop] = -1

        reach_up   = (armed == +1) & (zcode == 2)
        reach_down = (armed == -1) & (zcode == 0)
        success    = reach_up | reach_down
        if np.any(success):
            for i in np.where(success)[0]:
                dur = t_ps - t_start[i]
                if dur <= 0:
                    armed[i] = 0;  origin_zone[i] = -1;  continue

                frac_in = inside_frames[i] / max(1, total_frames[i])
                ok_frac = frac_in >= MIN_INSIDE_FRAC
                ok_mid  = middle_frames[i] >= MIN_MIDDLE_FRAMES
                ok_mono = True
                if REQUIRE_MONOTONIC_TREND:
                    dz_net  = z_cont[i] - z_arm0[i]
                    ok_mono = (dz_net > 0) if reach_up[i] else (dz_net < 0)

                if ok_frac and ok_mid and ok_mono:
                    events.append({
                        "monomer":     monomer_id,
                        "resid":       int(resids[i]),
                        "direction":   "UPWARD" if reach_up[i] else "DOWNWARD",
                        "t_start_ps":  float(t_start[i]),
                        "t_end_ps":    float(t_ps),
                        "duration_ps": float(dur),
                        "z_start_nm":  float(z_start[i]),
                        "z_end_nm":    float(z_cont[i]),
                        "inside_frac": float(frac_in),
                    })

                armed[i] = 0;  origin_zone[i] = -1

    cols = ["monomer","resid","direction","t_start_ps","t_end_ps",
            "duration_ps","z_start_nm","z_end_nm","inside_frac"]
    df = pd.DataFrame(events, columns=cols)
    df.sort_values(["t_start_ps","resid"], inplace=True)
    df.to_csv(os.path.join(outdir, f"permeations_monomer_{monomer_id}.csv"),
              index=False)
    print(f"    [{monomer_id}] {len(df)} events")
    return df

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
    if not os.path.isfile(topo):  missing.append(GRO_FILE)
    if not os.path.isfile(traj):  missing.append(XTC_FILE)
    for ch, xvg in XVG_FILES.items():
        if not os.path.isfile(os.path.join(folder, xvg)):
            missing.append(xvg)
    if missing:
        print(f"  [SKIP] missing: {', '.join(missing)}")
        return

    outdir = os.path.join(folder, "permeation_outputs")
    os.makedirs(outdir, exist_ok=True)

    try:
        u     = mda.Universe(topo, traj)
        water = u.select_atoms(WATER_SEL)
    except Exception as e:
        print(f"  [ERROR] MDAnalysis: {e}")
        return
    if water.n_atoms == 0:
        print(f"  [SKIP] WATER_SEL matched 0 atoms.")
        return

    mon_dfs  = []
    mon_rows = []
    for ch, xvg in XVG_FILES.items():
        xvg_path = os.path.join(folder, xvg)
        df = run_monomer(ch, xvg_path, u, outdir)
        mon_dfs.append(df)

        n_total    = len(df)
        n_upward   = int((df["direction"] == "UPWARD").sum())   if n_total else 0
        n_downward = int((df["direction"] == "DOWNWARD").sum()) if n_total else 0
        avg_dur    = float(df["duration_ps"].mean()) if n_total else float("nan")

        mon_rows.append({
            "aqp":       label,
            "monomer":   ch,
            "n_total":   n_total,
            "n_upward":  n_upward,
            "n_downward":n_downward,
            "avg_duration_ps": avg_dur,
        })
        global_rows.append(mon_rows[-1])

    if mon_dfs:
        combined = pd.concat(mon_dfs, ignore_index=True)
        combined.sort_values(["monomer","t_start_ps","resid"], inplace=True)
        combined.to_csv(os.path.join(outdir, "permeations_all.csv"), index=False)

    summary_df = pd.DataFrame(mon_rows)

    avg_row = {
        "aqp":            label,
        "monomer":        "AVERAGE",
        "n_total":        summary_df["n_total"].mean(),
        "n_upward":       summary_df["n_upward"].mean(),
        "n_downward":     summary_df["n_downward"].mean(),
        "avg_duration_ps":summary_df["avg_duration_ps"].mean(),
    }
    summary_df = pd.concat(
        [summary_df, pd.DataFrame([avg_row])], ignore_index=True
    )
    summary_df.to_csv(
        os.path.join(outdir, "permeation_summary.csv"), index=False
    )

    print(f"  avg permeations/monomer : {avg_row['n_total']:.2f}")
    print(f"  saved ??? {outdir}/")

    global_rows.append(avg_row)

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

    global_df = pd.DataFrame(global_rows)
    global_df.to_csv(GLOBAL_SUMMARY, index=False)

    print(f"\n{'='*60}")
    print(f"[DONE] Global summary => {GLOBAL_SUMMARY}")
    print(f"       {len(global_df)} rows written.")

    avg_df = global_df[global_df["monomer"] == "AVERAGE"][
        ["aqp", "n_total", "n_upward", "n_downward", "avg_duration_ps"]
    ].reset_index(drop=True)
    print("\n?????? Average permeations per monomer ???????????????????????????????????????????????????????????????")
    print(avg_df.to_string(index=False))

if __name__ == "__main__":
    main()
