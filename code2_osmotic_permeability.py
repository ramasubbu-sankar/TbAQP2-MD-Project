import os, traceback
import numpy as np
import MDAnalysis as mda
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from typing import Dict, Tuple

AQP_FOLDERS = [
    "/media/dsk_3/param/Aquaporin_water/tb_crys_50_ns/gromacs/50_ns",
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

Z_BELOW_NM        = 0.5
Z_ABOVE_NM        = 1.5
CHANNEL_LENGTH_NM = Z_BELOW_NM + Z_ABOVE_NM

CYLINDER_RADIUS_NM = 0.5

TIME_MIN_PS = 0.0
TIME_MAX_PS = 50000.0

MAX_LAG_PS     = 500.0
FIT_TAU_MIN_PS =  10.0

VW_CM3   = 2.99e-23
PLOT_DPI = 1200

OUTDIR_NAME = "pf_outputs_1ps_50ns_boundary_fix"

try:
    _script_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    _script_dir = os.getcwd()
GLOBAL_SUMMARY = os.path.join(_script_dir, "all_aqp_pf_summary_1ps_50ns_boundary_fix.csv")

def load_xvg_com(path: str) -> Tuple[np.ndarray, np.ndarray]:
    times, xyz = [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(("@", "#")):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            times.append(float(parts[0]))
            xyz.append((float(parts[1]), float(parts[2]), float(parts[3])))
    return np.array(times, dtype=np.float64), np.array(xyz, dtype=np.float64)

class COMByTime:
    def __init__(self, times_ps, centers_nm):
        self.t   = times_ps
        self.xyz = centers_nm

    def at(self, time_ps: float) -> np.ndarray:
        idx = int(np.clip(np.searchsorted(self.t, time_ps), 1, len(self.t) - 1))
        if abs(self.t[idx] - time_ps) < abs(self.t[idx - 1] - time_ps):
            return self.xyz[idx]
        return self.xyz[idx - 1]

def min_image(delta, box_len: float):
    if box_len <= 0.0:
        return delta
    return delta - box_len * np.round(delta / box_len)

def compute_msd(n_values: np.ndarray, dt_ps: float, max_lag_ps: float):
    N = n_values.size
    if N < 3:
        return np.array([]), np.array([])
    max_k = max(1, min(int(max_lag_ps / dt_ps), N - 1))
    taus  = np.arange(1, max_k + 1, dtype=np.int64) * dt_ps
    msd   = np.array([np.mean((n_values[k:] - n_values[:-k]) ** 2)
                      for k in range(1, max_k + 1)])
    return taus, msd

def linear_fit(x, y):
    if x.size < 2:
        return np.nan, np.nan, np.nan
    A    = np.vstack([x, np.ones_like(x)]).T
    m, b = np.linalg.lstsq(A, y, rcond=None)[0]
    y_p  = m * x + b
    ss_res = np.sum((y - y_p) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1.0 - (ss_res / ss_tot if ss_tot > 0 else np.nan)
    return m, b, r2

def run_folder(folder: str, label: str, global_rows: list):
    print(f"\n{'='*60}")
    print(f"[AQP] {label}  (1 ps / 50 ns run, boundary-fix)")
    print(f"      {folder}")
    print(f"{'='*60}")

    if not os.path.isdir(folder):
        print(f"  [SKIP] folder not found.")
        return

    topo = os.path.join(folder, GRO_FILE)
    traj = os.path.join(folder, XTC_FILE)
    missing = []
    if not os.path.isfile(topo):
        missing.append(GRO_FILE)
    if not os.path.isfile(traj):
        missing.append(XTC_FILE)
    for chain, xvg in XVG_FILES.items():
        if not os.path.isfile(os.path.join(folder, xvg)):
            missing.append(xvg)
    if missing:
        print(f"  [SKIP] missing files: {', '.join(missing)}")
        return

    outdir = os.path.join(folder, OUTDIR_NAME)
    os.makedirs(outdir, exist_ok=True)

    com_map: Dict[str, COMByTime] = {}
    for chain, xvg in XVG_FILES.items():
        t_ps, centers_nm = load_xvg_com(os.path.join(folder, xvg))
        com_map[chain]   = COMByTime(t_ps, centers_nm)
        print(f"  xvg [{chain}] {xvg}: {len(t_ps)} frames "
              f"(t = {t_ps[0]:.0f}-{t_ps[-1]:.0f} ps)")

    try:
        U     = mda.Universe(topo, traj)
        water = U.select_atoms(WATER_SEL)
    except Exception as e:
        print(f"  [ERROR] MDAnalysis: {e}")
        return

    if water.n_atoms == 0:
        print(f"  [SKIP] WATER_SEL '{WATER_SEL}' matched 0 atoms.")
        return

    mon_data: Dict[str, dict] = {
        chain: {
            "times_ps": [],
            "n_values": [],
            "n_curr":   0.0,
            "prev_mask": None,
            "prev_z_nm": None,
        }
        for chain in XVG_FILES
    }

    n_frames = 0
    for ts in U.trajectory:
        t_ps = float(ts.time)
        if TIME_MIN_PS is not None and t_ps < TIME_MIN_PS:
            continue
        if TIME_MAX_PS is not None and t_ps > TIME_MAX_PS:
            continue
        n_frames += 1

        Lx = float(ts.dimensions[0]) / 10.0
        Ly = float(ts.dimensions[1]) / 10.0
        Lz = float(ts.dimensions[2]) / 10.0

        pos  = water.positions.astype(np.float64) / 10.0
        x_nm = pos[:, 0]; y_nm = pos[:, 1]; z_nm = pos[:, 2]

        for chain, com in com_map.items():
            cx, cy, cz = com.at(t_ps)

            dx   = min_image(x_nm - cx, Lx)
            dy   = min_image(y_nm - cy, Ly)
            z_lo = cz - Z_BELOW_NM
            z_hi = cz + Z_ABOVE_NM
            r2   = dx*dx + dy*dy
            mask_now = (r2 <= CYLINDER_RADIUS_NM**2) & (z_nm >= z_lo) & (z_nm <= z_hi)

            D = mon_data[chain]
            prev_mask = D["prev_mask"]
            prev_z    = D["prev_z_nm"]

            if prev_mask is not None:
                both    = prev_mask & mask_now
                entered = (~prev_mask) & mask_now
                exited  = prev_mask & (~mask_now)

                contrib = 0.0

                if np.any(both):
                    dz_both = min_image(z_nm[both] - prev_z[both], Lz)
                    contrib += float(np.sum(dz_both))

                if np.any(entered):
                    zp = prev_z[entered]
                    zn = z_nm[entered]
                    from_below = zp <= z_lo
                    from_above = zp >= z_hi
                    dz_enter = np.zeros_like(zn)
                    dz_enter[from_below] = zn[from_below] - z_lo
                    dz_enter[from_above] = zn[from_above] - z_hi
                    contrib += float(np.sum(dz_enter))

                if np.any(exited):
                    zp = prev_z[exited]
                    zn = z_nm[exited]
                    to_below = zn <= z_lo
                    to_above = zn >= z_hi
                    dz_exit = np.zeros_like(zp)
                    dz_exit[to_below] = z_lo - zp[to_below]
                    dz_exit[to_above] = z_hi - zp[to_above]
                    contrib += float(np.sum(dz_exit))

                D["n_curr"] += contrib / CHANNEL_LENGTH_NM

            D["times_ps"].append(t_ps)
            D["n_values"].append(D["n_curr"])

            D["prev_mask"] = mask_now
            D["prev_z_nm"] = z_nm.copy()

    print(f"  frames processed: {n_frames}")
    if n_frames < 10:
        print(f"  [WARN] very few frames ??? check TIME_MIN/MAX_PS.")
        return

    t_arr = np.array(mon_data[next(iter(mon_data))]["times_ps"])
    if t_arr.size < 2:
        print(f"  [SKIP] not enough frames for MSD.")
        return
    dt_ps = float(np.median(np.diff(t_arr)))
    print(f"  detected dt = {dt_ps:.2f} ps")

    folder_rows = []
    for chain, D in mon_data.items():
        times_ps = np.array(D["times_ps"], dtype=np.float64)
        n_vals   = np.array(D["n_values"],  dtype=np.float64)

        np.savetxt(
            os.path.join(outdir, f"{chain}_n_values.txt"),
            np.column_stack([times_ps, n_vals]),
            fmt="%.6f", header="time_ps n"
        )

        eff_lag      = min(MAX_LAG_PS, (times_ps[-1] - times_ps[0]) * 0.25)
        taus_ps, msd = compute_msd(n_vals, dt_ps, eff_lag)

        np.savetxt(
            os.path.join(outdir, f"{chain}_msd.txt"),
            np.column_stack([taus_ps, msd]),
            fmt="%.6f", header="tau_ps msd"
        )

        fit_mask = (taus_ps >= FIT_TAU_MIN_PS) & (taus_ps <= eff_lag)
        slope, intercept, r2 = linear_fit(taus_ps[fit_mask], msd[fit_mask])

        Dn      = 0.5 * (slope / 1.0e-12)
        pf_cm3  = VW_CM3 * Dn
        pf_1e14 = pf_cm3 / 1.0e-14

        plt.figure()
        plt.plot(times_ps, n_vals, lw=1.2)
        plt.xlabel("Time (ps)")
        plt.ylabel("n (cumulative)")
        plt.title(f"{label} ??? monomer {chain}  n(t)  [1ps/50ns, boundary-fix]")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"{chain}_n_vs_t.jpeg"), dpi=PLOT_DPI)
        plt.close()

        plt.figure()
        plt.scatter(taus_ps, msd, s=8, label="MSD")
        if np.isfinite(slope):
            plt.plot(taus_ps, slope * taus_ps + intercept, lw=1.5,
                     label=f"slope={slope:.4e}  R??={r2:.3f}")
        plt.xlabel("?? (ps)")
        plt.ylabel("MSD of n")
        plt.title(f"{label} ??? monomer {chain}  |  lag ??? {eff_lag:.0f} ps  [1ps/50ns, boundary-fix]")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"{chain}_msd_fit.jpeg"), dpi=PLOT_DPI)
        plt.close()

        print(f"  [{chain}] R??={r2:.3f}  Dn={Dn:.3e} s?????  "
              f"pf={pf_cm3:.3e} cm??/s  ({pf_1e14:.3f} ??10????????)")

        row = [label, chain, len(times_ps), f"{dt_ps:.2f}", f"{eff_lag:.1f}",
               f"{slope:.6e}", f"{r2:.4f}", f"{Dn:.6e}",
               f"{pf_cm3:.6e}", f"{pf_1e14:.4f}"]
        folder_rows.append(row)
        global_rows.append(row)

    if not folder_rows:
        print("  [SKIP] no monomers produced a valid fit.")
        return

    hdr = ("monomer,frames,dt_ps,max_lag_ps,slope_dndn_per_ps,"
           "R2,Dn_per_s,pf_cm3_per_s,pf_(x1e-14)")
    np.savetxt(
        os.path.join(outdir, "pf_summary.csv"),
        np.array([r[1:] for r in folder_rows], dtype=object),
        fmt="%s", delimiter=",", header=hdr, comments=""
    )
    print(f"  saved ??? {outdir}/pf_summary.csv")

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

    hdr = ("aqp,monomer,frames,dt_ps,max_lag_ps,"
           "slope_dndn_per_ps,R2,Dn_per_s,pf_cm3_per_s,pf_(x1e-14)")
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
