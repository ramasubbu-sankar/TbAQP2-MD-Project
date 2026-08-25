import os
import numpy as np
import pandas as pd
import MDAnalysis as mda

SYSTEMS = [
    {
        "name": "AQP1_rep2",
        "base_dir": "/media/maajid/dsk_4/mj/Paper_review_simulations/aqp_1_replicate_2/charmm-gui-8124847806/gromacs",
        "topology": "water_oh_pro.gro",
        "trajectory": "water_oh_pro_last_500.xtc",
        "perm_file": "permeations_all.csv",

        "sf_resids": [58, 182, 191, 197],
    },
    {
        "name": "TbAQP2_crystal",
        "base_dir": "/media/maajid/dsk_4/mj/Paper_review_simulations/tbaqp_2_8jy7/charmm-gui-8120238087/gromacs",
        "topology": "water_oh_pro.gro",
        "trajectory": "water_oh_pro_last_500.xtc",
        "perm_file": "permeations_all.csv",

        "sf_resids": [110, 249, 258, 264],
    },
    {
        "name": "TbAQP2_rep2",
        "base_dir": "/media/maajid/dsk_4/mj/Paper_review_simulations/rep_2_tryp_aqp_2/charmm-gui-1057003267/gromacs",
        "topology": "water_oh_pro.gro",
        "trajectory": "water_oh_pro_last_500.xtc",
        "perm_file": "permeations_all.csv",
        "sf_resids": [110, 249, 258, 264],
    },
]

WATER_SEL = "resname TIP3 and name OH2"
DIST_CUTOFF_NM = 0.40

def analyze_system(sysinfo):
    name       = sysinfo["name"]
    base_dir   = sysinfo["base_dir"]
    topology   = os.path.join(base_dir, sysinfo["topology"])
    trajectory = os.path.join(base_dir, sysinfo["trajectory"])
    perm_file  = os.path.join(base_dir, sysinfo["perm_file"])
    sf_resids  = sysinfo["sf_resids"]

    outdir = os.path.join(base_dir, "SF_residue_contacts")
    os.makedirs(outdir, exist_ok=True)

    print(f"\n{'='*70}\n  SYSTEM: {name}\n{'='*70}")

    for p in (topology, trajectory, perm_file):
        if not os.path.exists(p):
            raise FileNotFoundError(f"[{name}] Missing file: {p}")

    df = pd.read_csv(perm_file)
    if df.empty:
        raise ValueError(f"[{name}] No permeation events found!")

    print(f"[{name}] Total permeation events: {len(df)}")

    df["contact_any_ps"] = 0.0
    for r in sf_resids:
        df[f"contact_res{r}_ps"] = 0.0
    df["contact_fraction"] = 0.0

    t_min_global = df["t_start_ps"].min()
    t_max_global = df["t_end_ps"].max()
    print(f"[{name}] Global event window: {t_min_global:.1f} -> {t_max_global:.1f} ps")

    events = df.to_dict(orient="records")

    u = mda.Universe(topology, trajectory)

    water = u.select_atoms(WATER_SEL)
    if len(water) == 0:
        raise ValueError(f"[{name}] No water oxygens found!")

    resid_to_idx = {resid: i for i, resid in enumerate(water.resids)}

    sf_groups = {}
    for r in sf_resids:
        sel = f"protein and resid {r} and not name H*"
        atoms = u.select_atoms(sel)
        if len(atoms) == 0:
            raise ValueError(f"[{name}] No atoms found for SF residue {r}")
        sf_groups[r] = atoms
        print(f"[{name}] Residue {r}: selected {len(atoms)} heavy atoms")

    prev_time = None
    total_frames = 0

    for ts in u.trajectory:
        t_ps = float(ts.time)
        total_frames += 1

        if prev_time is None:
            prev_time = t_ps
            continue

        dt = t_ps - prev_time
        prev_time = t_ps

        if t_ps < t_min_global or t_ps > t_max_global:
            continue

        water_pos = water.positions * 0.1

        for ev in events:
            if not (ev["t_start_ps"] <= t_ps <= ev["t_end_ps"]):
                continue

            resid = int(ev["resid"])
            if resid not in resid_to_idx:
                continue

            w_idx = resid_to_idx[resid]
            w_pos = water_pos[w_idx]

            contacted_any = False
            for r in sf_resids:
                sf_atoms = sf_groups[r]
                sf_pos = sf_atoms.positions * 0.1
                dist = np.linalg.norm(sf_pos - w_pos, axis=1)
                min_dist = dist.min()

                if min_dist <= DIST_CUTOFF_NM:
                    ev[f"contact_res{r}_ps"] += dt
                    contacted_any = True

            if contacted_any:
                ev["contact_any_ps"] += dt

    print(f"[{name}] Total frames read: {total_frames}")

    df_out = pd.DataFrame(events)
    df_out["contact_fraction"] = np.where(
        df_out["duration_ps"] > 0,
        df_out["contact_any_ps"] / df_out["duration_ps"],
        0.0
    )

    outfile = os.path.join(outdir, "new_SF_contacts_all_events.csv")
    df_out.to_csv(outfile, index=False)
    print(f"[{name}] Saved: {outfile}")
    print(df_out.head())

    return outfile

def main():
    outfiles = {}
    for sysinfo in SYSTEMS:
        outfiles[sysinfo["name"]] = analyze_system(sysinfo)

    print(f"\n{'='*70}\n  DONE ??? all 3 systems processed\n{'='*70}")
    for name, path in outfiles.items():
        print(f"  {name}: {path}")

if __name__ == "__main__":
    main()
