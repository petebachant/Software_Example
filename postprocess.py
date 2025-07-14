import numpy as np
import model
import plotting_utils as pu
import analysis_utils as au
import os

# folder setup
folder = "two-dim"
if not os.path.exists(f"{folder}/analysis"):
    os.makedirs(f"{folder}/analysis")

# find all checkpoint files in data folder
checkpoints = os.listdir(f"{folder}/data")
checkpoints = sorted(checkpoints, key=lambda x: int(x))  # sort by time
if len(checkpoints) == 0:
    raise ValueError(f"No checkpoints found in {folder}/data")

# loop over checkpoints
for checkpoint in checkpoints:
    t = int(checkpoint)
    print(f"t = {t}")

    # load model state
    m = model.load(folder, t)

    # calculate state variables
    q, q_phys = au.get_q(m)
    ens = np.abs(q) ** 2
    psi, psi_phys = au.get_psi(m)
    u, u_phys = au.get_u(m, psi)
    v, v_phys = au.get_v(m, psi)

    # fit power law to kinetic energy spectrum
    ke = 1 / 2 * (np.abs(u) ** 2 + np.abs(v) ** 2)
    kmin = 8e-5
    kmax = 4e-4
    c, n = au.fit_power_law(m, ke, kmin=kmin, kmax=kmax, p0=[1e-9, -3])
    print(f"n = {n}")

    # make some plots

    pu.plot_2d_field(
        m, q_phys/m.qy/m.a, r"Potential vorticity $q / \beta L$", vmax=2, fname=f"{folder}/analysis/q_{t:015.0f}.png"
    )

    pu.plot_k_spectrum(
        m,
        ens,
        r"$|q|^2$ (s$^{-2}$)",
        fname=f"{folder}/analysis/ens_spec_{t:015.0f}.png",
    )

    pu.plot_2d_field(
        m,
        psi_phys,
        r"$\psi$ (m$^2$ s$^{-1}$)",
        fname=f"{folder}/analysis/psi_{t:015.0f}.png",
    )

    pu.plot_2d_field(
        m, u_phys, r"$u$ (m s$^{-1}$)", fname=f"{folder}/analysis/u_{t:015.0f}.png"
    )

    pu.plot_2d_field(
        m, v_phys, r"$v$ (m s$^{-1}$)", fname=f"{folder}/analysis/v_{t:015.0f}.png"
    )

    pu.plot_k_spectrum(
        m,
        ke,
        r"Energy spectrum $\mathcal{E}(k)$ (m$^3$ s$^{-2}$)",
        fname=f"{folder}/analysis/ke_spec_{checkpoint}.png",
        power_law_data=(c, n, kmin, kmax),
        xlims=(1e-5, 1e-3),
        ylims=(1e-3, 1e3),
    )