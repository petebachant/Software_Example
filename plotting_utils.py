import numpy as np
import matplotlib.pyplot as plt

plt.style.use("plots.mplstyle")


def plot_2d_field(m, field, label, vmax=None, ilev=0, fname="2d_field.png"):
    """
    Helper function to plot a 2D (physical) field.

    Parameters
    ----------
    m : model object
        The model object containing grid information.
    field : ndarray
        The field to be plotted, expected to have shape (n, n, nz).
    label : str
        The label for the colorbar.
    vmax : float, optional
        The maximum value for the color scale (default is None, which uses the 
        max of the absolute value of the field).
    ilev : int, optional
        The index of the vertical level to plot (default is 0).
    fname : str, optional
        The filename to save the plot (default is "2d_field.png").
    """
    if vmax is None:
        vmax = np.max(np.abs(field[:, :, ilev]))
    fig, ax = plt.subplots(1, figsize=(3.2, 2.4))
    im = ax.pcolormesh(
        m.x[0, :, ilev] / 1e3,
        m.y[:, 0, ilev] / 1e3,
        field[:, :, ilev],
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=+vmax,
    )
    cb = plt.colorbar(im, ax=ax, label=label, shrink=0.6, pad=0.1)
    cb.ax.ticklabel_format(style="sci", useMathText=True, scilimits=(0, 0))
    ax.set_xlabel(r"$x$ (km)")
    ax.set_ylabel(r"$y$ (km)")
    ax.axis("equal")
    ax.set_xlim(0, m.a / 1e3)
    ax.set_ylim(0, m.a / 1e3)
    ax.set_xticks(np.array([0, m.a / 2, m.a]) / 1e3)
    ax.set_yticks(np.array([0, m.a / 2, m.a]) / 1e3)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    plt.savefig(fname)
    print(f"Saved '{fname}'")
    plt.close()


def plot_k_spectrum(
    m,
    field,
    label,
    ilev=0,
    fname="k_spectrum.png",
    power_law_data=None,
    xlims=None,
    ylims=None,
):
    """
    Helper function to plot the (k) spectrum of a field.

    This function computes the mean of the field over l and plots the resulting
    spectrum against the wavenumber `k`. It can also include a power law line
    if provided.

    Parameters
    ----------
    m : model object
        The model object containing grid information.
    field : ndarray
        The (spectral) field to be analyzed, expected to have shape (nl, nk, nz).
    label : str
        The label for the y-axis of the spectrum plot.
    ilev : int, optional
        The index of the vertical level to analyze (default is 0).
    fname : str, optional
        The filename to save the plot (default is "k_spectrum.png").
    power_law_data : tuple, optional
        A tuple containing (c, n, kmin, kmax) for plotting a power law line.
    xlims : tuple, optional
        The limits for the x-axis (default is None).
    ylims : tuple, optional
        The limits for the y-axis (default is None).
    """
    spectrum = np.mean(field, axis=0)  # take mean over l
    fig, ax = plt.subplots(1)
    ax.plot(m.k[0, :, ilev], spectrum, "C0-")
    if power_law_data is not None:
        c, n, kmin, kmax = power_law_data
        ax.plot(
            [kmin, kmax], [c * kmin**n, c * kmax**n], "k--", label=f"$k^{{{n:.1f}}}$"
        )
        ax.legend()
    ax.set_xscale("log")
    ax.set_yscale("log")
    if xlims is not None:
        ax.set_xlim(xlims)
    if ylims is not None:
        ax.set_ylim(ylims)
    ax.set_xlabel(r"$k$ (m$^{-1}$)")
    ax.set_ylabel(label)
    plt.savefig(fname)
    print(f"Saved '{fname}'")
    plt.close()
