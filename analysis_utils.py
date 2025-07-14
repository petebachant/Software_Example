import numpy as np
from scipy.optimize import curve_fit


def get_q(m):
    """
    Calculate the potential vorticity in spectral and physical space.

    To calculate the PV in physical space, we use the inverse FFT of the
    PV in spectral space and add the linear terms from the mean PV gradient.

    Parameters
    ----------
    m : model object
        The model object containing grid information and mean PV gradients.

    Returns
    -------
    q : ndarray
        The potential vorticity in spectral space, shape (nl, nk, nz).
    q_phys : ndarray
        The potential vorticity in physical space, shape (n, n, nz).
    """
    q_phys = m.irfft2(m.q, axes=(0, 1))
    q_phys += m.qx * (m.x - m.a / 2)
    q_phys += m.qy * (m.y - m.a / 2)
    return m.q, q_phys


def get_psi(m):
    """
    Calculate the streamfunction in spectral and physical space.

    The streamfunction is calculated by solving the linear system
    L * psi = q, where L is the linear operator defined by the model.

    Parameters
    ----------
    m : model object
        The model object containing grid information and the linear operator.

    Returns
    -------
    psi : ndarray
        The streamfunction in spectral space, shape (nl, nk, nz).
    psi_phys : ndarray
        The streamfunction in physical space, shape (n, n, nz).
    """
    psi = np.empty((m.l.size, m.k.size, m.nz), dtype=complex)
    for i in range(m.l.size):
        for j in range(m.k.size):
            psi[i, j, :] = np.linalg.solve(m.L[i, j, :, :], m.q[i, j, :])
    psi_phys = m.irfft2(psi, axes=(0, 1))
    return psi, psi_phys


def get_u(m, psi):
    """
    Calculate the zonal velocity in spectral and physical space.

    In physical space, the zonal velocity is u_phys = -d(psi_phys)/dy, so we can
    compute u in spectral space as u = i*l*psi, where l is the zonal wavenumber.

    Parameters
    ----------
    m : model object
        The model object containing grid information and the zonal wavenumber.
    psi : ndarray
        The streamfunction in spectral space, shape (nl, nk, nz).

    Returns
    -------
    u : ndarray
        The zonal velocity in spectral space, shape (nl, nk, nz).
    u_phys : ndarray
        The zonal velocity in physical space, shape (n, n, nz).
    """
    u = 1j * m.l * psi
    u_phys = m.u[0] + m.irfft2(u, axes=(0, 1))
    return u, u_phys


def get_v(m, psi):
    """
    Calculate the meridional velocity in spectral and physical space.

    In physical space, the meridional velocity is v_phys = d(psi_phys)/dx, so we can
    compute v in spectral space as v = -i*k*psi, where k is the meridional wavenumber.

    Parameters
    ----------
    m : model object
        The model object containing grid information and the meridional wavenumber.
    psi : ndarray
        The streamfunction in spectral space, shape (nl, nk, nz).

    Returns
    -------
    v : ndarray
        The meridional velocity in spectral space, shape (nl, nk, nz).
    v_phys : ndarray
        The meridional velocity in physical space, shape (n, n, nz).
    """
    v = -1j * m.k * psi
    v_phys = m.v[0] + m.irfft2(v, axes=(0, 1))
    return v, v_phys


def fit_power_law(m, ke, kmin=None, kmax=None, p0=None):
    """
    Fit a power law to the kinetic energy spectrum.

    The kinetic energy spectrum is averaged over the zonal wavenumber `l`
    and fitted to a linear function in log-log space using scipy's `curve_fit`.

    Parameters
    ----------
    m : model object
        The model object containing grid information.
    ke : ndarray
        The kinetic energy spectrum in spectral space, shape (nl, nk, nz).
    kmin : float, optional
        The minimum wavenumber to consider for the fit (default is None, uses m.k[0, 0, 0]).
    kmax : float, optional
        The maximum wavenumber to consider for the fit (default is None, uses m.k[0, -1, 0]).
    p0 : list, optional
        Initial guess for the fit parameters, expected to be a list of two elements:
        [c, n], where c is the coefficient and n is the exponent of the power law (default is [1e-9, -3]).
    Returns
    -------
    tuple
        A tuple containing the fitted coefficient `c` and exponent `n` of the power law
        in the form (c, n). The fit is performed in log-log space, so
        the returned `c` is the exponential of the intercept of the fit line,
        and `n` is the slope of the line.
    """
    ke_spec = np.mean(ke, axis=0)  # take mean over l

    def f(lk, lc, n):
        """
        Linear function in log k space.

        y = c*k**n -> f = log(y) = log(c) + n*log(k)
        """
        return lc + n * lk

    if kmin is None:
        kmin = m.k[0, 0, 0]
    if kmax is None:
        kmax = m.k[0, -1, 0]

    krange = np.where(np.logical_and(m.k[0, :, 0] >= kmin, m.k[0, :, 0] <= kmax))
    popt, pcov = curve_fit(
        f,
        np.log(m.k[0, krange, 0][0]),
        np.log(ke_spec[krange][:, 0]),
        p0=[np.log(p0[0]), p0[1]],
    )

    return np.exp(popt[0]), popt[1]
