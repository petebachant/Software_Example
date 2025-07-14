import numpy as np
import model

# 2D turbulence
folder = "two-dim"
m = model.TwoDim()

# mean PV gradients
qx = 0
qy = 2e-11  # this is beta (m s-2)
m.initmean(qx, qy)

# initialize numerics
domain_size = 5e5  # m
n = 128  # grid points
dt = 5000.0  # time step (s)
m.initnum(domain_size, n, dt)

# initial condition
amp = 1e-4  # amplitude of initial noise
q0 = amp * np.random.rand(n, n, m.nz)
m.initq(q0)
m.snapshot(folder)  # save an initial snapshot

# load a previous model state
# m = model.load(folder, 50000000)

# hyper- and hypo-viscosity parameters
m.nu = 2.5e46 * 4.0**20  # m20 s-1
m.diffexp = 20
m.hypodiff = 1e-16  # m-2 s-1

# run
nsteps = 20000  # number of time steps
tsave = 250000.0  # save every tsave seconds
for i in range(nsteps):
    m.timestep()
    m.screenlog()
    if m.clock % tsave < m.dt / 10 or m.clock % tsave - tsave > -m.dt / 10:
        m.save(folder)
        m.snapshot(folder)