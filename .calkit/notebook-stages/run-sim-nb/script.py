# This script was automatically generated by Calkit


import os
import sys

sys.path.insert(
    0,
    os.path.abspath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../")
    ),
)


import calkit


import model
import numpy as np

folder = "two-dim-nb"
m = model.TwoDim()
m.initmean(0, 2e-11)

m.initnum(5e5, 128, 5000.0)
m.initq(1e-4 * np.random.rand(128, 128, m.nz))
m.snapshot(folder)

m.nu = 2.5e46 * 4.0**20
m.diffexp = 20
m.hypodiff = 1e-16

for i in range(100):
    m.timestep()
    m.screenlog()
    if (
        m.clock % 250000.0 < m.dt / 10
        or m.clock % 250000.0 - 250000.0 > -m.dt / 10
    ):
        m.save(folder)
        m.snapshot(folder)
m
calkit.save_notebook_stage_out(m, stage_name='run-sim-nb', out_name='m')
