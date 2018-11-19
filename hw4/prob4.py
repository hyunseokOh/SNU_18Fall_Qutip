import matplotlib.pyplot as plt
import time
import numpy as np
from numpy import pi

from qutip import *
from qutip.control import *

T = 3
times = np.linspace(0, T, 100)
theta, phi = np.random.rand(2)

# target unitary transformation (random single qubit rotation)
U = rz(phi) * rx(theta);
print(U)
R = 200
H_ops = [sigmax(), sigmay(), sigmaz()]

H_labels = [r'$u_{x}$',r'$u_{y}$',r'$u_{z}$',]
H0 = 0 * pi * sigmaz()

from qutip.control.grape import plot_grape_control_fields, _overlap
from qutip.control.cy_grape import cy_overlap
from qutip.control.grape import cy_grape_unitary, grape_unitary_adaptive

from scipy.interpolate import interp1d
from qutip.ui.progressbar import TextProgressBar

u0 = np.array([np.random.rand(len(times)) * 2 * pi * 0.005 for _ in range(len(H_ops))])
u0 = [np.convolve(np.ones(10)/10, u0[idx,:], mode='same') for idx in range(len(H_ops))]

result = cy_grape_unitary(U, H0, H_ops, R, times, u_start=u0, eps=2*pi/T, phase_sensitive=False,progress_bar=TextProgressBar())


plot_grape_control_fields(times, result.u[:,:,:] / (2 * pi), H_labels, uniform_axes=True);


# target unitary
print(U)

# unitary from grape pulse
print(result.U_f)

# target / result overlap
_overlap(U, result.U_f).real, abs(_overlap(U, result.U_f))**2


c_ops = []
U_f_numerical = propagator(result.H_t, times[-1], c_ops, args={})
print(U_f_numerical)
_overlap(U, U_f_numerical)


psi0 = basis(2, 0)
e_ops = [sigmax(), sigmay(), sigmaz()]
me_result = mesolve(result.H_t, psi0, times, c_ops, e_ops)

b = Bloch()

b.add_points(me_result.expect)

b.add_states(psi0)
b.add_states(U * psi0)
b.render()

op_basis = [[qeye(2), sigmax(), sigmay(), sigmaz()]]
op_label = [["i", "x", "y", "z"]]

fig = plt.figure(figsize=(8,6))
U_ideal = spre(U) * spost(U.dag())
chi = qpt(U_ideal, op_basis)
fig = qpt_plot_combined(chi, op_label, fig=fig, threshold=0.001)

plt.show()
