import os, errno

import time

import matplotlib

import matplotlib.pylab as plt

import matplotlib.ticker

import numpy as np



figdir = 'figs'

try:

    os.makedirs(figdir)

except OSError as e:

    if e.errno != errno.EEXIST:

        raise

'''

Returns drag of a spherical particle given it Reynolds number based on Clift

 and Guavin 1971 model (-4%+6% deviation at Re<3e5)



Inputs:

    Re       :      Particle Reynolds number den_fluid * Particle_diameter *

                        Relative_velocity / Dynamic_viscosity



Returns:

    cd       :      Particle drag coefficient

'''





def sph_drag(Re):

    return (24.0 / Re) * (1 + 0.15 * Re ** 0.687) + 0.42 / (

            1.0 + 42500.0 / Re ** 1.16)





def calc_lag_track(dp, fl, el, partilce_density, fluid_density,

                   fluid_viscosity, x0, y0, z0, u0, v0, w0, uf0, vf0, wf0,

                   _dt, _solution_mode, drag_model):

    ''' tracks the particle until it reaches z =0.

        _solution_mode should be either 'euler' or 'analytical'

        if drag_model is set to 'sph' the spherical model will be

        used, otherwise, non-spherical model of BB2016 will be used.

        returns:

            components of particle position: x, y, z

            components of particle velocity: u, v, w

            time stamp of particle trajectory starting from the release

            moment: t

    '''

    # particle data

    flat = fl

    elon = el

    Fs = flat * elon ** 1.3  # Shape descriptor - Stoke's regime

    Fn = flat ** 2 * elon  # Shape descriptor - Newton regime

    Ks = .5 * (Fs ** (1 / 3) + Fs ** (-1 / 3))  # Stoke's drag correction

    Kn = 10 ** (.45 * (-np.log10(Fn)) ** .99)  # Newton's drag correction

    # particle diamter [m]

    diam = dp

    # time step [s]

    t = [0]

    dt = _dt

    # particle initial velocity [m/s]

    u = [u0]

    v = [v0]

    w = [w0]

    # particle initial position [m]

    x = [x0]

    y = [y0]

    z = [z0 + 1e-4]

    # fluid velocity in the whole domain [m/s]

    uf = uf0

    vf = vf0

    wf = wf0

    dens = partilce_density  # particle density

    denf = fluid_density  # air density

    visf = fluid_viscosity  # air viscosity [Pa.s]

    g = 9.8  # gravitational acceleration



    tau = dens * diam ** 2 / (18 * visf)  # Particle relaxation time (s)



    # iterate until partilce hit ground, i.e. z<0

    while z[-1] > 0:

        velr = np.sqrt((u[-1] - uf) ** 2 + (v[-1] - vf) ** 2 + (

                w[-1] - wf) ** 2)  # Particle relative velocity

        if u[-1] - uf == 0:

            theta = np.pi / 2.0

        else:

            # to check with Seb

            # direction of particle in the horizontal plane (x-y)

            theta = np.arctan((v[-1] - vf) / (u[-1] - uf))

        if np.sqrt((u[-1] - uf) ** 2 + (v[-1] - vf) ** 2) == 0:

            Beta = np.pi / 2.0

        else:

            # to check with Seb

            # direction of particle in the vertical plane (x-z or y-z)

            Beta = np.arctan((w[-1] - wf) / np.sqrt((u[-1] - uf) ** 2 +

                                                    (v[-1] - vf) ** 2))

        # to avoid possible divisions by zero, #to check with Seb

        if Beta in [0., np.pi / 2, np.pi, -np.pi / 2]:

            Beta += 1e-6

        if theta in [0., np.pi / 2, np.pi, -np.pi / 2]:

            theta += 1e-6

        Re_i = denf * velr * diam / visf  # Reynolds

        Re_S_i = Re_i * Kn / Ks  # Ganser Re



        # Drag coefficient,

        if Re_S_i > 3e5:  # to check with Seb

            # drag crisis

            Cd_i = 0.2

        else:

            if drag_model == 'sph':

                # to check with Seb

                Cd_i = sph_drag(Re_i)

            else:

                Cd_i = Kn * (24 * (

                        1 + .125 * Re_S_i ** (2 / 3)) / Re_S_i + .46 / (

                                     1 + 5330 / Re_S_i))  # Drag coef (eq. 34)



        ## If particle within region of reduced drag

        # if part.dis < P.adv.drag:

        #  part.Cd(i)   = 0;



        Fd = Cd_i * Re_i / (24.0 * tau)  # Total Drag force

        Fd_u = abs(

            Fd * np.cos(Beta) * np.cos(theta))  # Drag force in x direction

        G_u = 0.  # Other forces in x direction

        Fd_v = abs(

            Fd * np.cos(Beta) * np.sin(theta))  # Drag force in y direction

        G_v = 0.  # Other forces in y direction

        Fd_w = abs(Fd * np.sin(Beta))  # Drag force in z direction

        G_w = (1 - denf / dens) * -g  # Gravity force



        if _solution_mode == 'analytical':

            # to check with Seb

            u.append(

                uf + np.exp(-Fd_u * dt) * (u[-1] - uf) - G_u * (1 / Fd_u) * (

                        np.exp(-dt * Fd_u) - 1))

            x.append(x[-1] + (G_u * (1 / Fd_u) + uf) * dt + (1 / Fd_u) * (

                    1 - np.exp(-dt * Fd_u)) * (u[-1] - uf - G_u / Fd_u))



            v.append(

                vf + np.exp(-Fd_v * dt) * (v[-1] - vf) - G_v * (1 / Fd_v) * (

                        np.exp(-dt * Fd_v) - 1))

            y.append(y[-1] + (G_v * (1 / Fd_v) + vf) * dt + (1 / Fd_v) * (

                    1 - np.exp(-dt * Fd_v)) * (v[-1] - vf - G_v / Fd_v))



            w.append(

                wf + np.exp(-Fd_w * dt) * (w[-1] - wf) - G_w * (1 / Fd_w) * (

                        np.exp(-dt * Fd_w) - 1))

            z.append(z[-1] + (G_w * (1 / Fd_w) + wf) * dt + (1 / Fd_w) * (

                    1 - np.exp(-dt * Fd_w)) * (w[-1] - wf - G_w / Fd_w))



        elif _solution_mode == 'euler':

            # Euler semi-implicit

            u.append(((G_u + Fd_u * uf) * dt + u[-1]) / (1 + Fd_u * dt))

            x.append(x[-1] + .5 * dt * (u[-1] + u[-2]))



            v.append(((G_v + Fd_v * vf) * dt + v[-1]) / (1 + Fd_v * dt))

            y.append(y[-1] + .5 * dt * (v[-1] + v[-2]))



            w.append(((G_w + Fd_w * wf) * dt + w[-1]) / (1 + Fd_w * dt))

            z.append(z[-1] + .5 * dt * (w[-1] + w[-2]));

        else:

            raise ('solution mode not implemented')



        t.append(t[-1] + dt)

    return x, y, z, u, v, w, t





# start of the main code

# initialization of run paramters, all in SI units

dp = 1e-3  # particle diameter

fl = 1.0  # flatness

el = 1.0  # elongation

partilce_density = 1000.

fluid_density = 1.1615996268298834

fluid_viscosity = 1.8537151855672458e-05

# initial position of the particle

x0 = 0

y0 = 0

z0 = 10.627699

# initial release velocity of the particle

u0 = 0.

v0 = 0.0001

w0 = 0

# fluid velocity, constant-value throughout the whole domain

uf0 = 0

vf0 = 0.5

wf0 = 0

# creating a dict for storing position, velocity and time info of both

# simulation methods

(x, y, z, u, v, w, t) = (None, None, None, None, None, None, None)

res = {'euler': (x, y, z, u, v, w, t),

       'analytical': (x, y, z, u, v, w, t)}

# assigning dt

dt = 0.001

# running for both methods

for key in res:

    res[key] = calc_lag_track(dp, fl, el, partilce_density, fluid_density,

                              fluid_viscosity, x0, y0, z0, u0, v0, w0, uf0,

                              vf0, wf0, dt, key, 'sph')



# plotting results

fig = plt.figure(figsize=(7, 5))

ax = fig.add_subplot(1, 1, 1)

ax.plot(res['euler'][4], res['euler'][2], 'r*', label='euler')

ax.plot(res['analytical'][4], res['analytical'][2], 'b-', label='analytical')

ax.legend()

plt.tight_layout()

plt.show(fig)

