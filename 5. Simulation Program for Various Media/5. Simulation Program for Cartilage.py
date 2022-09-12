import numpy as np
from matplotlib import pyplot as plt
from math import pi, exp, cos, sin, sqrt, atan2

ke = 200
ex = np.zeros(ke)
dx = np.zeros(ke)
ix = np.zeros(ke)
sx = np.zeros(ke)
hy = np.zeros(ke)
# creating the geometry
ddx = 0.01  # Cell size
dt = ddx / 6e8  # Time step size
number_of_frequencies = 9
# required input frequency
freq_in = np.array((10e6, 100e6, 200e6, 300e6, 400e6,
                    500e6, 600e6, 700e6, 800e6))
t0 = 50
spread = 10


# information about the dielectric mediums
epsz = 8.854e-12
epsr = 43            # at 800MHz
sigma = 3.69E-1         # at 10MHz
tau = 0.001 * 1e-6
chi = 0
k_start = 100
# initializing boundary Conditions
boundary_low = [0, 0]
boundary_high = [0*i for i in range(2*int(sqrt(epsr)))]

# iterative variables
gax = np.ones(ke)
gbx = np.zeros(ke)
gcx = np.zeros(ke)

gax[k_start:] = 1 / (epsr + (sigma * dt / epsz) + chi * dt / tau)
gbx[k_start:] = sigma * dt / epsz
gcx[k_start:] = chi * dt / tau

del_exp = exp(-dt / tau)

# variables for Fourier transform:
arg = 2 * np.pi * freq_in * dt
real_pt = np.zeros((number_of_frequencies, ke))
imag_pt = np.zeros((number_of_frequencies, ke))
real_in = np.zeros(number_of_frequencies)
imag_in = np.zeros(number_of_frequencies)
amp_in = np.zeros(number_of_frequencies)
phase_in = np.zeros(number_of_frequencies)
amp = np.zeros((number_of_frequencies, ke))
phase = np.zeros((number_of_frequencies, ke))

# REF = np.zeros((number_of_frequencies,ke))

TRN = np.zeros((number_of_frequencies, ke))

nsteps = 1000

# Dictionary for the timestep
plotting_points = [{'num_steps': 350, 'ex': None, 'scaling_factor': 1,
                    'gb_scaling_factor': 1,
                    'y_ticks': (np.arange(-1, 1, step=0.5)),
                    'y_min': -1.3, 'y_max': 1.2, 'y_text_loc': 0.3,
                    'label': '(a)',
                    'label_loc': 1},
                   {'num_steps': 1000, 'ex': None, 'scaling_factor': 1,
                    'gb_scaling_factor': 1,
                    'y_ticks': (np.arange(-1, 1, step=0.5)),
                    'y_min': -1.3, 'y_max': 1.2,
                    'y_text_loc': 0.3, 'label': '(b)', 'label_loc': 1}]


# Main FDTD Loop
for time_step in range(1, nsteps + 1):

    # Calculate Dx in time domain
    for k in range(1, ke):
        dx[k] = dx[k] + 0.5 * (hy[k - 1] - hy[k])

    # Put a pulse wave in time domain
    pulse = exp(-0.5 * ((t0 - time_step) / spread) ** 2)
    dx[5] = pulse + dx[5]

    # Calculate the Ex field from Dx in time domain
    for k in range(1, ke):
        ex[k] = gax[k] * (dx[k] - ix[k] - del_exp * sx[k])
        ix[k] = ix[k] + gbx[k] * ex[k]
        sx[k] = del_exp * sx[k] + gcx[k] * ex[k]

    # Fourier Transform of the input pulse
    for k in range(ke):
        for m in range(number_of_frequencies):
            real_pt[m, k] = real_pt[m, k] + cos(arg[m] * time_step) * ex[k]
            imag_pt[m, k] = imag_pt[m, k] - sin(arg[m] * time_step) * ex[k]

    # Fourier Transform of the input pulse
    if time_step < (190):
        for m in range(number_of_frequencies):
            real_in[m] = real_in[m] + cos(arg[m] * time_step) * ex[10]
            imag_in[m] = imag_in[m] - sin(arg[m] * time_step) * ex[10]

    # Boundary Conditions
    ex[0] = boundary_low.pop(0)
    boundary_low.append(ex[1])
    ex[ke - 1] = boundary_high.pop(0)
    boundary_high.append(ex[ke - 2])

    #  Hy field calculations
    for k in range(ke - 1):
        hy[k] = hy[k] + 0.5 * (ex[k] - ex[k + 1])

# Saving the data plotting_points
    for plotting_point in plotting_points:
        if time_step == plotting_point['num_steps']:
            plotting_point['ex'] = np.copy(ex)

# Calculate the amplitude and phase at each frequency
            for m in range(number_of_frequencies):
                amp_in[m] = sqrt(imag_in[m] ** 2 + real_in[m] ** 2)
                phase_in[m] = atan2(imag_in[m], real_in[m])
                for k in range(ke):
                    amp[m, k] = (1 / amp_in[m]) * sqrt((real_pt[m, k]) ** 2
                                                       + imag_pt[m, k] ** 2)
                    phase[m, k] = atan2(imag_pt[m, k],
                                        real_pt[m, k]) - phase_in[m]
                for k in range(100, ke):
                    TRN[m, k] = ((amp[m, k])**2)*sqrt(epsr)

# Dictionary for transmission
fig = plt.figure(figsize=(8, 16))
plotting_trns = [{'freq': freq_in[0], 'TRN': TRN[0, ],
                  'label': '', 'x_label': ''},
                 {'freq': freq_in[1], 'TRN': TRN[1, ],
                  'label': '', 'x_label': ''},
                 {'freq': freq_in[2], 'TRN': TRN[2, ],
                  'label': '', 'x_label': ''},
                 {'freq': freq_in[3], 'TRN': TRN[3, ],
                  'label': '', 'x_label': ''},
                 {'freq': freq_in[4], 'TRN': TRN[4, ],
                  'label': '', 'x_label': ''},
                 {'freq': freq_in[5], 'TRN': TRN[5, ],
                  'label': '', 'x_label': ''},
                 {'freq': freq_in[6], 'TRN': TRN[6, ],
                  'label': '', 'x_label': ''},
                 {'freq': freq_in[7], 'TRN':TRN[7, ],
                  'label': '', 'x_label': ''},
                 {'freq': freq_in[8], 'TRN': TRN[8, ],
                  'label': '', 'x_label': 'FDTD Cells'}]

# plotting transmission
def plot_trn(data1, freq, label, x_label):
    """Plot of amplitude at one frequency"""
    plt.plot(data1, color='k', linewidth=1)
    plt.ylabel('TRNS')
    plt.xticks(np.arange(100, 199, step=20))
    plt.xlim(100, 198)
    # plt.yticks(np.arange(0, 2.1))
    # plt.ylim(-0.2, 2.0)
    plt.text(170, 0.015, 'Freq. at {} MHz'.format(int(round(freq / 1e6))),
             horizontalalignment='center')
    # plt.plot(gb * 1 / gb[120], 'k--',linewidth=0.75)
    # gb is just for scaling
    plt.text(-25, 0.6, label, horizontalalignment='center')
    plt.xlabel(x_label)
    return


for subplot_num1, plotting_trn in enumerate(plotting_trns):
    ax = fig.add_subplot(9, 1, subplot_num1+1)
    plot_trn(plotting_trn['TRN'], plotting_trn['freq'],
             plotting_trn['label'], plotting_trn['x_label'])

plt.subplots_adjust(bottom=0.1, hspace=0.45)
# saving the figure
fig.savefig('cartilagetran.png')

plt.rcParams['font.size'] = 12
fig = plt.figure(figsize=(8, 16))


def plot_e_field(data, gb, timestep, scaling_factor, gb_scaling_factor,
                 y_ticks, y_min, y_max, y_text_loc, label, label_loc):
    """Plot of E field at a single time step"""
    plt.plot(data * scaling_factor, color='k', linewidth=1)
    plt.ylabel('E$_x$(V/m)', fontsize='12')
    plt.xticks(np.arange(0, 199, step=20))
    plt.xlim(0, 199)
    plt.yticks(y_ticks)
    plt.ylim(y_min, y_max)
    plt.text(150, y_text_loc, 'Time Domain, T = {}'.format(timestep),
             horizontalalignment='center')
    plt.plot(gb * gb_scaling_factor / gb[120], 'k--', linewidth=0.75)
    # The math on gb is just for scaling
    plt.text(-25, label_loc, label, horizontalalignment='center')
    return


# Plotting e field
for subplot_num, plotting_point in enumerate(plotting_points):
    ax = fig.add_subplot(11, 1, subplot_num + 1)
    plot_e_field(plotting_point['ex'], gbx, plotting_point['num_steps'],
                 plotting_point['scaling_factor'],
                 plotting_point['gb_scaling_factor'],
                 plotting_point['y_ticks'],
                 plotting_point['y_min'],
                 plotting_point['y_max'], plotting_point['y_text_loc'],
                 plotting_point['label'],
                 plotting_point['label_loc'])

# Dictionary  for the amplitudes
plotting_freqs = [{'freq': freq_in[0], 'amp': amp[0],
                   'label': '(c)', 'x_label': ''},
                  {'freq': freq_in[1], 'amp': amp[1],
                   'label': '', 'x_label': ''},
                  {'freq': freq_in[2], 'amp': amp[2],
                   'label': '', 'x_label': ''},
                  {'freq': freq_in[3], 'amp': amp[3],
                   'label': '', 'x_label': ''},
                  {'freq': freq_in[4], 'amp': amp[4],
                   'label': '', 'x_label': ''},
                  {'freq': freq_in[5], 'amp': amp[5],
                   'label': '', 'x_label': ''},
                  {'freq': freq_in[6], 'amp': amp[6],
                   'label': '', 'x_label': ''},
                  {'freq': freq_in[7], 'amp': amp[7],
                   'label': '', 'x_label': ''},
                  {'freq': freq_in[8], 'amp': amp[8],
                   'label': '', 'x_label': 'FDTD Cells'}]


def plot_amp(data, gb, freq, label, x_label):
    """Plot of amplitude at one frequency"""
    plt.plot(data, color='k', linewidth=1)
    plt.ylabel('Amp')
    plt.xticks(np.arange(0, 199, step=20))
    plt.xlim(0, 198)
    plt.yticks(np.arange(0, 2.1, step=1))
    plt.ylim(-0.2, 2.0)
    plt.text(150, 1.2, 'Freq. Domain at {} MHz'.format(int(round(freq / 1e6))),
             horizontalalignment='center')
    plt.plot(gb * 1 / gb[120], 'k--',
             linewidth=0.75)
    # The math on gb is just for scaling
    plt.text(-25, 0.6, label, horizontalalignment='center')
    plt.xlabel(x_label)
    return


for subplot_num, plotting_freq in enumerate(plotting_freqs):
    ax = fig.add_subplot(11, 1, subplot_num+3)
    plot_amp(plotting_freq['amp'], gbx, plotting_freq['freq'],
             plotting_freq['label'], plotting_freq['x_label'])

plt.subplots_adjust(bottom=0.1, hspace=0.45)
plt.show()
fig.savefig('cartilagesimu.png')
