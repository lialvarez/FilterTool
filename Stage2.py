import numpy as np
from scipy import signal


class Stage2(object):
    """description of class"""

    def __init__(self, id):
        self.pole_pairs = []
        self.zero_pairs = []
        # Default gain
        self.gain = 1

        self.id = id + 1
        self.label = self.__generate_label()

        self.canvas = None
        self.figure = None

    def change_id(self, id):
        self.id = id + 1
        self.label = self.__generate_label()

    def clear(self):
        for pair in self.pole_pairs:
            pair.used = False
        for pair in self.zero_pairs:
            pair.used = False

    def remove_pole_pair(self, pair):
        for my_pair in self.pole_pairs:
            if pair == my_pair:
                self.pole_pairs.remove(pair)
                pair.used = False

    def remove_zero_pair(self, pair):
        for my_pair in self.zero_pairs:
            if pair == my_pair:
                self.zero_pairs.remove(pair)
                pair.used = False

    def add_pole_pair(self, pair):
        self.pole_pairs.append(pair)
        pair.used = True

    def add_zero_pair(self, pair):
        self.zero_pairs.append(pair)
        pair.used = True

    def set_gain(self, gain):
        self.gain = gain

    def link_canvas(self, canvas):
        self.canvas = canvas

    def link_figure(self, figure):
        self.figure = figure

    def plot_s_plane_to_figure(self, figure):
        # TODO: pasarle el axes como argumento.
        # remove axes
        for axes in figure.axes:
            figure.delaxes(axes)
        # create axes
        axes = figure.add_subplot(111, projection='polar')
        axes.set_title('{0}: Poles \& Zeros'.format(self.label), pad=10)
        marker_poles = 'x'
        marker_zeros = 'o'

        r_poles = []
        theta_poles = []
        r_zeros = []
        theta_zeros = []

        for pole_pair in self.pole_pairs:
            pole_theta = []
            pole_r = []
            for complex in pole_pair.points:
                r = np.real(complex)
                i = np.imag(complex)
                rho = np.sqrt(r ** 2 + i ** 2)
                theta = np.arctan2(i, r)
                pole_r.append(rho)
                pole_theta.append(theta)
            axes.scatter(pole_theta, pole_r, label='Pole {0}'.format(pole_pair.id), marker=marker_poles)
        for zero_pair in self.zero_pairs:
            zero_theta = []
            zero_r = []
            for complex in zero_pair.points:
                r = np.real(complex)
                i = np.imag(complex)
                rho = np.sqrt(r ** 2 + i ** 2)
                theta = np.arctan2(i, r)
                zero_r.append(rho)
                zero_theta.append(theta)
            axes.scatter(zero_theta, zero_r, label='Zero {0}'.format(zero_pair.id), marker=marker_zeros)
        axes.legend()

    def plot_transfer_to_figure(self, figure):
        n_points = 1000
        z = []
        p = []
        # get a list of all the zeros:
        for point_pair in self.zero_pairs:
            for point in point_pair.points:
                z.append(point)

        # get a list of all the poles
        for point_pair in self.pole_pairs:
            for point in point_pair.points:
                p.append(point)

        sys = signal.TransferFunction(signal.ZerosPolesGain(z, p, self.gain))
        w, mag, phase = signal.bode(sys)

        # create and plot
        title = 'Stage {0} - Frequency Response[dB]'.format(self.id)
        legend = 'Frequency Response[dB]'
        # remove axes
        for axes in figure.axes:
            figure.delaxes(axes)
        # create axes
        axes = figure.add_subplot(111, projection='rectilinear')
        axes.semilogx(w, mag, label=legend)
        axes.set_title(title, pad=10)
        axes.legend()
        axes.set_xlabel(r'$\omega [rad/seg]$')
        axes.set_ylabel(r'Mag')

    def __generate_label(self):
        label = 'Stage {0}'.format(self.id)
        return label
