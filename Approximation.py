from PyQt5 import QtGui, QtWidgets, QtCore
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
from Template import Template
from scipy import signal
import math
import numpy as np
import cmath
# TODO: sacar solo debug
import matplotlib.pyplot as mpl
from Etapa import *
from PolePair2 import PolePair2
from ZeroPair2 import ZeroPair2
from Stage2 import Stage2


class GroupBoxPar(QtWidgets.QGroupBox):
    def __init__(self, tab):
        QtWidgets.QGroupBox.__init__(self, tab)
        self.radioButtons = []
        self.par = None


class Approximation(object):
    """description of class"""

    def __init__(self, program, template, **kwargs):
        self.program = program
        self.template = template
        self.restriction = kwargs['restriction']
        self.min_order = kwargs['min_order']
        self.max_order = kwargs['max_order']
        self.max_q = kwargs['max_q']
        self.custom_order = kwargs['custom_order']
        self.approx_type = kwargs['approx_type']
        self.den = []
        self.num = []
        self.den_norm = []
        self.num_norm = []
        self.points = []
        self.stages = []
        self.stages_2 = []
        self.stages_canvases = []
        self.index = len(program.approximations)
        self.mainStageFigure = None
        self.mainStageFigureCanvas = None
        self.ax = None

        if self.approx_type == 'butterworth':
            self.approx_type_pretty = 'Butterworth'
        elif self.approx_type == 'bessel':
            self.approx_type_pretty = 'Bessel'
        elif self.approx_type == 'cheby_1':
            self.approx_type_pretty = 'Chebyshev I'
        elif self.approx_type == 'cheby_2':
            self.approx_type_pretty = 'Chebyshev II'
        elif self.approx_type == 'legendre':
            self.approx_type_pretty = 'Legendre'
        elif self.approx_type == 'gauss':
            self.approx_type_pretty = 'Gauss'
        elif self.approx_type == 'cauer':
            self.approx_type_pretty = 'Cauer'

        self.__compute_parameters()
        self.compute_approximation()

        # Generar las etapas automaticamente
        self.__init_stages()

        # # Esto es lo de mati
        # pointsQuantity = len(self.zeros) + len(self.poles)
        # cm = mpl.get_cmap('gist_rainbow', pointsQuantity)
        # self.colors = [cm(i) for i in range(pointsQuantity)]
        # self.initStages()
        # self.initPoints()
        # self.updateStages()
        # self.initFigures()

        # for pointList in [self.zeros, self.poles]:
        #     for pair in pointList:
        #         for radioButton in pair.groupBox.radioButtons:
        #             radioButton.clicked.connect(self.updateStages)

    def __compute_parameters(self):
        if self.template.filter_type == 'Low-pass':
            self.wp = self.template.omega_p1
            self.ws = self.template.omega_s1
            self.filter_t = 'lowpass'
        elif self.template.filter_type == 'High-pass':
            self.wp = self.template.omega_p1
            self.ws = self.template.omega_s1
            self.filter_t = 'highpass'
        elif self.template.filter_type == 'Band-pass':
            self.filter_t = 'bandpass'
            self.wp = [self.template.omega_p1, self.template.omega_p2]
            self.ws = [self.template.omega_s1, self.template.omega_s2]
        elif self.template.filter_type == 'Band-reject':
            self.filter_t = 'bandstop'
            self.wp = [self.template.omega_p1, self.template.omega_p2]
            self.ws = [self.template.omega_s1, self.template.omega_s2]

        if self.restriction == 'max_q':
            self.__parameters_max_q()
        elif self.restriction == 'custom_order':
            self.order = self.__parameters_custom_order()
        elif self.restriction == 'min_max_order':
            self.order = self.__parameters_min_max_order()

        return

    def __parameters_custom_order(self):
        order = self.custom_order
        if self.approx_type == 'butterworth':
            N, self.wn = signal.buttord(self.wp, self.ws, self.template.att_p, self.template.att_s, analog=True)
            N_norm, self.wn_N = signal.buttord(self.template.omega_pN, self.template.omega_sN, self.template.att_p,
                                               self.template.att_s, analog=True)
        elif self.approx_type == 'bessel':
            pass
        elif self.approx_type == 'cheby_1':
            N, self.wn = signal.cheb1ord(self.wp, self.ws, self.template.att_p, self.template.att_s, analog=True)
            N_norm, self.wn_N = signal.cheb2ord(self.template.omega_pN, self.template.omega_sN, self.template.att_p,
                                                self.template.att_s, analog=True)
        elif self.approx_type == 'cheby_2':
            N, self.wn = signal.cheb2ord(self.wp, self.ws, self.template.att_p, self.template.att_s, analog=True)
            N_norm, self.wn_N = signal.cheb2ord(self.template.omega_pN, self.template.omega_sN, self.template.att_p,
                                                self.template.att_s, analog=True)
        elif self.approx_type == 'legendre':
            pass
        elif self.approx_type == 'gauss':
            pass
        elif self.approx_type == 'cauer':
            N, self.wn = signal.ellipord(self.wp, self.ws, self.template.att_p, self.template.att_s, analog=True)
            N_norm, self.wn_N = signal.ellipord(self.template.omega_pN, self.template.omega_sN, self.template.att_p,
                                                self.template.att_s, analog=True)

        return order

    def __order_min_max_butter(self):
        N, self.wn = signal.buttord(self.wp, self.ws, self.template.att_p, self.template.att_s, analog=True)
        N_norm, self.wn_N = signal.buttord(self.template.omega_pN, self.template.omega_sN, self.template.att_p,
                                           self.template.att_s, analog=True)

        if N < self.min_order:
            order = self.min_order
        elif N > self.max_order:
            order = self.max_order
        else:
            order = N

        return order

    def __order_min_max_cheby_1(self):
        N, self.wn = signal.cheb1ord(self.wp, self.ws, self.template.att_p, self.template.att_s, analog=True)
        N_norm, self.wn_N = signal.cheb1ord(self.template.omega_pN, self.template.omega_sN, self.template.att_p,
                                            self.template.att_s, analog=True)

        if N < self.min_order:
            order = self.min_order
        elif N > self.max_order:
            order = self.max_order
        else:
            order = N

        return order

    def __order_min_max_cheby_2(self):
        N, self.wn = signal.cheb2ord(self.wp, self.ws, self.template.att_p, self.template.att_s, analog=True)
        N_norm, self.wn_N = signal.cheb2ord(self.template.omega_pN, self.template.omega_sN, self.template.att_p,
                                            self.template.att_s, analog=True)

        if N < self.min_order:
            order = self.min_order
        elif N > self.max_order:
            order = self.max_order
        else:
            order = N

        return order

    def __order_min_max_cauer(self):
        N, self.wn = signal.ellipord(self.wp, self.ws, self.template.att_p, self.template.att_s, analog=True)
        N_norm, self.wn_N = signal.ellipord(self.template.omega_pN, self.template.omega_sN, self.template.att_p,
                                            self.template.att_s, analog=True)

        if N < self.min_order:
            order = self.min_order
        elif N > self.max_order:
            order = self.max_order
        else:
            order = N

        return order

    def __order_min_max_bessel(self):
        order = self.min_order - 1
        founded = False
        while order <= self.max_order and not founded:
            order += 1
            num, den = signal.bessel(order, 1, 'low', analog=True, output='ba')
            # switched num and den to get attenuation 'points'
            sys = signal.TransferFunction(den, num)
            w, mag, phase = signal.bode(sys)
            mag_db = np.log10(mag)
            i = 0
            compatible = True
            for i in range(0, len(w)):
                if w[i] < 1:
                    if mag_db[i] > self.template.att_p:
                        compatible = False
                        break
                elif w[i] > self.template.omega_sN:
                    if mag_db[i] < self.template.att_s:
                        compatible = False
                        break
            if compatible:
                founded = True

        if order < self.min_order:
            order = self.min_order
        elif order > self.max_order:
            order = self.max_order

        return order

    def __parameters_max_q(self):
        if self.approx_type == 'butterworth':
            order_q = self.__order_max_q_butter()
            N, self.wn = signal.buttord(self.wp, self.ws, self.template.att_p, self.template.att_s, analog=True)
            N_norm, self.wn_N = signal.buttord(self.template.omega_pN, self.template.omega_sN, self.template.att_p,
                                               self.template.att_s, analog=True)
        elif self.approx_type == 'bessel':
            order_q = self.__order_max_q_bessel()
            N = order_q
        elif self.approx_type == 'cheby_1':
            order_q = self.__order_max_q_cheby_1()
            N, self.wn = signal.cheb1ord(self.wp, self.ws, self.template.att_p, self.template.att_s, analog=True)
            N_norm, self.wn_N = signal.cheb2ord(self.template.omega_pN, self.template.omega_sN, self.template.att_p,
                                                self.template.att_s, analog=True)
        elif self.approx_type == 'cheby_2':
            order_q = self.__order_max_q_cheby_2()
            N, self.wn = signal.cheb2ord(self.wp, self.ws, self.template.att_p, self.template.att_s, analog=True)
            N_norm, self.wn_N = signal.cheb1ord(self.template.omega_pN, self.template.omega_sN, self.template.att_p,
                                                self.template.att_s, analog=True)
        elif self.approx_type == 'legendre':
            pass
        elif self.approx_type == 'gauss':
            pass
        elif self.approx_type == 'cauer':
            order_q = self.__order_max_q_cauer()
            N, self.wn = signal.ellipord(self.wp, self.ws, self.template.att_p, self.template.att_s, analog=True)
            N_norm, self.wn_N = signal.ellipord(self.template.omega_pN, self.template.omega_sN, self.template.att_p,
                                                self.template.att_s, analog=True)

        if order_q > N:
            self.order = N
        elif order_q <= N:
            self.order = order_q
        return

    def __parameters_min_max_order(self):
        if self.approx_type == 'butterworth':
            order = self.__order_min_max_butter()
        elif self.approx_type == 'cheby_1':
            order = self.__order_min_max_cheby_1()
        elif self.approx_type == 'cheby_2':
            order = self.__order_min_max_cheby_2()
        elif self.approx_type == 'cauer':
            order = self.__order_min_max_cauer()
        elif self.approx_type == 'bessel':
            order = self.__order_min_max_bessel()
        return order

    def __order_max_q_butter(self):
        found = False
        order = 0
        wn = self.wp
        q = []
        while not found:
            order = order + 1
            zeros, poles, gain = signal.butter(order, wn, self.filter_t, analog=True, output='zpk')
            for p in poles:
                r, phi = cmath.polar(p)
                q.append(r / (2 * abs(r * math.cos(phi))))
            max_q = max(q)
            if (q and max_q >= self.max_q) or order > 15:
                found = True
        return order

    def __order_max_q_bessel(self):
        found = False
        order = 0
        wn = self.wp
        q = []
        while not found:
            order = order + 1
            zeros, poles, gain = signal.bessel(order, wn, self.filter_t, analog=True, output='zpk')
            for p in poles:
                r, phi = cmath.polar(p)
                q.append(r / (2 * abs(r * math.cos(phi))))
            max_q = max(q)
            if (q and max_q >= self.max_q) or order > 15:
                found = True
        return order

    def __order_max_q_cheby_1(self):
        found = False
        order = 0
        wn = self.wp
        q = []
        while not found:
            order = order + 1
            zeros, poles, gain = signal.cheby1(order, self.template.att_p, wn, self.filter_t, analog=True, output='zpk')
            for p in poles:
                r, phi = cmath.polar(p)
                q.append(r / (2 * abs(r * math.cos(phi))))
            max_q = max(q)
            if (q and max_q >= self.max_q) or order > 15:
                found = True
        return order

    def __order_max_q_cheby_2(self):
        found = False
        order = 0
        wn = self.wp
        q = []
        while not found:
            order = order + 1
            zeros, poles, gain = signal.cheby2(order, self.template.att_s, wn, self.filter_t, analog=True, output='zpk')
            for p in poles:
                r, phi = cmath.polar(p)
                q.append(r / (2 * abs(r * math.cos(phi))))
            max_q = max(q)
            if (q and max_q >= self.max_q) or order > 15:
                found = True
        return order

    def __order_max_q_cauer(self):
        found = False
        order = 0
        wn = self.wp
        q = []
        while not found:
            order = order + 1
            zeros, poles, gain = signal.ellip(order, self.template.att_p, self.template.att_s, wn, self.filter_t,
                                              analog=True, output='zpk')
            for p in poles:
                r, phi = cmath.polar(p)
                q.append(r / (2 * abs(r * math.cos(phi))))
            max_q = max(q)
            if (q and max_q >= self.max_q) or order > 15:
                found = True
        return order

    def compute_approximation(self):
        if self.approx_type == 'butterworth':
            self.__compute_approximation_denorm_butter()
            self.__compute_approximation_norm_butter()
        elif self.approx_type == 'bessel':
            self.__compute_approximation_denorm_bessel()
            self.__compute_approximation_norm_bessel()
        elif self.approx_type == 'cheby_1':
            self.__compute_approximation_denorm_cheby_1()
            self.__compute_approximation_norm_cheby_1()
        elif self.approx_type == 'cheby_2':
            self.__compute_approximation_denorm_cheby_2()
            self.__compute_approximation_norm_cheby_2()
        elif self.approx_type == 'cauer':
            self.__compute_approximation_denorm_cauer()
            self.__compute_approximation_norm_cauer()
            pass
        return

    def __compute_approximation_denorm_butter(self):
        self.num, self.den = signal.butter(self.order, self.wn, self.filter_t, analog=True, output='ba')
        self.zeros, self.poles, self.gain = signal.butter(self.order, self.wn, self.filter_t, analog=True, output='zpk')
        self.sos = signal.butter(self.order, self.wn, self.filter_t, analog=True, output='sos')

    def __compute_approximation_norm_butter(self):
        self.num_norm, self.den_norm = signal.butter(self.order, 1, 'lowpass', analog=True, output='ba')
        self.zeros_norm, self.poles_norm, self.gain_norm = signal.butter(self.order, self.wn_N, 'lowpass', analog=True,
                                                                         output='zpk')
        self.sos_norm = signal.butter(self.order, self.wn_N, 'lowpass', analog=True, output='sos')

    def __compute_approximation_denorm_bessel(self):
        self.num, self.den = signal.bessel(self.order, self.wp, self.filter_t, analog=True, output='ba')
        self.zeros, self.poles, self.gain = signal.bessel(self.order, self.wp, self.filter_t, analog=True, output='zpk')
        self.sos = signal.bessel(self.order, self.wp, self.filter_t, analog=True, output='sos')

    def __compute_approximation_norm_bessel(self):
        self.num_norm, self.den_norm = signal.bessel(self.order, 1, 'lowpass', analog=True, output='ba')
        self.zeros_norm, self.poles_norm, self.gain_norm = signal.bessel(self.order, self.template.omega_pN, 'lowpass',
                                                                         analog=True, output='zpk')
        self.sos_norm = signal.bessel(self.order, self.template.omega_pN, 'lowpass', analog=True, output='sos')

    def __compute_approximation_denorm_cheby_1(self):
        self.num, self.den = signal.cheby1(self.order, self.template.att_p, self.wn, self.filter_t, analog=True,
                                           output='ba')
        self.zeros, self.poles, self.gain = signal.cheby1(self.order, self.template.att_p, self.wn, self.filter_t,
                                                          analog=True, output='zpk')
        self.sos = signal.cheby1(self.order, self.template.att_p, self.wn, self.filter_t, analog=True, output='sos')

    def __compute_approximation_norm_cheby_1(self):
        self.num_norm, self.den_norm = signal.cheby1(self.order, self.template.att_p, 1, 'lowpass', analog=True,
                                                     output='ba')
        self.zeros_norm, self.poles_norm, self.gain_norm = signal.cheby1(self.order, self.template.att_p, self.wn_N,
                                                                         'lowpass', analog=True, output='zpk')
        self.sos_norm = signal.cheby1(self.order, self.template.att_p, self.wn_N, 'lowpass', analog=True, output='sos')

    def __compute_approximation_denorm_cheby_2(self):
        self.num, self.den = signal.cheby2(self.order, self.template.att_s, self.wn, self.filter_t, analog=True,
                                           output='ba')
        self.zeros, self.poles, self.gain = signal.cheby2(self.order, self.template.att_s, self.wn, self.filter_t,
                                                          analog=True, output='zpk')
        self.sos = signal.cheby2(self.order, self.template.att_s, self.wn, self.filter_t, analog=True, output='sos')

    def __compute_approximation_norm_cheby_2(self):
        self.num_norm, self.den_norm = signal.cheby2(self.order, self.template.att_s, self.wn_N, 'lowpass', analog=True,
                                                     output='ba')
        self.zeros_norm, self.poles_norm, self.gain_norm = signal.cheby2(self.order, self.template.att_s, self.wn_N,
                                                                         'lowpass', analog=True,
                                                                         output='zpk')  ##Aca cambie un 1 por self.wn_N
        self.sos_norm = signal.cheby2(self.order, self.template.att_s, self.wn_N, 'lowpass', analog=True, output='sos')

    def __compute_approximation_denorm_cauer(self):
        self.num, self.den = signal.ellip(self.order, self.template.att_p, self.template.att_s, self.wn, self.filter_t,
                                          analog=True, output='ba')
        self.zeros, self.poles, self.gain = signal.ellip(self.order, self.template.att_p, self.template.att_s, self.wn,
                                                         self.filter_t, analog=True, output='zpk')
        self.sos = signal.ellip(self.order, self.template.att_p, self.template.att_s, self.wn, self.filter_t,
                                analog=True, output='sos')

    def __compute_approximation_norm_cauer(self):
        self.num_norm, self.den_norm = signal.ellip(self.order, self.template.att_p, self.template.att_s, self.wn_N,
                                                    'lowpass', analog=True, output='ba')
        self.zeros_norm, self.poles_norm, self.gain_norm = signal.ellip(self.order, self.template.att_p,
                                                                        self.template.att_s, self.wn_N, 'lowpass',
                                                                        analog=True, output='zpk')
        self.sos_norm = signal.ellip(self.order, self.template.att_p, self.template.att_s, self.wn_N, 'lowpass',
                                     analog=True, output='sos')

    def compute_factorization(self):
        self.sos = signal.zpk2sos(selg.zeros, self.poles, self.gain)
        return

    def plot_preview_to_axes(self, axes, axes_N, limits):
        n_str = str(self.order)
        legend = self.approx_type_pretty + ' - Order: ' + n_str
        self.legend = legend
        left = limits['left']
        right = limits['right']
        left_N = limits['left_N']
        right_N = limits['right_N']
        line = None
        line_N = None
        if self.approx_type != 'legendre' and self.approx_type != 'gauss':
            w, h = signal.freqs(self.den, self.num, worN=np.logspace(math.log10(left), math.log10(right), 1000))
            line = axes.semilogx(w, 20 * np.log10(abs(h)), label=legend)
            max_value = np.max(20 * np.log10(abs(h)))
            axes.legend(loc='best')
            bottom, top = axes.get_ylim()
            if max_value > top:
                axes.set_ylim((bottom, max_value))
            w_n, h_n = signal.freqs(self.den_norm, self.num_norm,
                                    np.logspace(math.log10(left_N), math.log10(right_N), 1000))
            line_N = axes_N.semilogx(w_n, 20 * np.log10(abs(h_n)), label=legend)
            axes_N.legend(loc='best')
            max_value = np.max(20 * np.log10(abs(h_n)))
            axes.legend(loc='best')
            bottom, top = axes_N.get_ylim()
            if max_value > top:
                axes_N.set_ylim((bottom, max_value))
        return [line, line_N]

    def plot_attenuation_to_axes(self, axes):
        legend = 'Attenuation'
        n_points = 10000
        w, h = signal.freqs(self.den, self.num, n_points)
        line = axes.semilogx(w, 20 * np.log10(abs(h)), label=legend)
        axes.set_xlabel(r'$\omega [rad/seg]$')
        axes.set_ylabel(r'$|A(\omega)| [dB]$')
        axes.legend(loc='best')
        first_w = w.item(0)
        last_w = w.item(n_points - 1)

        return first_w, last_w

    def plot_phase_to_axes(self, axes, limits):
        legend = 'Phase'
        n_points = 1000
        sys = signal.TransferFunction(self.num, self.den)
        left = limits['left']
        right = limits['right']
        w_in = np.logspace(np.log10(left), np.log10(right), num=n_points)
        w, mag, phase = signal.bode(sys, w_in)
        line = axes.semilogx(w, phase, label=legend)
        axes.set_xlabel(r'$\omega [rad/seg]$')
        axes.set_ylabel(r'Phase')
        first_w = w.item(0)
        last_w = w[-1]

        return first_w, last_w

    def plot_group_delay_to_axes(self, axes, limits):
        legend = 'Group Delay'
        n_points = 1000
        w, gd = signal.group_delay((self.num, self.den))
        axes.plot(w, gd, label=legend)
        axes.set_xlabel(r'Frequency [rad/sample]')
        axes.set_ylabel(r'Group Delay [samples]')
        first_w = w.item(0)
        last_w = w[-1]

        return first_w, last_w

    def plot_step_response_to_axes(self, axes):
        legend = 'Step Response'
        n_points = 1000
        t, y = signal.step((self.num, self.den), N=n_points)
        axes.plot(t, y, label=legend)
        axes.set_xlabel(r'Time(seg)')
        axes.set_ylabel(r'V[Volts]')

    def plot_s_plane_to_axes(self, axes):
        legend_poles = 'Poles'
        marker_poles = 'x'
        legend_zeros = 'Zeros'
        marker_zeros = 'o'

        r_poles = []
        theta_poles = []
        r_zeros = []
        theta_zeros = []

        for complex in self.poles_norm:
            r = np.real(complex)
            i = np.imag(complex)
            rho = np.sqrt(r ** 2 + i ** 2)
            theta = np.arctan2(i, r)
            r_poles.append(rho)
            theta_poles.append(theta)

        for complex in self.zeros_norm:
            r = np.real(complex)
            i = np.imag(complex)
            rho = np.sqrt(r ** 2 + i ** 2)
            theta = np.arctan2(i, r)
            r_zeros.append(rho)
            theta_zeros.append(theta)

        if theta_poles:
            axes.scatter(theta_poles, r_poles, label=legend_poles, marker=marker_poles)
        if theta_zeros:
            axes.scatter(theta_zeros, r_zeros, label=legend_zeros, marker=marker_zeros)
        return

    def plot_freq_response_to_axes(self, axes, limits):
        legend = 'Frequency Response[dB]'
        n_points = 1000
        sys = signal.TransferFunction(self.num, self.den)
        left = limits['left']
        right = limits['right']
        w_in = np.logspace(np.log10(left), np.log10(right), num=n_points)
        w, mag, phase = signal.bode(sys, w_in)
        line = axes.semilogx(w, mag, label=legend)
        axes.set_xlabel(r'$\omega [rad/seg]$')
        axes.set_ylabel(r'Mag')
        first_w = w.item(0)
        last_w = w[-1]

        return first_w, last_w

    def __init_stages(self):
        nums = []
        dens = []
        # get numerator
        for num in self.sos[:, :3]:
            nums.append(num)
        for den in self.sos[:, 3:]:
            dens.append(den)

        n = 2 * len(dens)
        cm = mpl.get_cmap('gist_rainbow', n)
        self.colors = [cm(i) for i in range(n)]

        all_zeros = []
        all_poles = []
        all_gains = []

        self.poles_pairs = []
        self.zeros_pairs = []

        pole_pair_count = 0
        zero_pair_count = 0
        self.stages_2 = []
        for i in range(0, len(nums)):

            zeros, poles, gain = signal.tf2zpk(nums[i], dens[i])
            # Aca podria crear un par de polos y cers y una etapa y asignarselos
            zeros_empty = all(z == 0 for z in zeros)
            poles_empty = all(p == 0 for p in poles)
            temp_stage = Stage2(i)
            temp_stage.set_gain(gain)
            # Aca agrego los pares a su stage correspondiente
            if not poles_empty:
                pole_pair_count += 1
                temp_pair = PolePair2(poles, self.colors[i], pole_pair_count)
                self.poles_pairs.append(temp_pair)
                temp_stage.add_pole_pair(temp_pair)

            if not zeros_empty:
                zero_pair_count += 1
                temp_pair = ZeroPair2(zeros, self.colors[i], zero_pair_count)
                self.zeros_pairs.append(temp_pair)
                temp_stage.add_zero_pair(temp_pair)

            self.stages_2.append(temp_stage)
            self.stages_canvases = []

        # Creo las figuras y los canvas para los plots de cada stage y los linkeo
        for stage in self.stages_2:
            temp_figure = Figure()
            temp_canvas = FigureCanvas(temp_figure)
            stage.link_canvas(temp_canvas)
            stage.link_figure(temp_figure)

    def getZPKList(self):
        ##Recibe un nominador y denominador de una transferencia y devuelve:
        ##  -Zeros: Una lista de listas de zeros agrupados de a par con parte real igual
        ##      Pej: Si los zeros son   Z1 = 1+1j
        ##                              Z2 = 1-1j
        ##                              Z3 = 3
        ##              -> Zeros = [ [1+1j , 1-1j], [3] ]
        ##  -Poles: una lista de listas de polos agrupados de la misma manera
        ##  -K: Ganancia total de la transferenci
        z, p, k = signal.tf2zpk(self.num, self.den)
        zeros = []
        poles = []
        used = []
        colors = self.colors
        self.polar_legends = []
        for i in range(0, len(z)):
            if i not in used:
                grouped = False
                used.append(i)
                for u in range(0, len(z)):
                    if u != i:
                        if np.real(z[i]) == np.real(z[u]) and np.real(z[i]) != 0 and u not in used:
                            temp = [z[i], z[u]]
                            used.append(u)
                            zeros.append(temp)
                            grouped = True
                        elif np.imag(z[i]) ** 2 == np.imag(z[u]) ** 2 and np.imag(z[i]) != 0 and u not in used:
                            temp = [z[i], z[u]]
                            used.append(u)
                            zeros.append(temp)
                            grouped = True
                if grouped == False:
                    temp = [z[i]]
                    zeros.append(temp)
        used = []
        for i in range(0, len(p)):
            if i not in used:
                grouped = False
                used.append(i)
                for u in range(0, len(p)):
                    if u != i:
                        if np.real(p[i]) == np.real(p[u]) and u not in used:
                            used.append(u)
                            temp = [p[i], p[u]]
                            poles.append(temp)
                            grouped = True
                        if np.imag(p[i]) ** 2 == np.imag(p[u]) ** 2 and u not in used:
                            used.append(u)
                            temp = [p[i], p[u]]
                            poles.append(temp)
                            grouped = True
                if grouped == False:
                    temp = [p[i]]
                    poles.append(temp)
        returnPoles = []
        tempPair = []
        returnZeros = []
        for index, pair in enumerate(zeros):
            tempPair.clear()
            for zero in pair:
                zero = Zero(zero, colors[index])
                tempPair.append(zero)
            name = 'Zero Pair ' + str(len(returnZeros) + 1)
            self.polar_legends.append(name)
            returnZeros.append(PointPair(self, tempPair, "zero", name))
        for index, pair in enumerate(poles):
            tempPair.clear()
            for pole in pair:
                pole = Pole(pole, colors[index])
                tempPair.append(pole)
            name = 'Pole Pair ' + str(len(returnPoles) + 1)
            self.polar_legends.append(name)

            returnPoles.append(PointPair(self, tempPair, "pole", name))
        return returnZeros, returnPoles, k

    class stageRadioButton(QtWidgets.QRadioButton):
        stage = None

    def initPoints(self):
        appIndex = self.index
        zeros, poles, k = self.getZPKList()
        _translate = QtCore.QCoreApplication.translate
        polesGroupBoxes = [GroupBoxPar(self.program.tab_2) for i in range(len(poles))]
        zerosGroupBoxes = [GroupBoxPar(self.program.tab_2) for i in range(len(zeros))]
        for index, par in enumerate(zeros):
            groupBox = zerosGroupBoxes[index]

            groupBox.setObjectName("GroupBoxZeros" + str(index) + "approx" + str(appIndex))
            gridLayout = QtWidgets.QGridLayout(groupBox)
            gridLayout.setObjectName("gridLayoutZeros" + str(index) + "approx" + str(appIndex))
            radioButtonLayoutPar = QtWidgets.QVBoxLayout()
            radioButtonLayoutPar.setObjectName("radioButtonLayoutZeros" + str(index) + "approx" + str(appIndex))

            for indexetapa, etapa in enumerate(self.stages):
                radioButtonEtapa = self.stageRadioButton(groupBox)
                radioButtonEtapa.stage = etapa
                radioButtonEtapa.setObjectName(
                    "radioButtonEtapa" + str(indexetapa) + "_Cero" + str(index) + "approx" + str(appIndex))
                radioButtonLayoutPar.addWidget(radioButtonEtapa)
                groupBox.radioButtons.append(radioButtonEtapa)
                radioButtonEtapa.setText(_translate("MainWindow", "Stage " + str(indexetapa + 1)))

            gridLayout.addLayout(radioButtonLayoutPar, 0, 0, 1, 1)
            # self.program.paresGroupBoxLayout.addWidget(groupBox)
            groupBox.setTitle(_translate("MainWindow", "Zeros Pair " + str(index + 1)))
            par.groupBox = groupBox

        for index, par in enumerate(poles):
            groupBox = polesGroupBoxes[index]

            groupBox.setObjectName("GroupBoxPoles" + str(index) + "approx" + str(appIndex))
            gridLayout = QtWidgets.QGridLayout(groupBox)
            gridLayout.setObjectName("gridLayoutPoles" + str(index) + "approx" + str(appIndex))
            radioButtonLayoutPar = QtWidgets.QVBoxLayout()
            radioButtonLayoutPar.setObjectName("radioButtonLayoutPoles" + str(index) + "approx" + str(appIndex))

            for indexetapa, etapa in enumerate(self.stages):
                radioButtonEtapa = self.stageRadioButton(groupBox)
                radioButtonEtapa.stage = etapa
                radioButtonEtapa.setObjectName(
                    "radioButtonEtapa" + str(indexetapa) + "_Polo" + str(index + 1) + "approx" + str(appIndex))
                radioButtonLayoutPar.addWidget(radioButtonEtapa)
                groupBox.radioButtons.append(radioButtonEtapa)
                radioButtonEtapa.setText(_translate("MainWindow", "Stage " + str(indexetapa + 1)))

            gridLayout.addLayout(radioButtonLayoutPar, 0, 0, 1, 1)

            # self.program.paresGroupBoxLayout.addWidget(groupBox)
            groupBox.setTitle(_translate("MainWindow", "Pole Pair " + str(index + 1)))
            par.groupBox = groupBox

        self.zeros, self.poles, self.k = zeros, poles, k

    def initFigures(self):
        self.mainStageFigure = Figure()
        self.mainStageFigureCanvas = FigureCanvas(self.mainStageFigure)

    def getNewStageFigure(self):
        fig = Figure()
        figure_canvas = FigureCanvas(fig)
        return fig, figure_canvas

    def initStages(self):
        for i in range(math.ceil(self.order / 2)):
            name = 'Stage ' + str(i + 1)
            self.stages.append(Etapa(self.getNewStageFigure(), name))
        # self.stages.extend([Etapa(self.getNewStageFigure()) for i in range(math.ceil(self.order/2))])

    def plotToStagesTab(self):
        if self.mainStageFigure is not None:
            if self.ax is not None: self.mainStageFigure.delaxes(self.ax)

            self.ax = self.mainStageFigure.add_subplot(111, projection='polar')
            polarPlot(self.zeros, self.poles, self.ax)
            self.mainStageFigureCanvas.draw()

        else:
            print("Error: Set figure first")

    def updateStages(self):
        for polePair in self.poles:
            for radioButton in polePair.groupBox.radioButtons:
                if radioButton.isChecked() and radioButton.stage.canAddPole:
                    radioButton.stage.addPair(polePair)
                    # hacer not checkeable a los radiobuttons de la stage en los otroas groupboxes
                    selected_stage = radioButton.stage
                    selected_polePair = polePair
                    for polePair_2 in self.poles:
                        for radioButton_2 in polePair_2.groupBox.radioButtons:
                            if radioButton_2.stage == selected_stage and polePair_2 != selected_polePair:
                                radioButton_2.setEnabled(False)
                elif not radioButton.isChecked():
                    try:
                        radioButton.setChecked(False)
                        radioButton.stage.removePair(polePair)
                    except:
                        pass

        for zeroPair in self.zeros:
            for radioButton in zeroPair.groupBox.radioButtons:
                if radioButton.isChecked() and radioButton.stage.canAddZero:
                    radioButton.stage.addPair(zeroPair)
                    # hacer not checkeable a los radiobuttons de la stage en los otroas groupboxes
                    selected_stage = radioButton.stage
                    selected_zeroPair = zeroPair
                    for zeroPair_2 in self.zeros:
                        for radioButton_2 in zeroPair_2.groupBox.radioButtons:
                            if radioButton_2.stage == selected_stage and zeroPair_2 != selected_zeroPair:
                                radioButton_2.setEnabled(False)
                elif radioButton.isChecked() and not radioButton.stage.canAddZero:
                    radioButton.setChecked(False)
                elif not radioButton.isChecked():
                    try:
                        radioButton.setChecked(False)
                        radioButton.stage.removePair(zeroPair)
                    except:
                        pass

                for stage in self.stages:
                    for polePair in self.poles:
                        for radioButton in polePair.groupBox.radioButtons:
                            if radioButton.stage == stage:
                                # If stage can add poles, set stage buttons to enabled in all poles
                                # otherwise set to disable in all poles, except the pole contained in the stage
                                if polePair != stage.assignedPolePair:
                                    radioButton.setEnabled(stage.canAddPole)
                    for zeroPair in self.zeros:
                        for radioButton in zeroPair.groupBox.radioButtons:
                            if radioButton.stage == stage:
                                # If stage can add poles, set stage buttons to enabled in all poles
                                # otherwise set to disable in all poles, except the pole contained in the stage
                                if zeroPair != stage.assignedZeroPair:
                                    radioButton.setEnabled(stage.canAddZero)

        # TODO: si no funciona lo de arriba volver a esto
        # for pairList in [self.zeros, self.poles]:
        #    for pair in pairList:
        #        for radioButton in pair.groupBox.radioButtons:
        #            if radioButton.isChecked():
        #                try: 
        #                    radioButton.stage.addPair(pair)
        #                except: 
        #                    radioButton.setChecked(False)
        #                    print('checked out')
        #            else:
        #                try: 
        #                    radioButton.stage.removePair(pair)
        #                    print('checked out 2')
        #                    radioButton.setChecked(False)
        #                except: 
        #                    continue

        for stage in self.stages:
            stage.plotStage()
