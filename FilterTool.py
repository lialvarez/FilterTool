from PyQt5 import QtGui, QtWidgets, QtCore
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
import FilterToolDesign
import sys
import numpy as np
from Template import Template
from Approximation import Approximation
import matplotlib.pyplot as mpl
from distutils.spawn import find_executable
import math
from Etapa import *


def qt_message_handler(mode, context, message):
    if mode == QtCore.QtInfoMsg:
        mode = 'INFO'
    elif mode == QtCore.QtWarningMsg:
        mode = 'WARNING'
    elif mode == QtCore.QtCriticalMsg:
        mode = 'CRITICAL'
    elif mode == QtCore.QtFatalMsg:
        mode = 'FATAL'
    else:
        mode = 'DEBUG'
    file = open("log.txt", "a+")
    file.write('qt_message_handler: line: %d, func: %s(), file: %s' % (
        context.line, context.function, context.file))
    file.write('  %s: %s\n' % (mode, message))
    file.close()


## Class inherited from Ui_MainWindow to customize GUI design and behaviour
# noinspection PyAttributeOutsideInit
class FilterTool(QtWidgets.QMainWindow, FilterToolDesign.Ui_MainWindow):
    def __init__(self, parent=None):
        super(FilterTool, self).__init__(parent)
        self.setupUi(self)

        # configure to use latex interpreter

        if find_executable('latex'):
            mpl.rc('font', **{'family': 'serif', 'serif': ['Palatino']})
            mpl.rc('text', usetex=True)

        # array containing approximation types
        self.approx_types = ['butterworth', 'bessel', 'cheby_1', 'cheby_2', 'legendre', 'gauss', 'cauer']


        # connect 'preview approximation' button to handler
        self.previewApproxButton.clicked.connect(self.on_preview_clicked)
        # connect clear approximation button to handler
        self.clearApproxButton.clicked.connect(self.on_clear_previews_clicked)
        # connect approximation type combo box to handler
        self.approximationComboBox.currentIndexChanged.connect(self.validate_approximation)
        # connect 'Compute' button to handler.
        self.computePushButton.clicked.connect(self.on_compute_clicked)
        # connect 'Set Template' button to handler
        self.setTemplateButton.clicked.connect(self.on_set_template_clicked)
        # connect filter type combo box selection changed to handler
        self.filterTypeComboBox.currentIndexChanged.connect(self.on_filter_type_changed)
        # connect spinboxes valueChanged to handlers
        self.minAttnDoubleSpinBox.valueChanged.connect(self.validate_template)
        self.maxAttnDoubleSpinBox.valueChanged.connect(self.validate_template)
        self.passBand1DoubleSpinBox.valueChanged.connect(self.validate_template)
        self.passBand2DoubleSpinBox.valueChanged.connect(self.validate_template)
        self.stopBand1DoubleSpinBox.valueChanged.connect(self.validate_template)
        self.stopBand2DoubleSpinBox.valueChanged.connect(self.validate_template)
        # connect approximation input config parameters widgets to handler
        self.approximationComboBox.currentIndexChanged.connect(self.validate_approximation)
        self.minMaxRadioButton.clicked.connect(self.validate_approximation)
        self.minOrderSpinBox.valueChanged.connect(self.validate_approximation)
        self.maxOrderSpinBox.valueChanged.connect(self.validate_approximation)
        self.maxQRadioButton.clicked.connect(self.validate_approximation)
        self.maxQdoubleSpinBox.valueChanged.connect(self.validate_approximation)
        self.customOrderRadioButton.clicked.connect(self.validate_approximation)
        self.customOrderSpinBox.valueChanged.connect(self.validate_approximation)
        # connect plot settings to handlers
        self.attenuationCurvesCheckBox.clicked.connect(self.validate_approximation)
        self.phaseCheckBox.clicked.connect(self.validate_approximation)
        self.groupDelayCheckBox.clicked.connect(self.validate_approximation)
        self.stepResponseCheckBox.clicked.connect(self.validate_approximation)
        self.sPlaneCheckBox.clicked.connect(self.validate_approximation)
        self.freqResponseCheckBox.clicked.connect(self.validate_approximation)

        # connect current approximation in stages page changed
        self.current_approx_comboBox.currentIndexChanged.connect(self.on_current_approx_changed)

        # connect stage plot type to handler
        self.stage_plot_comboBox.currentIndexChanged.connect(self.on_current_stage_plot_type_changed)

        # connect current stage combobox to handler
        self.current_stage_comboBox.currentIndexChanged.connect(self.on_current_stage_changed)

        # connect appoximation selected pole changed
        self.poles_listWidget.currentItemChanged.connect(self.on_selected_approx_pole_changed)

        # connect approximation slected zero changed
        self.zeros_listWidget.currentItemChanged.connect(self.on_selected_approx_zero_changed)

        # connect remove pole from stage button to handler
        self.remove_poleButton.clicked.connect(self.on_remove_stage_pole)

        # connect remove zero from stage button to handler
        self.remove_zeroButton.clicked.connect(self.on_remove_stage_zero)

        # connect add pole to stage to handler
        self.add_poleButton.clicked.connect(self.on_add_pole)

        # connect add zero to stage to handler
        self.add_zeroButton.clicked.connect(self.on_add_zero)

        # connect stage gain spin box to handler
        self.gain_spinBox.valueChanged.connect(self.on_stage_gain_changed)

        # connect remove stage button to handler
        self.delete_stageButton.clicked.connect(self.on_remove_stage)

        # hide stopband2 & passband2 widgets
        self.hide_pass_stop_2()

        # validate template parameters
        self.validate_template()

        # set template_submited to False

        ## @var template_submited
        #  Bool variable indicatign whether the template has been submitted
        self.template_submited = False

        # validate approximation settings
        self.validate_approximation()

        ## @var validation_computed
        #  Bool variable indicating whether the approximation has been computed
        self.validation_computed = False

        # Disable save pushbutton
        self.saveResultsPushButton.setEnabled(False)

        # setup figures and canvases
        self.set_mpl()

        # array containing every preview ploted
        self.plotted_previews = []

        # initialize plot options
        self.any_plots_checked()

        ## @var current_stage
        #  vriable containing the current selected stage in stages tab
        self.current_stage = None

        # indica si alguna vez se ploteo esto
        self.att_plotted = False
        self.phase_plotted = False
        self.step_plotted = False
        self.s_plane_plotted = False
        self.freq_response_plotted = False
        self.group_delay_plotted = False

        self.preview_already_plotted = False
        self.currentApproximation = None
        self.approximations = []
        self.zerosGroupBoxes = []
        self.polesGroupBoxes = []
        self.stages = []
        self.mainStageFigure = None

    ## Helper to hide unused GUI widgets.
    #  Hides unused frequency input fields when 'Low-pass' or 'High-pass'
    #  @details filter type is selected.
    #  @param self The object pointer
    def hide_pass_stop_2(self):
        self.stopBand2_label.hide()
        self.stopBand2DoubleSpinBox.hide()
        self.passBand2_label.hide()
        self.passBand2DoubleSpinBox.hide()

    ## Helper to show hided GUI widgets
    #  @details Shows hided frequency input fields when needed.
    #  @param self The object pointer
    def show_pass_stop_2(self):
        self.stopBand2_label.show()
        self.stopBand2DoubleSpinBox.show()
        self.passBand2_label.show()
        self.passBand2DoubleSpinBox.show()

    ## Helper to instantiate figures and canvas
    #  @details Instantiates every figure and canvas to be used by de program
    #  @param self The object pointer
    def set_mpl(self):
        self.template_fig = Figure()
        self.template_canvas = FigureCanvas(self.template_fig)
        self.att_fig = Figure()
        self.att_canvas = FigureCanvas(self.att_fig)
        self.phase_fig = Figure()
        self.phase_canvas = FigureCanvas(self.phase_fig)
        self.group_fig = Figure()
        self.group_canvas = FigureCanvas(self.group_fig)
        self.step_fig = Figure()
        self.step_canvas = FigureCanvas(self.step_fig)
        self.s_plane_fig = Figure()
        self.s_plane_canvas = FigureCanvas(self.s_plane_fig)
        self.freq_fig = Figure()
        self.freq_canvas = FigureCanvas(self.freq_fig)

        self.composition_figure = Figure()
        self.composition_canvas = FigureCanvas(self.composition_figure)

        self.stage_figure = Figure()
        self.stage_canvas = FigureCanvas(self.stage_figure)
        self.stage_layout.addWidget(self.stage_canvas)

        self.template_toolbar = NavigationToolbar(self.template_canvas, self)
        self.plotsLayout.addWidget(self.template_canvas)
        self.plotsLayout.addWidget(self.template_toolbar)

    ## Handles a filter type combo box index change.
    #  @details Executes whenever the selected index in filter type combo box changes.
    #  @param self The object pointer
    def on_filter_type_changed(self):
        self.validate_template()
        selectedText = self.filterTypeComboBox.currentText()
        if selectedText == 'Low-pass':
            self.hide_pass_stop_2()
        if selectedText == 'High-pass':
            self.hide_pass_stop_2()
        if selectedText == 'Band-pass':
            self.show_pass_stop_2()
        if selectedText == 'Band-reject':
            self.show_pass_stop_2()
        return

    ## Handler that validates a new template parameter setting.
    #  @details Whenever a template parameter is changed, the new configuration is
    #  validated to determine if the user can submit it.
    #  @param self The object pointer
    def validate_template(self):
        selected_filter = self.filterTypeComboBox.currentText()
        template = Template(selected_filter,
                            omega_p1=self.passBand1DoubleSpinBox.value(),
                            omega_s1=self.stopBand1DoubleSpinBox.value(),
                            omega_p2=self.passBand2DoubleSpinBox.value(),
                            omega_s2=self.stopBand2DoubleSpinBox.value(),
                            att_p=self.maxAttnDoubleSpinBox.value(),
                            att_s=self.minAttnDoubleSpinBox.value(),
                            final=False)

        valid = template.is_valid()
        self.setTemplateButton.setEnabled(valid)

    ## Helper to clear all figures
    #  @details Deletes every canvas contained in the plotslayout
    #  @param self The object pointer
    def clear_plot_layout(self):
        for i in reversed(range(self.plotsLayout.count())):
            if isinstance(self.plotsLayout.itemAt(i), FigureCanvas):
                self.plotsLayout.itemAt(i).widget().figure.clf()
                self.plotsLayout.itemAt(i).widget().draw()
                self.plotsLayout.itemAt(i).widget().setParent(None)

    ## Set template handler
    #  @details Executes when a new template is submitted (previously validated)
    #  and plots it to the specified axes.
    #  @param self The object pointer
    def on_set_template_clicked(self):
        # Create the specified template
        ## @var selected_filter
        #  String value. Indicates the filter type selected for the template
        self.selected_filter = self.filterTypeComboBox.currentText()
        ## @var template
        #  Template insatance with the specified parameters.
        self.template = Template(self.selected_filter,
                                 omega_p1=self.passBand1DoubleSpinBox.value(),
                                 omega_s1=self.stopBand1DoubleSpinBox.value(),
                                 omega_p2=self.passBand2DoubleSpinBox.value(),
                                 omega_s2=self.stopBand2DoubleSpinBox.value(),
                                 att_p=self.maxAttnDoubleSpinBox.value(),
                                 att_s=self.minAttnDoubleSpinBox.value(),
                                 final=True)
        # remove all plots widgets. in plotLayout
        self.clear_plot_layout()
        # plot template in specified axes
        self.template_axes = self.template_fig.add_subplot(211)
        self.template_axes_norm = self.template_fig.add_subplot(212)
        self.template.plot_template_in_axes(self.template_axes, self.template_axes_norm)
        # add template_canas to plotLayout
        self.plotsLayout.addWidget(self.template_canvas)
        # show template_canvas
        self.template_canvas.show()
        # redraw template_canvas
        self.template_canvas.draw()
        # set template submitted to True
        self.template_submited = True
        self.validate_approximation()
        pass

    ## Validates the approximation settings
    #  @details Verifies if the approximation settings are valid and whether a template is 
    #  submitted.
    #  @param self The object pointer
    #  @return A boolean value. True if approximation is ready to compute. False otherwise
    def validate_approximation(self):
        valid = False

        approx_index = self.approximationComboBox.currentIndex()
        self.approx_type = self.approx_types[approx_index]

        if approx_index == 4 or approx_index == 5:  # Bessel, Legendre y Gauss
            self.minMaxRadioButton.setEnabled(False)
        else:
            self.minMaxRadioButton.setEnabled(True)

        if self.customOrderRadioButton.isChecked():
            self.restriction = 'custom_order'
            valid = self.customOrderSpinBox.value() > 0
        elif self.maxQRadioButton.isChecked():
            self.restriction = 'max_q'
            valid = self.maxQdoubleSpinBox.value() > 0
        elif self.minMaxRadioButton.isChecked():
            self.restriction = 'min_max_order'
            valid = self.minOrderSpinBox.value() < self.maxOrderSpinBox.value()
            valid = valid and self.minOrderSpinBox.value() > 0
            valid = valid and self.maxOrderSpinBox.value() > 0
        else:
            valid = False
        can_compute = valid and self.any_plots_checked()
        self.computePushButton.setEnabled(can_compute and self.template_submited)
        self.previewApproxButton.setEnabled(valid and self.template_submited)

    ## Checks if any plot option selected
    #  @details Verifies if any of the plots option check-boxes is checked
    #  @param self The object pointer
    #  @return A boolean value. True if any check-box checked. False otherwise
    def any_plots_checked(self):
        ret = False
        ## @var plot_attenuation
        #  Boolean variable that indicates whether the attenuation curve ust be plotted
        self.plot_attenuation = self.attenuationCurvesCheckBox.isChecked()
        ## @var plot_phase
        #  Boolean variable that indicates whether the phase response must be plotted
        self.plot_phase = self.phaseCheckBox.isChecked()
        ## @var plot_group_delay
        #  Boolean variable that indicates whether the group delay must be plotted
        self.plot_group_delay = self.groupDelayCheckBox.isChecked()
        ## @var plot_step_response
        #  Boolean variable that indicates whether the step response must be plotted
        self.plot_step_response = self.stepResponseCheckBox.isChecked()
        ## @var plot_s_plane
        #  Boolean variable that indicates whether the S plane must be plotted
        self.plot_s_plane = self.sPlaneCheckBox.isChecked()
        ## @var plot_freq_resp
        #  Boolean varialbe that indicates whether the frequency response must be plotted
        self.plot_freq_resp = self.freqResponseCheckBox.isChecked()

        ret = self.plot_attenuation or self.plot_phase or self.plot_group_delay or self.plot_step_response or self.plot_s_plane or self.plot_freq_resp
        return ret

    ## Active approx selection changed handler
    #  @details Handles the selection of a different active approximation
    #  in the stages tab.
    def on_current_approx_changed(self):
        selection = str(self.current_approx_comboBox.currentText())
        if selection == 'None':
            self.currentApproximation = None
        else:
            for approximation in self.approximations:
                if selection == approximation.legend:
                    self.currentApproximation = approximation

        # Clear Stage Figure
        self.stage_figure.clf()
        # clear poles list
        self.poles_listWidget.clear()
        # clear zeros list
        self.zeros_listWidget.clear()
        # clear stage combobox
        self.current_stage_comboBox.clear()

        # if current approx is not None:
        if self.currentApproximation != None:
            # clear poles list
            self.poles_listWidget.clear()
            # Add all poles to list
            for pole_pair in self.currentApproximation.poles_pairs:
                self.poles_listWidget.addItem(QtWidgets.QListWidgetItem(pole_pair.label))
            # Add all zeros to list
            for zero_pair in self.currentApproximation.zeros_pairs:
                self.zeros_listWidget.addItem(QtWidgets.QListWidgetItem(zero_pair.label))

            # disable add pole and zeros buttons if list are empty
            self.add_poleButton.setEnabled(self.poles_listWidget.count() > 0)
            self.add_zeroButton.setEnabled(self.zeros_listWidget.count() > 0)

            selected_stage_plot = str(self.stage_plot_comboBox.currentText())

            for stage in self.currentApproximation.stages_2:
                self.current_stage_comboBox.addItem('Stage {}'.format(str(stage.id)))

            # set first stage by defult
            self.current_stage_comboBox.setCurrentIndex(0)

            # Stages poles and zeros
            if selected_stage_plot == 'Poles & Zeros':
                self.currentApproximation.stages_2[0].plot_s_plane_to_figure(self.stage_figure)
                self.stage_polesList.setVisible(True)
                self.stage_zerosList.setVisible(True)
                # list stage poles and zeros
                self.stage_zerosList.clear()
                for zero_pair in self.current_stage.zero_pairs:
                    self.stage_zerosList.addItem(QtWidgets.QListWidgetItem(zero_pair.label))
                self.stage_polesList.clear()
                for pole_pair in self.current_stage.pole_pairs:
                    self.stage_polesList.addItem(QtWidgets.QListWidgetItem(pole_pair.label))
            else:  # else hide
                self.stage_polesList.setVisible(False)
                self.stage_zerosList.setVisible(False)
        self.update_composition()
        self.stage_canvas.draw()

    ## New stage plot type handler
    #  @details Handles the selection of a different plot for the stage
    def on_current_stage_plot_type_changed(self):
        selected_stage_plot = str(self.stage_plot_comboBox.currentText())
        if selected_stage_plot == 'Transfer':
            self.current_stage.plot_transfer_to_figure(self.stage_figure)
        elif selected_stage_plot == 'Poles & Zeros':
            self.current_stage.plot_s_plane_to_figure(self.stage_figure)
            # list stage poles and zeros
            self.stage_zerosList.clear()
            for zero_pair in self.current_stage.zero_pairs:
                self.stage_zerosList.addItem(QtWidgets.QListWidgetItem(zero_pair.label))
            self.stage_polesList.clear()
            for pole_pair in self.current_stage.pole_pairs:
                self.stage_polesList.addItem(QtWidgets.QListWidgetItem(pole_pair.label))
        self.stage_canvas.draw()

    ## Stage selector handler
    #  @details Handles the selection of a different 'active' stage 
    def on_current_stage_changed(self):

        # TODO: en vez de limpiar el canvas, borrar la figure.
        # Clear figure
        self.stage_figure.clf()
        # get the stage index
        index = self.current_stage_comboBox.currentIndex()
        self.current_stage = self.currentApproximation.stages_2[index]

        selected_stage_plot = str(self.stage_plot_comboBox.currentText())
        # Stages poles and zeros
        if selected_stage_plot == 'Poles & Zeros':
            # list stage poles and zeros
            self.stage_zerosList.clear()
            self.current_stage.plot_s_plane_to_figure(self.stage_figure)
            for zero_pair in self.current_stage.zero_pairs:
                self.stage_zerosList.addItem(QtWidgets.QListWidgetItem(zero_pair.label))
            self.stage_polesList.clear()
            for pole_pair in self.current_stage.pole_pairs:
                self.stage_polesList.addItem(QtWidgets.QListWidgetItem(pole_pair.label))
        elif selected_stage_plot == 'Transfer':
            self.current_stage.plot_transfer_to_figure(self.stage_figure)
        self.stage_canvas.draw()
        empty_stages = self.current_stage_comboBox.count() == 0
        if empty_stages:
            self.stage_figure.clf()
            self.composition_figure.clear()

    ## Handles te selection of a new approx pole
    def on_selected_approx_pole_changed(self):
        # if pair is already used disable add pole button.
        index = self.poles_listWidget.currentRow()
        if self.currentApproximation.poles_pairs[index].used or self.current_stage.pole_pairs:
            self.add_poleButton.setEnabled(False)
        else:
            self.add_poleButton.setEnabled(True)

    ## Handles te selection of a new approx zero
    def on_selected_approx_zero_changed(self):
        # if pair is already used disable add pole button.
        index = self.zeros_listWidget.currentRow()
        if self.currentApproximation.zeros_pairs[index].used or self.current_stage.zero_pairs:
            self.add_zeroButton.setEnabled(False)
        else:
            self.add_zeroButton.setEnabled(True)

    ## Handles the click on remove stage pole
    def on_remove_stage_pole(self):
        index = self.stage_polesList.currentRow()
        pole_pair = self.current_stage.pole_pairs[index]
        self.current_stage.remove_pole_pair(pole_pair)
        # update stage poles and zeros list
        self.stage_zerosList.clear()
        for zero_pair in self.current_stage.zero_pairs:
            self.stage_zerosList.addItem(QtWidgets.QListWidgetItem(zero_pair.label))
        self.stage_polesList.clear()
        for pole_pair in self.current_stage.pole_pairs:
            self.stage_polesList.addItem(QtWidgets.QListWidgetItem(pole_pair.label))
        self.current_stage.plot_s_plane_to_figure(self.stage_figure)
        self.stage_canvas.draw()
        self.update_composition()
        self.on_selected_approx_pole_changed()

    ## Handles the click on remove stage zero
    def on_remove_stage_zero(self):
        index = self.stage_zerosList.currentRow()
        zero_pair = self.current_stage.zero_pairs[index]
        self.current_stage.remove_zero_pair(zero_pair)
        # update stage poles and zeros list
        self.stage_zerosList.clear()
        for zero_pair in self.current_stage.zero_pairs:
            self.stage_zerosList.addItem(QtWidgets.QListWidgetItem(zero_pair.label))
        self.stage_polesList.clear()
        for pole_pair in self.current_stage.pole_pairs:
            self.stage_polesList.addItem(QtWidgets.QListWidgetItem(pole_pair.label))
        self.current_stage.plot_s_plane_to_figure(self.stage_figure)
        self.stage_canvas.draw()
        self.update_composition()
        self.on_selected_approx_zero_changed()

    ## Handles the click on add pole to stage
    def on_add_pole(self):
        index = self.poles_listWidget.currentRow()
        pair = self.currentApproximation.poles_pairs[index]
        self.current_stage.add_pole_pair(pair)
        selected_plot = str(self.stage_plot_comboBox.currentText())
        if selected_plot == 'Poles & Zeros':
            self.current_stage.plot_s_plane_to_figure(self.stage_figure)
            self.stage_canvas.draw()
            # list stage poles and zeros
            self.stage_zerosList.clear()
            for zero_pair in self.current_stage.zero_pairs:
                self.stage_zerosList.addItem(QtWidgets.QListWidgetItem(zero_pair.label))
            self.stage_polesList.clear()
            for pole_pair in self.current_stage.pole_pairs:
                self.stage_polesList.addItem(QtWidgets.QListWidgetItem(pole_pair.label))
            self.on_selected_approx_pole_changed()
        self.update_composition()

    ## Handles the click on add zero to stage
    def on_add_zero(self):
        index = self.zeros_listWidget.currentRow()
        pair = self.currentApproximation.zeros_pairs[index]
        self.current_stage.add_zero_pair(pair)
        selected_plot = str(self.stage_plot_comboBox.currentText())
        if selected_plot == 'Poles & Zeros':
            self.current_stage.plot_s_plane_to_figure(self.stage_figure)
            self.stage_canvas.draw()
            # list stage poles and zeros
            self.stage_zerosList.clear()
            for zero_pair in self.current_stage.zero_pairs:
                self.stage_zerosList.addItem(QtWidgets.QListWidgetItem(zero_pair.label))
            self.stage_polesList.clear()
            for pole_pair in self.current_stage.pole_pairs:
                self.stage_polesList.addItem(QtWidgets.QListWidgetItem(pole_pair.label))
        self.on_selected_approx_zero_changed()
        self.update_composition()

    ## Handles the gain change
    def on_stage_gain_changed(self):
        self.current_stage.gain = self.gain_spinBox.value()

        self.update_composition()

    ## Updates the transfer generated by the stages
    def update_composition(self):
        zeros = []
        poles = []
        gain = 1
        for stage in self.currentApproximation.stages_2:
            for zero_pair in stage.zero_pairs:
                for zero in zero_pair.points:
                    zeros.append(zero)
            for pole_pair in stage.pole_pairs:
                for pole in pole_pair.points:
                    poles.append(pole)
            gain *= stage.gain
        title = '{0} Accumulated Stages'.format(self.currentApproximation.legend)
        legend = 'Stage Composition Freq. Response'
        self.composition_figure.clear()
        axes = self.composition_figure.add_subplot(111)
        sys = signal.TransferFunction(signal.ZerosPolesGain(zeros, poles, gain))
        w, mag, phase = signal.bode(sys)
        axes.semilogx(w, mag, label=legend)
        axes.set_xlabel(r'$\omega [rad/seg]$')
        axes.set_ylabel(r'Mag')
        axes.set_title('')
        axes.legend()
        axes.grid(True, which='minor', axis='both')
        axes.yaxis.grid()
        self.clearLayout(self.composition_layout)
        self.composition_layout.addWidget(self.composition_canvas)
        self.composition_canvas.draw()

    ## Handles the deletion of a stage
    def on_remove_stage(self):
        if self.current_stage_comboBox.count() > 1:
            # remove stage from stage combo box
            # should generate currentSelected changed
            index = self.current_stage_comboBox.currentIndex()
            self.current_stage_comboBox.removeItem(index)
            # eliminarle los polos y ceros
            self.current_stage.clear()
            # remove from list
            if self.current_stage is not None:
                self.currentApproximation.stages_2.remove(self.current_stage)
            # rearrange ids
            for i in range(0, len(self.currentApproximation.stages_2)):
                self.currentApproximation.stages_2[i].change_id(i)
            # update stage names in combobox
            for i in range(0, self.current_stage_comboBox.count()):
                self.current_stage_comboBox.setItemText(i, self.currentApproximation.stages_2[i].label)
            self.update_composition()
            self.stage_canvas.draw()
            self.composition_canvas.draw()

    def on_preview_clicked(self):
        index = self.approximationComboBox.currentIndex()
        approximation = Approximation(self, self.template,
                                      restriction=self.restriction,
                                      min_order=self.minOrderSpinBox.value(),
                                      max_order=self.maxOrderSpinBox.value(),
                                      max_q=self.maxQdoubleSpinBox.value(),
                                      custom_order=self.customOrderSpinBox.value(),
                                      approx_type=self.approx_type)

        left, right = self.template_axes.get_xlim()
        left_N, right_N = self.template_axes_norm.get_xlim()
        limits = {'left': left, 'right': right, 'left_N': left_N, 'right_N': right_N}
        preview_lines = approximation.plot_preview_to_axes(self.template_axes, self.template_axes_norm, limits)
        # Add approx to combo box (LICHA)
        self.current_approx_comboBox.addItem(approximation.legend)

        self.approximations.append(approximation)

        for line in preview_lines:
            if line != None:
                self.plotted_previews.append(line)
        # show template_canvas
        self.template_canvas.show()
        # redraw template_canvas
        self.template_canvas.draw()
        self.preview_already_plotted = True
        # self.initStages()

    def on_clear_previews_clicked(self):
        if self.preview_already_plotted:
            for line in self.plotted_previews:
                line.pop(0).remove()
            self.template_canvas.draw()
            self.plotted_previews = []
            self.template_axes.get_legend().remove()
            self.template_axes_norm.get_legend().remove()
            self.template_canvas.draw()
            self.preview_already_plotted = False
        # redraw template to sclae properly
        # remove all plots widgets. in plotLayout
        self.clear_plot_layout()
        # plot template in specified axes
        self.template_axes = self.template_fig.add_subplot(211)
        self.template_axes_norm = self.template_fig.add_subplot(212)
        self.template.plot_template_in_axes(self.template_axes, self.template_axes_norm)
        # add template_canvas to plotLayout
        self.plotsLayout.addWidget(self.template_canvas)
        # show template_canvas
        self.template_canvas.show()
        # redraw template_canvas
        self.template_canvas.draw()
        # set template submitted to True
        self.template_submited = True
        self.approximations.clear()
        # empty currentComboBox list
        self.currentApproxComboBox.setCurrentIndex(0)
        self.currentApproxComboBox.clear()
        self.currentApproxComboBox.addItem("None")
        self.currentApproxComboBox.setCurrentIndex(0)
        # set current approx to None
        self.currentApproximation = None
        # validate approximation for GUI porpuses
        self.validate_approximation()
        self.currentApproxComboBox.update()

        return

    def on_compute_clicked(self):
        # eliminar los plots de plotsLayout
        compute_options = [self.plot_attenuation, self.plot_phase, self.plot_group_delay, self.plot_step_response,
                           self.plot_s_plane, self.plot_freq_resp]
        count = 0

        for opt in compute_options:
            if opt:
                count = count + 1
        self.plots_axes = []
        if self.plot_attenuation:
            # if previously plotted, clear
            if self.att_plotted:
                self.att_axes.clear()
            # create axes in figure
            self.att_axes = self.att_fig.add_subplot(111)
            for approximation in self.approximations:
                left, right = approximation.plot_attenuation_to_axes(self.att_axes)

            t_left, t_right = self.template_axes.get_xlim()
            if t_left > left:
                left = t_left

            if t_right < right:
                right = t_right

            self.att_axes.set_xlim(left, right)
            self.att_axes.set_ylim(auto=True)
            self.att_axes.grid(True, which='minor', axis='both')
            self.att_axes.yaxis.grid()
            self.att_plotted = True
            self.att_canvas.draw()
            self.att_canvas.show()
        if self.plot_phase:
            # if previously plotted, clear
            if self.phase_plotted:
                self.phase_axes.clear()
            # create axes in figure
            self.phase_axes = self.phase_fig.add_subplot(111)
            left, right = self.template_axes.get_xlim()
            limits = {'left': left, 'right': right}
            for approximation in self.approximations:
                p_left, p_right = approximation.plot_phase_to_axes(self.phase_axes, limits)

            if p_left > left:
                left = p_left

            if p_right < right:
                right = p_right
            self.phase_axes.set_xlim(left, right)
            self.phase_axes.set_ylim(auto=True)
            self.phase_axes.legend()
            self.phase_axes.grid(True, which='minor', axis='both')
            self.phase_axes.yaxis.grid()
            self.phase_plotted = True
            self.phase_canvas.draw()
            self.phase_canvas.show()
        if self.plot_group_delay:
            # if previously plotted, clear
            if self.group_delay_plotted:
                self.group_axes.clear()
            # create axes in figure
            self.group_axes = self.group_fig.add_subplot(111)
            left, right = self.template_axes.get_xlim()
            limits = {'left': left, 'right': right}
            for approximation in self.approximations:
                p_left, p_right = approximation.plot_group_delay_to_axes(self.group_axes, limits)

            self.group_axes.set_ylim(auto=True)
            self.group_axes.legend()
            self.group_axes.grid(True, which='minor', axis='both')
            self.group_axes.xaxis.grid()
            self.group_axes.yaxis.grid()
            self.group_delay_plotted = True
            self.group_canvas.draw()
            self.group_canvas.show()
        if self.plot_step_response:
            # if previously plotted, clear
            if self.step_plotted:
                self.step_axes.clear()
            # create axes in figure
            self.step_axes = self.step_fig.add_subplot(111)
            left, right = self.template_axes.get_xlim()
            limits = {'left': left, 'right': right}
            for approximation in self.approximations:
                approximation.plot_step_response_to_axes(self.step_axes)

            self.step_axes.set_ylim(auto=True)
            self.step_axes.legend()
            self.step_axes.grid(True, which='minor', axis='both')
            self.step_axes.xaxis.grid()
            self.step_axes.yaxis.grid()
            self.step_plotted = True
            self.step_canvas.draw()
            self.step_canvas.show()
        if self.plot_s_plane:
            # if previously plotted, clear
            if self.s_plane_plotted:
                self.s_plane_axes.clear()
            # create axes in figure
            self.s_plane_axes = self.s_plane_fig.add_subplot(111, projection='polar')
            left, right = self.template_axes.get_xlim()
            limits = {'left': left, 'right': right}
            for approximation in self.approximations:
                approximation.plot_s_plane_to_axes(self.s_plane_axes)

            self.s_plane_axes.set_ylim(auto=True)
            self.s_plane_axes.legend()
            self.s_plane_plotted = True
            self.s_plane_canvas.draw()
            self.s_plane_canvas.show()
        if self.plot_freq_resp:
            # if previously plotted, clear
            if self.freq_response_plotted:
                self.freq_axes.clear()
            # create axes in figure
            self.freq_axes = self.freq_fig.add_subplot(111)
            left, right = self.template_axes.get_xlim()
            limits = {'left': left, 'right': right}
            for approximation in self.approximations:
                p_left, p_right = approximation.plot_freq_response_to_axes(self.freq_axes, limits)

            if p_left > left:
                left = p_left

            if p_right < right:
                right = p_right
            self.freq_axes.set_xlim(left, right)
            self.freq_axes.set_ylim(auto=True)
            self.freq_axes.legend()
            self.freq_axes.grid(True, which='minor', axis='both')
            self.freq_axes.yaxis.grid()
            self.freq_response_plotted = True
            self.freq_canvas.draw()
            self.freq_canvas.show()
        return

    def on_set_current_approx_clicked(self):
        selectedLegend = str(self.currentApproxComboBox.currentText())
        if selectedLegend == '':
            return
        if selectedLegend != "None":
            for approximation in self.approximations:
                if approximation.legend == selectedLegend: selectedApproximation = approximation
            self.setCurrentApproximation(selectedApproximation)
        else:
            self.setCurrentApproximation(None)

    class groupBoxPar(QtWidgets.QGroupBox):
        radioButtons = []
        par = None

    class stageRadioButton(QtWidgets.QRadioButton):
        stage = None

    def updateStages(self):
        self.currentApproximation.updateStages()

    def clearLayout(self, layout):
        if layout.count() != 0:
            for i in reversed(range(layout.count())):
                tempWidget = layout.itemAt(i).widget()
                if tempWidget is not None and str(tempWidget.accessibleName()) != "scroll": tempWidget.setParent(None)
        return

    def setCurrentApproximation(self, approximation):
        self.clearLayout(self.horizontalLayout_30)
        self.clearLayout(self.horizontalLayout_29)
        self.clearLayout(self.paresGroupBoxLayout)
        self.currentApproximation = approximation
        if approximation is not None:
            self.setStagesPushButton.show()
            self.horizontalLayout_30.addWidget(approximation.mainStageFigureCanvas)
            approximation.plotToStagesTab()
            for stage in approximation.stages:
                self.horizontalLayout_29.addWidget(stage.figurecanvas)
                stage.plotStage()
            for pairList in [approximation.poles, approximation.zeros]:
                for pair in pairList:
                    self.paresGroupBoxLayout.addWidget(pair.groupBox)
            # currentItem=[self.currentApproxComboBox.itemText(count) for count in range(self.currentApproxComboBox.count())]
        else:
            self.setStagesPushButton.hide()


def main():
    QtCore.qInstallMessageHandler(qt_message_handler)
    app = QtWidgets.QApplication(sys.argv)
    app.setStyle('Fusion')
    filter_tool = FilterTool()
    filter_tool.showMaximized()
    app.exec_()


if __name__ == "__main__":
    main()
