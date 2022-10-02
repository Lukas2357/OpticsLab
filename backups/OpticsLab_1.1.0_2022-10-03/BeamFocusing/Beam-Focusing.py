#!/usr/bin/env python3.6

import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import Fcts as Fcts
import tkinter as tk
import os
import shutil


class Focus2d:
    def __init__(self, n1, n2, num_ap, wd, f0, wl, power, rho_min, rho_max, rho_steps, phi, z_min, z_max, z_steps, z0):
        self.n1 = n1
        self.n2 = n2
        self.NA = num_ap
        self.WD = wd
        self.f0 = f0
        self.WL = wl
        self.P = power
        self.f = np.sqrt(1 - self.NA ** 2) * self.WD  # Objective focal length
        self.k1 = 2 * np.pi / self.WL * self.n1  # Absolute value of the k-vector in medium 1
        self.k2 = 2 * np.pi / self.WL * self.n2  # Absolute value of the k-vector in medium 2
        self.rho = np.linspace(rho_min, rho_max * self.WL, rho_steps)  # Simulated radii
        self.phi = np.pi * phi  # Simulated azimuthal angles
        self.z = np.linspace(z_min * self.WL, z_max * self.WL, z_steps)  # Simulated axial distances
        self.z0 = z0 * self.WL  # Simulated interface positions
        self.xyAxis = np.append(np.flipud(-self.rho), self.rho) / self.WL  # axis in x-y-plane in units of the WL
        self.zAxis = self.z / self.WL  # z-axis (for visualization)

    def calc_focus(self):
        self.intensities = Fcts.calc_focus(self.k1, self.k2, self.n1, self.n2, self.f0, self.NA, self.f, self.P,
                                           self.z0, self.rho, self.phi, self.z)

    def get_plots_radial(self):
        return Fcts.plot_radial(self.WL, self.z0, self.xyAxis, self.intensities)

    def get_plots_2d(self):
        return Fcts.plot_2d(self.WL, self.z0, self.xyAxis, self.zAxis, self.intensities)


class Focus1d:
    def __init__(self, n1, n2, num_ap, wd, f0, wl, power, z_min, z_max, z_steps, z0_min, z0_max, z0_steps, n_lines):
        self.n1 = n1
        self.n2 = n2
        self.NA = num_ap
        self.WD = wd
        self.f0 = f0
        self.WL = wl
        self.P = power
        self.f = np.sqrt(1 - self.NA ** 2) * self.WD  # Objective focal length
        self.k1 = 2 * np.pi / self.WL * self.n1  # Absolute value of the k-vector in medium 1
        self.k2 = 2 * np.pi / self.WL * self.n2  # Absolute value of the k-vector in medium 2
        self.rho = [0]  # Simulated radii for simulation along z
        self.phi = np.pi / 4  # Simulated azimuthal angles for simulation along z
        self.z = np.linspace(z_min * self.WL, z_max * self.WL, z_steps)  # Sim axial distances for simulation along z
        self.z0 = np.linspace(z0_min * self.WL, z0_max * self.WL, z0_steps)  # Sim interface positions for sim along z
        self.zAxis = self.z / self.WL  # z-axis (for visualization)
        self.z0Axis = self.z0 / self.WL  # z0-axis (for visualization)
        self.n_lines = n_lines

    def calc_focus(self):
        self.params = Fcts.calc_focus_axial(self.WL, self.k1, self.k2, self.n1, self.n2, self.f0, self.NA, self.f,
                                            self.P, self.z0, self.rho, self.phi, self.z, self.zAxis)
        self.intensities = self.params['Intensity']
        self.zPeak = self.params['zPeak']
        self.FWHM = self.params['FWHM']

    def get_plots_axial(self):
        return Fcts.plot_axial(self.z0, self.zAxis, self.intensities, self.n_lines)

    def get_peak_pos_fit(self):
        global fit_peak_params
        fit_peak = Fcts.fit_peak(self.z0Axis, self.zPeak)
        fit_peak_params = fit_peak[0]
        return fit_peak[1]

    def get_peak_fwhm_fit(self):
        global fit_fwhm_params
        fit_fwhm = Fcts.fit_fwhm(self.z0Axis, self.FWHM)
        fit_fwhm_params = fit_fwhm[0]
        return fit_fwhm[1]


def calc_2d(n1, n2, num_ap, wd, f0, wl, power, ready, button2d, button_radial, rho_min, rho_max, rho_steps, phi, z_min,
            z_max, z_steps, z0):

    global focus2d
    global Int2d
    focus2d = Focus2d(n1, n2, num_ap, wd, f0, wl, power, rho_min, rho_max, rho_steps, phi, z_min, z_max, z_steps, z0)
    focus2d.calc_focus()
    Int2d = focus2d.intensities

    button2d.bind("<Button-1>", lambda event: fill_canvas(focus2d.get_plots_2d(), plot_window, 0, 2))
    button_radial.bind("<Button-1>", lambda event: fill_canvas(focus2d.get_plots_radial(), plot_window, 0, 2))

    ready.insert('0', 'I')


def calc_1d(n1, n2, num_ap, wd, f0, wl, power, ready, button_axial, button_fit1, button_fit2, z_min, z_max, z_steps,
            z0_min, z0_max, z0_steps, n_lines):

    global focus1d
    global Int1d
    focus1d = Focus1d(n1, n2, num_ap, wd, f0, wl, power, z_min, z_max, z_steps, z0_min, z0_max, z0_steps, n_lines)
    focus1d.calc_focus()
    Int1d = focus1d.intensities

    button_axial.bind("<Button-1>", lambda event: fill_canvas(focus1d.get_plots_axial(), plot_window, 0, 2))
    button_fit1.bind("<Button-1>", lambda event: fill_canvas(focus1d.get_peak_pos_fit(), plot_window, 0, 2))
    button_fit2.bind("<Button-1>", lambda event: fill_canvas(focus1d.get_peak_fwhm_fit(), plot_window, 0, 2))

    ready.insert('0', 'I')


def fill_canvas(figure, plot_window, column, row):
    canvas = FigureCanvasTkAgg(figure, master=plot_window)
    canvas.draw()
    canvas.get_tk_widget().grid(column=column, row=row)


def create_parameter(column, row, label, unit, default):
    label = tk.Label(text=label, anchor="e", width=16)
    label.grid(column=column, row=row)
    entry = tk.Entry(width=8)
    entry.grid(column=column+1, row=row)
    entry.insert(tk.END, default)
    unit = tk.Label(text=unit, width=4)
    unit.grid(column=column+2, row=row)
    return entry


class Parameter:
    def __init__(self):
        param_window = tk.Tk()
        param_window.wm_title("Parameter Window")
        param_window.geometry('900x200')
        frame = tk.Frame(param_window)
        frame.grid()

        label = tk.Label(text='          Initial Parameters:')
        label.grid(column=0, row=0)
        label = tk.Label(text='          Params for 2D Calc:')
        label.grid(column=3, row=0)
        label = tk.Label(text='          Params for 1D Calc:')
        label.grid(column=6, row=0)

        self.n1 = create_parameter(0, 1, 'IOR 1', '', 1)
        self.n2 = create_parameter(0, 2, 'IOR 2', '', 2.409)
        self.NA = create_parameter(0, 3, 'NA', '', 0.95)
        self.WD = create_parameter(0, 4, 'WD', 'mm', 0.35)
        self.f0 = create_parameter(0, 5, 'f0', '', 1)
        self.WL = create_parameter(0, 6, 'WL', 'nm', 656)
        self.P = create_parameter(0, 7, 'P', 'mW', 1)

        self.max_r = create_parameter(3, 1, 'r max', 'WL', 3)
        self.min_z_2d = create_parameter(3, 2, 'z min', 'WL', -10)
        self.max_z_2d = create_parameter(3, 3, 'z max', 'WL', 10)
        self.steps_r = create_parameter(3, 4, '# steps r', '', 21)
        self.steps_z_2d = create_parameter(3, 5, '# steps z', '', 21)
        self.phi = create_parameter(3, 6, 'phi', 'rad', 0.25)
        self.z0 = create_parameter(3, 7, 'z0', 'WL', 0)

        self.min_z_1d = create_parameter(6, 1, 'z min', 'WL', -5)
        self.max_z_1d = create_parameter(6, 2, 'z max', 'WL', 50)
        self.steps_z_1d = create_parameter(6, 3, '# steps z', '', 101)
        self.min_z0 = create_parameter(6, 4, 'z0 min', 'WL', -10)
        self.max_z0 = create_parameter(6, 5, 'z0 max', 'WL', 0)
        self.steps_z0 = create_parameter(6, 6, '# steps z0', '', 10)


def create_plot_window(Parameter):

    params = Parameter()

    plot_window = tk.Tk()
    plot_window.wm_title("Simulation of focus inside medium")
    plot_window.geometry('2000x600')
    Frame = tk.Frame(plot_window)
    Frame.grid()

    def _save_results():
        if not os.path.exists(os.getcwd()+'/Results'):
            os.makedirs(os.getcwd()+'/Results')
        try:
            for i in range(0, 3):
                label = 'Int_x' if i == 0 else 'Int_y' if i == 1 else 'Int_z'
                np.savetxt(label, Int2d[i])
                shutil.move(os.getcwd()+'/'+label, os.getcwd()+'/Results'+'/'+label)
        except NameError:
            pass
        try:
            np.savetxt('PeakPosFitParams', fit_peak_params)
            shutil.move(os.getcwd()+'/PeakPosFitParams', os.getcwd()+'/Results'+'/PeakPosFitParams')
        except NameError:
            pass
        try:
            np.savetxt('PeakFWHMFitParams', fit_fwhm_params)
            shutil.move(os.getcwd()+'/PeakFWHMFitParams', os.getcwd()+'/Results'+'/PeakFWHMFitParams')
        except NameError:
            pass

    def _help():
        with open(os.getcwd() + "/Help.txt") as file:
            help_file = file.read()
            file.close()
        help_window = tk.Tk()
        help_window.geometry('650x800')
        help_window.grid()
        help_window.wm_title("Help")
        help_text = tk.Text(help_window, height=50, width=80)
        help_text.insert(tk.END, help_file)
        help_text.pack()

    def _quit():
        plot_window.quit()
        plot_window.destroy()

    Frame.grid_rowconfigure(0, minsize=70)
    Frame.grid_columnconfigure(0, minsize=200)
    Frame.grid_columnconfigure(3, minsize=200)
    Frame.grid_columnconfigure(7, minsize=200)
    Frame.grid_columnconfigure(9, minsize=200)

    button2d = tk.Button(Frame, text='Show 2D Plots')
    button2d.grid(column=1, row=0)

    buttonRadial = tk.Button(Frame, text='Show Radial Plot')
    buttonRadial.grid(column=2, row=0)

    entry_ready = tk.Label(Frame, text='Calculations done:')
    entry_ready.grid(column=7, row=0)
    entry_ready = tk.Entry(Frame)
    entry_ready.grid(column=8, row=0)

    button2d_calc = tk.Button(Frame, text='Calculate 2D Focus')
    button2d_calc.bind("<Button-1>", lambda event: calc_2d(float(params.n1.get()), float(params.n2.get()),
                                                           float(params.NA.get()), float(params.WD.get()) * 10 ** -3,
                                                           float(params.f0.get()), float(params.WL.get()) * 10 ** -9,
                                                           float(params.P.get()) * 10 ** -3,
                                                           entry_ready, button2d, buttonRadial, 0,
                                                           float(params.max_r.get()), int(params.steps_r.get()),
                                                           float(params.phi.get()), float(params.min_z_2d.get()),
                                                           float(params.max_z_2d.get()), int(params.steps_z_2d.get()),
                                                           float(params.z0.get())))
    button2d_calc.grid(column=0, row=0)

    buttonAxial = tk.Button(Frame, text='Show Axial Plot')
    buttonAxial.grid(column=4, row=0)

    buttonFit1 = tk.Button(Frame, text='Show Peak Pos Fit')
    buttonFit1.grid(column=5, row=0)

    buttonFit2 = tk.Button(Frame, text='Show Peak FWHM Fit')
    buttonFit2.grid(column=6, row=0)

    button = tk.Button(Frame, text="Save Results", command=_save_results)
    button.grid(column=9, row=0)

    button = tk.Button(Frame, text="Help", command=_help)
    button.grid(column=10, row=0)

    button = tk.Button(Frame, text="Quit", command=_quit)
    button.grid(column=11, row=0)

    button1d_calc = tk.Button(Frame, text='Calculate Axial Focus')
    button1d_calc.bind("<Button-1>", lambda event: calc_1d(float(params.n1.get()), float(params.n2.get()),
                                                           float(params.NA.get()), float(params.WD.get()) * 10 ** -3,
                                                           float(params.f0.get()), float(params.WL.get()) * 10 ** -9,
                                                           float(params.P.get()) * 10 ** -3,
                                                           entry_ready, buttonAxial, buttonFit1, buttonFit2,
                                                           float(params.min_z_1d.get()), float(params.max_z_1d.get()),
                                                           int(params.steps_z_1d.get()), float(params.min_z0.get()),
                                                           float(params.max_z0.get()), int(params.steps_z0.get()),
                                                           int(params.steps_z0.get())))
    button1d_calc.grid(column=3, row=0)

    return plot_window


plot_window = create_plot_window(Parameter)

plot_window.mainloop()
