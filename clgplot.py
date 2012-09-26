#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, re, pylab
import Tkinter as Tki
from os.path import basename
from numpy import pi, exp, log10, array, arange
from tkFileDialog import askopenfilename
from optparse import OptionParser

def gradient(xs, ys):
    '''Return a new array containing the gradients at all the x points'''
    result = []
    for i in range(0, len(xs)):
        grad1, grad2 = None, None
        if i<len(xs)-1: grad1 = (ys[i+1]-ys[i]) / (xs[i+1]-xs[i])
        if i>0: grad2 = (ys[i]-ys[i-1]) / (xs[i]-xs[i-1])
        if grad1 == None: grad1 = grad2
        if grad2 == None: grad2 = grad1
        result.append((grad1+grad2) / 2)
    return array(result)

class DataSeries:
    '''A lightly wrapped 2-column matrix with a method for reading
    it from a file.'''

    def __init__(self, data, name = None, filename = None):
        self.data = data
        self.filename = filename
        self.name = name
        if name != None:
            self.name = name
        else:
            if filename != None:
                self.name = os.path.basename(filename)
            else:
                self.name = None

    @staticmethod
    def read(filename, col1 = 0, col2 = 1, name = None):
        '''Reads a series from a two-column whitespace-delimited text file. If
        there are more than two columns, the extra ones are ignored. If there
        is a header line (or any other non-numeric line), it is ignored.'''
        
        rows = []
        # 'U' for universal newlines
        with open(filename, 'U') as fh:
            for line in fh.readlines():
                parts = line.split()
                try:
                    position = float(parts[col1])
                    if len(parts) > col2:
                        value = float(parts[col2])
                    else:
                        print 'WARNING: missing data at '+str(position)
                        value = 0
                    rows.append([position, value])
                except ValueError:
                    pass # ignore non-numeric lines
        data = array(rows).transpose()
        return DataSeries(data, name=name, filename=filename)

class Gaussian:

    def __init__(self, m_abs, m, bhalf, dp):
        self.m_abs = m_abs # absolute contribution (not used)
        self.m = m # relative contribution (size of peak)
        self.a = m / (dp * (2*pi)**0.5) # corrected for dispersion
        self.bhalf = bhalf # mean log of field (position of peak)
        self.dp = dp # dispersion parameter (width of peak)

    def evaluate(self, x):
        (a,b,c) = (self.a, self.bhalf, self.dp)
        return a*exp(-(((x-b)**2) / (2*(c**2))))

    def cdf(self, x):
        (a,b,c) = (self.a, self.bhalf, self.dp)
        return 0.5*(1+erf((x-b)/(sqrt(2*c**2))))

    def to_csv_line(self):
        return '%.2f,%.2f,%.2f,%.2f,%.2f' % \
            (self.m_abs, self.m, self.a, self.bhalf, self.dp)

    @staticmethod
    def csv_header():
        return 'M_abs,m_rel,a,Bhalf,DP'
        
class IrmCurves:
    '''A collection of cumulative log-gaussian functions
    which can be evaluated to give a modelled IRM remanence
    for a specified applied field.'''

    def __init__(self, name, sirm, params):
        self.name = name
        self.sirm = sirm
        self.components = [Gaussian(*p) for p in params]

    def evaluate(self, x, normalize = False):
        result = sum([g.evaluate(x) for g in self.components])
        if not normalize: result = result * self.sirm
        return result

    @staticmethod
    def read_file(filename):
        re1 = re.compile(r'^ True SIRM= +([0-9.E-]+)')
        re2 = re.compile(r'^ Abs Cont= +([0-9.E-]+)')
        re3 = re.compile(r'^ Rel Cont= +([.0-9]+) +Mean= +([.0-9]+) +' +
                         r'DP= +([.0-9]+)\s+$')
        infile = open(filename)
        sirm = float(re1.search(infile.readline()).groups()[0])
        params = []
        infile.readline()
        while True:
            comp = infile.readline()
            if not comp.startswith(' Component'): break
            param = [float(re2.search(infile.readline()).groups()[0])]
            param += map(float,re3.search(infile.readline()).groups())
            params.append(param)
            infile.readline() # skip blank line
        return IrmCurves(basename(filename), sirm, params)

    def to_csv_line(self):
        result = self.name + ','
        result += ','.join([c.to_csv_line() for c in self.components])
        return result

    def csv_header(self):
        return 'Sample,' + \
            ','.join([Gaussian.csv_header()] * len(self.components))

def plot_clg_fit(series, curves):

    xs = map(log10, series.data[0][1:])
    ys = series.data[1][1:]
    pylab.plot(xs, gradient(xs, ys) / curves.sirm, marker='o',
             ls='', color='black', markerfacecolor='none', markersize=6)

    xs = arange(0.5, 3, 0.02)
    ys = [curves.evaluate(x, True) for x in xs]
    pylab.plot(xs, ys, linewidth=1.0, color='black')
    j = 0
    for curve in curves.components:
        ys2 = [curve.evaluate(x) for x in xs]
        pylab.plot(xs, ys2, linewidth=0.5, color='black')
        j += 1
    pylab.ylim(ymin=0)
    # fettle_subplot(ax, i,
    #                'log10(Applied field (mT))',
    #                'Gradient of normalized magnetization',
    #                5, 2, 2)
    pylab.ion()
    pylab.show()

class App:

    def __init__(self, master, options):

        self.series = None
        self.curves = None

        frame = Tki.Frame(master)
        frame.pack()

        self.data_button = \
            Tki.Button(frame, text="Choose Data file",
                           command=self.choose_data_file)
        self.data_button.pack(side=Tki.TOP)

        self.irmunmix_button = \
            Tki.Button(frame, text="Choose IrmUnmix file",
                           command=self.choose_curves_file)
        self.irmunmix_button.pack(side=Tki.TOP)

        self.plot_button = \
            Tki.Button(frame, text="Plot data", command=self.plot)
        self.plot_button.pack(side=Tki.TOP)

        self.quit_button = Tki.Button(frame, text="Quit",
                                     command=frame.quit)
        self.quit_button.pack(side=Tki.TOP)

        master.update_idletasks()

        w = master.winfo_screenwidth()
        h = master.winfo_screenheight()
        mastersize = tuple(int(_) for _ in master.geometry().split('+')[0].split('x'))
        x = w/2 - mastersize[0]/2
        y = h/2 - mastersize[1]/2
        # master.geometry("%dx%d+%d+%d" % (mastersize + (x, y)))

        if options.data_file: self.read_data_file(options.data_file)
        if options.curves_file: self.read_curves_file(options.curves_file)
        if options.plot_now: self.plot()
    
    def choose_curves_file(self):
        input_file = \
            askopenfilename(title = 'Select IrmUnmix parameter file')
        read_curves_file(self, input_file)

    def read_curves_file(self, input_file):
        self.curves = IrmCurves.read_file(input_file)

    def choose_data_file(self):
        input_file = \
            askopenfilename(title = 'Select IRM data file')
        read_data_file(self, input_file)

    def read_data_file(self, input_file):
        self.series = DataSeries.read(input_file)

    def plot(self):
        plot_clg_fit(self.series, self.curves)

def main():
    parser = OptionParser(usage = 'usage: clgplot [options]')
    parser.add_option('-d', '--data', dest='data_file',
                      help='Read IRM intensities from FILE.', metavar='FILE')
    parser.add_option('-c', '--curves', dest='curves_file',
                      help='Read curve parameters from FILE.', metavar='FILE')
    parser.add_option('-p', '--plot', action='store_true', dest='plot_now',
                      default=False, help='Plot at once.')
    (options, args) = parser.parse_args()
    
    root = Tki.Tk()
    app = App(root, options)
    root.mainloop()

main()
