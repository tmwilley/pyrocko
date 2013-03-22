
import numpy as num

from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4.QtOpenGL import *
from PyQt4.QtSvg import *

from pyrocko.drum.state import State
from pyrocko.gui_util import make_QPolygonF
from pyrocko import trace

class Project(object):
    def __init__(self, xyrange=((0.,1.),(0.,1.)), uvrange=((0.,1.),(0.,1.))):
        (self.xmin, self.xmax), (self.ymin, self.ymax) = xyrange
        (self.umin, self.umax), (self.vmin, self.vmax) = uvrange
        self.cxu = (self.umax - self.umin) / (self.xmax - self.xmin)
        self.cyv = (self.vmax - self.vmin) / (self.ymax - self.ymin)

    def u(self, x):
        return self.umin + (x-self.xmin) * self.cxu

    def v(self, y):
        return self.vmin + (y-self.ymin) * self.cyv


class PlotEnv(QPainter):
    def __init__(self, *args):
        QPainter.__init__(self, *args)
        self.umin = 0.
        self.umax = 1.
        self.vmin = 0.
        self.vmax = 1.

    @property
    def uvrange(self):
        return (self.umin, self.umax), (self.vmin, self.vmax)

def project_trace_data(tr, project):
    vdata = project.v(tr.get_ydata())
    udata_min = float(project.u(tr.tmin))
    udata_max = float(project.u(tr.tmin+tr.deltat*(vdata.size-1)))
    udata = num.linspace(udata_min, udata_max, vdata.size)
    return udata, vdata

class DrumLine(QObject):
    def __init__(self, iline, tmin, tmax, traces):
        self.traces = traces
        self.tmin = tmin
        self.tmax = tmax
        self.ymin = 0.
        self.ymax = 1.
        self.ymin_data = 0.
        self.ymax_data = 1.
        self.iline = iline
        self._last_mode = None

    def data_range(self, mode='min-max'):
        modemap = {
                'min-max': 'minmax',
                'mean-2std': 2.,
                'mean-4std': 4
        }
                
        if self._last_mode != mode:
            mi, ma = trace.minmax(self.traces, key=lambda tr: None, mode=modemap[mode])[None]
            self.ymin_data = mi
            self.ymax_data = ma
            self._last_mode = mode

        return self.ymin_data, self.ymax_data

    def set_yrange(self, ymin, ymax):
        self.ymin = ymin
        self.ymax = ymax

    @property
    def tyrange(self):
        return (self.tmin, self.tmax), (self.ymin, self.ymax)

    def draw(self, plotenv):
        project = Project(self.tyrange, plotenv.uvrange)
        for tr in self.traces:
            udata, vdata = project_trace_data(tr, project)
            qpoints = make_QPolygonF( udata, vdata )
            plotenv.drawPolyline(qpoints)
    
class DrumViewMain(QWidget):

    def __init__(self, pile, *args):
        QWidget.__init__(self, *args)

        self.state = State()
        self.pile = pile
        self.drumlines = {}

    def draw(self, plotenv):
        plotenv.umin = 0.
        plotenv.umax = self.width()
        self.draw_time_axis(plotenv)
        self.draw_lines(plotenv)

    def draw_time_axis(self, plotenv):
        pass

    def draw_lines(self, plotenv):
        st = self.state
        drumlines_seen = []
        for iline in range(st.iline, st.iline+st.nlines):
            self.update_line(iline)
            drumline = self.drumlines.get(iline, None)
            if drumline:
                drumlines_seen.append(drumline)

        self.autoscale(drumlines_seen)
        
        for drumline in drumlines_seen:
            plotenv.vmin = float(drumline.iline-st.iline+1)/st.nlines * self.height()
            plotenv.vmax = float(drumline.iline-st.iline)/st.nlines * self.height()
            drumline.draw(plotenv)

    def autoscale(self, drumlines):
        st = self.state

        data = []
        for drumline in drumlines:
            data.append(drumline.data_range(st.scaling.base))

        mi, ma = num.array(data, dtype=num.float).T
        
        if st.scaling.mode == 'same':
            ymin, ymax = mi.min(), ma.max()
            for drumline in drumlines:
                drumline.set_yrange(ymin, ymax)
        elif st.scaling.mode == 'individual':
            for drumline, ymin, ymax in zip(drumlines, mi, ma):
                drumline.set_yrange(ymin, ymax)
        elif st.scaling.mode == 'fixed':
            for drumline in drumlines:
                drumline.set_yrange(st.scaling.min, st.scaling.max)

    def update_line(self, iline):

        if iline not in self.drumlines:
            st = self.state
            tmin = iline*st.tline
            tmax = (iline+1)*st.tline
            if st.filters:
                tpad = max(x.tpad(factor=1.0) for x in st.filters)
            else:
                tpad = 0.0

            traces = self.pile.all(
                    tmin = iline*st.tline,
                    tmax = (iline+1)*st.tline,
                    tpad = tpad,
                    trace_selector = lambda tr: tr.nslc_id == st.nslc )

            for tr in traces:
                for filter in st.filters:
                    filter.apply(tr)

            if traces:
                dl = self.drumlines[iline] = DrumLine(iline, tmin, tmax, traces)

    def paintEvent(self, paint_ev ):
        '''Called by QT whenever widget needs to be painted'''
        plotenv = PlotEnv(self)
    
        if self.state.drawing.antialiasing:
            plotenv.setRenderHint( QPainter.Antialiasing )
        
        self.draw( plotenv )
        
