
import math
import numpy as num

from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4.QtOpenGL import *
from PyQt4.QtSvg import *

from pyrocko.drum.state import State, Listener
from pyrocko.gui_util import make_QPolygonF, Label
from pyrocko import trace, util


class ProjectXU(object):
    def __init__(self, xr=(0.,1.), ur=(0.,1.)):
        (self.xmin, self.xmax) = xr
        (self.umin, self.umax) = ur
        self.cxu = (self.umax - self.umin) / (self.xmax - self.xmin)

    def u(self, x):
        return self.umin + (x-self.xmin) * self.cxu

    def x(self, u):
        return self.xmin + (u-self.umin) / self.cxu

class ProjectYV(object):
    def __init__(self, yr=(0.,1.), vr=(0.,1.)):
        (self.ymin, self.ymax) = yr
        (self.vmin, self.vmax) = vr
        self.cyv = (self.vmax - self.vmin) / (self.ymax - self.ymin)

    def v(self, y):
        return self.vmin + (y-self.ymin) * self.cyv

    def y(self, v):
        return self.ymin + (v-self.vmin) / self.cyv

class ProjectXYUV(object):
    def __init__(self, xyr=((0.,1.),(0.,1.)), uvr=((0.,1.),(0.,1.))):
        (self.xmin, self.xmax), (self.ymin, self.ymax) = xyr
        (self.umin, self.umax), (self.vmin, self.vmax) = uvr
        self.cxu = (self.umax - self.umin) / (self.xmax - self.xmin)
        self.cyv = (self.vmax - self.vmin) / (self.ymax - self.ymin)

    def u(self, x):
        return self.umin + (x-self.xmin) * self.cxu

    def v(self, y):
        return self.vmin + (y-self.ymin) * self.cyv

    def x(self, u):
        return self.xmin + (u-self.umin) / self.cxu

    def y(self, v):
        return self.ymin + (v-self.vmin) / self.cyv

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


def time_fmt_drumline(t):
    ti = int(t)
    if ti % 60 == 0:
        fmt = '%H:%M'
    else:
        fmt = '%H:%M:%S'

    if ti % (3600*24) == 0:
        fmt = fmt + ' %Y-%m-%d'

    return fmt

class Empty(Exception):
    pass

class DrumLine(QObject, Listener):
    def __init__(self, iline, tmin, tmax, traces, state):
        QObject.__init__(self)
        self.traces = traces
        self.tmin = tmin
        self.tmax = tmax
        self.ymin = 0.
        self.ymax = 1.
        self.ymin_data = 0.
        self.ymax_data = 1.
        self.iline = iline
        self._last_mode = None
        self._ydata_cache = {}
        self._time_per_pixel = None
        self.access_counter = 0

        state.add_listener(self.listener_no_args(self.empty_cache), path='style.trace_resolution')

    def data_range(self, mode='min-max'):
        if not self.traces:
            raise Empty()

        modemap = {
                'min-max': 'minmax',
                'mean-plusminus-1-sigma': 1.,
                'mean-plusminus-2-sigma': 2.,
                'mean-plusminus-4-sigma': 4.,
        }
                
        if self._last_mode != mode:
            mi, ma = trace.minmax(self.traces, key=lambda tr: None, mode=modemap[mode])[None]
            self.ymin_data = mi
            self.ymax_data = ma
            self._last_mode = mode

        return self.ymin_data, self.ymax_data

    def empty_cache(self):
        self._ydata_cache = {}

    def set_yrange(self, ymin, ymax):
        if ymax == ymin:
            ymax = 0.5*(ymin + ymax) + 1.0
            ymin = 0.5*(ymin + ymax) - 1.0

        self.ymin = ymin
        self.ymax = ymax

    @property
    def tyrange(self):
        return (self.tmin, self.tmax), (self.ymin, self.ymax)

    def draw(self, plotenv):
        self.draw_traces(plotenv)
        self.draw_time_label(plotenv)

    def draw_time_label(self, plotenv):
        style = plotenv.style
        font = style.label_textstyle.qt_font
        c = QColor()
        text =  util.time_to_str(self.tmin, format=time_fmt_drumline(self.tmin))
        lab = Label(plotenv, font.pointSize(), plotenv.vmin+font.pointSize(), 
                text,
                label_bg=None, anchor='ML', 
                font = font,
                color=c)
        lab.draw()

    def projected_trace_data(self, tr, project, trace_resolution):
        n = tr.data_len()
        if trace_resolution > 0 and n > 2 and tr.deltat < 0.5/trace_resolution*self._time_per_pixel:
            spp = int(self._time_per_pixel/tr.deltat/trace_resolution)
            if tr not in self._ydata_cache:
                nok = (tr.data_len()/spp) * spp
                ydata_rs = tr.ydata[:nok].reshape((-1, spp))
                ydata = num.empty((nok/spp)*2)
                ydata[::2] = num.min(ydata_rs, axis=1)
                ydata[1::2] = num.max(ydata_rs, axis=1)
                self._ydata_cache[tr] = ydata
            else:
                ydata = self._ydata_cache[tr]

            udata_min = float(project.u(tr.tmin))
            udata_max = float(project.u(tr.tmin+0.5*tr.deltat*spp*(ydata.size-1)))
        else:
            ydata = tr.ydata
            udata_min = float(project.u(tr.tmin))
            udata_max = float(project.u(tr.tmin+tr.deltat*(n-1)))

        vdata = project.v(ydata)
        udata = num.linspace(udata_min, udata_max, vdata.size)
        return udata, vdata

    def draw_traces(self, plotenv):
        project = ProjectXYUV(self.tyrange, plotenv.uvrange)
        tpp = (project.xmax - project.xmin) / (project.umax - project.umin)
        if self._time_per_pixel != tpp:
            self.empty_cache()
            self._time_per_pixel = tpp

        for tr in self.traces:
            udata, vdata = self.projected_trace_data(tr, project, plotenv.style.trace_resolution)
            qpoints = make_QPolygonF( udata, vdata )
            plotenv.setPen(plotenv.style.trace_color.qt_color)
            plotenv.drawPolyline(qpoints)
    
class DrumViewMain(QWidget, Listener):

    def __init__(self, pile, *args):
        QWidget.__init__(self, *args)

        self.state = State()
        self.pile = pile

        st = self.state

        self._drumlines = {}
        self._wheel_pos = 0
        self._iline_float = None
        self._project_iline_to_screen = ProjectYV((st.iline-0.5, st.iline+st.nlines+0.5), (0.,self.height()))
        self._access_counter = 0
        self.state.add_listener(self.listener_no_args(self.state_changed))
        self.state.add_listener(self.listener_no_args(self.drop_cached_drumlines), 'filters')
        self.state.add_listener(self.listener_no_args(self.drop_cached_drumlines), 'nslc')
        self.state.add_listener(self.listener_no_args(self.drop_cached_drumlines), 'tline')

        self.next_nslc()
        if self.pile.tmax:
            self.state.iline = int(math.ceil(self.pile.tmax / self.state.tline))-self.state.nlines

    def drop_cached_drumlines(self):
        self._drumlines = {}

    def next_nslc(self):
        nslc_ids = sorted(self.pile.nslc_ids.keys())
        if nslc_ids:
            try:
                i = nslc_ids.index(self.state.nslc)
            except ValueError:
                i = -1

            self.state.nslc = nslc_ids[(i+1)%len(nslc_ids)]

    def state_changed(self):
        self.update()

    def draw(self, plotenv):
        plotenv.umin = 0.
        plotenv.umax = self.width()
        self.draw_title(plotenv)
        self.draw_time_axis(plotenv)
        self.draw_lines(plotenv)

    def title(self):
        return '.'.join(x for x in self.state.nslc if x)

    def draw_title(self, plotenv):
        style = plotenv.style
        font = style.title_textstyle.qt_font
        c = QColor()
        lab = Label(plotenv, 0.5*(plotenv.umin + plotenv.umax), font.pointSize(),
                self.title(),
                label_bg=None, anchor='TC', 
                font = font,
                color=c)

        lab.draw()

    def draw_time_axis(self, plotenv):
        pass

    def draw_lines(self, plotenv):
        st = self.state
        drumlines_seen = []
        for iline in range(st.iline, st.iline+st.nlines):
            self.update_line(iline)
            drumline = self._drumlines.get(iline, None)
            if drumline:
                drumlines_seen.append(drumline)

        self.autoscale(drumlines_seen)

        top_margin = 50.
        bottom_margin = 50.

        self._project_iline_to_screen = ProjectYV((st.iline-0.5, st.iline+st.nlines-0.5), (top_margin,self.height()-bottom_margin))
        
        for drumline in drumlines_seen:
            plotenv.vmin = self._project_iline_to_screen.v(drumline.iline-0.5)
            plotenv.vmax = self._project_iline_to_screen.v(drumline.iline+0.5)
            drumline.draw(plotenv)
            drumline.access_counter = self._access_counter
            self._access_counter += 1 

        drumlines_by_access = sorted(self._drumlines.values(), key=lambda dl: dl.access_counter)
        for drumline in drumlines_by_access[:-st.npages_cache*st.nlines]:
            del self._drumlines[drumline.iline]

    def autoscale(self, drumlines):
        if not drumlines:
            return

        st = self.state

        data = []
        for drumline in drumlines:
            try:
                data.append(drumline.data_range(st.scaling.base))
            except Empty:
                pass

        if not data:
            data = [ [ 0, 0 ] ]

        mi, ma = num.array(data, dtype=num.float).T
        gain = st.scaling.gain 
        if st.scaling.mode == 'same':
            ymin, ymax = mi.min(), ma.max()
            for drumline in drumlines:
                drumline.set_yrange(ymin/gain, ymax/gain)
        elif st.scaling.mode == 'individual':
            for drumline, ymin, ymax in zip(drumlines, mi, ma):
                drumline.set_yrange(ymin/gain, ymax/gain)
        elif st.scaling.mode == 'fixed':
            for drumline in drumlines:
                drumline.set_yrange(st.scaling.min/gain, st.scaling.max/gain)

    def update_line(self, iline):

        if iline not in self._drumlines:
            st = self.state
            tmin = iline*st.tline
            tmax = (iline+1)*st.tline
            if st.filters:
                tpad = max(x.tpad() for x in st.filters)
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

            dl = self._drumlines[iline] = DrumLine(iline, tmin, tmax, traces, self.state)

    def paintEvent(self, paint_ev ):
        plotenv = PlotEnv(self)
        plotenv.style = self.state.style
    
        if plotenv.style.antialiasing:
            plotenv.setRenderHint( QPainter.Antialiasing )
        
        self.draw( plotenv )
        
    def wheelEvent(self, wheel_event):
        self._wheel_pos += wheel_event.delta()
        n = self._wheel_pos / 120
        self._wheel_pos = self._wheel_pos % 120
        if n == 0:
            return

        amount = max(1.,self.state.nlines/24.)
        wdelta = amount * n
        
        if wheel_event.modifiers() & Qt.ControlModifier:
            proj = self._project_iline_to_screen

            anchor = (proj.y(wheel_event.y()) - proj.ymin) / (proj.ymax - proj.ymin)

            nlines = max(1, self.state.nlines + int(round(wdelta)))

            if self._iline_float is None:
                iline_float = float(self.state.iline)
            else:
                iline_float = self._iline_float

            self._iline_float = iline_float-anchor*wdelta

            self.state.iline = int(round(iline_float))
            self.state.nlines = nlines

        else:
            self.state.iline -= int(wdelta)
            self._iline_float = None

    def keyPressEvent(self, key_event):
       
        keytext = str(key_event.text())

        if key_event.key() == Qt.Key_PageDown:
            self.state.iline += self.state.nlines
            
        elif key_event.key() == Qt.Key_PageUp:
            self.state.iline -= self.state.nlines
            
        elif keytext == '+':
            self.state.scaling.gain *= 1.5
        
        elif keytext == '-':
            self.state.scaling.gain *= 1.0/1.5

        elif keytext == ' ':
            self.next_nslc()

