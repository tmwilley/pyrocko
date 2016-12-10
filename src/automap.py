#!/usr/bin/env python

import math
import random
import logging

from cStringIO import StringIO

import numpy as num

from pyrocko.guts import Object, Float, Bool, Int, Tuple, String, List
from pyrocko.guts import Unicode, Dict
from pyrocko.guts_array import Array
from pyrocko import orthodrome as od
from pyrocko import gmtpy, topo

logger = logging.getLogger('pyrocko.automap')

earthradius = 6371000.0
r2d = 180./math.pi
d2r = 1./r2d
km = 1000.
d2m = d2r*earthradius
m2d = 1./d2m
cm = gmtpy.cm


def point_in_region(p, r):
    p = [num.mod(x, 360.) for x in p]
    r = [num.mod(x, 360.) for x in r]
    if r[0] <= r[1]:
        blon = r[0] <= p[0] <= r[1]
    else:
        blon = not (r[1] < p[0] < r[0])
    if r[2] <= r[3]:
        blat = r[2] <= p[1] <= r[3]
    else:
        blat = not (r[3] < p[1] < r[2])

    return blon and blat


def darken(c, f=0.7):
    return (c[0]*f, c[1]*f, c[2]*f)


def corners(lon, lat, w, h):
    ll_lat, ll_lon = od.ne_to_latlon(lat, lon, -0.5*h, -0.5*w)
    ur_lat, ur_lon = od.ne_to_latlon(lat, lon, 0.5*h, 0.5*w)
    return ll_lon, ll_lat, ur_lon, ur_lat


def extent(lon, lat, w, h, n):
    x = num.linspace(-0.5*w, 0.5*w, n)
    y = num.linspace(-0.5*h, 0.5*h, n)
    slats, slons = od.ne_to_latlon(lat, lon, y[0], x)
    nlats, nlons = od.ne_to_latlon(lat, lon, y[-1], x)
    south = slats.min()
    north = nlats.max()

    wlats, wlons = od.ne_to_latlon(lat, lon, y, x[0])
    elats, elons = od.ne_to_latlon(lat, lon, y, x[-1])
    elons = num.where(elons < wlons, elons + 360., elons)

    if elons.max() - elons.min() > 180 or wlons.max() - wlons.min() > 180.:
        west = -180.
        east = 180.
    else:
        west = wlons.min()
        east = elons.max()

    return topo.positive_region((west, east, south, north))


class NoTopo(Exception):
    pass


class OutOfBounds(Exception):
    pass


class FloatTile(Object):
    xmin = Float.T()
    ymin = Float.T()
    dx = Float.T()
    dy = Float.T()
    data = Array.T(shape=(None, None), dtype=num.float, serialize_as='table')

    def __init__(self, xmin, ymin, dx, dy, data):
        Object.__init__(self, init_props=False)
        self.xmin = float(xmin)
        self.ymin = float(ymin)
        self.dx = float(dx)
        self.dy = float(dy)
        self.data = data
        self._set_maxes()

    def _set_maxes(self):
        self.ny, self.nx = self.data.shape
        self.xmax = self.xmin + (self.nx-1) * self.dx
        self.ymax = self.ymin + (self.ny-1) * self.dy

    def x(self):
        return self.xmin + num.arange(self.nx) * self.dx

    def y(self):
        return self.ymin + num.arange(self.ny) * self.dy

    def get(self, x, y):
        ix = int(round((x - self.xmin) / self.dx))
        iy = int(round((y - self.ymin) / self.dy))
        if 0 <= ix < self.nx and 0 <= iy < self.ny:
            return self.data[iy, ix]
        else:
            raise OutOfBounds()


class City(Object):
    def __init__(self, name, lat, lon, population=None, asciiname=None):
        name = unicode(name)
        lat = float(lat)
        lon = float(lon)
        if asciiname is None:
            asciiname = name.encode('ascii', errors='replace')

        if population is None:
            population = 0
        else:
            population = int(population)

        Object.__init__(self, name=name, lat=lat, lon=lon,
                        population=population, asciiname=asciiname)

    name = Unicode.T()
    lat = Float.T()
    lon = Float.T()
    population = Int.T()
    asciiname = String.T()


class Map(Object):
    lat = Float.T(optional=True)
    lon = Float.T(optional=True)
    radius = Float.T(optional=True)
    width = Float.T(default=20.)
    height = Float.T(default=14.)
    margins = List.T(Float.T())
    illuminate = Bool.T(default=True)
    skip_feature_factor = Float.T(default=0.02)
    show_grid = Bool.T(default=False)
    show_topo = Bool.T(default=True)
    show_topo_scale = Bool.T(default=False)
    show_center_mark = Bool.T(default=False)
    show_rivers = Bool.T(default=True)
    illuminate_factor_land = Float.T(default=0.5)
    illuminate_factor_ocean = Float.T(default=0.25)
    color_wet = Tuple.T(3, Int.T(), default=(216, 242, 254))
    color_dry = Tuple.T(3, Int.T(), default=(172, 208, 165))
    topo_resolution_min = Float.T(
        default=40.,
        help='minimum resolution of topography [dpi]')
    topo_resolution_max = Float.T(
        default=200.,
        help='maximum resolution of topography [dpi]')
    replace_topo_color_only = FloatTile.T(
        optional=True,
        help='replace topo color while keeping topographic shading')
    topo_cpt_wet = String.T(default='light_sea')
    topo_cpt_dry = String.T(default='light_land')
    axes_layout = String.T(optional=True)
    custom_cities = List.T(City.T())
    gmt_config = Dict.T(String.T(), String.T())
    comment = String.T(optional=True)

    def __init__(self, gmtversion='newest', **kwargs):
        Object.__init__(self, **kwargs)
        self._gmt = None
        self._scaler = None
        self._widget = None
        self._corners = None
        self._wesn = None
        self._minarea = None
        self._coastline_resolution = None
        self._rivers = None
        self._dems = None
        self._have_topo_land = None
        self._have_topo_ocean = None
        self._jxyr = None
        self._prep_topo_have = None
        self._labels = []
        self._gmtversion = gmtversion

    def save(self, outpath, resolution=75., oversample=2., size=None,
             width=None, height=None):

        '''Save the image.

        Save the image to *outpath*. The format is determined by the filename
        extension. Formats are handled as follows: ``'.eps'`` and ``'.ps'``
        produce EPS and PS, respectively, directly with GMT. If the file name
        ends with ``'.pdf'``, GMT output is fed through ``gmtpy-epstopdf`` to
        create a PDF file. For any other filename extension, output is first
        converted to PDF with ``gmtpy-epstopdf``, then with ``pdftocairo`` to
        PNG with a resolution oversampled by the factor *oversample* and
        finally the PNG is downsampled and converted to the target format with
        ``convert``. The resolution of rasterized target image can be
        controlled either by *resolution* in DPI or by specifying *width* or
        *height* or *size*, where the latter fits the image into a square with
        given side length.'''

        gmt = self.gmt
        self.draw_labels()
        self.draw_axes()
        if self.show_topo and self.show_topo_scale:
            self._draw_topo_scale()

        gmt.save(outpath, resolution=resolution, oversample=oversample,
                 size=size, width=width, height=height)

    @property
    def scaler(self):
        if self._scaler is None:
            self._setup_geometry()

        return self._scaler

    @property
    def wesn(self):
        if self._wesn is None:
            self._setup_geometry()

        return self._wesn

    @property
    def widget(self):
        if self._widget is None:
            self._setup()

        return self._widget

    @property
    def layout(self):
        if self._layout is None:
            self._setup()

        return self._layout

    @property
    def jxyr(self):
        if self._jxyr is None:
            self._setup()

        return self._jxyr

    @property
    def pxyr(self):
        if self._pxyr is None:
            self._setup()

        return self._pxyr

    @property
    def gmt(self):
        if self._gmt is None:
            self._setup()

        if self._have_topo_ocean is None:
            self._draw_background()

        return self._gmt

    def _setup(self):
        if not self._widget:
            self._setup_geometry()

        self._setup_lod()
        self._setup_gmt()

    def _setup_geometry(self):
        wpage, hpage = self.width, self.height
        ml, mr, mt, mb = self._expand_margins()
        wpage -= ml + mr
        hpage -= mt + mb

        wreg = self.radius * 2.0
        hreg = self.radius * 2.0
        if wpage >= hpage:
            wreg *= wpage/hpage
        else:
            hreg *= hpage/wpage

        self._corners = corners(self.lon, self.lat, wreg, hreg)
        west, east, south, north = extent(self.lon, self.lat, wreg, hreg, 10)

        x, y, z = ((west, east), (south, north), (-6000., 4500.))

        xax = gmtpy.Ax(mode='min-max', approx_ticks=4.)
        yax = gmtpy.Ax(mode='min-max', approx_ticks=4.)
        zax = gmtpy.Ax(mode='min-max', inc=1000., label='Height',
                       scaled_unit='km', scaled_unit_factor=0.001)

        scaler = gmtpy.ScaleGuru(data_tuples=[(x, y, z)], axes=(xax, yax, zax))

        par = scaler.get_params()

        west = par['xmin']
        east = par['xmax']
        south = par['ymin']
        north = par['ymax']

        self._wesn = west, east, south, north
        self._scaler = scaler

    def _setup_lod(self):
        w, e, s, n = self._wesn
        if self.radius > 1500.*km:
            coastline_resolution = 'i'
            rivers = False
        else:
            coastline_resolution = 'f'
            rivers = True

        self._minarea = (self.skip_feature_factor * self.radius/km)**2

        self._coastline_resolution = coastline_resolution
        self._rivers = rivers

        self._prep_topo_have = {}
        self._dems = {}

        cm2inch = gmtpy.cm/gmtpy.inch

        dmin = 2.0 * self.radius * m2d / (self.topo_resolution_max *
                                          (self.height * cm2inch))
        dmax = 2.0 * self.radius * m2d / (self.topo_resolution_min *
                                          (self.height * cm2inch))

        for k in ['ocean', 'land']:
            self._dems[k] = topo.select_dem_names(k, dmin, dmax, self._wesn)
            if self._dems[k]:
                logger.debug('using topography dataset %s for %s'
                             % (','.join(self._dems[k]), k))

    def _expand_margins(self):
        if len(self.margins) == 0 or len(self.margins) > 4:
            ml = mr = mt = mb = 2.0
        elif len(self.margins) == 1:
            ml = mr = mt = mb = self.margins[0]
        elif len(self.margins) == 2:
            ml = mr = self.margins[0]
            mt = mb = self.margins[1]
        elif len(self.margins) == 4:
            ml, mr, mt, mb = self.margins

        return ml, mr, mt, mb

    def _setup_gmt(self):
        w, h = self.width, self.height
        scaler = self._scaler

        if gmtpy.is_gmt5(self._gmtversion):
            gmtconf = dict(
                MAP_TICK_PEN_PRIMARY='1.25p',
                MAP_TICK_PEN_SECONDARY='1.25p',
                MAP_TICK_LENGTH_PRIMARY='0.2c',
                MAP_TICK_LENGTH_SECONDARY='0.6c',
                FONT_ANNOT_PRIMARY='12p,1,black',
                FONT_LABEL='12p,1,black',
                PS_CHAR_ENCODING='ISOLatin1+',
                MAP_FRAME_TYPE='fancy',
                FORMAT_GEO_MAP='D',
                PS_MEDIA='Custom_%ix%i' % (
                    w*gmtpy.cm,
                    h*gmtpy.cm),
                PS_PAGE_ORIENTATION='portrait',
                MAP_GRID_PEN_PRIMARY='thinnest,0/50/0',
                MAP_ANNOT_OBLIQUE='6')
        else:
            gmtconf = dict(
                TICK_PEN='1.25p',
                TICK_LENGTH='0.2c',
                ANNOT_FONT_PRIMARY='1',
                ANNOT_FONT_SIZE_PRIMARY='12p',
                LABEL_FONT='1',
                LABEL_FONT_SIZE='12p',
                CHAR_ENCODING='ISOLatin1+',
                BASEMAP_TYPE='fancy',
                PLOT_DEGREE_FORMAT='D',
                PAPER_MEDIA='Custom_%ix%i' % (
                    w*gmtpy.cm,
                    h*gmtpy.cm),
                GRID_PEN_PRIMARY='thinnest/0/50/0',
                DOTS_PR_INCH='1200',
                OBLIQUE_ANNOTATION='6')

        gmtconf.update(
            (k.upper(), v) for (k, v) in self.gmt_config.iteritems())

        gmt = gmtpy.GMT(config=gmtconf, version=self._gmtversion)

        layout = gmt.default_layout()

        layout.set_fixed_margins(*[x*cm for x in self._expand_margins()])

        widget = layout.get_widget()
        widget['P'] = widget['J']
        widget['J'] = ('-JA%g/%g' % (self.lon, self.lat)) + '/%(width)gp'
        scaler['R'] = '-R%g/%g/%g/%gr' % self._corners

        # aspect = gmtpy.aspect_for_projection(
        #     gmt.installation['version'], *(widget.J() + scaler.R()))

        aspect = self._map_aspect(jr=widget.J() + scaler.R())
        widget.set_aspect(aspect)

        self._gmt = gmt
        self._layout = layout
        self._widget = widget
        self._jxyr = self._widget.JXY() + self._scaler.R()
        self._pxyr = self._widget.PXY() + [
            '-R%g/%g/%g/%g' % (0, widget.width(), 0, widget.height())]
        self._have_drawn_axes = False
        self._have_drawn_labels = False

    def _draw_background(self):
        self._have_topo_land = False
        self._have_topo_ocean = False
        if self.show_topo:
            self._have_topo = self._draw_topo()

        self._draw_basefeatures()

    def _get_topo_tile(self, k):
        t = None
        demname = None
        for dem in self._dems[k]:
            t = topo.get(dem, self._wesn)
            demname = dem
            if t is not None:
                break

        if not t:
            raise NoTopo()

        return t, demname

    def _prep_topo(self, k):
        gmt = self._gmt
        t, demname = self._get_topo_tile(k)

        if demname not in self._prep_topo_have:

            grdfile = gmt.tempfilename()
            gmtpy.savegrd(
                t.x(), t.y(), t.data, filename=grdfile, naming='lonlat')

            if self.illuminate:
                if k == 'ocean':
                    factor = self.illuminate_factor_ocean
                else:
                    factor = self.illuminate_factor_land

                ilumfn = gmt.tempfilename()
                gmt.grdgradient(
                    grdfile,
                    N='e%g' % factor,
                    A=-45,
                    G=ilumfn,
                    out_discard=True)

                ilumargs = ['-I%s' % ilumfn]
            else:
                ilumargs = []

            if self.replace_topo_color_only:
                t2 = self.replace_topo_color_only
                grdfile2 = gmt.tempfilename()

                gmtpy.savegrd(
                    t2.x(), t2.y(), t2.data, filename=grdfile2,
                    naming='lonlat')

                if gmt.is_gmt5():
                    gmt.grdsample(
                        grdfile2,
                        G=grdfile,
                        n='l',
                        I='%g/%g' % (t.dx, t.dy),
                        R=grdfile,
                        out_discard=True)
                else:
                    gmt.grdsample(
                        grdfile2,
                        G=grdfile,
                        Q='l',
                        I='%g/%g' % (t.dx, t.dy),
                        R=grdfile,
                        out_discard=True)

                gmt.grdmath(
                    grdfile, '0.0', 'AND', '=', grdfile2,
                    out_discard=True)

                grdfile = grdfile2

            self._prep_topo_have[demname] = grdfile, ilumargs

        return self._prep_topo_have[demname]

    def _draw_topo(self):
        widget = self._widget
        scaler = self._scaler
        gmt = self._gmt
        cres = self._coastline_resolution
        minarea = self._minarea

        JXY = widget.JXY()
        R = scaler.R()

        try:
            grdfile, ilumargs = self._prep_topo('ocean')
            gmt.pscoast(D=cres, S='c', A=minarea, *(JXY+R))
            gmt.grdimage(grdfile, C=topo.cpt(self.topo_cpt_wet),
                         *(ilumargs+JXY+R))
            gmt.pscoast(Q=True, *(JXY+R))
            self._have_topo_ocean = True
        except NoTopo:
            self._have_topo_ocean = False

        try:
            grdfile, ilumargs = self._prep_topo('land')
            gmt.pscoast(D=cres, G='c', A=minarea, *(JXY+R))
            gmt.grdimage(grdfile, C=topo.cpt(self.topo_cpt_dry),
                         *(ilumargs+JXY+R))
            gmt.pscoast(Q=True, *(JXY+R))
            self._have_topo_land = True
        except NoTopo:
            self._have_topo_land = False

    def _draw_topo_scale(self, label='Elevation [km]'):
        dry = read_cpt(topo.cpt(self.topo_cpt_dry))
        wet = read_cpt(topo.cpt(self.topo_cpt_wet))
        combi = cpt_merge_wet_dry(wet, dry)
        for level in combi.levels:
            level.vmin /= km
            level.vmax /= km

        topo_cpt = self.gmt.tempfilename() + '.cpt'
        write_cpt(combi, topo_cpt)

        (w, h), (xo, yo) = self.widget.get_size()
        self.gmt.psscale(
            D='%gp/%gp/%gp/%gph' % (xo + 0.5*w, yo - 2.0*gmtpy.cm, w,
                                    0.5*gmtpy.cm),
            C=topo_cpt,
            B='1:%s:' % label)

    def _draw_basefeatures(self):
        gmt = self._gmt
        cres = self._coastline_resolution
        rivers = self._rivers
        minarea = self._minarea

        color_wet = self.color_wet
        color_dry = self.color_dry

        if self.show_rivers and rivers:
            rivers = ['-Ir/0.25p,%s' % gmtpy.color(self.color_wet)]
        else:
            rivers = []

        fill = {}
        if not self._have_topo_land:
            fill['G'] = color_dry

        if not self._have_topo_ocean:
            fill['S'] = color_wet

        gmt.pscoast(
            D=cres,
            W='thinnest,%s' % gmtpy.color(darken(gmtpy.color_tup(color_dry))),
            A=minarea,
            *(rivers+self._jxyr), **fill)

    def _draw_axes(self):
        gmt = self._gmt
        scaler = self._scaler
        widget = self._widget

        if self.axes_layout is None:
            if self.lat > 0.0:
                axes_layout = 'WSen'
            else:
                axes_layout = 'WseN'
        else:
            axes_layout = self.axes_layout

        scale_km = gmtpy.nice_value(self.radius/5.) / 1000.

        if self.show_center_mark:
            gmt.psxy(
                in_rows=[[self.lon, self.lat]],
                S='c20p', W='2p,black',
                *self._jxyr)

        if self.show_grid:
            btmpl = ('%(xinc)gg%(xinc)g:%(xlabel)s:/'
                     '%(yinc)gg%(yinc)g:%(ylabel)s:')
        else:
            btmpl = '%(xinc)g:%(xlabel)s:/%(yinc)g:%(ylabel)s:'

        gmt.psbasemap(
            B=(btmpl % scaler.get_params())+axes_layout,
            L=('x%gp/%gp/%g/%g/%gk' % (
                6./7*widget.width(),
                widget.height()/7.,
                self.lon,
                self.lat,
                scale_km)),
            *self._jxyr)

        if self.comment:
            fontsize = self.gmt.label_font_size()

            _, east, south, _ = self._wesn
            if gmt.is_gmt5():
                row = [
                    1, 0,
                    '%gp,%s,%s' % (fontsize, 0, 'black'), 'BR',
                    self.comment]

                farg = ['-F+f+j']
            else:
                row = [1, 0, fontsize, 0, 0, 'BR', self.comment]
                farg = []

            gmt.pstext(
                in_rows=[row],
                N=True,
                R=(0, 1, 0, 1),
                D='%gp/%gp' % (-fontsize*0.2, fontsize*0.3),
                *(widget.PXY() + farg))

    def draw_axes(self):
        if not self._have_drawn_axes:
            self._draw_axes()
            self._have_drawn_axes = True

    def _have_coastlines(self):
        gmt = self._gmt
        cres = self._coastline_resolution
        minarea = self._minarea

        checkfile = gmt.tempfilename()

        gmt.pscoast(
            M=True,
            D=cres,
            W='thinnest,black',
            A=minarea,
            out_filename=checkfile,
            *self._jxyr)

        with open(checkfile, 'r') as f:
            for line in f:
                ls = line.strip()
                if ls.startswith('#') or ls.startswith('>') or ls == '':
                    continue
                plon, plat = [float(x) for x in ls.split()]
                if point_in_region((plon, plat), self._wesn):
                    return True

        return False

    def have_coastlines(self):
        self.gmt
        return self._have_coastlines()

    def project(self, lats, lons, jr=None):
        onepoint = False
        if isinstance(lats, float) and isinstance(lons, float):
            lats = [lats]
            lons = [lons]
            onepoint = True

        if jr is not None:
            j, r = jr
            gmt = gmtpy.GMT(version=self._gmtversion)
        else:
            j, _, _, r = self.jxyr
            gmt = self.gmt

        f = StringIO()
        gmt.mapproject(j, r, in_columns=(lons, lats), out_stream=f, D='p')
        f.seek(0)
        data = num.loadtxt(f, ndmin=2)
        xs, ys = data.T
        if onepoint:
            xs = xs[0]
            ys = ys[0]
        return xs, ys

    def _map_box(self, jr=None):
        ll_lon, ll_lat, ur_lon, ur_lat = self._corners

        xs_corner, ys_corner = self.project(
            (ll_lat, ur_lat), (ll_lon, ur_lon), jr=jr)

        w = xs_corner[1] - xs_corner[0]
        h = ys_corner[1] - ys_corner[0]

        return w, h

    def _map_aspect(self, jr=None):
        w, h = self._map_box(jr=jr)
        return h/w

    def _draw_labels(self):
        if self._labels:
            fontsize = self.gmt.label_font_size()

            n = len(self._labels)

            lons, lats, texts, sx, sy, styles = zip(*self._labels)

            sx = num.array(sx, dtype=num.float)
            sy = num.array(sy, dtype=num.float)

            xs, ys = self.project(lats, lons)

            w, h = self._map_box()

            dxs = num.zeros(n)
            dys = num.zeros(n)

            for i in xrange(n):
                dx, dy = gmtpy.text_box(
                    texts[i], font=1, fontsize=fontsize, **styles[i])

                dxs[i] = dx
                dys[i] = dy

            la = num.logical_and
            anchors_ok = (
                la(xs + sx + dxs < w, ys + sy + dys < h),
                la(xs - sx - dxs > 0., ys - sy - dys > 0.),
                la(xs + sx + dxs < w, ys - sy - dys > 0.),
                la(xs - sx - dxs > 0., ys + sy + dys < h),
            )

            arects = [
                (xs, ys, xs + sx + dxs, ys + sy + dys),
                (xs - sx - dxs, ys - sy - dys, xs, ys),
                (xs, ys - sy - dys, xs + sx + dxs, ys),
                (xs - sx - dxs, ys, xs, ys + sy + dys)]

            def no_points_in_rect(xs, ys, xmin, ymin, xmax, ymax):
                xx = not num.any(la(la(xmin < xs, xs < xmax),
                                    la(ymin < ys, ys < ymax)))
                return xx

            for i in xrange(n):
                for ianch in xrange(4):
                    anchors_ok[ianch][i] &= no_points_in_rect(
                        xs, ys, *[xxx[i] for xxx in arects[ianch]])

            anchor_choices = []
            anchor_take = []
            for i in xrange(n):
                choices = [ianch for ianch in xrange(4)
                           if anchors_ok[ianch][i]]
                anchor_choices.append(choices)
                if choices:
                    anchor_take.append(choices[0])
                else:
                    anchor_take.append(None)

            def roverlaps(a, b):
                return (a[0] < b[2] and b[0] < a[2] and
                        a[1] < b[3] and b[1] < a[3])

            def cost(anchor_take):
                noverlaps = 0
                for i in xrange(n):
                    for j in xrange(n):
                        if i != j:
                            i_take = anchor_take[i]
                            j_take = anchor_take[j]
                            if i_take is None or j_take is None:
                                continue
                            r_i = [xxx[i] for xxx in arects[i_take]]
                            r_j = [xxx[j] for xxx in arects[j_take]]
                            if roverlaps(r_i, r_j):
                                noverlaps += 1

                return noverlaps

            cur_cost = cost(anchor_take)
            imax = 30
            while cur_cost != 0 and imax > 0:
                for i in xrange(n):
                    for t in anchor_choices[i]:
                        anchor_take_new = list(anchor_take)
                        anchor_take_new[i] = t
                        new_cost = cost(anchor_take_new)
                        if new_cost < cur_cost:
                            anchor_take = anchor_take_new
                            cur_cost = new_cost

                imax -= 1

            while cur_cost != 0:
                for i in xrange(n):
                    anchor_take_new = list(anchor_take)
                    anchor_take_new[i] = None
                    new_cost = cost(anchor_take_new)
                    if new_cost < cur_cost:
                        anchor_take = anchor_take_new
                        cur_cost = new_cost
                        break

            anchor_strs = ['BL', 'TR', 'TL', 'BR']

            for i in xrange(n):
                ianchor = anchor_take[i]
                if ianchor is not None:
                    anchor = anchor_strs[ianchor]
                    yoff = [-sy[i], sy[i]][anchor[0] == 'B']
                    xoff = [-sx[i], sx[i]][anchor[1] == 'L']
                    if self.gmt.is_gmt5():
                        row = (
                            lons[i], lats[i],
                            '%i,%s,%s' % (fontsize, 1, 'black'),
                            anchor,
                            texts[i])

                        farg = ['-F+f+j']
                    else:
                        row = (
                            lons[i], lats[i],
                            fontsize, 0, 1, anchor,
                            texts[i])
                        farg = []

                    self.gmt.pstext(
                        in_rows=[row],
                        D='%gp/%gp' % (xoff, yoff),
                        *(self.jxyr + farg),
                        **styles[i])

    def draw_labels(self):
        self.gmt
        if not self._have_drawn_labels:
            self._draw_labels()
            self._have_drawn_labels = True

    def add_label(self, lat, lon, text, offset_x=5., offset_y=5., style={}):
        self._labels.append((lon, lat, text, offset_x, offset_y, style))

    def cities_in_region(self):
        from pyrocko import geonames
        cities = geonames.get_cities_region(region=self.wesn, minpop=0)
        cities.extend(self.custom_cities)
        cities.sort(key=lambda x: x.population)
        return cities

    def draw_cities(self,
                    exact=None,
                    include=[],
                    exclude=[],
                    nmax_soft=10,
                    psxy_style=dict(S='s5p', G='black')):

        cities = self.cities_in_region()

        if exact is not None:
            cities = [c for c in cities if c.name in exact]
            minpop = None
        else:
            cities = [c for c in cities if c.name not in exclude]
            minpop = 10**3
            for minpop_new in [1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6, 3e6, 1e7]:
                cities_new = [
                    c for c in cities
                    if c.population > minpop_new or c.name in include]

                if len(cities_new) == 0 or (
                        len(cities_new) < 3 and len(cities) < nmax_soft*2):
                    break

                cities = cities_new
                minpop = minpop_new
                if len(cities) <= nmax_soft:
                    break

        if cities:
            lats = [c.lat for c in cities]
            lons = [c.lon for c in cities]

            self.gmt.psxy(
                in_columns=(lons, lats),
                *self.jxyr, **psxy_style)

            for c in cities:
                try:
                    text = c.name.encode('iso-8859-1')
                except UnicodeEncodeError:
                    text = c.asciiname

                self.add_label(c.lat, c.lon, text)

        self._cities_minpop = minpop


def rand(mi, ma):
    mi = float(mi)
    ma = float(ma)
    return random.random() * (ma-mi) + mi


def split_region(region):
    west, east, south, north = topo.positive_region(region)
    if east > 180:
        return [(west, 180., south, north),
                (-180., east-360., south, north)]
    else:
        return [region]


class CPTLevel(Object):
    vmin = Float.T()
    vmax = Float.T()
    color_min = Tuple.T(3, Float.T())
    color_max = Tuple.T(3, Float.T())


class CPT(Object):
    color_below = Tuple.T(3, Float.T(), optional=True)
    color_above = Tuple.T(3, Float.T(), optional=True)
    color_nan = Tuple.T(3, Float.T(), optional=True)
    levels = List.T(CPTLevel.T())

    def scale(self, vmin, vmax):
        vmin_old, vmax_old = self.levels[0].vmin, self.levels[-1].vmax
        for level in self.levels:
            level.vmin = (level.vmin - vmin_old) / (vmax_old - vmin_old) * \
                (vmax - vmin) + vmin
            level.vmax = (level.vmax - vmin_old) / (vmax_old - vmin_old) * \
                (vmax - vmin) + vmin


class CPTParseError(Exception):
    pass


def read_cpt(filename):
    with open(filename) as f:
        color_below = None
        color_above = None
        color_nan = None
        levels = []
        try:
            for line in f:
                line = line.strip()
                toks = line.split()

                if line.startswith('#'):
                    continue

                elif line.startswith('B'):
                    color_below = tuple(map(float, toks[1:4]))

                elif line.startswith('F'):
                    color_above = tuple(map(float, toks[1:4]))

                elif line.startswith('N'):
                    color_nan = tuple(map(float, toks[1:4]))

                else:
                    values = map(float, line.split())
                    vmin = values[0]
                    color_min = tuple(values[1:4])
                    vmax = values[4]
                    color_max = tuple(values[5:8])
                    levels.append(CPTLevel(
                        vmin=vmin,
                        vmax=vmax,
                        color_min=color_min,
                        color_max=color_max))

        except:
            raise CPTParseError()

    return CPT(
        color_below=color_below,
        color_above=color_above,
        color_nan=color_nan,
        levels=levels)


def color_to_int(color):
    return tuple(max(0, min(255, int(round(x)))) for x in color)


def write_cpt(cpt, filename):
    with open(filename, 'w') as f:
        for level in cpt.levels:
            f.write(
                '%e %i %i %i %e %i %i %i\n' %
                ((level.vmin, ) + color_to_int(level.color_min) +
                 (level.vmax, ) + color_to_int(level.color_max)))

        if cpt.color_below:
            f.write('B %i %i %i\n' % color_to_int(cpt.color_below))

        if cpt.color_above:
            f.write('F %i %i %i\n' % color_to_int(cpt.color_above))

        if cpt.color_nan:
            f.write('N %i %i %i\n' % color_to_int(cpt.color_nan))


def cpt_merge_wet_dry(wet, dry):
    levels = []
    for level in wet.levels:
        if level.vmin < 0.:
            if level.vmax > 0.:
                level.vmax = 0.

            levels.append(level)

    for level in dry.levels:
        if level.vmax > 0.:
            if level.vmin < 0.:
                level.vmin = 0.

            levels.append(level)

    combi = CPT(
        color_below=wet.color_below,
        color_above=dry.color_above,
        color_nan=dry.color_nan,
        levels=levels)

    return combi


if __name__ == '__main__':
    from pyrocko import util
    util.setup_logging('pyrocko.automap', 'info')

    m = Map(
        lat=rand(40, 70.),
        lon=rand(0., 20.),
        radius=math.exp(rand(math.log(10*km), math.log(2000*km))),
        width=rand(10, 20), height=rand(10, 20),
        show_grid=True,
        show_topo=True,
        illuminate=True)

    m.draw_cities()
    print m
    m.save('map.pdf')
