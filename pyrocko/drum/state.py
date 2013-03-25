
from guts import *
from weakref import ref

class listdict(dict):
    def __missing__(self, k):
        self[k] = []
        return self[k]


class Listener(object):
    
    def listener(self, listener):
        if not hasattr(self, '_strong_refs'):
            self._strong_refs = []

        self._strong_refs.append(listener)
        return listener

    def listener_no_args(self, listener):
        def listener1(k,v):
            listener()

        return self.listener(listener1)


class ScalingMode(StringChoice):
    choices = ['same', 'individual', 'fixed']


class ScalingBase(StringChoice):
    choices = ['min-max', 'mean-plusminus-1-sigma', 
            'mean-plusminus-2-sigma', 'mean-plusminus-4-sigma']


class Filter(Object):

    def apply(self, tr):
        pass

    def tpad(self):
        return 0.0


class ButterLowpass(Filter):
    order = Int.T(default=4)
    corner = Float.T()
    pad_factor = Float.T(default=1.0, optional=True)

    def apply(self, tr):
        tr.lowpass(self.order, self.corner)

    def tpad(self):
        return self.pad_factor/self.corner


class ButterHighpass(Filter):
    order = Int.T(default=4)
    corner = Float.T()
    pad_factor = Float.T(default=1.0, optional=True)

    def apply(self, tr):
        tr.highpass(self.order, self.corner)

    def tpad(self):
        return self.pad_factor/self.corner


class Downsample(Filter):
    deltat = Float.T()

    def apply(self, tr):
        tr.downsample_to(self.deltat)


class Talkie(Object):

    def __setattr__(self, name, value):
        try:
            self.T.get_property(name)
        except ValueError:
            Object.__setattr__(self, name, value)
            return

        root = self
        path = [ name ]
        while True:
            try:
                root, name_at_parent = root.parent
                path.append(name_at_parent)
            except AttributeError:
                break

        oldvalue = getattr(self, name, None)
        if oldvalue:
            if isinstance(oldvalue, Object):
                Object.__delattr__(oldvalue, 'parent')

        if isinstance(value, Object):
            Object.__setattr__(value, 'parent', (self, name))

        Object.__setattr__(self, name, value)
        root.set_event('.'.join(path[::-1]), value)

    def set_event(self, path, value):
        pass


class TextStyle(Talkie):
    family = String.T(default='default', optional=True)
    size = Float.T(default=9.0, optional=True)
    bold = Bool.T(default=False, optional=True)
    italic = Bool.T(default=False, optional=True)

    @property
    def qt_font(self):
        from PyQt4.QtGui import QFont
        font = QFont(self.family)
        font.setPointSizeF(self.size)
        font.setBold(self.bold)
        font.setItalic(self.italic)
        return font


class Color(Object):
    r = Float.T(default=0.0)
    g = Float.T(default=0.0)
    b = Float.T(default=0.0)
    a = Float.T(default=1.0)

    @property
    def qt_color(self):
        from PyQt4.QtGui import QColor
        color = QColor(*(int(round(x*255)) for x in (self.r, self.g, self.b, self.a)))
        return color

class Style(Talkie):
    antialiasing = Bool.T(default=False, optional=True)
    label_textstyle = TextStyle.T(default=TextStyle.D(bold=True))
    title_textstyle = TextStyle.T(default=TextStyle.D(bold=True))
    trace_resolution = Float.T(default=2.0, optional=True)
    trace_color = Color.T(default=Color.D())


class Scaling(Talkie):
    mode = ScalingMode.T(default='same')
    base = ScalingBase.T(default='min-max')
    min = Float.T(default=-1.0, optional=True)
    max = Float.T(default=1.0, optional=True)
    gain = Float.T(default=1.0, optional=True)


class State(Talkie):
    nslc = Tuple.T(4, String.T(default=''))
    tline = Float.T(default=60*60)
    nlines = Int.T(default=24)
    iline = Int.T(default=0)

    style = Style.T(default=Style.D())
    filters = List.T(Filter.T())
    scaling = Scaling.T(default=Scaling.D())

    npages_cache = Int.T(default=10, optional=True)

    def __init__(self, **kwargs):
        self.listeners = listdict()
        Talkie.__init__(self, **kwargs)

    def add_listener(self, listener, path=''):
        self.listeners[path].append(ref(listener))

    def set_event(self, path, value):
        parts = path.split('.')
        for i in xrange(len(parts)+1): 
            path = '.'.join(parts[:i])
            target_refs = self.listeners[path]
            delete = []
            for target_ref in target_refs:
                target = target_ref()
                if target:
                    target(path, value)
                else:
                    delete.append(target_ref)

            for target_ref in delete:
                target_refs.remove(target_ref)


