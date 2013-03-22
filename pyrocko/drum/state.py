
from guts import *
from weakref import ref

class ScalingMode(StringChoice):
    choices = ['same', 'individual', 'fixed']

class ScalingBase(StringChoice):
    choices = ['min-max', 'mean-2std', 'mean-4std']

class Filter(Object):
    pass

class ButterLowpass(Filter):
    order = Int.T(default=4)
    corner = Float.T()

    def apply(self, tr):
        tr.lowpass(self.order, self.corner)

    def tpad(self, factor=1.0):
        return factor/self.corner

class ButterHighpass(Filter):
    order = Int.T(default=4)
    corner = Float.T()

    def apply(self, tr):
        tr.highpass(self.order, self.corner)

    def tpad(self, factor=1.0):
        return factor/self.corner

class Drawing(Object):
    antialiasing = Bool.T(default=False)


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

class Scaling(Talkie):
    mode = ScalingMode.T(default='same')
    base = ScalingBase.T(default='min-max')
    min = Float.T(default=-1.0, optional=True)
    max = Float.T(default=1.0, optional=True)

class listdict(dict):
    def __missing__(self, k):
        self[k] = []
        return self[k]

class State(Talkie):
    nslc = Tuple.T(4, String.T(default=''))
    tline = Float.T(default=60*60)
    nlines = Int.T(default=24)
    iline = Int.T(default=0)

    drawing = Drawing.T(default=Drawing())
    filters = List.T(Filter.T())
    scaling = Scaling.T(default=Scaling())

    def __init__(self, **kwargs):
        Talkie.__init__(self, **kwargs)
        self.listeners = listdict()

    def add_listener(self, path, listener):
        self.listeners[path].append(ref(listener))

    def set_event(self, path, value):
        if hasattr(self, 'listeners'):
            target_refs = self.listeners[path]
            for target_ref in target_refs:
                target = target_ref()
                if target:
                    target.set_event(path, value)


class Listener:
    def set_event(self, path, value):
        print 'set', path, value

s = State()
s.scaling.min = 10.0

l = Listener()
s.add_listener('scaling.min', l)

s.scaling.min = 12.0




