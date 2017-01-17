import math
import glob
import numpy as num
from pyrocko import trace, io
from pyrocko import ahfullgreen_ext as ext
from subprocess import check_call


class AhfullgreenError(Exception):
    pass


def add_seismogram(
        vp, vs, density, qp, qs, x, f, m6,
        out_quantity, out_delta, out_offset,
        out_x, out_y, out_z, stf,
        want_far=1, want_intermediate=1, want_near=1):

    ns = [out.size for out in (out_x, out_y, out_z) if out is not None]

    if not all(n == ns[0] for n in ns):
        raise AhfullgreenError('length of component arrays must be identical')

    n = ns[0]

    nout = n // 2 + 1

    specs = []
    for out in (out_x, out_y, out_z):
        if out is not None:
            specs.append(num.zeros(nout, dtype=num.complex))
        else:
            specs.append(None)

    x = num.asarray(x, num.float)
    f = num.asarray(f, num.float)
    m6 = num.asarray(m6, num.float)

    oc_c = {
        'displacement': 1,  # treated externally
        'velocity': 1,
        'acceleration': 2}[out_quantity]

    out_spec_delta = float(2.0 * math.pi / (n*out_delta))
    out_spec_offset = 0.0

    omega = out_spec_offset + out_spec_delta * num.arange(nout)

    coeffs_stf = stf(omega/(2.*math.pi)).astype(num.complex)
    if out_offset != 0.0:
        coeffs_stf *= num.exp(1.0j * omega * out_offset)

    r = math.sqrt(num.sum(x**2))

    ext.add_seismogram(
        float(vp), float(vs), float(density), float(qp), float(qs),
        x, f, m6, oc_c, out_spec_delta, out_spec_offset,
        specs[0], specs[1], specs[2], want_far, want_intermediate, want_near)

    tp = r / vp
    ts = r / vs

    tpad = stf.t_cutoff()
    if tpad is None:
        tpad = out_delta*5

    icut1 = max(0, int(num.floor((tp - tpad - out_offset) / out_delta)))
    icut2 = min(n, int(num.ceil((ts + tpad - out_offset) / out_delta)))

    for i, out in enumerate((out_x, out_y, out_z)):
        if out is None:
            continue

        temp = num.fft.irfft(coeffs_stf * specs[i], n)
        temp /= out_delta
        assert temp.size // 2 + 1 == specs[i].size
        assert temp.size == out.size

        m1 = num.mean(temp[:icut1])
        m2 = num.mean(temp[icut2:])

        print m1, m2

        temp -= m1

        if out_quantity == 'displacement':
            out[:] += num.cumsum(temp) * out_delta
        else:
            out[:] += temp


class Impulse(object):
    def __init__(self):
        pass

    def t_cutoff(self):
        return None

    def __call__(self, f):
        return num.ones(f.size, dtype=num.complex)


class Gauss(object):
    def __init__(self, tau):
        self._tau = tau

    def t_cutoff(self):
        return self._tau * 2.

    def __call__(self, f):
        omega = f * 2. * math.pi

        return num.exp(-(omega**2 * self._tau**2 / 8.))
