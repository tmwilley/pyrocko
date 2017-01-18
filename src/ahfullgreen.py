import math
import numpy as num
from pyrocko import trace
from pyrocko import ahfullgreen_ext as ext


class AhfullgreenError(Exception):
    pass


def add_seismogram(
        vp, vs, density, qp, qs, x, f, m6,
        out_quantity, out_delta, out_offset,
        out_x, out_y, out_z, stf,
        want_far=1, want_intermediate=1, want_near=1):

    npad_levelling = 40
    ns = [out.size for out in (out_x, out_y, out_z) if out is not None]

    if not all(n == ns[0] for n in ns):
        raise AhfullgreenError('length of component arrays must be identical')

    n = ns[0]

    x = num.asarray(x, num.float)
    f = num.asarray(f, num.float)
    m6 = num.asarray(m6, num.float)

    r = math.sqrt(num.sum(x**2))

    tp = r / vp
    ts = r / vs

    if ts <= tp:
        raise AhfullgreenError('unsupported material properties')

    tpad = stf.t_cutoff() or out_delta * 10.

    tpad *= 10.

    tstart = tp - tpad - npad_levelling * out_delta
    tstart = out_offset + round((tstart - out_offset) / out_delta) * out_delta

    ntemp = trace.nextpow2(int(math.ceil(
        (ts - tp + 2 * tpad + 2*npad_levelling * out_delta) / out_delta)))

    nspec = ntemp // 2 + 1

    specs = []
    for out in (out_x, out_y, out_z):
        if out is not None:
            specs.append(num.zeros(nspec, dtype=num.complex))
        else:
            specs.append(None)

    oc_c = {
        'displacement': 1,  # treated in post processing
        'velocity': 1,
        'acceleration': 2}[out_quantity]

    out_spec_delta = float(2.0 * math.pi / (ntemp*out_delta))
    out_spec_offset = 0.0

    omega = out_spec_offset + out_spec_delta * num.arange(nspec)

    coeffs_stf = stf(omega/(2.*math.pi)).astype(num.complex)
    coeffs_stf *= num.exp(1.0j * omega * tstart)

    ext.add_seismogram(
        float(vp), float(vs), float(density), float(qp), float(qs),
        x, f, m6, oc_c, out_spec_delta, out_spec_offset,
        specs[0], specs[1], specs[2], want_far, want_intermediate, want_near)

    for i, out in enumerate((out_x, out_y, out_z)):
        if out is None:
            continue

        temp = num.fft.irfft(coeffs_stf * specs[i], ntemp)
        temp /= out_delta
        assert temp.size // 2 + 1 == specs[i].size

        m1 = num.mean(temp[:npad_levelling] * num.linspace(1., 0., npad_levelling))
        m2 = num.mean(temp[-npad_levelling:] * num.linspace(0., 1., npad_levelling))

        print m1, m2

        temp -= m1 * 2.

        if out_quantity == 'displacement':
            temp = num.cumsum(temp) * out_delta

        tmin = max(out_offset, tstart)
        tmax = min(
            out_offset + (n-1) * out_delta,
            tstart + (ntemp-1) * out_delta)

        def ind(t, t0):
            return int(round((t-t0)/out_delta))

        out[ind(tmin, out_offset):ind(tmax, out_offset)+1] \
            = temp[ind(tmin, tstart):ind(tmax, tstart)+1]

        out[:ind(tmin, out_offset)] = 0.
        out[ind(tmax, out_offset)+1:] = temp[ind(tmax, tstart)]


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
