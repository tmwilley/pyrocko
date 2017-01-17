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
    print n, nout

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
        'displacement': 0,  # treated externally
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
    tpad = None

    if tpad is not None:
        icut1 = max(0, int(num.floor((tp - tpad - out_offset) / out_delta)))
        icut2 = min(n, int(num.ceil((ts + tpad - out_offset) / out_delta)))
    else:
        icut1 = 0
        icut2 = n

    for i, out in enumerate((out_x, out_y, out_z)):
        if out is None:
            continue

        temp = num.fft.irfft(coeffs_stf * specs[i], n)
        temp /= out_delta
        assert temp.size // 2 + 1 == specs[i].size
        assert temp.size == out.size

        temp[:icut1] = 0.0
        temp[icut2:] = 0.0

        if out_quantity == 'displacement_extern':
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


class Step(object):
    def __init__(self):
        pass

    def t_cutoff(self):
        return None

    def __call__(self, f):
        tf = num.ones(f.size, dtype=num.complex)
        tf[0] = 0.
        tf[1:] /= (1.0j * 2. * num.pi*f[1:])
        return tf


class Gauss(object):
    def __init__(self, tau):
        self._tau = tau

    def t_cutoff(self):
        return self._tau * 2.

    def __call__(self, f):
        omega = f * 2. * math.pi

        return num.exp(-(omega**2 * self._tau**2 / 8.))


if __name__ == '__main__':

    x = (1000., 0., 0.)
    f = (0., 0., 0.)
    m6 = (1., 0., 0., 0., 0., 0.)
    #m6 = (1., 1., 1., 0., 0., 0.)

    vp = 3600.
    vs = 2000.
    density = 2800.

    d3d = math.sqrt(x[0]**2 + x[1]**2 + x[2]**2)

    tlen = d3d / vs * 5.

    deltat = 0.001

    n = int(num.round(tlen / deltat))

    out_x = num.zeros(n)
    out_y = num.zeros(n)
    out_z = num.zeros(n)

    import pylab as lab

    tau = 0.01

    nstf = int(round(tau * 5. / deltat))
    t = num.arange(nstf) * deltat
    t0 = nstf * deltat / 2.
    stf = num.exp(-(t-t0)**2/(tau/math.sqrt(2.))**2)
    lab.plot(t, stf)

    nf = True

    add_seismogram(
        vp, vs, density, 10000.0, 10000.0, x, f, m6, 'displacement', deltat, 0.0,
        out_x, out_y, out_z, Gauss(tau), want_far=True, want_near=nf, want_intermediate=nf)

    trs = []
    for out, comp in zip([out_x, out_y, out_z], 'NED'):
        tr = trace.Trace(
            '', '1', 'P', comp, deltat=deltat, tmin=0.0, ydata=out)
        tmin = d3d / vp - t0 * 2.
        tmax = d3d / vp - t0 
        tr2 = tr.chop(tmin, tmax, inplace=False)
        tr.ydata -= num.mean(tr2.ydata)
        trs.append(tr)
    
    def dump(stuff, fn):
        with open(fn, 'w') as f:
            f.write(' '.join('%s' % x for x in stuff))
            f.write('\n')
    
    dump((0., 0., 0., 0.) + m6 + f, 'sources.txt')
    dump(x + (int(nf), 1), 'receivers.txt')
    dump((density, vp, vs), 'material.txt')
    
    stf = num.cumsum(stf)
    stf /= stf[-1]
    stf[0] = 0.0

    data = num.vstack((t, stf)).T
    num.savetxt('stf.txt', data)

    check_call(['ahfull', 'sources.txt', 'receivers.txt', 'material.txt', 'stf.txt', '%g' % deltat, 'ahfull', 'mseed'])


    fns = glob.glob('ahfull-1-*.mseed')

    trs2 = []
    for fn in fns:
        trs2.extend(io.load(fn))

    for tr in trs2:
        tr.set_codes(
            channel={'x': 'N', 'y': 'E', 'z': 'D'}[tr.channel])
        tr.shift(-t0)

    trace.snuffle(trs + trs2)
