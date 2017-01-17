
import random
import math
import unittest
import logging
import os.path as op
from tempfile import mkdtemp

from subprocess import check_call
import numpy as num

from pyrocko import util, ahfullgreen, trace, io


logger = logging.getLogger('test_gf_ahfull')

km = 1000.


def rand(mi, ma):
    return mi + random.random() * (ma-mi)


class AhfullTestCase(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        unittest.TestCase.__init__(self, *args, **kwargs)
        self.tempdirs = []

    def __del__(self):
        import shutil

        for d in self.tempdirs:
            shutil.rmtree(d)

    def test_ahfull_kiwi(self):
        x = (rand(100., 1000.), rand(100., 1000.), rand(100., 1000.))
        f = (rand(-1., 1.), rand(-1., 1.), rand(-1., 1.))
        m6 = tuple(rand(-1., 1.) for _ in xrange(6))

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

        tau = 0.005

        nstf = int(round(tau * 5. / deltat))
        t = num.arange(nstf) * deltat
        t0 = nstf * deltat / 2.
        stf = num.exp(-(t-t0)**2/(tau/math.sqrt(2.))**2)

        nf = True

        ahfullgreen.add_seismogram(
            vp, vs, density, 1000000.0, 1000000.0, x, f, m6, 'displacement',
            deltat, 0.,
            out_x, out_y, out_z,
            ahfullgreen.Gauss(tau),
            want_far=True, want_near=nf, want_intermediate=nf)

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

        dn = mkdtemp(prefix='test-ahfull-')
        fn_sources = op.join(dn, 'sources.txt')
        fn_receivers = op.join(dn, 'receivers.txt')
        fn_material = op.join(dn, 'material.txt')
        fn_stf = op.join(dn, 'stf.txt')

        dump((0., 0., 0., 0.) + m6 + f, fn_sources)
        dump(x + (int(nf), 1), fn_receivers)
        dump((density, vp, vs), fn_material)

        stf = num.cumsum(stf)
        stf /= stf[-1]
        stf[0] = 0.0

        data = num.vstack((t, stf)).T
        num.savetxt(fn_stf, data)

        check_call(
            ['ahfull', fn_sources, fn_receivers, fn_material, fn_stf,
             '%g' % deltat, op.join(dn, 'ahfull'), 'mseed', '0'],
            stdout=open('/dev/null', 'w'))

        fns = [op.join(dn, 'ahfull-1-%s-1.mseed' % c) for c in 'xyz']

        trs2 = []
        for fn in fns:
            trs2.extend(io.load(fn))

        for tr in trs2:
            tr.set_codes(
                channel={'x': 'N', 'y': 'E', 'z': 'D'}[tr.channel])
            tr.shift(-t0)

        trace.snuffle(trs + trs2)


if __name__ == '__main__':
    util.setup_logging('test_ahfull', 'warning')
    unittest.main()
