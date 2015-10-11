import unittest
import numpy as num
from cStringIO import StringIO as StringIO

from pyrocko import util, moment_tensor as mtm, beachball

from math import pi as PI
from random import random, choice


def fuzz_angle(mi, ma):
    crits = [mi, ma, 0., 1., -1., 90., -90., 180., -180., 270., -270.]
    crits = [x for x in crits if mi <= x and x <= ma]
    ran = [
        0.,
        0.,
        0.,
        +random()*(ma-mi)*1.0e-8,
        -random()*(ma-mi)*1.0e-8,
        +random()*(ma-mi)*1.0e-4,
        -random()*(ma-mi)*1.0e-4,
        random()*(ma-mi),
        -random()*(ma-mi)]

    while True:
        v = choice(crits) + choice(ran)
        if mi <= v and v <= ma:
            return v


class BeachballTestCase(unittest.TestCase):

    def test_beachball(self):
        from matplotlib import pyplot as plt
        from matplotlib import image

        nx = 100

        for x in range(nx):
            # m6 = num.random.random(6)*2.-1.
            # m6 = [0., 0., 0., 0., 0., 1.]
            # m = mtm.symmat6(*m6)

            # mt = mtm.MomentTensor(m=m)
            strike = fuzz_angle(0., 360.)
            dip = fuzz_angle(0., 90.)
            rake = fuzz_angle(-180., 180.)

            #strike = 270.
            #dip = 0.0
            #rake = 0.01

            mt = mtm.MomentTensor(
                strike=strike,
                dip=dip,
                rake=rake)

            # mt = mt.deviatoric()

            imgs = []
            for iplot, plot in enumerate([
                    beachball.plot_beachball_mpl,
                    beachball.plot_beachball_mpl_pixmap]):

                fig = plt.figure(figsize=(3, 3), dpi=100)
                axes = fig.add_subplot(1, 1, 1, aspect=1.)
                axes.cla()
                axes.axison = False
                axes.set_xlim(-1.05, 1.05)
                axes.set_ylim(-1.05, 1.05)

                plot(mt, axes)

                f = StringIO()
                fig.savefig(f, format='png')
                f.seek(0)
                imgs.append(image.imread(f, format='png'))
                fig.clear()
                plt.close(fig)

            a, b = imgs

            d = num.abs(a-b)
            d[:, :, 3] = 1.
            dsum = num.sum(d[:, :, :3])
            print dsum
            if nx == 1 or dsum > 1500:
                print repr(strike), repr(dip), repr(rake)
                print mt
                fig = plt.figure()
                axes1 = fig.add_subplot(1, 3, 1, aspect=1.)
                axes2 = fig.add_subplot(1, 3, 2, aspect=1.)
                axes3 = fig.add_subplot(1, 3, 3, aspect=1.)
                axes1.imshow(a)
                axes2.imshow(b)
                axes3.imshow(d)
                plt.show()
                plt.close(fig)


if __name__ == "__main__":
    util.setup_logging('test_moment_tensor', 'warning')
    unittest.main()
