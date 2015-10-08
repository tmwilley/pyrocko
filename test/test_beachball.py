import unittest
import numpy as num
from cStringIO import StringIO as StringIO

from pyrocko import util, moment_tensor as mtm, beachball


class BeachballTestCase(unittest.TestCase):

    def test_beachball(self):
        from matplotlib import pyplot as plt
        from matplotlib import image

        nx = 10

        for x in range(nx):
            m6 = num.random.random(6)*2.-1.
            m = mtm.symmat6(*m6)

            mt = mtm.MomentTensor(m=m)
            mt = mt.deviatoric()

            imgs = []
            for iplot, plot in enumerate((beachball.plot_beachball_mpl, 
                    beachball.plot_beachball_mpl_pixmap)):

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
            d[:,:,3] = 1.
            dsum = num.sum(d[:, :, :3])
            print dsum
            if dsum > 1000:
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
