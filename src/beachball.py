from pyrocko import moment_tensor as mtm
from random import random, choice
from math import pi as PI
import numpy as num
import sys

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa


def vnorm(points):
    return num.sqrt(num.sum(points**2, axis=1))


def clean_poly(points):
    if not num.all(points[0, :] == points[-1, :]):
        points = num.vstack((points, points[0:1, :]))

    dupl = num.concatenate(
        (num.all(points[1:, :] == points[:-1, :], axis=1), [False]))
    points = points[num.logical_not(dupl)]
    return points


def circulation(points, axis):
    points2 = points[:, ((axis+2) % 3, (axis+1) % 3)]
    phi1 = num.arctan2(points2[:, 1], points2[:, 0])
    points2[:, 0] += num.cos(phi1) * (1.0-abs(points[:, axis]))
    points2[:, 1] += num.sin(phi1) * (1.0-abs(points[:, axis]))
    points2 = clean_poly(points2)
    vecs = points2[1:] - points2[:-1]
    vecs = num.vstack((vecs, vecs[0:1, :]))
    av = vecs[:-1, :]
    bv = vecs[1:, :]
    zs = num.cross(av, bv)
    phi = num.arcsin(zs / (vnorm(av) * vnorm(bv)))
    flip = num.sum(av*bv, axis=1) < 0.0
    phi[flip] = num.sign(phi[flip])*(PI - num.abs(phi[flip]))
    if num.any(phi == PI) or num.any(phi == -PI):
        sys.exit('ambiguous circulation!!!')

    result = num.sum(phi) / (2.0*PI)
    if round(result*100.) not in [100, -100]:
        sys.exit('circulation error')

    return result


def spoly_cut(l_points, axis=0, nonsimple=True):
    dphi = 2.*PI / 360.

    # cut sub-polygons and gather crossing point information
    crossings = []
    snippets = {}
    for ipath, points in enumerate(l_points):
        if not num.all(points[0, :] == points[-1, :]):
            points = num.vstack((points, points[0:1, :]))

        # get upward crossing points
        iup = num.where(num.logical_and(points[:-1, axis] <= 0.,
                                        points[1:, axis] > 0.))[0]
        aup = - points[iup, axis] / (points[iup+1, axis] - points[iup, axis])
        pup = points[iup, :] + aup[:, num.newaxis] * (points[iup+1, :] -
                                                      points[iup, :])
        phiup = num.arctan2(pup[:, (axis+2) % 3], pup[:, (axis+1) % 3])

        for i in xrange(len(iup)):
            crossings.append((phiup[i], ipath, iup[i], 1, pup[i], [1, -1]))

        # get downward crossing points
        idown = num.where(num.logical_and(points[:-1, axis] > 0.,
                                          points[1:, axis] <= 0.))[0]
        adown = - points[idown+1, axis] / (points[idown, axis] -
                                           points[idown+1, axis])
        pdown = points[idown+1, :] + adown[:, num.newaxis] * (
            points[idown, :] - points[idown+1, :])
        phidown = num.arctan2(pdown[:, (axis+2) % 3], pdown[:, (axis+1) % 3])

        for i in xrange(idown.size):
            crossings.append(
                (phidown[i], ipath, idown[i], -1, pdown[i], [1, -1]))

        icuts = num.sort(num.concatenate((iup, idown)))

        for i in xrange(icuts.size-1):
            snippets[ipath, icuts[i]] = (
                ipath, icuts[i+1], points[icuts[i]+1:icuts[i+1]+1])

        if icuts.size:
            points_last = num.concatenate((
                points[icuts[-1]+1:],
                points[:icuts[0]+1]))

            snippets[ipath, icuts[-1]] = (ipath, icuts[0], points_last)
        else:
            snippets[ipath, 0] = (ipath, 0, points)

    crossings.sort()

    # assemble new sub-polygons
    current = snippets.pop(snippets.keys()[0])
    outs = [[]]
    while True:
        outs[-1].append(current[2])
        for i, c1 in enumerate(crossings):
            if c1[1:3] == current[:2]:
                direction = -1 * c1[3]
                break
        else:
            if not snippets:
                break
            current = snippets.pop(snippets.keys()[0])
            outs.append([])
            continue

        while True:
            i = (i + direction) % len(crossings)
            if crossings[i][3] == direction and direction in crossings[i][-1]:
                break

        c2 = crossings[i]
        c2[-1].remove(direction)

        phi1 = c1[0]
        phi2 = c2[0]
        if direction == 1:
            if phi1 > phi2:
                phi2 += PI * 2.

        if direction == -1:
            if phi1 < phi2:
                phi2 -= PI * 2.

        n = int(abs(phi2 - phi1) / dphi) + 2

        phis = num.linspace(phi1, phi2, n)
        cpoints = num.zeros((n, 3))
        cpoints[:, (axis+1) % 3] = num.cos(phis)
        cpoints[:, (axis+2) % 3] = num.sin(phis)
        cpoints[:, axis] = 0.0

        outs[-1].append(cpoints)

        try:
            current = snippets[c2[1:3]]
            del snippets[c2[1:3]]

        except KeyError:
            if not snippets:
                break

            current = snippets.pop(snippets.keys()[0])
            outs.append([])

    # separate hemispheres, force polygons closed, remove duplicate points
    # remove polygons with less than 3 points (4, when counting repeated
    # endpoint)

    outs_upper = []
    outs_lower = []
    for out in outs:
        if out:
            out = clean_poly(num.vstack(out))
            if out.shape[0] >= 4:
                if num.sum(out[:, axis]) > 0.0:
                    outs_upper.append(out)
                else:
                    outs_lower.append(out)

    if nonsimple and (len(outs_upper) == 0 or len(outs_lower) == 0):
        # check if we are cutting between holes
        need_divider = False
        if outs_upper:
            candi = sorted(
                outs_upper, key=lambda out: num.min(out[:, axis]))[0]

            if circulation(candi, axis) > 0.0:
                need_divider = True

        if outs_lower:
            candi = sorted(
                outs_lower, key=lambda out: num.max(out[:, axis]))[0]

            if circulation(candi, axis) < 0.0:
                need_divider = True

        if need_divider:
            phi1 = 0.
            phi2 = PI*2.
            n = int(abs(phi2 - phi1) / dphi) + 2

            phis = num.linspace(phi1, phi2, n)
            cpoints = num.zeros((n, 3))
            cpoints[:, (axis+1) % 3] = num.cos(phis)
            cpoints[:, (axis+2) % 3] = num.sin(phis)
            cpoints[:, axis] = 0.0

            outs_upper.append(cpoints)
            outs_lower.append(cpoints[::-1, :])

    return outs_lower, outs_upper


def fuzz(mi, ma):

    vals = [
        0., 1., PI, 2.*PI, 0.5*PI/2., 0.25*PI,
        mi+random()*(ma-mi)*1.0e-8,
        ma-random()*(ma-mi)*1.0e-8,
        mi+random()*(ma-mi)*1.0e-4,
        ma-random()*(ma-mi)*1.0e-4,
        mi, ma,
        mi+random()*(ma-mi)]

    while True:
        v = choice(vals)
        if mi <= v and v <= ma:
            return v


def numpy_rtp2xyz(rtp):
    r = rtp[:, 0]
    theta = rtp[:, 1]
    phi = rtp[:, 2]
    vecs = num.empty(rtp.shape, dtype=num.float)
    vecs[:, 0] = r*num.sin(theta)*num.cos(phi)
    vecs[:, 1] = r*num.sin(theta)*num.sin(phi)
    vecs[:, 2] = r*num.cos(theta)
    return vecs


def numpy_xyz2rtp(xyz):
    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    vecs = num.empty(xyz.shape, dtype=num.float)
    vecs[:, 0] = num.sqrt(x**2+y**2+z**2)
    vecs[:, 1] = num.arccos(z/vecs[:, 0])
    vecs[:, 2] = num.arctan2(y/vecs[:, 0], x/vecs[:, 0])
    return vecs


def eig2gx(eig):
    aphi = num.linspace(0., 2.*PI, 181)
    ep, en, et, vp, vn, vt= eig
    groups = []
    for (pt_name, pt_sign) in [('P', -1.), ('T', 1.)]:
        patches = []
        patches_lower = []
        patches_upper = []
        lines = []
        lines_lower = []
        lines_upper = []
        for (va, vb, vc, ea, eb, ec) in [
                (vp, vn, vt, ep, en, et),
                (vt, vp, vn, et, ep, en),
                (vn, vt, vp, en, et, ep)]:

            to_e = num.vstack((vb, vc, va))
            from_e = to_e.T

            poly_es = []
            polys = []
            for sign in (-1., 1.):
                xphi = pt_sign*sign*aphi
                denom = eb*num.cos(xphi)**2 + ec*num.sin(xphi)**2
                if num.any(denom == 0.):
                    continue

                Y = -ea/denom
                if num.any(Y < 0.):
                    continue

                xtheta = num.arctan(num.sqrt(Y))
                rtp = num.empty(xphi.shape+(3,), dtype=num.float)
                rtp[:, 0] = 1.
                if sign > 0:
                    rtp[:, 1] = xtheta
                else:
                    rtp[:, 1] = PI - xtheta

                rtp[:, 2] = xphi
                poly_e = numpy_rtp2xyz(rtp)
                poly = num.dot(from_e, poly_e.T).T

                poly_es.append(poly_e)
                polys.append(poly)

            if poly_es:
                for aa in spoly_cut(poly_es, 0):
                    for bb in spoly_cut(aa, 1):
                        for cc in spoly_cut(bb, 2):
                            for poly_e in cc:
                                poly = num.dot(from_e, poly_e.T).T
                                polys_lower, polys_upper = spoly_cut(
                                    [poly], 2, nonsimple=False)

                                patches.append(poly)
                                patches_lower.extend(polys_lower)
                                patches_upper.extend(polys_upper)

        groups.append((
            pt_name,
            patches, patches_lower, patches_upper,
            lines, lines_lower, lines_upper))

    return groups


def extr(points):
    pmean = num.mean(points, axis=0)
    return points + pmean*0.05


def plot_beachball_mpl(mt, axes):
    mt_devi = mt.deviatoric()
    eig = mt_devi.eigensystem()

    for (group, patches, patches_lower, patches_upper,
            lines, lines_lower, lines_upper) in eig2gx(eig):

        if group == 'P':
            color = 'red'
        else:
            color = 'blue'

        for poly in patches_lower:
            px, py, pz = poly.T
            axes.fill(py, px, alpha=0.4, lw=1, color=color)

        for poly in lines_lower:
            px, py, pz = poly.T
            axes.plot(py, px, color='black')


if __name__ == '__main__':
    plt.ion()
    plt.show()
    fig = plt.figure()
    axes = fig.add_subplot(1, 1, 1)

    for x in range(100):
        m = mtm.symmat6(*(num.random.random(6)*2.-1.))
        mt = mtm.MomentTensor(m=m)
        plot_beachball_mpl(mt, axes)
        fig.canvas.draw()
