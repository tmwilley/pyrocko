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
        raise Exception('ambiguous circulation')

    result = num.sum(phi) / (2.0*PI)
    if int(round(result*100.)) not in [100, -100]:
        raise Exception('circulation error')

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
    vecs[:, 1] = num.arctan2(num.sqrt(x**2+y**2), z)
    vecs[:, 2] = num.arctan2(y, x)
    return vecs


def eig2gx(eig):
    aphi = num.linspace(0., 2.*PI, 181)
    ep, en, et, vp, vn, vt = eig
    groups = []
    for (pt_name, pt_sign) in [('P', -1.), ('T', 1.)]:
        patches = []
        patches_lower = []
        patches_upper = []
        lines = []
        lines_lower = []
        lines_upper = []
        for iperm, (va, vb, vc, ea, eb, ec) in enumerate([
                (vp, vn, vt, ep, en, et),
                (vt, vp, vn, et, ep, en)]):
                #(vn, vt, vp, en, et, ep)]):

            perm_sign = [-1.0, 1.0][iperm]
            to_e = num.vstack((vb, vc, va))
            from_e = to_e.T

            poly_es = []
            polys = []
            for sign in (-1., 1.):
                xphi = perm_sign*pt_sign*sign*aphi
                denom = eb*num.cos(xphi)**2 + ec*num.sin(xphi)**2
                if num.any(denom == 0.):
                    continue

                Y = -ea/denom
                if num.any(Y < 0.):
                    continue

                print iperm, pt_name, sign

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

            if polys:
                polys_lower, polys_upper = spoly_cut(polys, 2)
                lines.extend(polys)
                lines_lower.extend(polys_lower)
                lines_upper.extend(polys_upper)

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
    mt_devi = mt#.deviatoric()
    eig = mt_devi.eigensystem()

    for (group, patches, patches_lower, patches_upper,
            lines, lines_lower, lines_upper) in eig2gx(eig):

        if group == 'P':
            color = 'white'
        else:
            color = 'red'

        for poly in patches_lower:
            px, py, pz = poly.T
            axes.fill(-py, -px, lw=1, color=color, fc=color)

        for poly in lines_lower:
            px, py, pz = poly.T
            axes.plot(-py, -px, lw=2, color='black')


def plot_beachball_mpl_pixmap(mt, axes):
    mt_devi = mt#.deviatoric()
    ep, en, et, vp, vn, vt = mt_devi.eigensystem()

    nx = 400
    ny = 400

    x = num.linspace(-1., 1., nx)
    y = num.linspace(-1., 1., ny)

    vecs = num.zeros((nx*ny, 3), dtype=num.float)
    vecs[:, 0] = num.tile(x, ny)
    vecs[:, 1] = num.repeat(y, nx)
    ii_ok = vecs[:, 0]**2 + vecs[:, 1]**2 <= 1.0
    vecs[ii_ok, 2] = num.sqrt(1.0 - (vecs[ii_ok, 0]**2 + vecs[ii_ok, 1]**2))
    vecs_ok = vecs[ii_ok, :]
    to_e = num.vstack((vn, vt, vp))

    vecs_e = num.dot(to_e, vecs_ok.T).T
    rtp = numpy_xyz2rtp(vecs_e)

    atheta, aphi = rtp[:, 1], rtp[:, 2]
    amps_ok = ep * num.cos(atheta)**2 + (
        en * num.cos(aphi)**2 + et * num.sin(aphi)**2) * num.sin(atheta)**2

    amps = num.zeros(nx*ny, dtype=num.float)
    amps[:] = num.nan
    amps[ii_ok] = amps_ok

    amps = num.reshape(amps, (ny, nx))


    axes.contourf(y, x, amps.T, levels=[-2., 0., 2.], colors=['white', 'red'])
    axes.contour(y, x, amps.T, levels=[0.], colors=['black'], linewidths=2)

    phi = num.linspace(0., 2*PI, 361)
    x = num.cos(phi)
    y = num.sin(phi)
    axes.plot(x, y, lw=2, color='black')


if __name__ == '__main__':
    nx = 100
    if nx > 1:
        plt.ion()
        plt.show()

    fig = plt.figure()
    axes1 = fig.add_subplot(1, 3, 1, aspect=1.)
    axes2 = fig.add_subplot(1, 3, 2, aspect=1.)
    axes3 = fig.add_subplot(1, 3, 3, aspect=1.)

    import mopad
    import time

    for x in range(nx):
        m6 = num.random.random(6)*2.-1.
        m = mtm.symmat6(*m6)

        mt = mtm.MomentTensor(m=m)
        mt = mt.deviatoric()
        for axes in (axes1, axes2, axes3):
            axes.cla()
            axes.axison = False
            axes.set_xlim(-1.05, 1.05)
            axes.set_ylim(-1.05, 1.05)

        axes1.set_title('Copacabana')
        axes2.set_title('Contour')
        axes3.set_title('MoPaD')

        plot_beachball_mpl(mt, axes1)
        plot_beachball_mpl_pixmap(mt, axes2)

        mop_mt = mopad.MomentTensor(M=mt.m6())
        mop_beach = mopad.BeachBall(mop_mt)
        kwargs = dict(
            plot_projection='ortho',
            plot_nodalline_width=2,
            plot_faultplane_width=2,
            plot_outerline_width=2)

        mop_beach.ploBB(kwargs, ax=axes3)

        fig.canvas.draw()

    if nx == 1:
        plt.show()
