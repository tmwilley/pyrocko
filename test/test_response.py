import unittest, os
import numpy as num
from pyrocko import util, evalresp, pz, trace, guts

import common

def numeq(a,b, eps):
    return num.all(num.abs(num.array(a) - num.array(b)) < eps)

class ResponseTestCase(unittest.TestCase):
    
    def test_evalresp(self, plot=False):

        resp_fpath = common.test_data_file('test2.resp')

        freqs = num.logspace(num.log10(0.001), num.log10(10.), num=1000)

        transfer = evalresp.evalresp(sta_list='BSEG',
                          cha_list='BHZ',
                          net_code='GR',
                          locid='',
                          instant=util.str_to_time('2012-01-01 00:00:00'),
                          freqs=freqs,
                          units='DIS',
                          file=resp_fpath,
                          rtype='CS')[0][4]

        pz_fpath = common.test_data_file('test2.sacpz')

        zeros, poles, constant = pz.read_sac_zpk(pz_fpath)
        
        resp = trace.PoleZeroResponse(zeros, poles, constant)

        transfer2 = resp.evaluate(freqs)

        if plot:
            import pylab as lab
            lab.plot(freqs, num.imag(transfer))
            lab.plot(freqs, num.imag(transfer2))
            lab.gca().loglog() 
            lab.show()

        assert numeq(transfer, transfer2, 1e-4)

    def test_dump_load(self):
        
        r = trace.FrequencyResponse()

        r = trace.PoleZeroResponse([0j, 0j], [1j, 2j, 1+3j, 1-3j], 1.0)
        r.regularize()
        r2 = guts.load_string(r.dump())
        assert cnumeq(r.poles, r2.poles, 1e-6)
        assert cnumeq(r.zeros, r2.zeros, 1e-6)
        assert numeq(r.constant, r2.constant)

        r = trace.SampledResponse(
                [0., 1., 5., 10.],
                [0., 1., 1., 0.])

        r.regularize()
        r2 = guts.load_string(r.dump())
        assert numeq(r.frequencies, r2.frequencies, 1e-6)
        assert cnumeq(r.values, r2.values, 1e-6)

        r = trace.IntegrationResponse(2, 5.0)
        r2 = guts.load_string(r.dump())
        assert numeq(r.n, r2.n)
        assert numeq(r.gain, r2.gain, 1e-6)

        r = trace.DifferentiationResponse(2, 5.0)
        r2 = guts.load_string(r.dump())
        assert numeq(r.n, r2.n)
        assert numeq(r.gain, r2.gain, 1e-6)

        r = trace.AnalogFilterResponse(
            a=[1.0, 2.0, 3.0],
            b=[2.0, 3.0])
        r2 = guts.load_string(r.dump())
        assert numeq(r.a, r2.a, 1e-6)
        assert numeq(r.b, r2.b, 1e-6)


def numeq(a, b, eps=0):
    a = num.asarray(a)
    b = num.asarray(b)
    return num.max(num.abs(a - b)) <= eps


def cnumeq(a, b, eps=0):
    a = num.asarray(a)
    b = num.asarray(b)
    return num.max(num.abs(a.real - b.real)) <= eps and \
           num.max(num.abs(a.imag - b.imag)) <= eps


if __name__ == "__main__":
    util.setup_logging('test_response', 'warning')
    unittest.main()
