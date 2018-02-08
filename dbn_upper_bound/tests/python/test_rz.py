import unittest
from dbn_upper_bound.python.mputility import RSZ_plain, Ht_AFE_ABC
import mpmath as mpm


class TestRiemannZetaMethods(unittest.TestCase):

    def test_rsz_plain(self):
        print("Testing test_rsz_plain")
        self.assertEqual(RSZ_plain(7), mpm.ctx_mp_python.mpf('-1.86464680406036415287713643164771'))
        self.assertEqual(RSZ_plain(70), mpm.ctx_mp_python.mpf('0.74491432034723816079924542693526'))
        self.assertEqual(RSZ_plain(500), mpm.ctx_mp_python.mpf('1.72369135765247907290535354444898'))

    def test_ht_afe_abc(self):
        print("Testing test_ht_afe_abc")
        self.assertEqual(Ht_AFE_ABC(35+10j, 0),
                         mpm.ctx_mp_python.mpc(real='0.000313552332336939669577229294629972',
                                               imag='0.000076056232640656722945674013292791'))
        self.assertEqual(Ht_AFE_ABC(35+10j, 1),
                         mpm.ctx_mp_python.mpc(real='0.000325196711035042505456652179487761',
                                               imag='0.0000269811588736651272645671954530276'))


if __name__ == '__main__':
    unittest.main()