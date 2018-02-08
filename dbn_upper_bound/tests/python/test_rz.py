import unittest
from dbn_upper_bound.python.mputility import RSZ_plain
import mpmath as mpm


class TestRiemannZetaMethods(unittest.TestCase):

    def test_rsz_plain(self):
        self.assertEqual(RSZ_plain(7), mpm.ctx_mp_python.mpf('-1.86464680406036415287713643164771'))
        self.assertEqual(RSZ_plain(70), mpm.ctx_mp_python.mpf('0.74491432034723816079924542693526'))
        self.assertEqual(RSZ_plain(500), mpm.ctx_mp_python.mpf('1.72369135765247907290535354444898'))


if __name__ == '__main__':
    unittest.main()