import unittest
import mpmath as mpm
from dbn_upper_bound.python.mputility import RSZ_plain


mpm.mp.dps = 64


class TestRiemannZetaMethods(unittest.TestCase):

    def test_rsz_plain(self):
        print("Testing test_rsz_plain")
        self.assertEqual(RSZ_plain(7),
                         mpm.ctx_mp_python.mpf('-1.864646804060364152877136431647678106608335283311309510218784337516'))
        self.assertEqual(RSZ_plain(70),
                         mpm.ctx_mp_python.mpf('0.7449143203472381607992454269461905832009655243043250990620264020776'))
        self.assertEqual(RSZ_plain(500),
                         mpm.ctx_mp_python.mpf('1.723691357652479072905353544031637533254727488320737100506859803423'))


if __name__ == '__main__':
    unittest.main()
