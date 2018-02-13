import unittest
import mpmath as mpm
from dbn_upper_bound.python.mputility import Ht_AFE_ABC


mpm.mp.dps = 64


class TestHtMethods(unittest.TestCase):

    def test_ht_afe_abc(self):
        print("Testing test_ht_afe_abc")
        self.assertEqual(Ht_AFE_ABC(35+10j, 0),
                         mpm.ctx_mp_python.mpc(real='0.000313552332336939669577229294629866377113376343853824475809041362996',
                                               imag='0.00007605623264065672294567401329292513426458248255176130198721233252061'))
        self.assertEqual(Ht_AFE_ABC(35+10j, 1),
                         mpm.ctx_mp_python.mpc(real='0.0003251967110350425054566521794876102859158315098752036658198040231417',
                                               imag='0.00002698115887366512726456719545317965659234388554917279737689136923695'))

        self.assertEqual(Ht_AFE_ABC(1000000.0, 0.23),
                         mpm.ctx_mp_python.mpf('-2.07911712735971640170536766390318417562939296431867628827775220872e-170537'))
        
        self.assertEqual(Ht_AFE_ABC(3000000.0, 0.23),
                         mpm.ctx_mp_python.mpf('7.754307373196585578741651057258108684598351023017606312425025045854e-511631'))
        
if __name__ == '__main__':
    unittest.main()
