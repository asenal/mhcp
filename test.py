import unittest,os,sys
here=os.popen('pwd').readline().strip()
sys.path.append(here+'libs/libsvm-3.12/python')
sys.path.append(here+'libs/')
import hla

class HlaTestCase(unittest.TestCase):
    def test 
