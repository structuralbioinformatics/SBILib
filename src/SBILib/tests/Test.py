import unittest
from test_modules import * 


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(Test_structure))
    suite.addTest(unittest.makeSuite(Test_Blast))
    suite.addTest(unittest.makeSuite(Test_Complex))
    suite.addTest(unittest.makeSuite(Test_Loops))
    suite.addTest(unittest.makeSuite(Test_Grafting))
    runner = unittest.TextTestRunner()
    runner.run(suite)
