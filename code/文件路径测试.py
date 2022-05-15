import sys,os
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_DIR)
from module.module1 import tese_module1_1,tese_module1_2
from module.module2 import tese_module2_1,tese_module2_2

if __name__ == '__main__':
    tese_module1_1()
    tese_module1_2()
    tese_module2_1()
    tese_module2_2()
    print(os.path.dirname(os.path.abspath(__file__)))
