'''
Created on 05.12.2015

@author: Juri Bieler
'''

import time
import sys

class Tee(object):
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f in self.files:
            f.write(obj)


class Logger(object):

    def __init__(self):
        self.file = open('log/cout_' + str(time.strftime("%Y-%m-%d_%H-%M-%S")) + '.txt', 'w')
        self.original = sys.stdout
        sys.stdout = Tee(sys.stdout, self.file)

            
    def terminate(self):
        sys.stdout = self.original
        self.file.close()