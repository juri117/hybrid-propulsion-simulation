'''
Created on 24.12.2015

@author: Juri Bieler
'''

#from ConfigParser import SafeConfigParser
#from backports import configparser
import configparser
import generalUtils as generalUtils
import const as c

class Ini:

    ### CONSTRUCKTOR #####################################
    '''
    liest config.ini datei aus und speichert Werte in dict
    '''
    def __init__(self, mission, resPath=""):
        self.path = generalUtils.getMissionPath(mission) + resPath + "config.ini"
        self.path = generalUtils.fixPath(self.path)
        self.parser = configparser.SafeConfigParser()
        self.parser.read(self.path)

    def get(self, cat, key):
        if(self.keyExists(cat, key)):
            return self.parser.get(cat, key)
        return float("nan")
    
    def keyExists(self, cat, key):
        return self.parser.has_option(cat, key)
    
    def set(self, cat, key, val):
        if(not self.parser.has_section(cat)):
            self.parser.add_section(cat)
        
        self.parser.set(cat, key, str(val))
        # Writing our configuration file to 'example.cfg'
        with open(self.path, 'wt') as configfile:
            self.parser.write(configfile)
            
    def copyConfig(self, path):
        with open(path, 'wt') as configfile:
            self.parser.write(configfile)
            
    def getConfig(self, propulsionType):
        #keys = self.parser.get_keywords(propulsionType)
        config = dict(self.parser.items(propulsionType))
        for key, value in iter(config.items()):
            config[key] = float(config[key])
        config["mpl"] = c.M_PL
        if(self.keyExists("general", "mpl")):
            config["mpl"] = float(self.ini.get("general", "mpl"))
        return config