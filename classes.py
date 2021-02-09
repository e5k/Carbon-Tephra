# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 14:40:23 2021

@author: Ben
"""
import utm


class ZONE:
    def __init__(self, extended, num, let, orientation=None):
        if extended:
            if orientation == "N":
                self.num = num
                self.let = chr(ord(let) + 1)
            elif orientation == "S":
                self.num = num
                self.let = chr(ord(let) - 1)
            elif orientation == "E":
                self.let = let
                if not num == 1:
                    self.num = num + 1
                else:
                    self.num = 60
            elif orientation == "W":
                if not num == 60:
                    self.num = num - 1
                else:
                    self.num = 1
            
            self.name = chr(self.num) + self.let
        else:
            self.num = num
            self.let = let
            self.name = chr(self.num) + self.let
            
        self.years = []
        self.VEIs = []
        self.carbon = 0
        self.surface = 0
        
        
class VOLCANO:
    def __init__(self, name, lat, lon):
        coords = utm.from_latlon(lat,lon)
        self.zone_num = coords[2]
        self.zone_let = coords[3]
        self.easting = coords[0]
        self.northing = coords[1]
        proba = {}
        