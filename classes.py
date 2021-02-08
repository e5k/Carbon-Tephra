# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 14:40:23 2021

@author: Ben
"""

class ZONE:
    def __init__(self, new, num, let, orientation=None):
        if new:
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
            
        