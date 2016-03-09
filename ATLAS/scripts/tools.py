import math

"""Generates a formatter function for rounding a float to "digits" decimal
places"""
def round(digits):
    fstring = "{:." + str(digits) + "f}"
    
    def ret(x):
        if math.isnan(x):
            return ""
        else:
            return fstring.format(x)
    
    return ret
