import math
import matplotlib.ticker as tkr
import matplotlib.pyplot as plt

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

"""Latex comma formatter for plotting"""
def comma(x, pos):
    s = str(x)
    ind = s.index(".")
    return s[:ind] + "{,}" + s[ind+1:]

"""Set decimal comma"""
def set_mpl_comma():
    ax = plt.gca()
    formatter = tkr.FuncFormatter(comma)
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)
