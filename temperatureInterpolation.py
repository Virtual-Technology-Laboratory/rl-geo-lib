from __future__ import print_function

 # Copyright (c) 2015, Roger Lew (rogerlew.gmail.com)
 # Date: 6/16/2015
 # License: BSD (3-clause license)
 # 
 # The project described was supported by NSF award number IIA-1301792
 # from the NSF Idaho EPSCoR Program and by the National Science Foundation.
  

"""
This module provides functionality for reading daily temperature timeseries
obtained from the NOAA National Climatic Data Center (NCDC).

https://www.ncdc.noaa.gov/data-access/land-based-station-data

The module interpolates daily minimum and maximum temperature values using the
same algorithm as the interp.T package for R [1]. Some of the assumptions that
weren't explicitly listed in their paper [2] were taken from the original 
algorithm [3].

References:
[1] Eccel, E. & Cordano, C. (2013). Interpol.T: Hourly interpolation of multiple 
    temperature daily series.
    http://cran.r-project.org/web/packages/Interpol.T/index.html

[2] Eccel, E. (2010). What we can ask to hourly temperature recording. Part II:
    Hourly interpolation of temperatures for climatology and modelling. 
    http://www.agrometeorologia.it/documenti/Rivista2010_2/AIAM%202-2010_pag45.pdf

[3] Cesaraccio, C., Spano, D., Duce, P., Snyder, R. L. (2001). An improved model
    for determining degree-day values from daily temperature data. Int. J. 
    Biometeorol, 45, 161-169.
    http://download.springer.com/static/pdf/789/art%253A10.1007%252Fs004840100104.pdf?originUrl=http%3A%2F%2Flink.springer.com%2Farticle%2F10.1007%2Fs004840100104&token2=exp=1434487094~acl=%2Fstatic%2Fpdf%2F789%2Fart%25253A10.1007%25252Fs004840100104.pdf%3ForiginUrl%3Dhttp%253A%252F%252Flink.springer.com%252Farticle%252F10.1007%252Fs004840100104*~hmac=efc7b39f531ab587adda01b525aa189cf9b449a49b2a50b4f1b48d811cc72cc9
"""

import csv
from collections import OrderedDict
from datetime import datetime, timedelta
from math import pi, sin, sqrt

import astral
import pytz

def total_hours(tdelta):
    return tdelta.total_seconds() / 3600.0

_24hr = timedelta(hours=24)

class DailyTemperature:
    def __init__(self, tmin, tmax, date, solar,
                 c=0.39, z=0.5, HxHsInterval=4, HnHsInterval=1):
        """
        Represents daily temperature values and performs temperature 
        interpolation.
        
        Parameters
        ----------
        
        tmin : float
            the minimum temperature for the day
            
        tmax : float
            the maximum temperature for the day
            
        date : datetime
            datetime object representing the current day (at 00:00) in localtime
            
        solar : dict
            dictionary containing sunrise and sunset times (obtained from Astral)
            
        c : float (optional, default=0.39)
            Used to estimate the temperature at sunset. The 0.39 is a generic
            value provided by Cesaraccio, Spano,  Duce, & Snyder (2001).
            
        z : float (optional, default=0.5)
            during the sunset to sunrise periods used to parameterize whether
            the sky is clear (0.5) or cloudy (1.0). This was introduced into
            the revised model by Eccel (2010).
        """
        self.tmin = tmin
        self.tmax = tmax
        self.date = date
        self.solar = solar
        self.c = c
        self.z = z # clear day -> 0.5, cloudy -> 1
        self.HxHsInterval = HxHsInterval
        self.HnHsInterval = HnHsInterval
        self.nextDay = None
        self.prevDay = None
    
    def _Tp_getter(self):
        """
        returns the next days minimum temperature
        """
        if (self.nextDay == None):
            return None
        
        return self.nextDay.tmin

    Tp = property(_Tp_getter)
    
    def _Ts_getter(self):
        """
        returns estimate of the temperature at sunset
        """
        return self.tmax - self.c * (self.tmax - self.tmin)
    
    Ts = property(_Ts_getter)
    
    def _Hx_getter(self):
        """
        returns estimated hottest time of day
        """
        return self.solar['sunset'] - timedelta(hours=self.HxHsInterval)
    
    Hx = property(_Hx_getter)
    
    def _Hn_getter(self):
        """
        returns estimated hottest time of day
        """
        return self.solar['sunrise'] - timedelta(hours=self.HnHsInterval)
    
    Hn = property(_Hn_getter)
    
    def _Hp_getter(self):
        """
        return the time of the next day's sunrise 
        """
        if (self.nextDay == None):
            return None
        
        return self.nextDay.solar['sunrise']
        
    Hp = property(_Hp_getter)
        
    def _H0_getter(self):
        """
        returns the starting time of the represented day
        """
        return datetime(date.year, date.month, date.day)
        
    H0 = property(_Hp_getter)
    
    def _Hend_getter(self):
        """
        returns the starting time of the next day
        """
        global _24hr
        return H0 + _24hr
        
    Hend = property(_Hp_getter)
    
    def _Tn_getter(self):
        """
        returns the day's minimum temperature
        """
        return self.tmin
        
    Tn = property(_Tn_getter)
    
    def _Tx_getter(self):
        """
        returns the day's maximum temperature
        """
        return self.tmax
        
    Tx = property(_Tx_getter)
    
    def __call__(self, h):
        global _24hr
        
        c = self.c
        solar = self.solar
        HxHsInterval = self.HxHsInterval
        
        Tn = self.Tn
        Tx = self.tmax
        Tp = self.Tp
        Ts = self.Ts
            
        H0 = self.H0
        Hend = self.Hend
        
        Hn = self.Hn #solar['sunrise']
        Hs = solar['sunset']
        Hx = self.Hx #Hs - timedelta(hours=HxHsInterval)
        Hp = self.Hp #self.nextDay.solar['sunrise']
        z =  self.z 
        
            
        if (Hn < h and h <= Hx): # sunrise to max temp (4 hours before sunset)
            t_frac = total_hours(h - Hn) / total_hours(Hx - Hn)
            return Tn + ((Tx - Tn) / 2.0) * (1 + sin(pi * t_frac - (pi/2)))
        
        elif (Hx < h and h <= Hs): # max temp to sunset
            t_frac = total_hours(h - Hx) / total_hours(Hs - Hx)
            return Ts + (Tx - Ts) * sin((pi/2) * (1 + t_frac))
        
        elif (Hs < h and h <= Hend): # sunset to midnight
            Tn_next = self.nextDay.Tn
            deltaII = (Tn_next - Ts) / total_hours(Hn + _24hr - Hs) ** z 
            return Ts + deltaII * total_hours(h - Hs) ** z
        
        else: # midnight previous day to sunrise
            Ts_prev = self.prevDay.Ts
            deltaI = (Tn - Ts_prev) / total_hours(Hn + _24hr - Hs) ** z
            return Ts_prev + deltaI * total_hours(h + _24hr - Hs) **z
        
    def __repr__(self):
        return 'DailyTemperature(%f, %f, %s, %s)' % \
        (self.tmin, self.tmax, repr(self.date), repr(self.solar))
    
class DailyTemperatureSeries:

    def __init__(self, fname, astral_location):
        series = OrderedDict()

        with open(fname) as csvfile:
            reader = csv.DictReader(csvfile)

            last = None
            for row in reader:
                date_str = row['DATE']
                year = int(date_str[:4])
                month = int(date_str[4:6])
                day = int(date_str[6:])
                date = datetime(year, month, day)

                tmin = float(row['TMIN']) * 0.1
                tmax = float(row['TMAX']) * 0.1

                solar = location.sun(date)
                series[date] = DailyTemperature(tmin, tmax, date, solar)
                
                if (len(series) > 1):
                    series[last].nextDay = series[date]
                    series[date].prevDay = series[last]
                    
                last = datetime(year, month, day)

        self.series = series
        self.location = astral_location
        
    def datekey(self, h):
        return datetime(h.year, h.month, h.day)
        
    def __call__(self, h):
        date = self.datekey(h) 
        
        assert date in self.series     
        return self.series[date](h)
        
if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt

                               # name, region, latitude, longitude, timezone, elevation
    location = astral.Location(['Kenai', 'Alaska', 60.567402, -151.246719, 'US/Alaska', 26])
    dailyTemperatureSeries = DailyTemperatureSeries('KenaiTempData.csv', location)
    
    
    t0 = datetime(1900, 1, 12, 9, 23, 55, tzinfo=pytz.timezone('US/Alaska'))
    dates = [t0 + timedelta(days=x) for x in np.linspace(0,10,1000)]
    y = [dailyTemperatureSeries(date) for date in dates]
    
    
    plt.figure(figsize=(12,4))
    plt.plot(dates,y)

    for i in xrange(10):
        date = dailyTemperatureSeries.datekey(t0 + i * _24hr)
        
        sunrise = dailyTemperatureSeries.series[date].solar['sunrise']
        sunset = dailyTemperatureSeries.series[date].solar['sunset']
        plt.axvline(sunrise, c=(1, 0.5, 0), ls=':')
        plt.axvline(sunset,  c='k', ls=':')
        
    plt.show()