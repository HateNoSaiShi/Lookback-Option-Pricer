# -*- coding: utf-8 -*-
"""
Created on Wed Aug 02 10:23:31 2017

@author: HateNoSaiSHI
"""

import numpy as np
from scipy.stats import norm

class LookbackOption(object):
    
    def __init__(self, spot, r, sigma, T, OptionType, **kwargs):
        """
        Attributes
        ----------
        Smin/Smax:
            Running min/max
        
        T:
            Time to maturity. In year conventiom.
            
        OptionType:
            String only.
            Invalid input:
                FloatingStrike_call
                FloatingStrike_put
                FixedStrike_call
                FixedStrike_put
                PartialTimeFloatingStrike_call
                PartialTimeFloatingStrike_put
                PartialTimeFixedStrike_call
                PartialTimeFixedStrike_put
                
        kwargs:
            Dictionary. Inputs should be consistent with option type.
            Invalid key value:
                Smax
                Smin
                strike
                partial_T
                dividend
        """
        
        self.spot = float(spot)
        self.r = float(r)
        self.sigma = float(sigma)
        self.T = float(T)
        self.OptionType = OptionType
        if 'dividend' in kwargs.keys():
            self.b = self.r - kwargs['dividend']
        else:
            self.b = self.r
        
        try:
            if OptionType == 'FloatingStrike_call':
                self.Smin = kwargs['Smin']
            elif OptionType == 'FloatingStrike_put':
                self.Smax = kwargs['Smax']
            elif OptionType == 'FixedStrike_call':
                self.Smax = kwargs['Smax']
                self.strike = kwargs['strike']
            elif OptionType == 'FixedStrike_put':
                self.Smin = kwargs['Smin']
                self.strike = kwargs['strike']
            else:
                self.Partial_T = kwargs['partial_T']
                if self.Partial_T >= self.T:
                    raise ValueError('Partial time should be smaller than T')
                if OptionType == 'PartialTimeFloatingStrike_call':
                    self.Smin = kwargs['Smin']
                elif OptionType == 'PartialTimeFloatingStrike_put':
                    self.Smax = kwargs['Smax']
                elif OptionType == 'PartialTimeFixedStrike_call':
                    self.strike = kwargs['strike']
                elif OptionType == 'PartialTimeFixedStrike_put':
                    self.strike = kwargs['strike']
                else:
                    raise TypeError('Invalid option type')
        except KeyError:
            raise KeyError('Invalid keys in kwargs.')
            

    def non_partial(self):
        spot = self.spot
        r = self.r
        sigma = self.sigma
        T = self.T
        b = self.b
        if self.OptionType == 'FloatingStrike_call':
            if not b == 0:
                a1 = (np.log(spot / self.Smin) + (b + sigma ** 2 / 2) * T) / sigma / np.sqrt(T)
                a2 = a1 - sigma * np.sqrt(T)
                a3 = (spot/self.Smin) ** (-2 * b / sigma ** 2) * norm.cdf(-a1 + 2 * b / sigma * np.sqrt(T)) - np.exp(b * T) * norm.cdf(-a1)
                c = spot * np.exp((b-r) * T) * norm.cdf(a1) - self.Smin * np.exp(-r * T) * norm.cdf(a2)
                c += spot * np.exp(-r * T) * sigma ** 2 / 2 / b * a3
                return c
            else:
                c = spot * np.exp((-r) * T) * norm.cdf(a1) - self.Smin * np.exp(-r * T) * norm.cdf(a2)
                c += spot * np.exp(-r * T) * sigma * np.sqrt(T) * (norm.pdf(a1) + a1 * (norm.cdf(a1) - 1))
                return c
        elif self.OptionType == 'FloatingStrike_put':
            if not b == 0:
                b1 = (np.log(spot / self.Smax) + (b + sigma ** 2 / 2) * T) / sigma / np.sqrt(T)
                b2 = a1 - sigma * np.sqrt(T)
                b3 = -(spot/self.Smin) ** (-2 * b / sigma ** 2) * norm.cdf(b1 - 2 * b / sigma * np.sqrt(T)) + np.exp(b * T) * norm.cdf(a1)
                p = -spot * np.exp((b-r) * T) * norm.cdf(-b1) + self.Smax * np.exp(-r * T) * norm.cdf(-b2)
                p += spot * np.exp(-r * T) * sigma ** 2 / 2 / b * b3
                return p
            else:
                p = -spot * np.exp((-r) * T) * norm.cdf(-b1) + self.Smax * np.exp(-r * T) * norm.cdf(-b2)
                p += spot * np.exp(-r * T) * sigma * np.sqrt(T) * (norm.pdf(b1) + b1 * norm.cdf(b1))
                return p
        elif self.OptionType == 'FixedStrike_call':
            if self.strike <= self.Smax:
                e1 = (np.log(spot / self.Smax) + (b + sigma ** 2 / 2) * T) / sigma / np.sqrt(T)
                e2 = e1 - sigma * np.sqrt(T)
                e3 = -(spot/self.Smax) ** (-2 * b / sigma ** 2) * norm.cdf(e1 - 2 * b / sigma * np.sqrt(T)) + np.exp(b * T) * norm.cdf(e1)
                c = np.exp(-r * T) * (self.Smax - self.strike) + spot * np.exp((b - r) * T) * norm.cdf(e1) - self.Smax * np.exp(-r * T) * norm.cdf(e2)
                c += spot * np.exp(-r * T) * sigma ** 2 / 2 / b * e3
                return c
            else:
                d1 = (np.log(spot / self.strike) + (b + sigma ** 2 / 2) * T) / sigma / np.sqrt(T)
                d2 = d1 - sigma * np.sqrt(T)
                d3 = -(spot/self.strike) ** (-2 * b / sigma ** 2) * norm.cdf(d1 - 2 * b / sigma * np.sqrt(T)) + np.exp(b * T) * norm.cdf(d1)
                c = spot * np.exp((b - r) * T) * norm.cdf(d1) - self.strike * np.exp(-r * T) * norm.cdf(d2)
                c += spot * np.exp(-r * T) * sigma ** 2 / 2 / b * (d3)
                return c
        elif self.OptionType == 'FixedStrike_put':
            if self.strike >= self.Smin:
                f1 = (np.log(spot / self.Smin) + (b + sigma ** 2 / 2) * T) / sigma / np.sqrt(T)
                f2 = f1 - sigma * np.sqrt(T)
                f3 = (spot/self.Smin) ** (-2 * b / sigma ** 2) * norm.cdf(-f1 + 2 * b / sigma * np.sqrt(T)) - np.exp(b * T) * norm.cdf(-f1)
                p = np.exp(-r * T) * (-self.Smin + self.strike) - spot * np.exp((b - r) * T) * norm.cdf(-f1) + self.Smin * np.exp(-r * T) * norm.cdf(-f2)
                p +=  spot * np.exp(-r * T) * sigma ** 2 / 2 / b * f3
                return p
            else:
                d1 = (np.log(spot / self.strike) + (b + sigma ** 2 / 2) * T) / sigma / np.sqrt(T)
                d2 = d1 - sigma * np.sqrt(T)
                d3 = (spot/self.strike) ** (-2 * b / sigma ** 2) * norm.cdf(-d1 + 2 * b / sigma * np.sqrt(T)) - np.exp(b * T) * norm.cdf(-d1)
                c = -spot * np.exp((b - r) * T) * norm.cdf(-d1) + self.strike * np.exp(-r * T) * norm.cdf(-d2)
                c += spot * np.exp(-r * T) * sigma ** 2 / 2 / b * (d3)
                return c
                      
    def f(self, x, y, a, b, pho):
        a1 = a / np.sqrt(2 * (1 - pho ** 2))
        b1 = b / np.sqrt(2 * (1 - pho ** 2))
        ff = np.exp(a1 * (2 * x - a1) + b1 * (2 * y - b1) + 2 * pho * (x - a1) * (y - b1))
        return ff
        
    def phi(self, a, b, pho):
        x1 = 0.24840615
        x2 = 0.39233107
        x3 = 0.21141819
        x4 = 0.033246660
        x5 = 0.00082485334
        y1 = 0.10024215
        y2 = 0.48281397
        y3 = 1.0609498
        y4 = 1.7797294
        y5 = 2.6697604
        x = [x1,x2,x3,x4,x5]
        y = [y1,y2,y3,y4,y5]
        phii = 0
        for i in range(5):
            for j in range(5):
                phii += x[i] * x[j] * self.f(y[i], y[j], a, b, pho)
        phii *= np.sqrt(1 - pho ** 2) / np.pi
        return phii
        
    def M(self, a, b, pho):
        if a <= 0 and b <= 0 and pho <= 0:
            return self.phi(a, b, pho)
        elif a <= 0 and b >= 0 and pho >= 0:
            return -self.phi(a, -b, -pho) + norm.cdf(a)
        elif a >= 0 and b <= 0 and pho >= 0: 
            return norm.cdf(b) - self.phi(-a, b, -pho)
        elif a >= 0 and b >= 0 and pho <= 0:
            return norm.cdf(a) + norm.cdf(b) - 1 + self.phi(-a ,-b ,pho)
        else:
            sign_a = np.sign(a)
            sign_b = np.sign(b)
            if sign_a == 0:
                sign_a =1
            if sign_b ==0:
                sign_b =1          
            p1 = ((pho * a - b) * sign_a) / np.sqrt(a ** 2 - 2 * pho * a * b + b ** 2)
            p2 = ((pho * b - a) * sign_b) / np.sqrt(a ** 2 - 2 * pho * a * b + b ** 2)
            delta = (1 - sign_a * sign_b) / 4
            return self.M(a, 0 , p1) + self.M(b, 0, p2) - delta
            
    def partial(self):
        spot = self.spot
        r = self.r
        sigma = self.sigma
        t2 = self.T
        b = self.b
        t1 = self.Partial_T
        if self.OptionType == 'PartialTimeFloatingStrike_call':
            M0 = self.Smin
        elif self.OptionType == 'PartialTimeFloatingStrike_put':
            M0 = self.Smax
        elif self.OptionType == 'PartialTimeFixedStrike_call':
            M0 = self.strike
        elif self.OptionType == 'PartialTimeFixedStrike_put':
            M0 = self.strike
            
        d1 = (np.log(spot / M0) + (b + sigma ** 2 / 2) * t2) / sigma / np.sqrt(t2)
        d2 = d1 - sigma * np.sqrt(t2)
        e1 = (b + sigma ** 2 / 2) * (t2- t1) / sigma / np.sqrt(t2 - t1)
        e2 = e1 - sigma * np.sqrt(t2 - t1)
        f1 = (np.log(spot / M0) + (b + sigma ** 2 / 2) * t1) / sigma / np.sqrt(t1)
        f2 = f1 - sigma * np.sqrt(t1)
        g1 = 0
        g2 = 0
        if self.OptionType == 'PartialTimeFloatingStrike_call':       
            c = spot * np.exp((b - r) * t2) * norm.cdf(d1 - g1) - self.Smin * np.exp(-r * t2) * norm.cdf(d2 - g1)
            c += spot * np.exp(-r * t2) * sigma ** 2 / 2 / b * ((spot / self.Smin) ** (-2 * b / sigma ** 2) * self.M(
                -f1 + 2 * b * np.sqrt(t1) / sigma, - d1 + 2 * b * np.sqrt(t2) / sigma - g1, np.sqrt(t1 / t2)) - np.exp(b * t2) * self.M(
                -d1 - g1, e1 + g2, -np.sqrt(1 - t1 / t2)))
            c += spot * np.exp((b - r) * t2) * self.M(-d1 + g1, e1 - g2, -np.sqrt(1- t1 / t2)) + self.Smin * np.exp(
                -r * t2) * self.M(-f2, d2 - g1, - np.sqrt(t1 / t2)) - np.exp(-b * (t2 - t1)) * (1 + sigma ** 2 / 2 / b) * spot * np.exp(
                (b - r) * t2) * norm.cdf(e2 - g2) * norm.cdf(-f1)
            return c
        elif self.OptionType == 'PartialTimeFloatingStrike_put':
            c = -spot * np.exp((b - r) * t2) * norm.cdf(-d1 + g1) + self.Smax * np.exp(-r * t2) * norm.cdf(-d2 + g1)
            c += spot * np.exp(-r * t2) * sigma ** 2 / 2 / b * (-(spot / self.Smax) ** (-2 * b / sigma ** 2) * self.M(
                f1 - 2 * b * np.sqrt(t1) / sigma, d1 - 2 * b * np.sqrt(t2) / sigma + g1, np.sqrt(t1 / t2)) + np.exp(b * t2) * self.M(
                d1 + g1, -e1 - g2, -np.sqrt(1 - t1 / t2)))
            c += -spot * np.exp((b - r) * t2) * self.M(d1 - g1, -e1 + g2, -np.sqrt(1- t1 / t2)) - self.Smax * np.exp(
                -r * t2) * self.M(f2, -d2 + g1, - np.sqrt(t1 / t2)) + np.exp(-b * (t2 - t1)) * (1 + sigma ** 2 / 2 / b) * spot * np.exp(
                (b - r) * t2) * norm.cdf(-e2 + g2) * norm.cdf(f1)
            return c
        elif self.OptionType == 'PartialTimeFixedStrike_call':
            c = spot * np.exp((b - r) * t2) * norm.cdf(d1) - self.strike * np.exp(-r * t2) * norm.cdf(d2) + spot * np.exp(-r * t2) * sigma ** 2 / 2 / b *(
                -(spot / self.strike) ** (- 2 * b / sigma ** 2) * self.M(d1 - 2 * b * np.sqrt(t2) / sigma, - f1 + 2 * b * np.sqrt(
                t1) / sigma, -np.sqrt(t1 / t2)) + np.exp(b * t2) * self.M(e1, d1, np.sqrt(1 - t1 / t2)))
            c += -spot * np.exp((b - r) * t2) * self.M(-e1, d1, -np.sqrt(1- t1 / t2)) - self.strike * np.exp(-r * t2) * self.M(
                f2, -d2, -np.sqrt(t1 / t2)) + np.exp(-b * (t2 - t1)) * (1 - sigma ** 2 / 2 / b) * spot * np.exp((b - r) * t2) * norm.cdf(
                f1) * norm.cdf(-e2)
            return c
        elif self.OptionType == 'PartialTimeFixedStrike_put':
            c = -spot * np.exp((b - r) * t2) * norm.cdf(-d1) + self.strike * np.exp(-r * t2) * norm.cdf(-d2) + spot * np.exp(-r * t2) * sigma ** 2 / 2 / b *(
                (spot / self.strike) ** (- 2 * b / sigma ** 2) * self.M(-d1 + 2 * b * np.sqrt(t2) / sigma, f1 - 2 * b * np.sqrt(
                t1) / sigma, -np.sqrt(t1 / t2)) - np.exp(b * t2) * self.M(-e1, -d1, np.sqrt(1 - t1 / t2)))
            c += spot * np.exp((b - r) * t2) * self.M(e1, -d1, -np.sqrt(1- t1 / t2)) + self.strike * np.exp(-r * t2) * self.M(
                -f2, d2, -np.sqrt(t1 / t2)) - np.exp(-b * (t2 - t1)) * (1 - sigma ** 2 / 2 / b) * spot * np.exp((b - r) * t2) * norm.cdf(
                -f1) * norm.cdf(e2)
            return c            
            
    def execute(self):
        if self.OptionType == 'FloatingStrike_call':
            return self.non_partial()
        elif self.OptionType == 'FloatingStrike_put':
            return self.non_partial()
        elif self.OptionType == 'FixedStrike_call':
            return self.non_partial()
        elif self.OptionType == 'FixedStrike_put':
            return self.non_partial()
        else:
            return self.partial()
        
        
            
            
            
if __name__ == '__main__':
    obj = LookbackOption(120,0.1,0.3,0.5,'FloatingStrike_call', Smin=100, dividend=0.06)
    print (obj.execute())
    obj2 = LookbackOption(100,0.1,0.3,0.5,'FixedStrike_call', Smax=100,strike = 95)
    print (obj2.execute())
    obj3 = LookbackOption(100,0.1,0.2,0.5,'FixedStrike_call', Smax=100,strike = 105)
    print (obj3.execute())
    obj4 = LookbackOption(100,0.1,0.2,1,'FixedStrike_put', Smin=100,strike = 105)      
    print (obj4.execute())
    obj5 = LookbackOption(100,0.1,0.2,0.5,'FixedStrike_put', Smin=100,strike = 95)      
    print (obj5.execute())
    print ('')
    obj6= LookbackOption(90,0.06,0.3,1,'PartialTimeFloatingStrike_call',partial_T=0.5, Smin=90)
    print (obj6.execute())
    obj7= LookbackOption(100,0.06,0.3,1,'PartialTimeFixedStrike_call',partial_T=0.75, strike=90)
    print (obj7.execute())
    obj8= LookbackOption(100,0.06,0.3,1,'PartialTimeFixedStrike_put',partial_T=0.75, strike=110)
    print (obj8.execute())
    obj9= LookbackOption(90,0.06,0.3,1,'PartialTimeFloatingStrike_put',partial_T=0.5, Smax=90)
    print (obj9.execute())
    print ('')
    obj10 = LookbackOption(85,0.1,0.3,0.75/365,'FixedStrike_call', Smax=85, strike=90)
    v0=obj10.execute()
    obj11 = LookbackOption(85,0.1,0.3,0.75/365,'FixedStrike_call', Smax=85, strike=0.00000000000001)
    v1=obj11.execute()
    print (v0+np.exp(-0.1*0.75/365)*90 - v1)