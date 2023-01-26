from math import *
import numpy as np

class NS_atmosphere:
    def __init__ (self, model_name, g14, Tb, B):

        if model_name == 'Potekhin_2015':
            self.name = 1
            self.g14 = g14
            self.T9  = Tb / 1e9
            self.B12 = B / 1e12
        elif model_name == 'Potekhin_2003_iron':
            self.name = 2
            self.g14 = g14
            self.T9 = Tb / 1e9
            self.B12 = B / 1e12

        elif model_name == 'Potekhin_Yakovlev_2001_iron':
            self.name = 3
            self.g14 = g14
            self.T9 = Tb / 1e9
            self.B12 = B / 1e12

        elif model_name == 'Potekhin_2003_accr':
            self.name = 4
            self.g14 = g14
            self.T9 = Tb / 1e9
            self.B12 = B / 1e12


        else:
            print ('Model is unknown')
            self.name = 0
            sys.exit(1)

    ## Here we introduce functions to describe individual models

    def Potekhin_2015_T_polar_0 (self):

        a   = 0.337 / (1.0 + 0.02 * sqrt(self.B12))
        T1  = 1.13 * pow(self.B12, 0.119) * pow(self.T9, a)
        T0  = pow(15.7 * pow(self.T9, 3.0/2.0) + 1.36*self.T9, 0.3796)
        Tp  = pow(self.g14 * (T1**4 + (1.0 + 0.15 * sqrt(self.B12))*T0**4), 1.0 / 4.0) * 1e6

        self.Tp = Tp
        return Tp

    def Potekhin_2015_T_polar_max (self):
        res = (5.2 * pow(self.g14, 0.65) + 0.093*sqrt(self.g14*self.B12)) * 1e6
        return res

    def Potekhin_2015_T_p (self):
        Tp_0 = self.Potekhin_2015_T_polar_0 ()
        Tp_max = self.Potekhin_2015_T_polar_max ()
        res = Tp_0 * pow(1.0 + pow(Tp_0 / Tp_max,4.0), -1.0 / 4.0)
        return res

    def Potekhin_2015_T_eq (self):
        Tp_val = self.Potekhin_2015_T_p ()
        res = 1.0 + pow(1230*self.T9, 3.35) * self.B12 * sqrt(1.0 + 2*self.B12*self.B12) / pow(self.B12 + 450*self.T9 + 119*self.B12*self.T9,4.0)
        res = res + 0.0066*pow(self.B12, 5.0 / 2.0) / (sqrt(self.T9) + 0.00258*pow(self.B12, 5.0/2.0))
        res = Tp_val / res
        return res





    def Ts (self, theta):

        if self.name == 1:
            a2 = 10 * self.B12 / (sqrt(self.T9) + 0.1 * self.B12 * pow(self.T9, -1.0 / 4.0))
            a1 = a2 * sqrt(self.T9) / 3.0
            gf = (1.0 + a1 + a2) * np.cos(theta) * np.cos(theta) / (1.0 + a1 * np.abs(np.cos(theta)) + a2 * np.cos(theta) * np.cos(theta))
            Tp_val   = self.Potekhin_2015_T_p ()
            T_eq_val = self.Potekhin_2015_T_eq ()
            res = gf * (Tp_val - T_eq_val) + T_eq_val
            return res

        if self.name == 2:

             ## Numerical constants from Table 2 in Appendix
             a1 = 1.76e-4
             a2 = 0.038
             a3 = 1.5
             a4 = 0.0132
             a5 = 0.620
             a6 = 0.318
             a7 = 2.3e-9
             a8 = 3
             a9 = 0.160
             a10 = 21
             a11 = 4.7e5

             b1 = 159
             b2 = 270
             b3 = 172
             b4 = 110
             b5 = 0.363
             b6 = 0.181
             b7 = 0.50
             b8 = 0.619

             ksi = self.T9 - 1e-3 * pow(self.g14, 1./4.) * sqrt(7.0*self.T9)

             T0 = self.g14 * (pow(7*ksi, 2.25) + pow(0.33*ksi, 1.25))
             T0 = np.power(T0, 1.0 / 4.0)

             #print (T0)


             beta = 1.0 / (1.0 + b5 * pow(self.T9, b6))
             ksi_par = (1.0 + (a1 + a2 * pow(self.T9, a3))/(pow(self.T9, 2) + a4*pow(self.T9, a5)) * pow(self.B12, a6) / pow(1.0 + a7*self.B12 / pow(self.T9, a8), a9))
             ksi_par = ksi_par / (1.0 + 1.0 / (3.7 + (a10 + a11*pow(self.B12, -3.0 / 2.0))*self.T9*self.T9) )

             ksi_ort = sqrt(1.0 + b1*self.B12 / (1.0 + b2*pow(self.T9, b7))) / pow(1.0 + b3*self.B12 / (1.0 + b4*pow(self.T9, b8)),beta)

             alpha = 4.0 + sqrt(ksi_ort / ksi_par)

             res = np.power(pow(ksi_par, alpha)*np.power(np.cos(theta),2.0) + pow(ksi_ort, alpha) * np.power(np.sin(theta), 2.0)  , 1.0 / alpha)
             res = res * T0 * 1e6

             return res


        if self.name == 3:

             ksi = self.T9 - 0.001 * pow(self.g14, 1.0/4.0) * sqrt(7.0 * self.T9)

             Ts_0 = 1e6 * pow(self.g14, 1./4.) * pow(pow(7.0 * ksi, 2.25) + pow(ksi/3.0, 1.25), 1.0 / 4.0)

             ksi_par = 1.0 + 0.0492 * pow(self.B12, 0.292) / pow(self.T9, 0.240)

             ksi_ort = sqrt(1.0 + 0.1076 * self.B12*pow(0.03 + self.T9, -0.559)) / pow(1.0 + 0.819*self.B12 / (0.03 + self.T9), 0.6463)

             ksi_t = np.power(pow(ksi_par, 9.0 / 2.0) * np.power(np.cos(theta), 2.0) + pow(ksi_ort, 9.0 / 2.0) * np.power(np.sin(theta), 2.0), 2.0 / 9.0)

             res = Ts_0 * ksi_t

             return res

        if self.name == 4:

             ## Numerical constants from Table 2 in Appendix
             a1 = 4.5e-3
             a2 = 0.055
             a3 = 2.0
             a4 = 0.0595
             a5 = 0.328
             a6 = 0.237
             a7 = 6.8e-7
             a8 = 2
             a9 = 0.113
             a10 = 163
             a11 = 3.4e5

             b1 = 172
             b2 = 155
             b3 = 383
             b4 = 94
             b5 = 0.383
             b6 = 0.367
             b7 = 2.82
             b8 = 1.69

             ksi = self.T9 - 1e-3 * pow(self.g14, 1./4.) * sqrt(7.0*self.T9)

             TFe6 = self.g14 * (pow(7*ksi, 2.25) + pow(0.33*ksi, 1.25))
             TFe6 = np.power(TFe6, 1.0 / 4.0)


             T0 = (self.g14 * pow(18.1 * self.T9, 2.42) * (0.447 + 0.075 * log10(self.T9*1e9) / (1 + pow(6.2*self.T9, 4.0)) ) + 3.2*pow(self.T9, 1.67)*pow(TFe6, 4)) / (1.0 + 3.2*pow(self.T9, 1.67))  
             T0 = np.power(T0, 1.0 / 4.0)

             #print (T0)


             beta = 1.0 / (1.0 + b5 * pow(self.T9, b6))
             ksi_par = (1.0 + (a1 + a2 * pow(self.T9, a3))/(pow(self.T9, 2) + a4*pow(self.T9, a5)) * pow(self.B12, a6) / pow(1.0 + a7*self.B12 / pow(self.T9, a8), a9))
             ksi_par = ksi_par / (1.0 + 1.0 / (3.7 + (a10 + a11*pow(self.B12, -3.0 / 2.0))*self.T9*self.T9) )

             ksi_ort = sqrt(1.0 + b1*self.B12 / (1.0 + b2*pow(self.T9, b7))) / pow(1.0 + b3*self.B12 / (1.0 + b4*pow(self.T9, b8)),beta)

             alpha = 4.0 + sqrt(ksi_ort / ksi_par)

             res = np.power(pow(ksi_par, alpha)*np.power(np.cos(theta),2.0) + pow(ksi_ort, alpha) * np.power(np.sin(theta), 2.0)  , 1.0 / alpha)
             res = res * T0 * 1e6

             return res




    def describe (self):
        if self.name == 1:
            print ("Surface temperatures for magnetised envelope of NS")
            print ("Following fit from Potekhin, Pons & Page (2015) for iron envelope")
            print ("Space Science Reviews, Volume 191, Issue 1-4, pp. 239-291 2015SSRv..191..239P")
        elif self.name == 2:
            print ("Surface temperatures for magnetised envelope of NS")
            print ("Following fit from Potekhin, Yakovlev, Chabrier & Gnedin (2003) for iron envelope")
            print ("The Astrophysical Journal, Volume 594, Issue 1, pp. 404-418  2003ApJ...594..404P")
        elif self.name == 3:
            print ("Surface temperatures for magnetised envelope of NS")
            print ("Following fit from Potekhin & Yakovlev (2001) for iron envelope")
            print ("Astronomy and Astrophysics, v.374, p.213-226 (2001) 2001A&A...374..213P")
        elif self.name == 4:
            print ("Surface temperatures for magnetised envelope of NS")
            print ("Following fit from Potekhin, Yakovlev, Chabrier & Gnedin (2003) for fully accreted envelope")
            print ("The Astrophysical Journal, Volume 594, Issue 1, pp. 404-418  2003ApJ...594..404P")




