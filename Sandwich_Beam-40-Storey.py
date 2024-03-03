import numpy as np
import pandas as pd
from scipy import interpolate

num_f = 40             # Kat sayısı
num_h = 3.29           # (m) Kat yüksekliği
H = num_f * num_h      # (m) Bina yüksekliği
bay1 = 7               # (m) 2. açıklık
bay2 = 7               # (m) 3. açıklık
E = 200 * 10 ** 6      # (kN/m2) Elastisite modülü
G = 77 * 10 ** 6       # (kN/m2) Kayma modülü
mass = 250             # (t/kat)
Mt = 2 * mass * num_f  # (t) Toplam kütle (2 adet perde)
m = mass / num_h       # (t/m) Kat yüksekliğine yayılan kütle

"""### Kolon Özellikleri ###"""
# HD400x347
Ic = 0.1249 * 10**(-2)   # (m4) Atalet momenti
Ac = 442 * 10**(-4)      # (m2) Alan
d = 0.407                # (m)  Profil yüksekliği
Afl = 17655 * 10**(-6)   # (m2) Flanş alanı
Aweb = 11150 * 10**(-6)  # (m2) Gövde alanı
tweb = 27.2 * 10**(-3)   # (m)  Gövde kalınlığı

"""### Kiriş Özellikleri ###"""
# HEA400
Ib = 0.4507 * 10**(-3)  # (m4) Atalet momenti
Lb = 7                  # (m)  Kiriş uzunluğu

"""### Duvar Özellikleri ###"""
plw = 6      # (m) Levha genişliği
ptk = 0.006  # (m) Levha kalınlığı

"""### Deprem Parametreleri ###"""
Sa1 = 0.123654  # m/s2 (Spektral ivme)
Sa2 = 0.558422  # m/s2 (Spektral ivme)
Sa3 = 1.225113  # m/s2 (Spektral ivme)
Sd1 = 0.174414  # m (Spektral yer değiştirme)
Sd2 = 0.059739  # m (Spektral yer değiştirme)
Sd3 = 0.027230  # m (Spektral yer değiştirme)

correction = [[0.493], [0.653], [0.770], [0.812], [0.842], [0.863], [0.879],
              [0.892], [0.902], [0.911], [0.918], [0.924], [0.929], [0.934],
              [0.938], [0.941], [0.947], [0.952], [0.961], [0.967], [0.980]]

index_rf = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 25, 30, 50]
columns_rf = ["rf"]

df_rf = pd.DataFrame(correction, index=index_rf, columns=columns_rf)

x = False
for i in index_rf:
    if i == num_f:
        df_rf = df_rf.loc[num_f]
        x = True

if x != True:
    df_rf.loc[num_f] = np.nan  # Kat sayısına göre bilinmeyen satırı eklenmesi
    df_rf = df_rf.reindex(sorted(df_rf.index), axis=0)  # Listenin tekrar sıralanması
    df_rf = df_rf.interpolate(method="polynomial", order=1).round(4)  # 1. derece polinom interpolasyonu
    df_rf = df_rf.loc[num_f]

rf = df_rf.loc["rf"]  # Kat seviyelerinde yığılmış kütlelerin katkısını dikkate alan bir faktör

Q1 = Afl * (plw * 0.5 + d)
Q2 = Q1 + Aweb * 0.5 * (plw + d)
Q3 = Ac * 0.5 * (plw + d)
Q4 = Q3 + plw**2 * ptk / 8

beta1 = (Q1**2 + Q2**2) * d / tweb
beta2 = (Q3**2 + Q4**2) * plw / (2 * ptk)
beta = beta1 + beta2

Iw = (ptk * plw**3 / 12) + 2 * Ac * ((d + plw) / 2)**2 + 2 * Ic
# Iw: Yanal yüke dayanıklı sistemin bir parçası olan her duvar için, değiştirilmiş ikinci alan momenti
KGAw = (Iw**2 / beta) * G

fb = (rf * 0.5595 / H**2) * (E * Iw / m)**0.5
fs = (rf / (4 * H)) * (KGAw / m)**0.5
Tw = ((1 / fb**2) + (1 / fs**2))**0.5

Imw = (m * H**4) / (0.313 * rf**2 * Tw**2 * E)

I = Imw + 4 * Ic  # 4 kolon adedi

Ks1 = (12 * E) / (num_h * (1 / (4 * Ic / num_h) + 1 / (2 * Ib / Lb)))

r = plw / (2 * Lb)
n = 6 * Ic * Lb / (Ib * num_h)
s = (n - 3 * r - 1) / (n + 2)
Ks2 = (2 * (6 * E * Ib) / (Lb * num_h)) * ((1 + r) * (1 + 2 * r + s))

Ig = Ac*(2*(plw/2)**2 + 2*((plw/2)+bay1)**2 + 2*((plw/2)+bay1+bay2)**2)

ffs2 = (rf**2 * (Ks1 + Ks2)) / (m * (4 * H)**2)
ffb2 = (0.313 * rf**2 * E * Ig) / (m * H**4)
Ksi = ffb2 / (ffb2 + ffs2)
Ks = Ksi * (Ks1 + Ks2)

alpha = round(((Ks / (E * I))**0.5), 9)

k = round((alpha * H), 3)

cc = round(((num_f / (num_f + 2.06))**0.5), 7)  # Düzeltme katsayısı

data = [[1.788, 0.285, 0.102, 1.57, -0.87, 0.51, 0.61, 0.19, 0.07],
        [1.529, 0.276, 0.101, 1.55, -0.85, 0.50, 0.62, 0.18, 0.06],
        [1.160, 0.254, 0.098, 1.52, -0.82, 0.50, 0.65, 0.16, 0.06],
        [0.908, 0.227, 0.094, 1.47, -0.77, 0.48, 0.67, 0.14, 0.06],
        [0.744, 0.200, 0.089, 1.43, -0.70, 0.47, 0.69, 0.12, 0.06],
        [0.631, 0.178, 0.083, 1.39, -0.66, 0.47, 0.71, 0.11, 0.05],
        [0.547, 0.160, 0.078, 1.37, -0.62, 0.44, 0.72, 0.11, 0.05],
        [0.483, 0.144, 0.073, 1.35, -0.59, 0.42, 0.73, 0.10, 0.05],
        [0.432, 0.132, 0.068, 1.33, -0.58, 0.41, 0.73, 0.10, 0.04],
        [0.391, 0.121, 0.064, 1.32, -0.56, 0.38, 0.74, 0.10, 0.04],
        [0.357, 0.111, 0.060, 1.31, -0.52, 0.37, 0.75, 0.09, 0.04],
        [0.328, 0.103, 0.056, 1.31, -0.51, 0.37, 0.75, 0.09, 0.04],
        [0.304, 0.096, 0.053, 1.30, -0.50, 0.36, 0.76, 0.09, 0.04],
        [0.282, 0.090, 0.050, 1.30, -0.49, 0.35, 0.76, 0.09, 0.04],
        [0.264, 0.085, 0.047, 1.30, -0.48, 0.35, 0.76, 0.10, 0.04],
        [0.248, 0.080, 0.045, 1.29, -0.47, 0.34, 0.77, 0.09, 0.04],
        [0.234, 0.075, 0.043, 1.29, -0.47, 0.33, 0.77, 0.09, 0.04],
        [0.221, 0.072, 0.041, 1.29, -0.46, 0.32, 0.77, 0.10, 0.04],
        [0.209, 0.068, 0.039, 1.29, -0.45, 0.32, 0.77, 0.09, 0.04],
        [0.199, 0.065, 0.037, 1.29, -0.45, 0.31, 0.77, 0.09, 0.03],
        [0.190, 0.062, 0.036, 1.29, -0.44, 0.30, 0.78, 0.09, 0.04],
        [0.129, 0.042, 0.025, None, None, None, None, None, None]]

index_k = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
           11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 30.0]

columns_k = ["Z1", "Z2", "Z3", "mu1", "mu2", "mu3", "eko1", "eko2", "eko3"]

df = pd.DataFrame(data, index=index_k, columns=columns_k)

y = False
for i in index_k:
    if i == k:
        df = df.loc[k]
        y = True

if y != True:
    df.loc[k] = np.nan  # k değerine göre bilinmeyen satırı eklenmesi
    df = df.reindex(sorted(df.index), axis=0)  # Listenin k değeri ile tekrar sıralanması
    df = df.interpolate(method="polynomial", order=1).round(4)  # 1. derece polinom interpolasyonu
    df = df.loc[k]

Z1 = df.loc["Z1"]
Z2 = df.loc["Z2"]
Z3 = df.loc["Z3"]

mu1 = df.loc["mu1"]
mu2 = df.loc["mu2"]
mu3 = df.loc["mu3"]

eko1 = df.loc["eko1"]
eko2 = df.loc["eko2"]
eko3 = df.loc["eko3"]

T1 = (Z1 * H**2 / cc) * ((m / (E * I)))**0.5
T2 = (Z2 * H**2 / cc) * ((m / (E * I)))**0.5
T3 = (Z3 * H**2 / cc) * ((m / (E * I)))**0.5
print("T1 =", T1.round(3), "s")
print("T2 =", T2.round(3), "s")
print("T3 =", T3.round(3), "s")
print(" ")

Vt1 = eko1 * Mt * Sa1
Vt2 = eko2 * Mt * Sa2
Vt3 = eko3 * Mt * Sa3
print("Vt1 =", Vt1.round(3), "kN")
print("Vt2 =", Vt2.round(3), "kN")
print("Vt3 =", Vt3.round(3), "kN")
print(" ")

Vt = (Vt1**2 + Vt2**2 + Vt3**2)**0.5
print("Vt =", Vt.round(3), "kN")
print(" ")

dH1 = mu1 * Sd1
dH2 = mu2 * Sd2
dH3 = mu3 * Sd3
print("dH1 =", dH1.round(4), "m")
print("dH2 =", dH2.round(4), "m")
print("dH3 =", dH3.round(4), "m")
print(" ")

dH = (dH1**2 + dH2**2 + dH3**2)**0.5
print("dH =", dH.round(4), "m")
