import sympy as sym
import numpy as np

z=sym.Symbol('x')
z=sym.Symbol('z')

lambda0=756e-9 # laser central wavelength in meters
c=299792458
w0=2*np.pi*c/lambda0

f=(z-1)**2

# print(sym.diff(f,z).subs(z,0))

x=sym.Symbol('x')

n_air=1+0.05792105/(238.0185-(1e-12)*(x/(2*np.pi*c))**2)+0.00167917/(57.362-(1e-12)*(x/(2*np.pi*c))**2)
n_bk7=(1+1.03961212*(1e12)*(2*np.pi*c/x)**2/((1e12)*(2*np.pi*c/x)**2-0.00600069867)+0.231792344*(1e12)*(2*np.pi*c/x)**2/((1e12)*(2*np.pi*c/x)**2-0.0200179144)+1.01046945*(1e12)*(2*np.pi*c/x)**2/((1e12)*(2*np.pi*c/x)**2-103.560653))**0.5

n0_air=n_air.subs(x,w0)
np0_air=sym.diff(n_air,x).subs(x,w0)
npp0_air=sym.diff(sym.diff(n_air,x),x).subs(x,w0)
nppp0_air=sym.diff(sym.diff(sym.diff(n_air,x),x)).subs(x,w0)

k0_air=n0_air*w0/c
kp0_air=(n0_air+w0*np0_air)/c
kpp0_air=(2*np0_air+w0*npp0_air)/c
kppp0_air=(3*npp0_air+w0*nppp0_air)/c

print("n0_air=",n0_air)
print("np0_air=",np0_air)
print("npp0_air=",npp0_air)
print("npp0_air=",nppp0_air)

print("k0_air=",k0_air)
print("kp0_air=",kp0_air)
print("kpp0_air=",kpp0_air)
print("kppp0_air=",kppp0_air)

n0_bk7=n_bk7.subs(x,w0)
np0_bk7=sym.diff(n_bk7,x).subs(x,w0)
npp0_bk7=sym.diff(sym.diff(n_bk7,x),x).subs(x,w0)
nppp0_bk7=sym.diff(sym.diff(sym.diff(n_bk7,x),x)).subs(x,w0)

k0_bk7=n0_bk7*w0/c
kp0_bk7=(n0_bk7+w0*np0_bk7)/c
kpp0_bk7=(2*np0_bk7+w0*npp0_bk7)/c
kppp0_bk7=(3*npp0_bk7+w0*nppp0_bk7)/c

print("n0_bk7=",n0_bk7)
print("np0_bk7=",np0_bk7)
print("npp0_bk7=",npp0_bk7)
print("nppp0_bk7=",nppp0_bk7)

print("k0_bk7=",k0_bk7)
print("kp0_bk7=",kp0_bk7)
print("kpp0_bk7=",kpp0_bk7)
print("kppp0_bk7=",kppp0_bk7)