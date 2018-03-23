import csv
from mpmath import mp
mp.dps=30
mp.pretty=True

def read_data(filename):
    with open(filename, 'r') as csvfile:
        reader = csv.reader(csvfile)
        data = list(reader)    
    return data

def append_data(filename, rows):
    with open(filename, 'a') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerows(rows)

def periodic_saving(i,N,filename,rows): 
    if(i%N==0): 
       append_data(filename,rows)
       rows[:] = []    
  
def mpf(x): return mp.mpf(x);
def mpc(x,y): return mp.mpc(x,y);
def stringtonum(z): return mp.mpmathify(z.replace('*I','j'));
def floor(x): return mp.floor(x);
def ceil(x): return mp.ceil(x);
def log(n): return mp.log(n);
def exp(n): return mp.exp(n);
def sqrt(x): return mp.sqrt(x);
def power(a,x): return mp.power(a,x);
def gamma(z): return mp.gamma(z);
def cos(z): return mp.cos(z);
def sin(z): return mp.sin(z);
def conj(a): return a.conjugate();
def sum(n,N,summand): return mp.nsum(summand,[n,N]);
def intnum(u_lowlim,u_uplim,integrand): return mp.quad(integrand,[u_lowlim,u_uplim]);
Pi=mp.pi();
I = 1j;

def alpha1(s): return(1/(2*s) + 1/(s-1) + (1/2)*log(s/(2*Pi)));
def alpha1prime(s): return(-1/(2*s**2) - 1/(s-1)**2 + 1/(2*s));
def H01(s): return((1/2)*s*(s-1)*Pi**(-s/2)*sqrt(2*Pi)*exp((s/2-1/2)*log(s/2)-s/2));
def S(N,sigmavar,t): return(sum(1,N,lambda n: n**(sigmavar + (t/4.0)*log(N**2/n))));
def C0(p): return((exp(Pi*I*(p**2/2 + 3/8)) - I*sqrt(2)*cos(Pi*p/2))/(2*cos(Pi*p)));
def B0_eff(x,y=0.4,t=0.4): return((1/8)*exp((t/4)*alpha1((1+y-I*x)/2)**2)*H01((1+y-I*x)/2));

def thetafunc_new(x,y=0.4): return(Pi/8 - (1/4)*atan((9+y)/x));
def Jbound_t_th(t,th,X,a,beta): return(exp(-t*th**2)*intnum(0,X,lambda u: sqrt(th**2+u**2)*exp(t*u**2 - beta*exp(4*u)*cos(4*th) + a*u)));

def ddx_Ht_bound(x,y=0.4,t=0.4,n0=50,X=6):
   th = thetafunc_new(x,y);
   n0 = max(1,floor(sqrt(x/Pi)));
   main_est1 = 2*(Pi**2)*sum(1,n0,lambda n: (n**4)*(Jbound_t_th(t,th,X,(9-y),Pi*n**2) + Jbound_t_th(t,th,X,(9+y),Pi*n**2)));
   main_est2 = 3*Pi*sum(1,n0,lambda n: (n**2)*(Jbound_t_th(t,th,X,(5-y),Pi*n**2) + Jbound_t_th(t,th,X,(5+y),Pi*n**2)));
   main_est = (1/2)*(main_est1 + main_est2);
   return(main_est);   

def ddx_Ht_bound_optim(x,y=0.4,t=0.4,n0=50,X=6):
   th = thetafunc_new(x,y);
   n0 = max(1,floor(sqrt(x/Pi)));
   a1 = 9-y; a2 = 9+y; a3 = 5-y; a4 = 5+y; 
   pi_fac1 = 2*(Pi^2); pi_fac2 = 3*Pi;
   main_est = (1/2)*exp(-t*th**2)*intnum(0,X,lambda u: sum(1,n0,lambda n: sqrt(th**2+u**2)*exp(t*u**2 - (Pi*n**2)*exp(4*u)*cos(4*th))*(pi_fac1*(n**4)*(exp(a1*u)+exp(a2*u))+pi_fac2*(n**2)*(exp(a3*u)+exp(a4*u)))))
   return(main_est);   

def abceff_x(x,y=0.4,t=0.4):
    T = x/2;       
    Tdash = T + Pi*t/8;
    a=sqrt(Tdash/(2*Pi));
    N=floor(a);
    p = 1 - 2*(a-N);
    U = exp(-I*((Tdash/2)*log(Tdash/(2*Pi)) - Tdash/2 - Pi/8));
    sig = (1-y)/2;
    s = sig + I*T;
    sdash = sig + I*Tdash;
    alph1 = alpha1(s);
    alph2 = alpha1(1-s);
    A0 = exp((t/4)*alph1**2)*H01(s);
    B0 = exp((t/4)*alph2**2)*H01(1-s);
    A_sum = sum(1,N,lambda n: n**((t/4.0)*log(n) - (t/2.0)*alph1 - s));
    B_sum = sum(1,N,lambda n: n**((t/4.0)*log(n) - (t/2.0)*alph2 - (1-s)));
    A = A0 * A_sum;
    B = B0 * B_sum;
    termC1 = Pi**(-sdash/2)*gamma(sdash/2)*(a**(-sig))*C0(p)*U;
    termC2 = Pi**(-(1-sdash)/2)*gamma((1-sdash)/2)*(a**(-(1-sig)))*conj(C0(p))*conj(U);
    C = exp(t*Pi**2/64)*(sdash*(sdash-1)/2)*((-1)**N)*(termC1 + termC2);
    return((A+B-C)/8);

def newton_quot_abc(x,y=0.4,t=0.4,h=0.000001): return((abceff_x(x+h,y,t)-abceff_x(x,y,t))/h);

#fixed mesh data obtained from a different tool eg. Arb
in_filename = "Ht_data_1000_1600.txt"
out_filename = "Ht_data_1000_1600_verified_v4.txt"
data = read_data(in_filename)
data = [[stringtonum(x) for x in row] for row in data]


/*in_filename1 = "values2x.txt"
in_filename2 = "hvalues2.txt"
out_filename = "Ht_data_10_1000_verified_v4.txt"

xvals = read_data(in_filename1)[0]
xvals[0]=xvals[0][1:]; xvals[-1] = xvals[-1][:-1];
hvals = read_data(in_filename2)[0]
hvals[0]=hvals[0][1:]; hvals[-1] = hvals[-1][:-1];
data = [[mpf(xvals[i]),mpc(hvals[2*i],hvals[2*i+1])] for i in range(len(xvals))] 
*/

stepsize = 0.005
xc_step = 2.0
data_new=[]
i=0
x_next = 9.999
xc=100.0; thc = thetafunc_new(xc); D = ddx_Ht_bound_optim(xc);
for hrow in data:
    x=hrow[0];
    if x < x_next and abs(x-x_next)>0.00001: continue
    H = hrow[1];
    ABC = abceff_x(x); NQ_ABC = newton_quot_abc(x); B0 = B0_eff(x);
    xc_candidate = xc_step*floor(x/xc_step)
    if(xc < xc_candidate): 
       xc = xc_candidate; print(xc); thc = thetafunc_new(xc); D = ddx_Ht_bound_optim(xc);
    bnd_dH = D*exp(-thc*x);
    allowed_step = abs(H/bnd_dH); 
    if (stepsize <= allowed_step): interval_cleared = 1
    else: interval_cleared = 0
    new_hrow = [x,H,interval_cleared,abs(H/B0),abs(ABC/B0),thc,D,abs(bnd_dH/B0),abs(NQ_ABC/B0),allowed_step]
    #print(i, new_hrow)
    data_new.append(new_hrow)
    periodic_saving(i,100,out_filename,data_new)
    i+=1;
    x_next = x + max(stepsize,int(allowed_step/stepsize)/(1/stepsize))

append_data(out_filename,data_new)    
