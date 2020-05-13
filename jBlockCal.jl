module jBlockCal
#Calibrate a sensor block and write the block parameters and block error maps to files.

#include("txtLeasqr.jl")  #include as text file  #pasted in below 9/18/2018
using jConstants
using jCalInit     #added 2019/6/7 to set global tdustglass
using PyPlot       #added 3/24/2019
using jQuiverPlot  #added 3/24/2019
using Statistics
using LinearAlgebra
using jGeometry
using Printf
export P1blockcal, P1blockcal_test, leasqr_test, P1c2plane_e, P1modelS1e_test, P1blockcalS1S2S3, xP1c2plane_e_orig

#*******************************************************************

function dfdp(x,f,p,dp,func)
# numerical partial derivatives [Jacobian] df/dp for use with leasqr
# --------INPUT VARIABLES---------
# x=vec or matrix of indep var(used as arg to func) x=[x0 x1 ...]
# f=func[x,p] vector initialsed by user before each call to dfdp
# p= vec of current parameter values
# dp= fractional increment of p for numerical derivatives
#      dp[j]>0 central differences calculated
#      dp[j]<0 one sided differences calculated
#      dp[j]=0 sets corresponding partials to zero; i.e. holds p[j] fixed
# func=string naming the function (.m) file
#      e.g. to calc Jacobian for function expsum prt=dfdp[x,f,p,dp,"expsum"]
#----------OUTPUT VARIABLES-------
# prt= Jacobian Matrix prt[i,j]=df[i]/dp[j]
#---------------------------------

m=size(x,1)
if m==1
    m=size(x,2)
end  ## PAK: in case() #cols > #rows
n=length(p)      #dimensions
ps=p; prt=zeros(m,n); del=zeros(n,1);       # initialise Jacobian to Zero
for j=1:n
      del[j]=dp[j] .*p[j];    #cal delx=fract[dp]*param value[p]
      if p[j]==0
           del[j]=dp[j];     #if param=0 delx=fraction
      end
      p[j]=ps[j] + del[j]
      if del[j]!=0
           #f1=feval(func,x,p)
           f1=func(x,p)  #for Julia
           if dp[j] < 0 prt[:,j]=(f1-f)./del[j]
           else
                p[j]=ps[j]- del[j]
                #prt[:,j]=(f1-feval(func,x,p))./(2 .*del[j])
                prt[:,j]=(f1-func(x,p))./(2 .*del[j])  #for Julia
           end
      end
      p[j]=ps[j];     #restore p[j]
end
    #println("prt", size(prt))  #for debugging
return prt  #for Julia
    #return prt'  #for test 9/9/2018
end


#******************************************************************

#using PyPlot

# Copyright [C] 1992-1994 Richard Shrager, Arthur Jutan, Ray Muzic
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version().
#
# This program is distributed in the hope that it will be useful
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

function leasqr(x,y,pin,F,stol,niter,wt,dp,dFdp,options)
#function() [f,p,kvg,iter,corp,covp,covr,stdresid,Z,r2]=
#                   leasqr[x,y,pin,F,stol,niter,wt,dp,dFdp,options}]
#
# Levenberg-Marquardt nonlinear regression of f[x,p] to y[x].
#
# Version 3.beta()
# Optional parameters are in braces }.
# x = column vector or matrix of independent variables, 1 row per
#   observation: x = [x0 x1....xm].
# y = column vector of observed values, same number of rows as x.
# wt = column vector [dim=length(x)] of statistical weights.  These
#   should be set to be proportional to [sqrt of var(y)]^-1; (That is()
#   the covariance matrix of the data is assumed to be proportional to
#   diagonal with diagonal equal to [wt.^2]^-1.  The constant of
#   proportionality will be estimated.); default = ones(length(y),1).
# pin = column vec of initial parameters to be adjusted by leasqr.
# dp = fractional increment of p for numerical partial derivatives
#   default = .001*ones(size(pin))
#   dp[j] > 0 means central differences on j-th parameter p[j].
#   dp[j] < 0 means one-sided differences on j-th parameter p[j].
#   dp[j] = 0 holds p[j] fixed i.e. leasqr wont change initial guess: pin[j]
# F = name of function in quotes; the function shall be of the form y=f[x,p]
#   with y, x, p of the form y, x, pin as described above.
# dFdp = name of partial derivative function in quotes; default is() "dfdp", a
#   slow but general partial derivatives function(); the function shall be
#   of the form prt=dfdp[x,f,p,dp,F] (see dfdp.m).
# stol = scalar tolerance on fractional improvement in scalar sum of
#   squares = sum((wt.*(y-f))^2); default stol = .0001
# niter = scalar maximum number of iterations; default = 20
# options = matrix of n rows [same number of rows as pin] containing
#   column 1: desired fractional precision in parameter estimates.
#     Iterations are terminated if change in parameter vector [chg] on two
#     consecutive iterations is less than their corresponding elements
#     in options[:,1].  [ie. all(abs(chg*current parm est) < options[:,1])
#      on two consecutive iterations.], default = zeros().
#   column 2: maximum fractional step change in parameter vector.
#     Fractional change in elements of parameter vector is constrained to be
#     at most options[:,2] between sucessive iterations.
#     [ie. abs(chg[i])=abs(min([chg[i] options[i,2]*current param estimate])).]
#     default = Inf()*ones().
#
#          OUTPUT VARIABLES
# f = column vector of values computed: f = F[x,p].
# p = column vector trial or final parameters. i.e, the solution.
# kvg = scalar: = 1 if convergence, = 0 otherwise().
# iter = scalar number of iterations used.
# corp = correlation matrix for parameters.
# covp = covariance matrix of the parameters.
# covr = diagm(covariance matrix of the residuals).
# stdresid = standardized residuals.
# Z = matrix that defines confidence region [see comments in the source].
# r2 = coefficient of multiple determination.
#
# All Zero guesses not acceptable

# A modified version of Levenberg-Marquardt
# Non-Linear Regression program previously submitted by R.Schrager.
# This version corrects an error in that version and also provides
# an easier to use version with automatic numerical calculation of
# the Jacobian Matrix. In addition, this version calculates statistics
# such as correlation, etc....
#
# Version 3 Notes
# Errors in the original version submitted by Shrager [now called version 1]
# and the improved version of Jutan [now called version 2] have been corrected.
# Additional features, statistical tests, and documentation have also been
# included along with an example of usage.  BEWARE: Some the the input and
# output arguments were changed from the previous version().
#
#     Ray Muzic     <rfm2@ds2.uh.cwru.edu>
#     Arthur Jutan  <jutan@charon.engga.uwo.ca>

# Richard I. Shrager [301]-496-1122
# Modified by A.Jutan [519]-679-2111
# Modified by Ray Muzic 14-Jul-1992
#       1) add maxstep feature for limiting changes in parameter estimates
#          at each step.
#       2) remove forced columnization of x [x=x[:]] at beginning. x could be
#          a matrix with the ith row of containing values of the
#          independent variables at the ith observation.
#       3) add verbose option
#       4) add optional return arguments covp, stdresid, chi2
#       5) revise estimates of corp, stdev
# Modified by Ray Muzic 11-Oct-1992
#	1) revise estimate of Vy.  remove chi2, add Z as return values
# Modified by Ray Muzic 7-Jan-1994
#       1) Replace ones(x) with a construct that is compatible with versions
#          newer and older than v 4.1.
#       2) Added global declaration of verbose [needed for newer than v4.x]
#       3) Replace return value var(), the variance of the residuals with covr
#          the covariance matrix of the residuals.
#       4) Introduce options as 10th input argument.  Include
#          convergence criteria and maxstep in it.
#       5) Correct calculation of xtx which affects coveraince estimate.
#       6) Eliminate stdev (estimate of standard deviation of parameter
#          estimates) from the return values.  The covp is a much more()
#          meaningful expression of precision because it specifies a confidence
#          region in contrast to a confidence interval..  If needed, however
#          stdev may be calculated as stdev=sqrt(diagm(covp)).
#       7) Change the order of the return values to a more logical order.
#       8) Change to more efficent algorithm of Bard for selecting epsL.
#       9) Tighten up memory usage by making use of sparse matrices (if
#          MATLAB version() >= 4.0) in computation of covp, corp, stdresid.
# Modified by Francesco Potort�
#       for use in Octave
#
# References:
# Bard, Nonlinear Parameter Estimation, Academic Press, 1974.
# Draper and Smith, Applied Regression Analysis, John Wiley and Sons, 1981.
#
#set default args

# argument processing
#
#if [sscanf(version(),"#f") >= 4]
#=vernum= sscanf(version(),"#f")
if vernum[1] >= 4
  global verbose
  plotcmd="plot(x[:,1],y,''+'',x[:,1],f); figure(gcf())"
else
  plotcmd="plot(x[:,1],y,''+'',x[:,1],f); shg"
end
if [exist("OCTAVE_VERSION")]
  global verbose
  plotcmd="plot(x[:,1],y,"+;data;",x[:,1],f,";fit;");"
end

if[exist("verbose")!=1], #If verbose undefined, print nothing
	verbose=0;       #This will not tell them the results
end

verbose = 0; #added 12/22/2016 to avoid version error()
if [nargin() <= 8], dFdp="dfdp"; end
if [nargin() <= 7], dp=.001*(pin*0+1); end; #DT
if [nargin() <= 6], wt=ones(length(y),1); end;	# SMB modification
if [nargin() <= 5], niter=20; end
if [nargin == 4], stol=.0001; end
=#
    dFdp = dfdp  #for Julia
#
#sizey = size(y)  #for debugging
y=y[:]; wt=wt[:]; pin=pin[:]; dp=dp[:]; #change all vectors to columns
#x=x[:] #added 12/21/2016
# check data vectors- same length()?
m=size(y)[1]  #rows
n=length(pin)
p=pin;
#(m1,m2)=size(x)  #m2 is not used (from Octave)
m1 = size(x)[1]   #rows, for Julia

if m1!=m
  error("input(x)/output[y] data must have same number of rows ")
end
#if nargin() <= 9
if false  #do not allow variable number of arguments
  options=[zeros(n,1), Inf()*ones(n,1)]
  nor = n; noc = 2
else
  (nor, noc)=size(options)
  if nor != n
    error("options and parameter matrices must have same number of rows")
  end
  if noc != 2
    options=[options[:,1], Inf()*ones(nor,1)]
  end
end
pprec=options[:,1]
maxstep=options[:,2]
#

# set up for iterations
#
#f=feval(F,x,p);
    f = F(x,p)  #applying the function directly in Julia
    fbest=f; pbest=p
r=wt.*(y-f)
ss=r'*r
sbest=ss
nrm=zeros(n,1)
#chgprev=Inf()*ones(n,1)
    chgprev=Inf*ones(n,1)  #for Julia
kvg=0
epsLlast=1
epstab=[.1, 1, 1e2, 1e4, 1e6]

# do iterations
iter = 0  #define here to allow access outside of loop, 3/24/2019
for iter=1:niter
        #println("iter: ", iter)
  #global iter  #not needed?, already defined, 10/26/2018
  ###fflush[stdout[]]; #***********************************************************
  #pbest  #un-comment for debugging
println("iter: ", iter)
#println("pbest: ", pbest)
  pprev=pbest
  #prt=feval(dFdp,x,fbest,pprev,dp,F)
        #prt = dfdp(x,fbest,pprev,dp,F)  #for Julia
        prt = dFdp(x,fbest,pprev,dp,F)  #for Julia 9/9/2018
            #prt = dfdp(x,fbest,pprev,dp,fpoly)  #for Julia, test only
            prt = prt*2  #for Julia, test only -- factor of 2 needed in Julia, why?
        #println("prt: ",prt)
        #error("deliberate stop")
  r=wt.*(y-fbest)
  sprev=sbest
  global sgoal
  sgoal=(1-stol)*sprev
  for j=1:n
    if dp[j]==0
      nrm[j]=0
    else
      prt[:,j]=wt.*prt[:,j]
      nrm[j]=prt[:,j]'*prt[:,j]
      if nrm[j]>0
        nrm[j]=1/sqrt(nrm[j])
      end
    end
    prt[:,j]=nrm[j]*prt[:,j]
  end
# above loop could ? be replaced by:
# prt=prt.*wt[:,ones(1,n)]
# nrm=dp./sqrt(diagm(prt'*prt))
# prt=prt.*nrm[:,ones(1,m)]'
  #(prt,s,v)=svd(prt,0)  #with two arguments: "econ" version
        #println("prt: ",prt)
        #println(size(prt))
        (prt,s,v)=svd(prt)  #for Julia (default is "thin")
            # Note: svd() in Julia returns only the diagonal elements of s, unlike Octave
            #  which returns a diagonal matrix. This led to a bug that was hard to track down.
        #println("v: ",v)
        #println("prt: ",prt)
        #println(size(prt))
        #println("s: ",s)
        #println("r: ", r)
  g=prt'*r
        #println("g: ",g)
  for jjj=1:length(epstab)
    global epsL
    epsL = max.(epsLlast*epstab[jjj],1e-7) #max dotted for Julia
            #println("epsL: ",epsL)
    #se=sqrt((s.*s)+epsL)
            #se=sqrt.((s.*s)+epsL)  #for Julia
            se=sqrt.((s.*s).+epsL)  #for Julia 10/30/2018
            #println("se: ",se)
    gse=g./se
            #println("gse: ",gse)
            #println("nrm: ",nrm)
    global chg
    chg=((v*gse).*nrm)
#   check the change constraints and apply as necessary
    ochg=chg
            #println("chg: ",chg)
            #println("maxstep: ",maxstep)
    idx = .~isinf.(maxstep)
            #println("idx: ",idx)
            #println("pprev: ",pprev)
    limit = abs.(maxstep[idx].*pprev[idx])  #dotted abs for Julia
            #println("limit: ",limit)
    chg[idx] = min.(max.(chg[idx],-limit),limit) #min and max dotted for Julia

###verbose
###any(ochg ~= chg)
    #=
    if verbose && any(ochg != chg)
      display(["Change in parameter[s]: ",
         #sprintf("#d ',find(ochg != chg)), 'were constrained"])
         sprintf("%d ',findall(x->x!=chg, ochg)), 'were constrained"]) #11/14/2018
    end
    =#
    aprec=abs.(pprec.*pbest);       #---
# ss=scalar sum of squares=sum((wt.*(y-f))^2).
            #println("aprec: ",aprec)
            #println("chg :",chg)
    if any(abs.(chg') .> 0.1*aprec)   #---  # only worth evaluating function if
      p=chg+pprev;                       # there is some non-miniscule change
      #f=feval(F,x,p)                         #chg tranposed, dotted > in if test, for Julia
      f=F(x,p)  #for Julia
      r=wt.*(y-f)
      ss=r'*r
                #println("ss: ",ss)
                #println("sbest :",sbest)
      #if ss<sbest
      if ss[1]<sbest[1] #Converts 1x1 arrays to floats, for Julia
        pbest=p
        fbest=f
        sbest=ss
      end
      #if ss<=sgoal
      if ss[1]<=sgoal[1] #Converts 1x1 arrays to floats, for Julia
        break
      end
    end                          #---
  end
        #println("epsL: ",epsL)
  epsLlast = epsL
        verbose = false  #not expected to be used for Julia
  if verbose
    eval(plotcmd)
  end
  if ss[1]<eps()    #conversion of 1x1 array to a float, for Julia
    break
  end
  aprec=abs.(pprec.*pbest)  #dotted for Julia
#  [aprec, chg, chgprev]
  if all(abs.(chg) .< aprec) && all(abs.(chgprev) .< aprec)  #abs dotted for Julia
    kvg=1
    if [verbose]
      @sprintf("Parameter changes converged to specified precision\n")
    end
    break
  else
    chgprev=chg
  end
  if ss[1]>sgoal[1]  ##conversion of 1x1 arrays to a floats, for Julia
    break
  end
end
    #println("end iterations")  #for debugging
# set return values
#
p=pbest
f=fbest
ss=sbest
kvg=((sbest[1]>sgoal[1])|(sbest[1]<=eps())|kvg)  #conversion of 1x1 arrays to floats
if kvg != 1  display(" CONVERGENCE NOT ACHIEVED! ")
else
    println("Convergence achieved")
end
# CALC VARIANCE COV MATRIX AND CORRELATION MATRIX OF PARAMETERS
# re-evaluate the Jacobian at optimal values
#jac=feval(dFdp,x,f,p,dp,F)
## 9/9/2018 4 lines commented out below
jac = dFdp(x,f,p,dp,F)  #for Julia
msk = dp != 0
n = sum(msk);           # reduce n to equal number of estimated parameters
#jac = jac[:, msk];	# use only fitted parameters
jac = jac[:, Int(msk)];	# use only fitted parameters

    #println("pbest: ",pbest)
    #println("f :",f)

## following section is Ray Muzic's estimate for covariance and correlation
## assuming covariance of data is a diagonal matrix proportional to
## diagm(1/wt.^2).
## cov matrix of data est. from Bard Eq. 7-5-13, and Row 1 Table 5.1
#=
if false
#if exist("sparse()")  # save memory()
  Q=sparse(1:m,1:m,1./wt.^2)
  Qinv=sparse(1:m,1:m,wt.^2)
else
  Q=diagm((0*wt+1)./(wt.^2))
  Qinv=diagm(wt.*wt)
end
resid=y-f
#resid
sizeresid = size(resid)
sizey = size(y)
sizef = size(f)
                                    #un-weighted residuals
covr=resid'*Qinv*resid*Q/(m-n);                 #covariance of residuals
Vy=1/(1-n/m)*covr;  # Eq. 7-13-22, Bard         #covariance of the data

jtgjinv=inv(jac'*Qinv*jac);			#argument of inv may be singular
covp=jtgjinv*jac'*Qinv*Vy*Qinv*jac*jtgjinv; # Eq. 7-5-13, Bard #cov of parm est
d=sqrt(abs(diagm(covp)))
corp=covp./(d*d')

if false  #sparse exists in Julia, but let's first try without
#if exist("sparse()")
  covr=spdiags(covr,0)
  stdresid=resid./sqrt(spdiags(Vy,0))
else
  covr=diagm(covr);                 # convert returned values to compact storage
  stdresid=resid./sqrt(diagm(Vy));  # compute then convert for compact storage
end

  ###sizerresid = size(resid)
  ###sizeQinv = size(Qinv)
  ###sizejac = size(jac)
  ###size((m-n)*jac'*Qinv*jac)
  ###size(n*resid'*Qinv*resid)
Z=((m-n)*jac'*Qinv*jac)/(n*resid'*Qinv*resid)

### alt. est. of cov(). mat. of parm.:(Delforge, Circulation, 82:1494-1504, 1990
##display("Alternate estimate of cov(). of param. est.")
##acovp=resid'*Qinv*resid/(m-n)*jtgjinv

#Calculate R^2 [Ref Draper & Smith p.46]
#
   #removed 12/22/2016-DJT
   ### r=corrcoef([y[:],f[:]])
   ### r2=r[1,2].^2

# if someone has asked for it, let them have it
#
if [verbose]
  eval(plotcmd)
  display(" Least Squares Estimates of Parameters")
  display(p')
  display(" Correlation matrix of parameters estimated")
  display(corp)
  display(" Covariance matrix of Residuals" )
  display(covr)
  display(" Correlation Coefficient R^2")
  display(r2)
  sprintf(" 95## conf region: F[0.05](#.0f,#.0f)>= delta_pvec''*Z*delta_pvec",n,m-n)

#   runs test according to Bard. p 201.
  n1 = sum((f-y) < 0)
  n2 = sum((f-y) > 0)
  nrun=sum(abs(diff((f-y)<0)))+1
  if (n1>10)&(n2>10) # sufficent data for test?
    zed=(nrun-(2*n1*n2/(n1+n2)+1)+0.5)/(2*n1*n2*(2*n1*n2-n1-n2)
      /((n1+n2)^2*(n1+n2-1)))
    if [zed < 0]
      prob = erfc(-zed/sqrt(2))/2*100
      display([num2str(prob),"# chance of fewer than ',num2str(nrun),' runs."])
    else
      prob = erfc(zed/sqrt(2))/2*100
      display([num2str(prob),"# chance of greater than ',num2str(nrun),' runs."])
    end
  end
end
=#
#return (f,p,kvg,iter,corp,covp,covr,stdresid,Z,r2)
    return (f,p,kvg,iter)
end
#end

#************************************************

function leasqr_test()
#Test of leasqr

   p = [1.,2.,3.] #polynomial coefficients
   fpoly(x,p) = p[1].*x.^2 + p[2].*x + p[3]
   #x = [1 2 3 4 5]'*0.1 + 1
   #x = rand(1000,1)*0.1 + 1
   x = rand(1000,1)*0.1 .+ 1  #11/14/2018
   y = fpoly(x,p)

   #leasqr(x,y,pin,F,stol,niter,wt,dp,dFdp,options)
   pin = [.9,2.1,2.9]
   stol = .001;     #scalar tolerance on fractional improvement in least squares fit
   niter = 20;      #maximum iterations for least squares fit
   wts = ones(length(y),1);   #default
   ADJ = .001
   dp1 = [ADJ,ADJ,ADJ];
   #lim = 1000
   #options = [0 lim; 0 lim; 0 lim];
   options = [0 Inf; 0 Inf; 0 Inf];

   (f,p,kvg,iter) = leasqr(x,y,pin,fpoly,stol,niter,wts,dp1,dfdp,options)
   println("polynomial coefficients: ", p)
   println("should be: 1., 2., 3.")
end


#USAGE:
#   P1blockcal(sbdir,serialnum,runnum,sensornum,slitcen=0.)
#   P1blockcal_test(sensornum)

#*************************************************************

function P1readparameters(dirin,sensor)

    f = dirin * "_ParametersS" * string(sensor) * ".txt"
    fid = open(f, "r")
    P = readlines(fid)
    close(fid)

    (sx,sy,sz) = parse.(Float64, split(P[2])[2:4])
    (su,sv,sw) = parse.(Float64, split(P[3])[2:4])
    (gu,gv,gw) = parse.(Float64, split(P[4])[2:4])
    (ccdx,ccdy,ccdz) = parse.(Float64, split(P[5])[2:4])
    (ccdu,ccdv,ccdw) = parse.(Float64, split(P[6])[2:4])
    x0 = parse.(Float64, split(P[7])[2])
    tglass = parse.(Float64, split(P[8])[2])
    flen = parse.(Float64, split(P[9])[2])
    nglass = parse.(Float64, split(P[10])[2])
    #return (sx,sy,sz,su,sv,sw,gu,gv,gw,ccdx,ccdy,ccdz,ccdu,ccdv,ccdw,x0)
    return (sx,sy,sz,su,sv,sw,gu,gv,gw,ccdx,ccdy,ccdz,ccdu,ccdv,ccdw,flen,x0,tglass,nglass)
end

#*************************************************************

function P1getDataSX(sensor,dirin,fname="Cal32out2.dat",zoffset=0.)
#Get that portion of the calibration data(from Cal32out2.dat) seen by a given sensor.
#DJT--June 2017
#(similar to P1readCalDat except for line 20)
#returns an N by 4 matrix xyzc
#xyzc is a column vector with repeated sequence x,y,z,centroid

fid = open(dirin * fname)  #* is the operator for string concatenation
         #   V = sscanf(s,"c#d #f #f #f #f #f #f")  #   #if V[5]>0 && V[6]>0 && V[7]>0  #all seen, not -999;  #-999 rejection turned off--3 Feb 2017
M = readlines(fid)
    M = M[4:end]  #to skip first three lines
close(fid)

#idx = Array{Bool}(false,1,size(M)[1])
    idx = Array{Bool}(undef,size(M)[1])
    idx .= false
    #return idx
for i in eachindex(M[4:end])
	j = i + 1
	idx[i] = !contains(split(M[i])[sensor+4], "-999")
end
M = M[idx]

R = Meta.parse.(split(M[1])[2:end])
n = size(M)[1]

for k in 2:n
        S = split(M[k])
        P = Meta.parse.(Float64, S[2:end])  #ignoring first column
        R = [R P]
end
    R = R'
        #return R
    col = sensor + 3
    xyzc = [R[:,1:3] R[:,col]]
    return xyzc
end

#*************************************************************

## Author: Dennis <Dennis@DESKTOP-TRI7NMC>
## Created: 2017-06-08

function P1calcRBFcoeffs(dirin, sensor , miter, filter, XYZC, Tukey)
 ##Note: parameter XYZC is not used

  global P1c2plane_e  #handle to the function(), added 2/15/2018
  global mmPerPixel
  global CNORM
  global snell
  global rbfname
  global ep
  global omega

  tol = 0.0001
  RBF_on = 1

  #zoffset = 0
  #XYZC = P1getDataS1S2S3[sensor,dirin,"Cal32out2.dat",zoffset]
  XYZC = P1getDataSX(sensor,dirin)
  npts = size(XYZC)[1]
    println("XYZC: ", XYZC)

  (sx,sy,sz,su,sv,sw,gu,gv,gw,ccdx,ccdy,ccdz,ccdu,ccdv,ccdw,flen,x0,tglass,nglass) = P1readparameters(dirin,sensor)

   npts = size(XYZC)[1]  #rows
   s = zeros(npts)
   v = zeros(npts)
   for i in 1:npts

      x = XYZC[i,1]
      y = XYZC[i,2]
      z = XYZC[i,3]

      #Set the initial target point.
      ptarg = [x, y, z]

      #Find two points along the slit.
      mag = norm([su,sv,sw])
      p1 = [sx, sy, sz] - 10*[su, sv, sw]/mag
      p2 = [sx, sy, sz] + 10*[su, sv, sw]/mag

      kvg = 0;  #remains zero if convergence fails
      for k in 1:miter

         #Find the intersection of the plane defined by p1,p2,ptarg and the line defined by the CCD.
         (a,b,c,d) = P1planethru3pts( ptarg[1],ptarg[2],ptarg[3], p1[1],p1[2],p1[3],p2[1],p2[2],p2[3] )
         pn = [a,b,c]/norm([a,b,c]);  #normal to plane
         ppt = [p1[1],p1[2],p1[3]];  #a point on the plane
         ldir = [ccdu,ccdv,ccdw]/norm([ccdu,ccdv,ccdw]);  #direction of the CCD
         lpt = [ccdx,ccdy,ccdz];  #a point on the CCD
         ipt = P1intsectlineplane(ppt,pn,lpt,ldir)  #the intersection point, 9/13/2018
         global ipt  #to retain value after loop exit, 9/13/2018

         #Find the centroid corresponding to ipt.
         if sensor==1 || sensor==3
            centroid = ( ipt[1] - ccdx )/( mmPerPixel*ccdu/norm([ccdu,ccdv,ccdw]) ) + x0
         else  #sensor==2
            centroid = ( ipt[2] - ccdy )/( mmPerPixel*ccdv/norm([ccdu,ccdv,ccdw]) ) + x0
         end
         global centroid  #to retain value after loop exit, 9/13/2018

         #Find the exit plane corresponding to the targeting centroid.
         #[a,b,c,d] = P1c2plane_f[x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass,nglass, mmPerPixel]
         #Changed to version _e 4/28/2017 for test purposes
         (a,b,c,d) = P1c2plane_e(x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz,
                ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass,nglass, mmPerPixel)
         offset = P1distancefromplane( x,y,z, a,b,c,d )
         errvec = -1.0*offset*[a, b, c]/norm([a, b, c])

         if norm(errvec)<tol
            kvg = 1
            break
         end

         #Set a new target point.
         ptarg = ptarg - errvec

      end
      if kvg==1  #if mirror-image()-targeting converged
         #Calculate and save the centroid correction.
         s[i] = centroid - XYZC[i,4];  #centroid correction in pixels
         #Calculate and save the cosine of the ray along the direction of the slit.
         vray = ([x,y,z] - ipt)/norm([x,y,z] - ipt);  #direction from CCD intercept to CMM point
         #println("vray: ", vray)
         v[i] = CNORM*dot( vray, [su,sv,sw]/norm([su,sv,sw]) )
         #println("v: ", v[i])
         #global v  #to retain value after loop exit, 9/13/2018
         if RBF_on==0
             s[i] = 0;  #all corrections become zero
         end
         #global s  #to retain value after loop exit, 9/13/2018
      else
         s[i] = 0
         @printf("mirror image targeting failed to converge for point %s\n",string(i))
      end

   end
   println(s)

   #******************************************
   # Write the centroid corrections to a file.
   #******************************************
      fout = dirin * "_RBFmapS" * string(sensor) * ".txt"
      fid = open(fout,"w")
      @printf(fid,"%s\n", dirin * "_RBFmapS" * string(sensor) * ".txt");  #first line is a header
      Map = zeros(1,3)
      k = 1
      for i in 1:npts
         Map[k,1] = (XYZC[i,4] - x0)*mmPerPixel
         Map[k,2] = v[i]
         Map[k,3] = s[i]*mmPerPixel
         @printf(fid,"%12.4f %12.8f %12.8f\n", Map[k,1],Map[k,2],Map[k,3])
      end
      close(fid)

   RBF_plot_on = 0
   (ck,Ctrs,efitted) = P1createRBF(dirin, sensor, omega, rbfname, ep, RBF_plot_on, filter, Tukey)

   return(ck,Ctrs,efitted)
#endfunction

end

#*********************************************************

function P1cross(ux,uy,uz,vx,vy,vz)
#function s = P1cross(u,v)
#Find the cross product of two vectors.
#(faster than Octaves built-in function())

sx = uy*vz - uz*vy
sy = uz*vx - ux*vz
sz = ux*vy - uy*vx
return (sx,sy,sz)

end

#********************************************************

function P1Snell(vinc,vnorm,n1,n2)
#Calculate the direction of a refracted ray by Snell's law.
#DJT-January 2017

# n1 is the index of refraction of the external medium[generally air].
# n2 is the index of refraction of the refracting medium.
# vinc is unit vector giving the direction of the incident ray.
# vnorm is a unit vector normal to the refracting surface and directed towards the incoming ray.
# vrefrac is a unit vector giving the direction of the refracted ray.

r = n1/n2
c = -dot(vnorm,vinc)
#c = -c
#c = abs(c)
#vinc = -vinc
vrefrac = r*vinc + (r*c - sqrt(1 - r*r*(1-c*c)))*vnorm

return vrefrac
end

#********************************************************

function P1readCalDat(dirin,fname="Cal32out2.dat",zoffset=0.)
#Get that portion of the calibration data[from Cal32out.dat] seen by all three sensors.
#DJT--February 2017

   #xyzc is a column vector with repeated sequence x,y,z,centroid

   fid = open(dirin * fname)  #* is the operator for string concatenation
              #   V = sscanf(s,"c#d #f #f #f #f #f #f")    #   #if V[5]>0 && V[6]>0 && V[7]>0  #all seen, not -999;  #-999 rejection turned off--3 Feb 2017
   M = readlines(fid)
       M = M[4:end]  #to skip first three lines
   close(fid)
   #idx = Array{Bool}(false,1,size(M)[1])
   println("size(M) with -999s: ", size(M))
   #flush(stdout)  #removed 10/29/2018 DJT
       idx = Array{Bool}(undef, size(M)[1])
       idx .= false
       #return idx
   for i in eachindex(M[4:end])
   	   j = i + 1
   	   #idx[i] = !contains(M[i], "-999")
   	   idx[i] = occursin("-999", M[i])
   end
   #M = M[idx]  #11/20/2018
   #println("size(M) without -999s: ", size(M))

   #R = parse.(split(M[1])[2:end])
   R = Meta.parse.(split(M[1])[2:end])
   n = size(M)[1]

   for k in 2:n
       S = split(M[k])
       P = parse.(Float64, S[2:end])  #ignoring first column
       #P = Meta.parse.(Float64, S[2:end])  #ignoring first column
       R = [R P]
   end
       println("size(R'): ", size(R'))
       return R'
end

#************************************************************

function P1modelS1e( XYZC, p )
#Model of sensors--to be optimized with the Levenberg-Marquardt algorithm
#DJT--January 2017
#DJT--July 16, 2018

global P1c2plane_e #handle to the function(), added 2/15/2018
global sensor  #(1,2, or 3)

global sx
global sy
global sv
global su
global sv
global sw
global ccdx
global ccdy
global ccdz
global ccdu
global ccdv
global ccdw
global gu
global gv
global gw
global mmPerPixel
global flen
global x0
global tglass
global nglass
global snell


#adjustable parameters
if sensor==2
   sy = p[1]
   sz = p[2]
   sv = p[3]
   sw = p[4]
   gv = p[5]
   gu = -(sv*gv + sw) #Note: [1 sv sw].[gu gv 1] = |s||g|cos(theta) = 0
   flen = p[6]
   ccd = [sx,sy,sz] + flen*[gu,gv,gw]/norm([gu,gv,gw])
   ccdx = ccd[1]
   ccdy = ccd[2]
   ccdz = ccd[3]
   x0 = p[7]
   ccdw = p[8]
   ccdu = p[9]
   tglass = p[10]
   #par1 = [sy,sz,sv,sw,gv,flen,x0,ccdw,ccdu,tglass]
else #sensors 1 or 3
   sx = p[1]
   sz = p[2]
   su = p[3]
   sw = p[4]
   gu = p[5]
   gv = -(su*gu + sw) #Note: [1 sv sw].[gu gv 1] = |s||g|cos(theta) = 0
   flen = p[6]
   ccd = [sx,sy,sz] + flen*[gu,gv,gw]/norm([gu,gv,gw])
   ccdx = ccd[1]
   ccdy = ccd[2]
   ccdz = ccd[3]
   x0 = p[7]
   ccdw = p[8]
   ccdv = p[9]
   tglass = p[10]
   #par1 = [sx,sz,su,sw,gu,flen,x0,ccdw,ccdv,tglass]
end



n = size(XYZC)[1]  #rows
offset = zeros(n)  #added 9/8/2018
for i in 1:n
   x = XYZC[i,1]
   y = XYZC[i,2]
   z = XYZC[i,3]
   centroid = XYZC[i,4]
   #[a,b,c,d] = P1c2plane_f[ x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass,nglass, mmPerPixel ]
   #Changed to version e 4/28/2017 for test purposes
   (a,b,c,d) = P1c2plane_e( x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass,nglass, mmPerPixel )
   #[a,b,c,d] = P1c2plane_g[ x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass,nglass, mmPerPixel ]
      #Note: Version e replaced with version g, January 24, 2018 -- uses two exit rays to determine the plane.
   offset[i] = P1distancefromplane( x,y,z,a,b,c,d )
end
#ffset = offset'  #removed 9/10/2018
return offset
end

#*********************************************************

function P1planethru3pts(x1,y1,z1, x2,y2,z2, x3,y3,z3)
#Find the equation of a plane passing through three points[eq: a*x + b*y + c*z + d = 0]
#DJT--January 1017

###D = det([x1,x2,x3; y1,y2,y3; z1,z2,z3])
###d = 1
###a = (-d/D)*det([1,1,1; y1,y2,y3; z1,z2,z3])
###b = (-d/D)*det([x1,x2,x3; 1,1,1; z1,z2,z3])
###c = (-d/D)*det([x1,x2,x3; y1,y2,y3; 1,1,1])


#Revision DJT January 2018

### v1 = [x1,y1,z1] - [x2,y2,z2]
### v2 = [x1,y1,z1] - [x3,y3,z3]
### vn = cross(v1,v2)
### a = vn[1]
### b = vn[2]
### c = vn[3]

#[a,b,c] = P1cross[x1-x2,y1-y2,z1-z2, x1-x3,y1-y3,z1-z3]; #revised for speed 7/3/2018 DJT (Octave only)
(a,b,c) = P1cross(x1-x2,y1-y2,z1-z2, x1-x3,y1-y3,z1-z3)

d = -(a*x1 + b*y1 + c*z1)

##Check:  (uncomment for debugging)
#a*x1 + b*y1 + c*z1 + d
#a*x2 + b*y2 + c*z2 + d
#a*x3 + b*y3 + c*z3 + d
#d1 = P1distancefromplane[x1,y1,z1,a,b,c,d]
#d2 = P1distancefromplane[x2,y2,z2,a,b,c,d]
#d3 = P1distancefromplane[x3,y3,z3,a,b,c,d]

return (a,b,c,d)
end

#*************************************************************

function P1intsectlineplane(ppt,pn,lpt,ldir)
#Find the intersection of a line and a point.
#DJT--January 2017
#DJT--July 16, 2018, converted from Octave to Julia

#   ppt is a point on the plane.
#   pn is a unit vector normal to the plane.
#   lpt is a point on the line()
#   ldir a unit vector giving the direction of the line()

pn = pn/norm(pn)
ldir = ldir/norm(ldir)
if abs(dot(ldir,pn))<eps()
   error("the line and the plane do not intersect in a unique point")
else
   t = dot((ppt - lpt),pn)/dot(ldir,pn)
   point = lpt + t*ldir
end
return point
end

#*************************************************************

function P1c2plane_e( x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass,
        nglass, mmPerPixel )
#Given block paramters and a centroid calculate a refracted exit plane parallel to the incident plane.
# DJT-January 2017
# DJT-July 16, 2018, converted from Octave to Julia

#INPUT:
#   [x,y,z] are the coordinates of the emitter.
#   centroid is the center-of-gravity of the image on the CCD.
#   [sx,sy,sz] is the center of the slit.
#   [su,sv,sw] is a unit vector along the slit.
#   [ccdx,ccdy,ccdz] is a point on the CCD.
#   [ccdu,ccdv,ccdw] is a unit vector along the CCD.
#   snell is a flag to turn on the Snell effect corrections[1 for "on" and 0 for "off"].
#   [gu,gv,gw] is a unit vector normal to the glass
#   tglass is the thickness of the glass[in mm].

#Calculate the x,y,z[in mm] of the centroid along the CCD unit vector.
#mmPerPixel = 0.014
mag = norm([ccdu,ccdv,ccdw])
xc = ccdx + (centroid - x0)*mmPerPixel*ccdu/mag
yc = ccdy + (centroid - x0)*mmPerPixel*ccdv/mag
zc = ccdz + (centroid - x0)*mmPerPixel*ccdw/mag

#Find two points along the slit.
mag = norm([su,sv,sw])
p1x = sx - 10*su/mag
p1y = sy - 10*sv/mag
p1z = sz - 10*sw/mag
p2x = sx + 10*su/mag
p2y = sy + 10*sv/mag
p2z = sz + 10*sw/mag


#Place a third point on the slit. (varies with point location)
ys = 0
p3x = sx + ys*su/mag
p3y = sy + ys*sv/mag
p3z = sz + ys*sw/mag

if snell==1
   vnorm = [gu,gv,gw]/norm([gu,gv,gw])
#println("vnorm: ", vnorm)
   #find direction of the center refracted ray
#println("[p3x p3y p3z]: ", [p3x p3y p3z])
#println("[xc yc zc]: ", [xc yc zc])
   vinc1 = [p3x, p3y, p3z] - [xc, yc, zc]
   #(vx,vy,vz) = [p3x, p3y, p3z] - [xc, yc, zc]
   #vinc1 = [vx,vy,vz]/sqrt(vx*vx + vy*vy + vz*vz)
   #vinc1 = vinc1/norm(vinc1)
   #vinc1 = vinc1/P1norm[p3x-xc, p3y-yc, p3z-zc];  #revised for speed 7/3/2018 (Octave only)
   vinc1 = vinc1/sqrt( (p3x-xc)^2 + (p3y-yc)^2 + (p3z-zc)^2 )  #for Julia 11/1/2018
#println("vinc1: ", vinc1)
#println("nglass: ", nglass)
   vrefrac1 = P1Snell(vinc1,vnorm,1.000,nglass)  #(a unit vector)
#println("vrefrac1: ", vrefrac1)
   #find the exit point for the center ray
   #exitpoint = [p3x, p3y, p3z] + vrefrac1*tglass/abs(dot(vinc1,vnorm))
   exitpoint = [p3x, p3y, p3z] + vrefrac1*tglass/abs(dot(vrefrac1,vnorm))
#println("exitpoint: ", exitpoint)
   #calculate the plane of the incident rays
   (a,b,c,d) = P1planethru3pts( xc,yc,zc, p1x,p1y,p1z, p2x,p2y,p2z )
   #find a plane parallel to the incident plane and passing through the exit point
   d2 = -( a*exitpoint[1] + b*exitpoint[2] + c*exitpoint[3] )
   #returns [a,b,c,d2]
   return (a,b,c,d2)
else #(No Snell effect corrections)
   (a,b,c,d) = P1planethru3pts( xc,yc,zc, p1x,p1y,p1z, p2x,p2y,p2z )
   return (a,b,c,d)
end

end

#******************************************************************

function P1distancefromplane(x,y,z,a,b,c,d)
#Calculate the distance of a point from a plane.
#DJT--January 2017
#DJT--July 16, 2018, converted from Octave to julia

distance = (a*x + b*y + c*z + d)/norm([a,b,c])
#distance = (a*x + b*y + c*z + d)/P1norm[a,b,c];  #revised for speed 7/3/2018 (Octave only)
return distance
end

#*****************************************************************

function xP1c2plane_e_orig( x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass,
        nglass, mmPerPixel )
#Given block paramters and a centroid calculate a refracted exit plane parallel to the incident plane.
# DJT-January 2017
# DJT-July 16, 2018, converted from Octave to Julia

#INPUT:
#   [x,y,z] are the coordinates of the emitter.
#   centroid is the center-of-gravity of the image on the CCD.
#   [sx,sy,sz] is the center of the slit.
#   [su,sv,sw] is a unit vector along the slit.
#   [ccdx,ccdy,ccdz] is a point on the CCD.
#   [ccdu,ccdv,ccdw] is a unit vector along the CCD.
#   snell is a flag to turn on the Snell effect corrections[1 for "on" and 0 for "off"].
#   [gu,gv,gw] is a unit vector normal to the glass
#   tglass is the thickness of the glass[in mm].

#Calculate the x,y,z[in mm] of the centroid along the CCD unit vector.
#mmPerPixel = 0.014
mag = norm([ccdu,ccdv,ccdw])
xc = ccdx + (centroid - x0)*mmPerPixel*ccdu/mag
yc = ccdy + (centroid - x0)*mmPerPixel*ccdv/mag
zc = ccdz + (centroid - x0)*mmPerPixel*ccdw/mag

#Find two points along the slit.
mag = norm([su,sv,sw])
p1x = sx - 10*su/mag
p1y = sy - 10*sv/mag
p1z = sz - 10*sw/mag
p2x = sx + 10*su/mag
p2y = sy + 10*sv/mag
p2z = sz + 10*sw/mag

#Place a third point on the slit. (varies with point location)
ys = 0
p3x = sx + ys*su/mag
p3y = sy + ys*sv/mag
p3z = sz + ys*sw/mag

if snell==1
   vnorm = [gu,gv,gw]/norm([gu,gv,gw])
   #find direction of the center refracted ray
   vinc1 = [p3x, p3y, p3z] - [xc, yc, zc]
   vinc1 = vinc1/norm(vinc1)
   #vinc1 = vinc1/P1norm[p3x-xc, p3y-yc, p3z-zc];  #revised for speed 7/3/2018 (Octave only)

   vrefrac1 = P1Snell(vinc1,vnorm,1.000,nglass)  #(a unit vector)
   #find the exit point for the center ray
   #exitpoint = [p3x, p3y, p3z] + vrefrac1*tglass/abs(dot(vinc1,vnorm))
    exitpoint = [p3x, p3y, p3z] + vrefrac1*tglass/abs(dot(vrefrac1,vnorm))
   #calculate the plane of the incident rays
   (a,b,c,d) = P1planethru3pts( xc,yc,zc, p1x,p1y,p1z, p2x,p2y,p2z )
   #find a plane parallel to the incident plane and passing through the exit point
   d2 = ( a*exitpoint[1] + b*exitpoint[2] + c*exitpoint[3] )
   return (a,b,c,d2)
else() #(No Snell effect corrections)
   (a,b,c,d) = P1planethru3pts( xc,yc,zc, p1x,p1y,p1z, p2x,p2y,p2z )
   return (a,b,c,d)
end
end  #end P1c2plane_e_orig

#*******************************************************************

function P1blockcal(sbdir,serialnum,runnum,sensornum,slitcen=0.)
#Calibrate a sensor block.
#DJT--January 2017"
#DJT--Revised February 2018   Handle to function P1c2plane_e added
#DJT--6/14/2018 added 400mm sensor bar option

#INPUT:
#   sbdir is the directory holding the sensor bar (CMM) data.
#   serialnum is the serial number of the sensor bar().
#   runnum is the run number for the CMM data.
#   sensornum is the number of the sensor block [1, 2,or 3].

#FILE OUTPUT:
#   file "_ParametersSX.txt" (where X is 1, 2, or 3) holds the calculated parameters for this block.
#   file "_RBFmapSX.txt holds the centroid error corrections.

#USES:
#   XYZC = P1getDataSX[dirin,sensor] reads data from file "Call32out.dat" for sensor 1, 2 or 3 into matrix XYZC.
#   leasqr[...] is an Octave implementation of the Levenberg-Marquardt least squares fitting algorithm.
#   [a,b,c,d] = P1planethru3pts[...] constructs a plane though three given points.
#   ipt = P1intsectlineplane[...] finds the point of intersection of a line and a plane.
#   offset = P1distancefromplane[...] returns the [+/-] offset of a point from a plane.
#   [a,b,c,d] = P1c2plane_f[...] constructs an exit plane given a set of block parameters and a centroid.

#profile on  #added 7/2/2018

page_screen_output = 0
#more off
#warning("off", "All")
#dirin = strcat(sbdir,string(serialnum),"/Run',string(runnum),'/")
    dirin = sbdir * string(serialnum) * "/Run" * string(runnum) * "/"  #for julia
#if nargin()<5
   #slitcen = 0  #default argument value in Julia
#end

#**********************************************************************
#Global variables are used to pass parameters to the least squares fit.
#**********************************************************************
global P1c2plane_e      #handle to the function()
global P1calcRBFcoeffs  #handle to the function()
global sensorbartype    #900mm, 880mm, etc.
global zoffset          #shifts CMM z-coordinates by this amount
global CNORM            #normalization constant for vertical angles

#Control parmeters
global snell       #Snell corrections flag [1->"on",0->"off]
global stol        #scalar tolerance on fractional improvement in least squares fit
global niter       #maximum iterations for least squares fit
global tol         #convergence tolerance for mirror-imaging error vector [in mm]
global miter       #maximum iterations for mirror-imaging
global rbfname     #selects the radial basis function ("imq"->inverse multiquadrics)
global omega       #rbf smoothing paramters--larger gives less smoothing
global ep          #rbf shape paramter
global RBF_on      #to turn off RBF corrections
global RBF_plot_on #turn off for speed
global filter1     #filter level for base calibration outlier removal
global Tukey1      # ==1 for Tukey filtering of outliers in base calibration
global filter2     #filter level for RBF outlier removal [=0 to turn off filtering]
global Tukey2      # ==1 for Tukey filtering of outliers in RBF creation

global sensor  #    selects sensor 1,2 or 3
               #   -----------------------------
global sx      #    slit position and direction
global sy      #
global sv      #
global su      #
global sv      #
global sw      #
               #   -----------------------------
global ccdx    #    a point on the CCD and its
global ccdy    #    its direction
global ccdz    #
global ccdu    #
global ccdv    #
global ccdw    #
               #   -----------------------------
global gu      #    a vector normal to the glass
global gv      #
global gw      #
               #   -----------------------------
global mmPerPixel # pixel spacing
global flen       # effective focal length()
global x0         # principal pixel
global tglass     # slit glass thickness [in mm]
global tdustglass # dust glass thickness [in mm]
global nglass     # index of refraction
               #   -----------------------------

#Control parmeters [default values]
pars = "DT"
CNORM = 34;      #normalization constant for vertical angles
snell = 1;       #Snell corrections flag [1->"on",0->"off]
stol = .001;     #scalar tolerance on fractional improvement in least squares fit
niter = 20;      #maximum iterations for least squares fit
tol = .0001;     #convergence tolerance for mirror-imaging error vector [in mm]
miter = 10;      #maximum iterations for mirror-imaging
  miter = 20;      #maximum iterations for mirror-imaging  #9/19/2018
filter1 = 3.0    #filter level for base calibration outlier removal
filter2 = 5.0;   #filter level for RBF outlier removal [=0 to turn off filtering]
  filter2 = 3.0    ####setting RBF filter at the same level as the block cal filter()--3/16/2018################
#rbfname = "imq';#selects the radial basis function ('imq"->inverse multiquadrics)
#omega = 10000;  #rbf smoothing paramters--larger gives less smoothing
#ep = 8;         #rbf shape paramter
RBF_on = 1;      #to turn off RBF corrections
RBF_plot_on = 0; #turn off for speed
Tukey1 = 0
Tukey2 = 0
snell = 1  #added 9/20/2018

#Change some of the global control parameters from the default values.
#P1read_Ini_File[dirin]  #deactivated 9/6/2018

#*****************************************************************************
# Set starting values for all block parameters[some globals will be adjusted].
#*****************************************************************************
#nglass = 1.5098;  #index of refraction of N-BK7 at 850nm
nglass = N_SLIT_GLASS()
#tglass = 2.6
tglass = T_SLIT_GLASS()
tglass = tglass + tdustglass;

sensor = sensornum
println("\n***** Sensor " * string(sensor) * " *****\n")
   if sensorbartype == "900mm"   #(extrusion camera)
      theta = 0;        #yaw angle of end sensors [zero for extrusion sensor bars]
      barlength = 900;  #distance between end slits
      x0_13 = 6.0;      #principal pixel for sensors 1&3 [extrusion cameras]
   elseif sensorbartype == "800mm"   #(extrusion camera)
      theta = 0
      barlength = 800
      x0_13 = 6.0;
   elseif sensorbartype == "600mm"   #(extrusion camera)
      theta = 0
      barlength = 600
      x0_13 = 6.0;
   elseif sensorbartype == "880mm"   #(MIC-6 880)
      theta = 0;
      barlength = 880;
      x0_13 = 13.0  #11/21/2018 to match Octave version
      #x0_13 = 6.0  #11/13/2018 closer to converged value
   elseif sensorbartype == "880mm slanted"   #(MIC-6 880 with slanted CCDs in sensors 1&3)
      theta = 12.0;
      barlength = 880;
      x0_13 = 13.0
   elseif sensorbartype == "400mm"   #(extrusion camera) #added 6/14/2018
      theta = 0
      barlength = 400
      x0_13 = 6.0
    else
      error("sensor bar type not recognized")
    end

if sensor==2
   sx = slitcen;    #non-adjusting--pass in a starting value via argument 5
   sy = 0.01
   sz = 0.01

   su = 1;     #non-adjusting
   sv = 0.01
   sw = 0.01

   ccdx = 0;   #calculated
   ccdy = 0;   #calculated
   ccdz = 0;   #calculated

   ccdu = 0;
   ccdv = -1;  #non-adjusting
   ccdw = 0.01

   gu = 0;     #calculated
   gv = 0.01
   gw = 1;     #non-adjusting

   flen = 34
   mmPerPixel = 0.014
   x0 = 14/mmPerPixel
   #tglass = 2.5;

   #Designate the parameters to be adjusted.
   par1 = [sy,sz,sv,sw,gv,flen,x0,ccdw,ccdu,tglass]

else  #sensors 1&3
   sx = barlength/2
      if sensor==3
         sx = -sx
      end
   sy = slitcen;   #non-adjusting--pass in a starting value via argument 5
   sz = 0.01
println("barlength: ", barlength)
println("sx: ", sx)
#error("deliberate stop")

   su = cos((90 - theta)*pi/180) + 0.01;  #0.01 added to avoid a zero starting value
      if sensor==3
         su = -cos((90 - theta)*pi/180) + 0.01
      end
   sv = 1;   #non-adjusting
   sw = 0.01

   ccdx = 0; #calculated
   ccdy = 0; #calculated
   ccdz = 0; #calculated

   ccdu = 1; #non-adjusting
      if sensor==3
         ccdu = -1
      end
   ccdw = 0.01
   ccdv = 0

   gu = cos((90-theta)*pi/180) + 0.01;   #0.01 added to avoid a zero starting value
      if sensor==3
         gu = -cos((90 - theta)*pi/180) + 0.01
      end
   gv = 0;   #calculated
   gw = 1;   #non-adjusting

   flen = 34.1
   mmPerPixel = 0.014
   x0 = x0_13
   #x0 = x0_13/mmPerPixel  #for test 11/13/2018
   #tglass = 2.5

   #Designate the parameters to be adjusted.
   par1 = [sx,sz,su,sw,gu,flen,x0,ccdw,ccdv,tglass];
end

#**********************************************
# Adjust parameters for best least squares fit.
#**********************************************
#printf["Adjusting parameters...\n"]

#fid = fopen(strcat(dirin,"_RemovedS",string(sensor),".txt"), "w"); #to discard previous contents
    fid = open(dirin * "_RemovedS" * string(sensor) * ".txt", "w")  #for Julia 11/7/2018
close(fid)

#println("dirin: ", dirin)
#error("deliberate stop")

M = P1readCalDat(dirin,"Cal32out2.dat",zoffset);  #matrix with rows x,y,z,centroid
#Ix = find(M[:,3+sensornum]>0);  #to take out -999s
    Ix = findall(M[:,3+sensornum].>0);  #to take out -999s  #rev. 9/8/2018
XYZC = [M[Ix,1:3] M[Ix,3+sensornum]]

println("size(Ix): ", size(Ix))
println("size(XYZC): ", size(XYZC))  #for debugging 11/15/2018

#################################################################
for pass=1:3  #XYZC will be altered[filtered] after the first pass if outliers are present.
  global Ex
  global Ey
  global Ez

   n = size(XYZC)[1]  #rows
   y = zeros(n,1);  #column vector

   if pars=="DT"

      wts = ones(length(y),1);   #default
      pin = par1
      #options = [0 Inf(); 0 Inf(); 0 Inf(); 0 Inf(); 0 Inf(); 0 Inf(); 0 Inf(); 0 Inf(); 0 Inf(); 0 Inf()]
      options = [0 Inf; 0 Inf; 0 Inf; 0 Inf; 0 Inf; 0 Inf; 0 Inf; 0 Inf; 0 Inf; 0 Inf]  #for Julia 9/8/2018

      ADJ = .001
      FIX = 0
      if sensornum==2
        #println("[sy,sz,sv,sw,gv,flen,x0,ccdw,ccdu,tglass]")
        parnames = ("sy","sz","sv","sw","gv","flen","x0","ccdw","ccdu","tglass")  #added 9/10/2018
      else
        #println("[sx,sz,su,sw,gu,flen,x0,ccdw,ccdv,tglass]")
        parnames = (sx,sz,su,sw,gu,flen,x0,ccdw,ccdv,tglass)  #added 9/10/2018
        parnames = ("sx","sz","su","sw","gu","flen","x0","ccdw","ccdv","tglass")  #added 9/10/2018
      end
      #par1 = [sx,sz,su,sw,gu,flen,x0,ccdw,ccdv,tglass];  #template sensors 1&3
      #par1 = [sy,sz,sv,sw,gv,flen,x0,ccdw,ccdu,tglass];  #template sensor 2
      dp1 = [ADJ,ADJ,ADJ,ADJ,FIX,ADJ,ADJ,ADJ,FIX,FIX]
      dp2 = [ADJ,ADJ,ADJ,ADJ,ADJ,FIX,ADJ,ADJ,ADJ,FIX]
      dp3 = [FIX,FIX,FIX,FIX,FIX,ADJ,FIX,FIX,FIX,ADJ]
      dp4 = [FIX,FIX,FIX,FIX,ADJ,FIX,FIX,FIX,ADJ,FIX]
      dp5 = [ADJ,ADJ,ADJ,ADJ,ADJ,FIX,ADJ,ADJ,ADJ,FIX]

      #[f,p,kvg,iter,corp,covp,covr,stdresid,Z,r2] = leasqr[x,y,p,F,stol,niter,wt,dp,dFdp,options}]  #template--arguments in brackets are optional
      #print("\nAdjusting parameters: "),println(parnames[findall(dp1!=0)])  #01/30/2018
      print("\nAdjusting parameters: "),println(parnames[findall(x->x!=0, dp1)])  #11/13/2018
            #println("XYZC:",size(XYZC))  #for debugging
            #println("y:",size(y))
            #println("par1:",size(par1))
            #println("wts:",size(wts))

      #(f,pbest,kvg,iter,corp,covp,covr,stdresid,Z,r2) = leasqr( XYZC,y,par1,P1modelS1e,stol,niter,wts,dp1,
                #dfdp,options )
      println("par1:", par1,size(par1))
      (f,pbest,kvg,iter) = leasqr( XYZC,y,par1,P1modelS1e,stol,niter,wts,dp1,dfdp,options )
      println("iter: ", iter)
      #println("pbest: ", pbest)
      stol2 = 0.1*stol
      #print("\nAdjusting parameters: "),println(parnames[findall(dp2!=0)])
      print("\nAdjusting parameters: "),println(parnames[findall(x->x!=0, dp2)])  #11/13/2018
      (f,pbest,kvg,iter) = leasqr( XYZC,y,pbest,P1modelS1e,stol2,niter,wts,dp2,dfdp,options )
      #print("\nAdjusting parameters:"),println(parnames[findall(dp3)!=0])
      print("\nAdjusting parameters: "),println(parnames[findall(x->x!=0, dp3)])  #11/13/2018
      (f,pbest,kvg,iter) = leasqr( XYZC,y,pbest,P1modelS1e,stol2,niter,wts,dp3,dfdp,options )
      #print("\nAdjusting parameters:"),println(parnames[findall(dp4!=0)])
      print("\nAdjusting parameters: "),println(parnames[findall(x->x!=0, dp4)])  #11/13/2018
      (f,pbest,kvg,iter) = leasqr( XYZC,y,pbest,P1modelS1e,stol2,niter,wts,dp4,dfdp,options )
      #stol3 = stol
      stol3 = stol2
      print("\nAdjusting parameters:"),println(parnames[findall(dp5!=0)])
      print("\nAdjusting parameters: "),println(parnames[findall(x->x!=0, dp5)])  #11/13/2018
      (f,pbest,kvg,iter) = leasqr( XYZC,y,pbest,P1modelS1e,stol3,niter,wts,dp5,dfdp,options )


      kvg     # ==1 if convergence
      #iter    #number of iterations used
      #pbest   #the solution

      println("\nPlane offsets:")
      @printf( "max:%6.3f  min:%6.3f  std:%6.3f  mean:%6.3f\n", maximum(f),minimum(f),std(f),mean(f) )
      S_plane_offset = @sprintf("\nmax:%6.3f  min:%6.3f  std:%6.3f  mean:%6.3f\n",maximum(f),minimum(f),
                    std(f),mean(f) )
      fid = open( string(dirin,"Stats_plane_offsets_S",string(sensor),".txt"), "w" )
      #fprintf( fid,"\nmax:#6.3f  min():#6.3f  std():#6.3f  mean():#6.3f\n",maximum(f),minimum(f),std(f),mean(f) )
      @printf(fid,"%s", S_plane_offset)
      close(fid)

      if sensor==2
         sy = pbest[1]
         sz = pbest[2]
         sv = pbest[3]
         sw = pbest[4]
         gv = pbest[5]
         ccdu = pbest[9]

      else
         sx = pbest[1]
         sz = pbest[2]
         su = pbest[3]
         sw = pbest[4]
         gu = pbest[5]
         ccdv = pbest[9]
      end
      flen = pbest[6]
      x0 = pbest[7]
      ccdw = pbest[8]
      tglass = pbest[10]
   end

   if pars=="JR"
      # =============block 1===========
      if sensor==1
         cdir = [0.971507580 0.077863974 -0.223853129];
         sdir = [0.002846865 0.999994134 -0.001904655]
         gnorm = [0.164473,0.001410,0.986381]

         ccdpos = [425.314655890,0,31.682514443]
         glasspos = [430.589865053,0,-3.865713472]
         #tglass = 2.512399653;
      end

      # =============block 2===========
      if sensor==2
         cdir = [0.003173440 -0.999991372 0.002680596]
         sdir = [0.999998644 -0.001611621 0.000337768]
         gnorm = [-0.000379,-0.025646,0.999671]

         ccdpos = [0,4.976807809,28.870150146]
         glasspos = [0,-9.275043692,-4.367500317]
         #tglass = 2.549200184
      end

      # =============block 3===========
      if sensor==3
         cdir = [-0.974842983 0.000482230 -0.222892186]
         sdir = [0.003831187 0.999991876 -0.001252998]
         gnorm = [-0.243215,0.002147,0.969970]

         ccdpos = [-436.720745069,0,31.354116181]
         glasspos = [-442.037348878,0,-4.326149559]
         #tglass = 2.564262174
      end
      #Convert JR parameters to DT parameters.
      (sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, gu,gv,gw, tglass,nglass, mmPerPixel,
                flen, x0) = P1convertJR2DT[cdir,sdir,gnorm,ccdpos,glasspos,tglass,nglass,mmPerPixel]
   end

   svec = [su,sv,sw]/norm([su,sv,sw])
   gvec = [gu,gv,gw]/norm([gu,gv,gw])
   ccdvec = [ccdu,ccdv,ccdw]/norm([ccdu,ccdv,ccdw])

   @printf( "\nslit[x,y,z]   = [ %5.2f  %5.2f  %5.2f ]\n", sx, sy, sz)
   @printf( "slit[vx,vy,vz]   = [ %8.6f  %8.6f  %8.6f ]\n", svec[1],svec[2],svec[3] )
   @printf( "glass[vx,vy,vz]  = [ %8.6f  %8.6f  %8.6f ]\n", gvec[1],gvec[2],gvec[3] )
   @printf( "ccd[x,y,z]    = [ %5.2f  %5.2f  %5.2f ]\n", ccdx, ccdy, ccdz )
   @printf( "ccd[vx,vy,vz]    = [ %8.6f  %8.6f  %8.6f ]\n", ccdvec[1],ccdvec[2],ccdvec[3] )
   @printf( "flen   = %5.2f\n",flen)
   @printf( "x0     = %7.2f [%smm]\n",x0,string(x0*mmPerPixel)[1:6])
   @printf( "tglass = %6.4f\n",tglass)
   @printf( "gri  %7.4f\n",nglass)

   fout1 = dirin * "_ParametersS" * string(sensornum) * ".txt"
   fid1 = open(fout1,"w")
   @printf( fid1, "S/N %s Run%s Block%s\n", string(serialnum),string(runnum),string(sensornum) )
   @printf( fid1, "gl   %15.12f  %15.12f  %15.12f\n", sx, sy, sz)
   @printf( fid1, "sd   %14.12f  %14.12f  %14.12f\n", svec[1],svec[2],svec[3] )
   @printf( fid1, "gn   %14.12f  %14.12f  %14.12f ]\n", gvec[1],gvec[2],gvec[3] )
   @printf( fid1, "ccd  %15.12f  %15.12f  %15.12f ]\n", ccdx, ccdy, ccdz )
   @printf( fid1, "ccdd %14.12f  %14.12f  %14.12f ]\n", ccdvec[1],ccdvec[2],ccdvec[3] )
   @printf( fid1, "x0   %17.12f [%smm]\n",x0,string(x0*mmPerPixel))
   @printf( fid1, "gt   %14.12f\n",tglass)
   @printf( fid1, "flen %15.12f\n",flen)
   @printf( fid1, "gri  %7.4f\n",nglass)
   close(fid1)

   #************************************************************
   # Calculate centroid corrections using mirror image targeting
   #************************************************************
   @printf("\nCalculating centroid corrections...\n")

   #Find two points along the slit.
   mag = norm([su,sv,sw])
   p1 = [sx, sy, sz] - 10*[su, sv, sw]/mag
   p2 = [sx, sy, sz] + 10*[su, sv, sw]/mag

   npts = size(XYZC)[1]  #rows
   #npts = 2;  #for testing ***********************************

   s = zeros(npts)  #to initialize
   v = zeros(npts)  #to initialize
   for i=1:npts

      x = XYZC[i,1]
      y = XYZC[i,2]
      z = XYZC[i,3]

      #Set the initial target point.
      ptarg = [x, y, z]

      kvg = 0;  #remains zero if convergence fails
            centroid = 0  #added 9/10/2018
            ipt = zeros(3)  #added 9/10/2018
      for k=1:miter

         #Find the intersection of the plane defined by p1,p2,ptarg and the line defined by the CCD.
         (a,b,c,d) = P1planethru3pts( ptarg[1],ptarg[2],ptarg[3], p1[1],p1[2],p1[3],p2[1],p2[2],p2[3] )
         pn = [a,b,c]/norm([a,b,c]);  #normal to plane
         ppt = [p1[1],p1[2],p1[3]];  #a point on the plane
         ldir = [ccdu,ccdv,ccdw]/norm([ccdu,ccdv,ccdw]);  #direction of the CCD
         lpt = [ccdx,ccdy,ccdz];  #a point on the CCD
         ipt = P1intsectlineplane(ppt,pn,lpt,ldir);  #the intersection point
         #p4 = [ccdx,ccdy,ccdz]
         #p5 = p4 + 10*[ccdu,ccdv,ccdw]
         #ipt[1:3] = P1lineplane[p1,p2,ptarg,p4,p5]
         #if norm(ipt1 - ipt)>.000001
         #   i
         #   ipt1
         #   ipt
         #end

         #Find the centroid corresponding to ipt.
         #centroid = norm(ipt - lpt)/mmPerPixel + x0;
         #mag = norm([ccdu,ccdv,ccdw])
         #xc = ccdx + (centroid - x0)*mmPerPixel*ccdu/mag
         if sensor==1 || sensor==3
            centroid = ( ipt[1] - ccdx )/( mmPerPixel*ccdu/norm([ccdu,ccdv,ccdw]) ) + x0
         else  #sensor==2
            centroid = ( ipt[2] - ccdy )/( mmPerPixel*ccdv/norm([ccdu,ccdv,ccdw]) ) + x0
         end

         #Find the exit plane corresponding to the targeting centroid.
         #[a,b,c,d] = P1c2plane_f[x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass,nglass, mmPerPixel]
         #Changed to version _e 4/28/2017 for test purposes
         (a,b,c,d) = P1c2plane_e(x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, snell, gu,gv,gw,
                    tglass,nglass, mmPerPixel)
         offset = P1distancefromplane( x,y,z, a,b,c,d )
         errvec = -offset*[a, b, c]/norm([a, b, c])

   #### vnorm = [a,b,c]/norm([a,b,c])
   #### projpt = P1intsectlineplane[[p1[1],p1[2],p1[3]],vnorm,[x,y,z],vnorm]
   #### errvec = projpt - [x,y,z]


         #### #Find a ray , based on the target point and the centroid CCD-point, that passes through the slit [ignoring glass]
         #### vi = (ptarg - ipt)/norm(ptarg - ipt); #the direction from the centroid point to the target point
         #### vnorm = [gu,gv,gw]/norm([gu,gv,gw])
         #### #point = P1intsectlineplane[ppt,pn,lpt,ldir]
         #### pin = P1intsectlineplane[[sx,sy,sz],vnorm,ipt,vi];  #intersection of the incident ray with the glass
         ####
         #### #Find the location of the exit point on the glass.
         #### pback = [sx,sy,sz] + tglass*(-vnorm);  #a point on the back [exit()] surface of the glass
         #### vrefrac = P1Snell[vi,vnorm,1.000,nglass]
         #### pout = P1intsectlineplane[pback,-vnorm,pin,vrefrac]
         #### #pout = pin + vrefrac*tglass/abs(dot(vrefrac,vnorm))
         #### vout = vi
         #### errvec = P1pointlineoffset[[x,y,z],pout,vout]
         #### offset = norm(errvec)
         ####

         #if i==88 || i==89
         #   i
         #   ptarg
         #   x0
         #   centroid
         #   [a,b,c,d]
         #   offset
         #   errvec
         #end

         if norm(errvec)<tol
            kvg = 1
            break
         end

         #Set a new target point.
         ptarg = ptarg - errvec
      end
      if kvg==1  #if mirror-image()-targeting converged
         #Calculate and save the centroid correction.
         s[i] = centroid - XYZC[i,4];  #centroid correction in pixels
         #Calculate and save the cosine of the ray along the direction of the slit.
         vray = ([x,y,z] - ipt)/norm([x,y,z] - ipt);  #direction from CCD intercept to CMM point
         v[i] = CNORM*dot( vray, [su,sv,sw]/norm([su,sv,sw]) )
         if RBF_on==0
             s[i] = 0;  #all corrections become zero
         end
      else
         s[i] = 0
         printf["mirror image targeting failed to converge for point %s\n",string(i)]
         #XYZC[i,:]
      end
   end
   Ex = zeros(npts)
   Ey = zeros(npts)
   Ez = zeros(npts)
   xs = zeros(npts)
   ys = zeros(npts)
   zs = zeros(npts)
   for i=1:npts
         #Calculate the error after centroid correction.
         centroid = XYZC[i,4] + s[i]
         centroid = XYZC[i,4];  #added 5/22/2018 [to calculate planar offset before centroid correction]
         #[a,b,c,d] = P1c2plane_f[x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass,nglass, mmPerPixel]
         #Changed to version _e 4/28/2017 for test purposes:
         (a,b,c,d) = P1c2plane_e(x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, snell, gu,gv,gw,
                tglass,nglass, mmPerPixel)
         offset = P1distancefromplane( XYZC[i,1],XYZC[i,2],XYZC[i,3], a,b,c,d )
         errvec = -offset*[a, b, c]/norm([a, b, c])
   #printf[ "#4u #10.3f #10.3f #10.3f #10.3f #8.4f  #9.6f\n", i,XYZC[i,1],XYZC[i,2],XYZC[i,3],XYZC[i,4],s[i]*mmPerPixel,offset]
         Ex[i] = 0;       #added 5/22/2018
         Ey[i] = 0;       #added 5/22/2018
         Ez[i] = -offset; #added 5/22/2018
         xs[i] = XYZC[i,1];  #added 5/30/2018
         ys[i] = XYZC[i,2];  #added 5/30/2018
         zs[i] = XYZC[i,3];  #added 5/30/2018
   end

   #fflush[stdout[]]
   flush(stdout)  #for Julia 10/30/2018

   stats = [maximum(s) minimum(s) std(s) mean(s)]*mmPerPixel
println("stats: ", stats)
#println("s: ",s)
#error("deliberate stop")  #uncomment for debugging

   #******************************************
   # Write the centroid corrections to a file.
   #******************************************
      #tic()
      fout = dirin * "_RBFmapS" * string(sensornum) * ".txt"
      fid = open(fout,"w")
      @printf(fid,"S/N %s Run%s Block%s\n", string(serialnum),string(runnum),string(sensornum))
        #first line is a header
      npts
      Map = zeros(miter,3)  #added 9/10/2018
      k = 1  #added 9/10/2018
      for i=1:npts
         Map[k,1] = (XYZC[i,4] - x0)*mmPerPixel
         Map[k,2] = v[i]
         Map[k,3] = s[i]*mmPerPixel
         @printf(fid,"%12.4f %12.8f %12.8f\n", Map[k,1],Map[k,2],Map[k,3])
      end
      close(fid)
      #toc()

   #fflush[stdout[]]
   flush(stdout)  #for Julia 9/10/2018
   if pass==3
      break   #no further filtering
   end

   #********************************
   # Filter base calibration points.
   #********************************

   if Tukey1==1
     #Find outliers using Tukey's method.
     (Irej,Ikeep) = P1removeOutliers(s,filter1,Tukey1)
   else
     stddev = std(s)
     println("Filtering by standard deviation")
     #Re-calculate the standard deviation after removal of the worst outliers.
     #sf = s[find(s.<=stddev*filter1)]; #filtering out the worst outliers
     sf = s[findall(s.<=stddev*filter1)]; #filtering out the worst outliers  #11/14/2018
     stddev = std(sf);
     Irej = findall(s.>stddev*filter1)   #returns a row vector of indices  #11/14/2018
     Irej = findall(s.>stddev*filter1)   #returns a row vector of indices  #11/14/2018
     #Ikeep = find(s.<=stddev*filter1)
     Ikeep = findall(s.<=stddev*filter1)  #11/14/2018
   end


   stddev = std(s)
   sizes = size(s)
        println(size(Irej))
   nr = size(Irej)[1]  #nr = columns(Irej) 9/11/2018
   nk = size(Ikeep)[1]  #nk = columns(Ikeep) 9/11/2018
   if nr>0
      sigma = s[Irej]/stddev
      R = XYZC[Irej,1:4]
      #fid = fopen(strcat(dirin,"_RemovedS",string(sensor),".txt"), "a")
      fid = open(dirin * "_RemovedS" * string(sensor) * ".txt", "a")
      for j=1:nr
         @printf(fid, "point #%u %6.2f sigma\n", Irej[j],sigma[j])
         @printf(fid, "%12.6f %12.6f %12.6f %12.6f\n",R[j,1],R[j,2],R[j,3],R[j,4])
      end
      close(fid)
      ####Ikeep = find(s<=stddev*filter1)
      Tmp = XYZC[Ikeep,1:4]
      #clear XYZC
      XYZC = Tmp
      S = s[Ikeep]
      #clear s
      s = S
      println("Sensor #" * string(sensor) * ", starting pass " * string(pass+1) * " after filtering...")
   else
      break  #no need for another pass
   end
end   #end of pass

#for profiling the time used by the N top functions()
#profile off
#N = 20
#profshow(N)
#dirout = dirin;  #added 7/3/2018
#profexport[strcat(dirout,"profile()/")]
#profile clear()
#error("deliberate stop()")

#=
#Plot the plane offsets.
println("Plotting planar offset error for S" * string(sensor) * "...");
titl = "S" * string(sensor) * ": Plane Offsets\n";  #added 5/22/2018
#titl = titl * "\n" * S_plane_offset;               #added 5/24/2018
xlab = "";                                             #added 5/22/2018
figure(5+sensor);                                      #added 5/22/2018
#quiver_plot_planar_offsets[Ex,Ey,Ez,XYZC[:,1],XYZC[:,2],XYZC[:,3],titl,xlab];  #added 5/22/2018
quiver_plot_planar_offsets(Ex,Ey,Ez,xs,ys,zs,titl,xlab); #rev. 5/30/2018
#Create a .png file of the plane offsets.
println("Creating PNG file for planar offset error()...");
dirout = dirin;  #added 5/24/2018
fname = strcat("Plane_offsets_S',string(sensornum),'.png");
pfile = strcat(dirout,fname);
print(pfile,"-S400,300");
=#
  titl = "S" * string(sensor) * ": Plane Offsets\n";  #added 5/22/2018
  xlab = ""
  figure()
  quiver_plot_planar_offsets( Ex, Ey, Ez, XYZC[:,1], XYZC[:,2], XYZC[:,3], titl  )  #added 3/24/2019

  ### savefig( sbdir * string(serialnum) * "/" * string(runnum) * "/" * "S" * string(sensornum) * "_planar_offsets.png" )

  #savefig( "C:/MFG/" * "_planar_offsets.png" )
  return (sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, gu,gv,gw, tglass,nglass, mmPerPixel,
        flen, x0, #=ck, Ctrs=#) #removed ck, Ctrs, 9/13/2018
end

#*****************************************************************

function P1blockcal_test(serialnum,runnum,sensornum,slitcen=0.)

global P1c2plane_e
global P1calcRBFcoeffs
global sensorbartype
global zoffset
global tglass
global tdustglass

#RBF parameters
global rbfname
global ep
global omega

#P1c2plane_e = P1c2plane_e_orig
#P1calcRBFcoeffs = P1calcRBFcoeffs_orig
sensorbartype = "880mm"
zoffset = 0
tglass = 2.5
tdustglass = 1.5
rbfname = "imq"  #RBF type ('imq', 'tps', 'gauss")
ep = 8           #RBF shape parameter
omega = 0.2      #RBF smoothing parameter


sbdir = "C:/MFG/"
###serialnum = 627845
##
###runnum = 461
###sensornum = 2
###slitcen = 0.  #11/14/2018, added as argument
(sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, gu,gv,gw, tglass,nglass, mmPerPixel, flen, x0,
        #=ck, Ctrs=#) = P1blockcal(sbdir,serialnum,runnum,sensornum,slitcen) #removed ck, Ctrs, 9/13/2018

if sensornum==2
   return sy  #these returns added 11/21/2018
else
   return sx  #these returns added 11/21/2018
end

end

#********************************************************

function P1modelS1e_test( sensornum=1, snell_on=1 )
#Julia version 11/1/2018

	global mmPerPixel
	global nglass
	global snell

	global P1c2plane_e #handle to the function, added 2/15/2018
	global sensor  #(1,2, or 3)

	global sx
	global sy
	global sz

	global su
	global sv

	global ccdu
	global ccdv

	global gw

	mmPerPixel = 0.014;
	nglass = 1.5098;
	snell = snell_on;
	P1c2plane_e = P1c2plane_e
	sensor = sensornum;
	if sensornum == 2
		sx = 0.  #default for S2
		su = 1.0
		ccdv = 1.0
		gw = 1.0
	else
	    sy = -97.056  #for 627486 Run11 S1
	    sv = 0.99976
	    ccdu = 0.9999
	    gw = 0.9997
	end

	#627486 Run11 S1:

	p = [430.291 -97.056 0.00169 0.0216 0.00998 33.454 430.669 0.00780 -0.000180 4.100]
	XYZC = [494.999900  -215.999200  -646.002100   119.478700]

	offset = P1modelS1e( XYZC, p )
	println("offset: ", offset)
@printf(" %20.15f\n", offset[1])

end

#*******************************************************************

function P1blockcalS1S2S3(serialnum,runnum)
#Calibrate the three sensor blocks.
#DJT--11/21/2018

global sensorbartype
global tdustglass

   sensornum = 2
   xS2 = 0.
   yS2 = P1blockcal_test(serialnum,runnum,sensornum,xS2)

   sensornum = 1
   xS1 = P1blockcal_test(serialnum,runnum,sensornum,yS2)

   sensornum = 3
   xS2 = P1blockcal_test(serialnum,runnum,sensornum,yS2)

   sensornum = 2
   xS2 = (xS1 + xS2)/2
   P1blockcal_test(serialnum,runnum,sensornum,xS2)

end

#******************************************************************

end  #end of module jBlockCal
