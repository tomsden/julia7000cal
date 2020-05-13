module jSnell

using LinearAlgebra

export P1Snell, P1Snell_test

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
#println("r: ", r)
#println("c: ", c)
#c = abs(c)
vrefrac = r*vinc + (r*c - sqrt(1.0 - r*r*(1.0 - c*c)))*vnorm
vrefrac = vrefrac/norm(vrefrac)  #added 2019/6/25
return vrefrac
end
end #end module

### #=
### crossprod = cross(vnorm,vinc)
### vrefrac2 = r*(cross(vnorm,-crossprod)) - vnorm*sqrt( 1 - r*r*dot(crossprod,crossprod) )
###
### pivot = -cross(vinc,-vnorm)
###
### sini = norm(pivot);   #sine of the incident angle()
### sinr = (n1/n2)*sini;  #sine of the refracted angle()
###
### s = sinr
### c = sqrt(1 - s*s)
### pivot = pivot/norm(pivot)
### ux = pivot[1]
### uy = pivot[2]
### uz = pivot[3]
###
### #3D rotation matrix
### R[1,1] = c + ux*ux*(1 - c)
### R[1,2] = ux*uy*(1 - c) - uz*s
### R[1,3] = ux*uz*(1 - c) + uy*s
###
### R[2,1] = uy*ux*(1 - c) + uz*s
### R[2,2] = c + uy*uy*(1 - c)
### R[2,3] = uy*uz*(1 - c) - ux*s
###
### R[3,1] = uz*ux*(1 - c) - uy*s
### #R[3,1] = uy*ux*(1 - c) - uy*s;  #error in Mikahail,et. al.
### R[3,2] = uz*uy*(1 - c) + ux*s
### R[3,3] = c + uz*uz*(1 - c)
###
### vrefrac3 = R*(-vnorm')
###
### printf[ "Snell-1: #18.15f #18.15f #18.15f\n",vrefrac[:]];
### printf[ "Snell-2: #18.15f #18.15f #18.15f\n",vrefrac2[:]]
### printf[ "Snell-3: #18.15f #18.15f #18.15f\n",vrefrac3[:]]
###
###
### fid = fopen("C:/MFG/SnellAlgorithmTests.txt","w")
### fprintf(fid, "Snell-1: #18.15f #18.15f #18.15f\n",vrefrac[:])
### fprintf(fid, "Snell-2: #18.15f #18.15f #18.15f\n",vrefrac2[:])
### fprintf(fid, "Snell-3: #18.15f #18.15f #18.15f\n",vrefrac3[:])
### fclose(fid)
### =#
### end
###
### function P1Snell_test()
### #Test function P1Snell
###
### global tglass = 4.1
### n1 = 1.0000
### n2 = 1.5098
###
### vinc = [1. 1. 1.];
### vinc = vinc/norm(vinc)
### vnorm = [-1. 0. 0.]
### vrefrac = P1Snell(vinc,vnorm,n1,n2)
###
### end
###
### export P1Snell, P1Snell_test
### end #end module
