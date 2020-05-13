function P1Snell_test()
%Test function P1Snell

global tglass = 4.1
n1 = 1.0000
n2 = 1.5098

vinc = [1. 1. 1.];
vinc = vinc/norm(vinc)
vnorm = [-1. 0. 0.]
vrefrac = P1Snell(vinc,vnorm,n1,n2)

export P1Snell_test
end