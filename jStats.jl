module jStats

using Statistics, LinearAlgebra, Printf
export Stats, stats

mutable struct Stats
    xmean
    ymean
    zmean
    xrms
    yrms
    zrms
    xmax
    ymax
    zmax
    diagmean
    diagrms
    diagmax
end #endstruct

function Stats(errvec)  #returns Stats(object)

		n = size(errvec)[1]  #rows
		diag = zeros(n,1)
		for i=1:n
		    diag[i] = norm(errvec[i,1:3])
		end
		diagmax = maximum(diag)
		## diagstd = std(diag)
    diagrms = sqrt(sum(diag.^2)/n)
		diagmean = mean(diag)
		xmax = maximum(abs.(errvec[:,1]))
		ymax = maximum(abs.(errvec[:,2]))
		zmax = maximum(abs.(errvec[:,3]))
		## xstd = std(errvec[:,1])
		## ystd = std(errvec[:,2])
		## zstd = std(errvec[:,3])
		xmean = mean(errvec[:,1])
		ymean = mean(errvec[:,2])
		zmean = mean(errvec[:,3])

  	xsumsq = sum(errvec[:,1].^2)
  	xrms = sqrt(xsumsq/n)
  	ysumsq = sum(errvec[:,2].^2)
  	yrms = sqrt(ysumsq/n)
  	zsumsq = sum(errvec[:,3].^2)
  	zrms = sqrt(zsumsq/n)

  	return Stats(xmean,ymean,zmean,xrms,yrms,zrms,xmax,ymax,zmax,diagmean,diagrms,diagmax)
end #endfunction

function printon(s::jStats.Stats, heading="", fid=stdout)
   @printf("---------------------------------------------------------\n")
   @printf(fid, "%s\n", heading)
   @printf("---------------------------------------------------------\n")
	 @printf(fid, "  xmean =%7.3f   xrms =%7.3f   xmax =%7.3f\n", s.xmean, s.xrms, s.xmax)
   @printf(fid, "  ymean =%7.3f   yrms =%7.3f   ymax =%7.3f\n", s.ymean, s.yrms, s.ymax)
   @printf(fid, "  zmean =%7.3f   zrms =%7.3f   zmax =%7.3f\n", s.zmean, s.zrms, s.zmax)
   @printf(fid, "  diagmean =%7.3f   diagrms =%7.3f   diagmax =%7.3f\n", s.diagmean, s.diagrms, s.diagmax)
   @printf("---------------------------------------------------------\n")
end

function stats(s::jStats.Stats)  #returns Tuple
  means = [s.xmean, s.ymean, s.zmean, s.diagmean]
  rmss = [s.xrms, s.yrms, s.zrms, s.diagrms]
  maxs = [s.xmax, s.ymax, s.zmax, s.diagmax]
  return (means, rmss, maxs)
end

end #endmodule
