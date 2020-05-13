module jWriteHTML

using Printf
export writeAccTable2, write_planar_offsets

function writeAccTable2(dirin, dirout, means, rmss, maxs, title, seg, comment="")
#D. Toms, December 10, 2007
#July 2016--modification of "writeAccTable(...)" for RBF results presentation
#June 2019--conversion to julia
if seg==1	#write first segment
	#get the first line from Cal32Out2.dat
	fcal = open(dirin * "Cal32out.dat", "r");
	serialN_date = readline(fcal);
	close(fcal);
	#
	fid = open(dirout * "CalResults2.htm", "w");
	@printf(fid,"<html>\n<head>\n</head>\n");
	@printf(fid,"<body>\n");
	#@printf(fid,"<center><H1>3D Creator Calibration Results</H1></center>\n");
	#@printf(fid,"<center><H1>%s</H1>\n",serialN_date);
	#heading = sscanf(dirin,"[^1234567890]%s","C");
	@printf(fid,"<center><H1>%s</H1><H2>(%s)</H2>\n",dirin,comment);
	@printf(fid,"   <img alt=\"RBFsurface plot\" src=\"RBFsurface_S1.png\" width=400>");
	@printf(fid,"	<img alt=\"RBFsurface plot\" src=\"RBFsurface_S2.png\" width=400>");
	@printf(fid,"   <img alt=\"RBFsurface plot\" src=\"RBFsurface_S3.png\" width=400>\n");
	@printf(fid,"</center>\n>")

	@printf(fid,"<center><table><tr>\n")
	@printf(fid,"  <td><table bgcolor=\"white\" bgborder=1 cellpadding=5>\n");
	@printf(fid,"    <caption valign=\"top\"><strong>%s</strong></caption>\n",title);
	@printf(fid,"    <tr><th> </th><th>mean</th><th>rms</th><th>max</th></tr>\n");
	@printf(fid,"    <tr><th>x</th><td>%8.3f</td><td>%8.3f</td><td>%8.3f</td></tr>\n",means[1],rmss[1],maxs[1]);
	@printf(fid,"    <tr><th>y</th><td>%8.3f</td><td>%8.3f</td><td>%8.3f</td></tr>\n",means[2],rmss[2],maxs[2]);
	@printf(fid,"    <tr><th>z</th><td>%8.3f</td><td>%8.3f</td><td>%8.3f</td></tr>\n",means[3],rmss[3],maxs[3]);
	@printf(fid,"    <tr><th>Euclidean</th><td>%8.3f</td><td>%8.3f</td><td>%8.3f</td></tr>\n",means[4],rmss[4],maxs[4]);
	@printf(fid,"  </table></td>\n");
	close(fid);
end;
if seg==2 #append table
	fid = open(dirout * "CalResults2.htm", "a");
	@printf(fid,"  <td><table bgcolor=\"white\" bgborder=1 cellpadding=5>\n");
	@printf(fid,"    <caption valign=\"top\"><strong>%s</strong></caption>\n",title);
	@printf(fid,"    <tr><th> </th><th>mean</th><th>rms</th><th>max</th></tr>\n");
	@printf(fid,"    <tr><th>x</th><td>%8.3f</td><td>%8.3f</td><td>%8.3f</td></tr>\n",means[1],rmss[1],maxs[1]);
	@printf(fid,"    <tr><th>y</th><td>%8.3f</td><td>%8.3f</td><td>%8.3f</td></tr>\n",means[2],rmss[2],maxs[2]);
	@printf(fid,"    <tr><th>z</th><td>%8.3f</td><td>%8.3f</td><td>%8.3f</td></tr>\n",means[3],rmss[3],maxs[3]);
	@printf(fid,"    <tr><th>Euclidean</th><td>%8.3f</td><td>%8.3f</td><td>%8.3f</td></tr>\n",means[4],rmss[4],maxs[4]);
	@printf(fid,"  </table></center></td>\n");
	close(fid);
end;
if seg==3	#append end of HTML
	fid = open(dirout * "CalResults2.htm", "a");
	@printf(fid,"</tr></table></center>\n");
	@printf(fid,"<center>\n");
	@printf(fid,"<img alt=\"RBF histogram x\" src=\"RBFhist_x.png\"  width=400>\n");
	@printf(fid,"<img alt=\"RBF histogram y\" src=\"RBFhist_y.png\"  width=400>\n");
	@printf(fid,"<img alt=\"RBF histogram z\" src=\"RBFhist_z.png\"  width=400>\n");
	@printf(fid,"</center></body>\n</html>");
	close(fid);
end;
if seg==4   #omit histograms
  fid = open(dirout * "CalResults2.htm", "a");
	@printf(fid,"</tr></table></center>\n");
	@printf(fid,"<center>\n");
	#Write the length-error accuracy chart.
	@printf(fid,"   <img alt=\"AccuracyChart_L\" src=\"AccuracyChart_L.png\" width=400>");
	@printf(fid,"</center></body>\n</html>");
	close(fid);
end
end #endfunction

function write_planar_offsets(dirin, seg=1, comment="")
#Write planar offset plots to an HTML file
#DJT--June 2019

if seg==1	#write first segment
	#get the first line from Cal32Out2.dat
	fcal = open(dirin * "Cal32out.dat", "r");
	serialN_date = readline(fcal);
	close(fcal);
	#
	dirout = dirin
	fid = open(dirout * "Plane_offsets.htm", "w");
	@printf(fid,"<html>\n<head>\n</head>\n");
	@printf(fid,"<body>\n");
	#@printf(fid,"<center><H1>3D Creator Calibration Results</H1></center>\n");
	#@printf(fid,"<center><H1>%s</H1>\n",serialN_date);
	#heading = sscanf(dirin,"[^1234567890]%s","C");
	@printf(fid,"<center><H1>%s</H1><H2>(%s)</H2>\n",dirin,comment);
	@printf(fid,"   <img alt=\"Planar_offsets plot\" src=\"Plane_offsets_S1.png\" width=400>");
	@printf(fid,"	  <img alt=\"Planar_offsets plot\" src=\"Plane_offsets_S2.png\" width=400>");
	@printf(fid,"   <img alt=\"Planar_offsets plot\" src=\"Plane_offsets_S3.png\" width=400>\n");
	@printf(fid,"</center>\n>")
#=
	@printf(fid,"<center><table><tr>\n")
	@printf(fid,"  <td><table bgcolor=\"white\" bgborder=1 cellpadding=5>\n");
	@printf(fid,"    <caption valign=\"top\"><strong>%s</strong></caption>\n",title);
	@printf(fid,"    <tr><th> </th><th>mean</th><th>rms</th><th>max</th></tr>\n");
	@printf(fid,"    <tr><th>x</th><td>%8.3f</td><td>%8.3f</td><td>%8.3f</td></tr>\n",means[1],rmss[1],maxs[1]);
	@printf(fid,"    <tr><th>y</th><td>%8.3f</td><td>%8.3f</td><td>%8.3f</td></tr>\n",means[2],rmss[2],maxs[2]);
	@printf(fid,"    <tr><th>z</th><td>%8.3f</td><td>%8.3f</td><td>%8.3f</td></tr>\n",means[3],rmss[3],maxs[3]);
	@printf(fid,"    <tr><th>Euclidean</th><td>%8.3f</td><td>%8.3f</td><td>%8.3f</td></tr>\n",means[4],rmss[4],maxs[4]);
	@printf(fid,"  </table></td>\n");
=#
	close(fid);
end;
#=
if seg==2 #append table
	fid = open(dirout * "CalResults2.htm", "a");
	@printf(fid,"  <td><table bgcolor=\"white\" bgborder=1 cellpadding=5>\n");
	@printf(fid,"    <caption valign=\"top\"><strong>%s</strong></caption>\n",title);
	@printf(fid,"    <tr><th> </th><th>mean</th><th>rms</th><th>max</th></tr>\n");
	@printf(fid,"    <tr><th>x</th><td>%8.3f</td><td>%8.3f</td><td>%8.3f</td></tr>\n",means[1],rmss[1],maxs[1]);
	@printf(fid,"    <tr><th>y</th><td>%8.3f</td><td>%8.3f</td><td>%8.3f</td></tr>\n",means[2],rmss[2],maxs[2]);
	@printf(fid,"    <tr><th>z</th><td>%8.3f</td><td>%8.3f</td><td>%8.3f</td></tr>\n",means[3],rmss[3],maxs[3]);
	@printf(fid,"    <tr><th>Euclidean</th><td>%8.3f</td><td>%8.3f</td><td>%8.3f</td></tr>\n",means[4],rmss[4],maxs[4]);
	@printf(fid,"  </table></center></td>\n");
	close(fid);
end;
if seg==3	#append end of HTML
	fid = open(dirout * "CalResults2.htm", "a");
	@printf(fid,"</tr></table></center>\n");
	@printf(fid,"<center>\n");
	@printf(fid,"<img alt=\"RBF histogram x\" src=\"RBFhist_x.png\"  width=400>\n");
	@printf(fid,"<img alt=\"RBF histogram y\" src=\"RBFhist_y.png\"  width=400>\n");
	@printf(fid,"<img alt=\"RBF histogram z\" src=\"RBFhist_z.png\"  width=400>\n");
	@printf(fid,"</center></body>\n</html>");
	close(fid);
end;
if seg==4   #omit histograms
  fid = open(dirout * "CalResults2.htm", "a");
	@printf(fid,"</tr></table></center>\n");
	@printf(fid,"<center>\n");
	@printf(fid,"</center></body>\n</html>");
	close(fid);
end
=#
end #endfunction
end #endmodule
