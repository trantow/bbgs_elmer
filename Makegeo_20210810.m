% Matlab Script to create the mesh of Bering Glacier
% F. Gillet-Chaulet, O. Gagliardini, May 2011
% Edited for Bering Glacier by Thomas Trantow

% create bering.geo
clear;

lc_out=400; % size of element in the plane

A=dlmread('contour_BBGS_v12_UTM.dat');

fid1=fopen('BBGS_400_cv12.geo','w');
fprintf(fid1,'Mesh.Algorithm=5; \n'); % delaunay algorithm

As=size(A,1);

np=0;
for ii=1:As
    np=np+1;
    fprintf(fid1,'Point(%g)={%14.7e,%14.7e,0.0,%g}; \n',np,A(ii,1),A(ii,2),lc_out);
end

fprintf(fid1,'Spline(1)={');
for ii=1:As
  fprintf(fid1,'%g,',ii);
end
fprintf(fid1,'%g}; \n',1);

fprintf(fid1,'Line Loop(2)={1}; \n');
fprintf(fid1,'Plane Surface(3) = {2}; \n');
fprintf(fid1,'Physical Line(4) = {1}; \n');
fprintf(fid1,'Physical Surface(5) = {3}; \n');

fclose(fid1)

% create teterousse.msh using gmsh
%!gmsh bering.geo -1 -2 

% convert teterousse.gmsh in an Elmer type mesh
%!ElmerGrid 14 2 bering.msh -autoclean

% Extrude vertically the mesh (1m thick)
%!./ExtrudeMesh bering bering_ex 3 1 1 0 0 0 0 

% Deform vertically using the surface and bedrock DEM
% Input data are in mesh_input.dat
%!./MshGlacierDEM

% Make a .ep to visualize in ElmerPost the mesh
%!ElmerGrid 2 3 bering_ex
