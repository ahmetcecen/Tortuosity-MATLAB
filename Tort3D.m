function [TD, varargout]=Tort3D(Alpha)
% 3D Dimensional Tortuosity Calculation in the Z-Direction. A by-product of ongoing computational
% materials science research at MINED@Gatech. (http://mined.gatech.edu/)
%
% [TD]=Tort3D(Alpha) finds the tortuosity distribution of a properly
% weighted 3D dataset Alpha. Every voxel in the dataset Alpha must have a
% value representing the tortuosity of the phase that the voxel
% represents. TD is a list of calculated tortuosity values. 
%
% [TD Tormap]=Tort3D(Alpha) will also generate a non-unique visualization 
% of the calculated paths. 
%
% The current formulation searches shortest paths from every point at the entry
% to the entire exit surface (not every point), though this is easy enough to customize.
%
% Copyright (c) 2015, Ahmet Cecen  -  All rights reserved.

% --- SETUP TOOLS ---
check=dir('dijkstra.mexw64');
if size(check,1)==0
	try
		mex -largeArrayDims dijkstra.cpp
	catch
		fprintf('Please Install a Compatible C++ Compiler (Like Visual Studio Compilers for Windows)')
	end
end

% --- FRAMEWORK ---

% Capture Relevant Grid Properties
[nx,ny,nz] = size(Alpha);

% Generating Volume Boundary Indices
% Outer YZ Boundary Plane 
G=nx:nx:nx*ny*nz;
Gxplus=G';
% Outer XZ Boundary Plane 
G=1:nz;G=G-1;G=G.*nx*ny;G=repmat(G,nx,1);G1(1:nx*nz)=G(1:nx*nz);
G=1:nx;G=repmat(G,nz,1);G=G';G2(1:nx*nz)=G(1:nx*nz);G2=G2+(ny-1)*nx;
Gyplus=G1'+G2';
% Outer XY Boundary Plane 
G=(nz-1)*nx*ny+1:nx*ny*nz;
Gzplus=G';

% --- ADJACENCY TRANSFORM ---

% Assemble Diagonals
xplus=(Alpha+circshift(Alpha,[ -1 0 0 ]))/2.*Alpha.*circshift(Alpha,[ -1 0 0 ])./Alpha./circshift(Alpha,[ -1 0 0 ]);xplus=xplus(:);xplus=[xplus;0];xplus(isnan(xplus))=0;
yplus=(Alpha+circshift(Alpha,[ 0 -1 0 ]))/2.*Alpha.*circshift(Alpha,[ 0 -1 0 ])./Alpha./circshift(Alpha,[ 0 -1 0 ]);yplus=yplus(:);yplus=[yplus;0];yplus(isnan(yplus))=0;
zplus=(Alpha+circshift(Alpha,[ 0 0 -1 ]))/2.*Alpha.*circshift(Alpha,[ 0 0 -1 ])./Alpha./circshift(Alpha,[ 0 0 -1 ]);zplus=zplus(:);zplus=[zplus;0];zplus(isnan(zplus))=0;

% Applying Boundary Conditions
xplus(Gxplus)=0;
yplus(Gyplus)=0;
zplus(Gzplus)=0;

% Form a Pseudo-Banded Matrix
xplus=circshift(xplus,1);
yplus=circshift(yplus,nx);
zplus=circshift(zplus,nx*ny);
PB=[xplus,yplus,zplus];BI=[1,nx,nx*ny];

% Store the Pseudo-Banded Matrix as a Sparse Diagonal Matrix
Adj=spdiags(PB,BI,nx*ny*nz+1,nx*ny*nz+1);

% Modify The Adjacency Matrix for the Surface Node
Adj(nx*ny*(nz-1)+1:nx*ny*nz+1,nx*ny*nz+1)=1;

% Code Specific Arrangement
Adj=Adj+Adj';

%--- SOLVE FOR TORTUOSITY --- Uses Third Party Tool Dijkstra.cpp

[dist, path]=dijkstra(Adj,nx*ny*nz+1);

% --- PROCESS RESULTS ---

% Tortuosity Distribution
Tordata=(dist(1:nx*ny)-1)/(nz-1);
Tordata(isinf(Tordata))=[];
TD=Tordata;

if nargout==2

	% Visualize Shortest Paths (Paths Are Not Unique!)
	pathc=0;
	for i=1:nx*ny
		if path(i)~=1
			pathc=pathc+1;
			pp=i;
			pc=0;
			plist=[];
			while pp~=nx*ny*nz+1
				pp=path(pp);
				pc=pc+1;
				plist(pc)=path(pp);
			end
			Paths{pathc}=plist;
		end
	end

	Tormap=zeros(nx,ny,nz);
	for i=1:length(Paths)
		current=Paths{i};
		current(end-1:end)=[];
		Tormap(current)=1;
	end
	
	varargout{1}=Tormap;
	
end




